package stats

import (
	"context"
	"fmt"
	"math"
	"sort"

	"github.com/lechandonga/statistical-genetics-go/ingest"
	"github.com/lechandonga/statistical-genetics-go/numwrap"
	"github.com/lechandonga/statistical-genetics-go/variant"
)

// MetaModel 选择元分析模型。
type MetaModel uint8

const (
	MetaFixed      MetaModel = iota // 固定效应逆方差加权
	MetaRandomDL                    // DerSimonian–Laird 随机效应
	MetaRandomREML                  // 预留：当前后端不支持，显式拒绝
)

// OverlapPolicy 决定检测到样本重叠/重复队列时的行为。
type OverlapPolicy uint8

const (
	// OverlapFlag 仅标记（默认）：输出 Overlap 明细与 Flags；
	// 重叠队列对逐 SNP 贡献的权重照常给出，但 MetaResult 带 overlap 标志，
	// 不会把重叠误称为独立证据（解释层责任显式化）。
	OverlapFlag OverlapPolicy = iota
	// OverlapDrop 每组重叠中保留最先出现的队列，其余整体排除。
	OverlapDrop
	// OverlapReject 检测到任意重叠即整体拒绝。
	OverlapReject
)

// Cohort 是参与元分析的单队列输入。
type Cohort struct {
	Name      string
	Samples   int
	SampleIDs []string // 若提供，交集即重叠的硬证据
	Table     *ingest.SummaryTable
}

// MetaOptions 控制元分析。
type MetaOptions struct {
	Model   MetaModel
	Overlap OverlapPolicy
	Align   variant.AlignConfig
	Workers int
	// MinCohorts：每个变体至少需要的有效队列数（默认 2），否则该变体排除。
	MinCohorts int
	// Backend 选择本次元分析使用的数值后端；空 => numwrap 默认。
	// 后端在入口解析一次并贯穿全程，运行期间切换进程默认后端不影响本次任务。
	Backend numwrap.Backend
	// CrossCheck=true 时同一输入分别用 stdlib 与 gonum 各完整计算一遍：
	// 数值须落在 numwrap 声明容差内，离散判定（纳入/排除/标志）须完全一致；
	// 对不上即返回 *numwrap.BackendMismatchError，绝不挑选一侧结果。
	CrossCheck bool
}

// OverlapEvidence 描述一对队列的样本重叠证据。
type OverlapEvidence struct {
	CohortA, CohortB string
	SharedSamples    int
	Jaccard          float64
	Duplicate        bool // 样本集合完全一致
}

// MetaVariant 是单个变体的合并结论。
type MetaVariant struct {
	V             *variant.Variant
	Beta          float64
	SE            float64
	CI95          [2]float64
	P             float64
	WeightFixed   float64
	WeightRandom  float64
	Tau2          float64
	Q             float64
	Df            int
	PQ            float64
	I2            float64
	NCohorts      int
	Excluded      bool
	ExcludeReason variant.Reason
	PerCohort     []MetaCohortLine
	// IncludedCohorts/ExclusionLines 是 PerCohort 的便捷汇总，供输出直接展示。
	IncludedCohorts []string
	ExclusionLines  []MetaCohortLine
}

// MetaCohortLine 记录每个队列对该变体的贡献与对齐依据。
type MetaCohortLine struct {
	Cohort        string
	SourceKind    string // "beta_se"（完整效应）或 "z_only"（仅 Z+N 换算）
	BetaAligned   float64
	SE            float64
	Z             float64 // 该队列报告的原始 Z
	N             float64 // 换算所用样本量
	Weight        float64 // 实际用于合并的逆方差权重
	WeightBasis   string  // 权重/换算口径的人类可读说明
	Sign          float64
	Relation      variant.FlipRelation
	DirectionNote string // 方向依据（对齐解释；参考队列为 reference direction）
	Excluded      bool
	Reason        variant.Reason
}

// MetaResult 是元分析输出。
type MetaResult struct {
	Variants        []MetaVariant
	Overlap         []OverlapEvidence
	Flags           []string
	ExcludedCohorts []string
	Backend         numwrap.Backend
	CrossChecked    bool
	Diag            string
}

// metaError 是元分析层的可区分错误。
type metaError struct{ what, detail string }

func (e *metaError) Error() string { return "meta-analysis: " + e.what + ": " + e.detail }

// RunMeta 执行多队列元分析。
//
// 方向约定：每个变体以“第一个能成功对齐的有效队列”为参考方向，
// 后续队列 beta 按 AlignResult.Sign 变号；对齐失败的队列出现在
// PerCohort（Excluded=true + 原因），绝不静默丢弃。
func RunMeta(ctx context.Context, cohorts []Cohort, opts MetaOptions) (*MetaResult, error) {
	if len(cohorts) == 0 {
		return nil, &metaError{"no_cohorts", "at least one cohort required"}
	}
	if opts.Model == MetaRandomREML {
		return nil, &metaError{"unsupported_model",
			"REML random effects are not supported by the current numerical backend"}
	}
	if opts.Align == (variant.AlignConfig{}) {
		opts.Align = variant.DefaultAlignConfig()
	}
	if opts.MinCohorts == 0 {
		opts.MinCohorts = 2
	}
	for i := range cohorts {
		if cohorts[i].Table == nil {
			return nil, &metaError{"nil_table", "cohort " + cohorts[i].Name + " has nil table"}
		}
	}

	// 1) 重叠检测（仅基于显式样本 ID；不基于效应相关性猜测）
	overlap := detectOverlap(cohorts)
	excludedCohorts := map[string]bool{}
	var flags []string
	if len(overlap) > 0 {
		flags = append(flags, "sample_overlap_detected")
		switch opts.Overlap {
		case OverlapReject:
			return nil, &metaError{"overlap_rejected",
				fmt.Sprintf("%d overlapping cohort pair(s); first=%s/%s shared=%d",
					len(overlap), overlap[0].CohortA, overlap[0].CohortB, overlap[0].SharedSamples)}
		case OverlapDrop:
			flags = append(flags, "overlapping_cohorts_dropped")
			for _, ev := range overlap {
				// 保留对中先出现的，排除后者（已排除则跳过）
				if !excludedCohorts[ev.CohortA] {
					excludedCohorts[ev.CohortB] = true
				}
			}
		case OverlapFlag:
			if anyDuplicate(overlap) {
				flags = append(flags, "duplicate_cohort_detected")
			}
		}
	}
	// 2) 收集全部变体键（按 key 排序，保证确定性）
	keyOrder, perCohortIdx := collectKeys(cohorts, excludedCohorts)

	workers := normalizeWorkers(opts.Workers, len(keyOrder))
	backend := opts.Backend
	if backend == "" {
		backend = numwrap.DefaultBackend()
	}
	if !numwrap.Available(backend) {
		return nil, &numwrap.Error{Kind: numwrap.KindUnsupportedOption, Op: "RunMeta",
			Detail: "backend not compiled into this build: " + string(backend)}
	}

	// 第一阶段（后端无关）：逐变体对齐方向并换算贡献。
	prepared, err := mapShards(ctx, len(keyOrder), workers, func(start, end int) (prepareBlock, error) {
		blk := prepareBlock{lo: start, hi: end, items: make([]preparedVariant, 0, end-start)}
		for ki := start; ki < end; ki++ {
			blk.items = append(blk.items,
				prepareVariant(keyOrder[ki], cohorts, perCohortIdx, excludedCohorts, opts))
		}
		return blk, nil
	})
	if err != nil {
		return nil, err
	}

	flat := make([]preparedVariant, 0, len(keyOrder))
	for _, blk := range prepared {
		flat = append(flat, blk.items...)
	}

	// 第二阶段：选定后端做数值合并；CrossCheck 时两侧各算一遍并比较。
	primary, err := combineAll(ctx, len(keyOrder), workers, flat, opts, backend)
	if err != nil {
		return nil, err
	}
	if opts.CrossCheck {
		other := numwrap.BackendStd
		if backend == numwrap.BackendStd {
			other = numwrap.BackendGonum
		}
		secondary, err := combineAll(ctx, len(keyOrder), workers, flat, opts, other)
		if err != nil {
			return nil, err
		}
		if err := crossCheckMeta(primary, secondary); err != nil {
			return nil, err
		}
	}

	res := &MetaResult{Overlap: overlap, Backend: backend, CrossChecked: opts.CrossCheck}
	res.Variants = primary

	// 稳定排序的排除队列名
	var excl []string
	seen := map[string]bool{}
	for _, ev := range overlap {
		for _, name := range []string{ev.CohortB} {
			if excludedCohorts[name] && !seen[name] {
				seen[name] = true
				excl = append(excl, name)
			}
		}
	}
	sort.Strings(excl)
	res.ExcludedCohorts = excl
	res.Flags = flags
	nExcludedVariants := 0
	for _, v := range res.Variants {
		if v.Excluded {
			nExcludedVariants++
		}
	}
	res.Diag = fmt.Sprintf("cohorts=%d model=%s variantKeys=%d excludedVariants=%d overlapPairs=%d droppedCohorts=%d flags=%v",
		len(cohorts), modelName(opts.Model), len(keyOrder), nExcludedVariants,
		len(overlap), len(excl), flags)
	return res, nil
}

type resultBlock struct {
	lo, hi int
	vars   []MetaVariant
}

// contribution 是一条已对齐到参考方向、可参与数值合并的队列贡献。
type contribution struct {
	cohort      string
	sourceKind  string
	beta, se, z float64
	n, weight   float64
	weightBasis string
}

// preparedVariant 是后端无关的逐变体准备结果。
type preparedVariant struct {
	key      variant.MatchKey
	ref      *variant.Variant
	refCi    int // 方向参考队列下标；-1 表示无可参考队列
	lines    []MetaCohortLine
	contrib  []contribution
	excluded bool
	reason   variant.Reason
}

type prepareBlock struct {
	lo, hi int
	items  []preparedVariant
}

// combineAll 用指定后端对已准备变体做数值合并；
// 与准备阶段使用相同 (n,workers) 分片边界，下标一一对应。
func combineAll(ctx context.Context, n, workers int, prepared []preparedVariant,
	opts MetaOptions, backend numwrap.Backend) ([]MetaVariant, error) {
	blocks, err := mapShards(ctx, n, workers, func(start, end int) (resultBlock, error) {
		blk := resultBlock{lo: start, hi: end, vars: make([]MetaVariant, 0, end-start)}
		for bi := start; bi < end; bi++ {
			blk.vars = append(blk.vars, combineVariant(prepared[bi], opts, backend))
		}
		return blk, nil
	})
	if err != nil {
		return nil, err
	}
	var out []MetaVariant
	for _, blk := range blocks {
		out = append(out, blk.vars...)
	}
	return out, nil
}

func modelName(m MetaModel) string {
	switch m {
	case MetaFixed:
		return "fixed_IVW"
	case MetaRandomDL:
		return "random_DL"
	default:
		return "unknown"
	}
}

func anyDuplicate(evs []OverlapEvidence) bool {
	for _, e := range evs {
		if e.Duplicate {
			return true
		}
	}
	return false
}

// detectOverlap 基于显式样本 ID 集合计算两两交集与 Jaccard。
// 有任一队列未提供样本 ID 时，不对该队列下任何结论（禁止猜测）。
func detectOverlap(cohorts []Cohort) []OverlapEvidence {
	sets := make([]map[string]struct{}, len(cohorts))
	usable := make([]bool, len(cohorts))
	for i, c := range cohorts {
		if len(c.SampleIDs) == 0 {
			continue
		}
		m := make(map[string]struct{}, len(c.SampleIDs))
		for _, id := range c.SampleIDs {
			m[id] = struct{}{}
		}
		sets[i], usable[i] = m, true
	}
	var out []OverlapEvidence
	for i := 0; i < len(cohorts); i++ {
		for j := i + 1; j < len(cohorts); j++ {
			if !usable[i] || !usable[j] {
				continue
			}
			// 在较小集合上遍历
			a, b := sets[i], sets[j]
			shared := 0
			if len(a) > len(b) {
				a, b = b, a
			}
			for id := range a {
				if _, ok := b[id]; ok {
					shared++
				}
			}
			if shared == 0 {
				continue
			}
			union := len(sets[i]) + len(sets[j]) - shared
			jac := float64(shared) / float64(union)
			dup := shared == len(sets[i]) && shared == len(sets[j])
			out = append(out, OverlapEvidence{
				CohortA: cohorts[i].Name, CohortB: cohorts[j].Name,
				SharedSamples: shared, Jaccard: jac, Duplicate: dup})
		}
	}
	return out
}

// collectKeys 返回排序后的全部变体键，以及每队列 key->record 索引。
func collectKeys(cohorts []Cohort, excluded map[string]bool) ([]variant.MatchKey, []map[variant.MatchKey]int) {
	allKeys := map[variant.MatchKey]struct{}{}
	idx := make([]map[variant.MatchKey]int, len(cohorts))
	for ci, c := range cohorts {
		if excluded[c.Name] {
			idx[ci] = map[variant.MatchKey]int{}
			continue
		}
		m := make(map[variant.MatchKey]int, len(c.Table.Records))
		for ri, rec := range c.Table.Records {
			m[rec.V.Key()] = ri
			allKeys[rec.V.Key()] = struct{}{}
		}
		idx[ci] = m
	}
	keys := make([]variant.MatchKey, 0, len(allKeys))
	for k := range allKeys {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })
	return keys, idx
}

// prepareVariant 完成单个变体的后端无关工作：
// 逐队列对齐方向、把 Z-only 队列按统一口径换算为 (beta,se)、
// 记录每条纳入/排除审计行。任何不可用队列都给出可区分原因，
// 绝不静默丢弃，也绝不产生“参与了合并但数值不可解释且无标记”的行。
func prepareVariant(key variant.MatchKey, cohorts []Cohort,
	perCohortIdx []map[variant.MatchKey]int, excludedCohort map[string]bool,
	opts MetaOptions) preparedVariant {

	var ref *variant.Variant
	pv := preparedVariant{key: key}
	pv.refCi = -1

	// 第一遍：确定方向参考（第一个可对齐、非模糊的有效队列）。
	for ci, c := range cohorts {
		if excludedCohort[c.Name] {
			continue
		}
		ri, present := perCohortIdx[ci][key]
		if !present {
			continue
		}
		rec := c.Table.Records[ri]
		if rec.V.IsAmbiguous() && opts.Align.Ambiguous == variant.AmbiguousReject {
			continue
		}
		ref = rec.V
		pv.refCi = ci
		break
	}

	for ci, c := range cohorts {
		ri, present := perCohortIdx[ci][key]
		line := MetaCohortLine{Cohort: c.Name}
		if excludedCohort[c.Name] {
			line.Excluded = true
			line.Reason = variant.Reason("cohort_excluded_overlap")
			line.WeightBasis = "excluded before per-variant combination"
			pv.lines = append(pv.lines, line)
			continue
		}
		if !present {
			line.Excluded = true
			line.Reason = variant.Reason("variant_missing_in_cohort")
			pv.lines = append(pv.lines, line)
			continue
		}
		rec := c.Table.Records[ri]

		// 方向对齐（含模糊位点判定）。
		var sign float64 = 1
		var relation variant.FlipRelation = variant.RelationNone
		if ci == pv.refCi {
			line.DirectionNote = "reference direction cohort (first aligned cohort)"
		} else if ref != nil {
			al, err := variant.Align(ref, rec.V, opts.Align)
			if err != nil {
				r, _ := variant.AsReject(err)
				line.Excluded = true
				if r != nil {
					line.Reason = r.Reason
					line.DirectionNote = r.Detail
				} else {
					line.Reason = variant.ReasonAlleleMismatch
					line.DirectionNote = "alleles cannot be aligned"
				}
				pv.lines = append(pv.lines, line)
				continue
			}
			relation, sign = al.Relation, al.Sign
			line.DirectionNote = al.Explanation
		}
		if rec.V.IsAmbiguous() && opts.Align.Ambiguous == variant.AmbiguousReject {
			line.Relation = relation
			line.Sign = sign
			line.Excluded = true
			line.Reason = variant.ReasonAmbiguousStrand
			line.DirectionNote = "ambiguous A/T or C/G site: sign not identifiable"
			pv.lines = append(pv.lines, line)
			continue
		}

		// 数值来源：完整 (beta,se) 或 Z-only 换算。
		zonly := math.IsNaN(rec.SE) || rec.SE <= 0
		con, why := contributionOf(c.Name, rec, sign, zonly)
		if why != "" {
			line.Relation = relation
			line.Sign = sign
			line.Excluded = true
			line.Reason = why
			if zonly {
				line.SourceKind = "z_only"
				line.Z = rec.Z
				line.N = rec.N
			} else {
				line.SourceKind = "beta_se"
				line.BetaAligned = rec.Beta * sign
				line.SE = rec.SE
			}
			pv.lines = append(pv.lines, line)
			continue
		}

		line.SourceKind = con.sourceKind
		line.BetaAligned = con.beta
		line.SE = con.se
		line.Z = con.z
		line.N = con.n
		line.Weight = con.weight
		line.WeightBasis = con.weightBasis
		line.Relation = relation
		line.Sign = sign
		pv.lines = append(pv.lines, line)
		pv.contrib = append(pv.contrib, con)
	}

	pv.ref = ref
	if ref == nil {
		pv.excluded = true
		pv.reason = variant.Reason("no_alignable_cohort")
	} else if len(pv.contrib) < opts.MinCohorts {
		pv.excluded = true
		pv.reason = variant.Reason("insufficient_cohorts")
	}
	return pv
}

// contributionOf 把一条队列记录换算为已定向的 (beta,se,z,n,weight)。
//
// 统一口径（Z-only 队列默认真正参与，不一律排除）：
//
//	完整队列：直接 beta、SE；w = 1/SE^2。
//	Z-only 队列：p = 2*Phi(-|Z|)；SE = 1/sqrt(2p(1-p)N)；beta = |Z|*SE*sign；
//	             w = 1/SE^2 = 2p(1-p)N。
//
// 仅当缺换算所需的样本量（N<=0），或 Z 本身没法定向
// （NaN/Inf/恰好为 0：符号不存在）时才排除，并给出可区分原因。
func contributionOf(cohort string, rec ingest.SummaryRecord, sign float64, zonly bool) (contribution, variant.Reason) {
	if zonly {
		z := rec.Z
		if math.IsNaN(z) || math.IsInf(z, 0) {
			return contribution{}, variant.Reason("z_not_finite")
		}
		if z == 0 {
			return contribution{}, variant.Reason("z_undirected")
		}
		if !(rec.N > 0) {
			return contribution{}, variant.Reason("missing_sample_size")
		}
		p := 2 * normalSurvival(math.Abs(z))
		info := 2 * p * (1 - p) * rec.N
		if !(info > 0) || math.IsInf(info, 0) {
			return contribution{}, variant.Reason("z_conversion_non_finite")
		}
		se := 1 / math.Sqrt(info)
		// 方向：来源 Z 自带符号（beta/se），对齐时整体乘对齐符号 sign；
		// Z-only 与完整队列因此共用同一方向约定。
		beta := z * se * sign
		if math.IsNaN(beta) || math.IsInf(beta, 0) || se <= 0 {
			return contribution{}, variant.Reason("z_conversion_non_finite")
		}
		return contribution{
			cohort: cohort, sourceKind: "z_only", beta: beta, se: se,
			z: z, n: rec.N, weight: info,
			weightBasis: fmt.Sprintf("z-only: p=2*Phi(-|Z|), se=1/sqrt(2p(1-p)N), N=%.0f", rec.N),
		}, ""
	}
	if math.IsNaN(rec.Beta) || math.IsInf(rec.Beta, 0) ||
		math.IsNaN(rec.SE) || math.IsInf(rec.SE, 0) || rec.SE <= 0 {
		return contribution{}, variant.Reason("invalid_numeric")
	}
	beta := rec.Beta * sign
	w := 1 / (rec.SE * rec.SE)
	return contribution{
		cohort: cohort, sourceKind: "beta_se", beta: beta, se: rec.SE,
		z: rec.Z, n: rec.N, weight: w,
		weightBasis: "inverse-variance w=1/se^2",
	}, ""
}

// combineVariant 用指定后端对已准备变体做数值合并。
func combineVariant(pv preparedVariant, opts MetaOptions, backend numwrap.Backend) MetaVariant {
	mv := MetaVariant{V: pv.ref, PerCohort: pv.lines}
	for i := range pv.lines {
		if pv.lines[i].Excluded {
			mv.ExclusionLines = append(mv.ExclusionLines, pv.lines[i])
		} else {
			mv.IncludedCohorts = append(mv.IncludedCohorts, pv.lines[i].Cohort)
		}
	}
	if pv.excluded {
		mv.Excluded = true
		mv.ExcludeReason = pv.reason
		return mv
	}

	betas := make([]float64, len(pv.contrib))
	ses := make([]float64, len(pv.contrib))
	for i, c := range pv.contrib {
		betas[i], ses[i] = c.beta, c.se
	}

	fit, nerr := numwrap.MetaCombine(betas, ses, numwrap.MetaOptions{
		Backend: backend, Random: opts.Model == MetaRandomDL})
	if nerr != nil {
		// 输入已经过 contributionOf 校验；走到这里属于数值层异常，
		// 变体整体排除并显式标记，绝不留下不可解释的数值。
		mv.Excluded = true
		mv.ExcludeReason = variant.Reason("numerical_backend_error")
		mv.ExclusionLines = append(mv.ExclusionLines, MetaCohortLine{
			Cohort: "<backend>", Excluded: true,
			Reason:      variant.Reason("numerical_backend_error"),
			WeightBasis: nerr.Error()})
		return mv
	}

	mv.Beta, mv.SE, mv.WeightFixed = fit.FixedBeta, fit.FixedSE, fit.FixedWeight
	mv.Q, mv.Df, mv.I2 = fit.Q, fit.Df, fit.I2
	mv.PQ = chi2Survival(fit.Q, fit.Df)
	mv.NCohorts = len(pv.contrib)
	if opts.Model == MetaRandomDL {
		mv.Tau2 = fit.Tau2
		mv.Beta, mv.SE, mv.WeightRandom = fit.RandomBeta, fit.RandomSE, fit.RandomWeight
	} else {
		mv.WeightRandom = fit.FixedWeight
	}

	z95, _ := numwrap.NormalCDFInvB(0.975, backend)
	mv.CI95 = [2]float64{mv.Beta - z95*mv.SE, mv.Beta + z95*mv.SE}
	mv.P = 2 * normalSurvival(math.Abs(mv.Beta/mv.SE))

	// 最后防线：凡参与合并的结果，任何对外数值都必须有限，
	// 否则整体排除并标记（不允许“参与了却不可解释且无标记”）。
	if !finiteMeta(mv) {
		mv.Excluded = true
		mv.ExcludeReason = variant.Reason("non_finite_combined_result")
	}
	return mv
}

func finiteMeta(mv MetaVariant) bool {
	for _, v := range []float64{mv.Beta, mv.SE, mv.P, mv.CI95[0], mv.CI95[1],
		mv.Q, mv.I2, mv.WeightFixed, mv.WeightRandom} {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return false
		}
	}
	return true
}

// crossCheckMeta 比较两个后端的完整结果：
// 数值在声明容差内一致；离散判定（排除与否、原因、纳入/排除队列集合）
// 必须完全一致。任何不一致都是硬错误。
func crossCheckMeta(a, b []MetaVariant) error {
	const op = "meta"
	if len(a) != len(b) {
		return &numwrap.BackendMismatchError{Op: op,
			Detail: fmt.Sprintf("variant count %d != %d", len(a), len(b))}
	}
	for i := range a {
		id := "<excluded>"
		if a[i].V != nil {
			id = a[i].V.ID
		}
		tag := fmt.Sprintf("variant[%s]", id)
		if err := numwrap.CrossCheckBool(op, tag+".excluded", a[i].Excluded, b[i].Excluded); err != nil {
			return err
		}
		if a[i].ExcludeReason != b[i].ExcludeReason {
			return &numwrap.BackendMismatchError{Op: op,
				Detail: fmt.Sprintf("%s exclude_reason %q != %q", tag, a[i].ExcludeReason, b[i].ExcludeReason)}
		}
		if err := numwrap.CrossCheckInt(op, tag+".ncohorts", a[i].NCohorts, b[i].NCohorts); err != nil {
			return err
		}
		if a[i].Excluded {
			continue
		}
		checks := []struct {
			name string
			x, y float64
		}{
			{"beta", a[i].Beta, b[i].Beta},
			{"se", a[i].SE, b[i].SE},
			{"p", a[i].P, b[i].P},
			{"q", a[i].Q, b[i].Q},
			{"i2", a[i].I2, b[i].I2},
			{"pQ", a[i].PQ, b[i].PQ},
			{"tau2", a[i].Tau2, b[i].Tau2},
			{"ci95lo", a[i].CI95[0], b[i].CI95[0]},
			{"ci95hi", a[i].CI95[1], b[i].CI95[1]},
		}
		for _, ck := range checks {
			if err := numwrap.CrossCheckScalar(op, tag+"."+ck.name, ck.x, ck.y); err != nil {
				return err
			}
		}
		if err := crossCheckCohortLines(op, tag, a[i].PerCohort, b[i].PerCohort); err != nil {
			return err
		}
	}
	return nil
}

func crossCheckCohortLines(op, tag string, a, b []MetaCohortLine) error {
	if len(a) != len(b) {
		return &numwrap.BackendMismatchError{Op: op, Detail: tag + " per-cohort line count differs"}
	}
	for i := range a {
		p := tag + "." + a[i].Cohort
		if a[i].Cohort != b[i].Cohort || a[i].Excluded != b[i].Excluded ||
			a[i].Reason != b[i].Reason || a[i].SourceKind != b[i].SourceKind {
			return &numwrap.BackendMismatchError{Op: op,
				Detail: fmt.Sprintf("%s discrete audit differs: %+v vs %+v", p, a[i], b[i])}
		}
		if a[i].Excluded {
			continue
		}
		if err := numwrap.CrossCheckScalar(op, p+".weight", a[i].Weight, b[i].Weight); err != nil {
			return err
		}
		if err := numwrap.CrossCheckScalar(op, p+".betaAligned", a[i].BetaAligned, b[i].BetaAligned); err != nil {
			return err
		}
		if err := numwrap.CrossCheckScalar(op, p+".se", a[i].SE, b[i].SE); err != nil {
			return err
		}
	}
	return nil
}

// chi2Survival 返回卡方(df) 的上尾概率 P(Chi2_df > x)。
// 使用 Wilson–Hilferty 正态近似（df>=1），仅用于异质性 P 值展示；
// 文档中明确该近似，不宣称精确。
func chi2Survival(x float64, df int) float64 {
	if df <= 0 {
		return math.NaN()
	}
	k := float64(df)
	// WH: (x/k)^(1/3) approx Normal(1-2/(9k), 2/(9k))
	t := math.Pow(x/k, 1.0/3.0)
	mean := 1 - 2/(9*k)
	sd := math.Sqrt(2 / (9 * k))
	z := (t - mean) / sd
	return normalSurvival(z)
}

// normalSurvival 标准正态上尾 P(Z>z)（math.Erfc，确定性）。
func normalSurvival(z float64) float64 {
	return 0.5 * math.Erfc(z/math.Sqrt2)
}
