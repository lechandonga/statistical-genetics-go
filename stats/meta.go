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
	// Backend 指定本次调用使用的数值后端（空=默认；crosscheck=双后端验证）。
	Backend numwrap.Backend
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
}

// ContributionBasis 标明一个队列的效应/SE/权重来自哪种统一口径。
type ContributionBasis string

const (
	// BasisBetaSE：队列直接提供 beta 与 SE，SE 直接采用。
	BasisBetaSE ContributionBasis = "beta_se"
	// BasisZNEAF：队列仅提供 Z 与样本量 N，且 EAF 已知且在 (0,1)：
	// 按统一口径换算
	//
	//	SE(beta) = 1 / sqrt(2 * N * EAF * (1 - EAF))
	//	beta     = Z * SE
	//
	// （加性剂量 0/1/2 的方差 Var(X)=2p(1-p)）。方向即 Z 的符号。
	BasisZNEAF ContributionBasis = "z_n_eaf"
)

// 混合来源队列的可区分排除原因。
const (
	// ReasonZMissingN：Z-only 队列缺少换算所需样本量 N。
	ReasonZMissingN variant.Reason = "z_missing_n"
	// ReasonZUnorientableMAF：Z-only 队列没有可用 EAF（或 EAF 为 0/1，
	// 剂量方差为 0），beta 尺度无法确定，禁止参与合并。
	ReasonZUnorientableMAF variant.Reason = "z_unorientable_maf"
	// ReasonZInvalid：Z-only 队列的 Z 非有限数值。
	ReasonZInvalid variant.Reason = "z_invalid"
	// ReasonContributionNonFinite：换算后的贡献仍为非有限值
	// （防御性检查：绝不允许“参与了合并但结果无法解释且无标记”）。
	ReasonContributionNonFinite variant.Reason = "non_finite_contribution"
)

// MetaCohortLine 记录每个队列对该变体的贡献与对齐依据。
type MetaCohortLine struct {
	Cohort      string
	BetaAligned float64
	SE          float64
	Weight      float64
	Sign        float64
	Relation    variant.FlipRelation
	Excluded    bool
	Reason      variant.Reason
	// N 为该队列在该变体上的样本量（审计用）。
	N float64
	// Z 为原始 Z 分数（beta/SE 或直接提供），方向审计用。
	Z float64
	// EAF 为换算时使用的效应等位基因频率；beta/SE 直接提供时可为 NaN。
	EAF float64
	// Basis 标明贡献口径（beta_se / z_n_eaf）；排除的队列为空。
	Basis ContributionBasis
	// BasisDetail 为人类可读的换算依据。
	BasisDetail string
}

// MetaResult 是元分析输出。
type MetaResult struct {
	Variants        []MetaVariant
	Overlap         []OverlapEvidence
	Flags           []string
	ExcludedCohorts []string
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
	if opts.Backend == numwrap.BackendCrosscheck {
		return crosscheckRun("meta",
			func(o MetaOptions, b numwrap.Backend) MetaOptions { o.Backend = b; return o },
			compareMetaResults,
			func(o MetaOptions) (*MetaResult, error) { return RunMeta(ctx, cohorts, o) },
			opts)
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
	mv, err := mapShards(ctx, len(keyOrder), workers, func(start, end int) (resultBlock, error) {
		blk := resultBlock{lo: start, hi: end, vars: make([]MetaVariant, 0, end-start)}
		for ki := start; ki < end; ki++ {
			key := keyOrder[ki]
			mv := analyzeVariant(key, cohorts, perCohortIdx, excludedCohorts, opts)
			blk.vars = append(blk.vars, mv)
		}
		return blk, nil
	})
	if err != nil {
		return nil, err
	}
	res := &MetaResult{Overlap: overlap}
	for _, blk := range mv { // 严格按分片顺序展开
		res.Variants = append(res.Variants, blk.vars...)
	}

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

func analyzeVariant(key variant.MatchKey, cohorts []Cohort,
	perCohortIdx []map[variant.MatchKey]int, excludedCohort map[string]bool,
	opts MetaOptions) MetaVariant {

	var ref *variant.Variant
	lines := make([]MetaCohortLine, 0, len(cohorts))
	type contribution struct {
		beta, se, n, z, eaf float64
		basis               ContributionBasis
		detail              string
	}
	var contribs []contribution

	// acceptContribution 统一入口：任何参与合并的贡献都必须方向确定
	// 且 beta/se/权重全部有限，否则以可区分原因排除（绝不静默纳入）。
	acceptContribution := func(line MetaCohortLine, c contribution) {
		if c.beta != c.beta || c.se != c.se || c.se <= 0 {
			line.Excluded = true
			line.Reason = ReasonContributionNonFinite
			line.Basis = ""
			lines = append(lines, line)
			return
		}
		w := 1 / (c.se * c.se)
		if w != w || w <= 0 {
			line.Excluded = true
			line.Reason = ReasonContributionNonFinite
			line.Basis = ""
			lines = append(lines, line)
			return
		}
		line.BetaAligned = c.beta
		line.SE = c.se
		line.Weight = w
		line.N, line.Z, line.EAF = c.n, c.z, c.eaf
		line.Basis = c.basis
		line.BasisDetail = c.detail
		lines = append(lines, line)
		contribs = append(contribs, c)
	}

	for ci, c := range cohorts {
		ri, present := perCohortIdx[ci][key]
		line := MetaCohortLine{Cohort: c.Name}
		if excludedCohort[c.Name] {
			line.Excluded = true
			line.Reason = variant.Reason("cohort_excluded_overlap")
			lines = append(lines, line)
			continue
		}
		if !present {
			line.Excluded = true
			line.Reason = variant.Reason("variant_missing_in_cohort")
			lines = append(lines, line)
			continue
		}
		rec := c.Table.Records[ri]
		if ref == nil {
			// 第一个有效队列成为方向基准；但其自身若为链模糊位点，
			// 方向同样不可作为参考——整条变体按策略标记/排除。
			ref = rec.V
			line.Relation = variant.RelationNone
			line.Sign = 1
			if rec.V.IsAmbiguous() && opts.Align.Ambiguous == variant.AmbiguousReject {
				line.Excluded = true
				line.Reason = variant.ReasonAmbiguousStrand
				lines = append(lines, line)
				continue
			}
			beta, se, basis, detail, r := deriveContribution(rec, 1)
			if r != "" {
				line.Excluded = true
				line.Reason = r
				lines = append(lines, line)
				continue
			}
			acceptContribution(line, contribution{
				beta: beta, se: se, n: rec.N, z: rec.Z, eaf: rec.V.EAF,
				basis: basis, detail: detail})
			continue
		}
		al, err := variant.Align(ref, rec.V, opts.Align)
		if err != nil {
			r, _ := variant.AsReject(err)
			line.Excluded = true
			if r != nil {
				line.Reason = r.Reason
			} else {
				line.Reason = variant.ReasonAlleleMismatch
			}
			lines = append(lines, line)
			continue
		}
		line.Relation = al.Relation
		line.Sign = al.Sign
		beta, se, basis, detail, r := deriveContribution(rec, al.Sign)
		if r != "" {
			line.Excluded = true
			line.Reason = r
			lines = append(lines, line)
			continue
		}
		acceptContribution(line, contribution{
			beta: beta, se: se, n: rec.N, z: rec.Z * al.Sign, eaf: rec.V.EAF,
			basis: basis, detail: detail})
	}

	mv := MetaVariant{V: ref, PerCohort: lines}
	if ref == nil || len(contribs) < opts.MinCohorts {
		mv.Excluded = true
		mv.ExcludeReason = variant.Reason("insufficient_cohorts")
		return mv
	}

	betas := make([]float64, len(contribs))
	ses := make([]float64, len(contribs))
	for i, c := range contribs {
		betas[i] = c.beta
		ses[i] = c.se
	}

	fBeta, fSE, wFixed, q, df, i2 := fixedIVW(betas, ses)
	mv.Beta, mv.SE, mv.WeightFixed = fBeta, fSE, wFixed
	mv.Q, mv.Df, mv.I2 = q, df, i2
	mv.PQ = chi2Survival(q, df)
	mv.NCohorts = len(contribs)

	if opts.Model == MetaRandomDL {
		dl := analyzeDL(betas, ses, q, df, wFixed)
		mv.Tau2 = dl.tau2
		mv.Beta, mv.SE, mv.WeightRandom = dl.beta, dl.se, dl.wsum
	} else {
		mv.WeightRandom = wFixed
	}

	z95 := quantile975(opts.Backend)
	mv.CI95 = [2]float64{mv.Beta - z95*mv.SE, mv.Beta + z95*mv.SE}
	mv.P = 2 * normalSurvival(math.Abs(mv.Beta/mv.SE))
	return mv
}

// deriveContribution 把单条汇总统计换算为对齐到参考方向的 (beta, SE)，
// 并返回口径与人类可读依据；无法参与时返回非空 variant.Reason。
//
// 两种来源统一到同一逆方差口径：
//
//  1. beta + SE（EffectBeta/EffectOR）：直接采用（beta 已按 sign 变号）；
//  2. Z-only（EffectZ）：必须有 N；SE 按 2*N*p*(1-p) 反演（p=EAF）。
//     EAF 缺失或为 0/1 时剂量方差为 0、beta 尺度无法定向，返回
//     z_unorientable_maf；N 缺失返回 z_missing_n。
//
// sign 为对齐符号（±1）：对 Z-only 队列，方向直接随 Z 变号。
func deriveContribution(rec ingest.SummaryRecord, sign float64) (beta, se float64,
	basis ContributionBasis, detail string, reason variant.Reason) {

	zOnly := math.IsNaN(rec.SE)
	if !zOnly {
		if rec.SE <= 0 || math.IsInf(rec.SE, 0) {
			return 0, 0, "", "", ReasonContributionNonFinite
		}
		b := rec.Beta * sign
		if math.IsNaN(b) || math.IsInf(b, 0) {
			return 0, 0, "", "", ReasonContributionNonFinite
		}
		return b, rec.SE, BasisBetaSE,
			fmt.Sprintf("beta=%.6g se=%.6g provided directly; sign=%+v from allele alignment",
				rec.Beta, rec.SE, sign), ""
	}

	// Z-only 分支
	if math.IsNaN(rec.Z) || math.IsInf(rec.Z, 0) {
		return 0, 0, "", "", ReasonZInvalid
	}
	if !(rec.N > 0) {
		return 0, 0, "", "", ReasonZMissingN
	}
	p := rec.V.EAF
	if math.IsNaN(p) || p <= 0 || p >= 1 {
		return 0, 0, "", "", ReasonZUnorientableMAF
	}
	vx := 2 * rec.N * p * (1 - p)
	se = 1 / math.Sqrt(vx)
	zAligned := rec.Z * sign
	beta = zAligned * se
	if math.IsNaN(beta) || math.IsInf(beta, 0) || math.IsNaN(se) || se <= 0 {
		return 0, 0, "", "", ReasonContributionNonFinite
	}
	return beta, se, BasisZNEAF,
		fmt.Sprintf("Z=%.6g N=%.0f EAF=%.4f: SE=1/sqrt(2*N*p*(1-p))=%.6g, beta=Z*SE; sign=%+v",
			rec.Z, rec.N, p, se, sign), ""
}

// fixedIVW 返回固定效应 beta、SE、总权重、Cochran Q、自由度与 I2。
// 所有求和均为 Kahan 顺序求和（确定性）。
func fixedIVW(betas, ses []float64) (beta, se, wsum, q float64, df int, i2 float64) {
	type ka struct{ s, c float64 }
	add := func(a *ka, v float64) {
		y := v - a.c
		t := a.s + y
		a.c = (t - a.s) - y
		a.s = t
	}
	var sw, swb ka
	weights := make([]float64, len(betas))
	for i := range betas {
		w := 1 / (ses[i] * ses[i])
		weights[i] = w
		add(&sw, w)
		add(&swb, w*betas[i])
	}
	beta = swb.s / sw.s
	se = math.Sqrt(1 / sw.s)
	var sq ka
	for i := range betas {
		add(&sq, weights[i]*(betas[i]-beta)*(betas[i]-beta))
	}
	q = sq.s
	df = len(betas) - 1
	if q > float64(df) {
		i2 = (q - float64(df)) / q
		if i2 < 0 {
			i2 = 0
		}
	}
	return beta, se, sw.s, q, df, i2
}

// analyzeDL 使用权重数组完整计算 DL tau2 与随机效应合并量。
type dlResult struct {
	beta, se, wsum, tau2 float64
}

func analyzeDL(betas, ses []float64, q float64, df int, wFixed float64) dlResult {
	type ka struct{ s, c float64 }
	add := func(a *ka, v float64) {
		y := v - a.c
		t := a.s + y
		a.c = (t - a.s) - y
		a.s = t
	}
	var sw, sw2 ka
	for i := range betas {
		w := 1 / (ses[i] * ses[i])
		add(&sw, w)
		add(&sw2, w*w)
	}
	tau2 := 0.0
	if q > float64(df) {
		c := sw.s - sw2.s/sw.s
		if c > 0 {
			tau2 = (q - float64(df)) / c
		}
	}
	swR := ka{}
	swbR := ka{}
	for i := range betas {
		wr := 1 / (ses[i]*ses[i] + tau2)
		add(&swR, wr)
		add(&swbR, wr*betas[i])
	}
	_ = wFixed
	return dlResult{beta: swbR.s / swR.s, se: math.Sqrt(1 / swR.s), wsum: swR.s, tau2: tau2}
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
