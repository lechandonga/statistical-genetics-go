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
		beta, se float64
	}
	var contribs []contribution

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
			line.BetaAligned = rec.Beta
			line.SE = rec.SE
			line.Weight = weightOf(rec)
			if rec.V.IsAmbiguous() && opts.Align.Ambiguous == variant.AmbiguousReject {
				line.Excluded = true
				line.Reason = variant.ReasonAmbiguousStrand
				lines = append(lines, line)
				continue
			}
			lines = append(lines, line)
			contribs = append(contribs, contribution{rec.Beta, rec.SE})
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
		betaA := rec.Beta * al.Sign
		line.BetaAligned = betaA
		line.SE = rec.SE
		line.Weight = weightOf(rec)
		lines = append(lines, line)
		contribs = append(contribs, contribution{betaA, rec.SE})
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

	z95, _ := numwrap.NormalCDFInv(0.975)
	mv.CI95 = [2]float64{mv.Beta - z95*mv.SE, mv.Beta + z95*mv.SE}
	mv.P = 2 * normalSurvival(math.Abs(mv.Beta/mv.SE))
	return mv
}

func weightOf(rec ingest.SummaryRecord) float64 {
	if rec.SE == rec.SE && rec.SE > 0 {
		return 1 / (rec.SE * rec.SE)
	}
	// Z-only 队列：无法形成逆方差权重
	return math.NaN()
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
