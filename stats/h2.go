// Package stats 实现本流程的三类统计推断：
//   - h2.go：基于汇总统计与 LD 分数的截距校正遗传力（LDSC 风格）；
//   - meta.go：多队列固定/随机效应元分析与重叠检测；
//   - prs.go：个体多基因风险评分构建与评估。
//
// 所有面向用户的错误都来自 variant / numwrap / ingest 的可区分错误类型，
// 不直接透出第三方库异常。
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

// H2BoundaryPolicy 决定遗传力点估计落在可解释区间 [0, h2Max] 之外时的行为。
type H2BoundaryPolicy uint8

const (
	// H2ReportRaw 不改动点估计，仅在结果上置标志（推荐用于研究审计）。
	H2ReportRaw H2BoundaryPolicy = iota
	// H2Clamp 截断到 [0, h2Max]，并置 Clamped 标志；CI 仍基于未截断分布报告。
	H2Clamp
	// H2Reject 视为不可解释结果，返回带 Reason 的错误。
	H2Reject
)

// H2Options 控制遗传力估计。
type H2Options struct {
	// NBar：SNP 级 N 缺失时的平均样本量；0 表示必须逐 SNP 有 N。
	NBar float64
	// Weighted=true 时使用 LDSC 风格权重 1/(2*l2+const) 的简化形式
	// 1/max(2*l2, eps)（固定截距迭代权重留给后续版本，显式不暴露）。
	Weighted bool
	// FreeIntercept=true：自由估计截距（样本膨胀校正）；false：固定为 1。
	FreeIntercept bool
	// Boundary 控制越界策略；h2Max 缺省 1。
	Boundary H2BoundaryPolicy
	H2Max    float64
	// MinVariants：可估计的最小 SNP 数（默认 20）。
	MinVariants int
	// Workers：分片并行度；1 为串行。任意取值都产生逐位一致结果。
	Workers int
	// Align 控制与 LD 参考的等位基因对齐。
	Align variant.AlignConfig
}

// H2Flag 以可区分位标记描述估计过程中发生的情况。
type H2Flag uint8

const (
	H2FlagNone             H2Flag = 0
	H2FlagNegativeRaw      H2Flag = 1 << iota // 原始点估计为负
	H2FlagAboveMaxRaw                         // 原始点估计超过 h2Max
	H2FlagClamped                             // 已按策略截断
	H2FlagInterceptFixed                      // 截距被固定为 1
	H2FlagLowPower                            // CI 过宽（半宽 > h2Max），不可解释
	H2FlagAmbiguousDropped                    // 因链模糊被拒绝的位点数>0
)

func (f H2Flag) String() string {
	if f == H2FlagNone {
		return "none"
	}
	var parts []string
	if f&H2FlagNegativeRaw != 0 {
		parts = append(parts, "negative_raw")
	}
	if f&H2FlagAboveMaxRaw != 0 {
		parts = append(parts, "above_max_raw")
	}
	if f&H2FlagClamped != 0 {
		parts = append(parts, "clamped")
	}
	if f&H2FlagInterceptFixed != 0 {
		parts = append(parts, "intercept_fixed")
	}
	if f&H2FlagLowPower != 0 {
		parts = append(parts, "low_power")
	}
	if f&H2FlagAmbiguousDropped != 0 {
		parts = append(parts, "ambiguous_dropped")
	}
	return joinParts(parts)
}

func joinParts(parts []string) string {
	out := ""
	for i, p := range parts {
		if i > 0 {
			out += ","
		}
		out += p
	}
	return out
}

// H2Result 是遗传力估计输出。
type H2Result struct {
	H2Raw       float64 // 未截断点估计
	H2          float64 // 策略处理后的报告值
	SE          float64
	CI95        [2]float64
	Intercept   float64
	InterceptSE float64
	NSnps       int
	Flags       H2Flag
	Rejected    []ingest.RowRef // 对齐/质量拒绝明细
	// Components 按分片给出充分统计量，用于审计并行=串行。
	Components []H2Component
	Diag       string
}

// H2Component 是单个分片的可加性充分统计量。
type H2Component struct {
	Key   string
	N     int
	SumXX float64
	SumXY float64
	SumW  float64
}

// matchedSNP 是内连接 + 对齐后的逐 SNP 回归数据。
type matchedSNP struct {
	key   variant.MatchKey
	x     float64 // l2
	y     float64 // chi2
	w     float64
	n     float64
	chrom string
}

// EstimateH2 用汇总统计与 LD 分数估计 SNP 遗传力。
//
// 模型（bulit-SLDSc 简化版）：
//
//	E[chi2] = intercept + (Nbar * h2 / M) * l2
//
// 因而 slope = Nbar*h2/M，h2 = slope*M/Nbar（M=纳入 SNP 数）。
func EstimateH2(ctx context.Context, sum *ingest.SummaryTable, ld *ingest.LDScoreTable, opts H2Options) (*H2Result, error) {
	if sum == nil || ld == nil {
		return nil, fmt.Errorf("EstimateH2: nil input table")
	}
	if opts.H2Max == 0 {
		opts.H2Max = 1
	}
	if opts.MinVariants == 0 {
		opts.MinVariants = 20
	}
	if opts.Align == (variant.AlignConfig{}) {
		opts.Align = variant.DefaultAlignConfig()
	}
	if opts.Align.Ambiguous == 0 && opts.Align.Flip == 0 {
		opts.Align = variant.DefaultAlignConfig()
	}

	matched, rejected, nAmb, err := joinAndAlign(sum, ld, opts)
	if err != nil {
		return nil, err
	}
	if len(matched) < opts.MinVariants {
		return nil, &h2EstimationError{Reason: "insufficient_variants",
			Detail: fmt.Sprintf("matched=%d < min=%d", len(matched), opts.MinVariants)}
	}

	// 按 key 排序，保证输入顺序不影响逐 SNP 回归的累加次序。
	sort.Slice(matched, func(i, j int) bool { return matched[i].key < matched[j].key })

	workers := normalizeWorkers(opts.Workers, len(matched))
	comps, err := mapShards(ctx, len(matched), workers, func(start, end int) (H2Component, error) {
		c := H2Component{Key: fmt.Sprintf("shard[%d:%d]", start, end)}
		var kxx, kxy, kw float64 // Kahan 补偿
		var cxx, cxy, cw float64
		add := func(acc, comp *float64, v float64) {
			yv := v - *comp
			t := *acc + yv
			*comp = (t - *acc) - yv
			*acc = t
		}
		for i := start; i < end; i++ {
			m := matched[i]
			add(&kxx, &cxx, m.w*m.x*m.x)
			add(&kxy, &cxy, m.w*m.x*m.y)
			add(&kw, &cw, m.w)
			c.N++
		}
		c.SumXX, c.SumXY, c.SumW = kxx, kxy, kw
		return c, nil
	})
	if err != nil {
		return nil, err
	}

	// 归并严格按分片顺序做 Kahan 累加（与串行分片完全一致）。
	var sxx, sxy, sw float64
	var cxx, cxy, cw float64
	add := func(acc, comp *float64, v float64) {
		yv := v - *comp
		t := *acc + yv
		*comp = (t - *acc) - yv
		*acc = t
	}
	for _, c := range comps {
		add(&sxx, &cxx, c.SumXX)
		add(&sxy, &cxy, c.SumXY)
		add(&sw, &cw, c.SumW)
	}
	if sxx <= 0 {
		return nil, &h2EstimationError{Reason: "singular_design", Detail: "weighted sum of l2^2 is zero"}
	}

	// 全样本 WLS（复用 numwrap，确保方差/截距与单遍求解一致）。
	xs := make([]float64, len(matched))
	ys := make([]float64, len(matched))
	ws := make([]float64, len(matched))
	for i, m := range matched {
		xs[i], ys[i], ws[i] = m.x, m.y, m.w
	}
	wlsOpts := numwrap.RegressionOptions{Weights: ws, InterceptFree: opts.FreeIntercept}
	if !opts.FreeIntercept {
		wlsOpts.FixedIntercept = 1
	}
	fit, nerr := numwrap.WLS(xs, ys, wlsOpts)
	if nerr != nil {
		return nil, fmt.Errorf("h2 regression failed: %w", nerr)
	}

	var intercept, slope, slopeSE, interceptSE float64
	flags := H2FlagNone
	if opts.FreeIntercept {
		intercept, slope = fit.Coefficients[0], fit.Coefficients[1]
		interceptSE, slopeSE = fit.StdErr[0], fit.StdErr[1]
	} else {
		intercept, slope = 1, fit.Coefficients[0]
		interceptSE, slopeSE = 0, fit.StdErr[0]
		flags |= H2FlagInterceptFixed
	}

	// Nbar：优先逐 SNP N（加权平均），否则用 opts.NBar。
	nbar := opts.NBar
	var nSum, nW float64
	for _, m := range matched {
		nSum += m.w * m.n
		nW += m.w
	}
	if nW > 0 && nSum/nW > 0 {
		nbar = nSum / nW
	}
	if nbar <= 0 {
		return nil, &h2EstimationError{Reason: "missing_sample_size",
			Detail: "no per-SNP N and H2Options.NBar <= 0"}
	}

	M := float64(len(matched))
	h2raw := slope * M / nbar
	h2se := slopeSE * M / nbar

	z95, _ := numwrap.NormalCDFInv(0.975)
	res := &H2Result{
		H2Raw:       h2raw,
		H2:          h2raw,
		SE:          h2se,
		CI95:        [2]float64{h2raw - z95*h2se, h2raw + z95*h2se},
		Intercept:   intercept,
		InterceptSE: interceptSE,
		NSnps:       len(matched),
		Components:  comps,
		Rejected:    rejected,
	}

	if h2raw < 0 {
		flags |= H2FlagNegativeRaw
	}
	if h2raw > opts.H2Max {
		flags |= H2FlagAboveMaxRaw
	}
	// 低功效：CI 半宽超过 h2Max（区间宽到无法解释）
	if z95*h2se > opts.H2Max {
		flags |= H2FlagLowPower
	}
	if nAmb > 0 {
		flags |= H2FlagAmbiguousDropped
	}

	switch opts.Boundary {
	case H2ReportRaw:
		// 保留原值与标志
	case H2Clamp:
		if h2raw < 0 {
			res.H2 = 0
			flags |= H2FlagClamped
		} else if h2raw > opts.H2Max {
			res.H2 = opts.H2Max
			flags |= H2FlagClamped
		}
	case H2Reject:
		if h2raw < 0 {
			return nil, &h2EstimationError{Reason: "negative_heritability",
				Detail: fmt.Sprintf("raw h2=%.4g is negative and policy=reject", h2raw)}
		}
		if h2raw > opts.H2Max {
			return nil, &h2EstimationError{Reason: "heritability_above_max",
				Detail: fmt.Sprintf("raw h2=%.4g > h2Max=%.4g", h2raw, opts.H2Max)}
		}
	}
	res.Flags = flags
	res.Diag = fmt.Sprintf(
		"model=chi2~intercept+slope*l2; M=%d Nbar=%.1f intercept=%.4f(%.4f) slope=%.6f h2raw=%.6f se=%.6f ci=[%.4f,%.4f] flags=%s; rejected=%d (ambiguous=%d)",
		len(matched), nbar, intercept, interceptSE, slope, h2raw, h2se,
		res.CI95[0], res.CI95[1], flags, len(rejected), nAmb)
	_ = sxy
	_ = sw
	return res, nil
}

// h2EstimationError 是估计层的可区分错误（与输入拒绝、数值错误区分）。
type h2EstimationError struct {
	Reason string
	Detail string
}

func (e *h2EstimationError) Error() string {
	return "h2 estimation: " + e.Reason + ": " + e.Detail
}

// joinAndAlign 以 variant.Key 内连接汇总统计与 LD 分数，并对齐等位基因。
//
// LD 参考作为方向基准；返回的 y=chi2 与方向无关，x=l2 同理，
// 但对齐仍被强制执行，以便把模糊位点/无法对齐位点显式剔除并记录。
func joinAndAlign(sum *ingest.SummaryTable, ld *ingest.LDScoreTable, opts H2Options) ([]matchedSNP, []ingest.RowRef, int, error) {
	ldByKey := make(map[variant.MatchKey]*ingest.LDScoreRecord, len(ld.Records))
	for i := range ld.Records {
		ldByKey[ld.Records[i].V.Key()] = &ld.Records[i]
	}
	var out []matchedSNP
	var rejected []ingest.RowRef
	nAmb := 0

	// 汇总侧重复已在解析时拒绝；LD 侧同。sum 行顺序决定遍历顺序，
	// 最终结果在 EstimateH2 中按 key 排序，因此这里的追加顺序无关结论。
	for _, rec := range sum.Records {
		lr, ok := ldByKey[rec.V.Key()]
		if !ok {
			rejected = append(rejected, ingest.RowRef{
				ID: rec.V.ID, Reason: variant.ReasonAlleleMismatch,
				Detail: "no LD score for variant key " + string(rec.V.Key())})
			continue
		}
		al, err := variant.Align(lr.V, rec.V, opts.Align)
		if err != nil {
			r, _ := variant.AsReject(err)
			rr := ingest.RowRef{ID: rec.V.ID, Detail: err.Error()}
			if r != nil {
				rr.Reason = r.Reason
				rr.Line = r.Line
				if r.Reason == variant.ReasonAmbiguousStrand {
					nAmb++
				}
			} else {
				rr.Reason = variant.ReasonAlleleMismatch
			}
			rejected = append(rejected, rr)
			continue
		}
		if al.Ambiguous {
			nAmb++
			rejected = append(rejected, ingest.RowRef{
				ID: rec.V.ID, Reason: variant.ReasonAmbiguousStrand,
				Detail: "kept-policy ambiguous site excluded from h2 regression"})
			continue
		}
		chi2 := rec.Z * rec.Z
		if math.IsNaN(chi2) || math.IsInf(chi2, 0) {
			rejected = append(rejected, ingest.RowRef{
				ID: rec.V.ID, Reason: variant.ReasonInvalidNumeric, Detail: "chi2 non-finite"})
			continue
		}
		w := 1.0
		if opts.Weighted {
			w = 1 / math.Max(2*lr.L2, 1e-12)
		}
		nEff := rec.N
		if lr.N > 0 {
			nEff = lr.N
		}
		out = append(out, matchedSNP{
			key: rec.V.Key(), x: lr.L2, y: chi2, w: w, n: nEff,
			chrom: rec.V.Chrom,
		})
	}
	return out, rejected, nAmb, nil
}
