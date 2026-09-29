package numwrap

// 元分析数值原语（固定效应 IVW + DerSimonian–Laird 随机效应）。
//
// 这些原语在 stdlib 与 gonum 两个后端各有一份独立实现：
// 两侧使用同一组 float64 直接公式、但实现路径独立（stdlib 走
// Kahan 顺序补偿求和；gonum 后端走独立的累加路径，不共用代码），
// 用于发现单一实现里的系统性问题。

import "math"

// MetaResult 是固定/随机效应合并的数值结果。
type MetaResult struct {
	FixedBeta    float64
	FixedSE      float64
	FixedWeight  float64
	Q            float64
	I2           float64
	Tau2         float64
	RandomBeta   float64
	RandomSE     float64
	RandomWeight float64
	Df           int
	Weights      []float64 // 固定效应逆方差权重 w_i = 1/se_i^2（审计用）
}

// MetaOptions 选择后端；空 Backend 使用进程默认。
type MetaOptions struct {
	Backend Backend
	// Random=false 只计算固定效应；true 同时计算 DL 随机效应。
	Random bool
}

// MetaCombine 对已对齐到同一方向的 (beta, se) 序列做固定效应 IVW，
// 并按需计算 DL 随机效应。se 必须全部有限且 >0。
func MetaCombine(betas, ses []float64, opts MetaOptions) (*MetaResult, error) {
	b, err := resolve(opts.Backend, "MetaCombine")
	if err != nil {
		return nil, err
	}
	if len(betas) != len(ses) {
		return nil, &Error{Kind: KindInvalidArgument, Op: "MetaCombine",
			Detail: "beta/se length mismatch"}
	}
	if len(betas) < 2 {
		return nil, &Error{Kind: KindInvalidArgument, Op: "MetaCombine",
			Detail: "at least two contributions required"}
	}
	for i := range ses {
		if math.IsNaN(betas[i]) || math.IsNaN(ses[i]) {
			return nil, &Error{Kind: KindNaN, Op: "MetaCombine", Detail: "NaN in inputs"}
		}
		if math.IsInf(betas[i], 0) || math.IsInf(ses[i], 0) {
			return nil, &Error{Kind: KindInf, Op: "MetaCombine", Detail: "Inf in inputs"}
		}
		if ses[i] <= 0 {
			return nil, &Error{Kind: KindInvalidArgument, Op: "MetaCombine",
				Detail: "se must be positive"}
		}
	}
	switch b {
	case BackendStd:
		return metaCombineImpl(betas, ses, opts.Random), nil
	case BackendGonum:
		return metaCombineGonum(betas, ses, opts.Random), nil
	default:
		return nil, &Error{Kind: KindUnsupportedOption, Op: "MetaCombine",
			Detail: "backend not available: " + string(b)}
	}
}

// metaCombineImpl 是 stdlib 后端：固定顺序 Kahan 补偿求和。
func metaCombineImpl(betas, ses []float64, random bool) *MetaResult {
	type ka struct{ s, c float64 }
	add := func(a *ka, v float64) {
		y := v - a.c
		t := a.s + y
		a.c = (t - a.s) - y
		a.s = t
	}
	n := len(betas)
	weights := make([]float64, n)
	var sw, swb ka
	for i := range betas {
		w := 1 / (ses[i] * ses[i])
		weights[i] = w
		add(&sw, w)
		add(&swb, w*betas[i])
	}
	beta := swb.s / sw.s
	se := math.Sqrt(1 / sw.s)
	var sq ka
	for i := range betas {
		add(&sq, weights[i]*(betas[i]-beta)*(betas[i]-beta))
	}
	df := n - 1
	i2 := 0.0
	if sq.s > float64(df) {
		i2 = (sq.s - float64(df)) / sq.s
		if i2 < 0 {
			i2 = 0
		}
	}
	res := &MetaResult{
		FixedBeta: beta, FixedSE: se, FixedWeight: sw.s,
		Q: sq.s, I2: i2, Df: df, Weights: weights,
	}
	if random {
		var sw2 ka
		for i := range weights {
			add(&sw2, weights[i]*weights[i])
		}
		tau2 := 0.0
		if sq.s > float64(df) {
			denom := sw.s - sw2.s/sw.s
			if denom > 0 {
				tau2 = (sq.s - float64(df)) / denom
			}
		}
		var swr, swbr ka
		for i := range betas {
			wr := 1 / (ses[i]*ses[i] + tau2)
			add(&swr, wr)
			add(&swbr, wr*betas[i])
		}
		res.Tau2 = tau2
		res.RandomBeta = swbr.s / swr.s
		res.RandomSE = math.Sqrt(1 / swr.s)
		res.RandomWeight = swr.s
	}
	return res
}
