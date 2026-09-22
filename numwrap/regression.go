package numwrap

import (
	"math"
)

// RegressionOptions 控制加权线性回归。
type RegressionOptions struct {
	Backend Backend // 空 => 默认后端
	// InterceptFree=true 时自由估计截距；false 时截距固定为 FixedIntercept
	// （LDSC 约束截距模式：固定为 1）。
	InterceptFree  bool
	FixedIntercept float64
	// Weights 为空表示等权；长度必须与观测数一致，权重须有限且 >=0。
	Weights []float64
}

// RegressionResult 是加权最小二乘结果。
type RegressionResult struct {
	Coefficients []float64 // 自由截距: [intercept, slope]；固定截距: [slope]
	StdErr       []float64
	Cov          []float64 // 行优先 2x2；固定截距时长度 1（slope 方差）
	ResidualSE   float64
	N            int
}

// WLS 求解 y = intercept + slope*x（单自变量）的加权最小二乘。
//
// 结果只依赖输入数值与权重顺序：内部以固定顺序 Kahan 累加，
// 相同输入在任何机器/并行度下逐位一致。
func WLS(x, y []float64, opts RegressionOptions) (*RegressionResult, error) {
	if opts.Backend == "" {
		opts.Backend = DefaultBackend()
	}
	if opts.Backend != BackendStd {
		return nil, &Error{Kind: KindUnsupportedOption, Op: "WLS",
			Detail: "backend not available: " + string(opts.Backend)}
	}
	n := len(x)
	if n != len(y) {
		return nil, &Error{Kind: KindInvalidArgument, Op: "WLS",
			Detail: "x/y length mismatch"}
	}
	hasW := len(opts.Weights) > 0
	if hasW && len(opts.Weights) != n {
		return nil, &Error{Kind: KindInvalidArgument, Op: "WLS",
			Detail: "weights length mismatch"}
	}
	need := 2
	if !opts.InterceptFree {
		need = 1
	}
	if n < need+1 { // 自由截距至少需要 2 点（这里要求 n>=3 求残差）；固定至少 2
		return nil, &Error{Kind: KindInvalidArgument, Op: "WLS",
			Detail: "not enough observations for requested model"}
	}

	w := func(i int) float64 {
		if hasW {
			return opts.Weights[i]
		}
		return 1
	}
	for i := 0; i < n; i++ {
		wi := w(i)
		if math.IsNaN(x[i]) || math.IsNaN(y[i]) || math.IsNaN(wi) {
			return nil, &Error{Kind: KindNaN, Op: "WLS", Detail: "NaN in inputs"}
		}
		if math.IsInf(x[i], 0) || math.IsInf(y[i], 0) || math.IsInf(wi, 0) {
			return nil, &Error{Kind: KindInf, Op: "WLS", Detail: "Inf in inputs"}
		}
		if wi < 0 {
			return nil, &Error{Kind: KindInvalidArgument, Op: "WLS",
				Detail: "negative weight"}
		}
	}

	// 固定截距：对 y' = y - c 做过原点回归。
	c := 0.0
	if !opts.InterceptFree {
		c = opts.FixedIntercept
	}

	// 累加顺序固定为 i=0..n-1 的 Kahan 求和，保证逐位可复现。
	var sw, sxx, sxy, syy float64
	var csw, csxx, csxy, csyy float64
	add := func(acc *float64, comp *float64, v float64) {
		yv := v - *comp
		t := *acc + yv
		*comp = (t - *acc) - yv
		*acc = t
	}
	for i := 0; i < n; i++ {
		wi := w(i)
		yc := y[i] - c
		add(&sw, &csw, wi)
		add(&sxx, &csxx, wi*x[i]*x[i])
		add(&sxy, &csxy, wi*x[i]*yc)
		add(&syy, &csyy, wi*yc*yc)
	}
	if sw <= 0 {
		return nil, &Error{Kind: KindInvalidArgument, Op: "WLS", Detail: "total weight is zero"}
	}

	res := &RegressionResult{N: n}
	if opts.InterceptFree {
		// 需要额外的 Sx, Sy
		var sx, sy float64
		var csx, csy float64
		for i := 0; i < n; i++ {
			wi := w(i)
			add(&sx, &csx, wi*x[i])
			add(&sy, &csy, wi*(y[i]))
		}
		xbar := sx / sw
		ybar := sy / sw
		// Sxx_c = sxx - sx^2/sw ; Sxy_c = sxy'... 使用中心化公式：
		sxxC := sxx - sx*sx/sw
		// syy/sxy 上面是围绕固定截距 c 的，需要重新围绕均值计算
		var sxyC, syyC float64
		var k1, k2 float64
		for i := 0; i < n; i++ {
			wi := w(i)
			dx := x[i] - xbar
			dy := y[i] - ybar
			add(&sxyC, &k1, wi*dx*dy)
			add(&syyC, &k2, wi*dy*dy)
		}
		if sxxC <= 0 {
			return nil, &Error{Kind: KindSingular, Op: "WLS", Detail: "zero variance in x"}
		}
		slope := sxyC / sxxC
		intercept := ybar - slope*xbar
		rss := syyC - slope*sxyC
		df := float64(n - 2)
		residVar := rss / df
		if residVar < 0 {
			rss = 0
			residVar = 0
		}
		// 经典 WLS 协方差（条件于权重）：
		// Var(slope) = sigma2 / SxxC ; Var(intercept) = sigma2*(1/sw + xbar^2/SxxC)
		varSlope := residVar / sxxC
		varIntercept := residVar * (1/sw + xbar*xbar/sxxC)
		cov := residVar * (-xbar / sxxC)
		res.Coefficients = []float64{intercept, slope}
		res.StdErr = []float64{math.Sqrt(max0(varIntercept)), math.Sqrt(max0(varSlope))}
		res.Cov = []float64{varIntercept, cov, cov, varSlope}
		res.ResidualSE = math.Sqrt(max0(residVar))
		return res, nil
	}

	// 固定截距：slope = sxy / sxx
	if sxx <= 0 {
		return nil, &Error{Kind: KindSingular, Op: "WLS", Detail: "zero weighted sum of x^2"}
	}
	slope := sxy / sxx
	rss := syy - slope*sxy
	if rss < 0 {
		rss = 0
	}
	residVar := rss / float64(n-1)
	varSlope := residVar / sxx
	res.Coefficients = []float64{slope}
	res.StdErr = []float64{math.Sqrt(max0(varSlope))}
	res.Cov = []float64{varSlope}
	res.ResidualSE = math.Sqrt(max0(residVar))
	return res, nil
}

func max0(v float64) float64 {
	if v < 0 {
		return 0
	}
	return v
}

// KahanSum 以输入切片顺序执行补偿求和；遇 NaN 返回 KindNaN，遇 Inf 返回 KindInf。
func KahanSum(xs []float64) (float64, error) {
	var sum, comp float64
	for _, v := range xs {
		if math.IsNaN(v) {
			return 0, &Error{Kind: KindNaN, Op: "KahanSum"}
		}
		if math.IsInf(v, 0) {
			return 0, &Error{Kind: KindInf, Op: "KahanSum"}
		}
		y := v - comp
		t := sum + y
		comp = (t - sum) - y
		sum = t
	}
	return sum, nil
}

// MeanStd 返回均值与（N-1 分母）样本标准差，单遍确定性算法。
func MeanStd(xs []float64) (mean, std float64, err error) {
	n := len(xs)
	if n == 0 {
		return 0, 0, &Error{Kind: KindInvalidArgument, Op: "MeanStd", Detail: "empty"}
	}
	s, e := KahanSum(xs)
	if e != nil {
		return 0, 0, e
	}
	mean = s / float64(n)
	if n == 1 {
		return mean, 0, nil
	}
	var ss, comp float64
	for _, v := range xs {
		d := v - mean
		add := d*d - comp
		t := ss + add
		comp = (t - ss) - add
		ss = t
	}
	variance := ss / float64(n-1)
	if variance < 0 {
		variance = 0
	}
	return mean, math.Sqrt(variance), nil
}

// NormalCDFInv 标准正态分位数函数，基于 math.Erfcinv（stdlib，逐平台一致）。
func NormalCDFInv(p float64) (float64, error) {
	if math.IsNaN(p) {
		return 0, &Error{Kind: KindNaN, Op: "NormalCDFInv"}
	}
	if p <= 0 || p >= 1 {
		return 0, &Error{Kind: KindInvalidArgument, Op: "NormalCDFInv",
			Detail: "p must be in (0,1)"}
	}
	return -math.Sqrt2 * math.Erfcinv(2*p), nil
}
