package numwrap

// gonum 后端：基于 gonum/gonum 的 WLS 独立第二实现。
//
// 本文件只暴露 numwrap 自有类型；gonum 的错误/panic 一律经
// guard() 翻译为 *Error（见 errors 映射表与 docs/compatibility.md），
// 不向接口外透出任何第三方错误文本。

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/gonum/stat/distuv"
)

// gonumNormalQuantile 用 gonum 分布包计算标准正态分位数（独立第二实现）。
func gonumNormalQuantile(p float64) (q float64, err error) {
	defer func() {
		if r := recover(); r != nil {
			err = &Error{Kind: KindInvalidArgument, Op: "NormalCDFInv[gonum]",
				Detail: "quantile failed"}
		}
	}()
	return distuv.UnitNormal.Quantile(p), nil
}

// gonumMeanStd 用 gonum/stat 独立计算均值与（N-1 分母）样本标准差。
func gonumMeanStd(xs []float64) (mean, sd float64, err error) {
	defer func() {
		if r := recover(); r != nil {
			err = &Error{Kind: KindNaN, Op: "MeanStd[gonum]", Detail: "moment failed"}
		}
	}()
	mean = stat.Mean(xs, nil)
	v := stat.Variance(xs, nil)
	if math.IsNaN(mean) || math.IsNaN(v) {
		return 0, 0, &Error{Kind: KindNaN, Op: "MeanStd[gonum]", Detail: "NaN moment"}
	}
	if v < 0 {
		v = 0
	}
	return mean, math.Sqrt(v), nil
}

// wlsGonum 用 gonum/mat 的 QR 求解加权最小二乘。
//
// 加权方式：对设计矩阵 X 与响应 y 逐行乘 sqrt(w)，再求解
// Xw' * b = yw 的最小二乘解；残差方差与系数协方差按经典
// WLS 公式 b' cov = sigma2 * (Xw'Xw)^-1 计算。
func wlsGonum(x, y []float64, opts RegressionOptions) (result *RegressionResult, err error) {
	defer func() {
		if r := recover(); r != nil {
			// gonum 在维度/奇异问题上可能 panic：绝不原样漏出。
			err = &Error{Kind: KindSingular, Op: "WLS[gonum]",
				Detail: fmt.Sprintf("linear algebra failure (%T)", r)}
		}
	}()

	n := len(x)
	hasW := len(opts.Weights) > 0
	p := 1
	if opts.InterceptFree {
		p = 2
	}

	// 固定截距：平移响应后做过原点回归。
	c := 0.0
	if !opts.InterceptFree {
		c = opts.FixedIntercept
	}

	xw := mat.NewDense(n, p, nil)
	yw := mat.NewVecDense(n, nil)
	var sw float64
	for i := 0; i < n; i++ {
		wi := 1.0
		if hasW {
			wi = opts.Weights[i]
		}
		s := math.Sqrt(wi)
		sw += wi
		if opts.InterceptFree {
			xw.Set(i, 0, s)
			xw.Set(i, 1, s*x[i])
		} else {
			xw.Set(i, 0, s*x[i])
		}
		yw.SetVec(i, s*(y[i]-c))
	}

	// QR 最小二乘：SolveVec 对 Xw (n×p, n>=p) 求最小范数残差解。
	var qr mat.QR
	qr.Factorize(xw)
	var coef mat.VecDense
	if err := qr.SolveVecTo(&coef, false, yw); err != nil {
		// 设计矩阵秩亏（x 零方差等）。
		return nil, &Error{Kind: KindSingular, Op: "WLS[gonum]",
			Detail: "rank-deficient design matrix"}
	}

	// 残差与 RSS（在加权空间中）。
	var pred mat.VecDense
	pred.MulVec(xw, &coef)
	rss := 0.0
	for i := 0; i < n; i++ {
		r := yw.AtVec(i) - pred.AtVec(i)
		rss += r * r
	}
	if rss < 0 {
		rss = 0
	}

	// 协方差：sigma2 * (Xw'Xw)^-1。
	var xtx mat.Dense
	xtx.Mul(xw.T(), xw)
	var inv mat.Dense
	if err := inv.Inverse(&xtx); err != nil {
		return nil, &Error{Kind: KindSingular, Op: "WLS[gonum]",
			Detail: "weighted information matrix singular"}
	}

	res := &RegressionResult{N: n}
	df := n - p
	if df <= 0 {
		return nil, &Error{Kind: KindInvalidArgument, Op: "WLS[gonum]",
			Detail: "not enough observations for residual degrees of freedom"}
	}
	residVar := rss / float64(df)
	if residVar < 0 {
		residVar = 0
	}
	if opts.InterceptFree {
		intercept := coef.AtVec(0)
		slope := coef.AtVec(1)
		varI := residVar * inv.At(0, 0)
		varS := residVar * inv.At(1, 1)
		cov := residVar * inv.At(0, 1)
		res.Coefficients = []float64{intercept, slope}
		res.StdErr = []float64{math.Sqrt(max0(varI)), math.Sqrt(max0(varS))}
		res.Cov = []float64{varI, cov, cov, varS}
	} else {
		slope := coef.AtVec(0)
		varS := residVar * inv.At(0, 0)
		res.Coefficients = []float64{slope}
		res.StdErr = []float64{math.Sqrt(max0(varS))}
		res.Cov = []float64{varS}
	}
	_ = sw
	res.ResidualSE = math.Sqrt(max0(residVar))

	for _, v := range res.Coefficients {
		if math.IsNaN(v) {
			return nil, &Error{Kind: KindNaN, Op: "WLS[gonum]", Detail: "NaN coefficient"}
		}
		if math.IsInf(v, 0) {
			return nil, &Error{Kind: KindInf, Op: "WLS[gonum]", Detail: "Inf coefficient"}
		}
	}
	return res, nil
}
