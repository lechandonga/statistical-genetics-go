package numwrap

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

// wlsGonum 使用 gonum mat 的稠密 QR 求解 WLS，并把 gonum 的任何
// panic/错误翻译成本包 *Error，绝不向接口外透出第三方类型/字符串。
//
// 数值口径与 wlsStd 相同：
//   - 自由截距：设计矩阵 [1, x]，系数 [intercept, slope]；
//   - 固定截距 c：y' = y-c，过原点回归，系数 [slope]；
//   - 残差方差用无偏估计 rss/(n-p)（p=2 或 1）；
//   - 协方差 = sigma2 * (X'WX)^-1。
func wlsGonum(x, y []float64, opts RegressionOptions) (res *RegressionResult, errOut error) {
	defer func() {
		if r := recover(); r != nil {
			res = nil
			errOut = translateGonumPanic("WLS", r)
		}
	}()

	n := len(x)
	hasW := len(opts.Weights) > 0
	wi := func(i int) float64 {
		if hasW {
			return opts.Weights[i]
		}
		return 1
	}

	p := 2
	if !opts.InterceptFree {
		p = 1
	}

	// 对 sqrt(w) 加权：Xw = sqrt(W) X，yw = sqrt(W) y。
	xw := make([]float64, n*p)
	yw := make([]float64, n)
	var sw float64
	for i := 0; i < n; i++ {
		w := wi(i)
		sw += w
		s := math.Sqrt(w)
		if opts.InterceptFree {
			xw[i*2] = s
			xw[i*2+1] = s * x[i]
		} else {
			xw[i] = s * x[i]
		}
		yc := y[i]
		if !opts.InterceptFree {
			yc = y[i] - opts.FixedIntercept
		}
		yw[i] = s * yc
	}
	if sw <= 0 {
		return nil, &Error{Kind: KindInvalidArgument, Op: "WLS", Detail: "total weight is zero"}
	}

	X := mat.NewDense(n, p, xw)
	b := mat.NewVecDense(n, yw)

	var qr mat.QR
	qr.Factorize(X)
	var beta mat.VecDense
	if err := qr.SolveVecTo(&beta, false, b); err != nil {
		return nil, &Error{Kind: KindSingular, Op: "WLS",
			Detail: "design matrix rank-deficient"}
	}

	// 残差与 RSS。
	var resid mat.VecDense
	resid.MulVec(X, &beta)
	resid.SubVec(b, &resid)
	rss := mat.Dot(&resid, &resid)
	if rss < 0 {
		rss = 0
	}
	df := float64(n - p)
	if df <= 0 {
		return nil, &Error{Kind: KindInvalidArgument, Op: "WLS",
			Detail: "no residual degrees of freedom"}
	}
	residVar := rss / df

	// 协方差矩阵 sigma2 * (Xw'Xw)^-1。
	var XtX mat.Dense
	XtX.Mul(X.T(), X)
	var covM mat.Dense
	if err := covM.Inverse(&XtX); err != nil {
		// 第三方错误不透出：仅用于判定奇异。
		return nil, &Error{Kind: KindSingular, Op: "WLS",
			Detail: "normal-equation matrix not invertible"}
	}
	covM.Scale(residVar, &covM)

	for _, v := range beta.RawVector().Data {
		if math.IsNaN(v) {
			return nil, &Error{Kind: KindNaN, Op: "WLS", Detail: "NaN coefficient"}
		}
		if math.IsInf(v, 0) {
			return nil, &Error{Kind: KindInf, Op: "WLS", Detail: "Inf coefficient"}
		}
	}

	res = &RegressionResult{N: n, ResidualSE: math.Sqrt(max0(residVar))}
	if opts.InterceptFree {
		res.Coefficients = []float64{beta.AtVec(0), beta.AtVec(1)}
		res.StdErr = []float64{math.Sqrt(max0(covM.At(0, 0))), math.Sqrt(max0(covM.At(1, 1)))}
		res.Cov = []float64{covM.At(0, 0), covM.At(0, 1), covM.At(1, 0), covM.At(1, 1)}
	} else {
		res.Coefficients = []float64{beta.AtVec(0)}
		res.StdErr = []float64{math.Sqrt(max0(covM.At(0, 0)))}
		res.Cov = []float64{covM.At(0, 0)}
	}
	return res, nil
}

// translateGonumPanic 把 gonum 内部 panic 归一化：数值型 panic 映射为
// nan/inf，其余统一为 unknown（不含第三方库名与内部文本）。
func translateGonumPanic(op string, r any) *Error {
	switch v := r.(type) {
	case *Error:
		return v
	case float64:
		if math.IsNaN(v) {
			return &Error{Kind: KindNaN, Op: op, Detail: "NaN produced by backend"}
		}
		if math.IsInf(v, 0) {
			return &Error{Kind: KindInf, Op: op, Detail: "Inf produced by backend"}
		}
	}
	return &Error{Kind: KindUnknown, Op: op,
		Detail: fmt.Sprintf("third-party backend internal failure: %T", r)}
}

// withinTol 按 numwrap 声明的交叉验证容差比较两个标量。
func withinTol(a, b float64) bool {
	if math.IsNaN(a) || math.IsNaN(b) {
		// 两边同为 NaN 视为一致（例如固定截距模式没有 intercept 概念）
		return math.IsNaN(a) && math.IsNaN(b)
	}
	if math.IsInf(a, 0) || math.IsInf(b, 0) {
		return a == b
	}
	d := math.Abs(a - b)
	scale := math.Max(math.Abs(a), math.Abs(b))
	return d <= math.Max(CrosscheckAbsTol, CrosscheckRelTol*scale)
}

// WithinTol 导出 withinTol 供上层（stats 层交叉验证）复用同一容差口径。
func WithinTol(a, b float64) bool { return withinTol(a, b) }

// compareErrors 校验两个后端对同一输入的错误分类一致。
// 分类一致时返回该分类错误（不透出任何一边的内部细节差异）；
// 不一致时返回 KindBackendMismatch（一个成功一个失败也属于不一致）。
func compareErrors(op string, err1, err2 error) (*RegressionResult, error) {
	e1, ok1 := AsError(err1)
	e2, ok2 := AsError(err2)
	if ok1 && ok2 && e1.Kind == e2.Kind {
		return nil, e1
	}
	return nil, &Error{Kind: KindBackendMismatch, Op: op,
		Detail: fmt.Sprintf("error classification differs: stdlib=%v gonum=%v", err1, err2)}
}

// CompareWLSResults 比较两个后端的 WLS 结果；任何系数/方差分量超出
// 声明容差都返回 KindBackendMismatch。
func CompareWLSResults(a, b *RegressionResult) error {
	if a == nil || b == nil {
		return &Error{Kind: KindBackendMismatch, Op: "WLS", Detail: "nil result"}
	}
	if a.N != b.N {
		return &Error{Kind: KindBackendMismatch, Op: "WLS",
			Detail: fmt.Sprintf("N differs: %d vs %d", a.N, b.N)}
	}
	if len(a.Coefficients) != len(b.Coefficients) ||
		len(a.StdErr) != len(b.StdErr) || len(a.Cov) != len(b.Cov) {
		return &Error{Kind: KindBackendMismatch, Op: "WLS", Detail: "result shape differs"}
	}
	check := func(name string, u, v []float64) error {
		for i := range u {
			if !withinTol(u[i], v[i]) {
				return &Error{Kind: KindBackendMismatch, Op: "WLS",
					Detail: fmt.Sprintf("%s[%d] differs: %.17g vs %.17g", name, i, u[i], v[i])}
			}
		}
		return nil
	}
	if err := check("coefficients", a.Coefficients, b.Coefficients); err != nil {
		return err
	}
	if err := check("stderr", a.StdErr, b.StdErr); err != nil {
		return err
	}
	if err := check("cov", a.Cov, b.Cov); err != nil {
		return err
	}
	if !withinTol(a.ResidualSE, b.ResidualSE) {
		return &Error{Kind: KindBackendMismatch, Op: "WLS",
			Detail: fmt.Sprintf("residualSE differs: %.17g vs %.17g", a.ResidualSE, b.ResidualSE)}
	}
	return nil
}
