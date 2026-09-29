package numwrap

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/gonum/stat/distuv"
)

// meanStdGonum 用 gonum/stat 计算均值与样本标准差（N-1 分母），
// 与 stdlib 口径一致；任何内部 panic 翻译为本包错误。
func meanStdGonum(xs []float64) (mean, std float64, errOut error) {
	defer func() {
		if r := recover(); r != nil {
			errOut = translateGonumPanic("MeanStd", r)
			mean, std = 0, 0
		}
	}()
	if len(xs) == 0 {
		return 0, 0, &Error{Kind: KindInvalidArgument, Op: "MeanStd", Detail: "empty"}
	}
	for _, v := range xs {
		if math.IsNaN(v) {
			return 0, 0, &Error{Kind: KindNaN, Op: "MeanStd"}
		}
		if math.IsInf(v, 0) {
			return 0, 0, &Error{Kind: KindInf, Op: "MeanStd"}
		}
	}
	mean = stat.Mean(xs, nil)
	if len(xs) == 1 {
		return mean, 0, nil
	}
	variance := stat.Variance(xs, nil)
	if variance < 0 {
		variance = 0
	}
	return mean, math.Sqrt(variance), nil
}

// quantileGonum 用 gonum 标准正态分布的分位数函数。
func quantileGonum(p float64) (q float64, errOut error) {
	defer func() {
		if r := recover(); r != nil {
			errOut = translateGonumPanic("NormalCDFInv", r)
			q = 0
		}
	}()
	if math.IsNaN(p) {
		return 0, &Error{Kind: KindNaN, Op: "NormalCDFInv"}
	}
	if p <= 0 || p >= 1 {
		return 0, &Error{Kind: KindInvalidArgument, Op: "NormalCDFInv",
			Detail: "p must be in (0,1)"}
	}
	return distuv.UnitNormal.Quantile(p), nil
}

func compareSimpleErrors(op string, err1, err2 error) (error, bool) {
	e1, ok1 := AsError(err1)
	e2, ok2 := AsError(err2)
	if ok1 && ok2 && e1.Kind == e2.Kind {
		return e1, true
	}
	return nil, false
}

func fmtSprintMismatch(m1, s1, m2, s2 float64) string {
	return fmt.Sprintf("mean/std differs: (%.17g,%.17g) vs (%.17g,%.17g)", m1, s1, m2, s2)
}
