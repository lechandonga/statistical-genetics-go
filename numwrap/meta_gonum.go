package numwrap

// gonum 后端的元分析组合：与 metaCombineImpl 独立的实现路径。
//
// 这里刻意不调用 metaCombineImpl，也不使用 Kahan 补偿，
// 而是用最简单直观的顺序累加（与固定公式定义一一对应），
// 使两套代码在任何共享子例程出错时仍能互相暴露问题。

import "math"

func metaCombineGonum(betas, ses []float64, random bool) (res *MetaResult) {
	n := len(betas)
	weights := make([]float64, n)
	var sumW, sumWB float64
	for i := 0; i < n; i++ {
		weights[i] = 1.0 / (ses[i] * ses[i])
		sumW += weights[i]
		sumWB += weights[i] * betas[i]
	}
	beta := sumWB / sumW
	var q float64
	for i := 0; i < n; i++ {
		d := betas[i] - beta
		q += weights[i] * d * d
	}
	df := n - 1
	i2 := 0.0
	if q > float64(df) {
		i2 = (q - float64(df)) / q
	}
	out := &MetaResult{
		FixedBeta:   beta,
		FixedSE:     math.Sqrt(1.0 / sumW),
		FixedWeight: sumW,
		Q:           q,
		I2:          i2,
		Df:          df,
		Weights:     weights,
	}
	if random {
		var sumW2 float64
		for i := 0; i < n; i++ {
			sumW2 += weights[i] * weights[i]
		}
		tau2 := 0.0
		if q > float64(df) {
			if den := sumW - sumW2/sumW; den > 0 {
				tau2 = (q - float64(df)) / den
			}
		}
		var swr, swbr float64
		for i := 0; i < n; i++ {
			wr := 1.0 / (ses[i]*ses[i] + tau2)
			swr += wr
			swbr += wr * betas[i]
		}
		out.Tau2 = tau2
		out.RandomBeta = swbr / swr
		out.RandomSE = math.Sqrt(1.0 / swr)
		out.RandomWeight = swr
	}
	return out
}
