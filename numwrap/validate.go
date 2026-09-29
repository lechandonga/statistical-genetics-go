package numwrap

import "math"

// validateWLSInputs 是两个 WLS 实现共享的输入校验：
// 维度、观测数、权重非负、有限性。错误一律为本包 *Error。
func validateWLSInputs(x, y []float64, opts RegressionOptions) error {
	n := len(x)
	hasW := len(opts.Weights) > 0
	if hasW && len(opts.Weights) != n {
		return &Error{Kind: KindInvalidArgument, Op: "WLS",
			Detail: "weights length mismatch"}
	}
	need := 2
	if !opts.InterceptFree {
		need = 1
	}
	if n < need+1 {
		return &Error{Kind: KindInvalidArgument, Op: "WLS",
			Detail: "not enough observations for requested model"}
	}
	for i := 0; i < n; i++ {
		wi := 1.0
		if hasW {
			wi = opts.Weights[i]
		}
		if math.IsNaN(x[i]) || math.IsNaN(y[i]) || math.IsNaN(wi) {
			return &Error{Kind: KindNaN, Op: "WLS", Detail: "NaN in inputs"}
		}
		if math.IsInf(x[i], 0) || math.IsInf(y[i], 0) || math.IsInf(wi, 0) {
			return &Error{Kind: KindInf, Op: "WLS", Detail: "Inf in inputs"}
		}
		if wi < 0 {
			return &Error{Kind: KindInvalidArgument, Op: "WLS",
				Detail: "negative weight"}
		}
	}
	return nil
}
