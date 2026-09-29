package numwrap

import (
	"math"
	"testing"
)

func approx(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

func TestWLSDeterministicLine(t *testing.T) {
	x := []float64{1, 2, 3, 4, 5}
	y := []float64{2 + 3*1, 2 + 3*2, 2 + 3*3, 2 + 3*4, 2 + 3*5}
	r1, err := WLS(x, y, RegressionOptions{InterceptFree: true})
	if err != nil {
		t.Fatal(err)
	}
	// 打乱顺序的子集不行；相同输入重复调用必须逐位相等
	r2, _ := WLS(x, y, RegressionOptions{InterceptFree: true})
	if r1.Coefficients[0] != r2.Coefficients[0] || r1.Coefficients[1] != r2.Coefficients[1] {
		t.Fatal("repeated WLS results differ")
	}
	if !approx(r1.Coefficients[0], 2, 1e-12) || !approx(r1.Coefficients[1], 3, 1e-12) {
		t.Fatalf("coef=%v", r1.Coefficients)
	}
	if r1.ResidualSE != 0 {
		t.Fatalf("perfect fit residual SE=%v", r1.ResidualSE)
	}
	t.Logf("[input] y=2+3x n=5 [judgement] intercept=%.10f slope=%.10f residSE=%.2e",
		r1.Coefficients[0], r1.Coefficients[1], r1.ResidualSE)
}

func TestWLSFixedIntercept(t *testing.T) {
	// y = 1 + 2x + 噪声；固定截距 1
	x := []float64{1, 2, 3, 4, 5, 6}
	y := []float64{3.1, 4.9, 7.0, 9.2, 10.8, 13.1}
	r, err := WLS(x, y, RegressionOptions{InterceptFree: false, FixedIntercept: 1})
	if err != nil {
		t.Fatal(err)
	}
	if len(r.Coefficients) != 1 || !approx(r.Coefficients[0], 2, 0.05) {
		t.Fatalf("slope=%v", r.Coefficients)
	}
	t.Logf("[input] y~1+2x, fixed intercept=1 [judgement] slope=%.4f se=%.4f",
		r.Coefficients[0], r.StdErr[0])
}

func TestWLSErrors(t *testing.T) {
	if _, err := WLS([]float64{1}, []float64{1, 2}, RegressionOptions{}); err == nil {
		t.Fatal("length mismatch must error")
	}
	if _, err := WLS([]float64{1}, []float64{1}, RegressionOptions{}); err == nil {
		t.Fatal("too few obs must error")
	}
	_, err := WLS([]float64{1, 2}, []float64{1, math.NaN()}, RegressionOptions{})
	if e, ok := AsError(err); !ok || e.Kind != KindNaN {
		t.Fatalf("want KindNaN got %v", err)
	}
	_, err = WLS([]float64{1, 2}, []float64{1, 2},
		RegressionOptions{InterceptFree: false, FixedIntercept: 0, Weights: []float64{-1, 1}})
	if e, ok := AsError(err); !ok || e.Kind != KindInvalidArgument {
		t.Fatalf("want invalid_argument got %v", err)
	}
	// 未编译/未知后端必须显式拒绝，绝不静默回退
	_, err = WLS([]float64{1, 2}, []float64{1, 2},
		RegressionOptions{InterceptFree: false, Backend: Backend("blis")})
	if e, ok := AsError(err); !ok || e.Kind != KindUnsupportedOption {
		t.Fatalf("want unsupported_option got %v", err)
	}
	t.Logf("[input] mismatch/NaN/negative-weight/unknown-backend [judgement] all mapped to typed numwrap errors")
}

func TestKahanAndQuantile(t *testing.T) {
	xs := []float64{1e16, -1e16, 1}
	sum, err := KahanSum(xs)
	if err != nil || sum != 1 {
		t.Fatalf("kahan sum=%v err=%v", sum, err)
	}
	// 重复调用逐位一致
	sum2, _ := KahanSum(xs)
	if sum != sum2 {
		t.Fatal("KahanSum not reproducible")
	}
	m, sd, err := MeanStd([]float64{1, 2, 3})
	if err != nil || !approx(m, 2, 1e-12) || !approx(sd, 1, 1e-12) {
		t.Fatalf("mean=%v sd=%v", m, sd)
	}
	q, err := NormalCDFInv(0.975)
	if err != nil || !approx(q, 1.9599639845, 1e-6) {
		t.Fatalf("q97.5=%v", q)
	}
	if _, err := NormalCDFInv(0); err == nil {
		t.Fatal("p=0 must error")
	}
	t.Logf("[input] {1e16,1,-1e16} [judgement] KahanSum=%v exactly, reproducible; q(0.975)=%.6f", sum, q)
}
