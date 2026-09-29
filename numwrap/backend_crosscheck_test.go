package numwrap

import (
	"math"
	"testing"
)

// TestWLSBackendsAgree 验证同一输入下 stdlib 与 gonum 的 WLS
// 系数/SE/协方差在声明容差内一致。
func TestWLSBackendsAgree(t *testing.T) {
	x := []float64{0.5, 1.5, 2.1, 3.7, 4.2, 5.9, 6.3, 7.8, 8.1, 9.4}
	y := []float64{1.1, 1.9, 2.7, 4.4, 4.6, 6.6, 7.2, 8.5, 8.9, 10.4}
	weights := []float64{1, 1.2, 0.7, 1.5, 0.9, 1.1, 0.8, 1.3, 1.0, 0.6}
	cases := []RegressionOptions{
		{InterceptFree: true},
		{InterceptFree: false, FixedIntercept: 1},
		{InterceptFree: true, Weights: weights},
		{InterceptFree: false, FixedIntercept: 1, Weights: weights},
	}
	for ci, opts := range cases {
		rs, err := WLS(x, y, withBackend(opts, BackendStd))
		if err != nil {
			t.Fatalf("case %d stdlib: %v", ci, err)
		}
		rg, err := WLS(x, y, withBackend(opts, BackendGonum))
		if err != nil {
			t.Fatalf("case %d gonum: %v", ci, err)
		}
		if err := CompareWLSResults(rs, rg); err != nil {
			t.Fatalf("case %d: %v", ci, err)
		}
		t.Logf("[input] case=%d [judgement] stdlib=%v gonum=%v within rel=%g abs=%g",
			ci, rs.Coefficients, rg.Coefficients, CrosscheckRelTol, CrosscheckAbsTol)
	}
}

func withBackend(o RegressionOptions, b Backend) RegressionOptions {
	o.Backend = b
	return o
}

// TestWLSSameBackendBitExact 验证同一后端重复运行逐位一致。
func TestWLSSameBackendBitExact(t *testing.T) {
	x := []float64{0.5, 1.5, 2.1, 3.7, 4.2, 5.9, 6.3}
	y := []float64{1.1, 1.9, 2.7, 4.4, 4.6, 6.6, 7.2}
	for _, b := range []Backend{BackendStd, BackendGonum} {
		r1, err := WLS(x, y, RegressionOptions{InterceptFree: true, Backend: b})
		if err != nil {
			t.Fatal(err)
		}
		for k := 0; k < 5; k++ {
			r2, _ := WLS(x, y, RegressionOptions{InterceptFree: true, Backend: b})
			if !bitEqualSlices(r1.Coefficients, r2.Coefficients) ||
				!bitEqualSlices(r1.StdErr, r2.StdErr) {
				t.Fatalf("backend %s rerun not bit-exact: %v vs %v", b, r1, r2)
			}
		}
	}
}

func bitEqualSlices(a, b []float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// TestCrosscheckWLSAndPrimitives 验证 crosscheck 后端：
// 两边一致时返回结果；错误分类一致时返回同类错误。
func TestCrosscheckWLSAndPrimitives(t *testing.T) {
	x := []float64{1, 2, 3, 4, 5, 6}
	y := []float64{1.2, 2.1, 3.3, 3.9, 5.1, 5.8}
	r, err := WLS(x, y, RegressionOptions{InterceptFree: true, Backend: BackendCrosscheck})
	if err != nil {
		t.Fatal(err)
	}
	rs, _ := WLS(x, y, RegressionOptions{InterceptFree: true, Backend: BackendStd})
	if r.Coefficients[0] != rs.Coefficients[0] {
		t.Fatal("crosscheck must report the stdlib reference value")
	}

	// 非法输入在两个后端得到同一错误分类
	_, err = WLS([]float64{1, 2}, []float64{1, math.NaN()},
		RegressionOptions{Backend: BackendCrosscheck})
	if e, ok := AsError(err); !ok || e.Kind != KindNaN {
		t.Fatalf("want KindNaN from crosscheck, got %v", err)
	}
	_, err = WLS([]float64{1}, []float64{1, 2},
		RegressionOptions{Backend: BackendCrosscheck})
	if e, ok := AsError(err); !ok || e.Kind != KindInvalidArgument {
		t.Fatalf("want invalid_argument from crosscheck, got %v", err)
	}

	q, err := NormalCDFInvBackend(0.975, BackendCrosscheck)
	if err != nil || !withinTol(q, 1.959963984540054) {
		t.Fatalf("crosscheck quantile=%v err=%v", q, err)
	}
	m, s, err := MeanStdBackend(y, BackendCrosscheck)
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("[judgement] crosscheck WLS coef=%v q=%.10g mean/std=%.10g/%.10g",
		r.Coefficients, q, m, s)
}

// TestUnknownBackendRejected 验证未知后端显式拒绝、不回退默认实现。
func TestUnknownBackendRejected(t *testing.T) {
	for _, name := range []string{"", "mkl", "GONUM", " openblas"} {
		if _, err := ResolveBackend(name); err == nil {
			t.Fatalf("backend %q must be rejected", name)
		}
	}
	b, err := ResolveBackend("gonum")
	if err != nil || b != BackendGonum {
		t.Fatalf("gonum must resolve, got %v %v", b, err)
	}
	if _, err := WLS([]float64{1, 2, 3}, []float64{1, 2, 3},
		RegressionOptions{InterceptFree: true, Backend: Backend("cuda")}); err == nil {
		t.Fatal("unknown WLS backend must error")
	}
	got := AvailableBackends()
	want := []Backend{BackendStd, BackendGonum, BackendCrosscheck}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("available backends=%v want %v", got, want)
		}
	}
}

// TestConcurrentBackendSwitchRaceFree 在 -race 下并发切换进程级默认
// 后端，同时按调用显式指定后端求解：验证默认值切换与进行中的计算互不
// 干扰（计算只使用每调用 opts.Backend）。
func TestConcurrentBackendSwitchRaceFree(t *testing.T) {
	done := make(chan struct{})
	go func() {
		defer close(done)
		for i := 0; i < 200; i++ {
			_ = SetDefaultBackend(BackendGonum)
			_ = SetDefaultBackend(BackendStd)
		}
	}()
	x := []float64{1, 2, 3, 4, 5, 6, 7, 8}
	y := []float64{2.1, 2.9, 4.2, 4.8, 6.1, 7.0, 7.9, 9.2}
	for i := 0; i < 200; i++ {
		b := BackendStd
		if i%2 == 0 {
			b = BackendGonum
		}
		r, err := WLS(x, y, RegressionOptions{InterceptFree: true, Backend: b})
		if err != nil {
			t.Fatal(err)
		}
		// 无论默认值如何切换，显式后端的结果必须始终一致。
		r2, _ := WLS(x, y, RegressionOptions{InterceptFree: true, Backend: b})
		if r.Coefficients[1] != r2.Coefficients[1] {
			t.Fatalf("per-call backend result drifted during switch")
		}
	}
	<-done
}
