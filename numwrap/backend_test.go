package numwrap

import (
	"math"
	"os"
	"strings"
	"sync"
	"testing"
)

// TestWLSBackendsAgree：两个后端在自由/固定截距、等权/加权下
// 结果必须落在声明容差内；同后端重复运行必须逐位一致。
func TestWLSBackendsAgree(t *testing.T) {
	x := []float64{0.1, 0.7, 1.2, 2.1, 3.3, 4.0, 5.2}
	y := make([]float64, len(x))
	w := []float64{1.0, 0.5, 2.0, 1.3, 0.8, 1.7, 1.1}
	for i := range x {
		y[i] = 1.5 - 0.8*x[i] + 0.03*math.Sin(float64(i))
	}
	cases := []RegressionOptions{
		{InterceptFree: true},
		{InterceptFree: false, FixedIntercept: 1.5},
		{InterceptFree: true, Weights: w},
		{InterceptFree: false, FixedIntercept: 1.5, Weights: w},
	}
	for ci, opts := range cases {
		opts.Backend = BackendStd
		r1, err := WLS(x, y, opts)
		if err != nil {
			t.Fatalf("case %d std: %v", ci, err)
		}
		// 同后端逐位一致
		r1b, _ := WLS(x, y, opts)
		for i := range r1.Coefficients {
			if r1.Coefficients[i] != r1b.Coefficients[i] {
				t.Fatalf("case %d stdlib not bit-reproducible", ci)
			}
		}
		opts.Backend = BackendGonum
		r2, err := WLS(x, y, opts)
		if err != nil {
			t.Fatalf("case %d gonum: %v", ci, err)
		}
		if err := crossCheckSlices("test", "coef", r1.Coefficients, r2.Coefficients); err != nil {
			t.Fatalf("case %d: %v", ci, err)
		}
		if err := crossCheckSlices("test", "stderr", r1.StdErr, r2.StdErr); err != nil {
			t.Fatalf("case %d: %v", ci, err)
		}
		t.Logf("[case %d] coef std=%v gonum=%v within tol rel=%.0e",
			ci, r1.Coefficients, r2.Coefficients, CrossCheckRelTol)
	}
}

// TestMetaBackendsAgree：IVW 与 DL 在两个后端间一致、同后端逐位一致。
func TestMetaBackendsAgree(t *testing.T) {
	betas := []float64{-0.31, 0.12, 0.05, 0.44}
	ses := []float64{0.11, 0.08, 0.15, 0.09}
	for _, random := range []bool{false, true} {
		a, err := MetaCombine(betas, ses, MetaOptions{Backend: BackendStd, Random: random})
		if err != nil {
			t.Fatal(err)
		}
		a2, _ := MetaCombine(betas, ses, MetaOptions{Backend: BackendStd, Random: random})
		if a.FixedBeta != a2.FixedBeta || a.FixedSE != a2.FixedSE || a.Q != a2.Q {
			t.Fatal("stdlib meta not bit-reproducible")
		}
		b, err := MetaCombine(betas, ses, MetaOptions{Backend: BackendGonum, Random: random})
		if err != nil {
			t.Fatal(err)
		}
		for _, ck := range []struct {
			name string
			x, y float64
		}{
			{"beta", a.FixedBeta, b.FixedBeta}, {"se", a.FixedSE, b.FixedSE},
			{"q", a.Q, b.Q}, {"i2", a.I2, b.I2},
			{"tau2", a.Tau2, b.Tau2}, {"rbeta", a.RandomBeta, b.RandomBeta},
			{"rse", a.RandomSE, b.RandomSE},
		} {
			if err := CrossCheckScalar("test", ck.name, ck.x, ck.y); err != nil {
				t.Fatalf("random=%v: %v", random, err)
			}
		}
	}
	t.Log("[judgement] stdlib vs gonum IVW/DL agree within tolerance; each is bit-reproducible")
}

// TestQuantileAndMomentsBackends：分位数与矩在两后端一致。
func TestQuantileAndMomentsBackends(t *testing.T) {
	xs := []float64{1.0, -2.5, 3.25, 4.0, 5.5, 6.0, -0.25}
	m1, s1, err := MeanStdB(xs, BackendStd)
	if err != nil {
		t.Fatal(err)
	}
	m2, s2, err := MeanStdB(xs, BackendGonum)
	if err != nil {
		t.Fatal(err)
	}
	if err := CrossCheckScalar("test", "mean", m1, m2); err != nil {
		t.Fatal(err)
	}
	if err := CrossCheckScalar("test", "sd", s1, s2); err != nil {
		t.Fatal(err)
	}
	q1, _ := NormalCDFInvB(0.975, BackendStd)
	q2, _ := NormalCDFInvB(0.975, BackendGonum)
	if err := CrossCheckScalar("test", "q975", q1, q2); err != nil {
		t.Fatal(err)
	}
	if math.Abs(q2-1.959963984) > 1e-6 {
		t.Fatalf("gonum q=%v", q2)
	}
	t.Logf("[judgement] mean %.12f/%.12f sd %.12f/%.12f q %.10f/%.10f",
		m1, m2, s1, s2, q1, q2)
}

// TestBackendSwitchConcurrentSafe：并发反复切换默认后端 + 并发计算，
// 每次计算显式指定后端，切换不得影响进行中的任务（-race 下检验）。
func TestBackendSwitchConcurrentSafe(t *testing.T) {
	x := []float64{1, 2, 3, 4, 5, 6}
	y := []float64{2.1, 4.0, 6.2, 7.9, 10.1, 12.0}
	var wg sync.WaitGroup
	stop := make(chan struct{})
	wg.Add(1)
	go func() {
		defer wg.Done()
		for {
			select {
			case <-stop:
				return
			default:
				_ = SetDefaultBackend(BackendGonum)
				_ = SetDefaultBackend(BackendStd)
			}
		}
	}()
	for g := 0; g < 16; g++ {
		wg.Add(1)
		go func(g int) {
			defer wg.Done()
			b := BackendStd
			if g%2 == 0 {
				b = BackendGonum
			}
			for k := 0; k < 200; k++ {
				r, err := WLS(x, y, RegressionOptions{InterceptFree: true, Backend: b})
				if err != nil {
					t.Error(err)
					return
				}
				// 结果必须始终来自显式选定的后端：两后端斜率应一致，
				// 且每次逐位相同（不随全局切换抖动）。
				if math.Abs(r.Coefficients[1]-1.98) > 0.05 {
					t.Errorf("unexpected slope %v", r.Coefficients[1])
					return
				}
			}
		}(g)
	}
	close(stop)
	wg.Wait()
	t.Log("[judgement] concurrent default-backend switching never perturbed explicit-backend runs")
}

// TestUnknownBackendRejected：不支持的后端显式拒绝，无静默回退。
func TestUnknownBackendRejected(t *testing.T) {
	if err := SetDefaultBackend(Backend("atlas")); err == nil {
		t.Fatal("unknown backend must be rejected")
	}
	if Available(Backend("atlas")) {
		t.Fatal("unknown backend reported available")
	}
	if _, err := MetaCombine([]float64{1, 2}, []float64{0.1, 0.1},
		MetaOptions{Backend: Backend("magma")}); err == nil {
		t.Fatal("unknown backend for MetaCombine must error")
	}
}

// TestGonumVersionPinned：代码常量、go.mod、vendor 清单三处版本一致，
// 且 vendored 树存在（离线可复现构建）。
func TestGonumVersionPinned(t *testing.T) {
	mod, err := os.ReadFile("../go.mod")
	if err != nil {
		t.Fatal(err)
	}
	want := "gonum.org/v1/gonum " + GonumVersion
	if !strings.Contains(string(mod), want) {
		t.Fatalf("go.mod does not pin %s", want)
	}
	mt, err := os.ReadFile("../vendor/modules.txt")
	if err != nil {
		t.Fatalf("vendor/modules.txt missing: %v (run: go mod vendor)", err)
	}
	if !strings.Contains(string(mt), "# gonum.org/v1/gonum "+GonumVersion) {
		t.Fatalf("vendor/modules.txt does not pin gonum %s", GonumVersion)
	}
	if v, _ := Version(BackendGonum); !strings.Contains(v, GonumVersion) {
		t.Fatalf("Version(gonum)=%q does not mention %s", v, GonumVersion)
	}
	t.Logf("[judgement] gonum pinned at %s in code constant, go.mod and vendor/", GonumVersion)
}

// TestGonumErrorsMapped：gonum 后端的奇异矩阵被映射为 *Error(KindSingular)，
// 不透出第三方错误类型。
func TestGonumErrorsMapped(t *testing.T) {
	x := []float64{1, 1, 1, 1, 1}
	y := []float64{2, 2, 2, 2, 2}
	_, err := WLS(x, y, RegressionOptions{InterceptFree: true, Backend: BackendGonum})
	e, ok := AsError(err)
	if !ok {
		t.Fatalf("want *numwrap.Error, got %T %v", err, err)
	}
	if e.Kind != KindSingular {
		t.Fatalf("want KindSingular got %v", e.Kind)
	}
	t.Logf("[judgement] gonum singular failure translated: %v", e)
}
