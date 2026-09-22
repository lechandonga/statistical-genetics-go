package numlib

import (
	"math"
	"testing"
)

func TestWeightedOLSKnownLine(t *testing.T) {
	// y = 2 + 3x, noiseless, 8 points, weights 1..8.
	x := [][]float64{}
	var y, w []float64
	for i := 0; i < 8; i++ {
		xi := float64(i)
		x = append(x, []float64{1, xi})
		y = append(y, 2+3*xi)
		w = append(w, float64(i+1))
	}
	res, err := WeightedOLS(x, y, w, SolveOptions{Method: SolveWeightedOLS})
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if math.Abs(res.Coeffs[0]-2) > 1e-9 || math.Abs(res.Coeffs[1]-3) > 1e-9 {
		t.Fatalf("coeffs = %v, want [2 3]", res.Coeffs)
	}
	if math.Abs(res.R2-1) > 1e-9 {
		t.Fatalf("R2 = %v, want 1", res.R2)
	}
	if res.ResidualSE > 1e-9 {
		t.Fatalf("residualSE = %v, want ~0", res.ResidualSE)
	}
}

func TestWeightedOLSSingular(t *testing.T) {
	x := [][]float64{{1, 2}, {1, 2}, {1, 2}}
	y := []float64{1, 2, 3}
	if _, err := WeightedOLS(x, y, nil, SolveOptions{}); err == nil {
		t.Fatal("expected singular-system error")
	}
}

func TestUnsupportedMethod(t *testing.T) {
	_, err := WeightedOLS([][]float64{{1}}, []float64{1}, nil, SolveOptions{Method: "magic"})
	if err == nil {
		t.Fatal("expected unsupported option error")
	}
}

func TestKahanAndCDF(t *testing.T) {
	if got := KahanSum([]float64{1e16, -1e16, 1}); math.Abs(got-1) > 1e-6 {
		t.Fatalf("Kahan sum = %v, want 1", got)
	}
	if got := NormalCDF(0); math.Abs(got-0.5) > 1e-12 {
		t.Fatalf("Phi(0)=%v", got)
	}
	if got := ChiSquareCDF(0, 1); math.Abs(got) > 1e-12 {
		t.Fatalf("chi2cdf(0)=%v", got)
	}
}
