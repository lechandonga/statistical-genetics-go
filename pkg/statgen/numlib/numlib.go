// Package numlib is the single seam to third-party numerical/statistical
// libraries (gonum). Packages above must not import gonum directly: syntax
// differences and internal panics are normalized here into stgerr errors
// with distinguishable kinds, and unsupported options are rejected
// explicitly before a library call is made.
package numlib

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/stat/distuv"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
)

// Backend names the wrapped numerical backend.
const Backend = "gonum"

// Version is the pinned backend module version.
const Version = "v0.17.0"

// SolverMethod enumerates the supported least-squares solvers.
type SolverMethod string

const (
	// SolveWeightedOLS fits weighted regression through gonum QR/Solve on the
	// weighted design matrix, with explicit singular-system detection.
	SolveWeightedOLS SolverMethod = "weighted_ols"
	// SolveQR fits the unweighted system via QR through gonum/mat.
	SolveQR SolverMethod = "qr"
)

// OLSResult is the outcome of a (weighted) linear regression.
type OLSResult struct {
	Coeffs     []float64
	ResidualSE float64
	Cov        []float64 // len p*p, covariance estimate sigma^2 (X'WX)^-1
	R2         float64
	N          int
	P          int
}

// SolveOptions selects the method.
type SolveOptions struct {
	Method SolverMethod
}

// WeightedOLS performs y = X b with per-observation weights w.
// w may be nil (all weights 1); otherwise len(w)==len(y), all w>=0.
// Degenerate systems return stgerr errors instead of gonum panics.
func WeightedOLS(x [][]float64, y, w []float64, opts SolveOptions) (res OLSResult, err error) {
	defer func() {
		if r := recover(); r != nil {
			err = stgerr.New(stgerr.KindLibraryFailure, "numlib",
				fmt.Sprintf("backend panic in WeightedOLS: %v", r))
		}
	}()

	if len(x) == 0 || len(x[0]) == 0 {
		return OLSResult{}, stgerr.New(stgerr.KindInvalidValue, "numlib", "empty design matrix")
	}
	n, p := len(x), len(x[0])
	if len(y) != n {
		return OLSResult{}, stgerr.New(stgerr.KindInvalidValue, "numlib",
			fmt.Sprintf("y length %d != rows %d", len(y), n))
	}
	method := opts.Method
	if method == "" {
		method = SolveWeightedOLS
	}
	if method != SolveWeightedOLS && method != SolveQR {
		return OLSResult{}, stgerr.New(stgerr.KindUnsupportedOption, "numlib",
			"unsupported solver method: "+string(method))
	}
	if method == SolveQR && w != nil {
		return OLSResult{}, stgerr.New(stgerr.KindUnsupportedOption, "numlib",
			"weights are only supported with with weighted_ols")
	}
	if w != nil && len(w) != n {
		return OLSResult{}, stgerr.New(stgerr.KindInvalidValue, "numlib",
			fmt.Sprintf("weights length %d != rows %d", len(w), n))
	}
	for i := 1; i < n; i++ {
		if len(x[i]) != p {
			return OLSResult{}, stgerr.New(stgerr.KindMalformedInput, "numlib", "ragged design matrix")
		}
	}
	for i := range y {
		if math.IsNaN(y[i]) || math.IsInf(y[i], 0) {
			return OLSResult{}, stgerr.New(stgerr.KindInvalidValue, "numlib",
				fmt.Sprintf("y[%d] NaN/Inf", i))
		}
		if w != nil && (math.IsNaN(w[i]) || w[i] < 0) {
			return OLSResult{}, stgerr.New(stgerr.KindInvalidValue, "numlib",
				fmt.Sprintf("w[%d] negative/NaN", i))
		}
	}

	// Effective design/response after pre-multiplication by sqrt(W).
	sqrtw := make([]float64, n)
	for i := 0; i < n; i++ {
		if w != nil {
			sqrtw[i] = math.Sqrt(w[i])
		} else {
			sqrtw[i] = 1
		}
	}
	Xm := mat.NewDense(n, p, nil)
	yv := mat.NewVecDense(n, nil)
	for i := 0; i < n; i++ {
		for j := 0; j < p; j++ {
			v := x[i][j]
			if math.IsNaN(v) || math.IsInf(v, 0) {
				return OLSResult{}, stgerr.New(stgerr.KindInvalidValue, "numlib",
					fmt.Sprintf("X[%d][%d] NaN/Inf", i, j))
			}
			Xm.Set(i, j, sqrtw[i]*v)
		}
		yv.SetVec(i, sqrtw[i]*y[i])
	}

	var qr mat.QR
	qr.Factorize(Xm)

	coef := mat.NewVecDense(p, nil)
	if err := qr.SolveVecTo(coef, false, yv); err != nil {
		return OLSResult{}, stgerr.New(stgerr.KindOutOfBounds, "numlib",
			"singular regression system: "+err.Error())
	}

	var pred mat.VecDense
	pred.MulVec(Xm, coef)
	rss, weightSum, active := 0.0, 0.0, 0
	for i := 0; i < n; i++ {
		ri := yv.AtVec(i) - pred.AtVec(i)
		wi := wgt(sqrtw[i])
		rss += ri * ri
		weightSum += wi
		if wi > 0 {
			active++
		}
	}
	if active <= p {
		return OLSResult{}, stgerr.New(stgerr.KindOutOfBounds, "numlib",
			fmt.Sprintf("degrees of freedom exhausted: active=%d <= p=%d", active, p))
	}
	// Residual degrees of freedom count active observations, as in R/lm.
	denom := float64(active - p)
	if denom <= 0 {
		return OLSResult{}, stgerr.New(stgerr.KindOutOfBounds, "numlib",
			"non-positive residual degrees of freedom")
	}
	sigma2 := rss / denom

	// Covariance: sigma^2 (X'WX)^-1, with an explicit singularity check.
	var XtWX mat.Dense
	XtWX.Mul(Xm.T(), Xm)
	var chol mat.Cholesky
	if ok := chol.Factorize(symmetricAdapter{m: &XtWX}); !ok {
		return OLSResult{}, stgerr.New(stgerr.KindOutOfBounds, "numlib",
			"X'WX is not positive definite")
	}
	eye := mat.NewDense(p, p, nil)
	for i := 0; i < p; i++ {
		eye.Set(i, i, 1)
	}
	var inv mat.Dense
	if serr := chol.SolveTo(&inv, eye); serr != nil {
		return OLSResult{}, stgerr.New(stgerr.KindLibraryFailure, "numlib",
			"covariance solve failed: "+serr.Error())
	}
	cov := make([]float64, p*p)
	for i := 0; i < p; i++ {
		for j := 0; j < p; j++ {
			cov[i*p+j] = sigma2 * inv.At(i, j)
		}
	}

	// Weighted R2 relative to the weighted-mean null model.
	swy := 0.0
	for i := 0; i < n; i++ {
		swy += wgt(sqrtw[i]) * y[i]
	}
	yMean := 0.0
	if weightSum > 0 {
		yMean = swy / weightSum
	}
	tss := 0.0
	for i := 0; i < n; i++ {
		d := yv.AtVec(i) - sqrtw[i]*yMean
		tss += d * d
	}
	r2 := 0.0
	if tss > 0 {
		r2 = 1 - rss/tss
	}

	coeffs := make([]float64, p)
	for j := 0; j < p; j++ {
		coeffs[j] = coef.AtVec(j)
	}
	return OLSResult{Coeffs: coeffs, ResidualSE: math.Sqrt(math.Max(sigma2, 0)),
		Cov: cov, R2: r2, N: n, P: p}, nil
}

func wgt(sqrtw float64) float64 { return sqrtw * sqrtw }

// NormalCDF is the standard normal cumulative distribution function.
func NormalCDF(z float64) float64 {
	return distuv.UnitNormal.CDF(z)
}

// NormalQuantile is the inverse standard normal CDF for p in (0,1).
// Out-of-range inputs return an error-bearing NaN sentinel documented to
// callers; callers in this pipeline always pass validated probabilities.
func NormalQuantile(p float64) float64 {
	if p <= 0 || p >= 1 || math.IsNaN(p) {
		return math.NaN()
	}
	return distuv.UnitNormal.Quantile(p)
}

// ChiSquareCDF returns P(Chi2(df) <= x).
func ChiSquareCDF(x, df float64) float64 {
	if x < 0 || df <= 0 {
		return math.NaN()
	}
	return distuv.ChiSquared{K: df}.CDF(x)
}

// Safe guards a backend call, converting panics into KindLibraryFailure.
func Safe(fn func() error) (err error) {
	defer func() {
		if r := recover(); r != nil {
			err = stgerr.New(stgerr.KindLibraryFailure, "numlib",
				fmt.Sprintf("backend panic: %v", r))
		}
	}()
	return fn()
}

// symmetricAdapter presents a symmetric X'WX Dense matrix as mat.Symmetric.
// X'WX is symmetric by construction; we read through At and supply the
// triangular slice as stored (lower-triangle duplication is not needed for
// Cholesky factorization, which consults Dim/At/T only).
type symmetricAdapter struct{ m *mat.Dense }

func (s symmetricAdapter) Dims() (int, int)    { r, c := s.m.Dims(); return r, c }
func (s symmetricAdapter) At(i, j int) float64 { return s.m.At(i, j) }
func (s symmetricAdapter) T() mat.Matrix       { return s }
func (s symmetricAdapter) SymmetricDim() int {
	r, _ := s.m.Dims()
	return r
}
