// Package ldsc implements LD-score regression (Bulik-Sullivan et al. 2015):
//
//	E[chi2] = intercept + (N / M) * h2 * l2
//
// Implementation notes and model assumptions:
//   - chi2 = (beta/SE)^2 is invariant to beta sign, but beta is aligned
//     anyway so downstream consumers see one fixed allele convention.
//   - weights follow the LDSC-style two-step approximation
//     w_i = 1 / (2*(1 + b*l2)^2), recomputed once with b = N*h2/M.
//   - negative, >1 and CI-out-of-range estimates are never clipped
//     silently; they carry an explicit Status.
//   - all inputs are sorted by canonical key; parallel preparation is merged
//     back in that fixed order, so Workers>1 is bit-identical to serial.
package ldsc

import (
	"fmt"
	"math"
	"sort"
	"sync"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/ldscore"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/numlib"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/sumstats"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/variant"
)

// Status classifies the interpretability of an h2 estimate.
type Status string

const (
	// StatusValid: point estimate in [0,1] and CI within [0,1].
	StatusValid Status = "valid"
	// StatusNegative: point estimate below zero (data consistent with zero).
	StatusNegative Status = "negative"
	// StatusAboveOne: point estimate above one (model/data inconsistency).
	StatusAboveOne Status = "above_one"
	// StatusCIOutOfBounds: point within [0,1] but CI leaves the range.
	StatusCIOutOfBounds Status = "ci_out_of_bounds"
	// StatusUnidentified: regression could not estimate h2 (degenerate data).
	StatusUnidentified Status = "unidentified"
)

// NegativePolicy controls negative point estimates.
type NegativePolicy string

const (
	// NegativeReport keeps the negative value and flags StatusNegative.
	NegativeReport NegativePolicy = "report"
	// NegativeZero truncates to zero with the flag preserved.
	NegativeZero NegativePolicy = "truncate_zero"
)

// Config controls the regression.
type Config struct {
	Name           string
	M              int     // score-normalizing SNP count M; 0 -> matched count
	Alpha          float64 // two-sided significance level (default 0.05)
	NegPolicy      NegativePolicy
	ConstrainedInt bool    // fix intercept at 1
	MaxChi2        float64 // LDSC outlier floor (default 80)
	MinVariants    int     // minimum matched variants (default 2 unconstrained)
	MaxVariants    int     // input-size cap; 0 = default 1,000,000
	Workers        int     // 1 = serial; >1 deterministic parallel matching
}

// Result is the full heritability estimate with uncertainty and provenance.
type Result struct {
	Cohort      string
	H2          float64
	SE          float64
	CI          [2]float64
	Intercept   float64
	InterceptSE float64
	Status      Status
	NVariants   int
	R2          float64
	Weights     []float64 // aligned with Variants
	Variants    []variant.Key
	Report      variant.Report
}

// Estimate runs LD-score regression on matched cohort/score data.
func Estimate(cohort *sumstats.CohortSummary, scores *ldscore.Table,
	cfg Config) (*Result, error) {

	if cohort == nil || scores == nil {
		return nil, stgerr.New(stgerr.KindMalformedInput, "ldsc", "nil cohort or LD table")
	}
	if cfg.Alpha == 0 {
		cfg.Alpha = 0.05
	}
	if cfg.Alpha <= 0 || cfg.Alpha >= 1 {
		return nil, stgerr.New(stgerr.KindInvalidValue, "ldsc",
			fmt.Sprintf("alpha must be in (0,1), got %v", cfg.Alpha))
	}
	if cfg.MaxChi2 == 0 {
		cfg.MaxChi2 = 80
	}
	if cfg.MinVariants == 0 {
		if cfg.ConstrainedInt {
			cfg.MinVariants = 2
		} else {
			cfg.MinVariants = 3
		}
	}
	maxVariants := cfg.MaxVariants
	if maxVariants == 0 {
		maxVariants = 1_000_000
	}
	name := cfg.Name
	if name == "" {
		name = cohort.Name
	}
	aligner := variant.DefaultAligner()

	matched, report := matchParallel(cohort, scores, aligner, cfg.Workers)
	if len(matched) > maxVariants {
		return nil, stgerr.New(stgerr.KindLimitExceeded, name,
			fmt.Sprintf("matched variants %d exceed limit %d", len(matched), maxVariants))
	}

	res := &Result{Cohort: name, Report: report,
		Variants: make([]variant.Key, 0, len(matched))}

	// LDSC chi2 outlier rule: drop variants with chi2 > MaxChi2 (rare
	// association artifacts that destabilize the regression).
	kept := make([]MatchedVariant, 0, len(matched))
	for _, mv := range matched {
		if mv.Chi2 > cfg.MaxChi2 {
			report.Add(variant.Rejection{
				Source: name, Key: string(mv.Key),
				Kind:   stgerr.KindInvalidValue,
				Reason: fmt.Sprintf("chi2=%.2f exceeds max %.2f", mv.Chi2, cfg.MaxChi2),
			})
			continue
		}
		kept = append(kept, mv)
		res.Variants = append(res.Variants, mv.Key)
	}
	res.NVariants = len(kept)
	res.Report = report

	if len(kept) < cfg.MinVariants {
		res.Status = StatusUnidentified
		return res, stgerr.New(stgerr.KindOutOfBounds, name,
			fmt.Sprintf("only %d matched variants after filtering (need >= %d)",
				len(kept), cfg.MinVariants))
	}

	m := cfg.M
	if m <= 0 {
		m = len(kept)
	}
	n := cohort.DefaultN
	if n <= 0 {
		// Fall back to the first explicit per-variant N; require consistency.
		for _, mv := range kept {
			if mv.Stats.N > 0 {
				n = mv.Stats.N
				break
			}
		}
	}
	if n <= 0 {
		return nil, stgerr.New(stgerr.KindMissingField, name,
			"sample size N unavailable (neither cohort default nor per-variant N)")
	}

	// Effective per-variant sample size: prefer per-variant N when given.
	nEff := make([]float64, len(kept))
	for i, mv := range kept {
		if mv.Stats.N > 0 {
			nEff[i] = mv.Stats.N
		} else {
			nEff[i] = n
		}
	}

	y := make([]float64, len(kept))
	l2 := make([]float64, len(kept))
	for i, mv := range kept {
		y[i] = mv.Chi2
		l2[i] = mv.Scores.L2
	}

	// Two-step weights: initial b=0 (weights all 1/2), refit with b1.
	b1, b0 := regressOnce(kept, y, l2, nEff, m, nil, cfg.ConstrainedInt)
	weights := make([]float64, len(kept))
	for i := range kept {
		weights[i] = ldscWeight(nEff[i], l2[i], b1, b0)
	}
	slope, intercept := regressOnce(kept, y, l2, nEff, m, weights, cfg.ConstrainedInt)

	// Recompute standard errors from the weighted fit covariance.
	X := make([][]float64, len(kept))
	yy := make([]float64, len(kept))
	for i := range kept {
		x := nEff[i] / float64(m) * l2[i]
		if cfg.ConstrainedInt {
			X[i] = []float64{x}
			yy[i] = y[i] - 1
		} else {
			X[i] = []float64{1, x}
			yy[i] = y[i]
		}
	}
	fit, ferr := numlib.WeightedOLS(X, yy, weights,
		numlib.SolveOptions{Method: numlib.SolveWeightedOLS})
	if ferr != nil {
		res.Status = StatusUnidentified
		return res, ferr
	}

	var h2, h2SE, inter, interSE float64
	if cfg.ConstrainedInt {
		h2, h2SE = slope, math.Sqrt(fit.Cov[0])
		inter, interSE = 1, 0
	} else {
		inter, interSE = intercept, math.Sqrt(fit.Cov[0])
		h2 = slope
		h2SE = math.Sqrt(fit.Cov[1*fit.P+1])
	}
	res.H2 = h2
	res.SE = math.Max(h2SE, 0)
	res.Intercept = inter
	res.InterceptSE = interSE
	res.R2 = fit.R2
	res.Weights = weights

	zcrit := math.Abs(numlib.NormalQuantile(1 - cfg.Alpha/2))
	res.CI = [2]float64{h2 - zcrit*res.SE, h2 + zcrit*res.SE}
	res.Status = classify(h2, res.CI)

	if res.Status == StatusNegative && cfg.NegPolicy == NegativeZero {
		res.H2 = 0
		res.CI = [2]float64{math.Max(res.CI[0], 0), math.Max(res.CI[1], 0)}
	}
	return res, nil
}

// classify maps an estimate and its CI to the explicit interpretability status.
func classify(h2 float64, ci [2]float64) Status {
	switch {
	case math.IsNaN(h2) || math.IsInf(h2, 0):
		return StatusUnidentified
	case h2 < 0:
		return StatusNegative
	case h2 > 1:
		return StatusAboveOne
	case ci[0] < 0 || ci[1] > 1:
		return StatusCIOutOfBounds
	default:
		return StatusValid
	}
}

// ldscWeight approximates the inverse variance of chi2:
// Var(chi2) ~ 2*(intercept + slope*l2)^2 with slope = N*h2/M.
func ldscWeight(nEff, l2, slope, intercept float64) float64 {
	b := nEff * slope * l2
	v := 2 * (intercept + b) * (intercept + b)
	if !(v > 0) || math.IsInf(v, 0) {
		return 0
	}
	return 1 / v
}

// regressOnce performs the slope extraction; weights nil => equal weights.
// Returns (slope, intercept) where slope already equals h2 (the x is N*l2/M).
func regressOnce(kept []MatchedVariant, y, l2, nEff []float64, m int,
	weights []float64, constrained bool) (slope, intercept float64) {

	n := len(kept)
	X := make([][]float64, n)
	yy := make([]float64, n)
	for i := 0; i < n; i++ {
		x := nEff[i] / float64(m) * l2[i]
		if constrained {
			X[i] = []float64{x}
			yy[i] = y[i] - 1
		} else {
			X[i] = []float64{1, x}
			yy[i] = y[i]
		}
	}
	fit, err := numlib.WeightedOLS(X, yy, weights,
		numlib.SolveOptions{Method: numlib.SolveWeightedOLS})
	if err != nil {
		return 0, 0
	}
	if constrained {
		return fit.Coeffs[0], 1
	}
	return fit.Coeffs[1], fit.Coeffs[0]
}

// matchParallel aligns cohort variants to LD scores. With Workers>1 the
// sorted key space is split into contiguous shards; each shard's output is
// written to its own slot (no shared mutation), then concatenated in order,
// producing bit-identical output to the serial path.
func matchParallel(cohort *sumstats.CohortSummary, scores *ldscore.Table,
	a *variant.Aligner, workers int) ([]MatchedVariant, variant.Report) {

	keys := make([]variant.Key, 0, len(cohort.Variants))
	for k := range cohort.Variants {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })

	if workers <= 1 {
		return matchShard(keys, cohort, scores, a)
	}
	if workers > len(keys) {
		workers = len(keys)
	}
	type shard struct {
		mv []MatchedVariant
		r  variant.Report
	}
	slots := make([]shard, workers)
	var wg sync.WaitGroup
	for w := 0; w < workers; w++ {
		lo := w * len(keys) / workers
		hi := (w + 1) * len(keys) / workers
		wg.Add(1)
		go func(wid, lo, hi int) {
			defer wg.Done()
			mv, r := matchShard(keys[lo:hi], cohort, scores, a)
			slots[wid] = shard{mv: mv, r: r}
		}(w, lo, hi)
	}
	wg.Wait()

	var out []MatchedVariant
	merged := variant.Report{Source: cohort.Name + "~" + scores.Source,
		RejectByKind: map[stgerr.Kind]int{}}
	for _, s := range slots {
		out = append(out, s.mv...)
		merged.Merge(s.r)
	}
	return out, merged
}

// matchShard aligns a contiguous key slice.
func matchShard(keys []variant.Key, cohort *sumstats.CohortSummary,
	scores *ldscore.Table, a *variant.Aligner) ([]MatchedVariant, variant.Report) {

	report := variant.Report{Source: cohort.Name + "~" + scores.Source,
		RejectByKind: map[stgerr.Kind]int{}}
	out := make([]MatchedVariant, 0, len(keys))
	for _, k := range keys {
		st := cohort.Variants[k]
		sc, ok := scores.Scores[k]
		if !ok {
			report.Add(variant.Rejection{
				Source: cohort.Name, Key: string(k),
				Kind:   stgerr.KindVariantNotFound,
				Reason: "no LD score for cohort variant",
			})
			continue
		}
		al, err := a.Align(st.Var, sc.Var)
		if err != nil {
			report.AddError(0, string(k), err)
			continue
		}
		beta := st.Beta * al.Sign
		z := beta / st.SE
		out = append(out, MatchedVariant{
			Key: k, Stats: st, Scores: sc, Aligned: al,
			Beta: beta, Chi2: z * z,
		})
	}
	report.Accepted = len(out)
	return out, report
}
