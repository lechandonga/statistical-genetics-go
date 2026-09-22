package meta

import (
	"fmt"
	"math"
	"sort"
	"sync"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/numlib"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/variant"
)

// Prepare validates cohorts and applies the overlap policy.
func Prepare(cohorts []Cohort, cfg Config) (*Analysis, error) {
	return prepare(cohorts, cfg)
}

// unionKeys returns the sorted set of variant keys shared by >= minCohorts
// effective cohorts. Iteration always follows this order in serial and
// parallel modes.
func (a *Analysis) unionKeys(minCohorts int) []variant.Key {
	count := map[variant.Key]int{}
	for _, ci := range a.Effective {
		for k := range a.Cohorts[ci].Summary.Variants {
			count[k]++
		}
	}
	keys := make([]variant.Key, 0, len(count))
	for k, n := range count {
		if n >= minCohorts {
			keys = append(keys, k)
		}
	}
	sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })
	return keys
}

// PoolAll pools every variant present in >= MinCohorts effective cohorts.
func (a *Analysis) PoolAll(cfg Config) ([]VariantResult, error) {
	min := cfg.MinCohorts
	if min == 0 {
		min = 1
	}
	keys := a.unionKeys(min)
	out := make([]VariantResult, len(keys))
	workers := cfg.Workers
	if workers <= 1 || len(keys) < 64 {
		for i, k := range keys {
			r, err := a.Pool(k, cfg)
			if err != nil {
				return nil, err
			}
			out[i] = r
		}
		return out, nil
	}
	if workers > len(keys) {
		workers = len(keys)
	}
	var (
		wg       sync.WaitGroup
		firstErr error
		errMu    sync.Mutex
	)
	setErr := func(e error) {
		errMu.Lock()
		if firstErr == nil {
			firstErr = e
		}
		errMu.Unlock()
	}
	for w := 0; w < workers; w++ {
		lo := w * len(keys) / workers
		hi := (w + 1) * len(keys) / workers
		wg.Add(1)
		go func(lo, hi int) {
			defer wg.Done()
			for i := lo; i < hi; i++ {
				r, err := a.Pool(keys[i], cfg)
				if err != nil {
					setErr(err)
					return
				}
				out[i] = r // each slot written by exactly one goroutine
			}
		}(lo, hi)
	}
	wg.Wait()
	if firstErr != nil {
		return nil, firstErr
	}
	return out, nil
}

// Pool computes the pooled effect for one variant across effective cohorts.
// Effects are aligned to the first effective cohort's allele convention;
// every per-cohort contribution carries the alignment provenance.
func (a *Analysis) Pool(k variant.Key, cfg Config) (VariantResult, error) {
	model := cfg.Model
	if model == "" {
		model = FixedInverseVariance
	}
	if model != FixedInverseVariance && model != FixedZScore && model != RandomDerSimonianLaird {
		return VariantResult{}, stgerr.New(stgerr.KindUnsupportedOption, "meta",
			"unsupported model: "+string(model))
	}

	var base *variant.Variant
	var effects []CohortEffect
	for _, ci := range a.Effective {
		st, ok := a.Cohorts[ci].Summary.Variants[k]
		if !ok {
			continue
		}
		if base == nil {
			v := st.Var
			base = &v
			effects = append(effects, CohortEffect{
				Cohort: a.Cohorts[ci].Summary.Name,
				Beta:   st.Beta, SE: st.SE, Z: st.Beta / st.SE,
				Weight: 1 / (st.SE * st.SE), N: effectiveN(ci, a),
				Aligned: variant.Aligned{Target: st.Var, Orient: variant.OrientMatch, Sign: 1},
			})
			continue
		}
		al, err := variant.DefaultAligner().Align(*base, st.Var)
		if err != nil {
			// Strand-ambiguous or incompatible cohort evidence is excluded
			// for this variant with a recorded warning, not silently summed.
			effects = append(effects, CohortEffect{
				Cohort:  a.Cohorts[ci].Summary.Name,
				Aligned: variant.Aligned{Target: *base},
			})
			effects[len(effects)-1].Beta = math.NaN()
			_ = err
			continue
		}
		beta := st.Beta * al.Sign
		effects = append(effects, CohortEffect{
			Cohort: a.Cohorts[ci].Summary.Name,
			Beta:   beta, SE: st.SE, Z: beta / st.SE,
			Weight: 1 / (st.SE * st.SE), N: effectiveN(ci, a),
			Aligned: al,
		})
	}

	usable := make([]CohortEffect, 0, len(effects))
	var warnings []string
	for _, e := range effects {
		if math.IsNaN(e.Beta) {
			warnings = append(warnings,
				fmt.Sprintf("cohort %s excluded: alleles could not be aligned", e.Cohort))
			continue
		}
		usable = append(usable, e)
	}
	effects = usable

	for _, ov := range a.Overlaps {
		involved := map[string]bool{}
		for _, e := range effects {
			if e.Cohort == ov.A || e.Cohort == ov.B {
				involved[e.Cohort] = true
			}
		}
		if len(involved) == 2 {
			if ov.Kind == OverlapDuplicate {
				warnings = append(warnings,
					fmt.Sprintf("duplicate cohort evidence: %s and %s are not independent",
						ov.A, ov.B))
			} else {
				warnings = append(warnings,
					fmt.Sprintf("sample overlap %s~%s: %d shared samples; pooled SE ignores correlation",
						ov.A, ov.B, ov.Shared))
			}
		}
	}

	res := VariantResult{Key: k, Model: model, Effects: effects, Warnings: warnings}
	if len(effects) == 0 {
		return res, stgerr.New(stgerr.KindVariantNotFound, "meta",
			"no alignable cohort evidence for "+string(k))
	}
	if base != nil {
		// keep target for sign provenance
		for i := range res.Effects {
			if res.Effects[i].Aligned.Target.Chrom == "" {
				res.Effects[i].Aligned.Target = *base
			}
		}
	}

	// Cochran's Q on inverse-variance weights (shared by all models).
	feBeta, feSE := fixedIV(effects)
	q := 0.0
	for _, e := range effects {
		d := e.Beta - feBeta
		q += e.Weight * d * d
	}
	df := len(effects) - 1
	res.Q = q
	res.QDf = df
	if df > 0 {
		res.QPValue = 1 - numlib.ChiSquareCDF(q, float64(df))
	}
	// I2 and H.
	if df > 0 && q > float64(df) {
		res.I2 = (q - float64(df)) / q
	}
	if float64(df) > 0 {
		res.H = math.Sqrt(math.Max(q/float64(df), 0))
	}

	switch model {
	case FixedInverseVariance:
		res.Beta, res.SE = feBeta, feSE
	case FixedZScore:
		zSum, wSum := 0.0, 0.0
		for _, e := range effects {
			w := math.Sqrt(math.Max(e.N, 1))
			zSum += w * e.Z
			wSum += w
			e.Weight = w
		}
		// re-slice copy so Weight updates persist
		for i := range effects {
			res.Effects[i].Weight = math.Sqrt(math.Max(effects[i].N, 1))
		}
		z := zSum / wSum
		res.Z = z
		res.PValue = 2 * (1 - numlib.NormalCDF(math.Abs(z)))
		// Express pooled beta/SE consistent with z.
		res.SE = 1 / wSum
		res.Beta = z * res.SE
		res.Tau2 = 0
		return res, nil
	case RandomDerSimonianLaird:
		tau2 := dersimonianLaird(q, df, effects)
		res.Tau2 = tau2
		// Re-weight with tau2.
		var beta, wSum float64
		for _, e := range effects {
			w := 1 / (e.SE*e.SE + tau2)
			beta += w * e.Beta
			wSum += w
			for i := range res.Effects {
				if res.Effects[i].Cohort == e.Cohort {
					res.Effects[i].Weight = w
				}
			}
		}
		res.Beta = beta / wSum
		res.SE = 1 / math.Sqrt(wSum)
	}
	res.Z = res.Beta / res.SE
	res.PValue = 2 * (1 - numlib.NormalCDF(math.Abs(res.Z)))
	return res, nil
}

// fixedIV performs classic inverse-variance fixed-effect pooling.
func fixedIV(effects []CohortEffect) (beta, se float64) {
	var bSum, wSum float64
	for _, e := range effects {
		bSum += e.Weight * e.Beta
		wSum += e.Weight
	}
	if wSum <= 0 {
		return 0, math.Inf(1)
	}
	return bSum / wSum, 1 / math.Sqrt(wSum)
}

// dersimonianLaird estimates between-cohort variance tau2 = max(0,(Q-df)/C).
func dersimonianLaird(q float64, df int, effects []CohortEffect) float64 {
	if df <= 0 || q <= float64(df) {
		return 0
	}
	var wSum, wSqSum float64
	for _, e := range effects {
		wSum += e.Weight
		wSqSum += e.Weight * e.Weight
	}
	c := wSum - wSqSum/wSum
	if c <= 0 {
		return 0
	}
	return (q - float64(df)) / c
}

func effectiveN(ci int, a *Analysis) float64 {
	s := a.Cohorts[ci].Summary
	if s.DefaultN > 0 {
		return s.DefaultN
	}
	for _, st := range s.Variants {
		if st.N > 0 {
			return st.N
		}
	}
	return 0
}
