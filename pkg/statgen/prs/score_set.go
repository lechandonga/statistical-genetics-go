package prs

import (
	"fmt"
	"math"
	"sort"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/meta"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/sumstats"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/variant"
)

// signConvention is the single, documented convention of every ScoreSet:
// weight is the change per additional copy of Target.Effect at that locus.
const signConvention = "weight = per-copy effect of Target.Effect allele"

// NewScoreSet validates effects: duplicate keys, non-finite weights and
// strand-ambiguous loci are rejected (ambiguous must be pre-marked through
// the Aligner), then returns an immutable set sorted by key.
func NewScoreSet(effects []Effect) (*ScoreSet, error) {
	if len(effects) == 0 {
		return nil, stgerr.New(stgerr.KindMissingField, "prs", "no scoring effects supplied")
	}
	seen := map[variant.Key]bool{}
	out := make([]Effect, 0, len(effects))
	for _, e := range effects {
		if math.IsNaN(e.Weight) || math.IsInf(e.Weight, 0) {
			return nil, stgerr.New(stgerr.KindInvalidValue, "prs",
				fmt.Sprintf("non-finite weight for %s", e.Target.Key()))
		}
		if e.Target.IsAmbiguous() {
			return nil, stgerr.New(stgerr.KindStrandAmbiguous, "prs",
				"ambiguous locus must be resolved before scoring: "+string(e.Target.Key()))
		}
		k := e.Target.Key()
		if seen[k] {
			return nil, stgerr.New(stgerr.KindDuplicateVariant, "prs",
				"duplicate scoring locus "+string(k))
		}
		seen[k] = true
		out = append(out, e)
	}
	sort.Slice(out, func(i, j int) bool { return out[i].Target.Key() < out[j].Target.Key() })
	return &ScoreSet{Effects: out, signConvention: signConvention}, nil
}

// FromCohort builds effects aligned to target variants (typically the
// genotype matrix variants) from one cohort. A loci missing from the cohort
// are simply absent from the set — callers that need them present should
// first compare keys; FromCohort never silently invents zero effects.
func FromCohort(targets []variant.Variant, c *sumstats.CohortSummary,
	a *variant.Aligner) (*ScoreSet, error) {

	if a == nil {
		a = variant.DefaultAligner()
	}
	var effects []Effect
	for _, tv := range targets {
		st, ok := c.Variants[tv.Key()]
		if !ok {
			continue
		}
		al, err := a.Align(tv, st.Var)
		if err != nil {
			return nil, err
		}
		if al.Marked {
			return nil, stgerr.New(stgerr.KindStrandAmbiguous, "prs",
				"ambiguous locus marked but cannot enter a score set: "+string(tv.Key()))
		}
		effects = append(effects, Effect{
			Target: tv, Weight: st.Beta * al.Sign, Source: c.Name,
		})
	}
	if len(effects) == 0 {
		return nil, stgerr.New(stgerr.KindVariantNotFound, "prs",
			"no cohort variants align to the supplied targets")
	}
	return NewScoreSet(effects)
}

// FromPooled builds effects from pooled meta-analysis results.
// Targets are looked up by key; pooled betas already share one convention
// because Pool aligned every cohort to its first-cohort target — here we
// realign onto the explicitly supplied targets so changing inputs can never
// silently change the sign convention.
func FromPooled(results []meta.VariantResult, targets []variant.Variant,
	a *variant.Aligner) (*ScoreSet, error) {

	if a == nil {
		a = variant.DefaultAligner()
	}
	byKey := make(map[variant.Key]meta.VariantResult, len(results))
	for _, r := range results {
		byKey[r.Key] = r
	}
	var effects []Effect
	for _, tv := range targets {
		r, ok := byKey[tv.Key()]
		if !ok {
			continue
		}
		// Pooled result's convention base: recover from first effect target.
		var base variant.Variant
		if len(r.Effects) > 0 {
			base = r.Effects[0].Aligned.Target
		}
		sign := 1.0
		if base.Chrom != "" {
			al, err := a.Align(tv, base)
			if err != nil {
				return nil, err
			}
			sign = al.Sign
		}
		effects = append(effects, Effect{
			Target: tv, Weight: r.Beta * sign, Source: string(r.Model),
		})
	}
	if len(effects) == 0 {
		return nil, stgerr.New(stgerr.KindVariantNotFound, "prs",
			"no pooled results align to the supplied targets")
	}
	return NewScoreSet(effects)
}

// Convention returns the documented sign convention for test logs.
func (s *ScoreSet) Convention() string { return s.signConvention }
