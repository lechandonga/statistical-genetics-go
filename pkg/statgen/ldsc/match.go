package ldsc

import (
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/ldscore"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/sumstats"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/variant"
)

// MatchedVariant is one cohort/score pair after allele alignment.
type MatchedVariant struct {
	Key     variant.Key
	Stats   sumstats.VariantStats
	Scores  ldscore.Scores
	Aligned variant.Aligned
	Beta    float64 // orientation-corrected beta relative to the target alleles
	Chi2    float64 // orientation-corrected chi-square, sign-invariant
}

// MatchVariants aligns cohort and LD-score variants under a's policy using
// the deterministic serial path. (Estimate uses the parallel-capable path;
// both consume keys in sorted order and yield identical results.)
func MatchVariants(cohort *sumstats.CohortSummary, scores *ldscore.Table,
	a *variant.Aligner) ([]MatchedVariant, variant.Report) {
	return matchParallel(cohort, scores, a, 1)
}
