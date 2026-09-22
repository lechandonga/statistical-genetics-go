// Package prs constructs and evaluates individual polygenic risk scores:
// score(s) = sum_i w_i * dosage_i(s). Effect signs follow one fixed allele
// convention (the target variant's effect allele); changing the variant set
// never silently flips sign conventions or rescales scores. Missing
// genotypes and missing effect loci are handled by explicit policies.
package prs

import (
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/genotype"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/variant"
)

// Effect is one scoring locus with an effect defined relative to Target.
type Effect struct {
	Target variant.Variant
	Weight float64
	Source string // cohort or pooled-result provenance
}

// Scale controls score centering/scaling across the variant set.
type Scale string

const (
	// ScaleRaw: plain sum of w * dosage.
	ScaleRaw Scale = "raw"
	// ScaleMeanCenter: subtract per-locus weighted mean dosage.
	ScaleMeanCenter Scale = "mean_center"
	// ScaleStandardized: center then divide by score SD over samples.
	ScaleStandardized Scale = "standardized"
)

// Config controls score construction.
type Config struct {
	MissingGenotype genotype.DosagePolicy
	MissingLocus    MissingLocusPolicy
	ScoreScale      Scale
	Workers         int
}

// MissingLocusPolicy controls effect loci absent from the genotype matrix.
type MissingLocusPolicy string

const (
	// MissingLocusError refuses to score when any effect locus is absent.
	MissingLocusError MissingLocusPolicy = "error"
	// MissingLocusSkip drops absent loci and records them in the report.
	MissingLocusSkip MissingLocusPolicy = "skip"
)

// ScoreSet is an immutable, validated collection of scoring effects.
// Building it separately from scoring guarantees that changing the variant
// set yields an explicit new set rather than a silent convention drift.
type ScoreSet struct {
	Effects        []Effect
	signConvention string
}

// ScoreReport records what happened while computing scores.
type ScoreReport struct {
	UsedLoci       int
	SkippedLoci    []variant.Key
	MissingData    map[variant.Key]int // per-locus count of imputed calls
	SignConvention string
}

// Scores is per-sample output plus provenance.
type Scores struct {
	SampleIDs []string
	Values    []float64
	UsedLoci  int
	Report    ScoreReport
}
