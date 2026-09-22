// Package meta implements multi-cohort fixed-effect (inverse-variance /
// Z-score) and random-effects (DerSimonian-Laird) meta-analysis with
// interpretable weights and heterogeneity statistics (Cochran's Q, I2,
// tau2, H). Cohort sample overlap and duplicate cohorts are detected from
// declared sample identifiers and flagged: overlapping evidence is never
// stacked as if independent.
package meta

import (
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/sumstats"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/variant"
)

// Model selects the pooling model.
type Model string

const (
	// FixedInverseVariance is classic inverse-variance fixed-effect pooling.
	FixedInverseVariance Model = "fixed_iv"
	// FixedZScore pools per-cohort Z scores (weights sqrt(N)).
	FixedZScore Model = "fixed_z"
	// RandomDerSimonianLaird adds a tau2 between-cohort variance.
	RandomDerSimonianLaird Model = "random_dl"
)

// Cohort wraps parsed summary stats with overlap-declaration metadata.
type Cohort struct {
	Summary    *sumstats.CohortSummary
	SampleIDs  []string // hashed/de-identified sample identifiers, if known
	OverlapKey string   // cohort provenance key; equal keys => duplicate cohort
}

// CohortEffect is one cohort's contribution at a variant.
type CohortEffect struct {
	Cohort  string
	Beta    float64
	SE      float64
	Z       float64
	Weight  float64
	N       float64
	Aligned variant.Aligned
}

// VariantResult is the pooled estimate at one variant.
type VariantResult struct {
	Key      variant.Key
	Beta     float64
	SE       float64
	Z        float64
	PValue   float64
	Model    Model
	Q        float64
	QDf      int
	QPValue  float64
	I2       float64
	Tau2     float64
	H        float64
	Effects  []CohortEffect
	Warnings []string
}

// Heterogeneity summarizes cross-variant heterogeneity bookkeeping.
type Heterogeneity struct {
	VariantCount int
}

// Config controls pooling and overlap handling.
type Config struct {
	Model         Model
	OverlapPolicy OverlapPolicy
	MinCohorts    int
	Workers       int
}

// OverlapPolicy decides how overlapping/duplicate cohorts are treated.
type OverlapPolicy string

const (
	// OverlapReject refuses analysis when overlap is detected.
	OverlapReject OverlapPolicy = "reject"
	// OverlapFlag keeps cohorts but marks every affected result.
	OverlapFlag OverlapPolicy = "flag"
	// OverlapDeduplicate drops later duplicate/overlapping cohorts.
	OverlapDeduplicate OverlapPolicy = "deduplicate"
)

// Analysis holds validated cohorts and overlap findings.
type Analysis struct {
	Cohorts   []Cohort
	Overlaps  []Overlap
	Dropped   []string
	Effective []int // indices into Cohorts after policy application
}
