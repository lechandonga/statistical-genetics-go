package genotype

import (
	"fmt"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
)

// DosagePolicy controls how missing genotypes are handled in scoring.
type DosagePolicy int

const (
	// DosageMeanImpute fills missing calls with the observed per-variant mean.
	DosageMeanImpute DosagePolicy = iota
	// DosageZeroImpute fills missing calls with the homozygous reference code 0.
	DosageZeroImpute
	// DosageRejectSample returns an error if any scored locus is missing for a sample.
	DosageRejectSample
)

// ValidateDosage checks a raw dosage value for the legal {0,1,2} set.
// The missing sentinel (genotype.MissingCode) is accepted as missing.
func ValidateDosage(x float64, missing float64) error {
	if x == missing {
		return nil
	}
	if x != 0 && x != 1 && x != 2 {
		return stgerr.New(stgerr.KindInvalidValue, "genotype",
			fmt.Sprintf("dosage must be one of {0,1,2,%g}, got %v", missing, x))
	}
	return nil
}

// MeanDosage computes the mean over observed calls, ignoring missing.
// With zero observed calls it returns 0 (zero-imputation fallback).
func MeanDosage(x []float64, missing []bool) float64 {
	sum, n := 0.0, 0
	for i := range x {
		if missing != nil && missing[i] {
			continue
		}
		sum += x[i]
		n++
	}
	if n == 0 {
		return 0
	}
	return sum / float64(n)
}

// ImputedRow applies a dosage policy, returning the row and the number of
// imputed entries. Under DosageRejectSample any missing entry errors.
func ImputedRow(x []float64, missing []bool, policy DosagePolicy) ([]float64, int, error) {
	out := make([]float64, len(x))
	copy(out, x)
	imputed := 0
	for i := range out {
		if missing != nil && missing[i] {
			switch policy {
			case DosageRejectSample:
				return nil, imputed, stgerr.New(stgerr.KindInvalidValue, "genotype",
					fmt.Sprintf("missing genotype at index %d rejected by policy", i))
			case DosageZeroImpute:
				out[i] = 0
			case DosageMeanImpute:
				out[i] = MeanDosage(x, missing)
			default:
				return nil, 0, stgerr.New(stgerr.KindUnsupportedOption, "genotype",
					"unknown dosage policy")
			}
			imputed++
		}
	}
	return out, imputed, nil
}
