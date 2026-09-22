// Package stgerr defines the error taxonomy for the statistical genetics
// pipeline. Every error that crosses a package boundary is a *stgerr.Error
// carrying a stable Kind so callers can distinguish rejection reasons
// programmatically instead of matching on message text.
package stgerr

import (
	"errors"
	"fmt"
	"strings"
)

// Kind is the stable, machine-readable category of a pipeline error.
type Kind string

const (
	// KindMalformedInput: syntax problem (unparsable line, wrong column count).
	KindMalformedInput Kind = "malformed_input"
	// KindMissingField: a required field is absent or empty.
	KindMissingField Kind = "missing_field"
	// KindInvalidValue: field present but value illegal (NaN, Inf, negative SE, ...).
	KindInvalidValue Kind = "invalid_value"
	// KindDuplicateVariant: same variant key seen more than once in one source.
	KindDuplicateVariant Kind = "duplicate_variant"
	// KindUnknownAllele: allele contains a base other than A/C/G/T (and case forms).
	KindUnknownAllele Kind = "unknown_allele"
	// KindStrandAmbiguous: A/T or C/G SNPs where strand cannot be inferred.
	KindStrandAmbiguous Kind = "strand_ambiguous"
	// KindIncompatibleAlleles: alleles cannot be aligned (no swap/flip/complement match).
	KindIncompatibleAlleles Kind = "incompatible_alleles"
	// KindVariantNotFound: lookup in a variant-indexed source failed.
	KindVariantNotFound Kind = "variant_not_found"
	// KindOutOfBounds: estimate outside its interpretable range (e.g. h2 < 0 or > 1).
	KindOutOfBounds Kind = "out_of_bounds"
	// KindUninterpretableInterval: confidence interval not contained in the valid range.
	KindUninterpretableInterval Kind = "uninterpretable_interval"
	// KindCohortOverlap: cohorts share samples and must not be treated as independent.
	KindCohortOverlap Kind = "cohort_overlap"
	// KindDuplicateCohort: same cohort supplied more than once.
	KindDuplicateCohort Kind = "duplicate_cohort"
	// KindUnsupportedOption: caller requested an option the package does not support.
	KindUnsupportedOption Kind = "unsupported_option"
	// KindLimitExceeded: configured row/sample/variant/time budget exceeded.
	KindLimitExceeded Kind = "limit_exceeded"
	// KindTruncated: input accepted but cut off at a configured limit.
	KindTruncated Kind = "truncated"
	// KindLibraryFailure: wrapped third-party numerical/statistical library failure.
	KindLibraryFailure Kind = "library_failure"
	// KindConcurrency: worker pool / shared-state failure.
	KindConcurrency Kind = "concurrency"
	// KindInternal: invariant violation that should never happen.
	KindInternal Kind = "internal"
)

// Error is the single structured error type of the pipeline.
type Error struct {
	Kind   Kind   // stable category
	Source string // input name, cohort name or package component
	Reason string // human-readable explanation
	Cause  error  // wrapped underlying error, possibly another *Error
}

func (e *Error) Error() string {
	var b strings.Builder
	b.WriteString("statgen: ")
	b.WriteString(string(e.Kind))
	if e.Source != "" {
		fmt.Fprintf(&b, " [%s]", e.Source)
	}
	if e.Reason != "" {
		b.WriteString(": ")
		b.WriteString(e.Reason)
	}
	if e.Cause != nil {
		b.WriteString(": ")
		b.WriteString(e.Cause.Error())
	}
	return b.String()
}

// Unwrap exposes the cause for errors.Is/errors.As.
func (e *Error) Unwrap() error { return e.Cause }

// New builds an *Error without a cause.
func New(kind Kind, source, reason string) *Error {
	return &Error{Kind: kind, Source: source, Reason: reason}
}

// Wrap builds an *Error wrapping cause.
func Wrap(kind Kind, source, reason string, cause error) *Error {
	return &Error{Kind: kind, Source: source, Reason: reason, Cause: cause}
}

// Is reports whether err (or any wrapped error) has the given kind.
func Is(err error, kind Kind) bool {
	var se *Error
	if errors.As(err, &se) {
		if se.Kind == kind {
			return true
		}
	}
	return false
}

// AsError returns the first *Error in err's chain, or nil.
func AsError(err error) *Error {
	var se *Error
	if errors.As(err, &se) {
		return se
	}
	return nil
}
