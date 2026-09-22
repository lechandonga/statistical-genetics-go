package meta

import (
	"fmt"
	"sort"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
)

// OverlapKind distinguishes duplicates from partial sample overlap.
type OverlapKind string

const (
	OverlapDuplicate OverlapKind = "duplicate_cohort"
	OverlapSamples   OverlapKind = "sample_overlap"
)

// Overlap describes two cohorts sharing samples or a duplicate submission.
type Overlap struct {
	A, B   string
	Shared int
	Kind   OverlapKind
}

// prepare validates cohorts and applies the overlap policy. Duplicate cohorts
// are identified by equal OverlapKey; partial overlap by intersecting the
// (de-identified) SampleID sets. Findings are deterministic: pairs are
// reported in cohort-index order.
func prepare(cohorts []Cohort, cfg Config) (*Analysis, error) {
	if len(cohorts) == 0 {
		return nil, stgerr.New(stgerr.KindMissingField, "meta", "no cohorts supplied")
	}
	min := cfg.MinCohorts
	if min == 0 {
		min = 1
	}
	policy := cfg.OverlapPolicy
	if policy == "" {
		policy = OverlapFlag
	}
	if policy != OverlapReject && policy != OverlapFlag && policy != OverlapDeduplicate {
		return nil, stgerr.New(stgerr.KindUnsupportedOption, "meta",
			"unsupported overlap policy: "+string(policy))
	}

	an := &Analysis{Cohorts: cohorts}
	for i := range cohorts {
		if cohorts[i].Summary == nil {
			return nil, stgerr.New(stgerr.KindMissingField, "meta",
				fmt.Sprintf("cohort index %d has nil summary", i))
		}
		if cohorts[i].Summary.Name == "" {
			cohorts[i].Summary.Name = fmt.Sprintf("cohort_%d", i+1)
		}
	}

	seenKey := map[string]int{}
	sampleSets := make([]map[string]struct{}, len(cohorts))
	for i := range cohorts {
		if len(cohorts[i].SampleIDs) > 0 {
			sampleSets[i] = make(map[string]struct{}, len(cohorts[i].SampleIDs))
			for _, id := range cohorts[i].SampleIDs {
				sampleSets[i][id] = struct{}{}
			}
		}
	}

	dropped := map[int]bool{}
	for i := range cohorts {
		key := cohorts[i].OverlapKey
		if key == "" {
			key = cohorts[i].Summary.Name
		}
		if j, ok := seenKey[key]; ok {
			an.Overlaps = append(an.Overlaps, Overlap{
				A: cohorts[j].Summary.Name, B: cohorts[i].Summary.Name,
				Shared: -1, Kind: OverlapDuplicate})
			if policy == OverlapReject {
				return nil, stgerr.New(stgerr.KindDuplicateCohort, "meta",
					fmt.Sprintf("cohort %q duplicates %q",
						cohorts[i].Summary.Name, cohorts[j].Summary.Name))
			}
			if policy == OverlapDeduplicate {
				dropped[i] = true
				an.Dropped = append(an.Dropped, cohorts[i].Summary.Name)
			}
			continue
		}
		seenKey[key] = i
	}

	for i := 0; i < len(cohorts); i++ {
		for j := i + 1; j < len(cohorts); j++ {
			if dropped[i] || dropped[j] {
				continue
			}
			if sampleSets[i] == nil || sampleSets[j] == nil {
				continue
			}
			shared := intersectSize(sampleSets[i], sampleSets[j])
			if shared > 0 {
				an.Overlaps = append(an.Overlaps, Overlap{
					A: cohorts[i].Summary.Name, B: cohorts[j].Summary.Name,
					Shared: shared, Kind: OverlapSamples})
				if policy == OverlapReject {
					return nil, stgerr.New(stgerr.KindCohortOverlap, "meta",
						fmt.Sprintf("cohorts %q and %q share %d samples",
							cohorts[i].Summary.Name, cohorts[j].Summary.Name, shared))
				}
				if policy == OverlapDeduplicate {
					dropped[j] = true
					an.Dropped = append(an.Dropped, cohorts[j].Summary.Name)
				}
			}
		}
	}

	for i := range cohorts {
		if !dropped[i] {
			an.Effective = append(an.Effective, i)
		}
	}
	if len(an.Effective) < min {
		return nil, stgerr.New(stgerr.KindCohortOverlap, "meta",
			fmt.Sprintf("only %d effective cohorts after overlap policy (min %d)",
				len(an.Effective), min))
	}
	sort.Slice(an.Overlaps, func(i, j int) bool {
		if an.Overlaps[i].A != an.Overlaps[j].A {
			return an.Overlaps[i].A < an.Overlaps[j].A
		}
		return an.Overlaps[i].B < an.Overlaps[j].B
	})
	return an, nil
}

func intersectSize(a, b map[string]struct{}) int {
	if len(a) > len(b) {
		a, b = b, a
	}
	n := 0
	for k := range a {
		if _, ok := b[k]; ok {
			n++
		}
	}
	return n
}
