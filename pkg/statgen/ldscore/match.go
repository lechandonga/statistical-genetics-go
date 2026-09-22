package ldscore

import (
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/variant"
)

// MatchResult reports alignment of one LD-score row to a target variant.
type MatchResult struct {
	Key     variant.Key
	Score   Scores
	Aligned variant.Aligned
}

// Match aligns every table entry to target variants using the given aligner,
// returning matched entries in target order (deterministic) and a report of
// entries that could not be aligned.
func Match(t *Table, targets []variant.Variant, a *variant.Aligner) ([]MatchResult, variant.Report) {
	report := variant.Report{Source: t.Source, RejectByKind: map[stgerr.Kind]int{}}
	out := make([]MatchResult, 0, len(targets))
	for _, tv := range targets {
		sc, ok := t.Scores[tv.Key()]
		if !ok {
			report.Add(variant.Rejection{
				Source: t.Source, Key: string(tv.Key()),
				Kind:   stgerr.KindVariantNotFound,
				Reason: "LD score missing for variant",
			})
			continue
		}
		al, err := a.Align(tv, sc.Var)
		if err != nil {
			report.AddError(0, string(tv.Key()), err)
			continue
		}
		out = append(out, MatchResult{Key: tv.Key(), Score: sc, Aligned: al})
	}
	report.Accepted = len(out)
	return out, report
}
