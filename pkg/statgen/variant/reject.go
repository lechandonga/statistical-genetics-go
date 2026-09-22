package variant

import (
	"fmt"
	"sort"
	"strings"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
)

// Rejection records one input row that was refused, with a stable reason.
// Rejected rows are never silently dropped: every parser returns a Report.
type Rejection struct {
	Source string
	Line   int
	Key    string
	Kind   stgerr.Kind
	Reason string
}

// Report aggregates accepted counts and per-row rejections.
type Report struct {
	Source       string
	Accepted     int
	Rejected     []Rejection
	Truncated    bool
	RejectByKind map[stgerr.Kind]int
}

// Add records a rejection and keeps the by-kind tally.
func (r *Report) Add(rj Rejection) {
	if rj.Source == "" {
		rj.Source = r.Source
	}
	r.Rejected = append(r.Rejected, rj)
	if r.RejectByKind == nil {
		r.RejectByKind = make(map[stgerr.Kind]int)
	}
	r.RejectByKind[rj.Kind]++
}

// AddError extracts an *Error kind from err and records it against a line.
func (r *Report) AddError(line int, key string, err error) {
	se := stgerr.AsError(err)
	if se == nil {
		r.Add(Rejection{Line: line, Key: key, Kind: stgerr.KindInternal, Reason: err.Error()})
		return
	}
	r.Add(Rejection{Line: line, Key: key, Kind: se.Kind, Reason: se.Reason})
}

// Merge folds another report (same kind of data, different shard) in.
func (r *Report) Merge(other Report) {
	r.Accepted += other.Accepted
	r.Rejected = append(r.Rejected, other.Rejected...)
	if r.RejectByKind == nil {
		r.RejectByKind = make(map[stgerr.Kind]int)
	}
	for k, n := range other.RejectByKind {
		r.RejectByKind[k] += n
	}
	if other.Truncated {
		r.Truncated = true
	}
}

// Sorted returns a line-ordered copy of the rejections (reports are built
// in read order already, but Merge from parallel shards requires a sort).
func (r *Report) Sorted() []Rejection {
	out := make([]Rejection, len(r.Rejected))
	copy(out, r.Rejected)
	sort.SliceStable(out, func(i, j int) bool {
		if out[i].Line != out[j].Line {
			return out[i].Line < out[j].Line
		}
		return out[i].Key < out[j].Key
	})
	return out
}

// Summary returns a deterministic one-line input summary for test logs:
// source name, accepted count, rejection counts grouped by stable kind,
// and a truncation flag.
func (r *Report) Summary() string {
	kinds := make([]string, 0, len(r.RejectByKind))
	for k := range r.RejectByKind {
		kinds = append(kinds, string(k))
	}
	sort.Strings(kinds)
	parts := make([]string, 0, len(kinds))
	for _, k := range kinds {
		parts = append(parts, fmt.Sprintf("%s=%d", k, r.RejectByKind[stgerr.Kind(k)]))
	}
	trunc := ""
	if r.Truncated {
		trunc = " truncated=true"
	}
	return fmt.Sprintf("source=%s accepted=%d rejected=%d%s rejects{%s}",
		r.Source, r.Accepted, len(r.Rejected), trunc, strings.Join(parts, ","))
}
