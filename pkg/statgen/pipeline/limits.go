// Package pipeline wires parsing, alignment, LD-score heritability,
// meta-analysis and PRS into one reproducible local run with configurable
// size/time budgets. Budget exhaustion aborts cleanly: no half-updated
// result state and no leaked goroutines or files.
package pipeline

import (
	"fmt"
	"sync"
	"time"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
)

// Limits bounds input size and compute time.
type Limits struct {
	MaxVariants     int
	MaxSamples      int
	MaxCohorts      int
	MaxFileBytes    int64
	Timeout         time.Duration
	TruncateInstead bool // truncate at the limit instead of rejecting
}

// DefaultLimits used when fields are zero.
var DefaultLimits = Limits{
	MaxVariants:  1_000_000,
	MaxSamples:   100_000,
	MaxCohorts:   64,
	MaxFileBytes: 512 << 20,
	Timeout:      10 * time.Minute,
}

// BudgetKind selects the bounded resource.
type BudgetKind int

const (
	BudgetVariants BudgetKind = iota
	BudgetSamples
	BudgetCohorts
	BudgetBytes
)

// Budget is a live concurrency-safe consumption counter.
type Budget struct {
	mu       sync.Mutex
	limits   Limits
	used     [4]int64
	deadline time.Time
	timedOut bool
}

// NewBudget constructs a budget, filling zero fields from defaults.
func NewBudget(l Limits) *Budget {
	if l.MaxVariants == 0 {
		l.MaxVariants = DefaultLimits.MaxVariants
	}
	if l.MaxSamples == 0 {
		l.MaxSamples = DefaultLimits.MaxSamples
	}
	if l.MaxCohorts == 0 {
		l.MaxCohorts = DefaultLimits.MaxCohorts
	}
	if l.MaxFileBytes == 0 {
		l.MaxFileBytes = DefaultLimits.MaxFileBytes
	}
	if l.Timeout == 0 {
		l.Timeout = DefaultLimits.Timeout
	}
	return &Budget{limits: l, deadline: time.Now().Add(l.Timeout)}
}

// Consume records n units of kind k; on overflow it returns KindLimitExceeded.
// A consumed budget is never partially rolled back: callers treat the error
// as terminal for the run, which is why the pipeline builds results in
// local state and only publishes a Result after success.
func (b *Budget) Consume(k BudgetKind, n int64) error {
	b.mu.Lock()
	defer b.mu.Unlock()
	var capv int64
	switch k {
	case BudgetVariants:
		capv = int64(b.limits.MaxVariants)
	case BudgetSamples:
		capv = int64(b.limits.MaxSamples)
	case BudgetCohorts:
		capv = int64(b.limits.MaxCohorts)
	case BudgetBytes:
		capv = b.limits.MaxFileBytes
	default:
		return stgerr.New(stgerr.KindInternal, "pipeline", "unknown budget kind")
	}
	if b.used[k]+n > capv {
		if b.limits.TruncateInstead && k == BudgetVariants {
			b.used[k] = capv
			return stgerr.New(stgerr.KindTruncated, "pipeline",
				fmt.Sprintf("resource %d truncated at cap %d", k, capv))
		}
		return stgerr.New(stgerr.KindLimitExceeded, "pipeline",
			fmt.Sprintf("resource %d would exceed cap %d (used %d + %d)",
				k, capv, b.used[k], n))
	}
	b.used[k] += n
	if time.Now().After(b.deadline) {
		b.timedOut = true
		return stgerr.New(stgerr.KindLimitExceeded, "pipeline",
			"configured timeout exceeded")
	}
	return nil
}

// Remaining reports the unused allowance for a resource.
func (b *Budget) Remaining(k BudgetKind) int64 {
	b.mu.Lock()
	defer b.mu.Unlock()
	switch k {
	case BudgetVariants:
		return int64(b.limits.MaxVariants) - b.used[k]
	case BudgetSamples:
		return int64(b.limits.MaxSamples) - b.used[k]
	case BudgetCohorts:
		return int64(b.limits.MaxCohorts) - b.used[k]
	case BudgetBytes:
		return b.limits.MaxFileBytes - b.used[k]
	}
	return 0
}
