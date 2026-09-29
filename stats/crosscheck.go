package stats

// 双后端交叉验证的结果核对：数值在 numwrap 声明容差内，
// 离散判定（Flags、排除明细的原因集合）必须完全一致。

import (
	"sort"

	"github.com/lechandonga/statistical-genetics-go/ingest"
	"github.com/lechandonga/statistical-genetics-go/numwrap"
)

const crossCheckOp = "crosscheck"

func crossCheckH2(a, b *H2Result) error {
	if err := numwrap.CrossCheckInt(crossCheckOp, "h2.nsnps", a.NSnps, b.NSnps); err != nil {
		return err
	}
	if a.Flags&^H2FlagClamped != b.Flags&^H2FlagClamped {
		return &numwrap.BackendMismatchError{Op: crossCheckOp,
			Detail: "h2 discrete flags differ between backends"}
	}
	checks := []struct {
		name string
		x, y float64
	}{
		{"h2raw", a.H2Raw, b.H2Raw},
		{"se", a.SE, b.SE},
		{"ci95lo", a.CI95[0], b.CI95[0]},
		{"ci95hi", a.CI95[1], b.CI95[1]},
		{"intercept", a.Intercept, b.Intercept},
		{"intercept_se", a.InterceptSE, b.InterceptSE},
	}
	for _, c := range checks {
		if err := numwrap.CrossCheckScalar(crossCheckOp, c.name, c.x, c.y); err != nil {
			return err
		}
	}
	return crossCheckRejects("h2.rejected", a.Rejected, b.Rejected)
}

// crossCheckRejects 比较两侧拒绝明细：同样的行必须给出同样的原因。
// 顺序在 h2 中由 sum 行序固定；这里再按 (Line,ID) 排序后比较，避免顺序偶然。
func crossCheckRejects(field string, a, b []ingest.RowRef) error {
	key := func(r ingest.RowRef) string {
		return string(r.Reason) + "|" + r.ID
	}
	sa := make([]string, len(a))
	sb := make([]string, len(b))
	for i, r := range a {
		sa[i] = key(r)
	}
	for i, r := range b {
		sb[i] = key(r)
	}
	sort.Strings(sa)
	sort.Strings(sb)
	if len(sa) != len(sb) {
		return &numwrap.BackendMismatchError{Op: crossCheckOp,
			Detail: field + " count differs"}
	}
	for i := range sa {
		if sa[i] != sb[i] {
			return &numwrap.BackendMismatchError{Op: crossCheckOp,
				Detail: field + " reason sets differ: " + sa[i] + " != " + sb[i]}
		}
	}
	return nil
}

func crossCheckPRS(a, b *PRSResult) error {
	if err := numwrap.CrossCheckInt(crossCheckOp, "prs.nsamples", len(a.Scores), len(b.Scores)); err != nil {
		return err
	}
	for i := range a.Scores {
		if err := numwrap.CrossCheckScalar(crossCheckOp, "prs.score", a.Scores[i], b.Scores[i]); err != nil {
			return err
		}
		if a.NActive[i] != b.NActive[i] {
			return &numwrap.BackendMismatchError{Op: crossCheckOp,
				Detail: "prs discrete active-site count differs"}
		}
	}
	if len(a.Audit) != len(b.Audit) {
		return &numwrap.BackendMismatchError{Op: crossCheckOp, Detail: "prs audit length differs"}
	}
	for i := range a.Audit {
		x, y := a.Audit[i], b.Audit[i]
		if x.ID != y.ID || x.Used != y.Used || x.Reason != y.Reason ||
			x.MissingGT != y.MissingGT || x.Sign != y.Sign {
			return &numwrap.BackendMismatchError{Op: crossCheckOp,
				Detail: "prs discrete audit differs for " + x.ID}
		}
	}
	return nil
}
