package ingest

import (
	"fmt"
	"strings"

	"github.com/lechandonga/statistical-genetics-go/variant"
)

// Err 返回报告对应的聚合错误：存在被拒绝行时返回 *RejectSummaryError。
func (r *ParseReport) Err() error {
	if len(r.Rejected) == 0 {
		return nil
	}
	return &RejectSummaryError{Source: r.Source, Report: r}
}

// RejectSummaryError 聚合一次解析的全部拒绝原因，可与 *variant.RejectError 区分。
type RejectSummaryError struct {
	Source string
	Report *ParseReport
}

func (e *RejectSummaryError) Error() string {
	var b strings.Builder
	fmt.Fprintf(&b, "source=%s: %d row(s) rejected of %d (accepted=%d)",
		e.Source, len(e.Report.Rejected), e.Report.TotalRows, e.Report.Accepted)
	for i := 0; i < len(e.Report.Rejected) && i < 10; i++ {
		b.WriteString("\n  - ")
		b.WriteString(e.Report.Rejected[i].String())
	}
	if len(e.Report.Rejected) > 10 {
		fmt.Fprintf(&b, "\n  ... and %d more", len(e.Report.Rejected)-10)
	}
	return b.String()
}

// FirstReject 返回第一条拒绝明细。
func (e *RejectSummaryError) FirstReject() RowRef { return e.Report.Rejected[0] }

// summarize 生成输入摘要文本（供测试日志与 CLI 使用）。
func summarize(r ParseReport) string {
	reasons := map[variant.Reason]int{}
	for _, x := range r.Rejected {
		reasons[x.Reason]++
	}
	parts := make([]string, 0, len(reasons))
	for reason, n := range reasons {
		parts = append(parts, fmt.Sprintf("%s=%d", reason, n))
	}
	return fmt.Sprintf("rows=%d accepted=%d rejected=%d reasons={%s} truncated=%v",
		r.TotalRows, r.Accepted, len(r.Rejected), strings.Join(parts, ","), r.Truncated)
}
