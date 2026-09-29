package ingest

import (
	"strings"
	"testing"
	"time"
)

const sumInput = "SNP\tA1\tA2\tBETA\tSE\tN\tEAF\n" +
	"rs1\tA\tG\t0.1\t0.2\t100\t0.3\n" +
	"rs2\tC\tT\t0.2\t0.2\t100\t0.3\n" +
	"rs3\tA\tG\t0.3\t0.2\t100\t0.3\n" +
	"rs4\tG\tA\t0.4\t0.2\t100\t0.3\n" +
	"rs5\tA\tC\t0.5\t0.2\t100\t0.3\n"

func TestUsageReportedForAllInputTypes(t *testing.T) {
	sum, err := ParseSummary(strings.NewReader(sumInput), SummaryOptions{
		Source: "sum", Limits: DefaultLimits()})
	if err != nil {
		t.Fatal(err)
	}
	u := sum.Report.Usage
	if u.BytesRead != int64(len(sumInput)) {
		t.Fatalf("bytes=%d want %d", u.BytesRead, len(sumInput))
	}
	if u.LinesRead != 6 || u.DataRows != 5 || u.Elapsed <= 0 {
		t.Fatalf("usage=%+v", u)
	}

	ldIn := "SNP\tA1\tA2\tL2\nrs1\tA\tG\t1.0\nrs2\tC\tT\t2.0\n"
	ld, err := ParseLDScores(strings.NewReader(ldIn), LDSOptions{
		Source: "ld", Limits: DefaultLimits()})
	if err != nil {
		t.Fatal(err)
	}
	if ld.Report.Usage.BytesRead != int64(len(ldIn)) ||
		ld.Report.Usage.DataRows != 2 {
		t.Fatalf("ld usage=%+v", ld.Report.Usage)
	}

	gtIn := "SNP\tA1\tA2\ts1\ts2\ts3\nrs1\tA\tG\t0\t1\t2\n"
	gt, err := ParseGenotype(strings.NewReader(gtIn), GenotypeOptions{
		Source: "gt", Limits: DefaultLimits()})
	if err != nil {
		t.Fatal(err)
	}
	if gt.Report.Usage.Samples != 3 || gt.Report.Usage.SamplesOffered != 3 ||
		gt.Report.Usage.BytesRead != int64(len(gtIn)) {
		t.Fatalf("gt usage=%+v", gt.Report.Usage)
	}
	t.Logf("[input-summary] sum: %s", u.Summary())
	t.Logf("[input-summary] ld : %s", ld.Report.Usage.Summary())
	t.Logf("[input-summary] gt : %s", gt.Report.Usage.Summary())
}

// TestLimitsSmallEnforced 用极小上限验证每类上限都真实生效。
func TestLimitsSmallEnforced(t *testing.T) {
	// MaxVariants=2 + reject：整体失败，不返回表
	_, err := ParseSummary(strings.NewReader(sumInput), SummaryOptions{
		Source: "sum",
		Limits: Limits{MaxVariants: 2, MaxLineBytes: 1 << 20,
			MaxFileBytes: 1 << 30, MaxElapsed: time.Minute, Policy: LimitReject}})
	if err == nil || !strings.Contains(err.Error(), "limit exceeded") {
		t.Fatalf("MaxVariants reject must fail, got %v", err)
	}

	// MaxVariants=2 + truncate：恰好 2 行，标记截断，无半结果
	tbl, err := ParseSummary(strings.NewReader(sumInput), SummaryOptions{
		Source: "sum",
		Limits: Limits{MaxVariants: 2, MaxLineBytes: 1 << 20,
			MaxFileBytes: 1 << 30, MaxElapsed: time.Minute, Policy: LimitTruncate}})
	if err != nil {
		t.Fatal(err)
	}
	if len(tbl.Records) != 2 || !tbl.Report.Truncated ||
		tbl.Report.TruncateAt != 2 {
		t.Fatalf("truncate: records=%d report=%+v", len(tbl.Records), tbl.Report)
	}

	// MaxFileBytes 极小 + reject
	_, err = ParseSummary(strings.NewReader(sumInput), SummaryOptions{
		Source: "sum",
		Limits: Limits{MaxVariants: 100, MaxLineBytes: 1 << 20,
			MaxFileBytes: 30, MaxElapsed: time.Minute, Policy: LimitReject}})
	if err == nil || !strings.Contains(err.Error(), "MaxFileBytes") {
		t.Fatalf("MaxFileBytes reject must fail, got %v", err)
	}

	// MaxLineBytes 极小：畸形巨行在任何策略下都拒绝
	_, err = ParseSummary(strings.NewReader(sumInput), SummaryOptions{
		Source: "sum",
		Limits: Limits{MaxVariants: 100, MaxLineBytes: 10,
			MaxFileBytes: 1 << 30, MaxElapsed: time.Minute, Policy: LimitTruncate}})
	if err == nil || !strings.Contains(err.Error(), "MaxLineBytes") {
		t.Fatalf("oversized line must always reject, got %v", err)
	}
}

// TestMaxSamplesPolicies 验证基因型样本上限：拒绝不留结果，截断保留前 N 列。
func TestMaxSamplesPolicies(t *testing.T) {
	gtIn := "SNP\tA1\tA2\ts1\ts2\ts3\ts4\nrs1\tA\tG\t0\t1\t2\t1\n"
	_, err := ParseGenotype(strings.NewReader(gtIn), GenotypeOptions{
		Source: "gt",
		Limits: Limits{MaxSamples: 2, MaxLineBytes: 1 << 20,
			MaxFileBytes: 1 << 30, Policy: LimitReject}})
	if err == nil || !strings.Contains(err.Error(), "MaxSamples") {
		t.Fatalf("MaxSamples reject must fail, got %v", err)
	}
	m, err := ParseGenotype(strings.NewReader(gtIn), GenotypeOptions{
		Source: "gt",
		Limits: Limits{MaxSamples: 2, MaxLineBytes: 1 << 20,
			MaxFileBytes: 1 << 30, Policy: LimitTruncate}})
	if err != nil {
		t.Fatal(err)
	}
	if m.NSamples() != 2 || m.SampleIDs[0] != "s1" || m.SampleIDs[1] != "s2" ||
		!m.Report.Truncated {
		t.Fatalf("sample truncate wrong: ids=%v report=%+v", m.SampleIDs, m.Report)
	}
	if m.Report.Usage.Samples != 2 || m.Report.Usage.SamplesOffered != 4 {
		t.Fatalf("usage samples=%d/%d", m.Report.Usage.Samples, m.Report.Usage.SamplesOffered)
	}
	if d, miss := m.At(0, 1); miss || d != 1 {
		t.Fatalf("kept dosage wrong: %v %v", d, miss)
	}
	t.Logf("[judgement] MaxSamples=2 truncate kept %v of %d offered",
		m.SampleIDs, m.Report.Usage.SamplesOffered)
}

// TestMaxElapsedTruncate 验证极小墙钟预算在截断策略下结束于完整行边界，
// 且返回的是完整记录（无半行）。
func TestMaxElapsedTruncate(t *testing.T) {
	var b strings.Builder
	b.WriteString("SNP\tA1\tA2\tBETA\tSE\tN\tEAF\n")
	for i := 0; i < 5000; i++ {
		b.WriteString("rs_0000000000\tA\tG\t0.1\t0.2\t100\t0.3\n")
	}
	tbl, err := ParseSummary(strings.NewReader(b.String()), SummaryOptions{
		Source: "sum",
		Limits: Limits{MaxVariants: 1 << 30, MaxLineBytes: 1 << 20,
			MaxFileBytes: 1 << 40, MaxElapsed: 1, Policy: LimitTruncate}})
	if err != nil {
		t.Fatal(err)
	}
	if !tbl.Report.Truncated {
		t.Fatalf("tiny elapsed budget should truncate, got %+v", tbl.Report.Usage)
	}
	t.Logf("[judgement] MaxElapsed=1ns truncated at %d complete rows",
		len(tbl.Records))
}

// TestUsagePerCallIsolated 并发解析多个任务，各自的 Usage 不得串扰。
func TestUsagePerCallIsolated(t *testing.T) {
	done := make(chan error, 3)
	go func() {
		t1, e := ParseSummary(strings.NewReader(sumInput), SummaryOptions{
			Source: "a", Limits: DefaultLimits()})
		if e == nil && t1.Report.Usage.BytesRead != int64(len(sumInput)) {
			e = &usageError{t1.Report.Usage.BytesRead}
		}
		done <- e
	}()
	go func() {
		var bld strings.Builder
		for i := 0; i < 50; i++ {
			bld.WriteString("rs")
			bld.WriteString(itoa(i))
			bld.WriteString("\tA\tG\t0.1\t0.2\t100\t0.3\n")
		}
		hdr := "SNP\tA1\tA2\tBETA\tSE\tN\tEAF\n"
		t2, e := ParseSummary(strings.NewReader(hdr+bld.String()), SummaryOptions{
			Source: "b", Limits: DefaultLimits()})
		if e == nil && t2.Report.Usage.DataRows != 50 {
			e = &usageError{int64(t2.Report.Usage.DataRows)}
		}
		done <- e
	}()
	go func() {
		var bld strings.Builder
		for i := 0; i < 10; i++ {
			bld.WriteString("rs")
			bld.WriteString(itoa(i))
			bld.WriteString("\tA\tG\t1.0\n")
		}
		t3, e := ParseLDScores(strings.NewReader("SNP\tA1\tA2\tL2\n"+bld.String()),
			LDSOptions{Source: "c", Limits: DefaultLimits()})
		if e == nil && t3.Report.Usage.DataRows != 10 {
			e = &usageError{int64(t3.Report.Usage.DataRows)}
		}
		done <- e
	}()
	for i := 0; i < 3; i++ {
		if err := <-done; err != nil {
			t.Fatal(err)
		}
	}
}

type usageError struct{ v int64 }

func (e *usageError) Error() string { return "usage cross-talk detected" }

func itoa(i int) string {
	if i == 0 {
		return "0"
	}
	var b [12]byte
	k := len(b)
	for i > 0 {
		k--
		b[k] = byte('0' + i%10)
		i /= 10
	}
	return string(b[k:])
}
