package ingest

import (
	"strings"
	"testing"
)

func TestUsageReportPopulated(t *testing.T) {
	in := "SNP\tA1\tA2\tBETA\tSE\tN\n" +
		"rs1\tA\tG\t0.1\t0.1\t10\n" +
		"rs2\tA\tG\t0.2\t0.1\t10\n"
	op := optsBeta()
	op.FallbackN = 0
	tab, err := ParseSummary(strings.NewReader(in), op)
	if err != nil {
		t.Fatal(err)
	}
	u := tab.Report.Usage
	if u.BytesRead <= 0 || u.Lines != 2 || u.Variants != 2 {
		t.Fatalf("usage wrong: %+v", u)
	}
	if u.Elapsed <= 0 {
		t.Fatal("elapsed should be recorded")
	}
	if u.LimitKind != "" {
		t.Fatalf("no limit should be hit, got %s", u.LimitKind)
	}
	t.Logf("[input] 2-row summary [judgement] %s", u)
}

func TestMaxFileBytesEnforcedRejectNoPartial(t *testing.T) {
	in := "SNP\tA1\tA2\tBETA\tSE\tN\n"
	for i := 0; i < 50; i++ {
		in += "rs00000000000000000000000000000000000000000000\tA\tG\t0.1\t0.1\t10\n"
	}
	op := optsBeta()
	op.FallbackN = 0
	op.Limits.MaxFileBytes = 120 // 小于完整文件
	tab, err := ParseSummary(strings.NewReader(in), op)
	if err == nil {
		t.Fatalf("tiny MaxFileBytes must reject, got %d records", len(tab.Records))
	}
	if !strings.Contains(err.Error(), "MaxFileBytes") {
		t.Fatalf("error must name MaxFileBytes: %v", err)
	}
	if tab != nil {
		t.Fatal("reject policy must not return partial results")
	}
	t.Logf("[input] 50 rows MaxFileBytes=120 [judgement] REJECT %v", err)

	// 极小行上限也必须真实生效
	op.Limits = DefaultLimits()
	op.FallbackN = 0
	op.Limits.MaxLineBytes = 10
	if _, err := ParseSummary(strings.NewReader(in), op); err == nil {
		t.Fatal("MaxLineBytes=10 must reject")
	}

	// MaxElapsed 已过 => 立即拒绝
	op.Limits = DefaultLimits()
	op.FallbackN = 0
	op.Limits.MaxElapsed = 1
	if _, err := ParseSummary(strings.NewReader(in), op); err == nil {
		t.Fatal("expired MaxElapsed must reject")
	}
}

func TestMaxVariantsTruncateReportsUsage(t *testing.T) {
	in := "SNP\tA1\tA2\tBETA\tSE\tN\n"
	for i := 0; i < 10; i++ {
		in += "rs" + string(rune('0'+i)) + "\tA\tG\t0.1\t0.1\t10\n"
	}
	op := optsBeta()
	op.FallbackN = 0
	op.Limits.MaxVariants = 4
	op.Limits.Policy = LimitTruncate
	tab, err := ParseSummary(strings.NewReader(in), op)
	if err != nil {
		t.Fatal(err)
	}
	if len(tab.Records) != 4 || !tab.Report.Usage.Truncated ||
		tab.Report.Usage.LimitKind != "MaxVariants" || tab.Report.Usage.Variants != 4 {
		t.Fatalf("truncate usage wrong: %+v", tab.Report.Usage)
	}
	t.Logf("[input] 10 rows cap=4 truncate [judgement] %s", tab.Report.Usage)
}

func TestGenotypeSampleTruncateAndReject(t *testing.T) {
	in := "SNP\tA1\tA2\ts1\ts2\ts3\ts4\n" +
		"rs1\tA\tG\t0\t1\t2\t1\n" +
		"rs2\tC\tG\t2\t1\t0\t1\n"

	// reject：超样本上限整体失败，无部分结果
	lim := DefaultLimits()
	lim.MaxSamples = 2
	if _, err := ParseGenotype(strings.NewReader(in),
		GenotypeOptions{Source: "gt", Limits: lim}); err == nil {
		t.Fatal("MaxSamples reject policy must fail")
	}

	// truncate：确定性保留前两个样本列
	lim.Policy = LimitTruncate
	m, err := ParseGenotype(strings.NewReader(in),
		GenotypeOptions{Source: "gt", Limits: lim})
	if err != nil {
		t.Fatal(err)
	}
	if m.NSamples() != 2 || !m.Report.Usage.Truncated ||
		m.Report.Usage.Samples != 2 || m.Report.Usage.LimitKind != "MaxSamples" {
		t.Fatalf("sample truncate wrong: %+v samples=%d", m.Report.Usage, m.NSamples())
	}
	if got := m.SampleIDs[0]; got != "s1" || m.SampleIDs[1] != "s2" {
		t.Fatalf("must keep first two sample columns, got %v", m.SampleIDs)
	}
	if d, _ := m.At(0, 1); d != 1 {
		t.Fatalf("dosage after truncation wrong: %v", d)
	}
	t.Logf("[input] 4 samples cap=2 truncate [judgement] kept %v; %s",
		m.SampleIDs, m.Report.Usage)
}

func TestUsageIsolatedAcrossConcurrentParses(t *testing.T) {
	// 并发解析各自不同上限，用量小结绝不能互相串。
	errs := make(chan error, 8)
	for g := 0; g < 8; g++ {
		go func(g int) {
			rows := 3 + g
			var sb strings.Builder
			sb.WriteString("SNP\tA1\tA2\tBETA\tSE\tN\n")
			for i := 0; i < rows; i++ {
				sb.WriteString("rs")
				sb.WriteString(string(rune('a' + g)))
				sb.WriteString(string(rune('0' + i)))
				sb.WriteString("\tA\tG\t0.1\t0.1\t10\n")
			}
			op := optsBeta()
			op.FallbackN = 0
			tab, err := ParseSummary(strings.NewReader(sb.String()), op)
			if err != nil {
				errs <- err
				return
			}
			if tab.Report.Usage.Variants != rows || tab.Report.Usage.Lines != rows {
				errs <- &usageLeak{got: tab.Report.Usage.Variants, want: rows, g: g}
				return
			}
			errs <- nil
		}(g)
	}
	for g := 0; g < 8; g++ {
		if err := <-errs; err != nil {
			t.Fatal(err)
		}
	}
	t.Log("[judgement] concurrent parses report independent usage")
}

type usageLeak struct {
	got, want, g int
}

func (e *usageLeak) Error() string {
	return "usage leaked between concurrent parses"
}
