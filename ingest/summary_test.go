package ingest

import (
	"strings"
	"testing"

	"github.com/lechandonga/statistical-genetics-go/variant"
)

func optsBeta() SummaryOptions {
	return SummaryOptions{
		Source: "sumstat.tsv", Effect: EffectBeta,
		FallbackN: 1000, Limits: DefaultLimits(),
	}
}

func TestParseSummaryHappyAndAliases(t *testing.T) {
	in := "SNP\tCHR\tBP\tA1\tA2\tBETA\tSE\tP\tN\tFRQ\n" +
		"rs1\t1\t100\tA\tG\t0.2\t0.05\t0.001\t500\t0.3\n" +
		"rs2\tchr2\t200\tt\tc\t-0.1\t0.04\t0.02\t500\t0.7\n"
	tab, err := ParseSummary(strings.NewReader(in), optsBeta())
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("[input-summary] %s", tab.Report.Summary)
	if len(tab.Records) != 2 {
		t.Fatalf("n=%d", len(tab.Records))
	}
	r1 := tab.Records[0]
	if r1.V.Key() != "RS1" || string(r1.V.EA) != "A" {
		t.Fatalf("r1 %+v", r1.V)
	}
	if r1.Z != 0.2/0.05 {
		t.Fatalf("z=%v", r1.Z)
	}
	if r1.N != 500 {
		t.Fatalf("n=%v", r1.N)
	}
	// 第二行小写等位基因归一化
	if string(tab.Records[1].V.EA) != "T" || tab.Records[1].V.Chrom != "2" {
		t.Fatalf("r2 not normalized: %+v", tab.Records[1].V)
	}
	t.Logf("[judgement] accepted 2 rows; aliases SNP/A1/A2/BETA/FRQ resolved; z=beta/se")
}

func TestParseSummaryRejections(t *testing.T) {
	in := "SNP\tA1\tA2\tBETA\tSE\n" +
		"rs1\tA\tG\t0.2\t0.05\n" + // ok（FallbackN）
		"rs2\tA\tA\t0.2\t0.05\n" + // 相同等位基因
		"rs3\tA\tG\tNA\t0.05\n" + // beta 非法
		"rs4\tA\tG\t0.2\t-0.05\n" + // 负 SE
		"rs5\tI\tG\t0.2\t0.05\n" + // 非法等位基因
		"rs1\tG\tA\t0.9\t0.05\n" + // 重复
		"\tA\tG\t0.2\t0.05\n" + // 无身份
		"rs6\tA\tG\t0.2\t0.05\n" // ok
	tab, err := ParseSummary(strings.NewReader(in), optsBeta())
	if err == nil {
		t.Fatalf("strict mode must reject, got %d records", len(tab.Records))
	}
	sumErr, ok := err.(*RejectSummaryError)
	if !ok {
		t.Fatal(err)
	}
	rep := sumErr.Report
	t.Logf("[input-summary] %s", rep.Summary)
	if rep.Accepted != 2 || len(rep.Rejected) != 6 {
		t.Fatalf("accepted=%d rejected=%d", rep.Accepted, len(rep.Rejected))
	}
	got := map[variant.Reason]bool{}
	for _, rr := range rep.Rejected {
		t.Logf("[judgement] REJECT %s", rr.String())
		got[rr.Reason] = true
	}
	for _, want := range []variant.Reason{
		variant.ReasonInvalidAllele, variant.ReasonInvalidNumeric,
		variant.ReasonDuplicateVariant, variant.ReasonMissingID} {
		if !got[want] {
			t.Fatalf("missing expected reason %s", want)
		}
	}
}

func TestParseSummaryLenientKeepsAccepted(t *testing.T) {
	in := "SNP\tA1\tA2\tBETA\tSE\n" +
		"rs1\tA\tG\t0.2\t0.05\n" +
		"rs2\tA\tA\t0.2\t0.05\n"
	op := optsBeta()
	op.Lenient = true
	tab, err := ParseSummary(strings.NewReader(in), op)
	if err != nil {
		t.Fatal(err)
	}
	if len(tab.Records) != 1 || len(tab.Report.Rejected) != 1 {
		t.Fatalf("records=%d rejected=%d", len(tab.Records), len(tab.Report.Rejected))
	}
	t.Logf("[input-summary] lenient: %s", tab.Report.Summary)
}

func TestParseSummaryMissingColumns(t *testing.T) {
	in := "SNP\tBETA\tSE\nrs1\t0.2\t0.05\n"
	_, err := ParseSummary(strings.NewReader(in), optsBeta())
	r, ok := variant.AsReject(err)
	if !ok || r.Reason != variant.ReasonMissingField {
		t.Fatalf("want missing_field reject, got %v", err)
	}
	t.Logf("[input] header without allele columns [judgement] REJECT %s: %s", r.Reason, r.Detail)
}

func TestParseSummaryORandZ(t *testing.T) {
	in := "SNP\tA1\tA2\tOR\tSE\tN\nrs1\tA\tG\t1.22\t0.1\t100\n"
	op := optsBeta()
	op.Effect = EffectOR
	op.FallbackN = 0
	tab, err := ParseSummary(strings.NewReader(in), op)
	if err != nil {
		t.Fatal(err)
	}
	if tab.Records[0].Beta <= 0 {
		t.Fatal("OR must be converted to log(beta)>0")
	}
	t.Logf("[input] OR=1.22 [judgement] beta=log(OR)=%v", tab.Records[0].Beta)

	zin := "SNP\tA1\tA2\tZ\nrs1\tA\tG\t3.1\n"
	zop := optsBeta()
	zop.Effect = EffectZ
	ztab, err := ParseSummary(strings.NewReader(zin), zop)
	if err != nil {
		t.Fatal(err)
	}
	if !ztab.Records[0].HasZ || ztab.Records[0].Z != 3.1 {
		t.Fatal("z mode failed")
	}
	t.Log("[input] Z-only row [judgement] accepted in z-mode, SE left NaN")
}

func TestParseSummaryMissingNRejected(t *testing.T) {
	in := "SNP\tA1\tA2\tBETA\tSE\nrs1\tA\tG\t0.2\t0.05\n"
	op := optsBeta()
	op.FallbackN = 0
	_, err := ParseSummary(strings.NewReader(in), op)
	if err == nil {
		t.Fatal("missing N without fallback must reject")
	}
	t.Logf("[input] no N column and FallbackN=0 [judgement] REJECT %v", err)
}

func TestParseSummaryLimitTruncate(t *testing.T) {
	in := "SNP\tA1\tA2\tBETA\tSE\tN\n"
	for i := 0; i < 10; i++ {
		in += "rs" + string(rune('0'+i)) + "\tA\tG\t0.1\t0.1\t10\n"
	}
	op := optsBeta()
	op.FallbackN = 0
	op.Limits.MaxVariants = 3
	op.Limits.Policy = LimitTruncate
	op.Lenient = true
	tab, err := ParseSummary(strings.NewReader(in), op)
	if err != nil {
		t.Fatal(err)
	}
	if len(tab.Records) != 3 || !tab.Report.Truncated {
		t.Fatalf("records=%d truncated=%v", len(tab.Records), tab.Report.Truncated)
	}
	t.Logf("[input-summary] truncate: %s", tab.Report.Summary)

	op.Limits.Policy = LimitReject
	if _, err := ParseSummary(strings.NewReader(in), op); err == nil {
		t.Fatal("reject policy must error on overflow")
	} else {
		t.Logf("[input] 10 rows limit=3 policy=reject [judgement] %v", err)
	}
}
