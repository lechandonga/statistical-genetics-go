package ingest

import (
	"strings"
	"testing"

	"github.com/lechandonga/statistical-genetics-go/variant"
)

func TestParseLDScore(t *testing.T) {
	in := strings.Join([]string{
		"SNP\tCHR\tBP\tA1\tA2\tL2\tMAF",
		"rs1\t1\t10\tA\tG\t1.5\t0.25",
		"rs2\t1\t20\tC\tT\t2.0\t0.10",
		"rs3\t1\t30\tA\tT\t1.0\t0.30", // 模糊位点在解析层允许（对齐层处理）
		"rsdup\t1\t40\tA\tG\t-1\t0.1", // 负 L2
	}, "\n") + "\n"
	_, err := ParseLDScores(strings.NewReader(in), LDSOptions{Source: "ld", Limits: DefaultLimits()})
	if err == nil {
		t.Fatal("negative L2 must reject whole file in strict mode")
	}
	se, ok := err.(*RejectSummaryError)
	if !ok {
		t.Fatal(err)
	}
	t.Logf("[input-summary] %s", se.Report.Summary)
	if se.FirstReject().Reason != variant.ReasonInvalidNumeric {
		t.Fatalf("reason=%s", se.FirstReject().Reason)
	}
	t.Logf("[judgement] ambiguous A/T accepted at parse layer; negative L2 rejected")
}

func TestParseGenotypeDosage(t *testing.T) {
	in := "SNP\tA1\tA2\ts1\ts2\ts3\n" +
		"rs1\tA\tG\t0\t1\t2\n" +
		"rs2\tC\tT\tNA\t1\t0\n" +
		"rs3\tA\tG\t3\t1\t0\n" // 剂量越界
	_, err := ParseGenotype(strings.NewReader(in), GenotypeOptions{
		Source: "gt", Encoding: GTDosage012, Limits: DefaultLimits()})
	if err == nil {
		t.Fatal("dosage=3 must reject in strict mode")
	}
	se := err.(*RejectSummaryError)
	t.Logf("[input-summary] %s", se.Report.Summary)
	if se.FirstReject().Reason != variant.ReasonInvalidNumeric {
		t.Fatalf("reason=%s", se.FirstReject().Reason)
	}

	m2, err := ParseGenotype(strings.NewReader(in), GenotypeOptions{
		Source: "gt", Encoding: GTDosage012, Limits: DefaultLimits(), Lenient: true})
	if err != nil {
		t.Fatal(err)
	}
	if m2.NVariants() != 2 || m2.NSamples() != 3 {
		t.Fatalf("dim=%d x %d", m2.NVariants(), m2.NSamples())
	}
	if d, miss := m2.At(1, 0); !miss || d == d {
		t.Fatalf("s1/rs2 should be missing NaN, got %v miss=%v", d, miss)
	}
	if d, miss := m2.At(0, 2); miss || d != 2 {
		t.Fatalf("s3/rs1 should be dosage 2, got %v", d)
	}
	t.Logf("[input-summary] lenient: %s", m2.Report.Summary)
}

func TestParseGenotypeACGT(t *testing.T) {
	in := "SNP\tA1\tA2\ts1\ts2\ts3\nrs1\tA\tG\tAA\tAG\tGG\n"
	m, err := ParseGenotype(strings.NewReader(in), GenotypeOptions{
		Source: "gt", Encoding: GTACGT, Limits: DefaultLimits()})
	if err != nil {
		t.Fatal(err)
	}
	want := []float64{2, 1, 0}
	for j := 0; j < 3; j++ {
		if d, _ := m.At(0, j); d != want[j] {
			t.Fatalf("j=%d dose=%v want %v", j, d, want[j])
		}
	}
	t.Logf("[input] AA/AG/GG with EA=A [judgement] dosages 2/1/0")

	bad := "SNP\tA1\tA2\ts1\nrs1\tA\tG\tAT\n"
	if _, err := ParseGenotype(strings.NewReader(bad), GenotypeOptions{
		Source: "gt", Encoding: GTACGT, Limits: DefaultLimits()}); err == nil {
		t.Fatal("AT against A/G site must reject")
	} else {
		t.Logf("[input] AT at A/G site [judgement] REJECT %v", err)
	}
}

func TestParseGenotypeDuplicateAndSampleLimit(t *testing.T) {
	in := "SNP\tA1\tA2\ts1\ts2\nrs1\tA\tG\t0\t1\nrs1\tG\tA\t1\t0\n"
	_, err := ParseGenotype(strings.NewReader(in), GenotypeOptions{
		Source: "gt", Encoding: GTDosage012, Limits: DefaultLimits()})
	if err == nil {
		t.Fatal("duplicate variant must reject")
	}
	t.Logf("[input] duplicate rs1 rows [judgement] REJECT %v", err)

	lim := DefaultLimits()
	lim.MaxSamples = 1
	_, err = ParseGenotype(strings.NewReader(in), GenotypeOptions{
		Source: "gt", Encoding: GTDosage012, Limits: lim})
	if err == nil {
		t.Fatal("sample cap must trigger reject")
	}
	t.Logf("[input] 2 samples cap=1 [judgement] REJECT %v", err)
}
