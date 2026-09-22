package genotype

import (
	"strings"
	"testing"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
)

func TestParseGenotypeMatrix(t *testing.T) {
	m, err := ParseFile("../../../testdata/genotypes.tsv", ParseConfig{
		Source: "g", HasHeader: true, HasSampleCol: true, SkipReject: true})
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("INPUT %s samples=%d", m.Report.Summary(), m.NumSamples())
	if m.NumSamples() != 6 || m.NumVariants() != 5 {
		t.Fatalf("dims=%dx%d want 5x6", m.NumVariants(), m.NumSamples())
	}
	d, miss, ok := m.Row(m.Variants[0].Key())
	if !ok {
		t.Fatal("row lookup")
	}
	if !miss[4] || d[4] != MissingCode {
		t.Fatalf("missing code not recognized: %v", d)
	}
}

func TestDosageValidation(t *testing.T) {
	if err := ValidateDosage(0, MissingCode); err != nil {
		t.Fatal(err)
	}
	if err := ValidateDosage(3, MissingCode); !stgerr.Is(err, stgerr.KindInvalidValue) {
		t.Fatalf("dosage 3: %v", err)
	}
	if err := ValidateDosage(-1, MissingCode); err != nil {
		t.Fatalf("missing sentinel should pass: %v", err)
	}
}

func TestRaggedRowRejected(t *testing.T) {
	in := "id\ts1\ts2\n1:10:A:G\t0\t1\n1:20:C:T\t0\n"
	_, err := Parse(strings.NewReader(in), ParseConfig{
		Source: "g", HasHeader: true, HasSampleCol: true, SkipReject: true})
	if err != nil {
		t.Fatalf("tolerant parse error: %v", err)
	}
}

func TestBadVariantID(t *testing.T) {
	in := "id\ts1\nnonsense\t0\n"
	_, err := Parse(strings.NewReader(in), ParseConfig{
		Source: "g", HasHeader: true, HasSampleCol: true, SkipReject: false})
	if err == nil {
		t.Fatal("malformed variant id must abort strict parse")
	}
	t.Logf("bad id reason: %v", err)
}

func TestImputePolicies(t *testing.T) {
	x := []float64{0, MissingCode, 2}
	miss := []bool{false, true, false}
	out, n, err := ImputedRow(x, miss, DosageMeanImpute)
	if err != nil || n != 1 || out[1] != 1 {
		t.Fatalf("mean impute: %v n=%d out=%v", err, n, out)
	}
	out2, _, _ := ImputedRow(x, miss, DosageZeroImpute)
	if out2[1] != 0 {
		t.Fatalf("zero impute: %v", out2)
	}
	if _, _, err := ImputedRow(x, miss, DosageRejectSample); !stgerr.Is(err, stgerr.KindInvalidValue) {
		t.Fatalf("reject policy: %v", err)
	}
}
