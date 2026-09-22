package sumstats

import (
	"strings"
	"testing"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
)

func TestParseBadRowsDistinguishable(t *testing.T) {
	c, err := ParseFile("../../../testdata/cohort_bad.tsv", ParseConfig{
		CohortName: "bad", HasHeader: true, SkipReject: true})
	if err != nil {
		t.Fatalf("tolerant parse: %v", err)
	}
	t.Logf("INPUT %s", c.Report.Summary())
	got := map[stgerr.Kind]int{}
	for _, rj := range c.Report.Rejected {
		got[rj.Kind]++
		t.Logf("  reject line=%d kind=%s reason=%q", rj.Line, rj.Kind, rj.Reason)
	}
	want := map[stgerr.Kind]int{
		stgerr.KindInvalidValue:     7, // NaN, bad chrom, identical alleles, P>1, SE=0, N<0, SE<0
		stgerr.KindUnknownAllele:    1, // A/X
		stgerr.KindDuplicateVariant: 1, // duplicate 1:100
	}
	for k, n := range want {
		if got[k] != n {
			t.Errorf("kind %s: got %d want %d", k, got[k], n)
		}
	}
	// bad chromosome line (23) is also KindInvalidValue — verify total.
	totalBad := 0
	for _, n := range got {
		totalBad += n
	}
	if totalBad != 9 {
		t.Errorf("expected 9 rejected rows, got %d", totalBad)
	}
}

func TestParseStrictAborts(t *testing.T) {
	_, err := ParseFile("../../../testdata/cohort_bad.tsv", ParseConfig{
		CohortName: "bad", HasHeader: true, SkipReject: false})
	if err == nil {
		t.Fatal("strict mode must abort on first rejection")
	}
	t.Logf("strict abort reason: %v", err)
}

func TestParseGoodCohort(t *testing.T) {
	c, err := ParseFile("../../../testdata/cohort_A.tsv", ParseConfig{
		CohortName: "A", HasHeader: true, SkipReject: true})
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("INPUT %s", c.Report.Summary())
	if c.Report.Accepted != 10 {
		t.Fatalf("accepted=%d want 10", c.Report.Accepted)
	}
}

func TestMissingRequiredColumn(t *testing.T) {
	_, err := Parse(strings.NewReader("CHR\tPOS\n1\t10\n"), ParseConfig{
		CohortName: "m", HasHeader: true})
	// First data row must surface missing ref/effect/beta/se fields.
	if !stgerr.Is(err, stgerr.KindMissingField) {
		t.Fatalf("expected missing-field error, got %v", err)
	}
}
