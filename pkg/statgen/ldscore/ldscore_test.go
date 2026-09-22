package ldscore

import (
	"strings"
	"testing"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
)

func TestParseLDGood(t *testing.T) {
	tab, err := ParseFile("../../../testdata/ldscores.tsv", ParseConfig{
		Source: "ld", HasHeader: true, SkipReject: true})
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("INPUT %s", tab.Report.Summary())
	if tab.Report.Accepted != 10 {
		t.Fatalf("accepted=%d want 10", tab.Report.Accepted)
	}
}

func TestParseLDRejections(t *testing.T) {
	in := "CHR\tPOS\tREF\tALT\tL2\n1\t10\tA\tG\t1.1\n1\t11\tA\tG\t-0.2\n1\t12\tA\tN\t1\n1\t10\tA\tG\t2.2\n"
	tab, err := Parse(strings.NewReader(in), ParseConfig{
		Source: "ld", HasHeader: true, SkipReject: true})
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("INPUT %s", tab.Report.Summary())
	kinds := map[stgerr.Kind]int{}
	for _, rj := range tab.Report.Rejected {
		kinds[rj.Kind]++
		t.Logf("  reject line=%d kind=%s reason=%q", rj.Line, rj.Kind, rj.Reason)
	}
	if kinds[stgerr.KindInvalidValue] != 1 {
		t.Errorf("negative L2: %d", kinds[stgerr.KindInvalidValue])
	}
	if kinds[stgerr.KindUnknownAllele] != 1 {
		t.Errorf("unknown allele: %d", kinds[stgerr.KindUnknownAllele])
	}
	if kinds[stgerr.KindDuplicateVariant] != 1 {
		t.Errorf("duplicate: %d", kinds[stgerr.KindDuplicateVariant])
	}
	if tab.Report.Accepted != 1 {
		t.Errorf("accepted=%d want 1", tab.Report.Accepted)
	}
}

func TestParseLDLimit(t *testing.T) {
	in := "CHR\tPOS\tREF\tALT\tL2\n1\t10\tA\tG\t1.1\n1\t11\tC\tT\t1.2\n"
	_, err := Parse(strings.NewReader(in), ParseConfig{
		Source: "ld", HasHeader: true, MaxRows: 1})
	if !stgerr.Is(err, stgerr.KindLimitExceeded) {
		t.Fatalf("limit: %v", err)
	}
	t.Logf("limit reason: %v", err)
}

func TestParseLDTruncationMarked(t *testing.T) {
	in := "CHR\tPOS\tREF\tALT\tL2\n" +
		"1\t10\tA\tG\t1\n1\t11\tC\tT\t1\n1\t12\tG\tA\t1\n"
	tab, err := Parse(strings.NewReader(in), ParseConfig{
		Source: "ld", HasHeader: true, MaxRows: 2, SkipReject: true})
	if err != nil {
		t.Fatal(err)
	}
	if !tab.Truncated || !tab.Report.Truncated {
		t.Fatal("truncation must be explicitly marked")
	}
	if tab.Report.Accepted != 2 {
		t.Fatalf("accepted=%d want 2", tab.Report.Accepted)
	}
	t.Logf("JUDGE %s", tab.Report.Summary())
}
