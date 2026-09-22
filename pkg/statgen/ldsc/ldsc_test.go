package ldsc

import (
	"bytes"
	"math"
	"strconv"
	"testing"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/ldscore"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/sumstats"
)

func TestParseFixturesAndEstimate(t *testing.T) {
	c, err := sumstats.ParseFile("../../../testdata/cohort_A.tsv",
		sumstats.ParseConfig{CohortName: "A", HasHeader: true, SkipReject: true})
	if err != nil {
		t.Fatal(err)
	}
	tab, err := ldscore.ParseFile("../../../testdata/ldscores.tsv",
		ldscore.ParseConfig{Source: "ld", HasHeader: true, SkipReject: true})
	if err != nil {
		t.Fatal(err)
	}
	res, err := Estimate(c, tab, Config{Name: "A"})
	if err != nil {
		t.Fatalf("estimate: %v", err)
	}
	t.Logf("INPUT %s", res.Report.Summary())
	t.Logf("JUDGE h2=%.4f se=%.4f ci=[%.4f,%.4f] intercept=%.3f status=%s n=%d r2=%.3f",
		res.H2, res.SE, res.CI[0], res.CI[1], res.Intercept, res.Status, res.NVariants, res.R2)
	if res.NVariants != 7 {
		t.Fatalf("matched=%d want 7 (3 ambiguous rejected)", res.NVariants)
	}
	if !(res.SE > 0) {
		t.Fatalf("SE must be positive: %v", res.SE)
	}
}

func TestSerialEqualsParallel(t *testing.T) {
	c, _ := sumstats.ParseFile("../../../testdata/cohort_A.tsv",
		sumstats.ParseConfig{CohortName: "A", HasHeader: true, SkipReject: true})
	tab, _ := ldscore.ParseFile("../../../testdata/ldscores.tsv",
		ldscore.ParseConfig{Source: "ld", HasHeader: true, SkipReject: true})
	r1, err := Estimate(c, tab, Config{Workers: 1})
	if err != nil {
		t.Fatal(err)
	}
	r4, err := Estimate(c, tab, Config{Workers: 4})
	if err != nil {
		t.Fatal(err)
	}
	if math.Float64bits(r1.H2) != math.Float64bits(r4.H2) {
		t.Fatalf("h2 serial=%v parallel=%v", r1.H2, r4.H2)
	}
	if math.Float64bits(r1.SE) != math.Float64bits(r4.SE) {
		t.Fatalf("SE serial=%v parallel=%v", r1.SE, r4.SE)
	}
	if r1.NVariants != r4.NVariants {
		t.Fatal("variant count drift")
	}
}

func TestNegativeHeritabilityReported(t *testing.T) {
	// Construct chi2 values all below 1 => negative slope => h2 < 0.
	var buf bytes.Buffer
	buf.WriteString("CHR\tPOS\tREF\tALT\tBETA\tSE\tN\n")
	l2 := []float64{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}
	for i, l := range l2 {
		// chi2 declines with l2 (0.9 -> 0.0): clearly negative slope.
		chi2 := 1.0 - 0.1*l
		beta := math.Sqrt(chi2)
		_ = i
		buf.WriteString("1\t")
		buf.WriteString(itoa(100 + i))
		buf.WriteString("\tA\tG\t")
		buf.WriteString(strconv.FormatFloat(beta, 'f', 4, 64))
		buf.WriteString("\t1\t10000\n")
	}
	var lbuf bytes.Buffer
	lbuf.WriteString("CHR\tPOS\tREF\tALT\tL2\n")
	for i, l := range l2 {
		lbuf.WriteString("1\t" + itoa(100+i) + "\tA\tG\t" + strconv.FormatFloat(l, 'f', 1, 64) + "\n")
	}
	c, err := sumstats.Parse(&buf, sumstats.ParseConfig{
		CohortName: "neg", HasHeader: true, SkipReject: true})
	if err != nil {
		t.Fatal(err)
	}
	tab, err := ldscore.Parse(&lbuf, ldscore.ParseConfig{
		Source: "ld", HasHeader: true, SkipReject: true})
	if err != nil {
		t.Fatal(err)
	}
	res, err := Estimate(c, tab, Config{Name: "neg", NegPolicy: NegativeReport})
	if err != nil {
		t.Fatalf("negative h2 must be a Result, not an error: %v", err)
	}
	t.Logf("JUDGE negative-case h2=%.4f se=%.4f status=%s", res.H2, res.SE, res.Status)
	if res.Status != StatusNegative {
		t.Fatalf("status=%s want negative", res.Status)
	}
	if res.H2 >= 0 {
		t.Fatalf("h2=%v should be negative", res.H2)
	}
	// Same input again -> bit-identical.
	res2, _ := Estimate(c, tab, Config{Name: "neg", NegPolicy: NegativeReport})
	if math.Float64bits(res.H2) != math.Float64bits(res2.H2) {
		t.Fatal("non-reproducible estimate")
	}

	// Truncation policy clamps to zero but preserves status.
	res3, _ := Estimate(c, tab, Config{Name: "neg", NegPolicy: NegativeZero})
	if res3.H2 != 0 || res3.Status != StatusNegative {
		t.Fatalf("truncate: h2=%v status=%s", res3.H2, res3.Status)
	}
}

func TestTooFewVariantsUnidentified(t *testing.T) {
	in := "CHR\tPOS\tREF\tALT\tBETA\tSE\tN\n1\t10\tA\tG\t0.1\t0.05\t1000\n"
	ldin := "CHR\tPOS\tREF\tALT\tL2\n1\t10\tA\tG\t1.0\n"
	c, _ := sumstats.Parse(bytes.NewBufferString(in),
		sumstats.ParseConfig{CohortName: "x", HasHeader: true})
	tab, _ := ldscore.Parse(bytes.NewBufferString(ldin),
		ldscore.ParseConfig{Source: "ld", HasHeader: true})
	_, err := Estimate(c, tab, Config{})
	if !stgerr.Is(err, stgerr.KindOutOfBounds) {
		t.Fatalf("expected out-of-bounds/unidentified, got %v", err)
	}
	t.Logf("JUDGE degenerate: %v", err)
}

func itoa(i int) string {
	if i == 0 {
		return "0"
	}
	b := []byte{}
	for i > 0 {
		b = append([]byte{byte('0' + i%10)}, b...)
		i /= 10
	}
	return string(b)
}
