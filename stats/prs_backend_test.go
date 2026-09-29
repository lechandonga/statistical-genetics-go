package stats

import (
	"context"
	"math"
	"testing"

	"github.com/lechandonga/statistical-genetics-go/ingest"
	"github.com/lechandonga/statistical-genetics-go/numwrap"
	"github.com/lechandonga/statistical-genetics-go/variant"
)

func prsTables(t *testing.T) (*ingest.SummaryTable, *ingest.GenotypeMatrix) {
	t.Helper()
	nS, nV := 40, 30
	recs := make([]ingest.SummaryRecord, 0, nV)
	variants := make([]*variant.Variant, 0, nV)
	dosage := make([]float64, 0, nV*nS)
	missing := make([]bool, 0, nV*nS)
	for i := 0; i < nV; i++ {
		id := idFor(i)
		v, _ := variant.New(id, "1", int64(100+i), "A", "G", 0.3, variant.StrandUnknown)
		recs = append(recs, ingest.SummaryRecord{
			V: v, Beta: 0.02 * float64(1+i%3-1), SE: 0.01, N: 1000,
			Z: 2, P: 0.01, HasZ: true})
		variants = append(variants, v)
		for j := 0; j < nS; j++ {
			dosage = append(dosage, float64((i+j)%3))
			missing = append(missing, false)
		}
	}
	sampleIDs := make([]string, nS)
	for j := range sampleIDs {
		sampleIDs[j] = "s" + itoa(j)
	}
	return &ingest.SummaryTable{Records: recs},
		&ingest.GenotypeMatrix{Variants: variants, SampleIDs: sampleIDs,
			Dosage: dosage, Missing: missing}
}

func TestPRSBackendsAgree(t *testing.T) {
	model, gt := prsTables(t)
	base := PRSOptions{ScoreScale: "std", Model: ModelIgnore, Backend: numwrap.BackendStd}
	rs, err := BuildPRS(context.Background(), model, gt, base)
	if err != nil {
		t.Fatal(err)
	}
	rg, err := BuildPRS(context.Background(), model, gt, PRSOptions{
		ScoreScale: "std", Model: ModelIgnore, Backend: numwrap.BackendGonum})
	if err != nil {
		t.Fatal(err)
	}
	if err := comparePRSResults(rs, rg); err != nil {
		t.Fatal(err)
	}
	rc, err := BuildPRS(context.Background(), model, gt, PRSOptions{
		ScoreScale: "std", Model: ModelIgnore, Workers: 4, Backend: numwrap.BackendCrosscheck})
	if err != nil {
		t.Fatal(err)
	}
	if len(rc.Scores) != len(rs.Scores) {
		t.Fatal("crosscheck score count mismatch")
	}
	// std 尺度：均值必须为 0
	var sum float64
	for _, s := range rc.Scores {
		sum += s
	}
	if math.Abs(sum) > 1e-9 {
		t.Fatalf("standardized mean=%v", sum)
	}
	t.Logf("[judgement] PRS std scale: stdlib/gonum/crosscheck agree over %d samples",
		len(rc.Scores))

	// EvaluatePRS 也支持双后端
	pheno := make([]float64, len(rs.Scores))
	for i := range pheno {
		pheno[i] = 0.5*rs.Scores[i] + float64(i%5)*0.01
	}
	es, err := EvaluatePRS(context.Background(), rs, EvalOptions{Phenotype: pheno, Backend: numwrap.BackendStd})
	if err != nil {
		t.Fatal(err)
	}
	eg, err := EvaluatePRS(context.Background(), rg, EvalOptions{Phenotype: pheno, Backend: numwrap.BackendGonum})
	if err != nil {
		t.Fatal(err)
	}
	if !numwrap.WithinTol(es.Beta, eg.Beta) || !numwrap.WithinTol(es.SE, eg.SE) ||
		!numwrap.WithinTol(es.R2, eg.R2) {
		t.Fatalf("evaluate backends differ: %+v vs %+v", es, eg)
	}
}
