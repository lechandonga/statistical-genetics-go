package meta

import (
	"math"
	"testing"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/sumstats"
)

func parseTwo(t *testing.T) []Cohort {
	t.Helper()
	a, err := sumstats.ParseFile("../../../testdata/cohort_A.tsv",
		sumstats.ParseConfig{CohortName: "A", HasHeader: true, SkipReject: true})
	if err != nil {
		t.Fatal(err)
	}
	b, err := sumstats.ParseFile("../../../testdata/cohort_B.tsv",
		sumstats.ParseConfig{CohortName: "B", HasHeader: true, SkipReject: true})
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("INPUT A: %s", a.Report.Summary())
	t.Logf("INPUT B: %s", b.Report.Summary())
	return []Cohort{{Summary: a}, {Summary: b}}
}

func TestFixedEffectAlignmentAndPool(t *testing.T) {
	cohorts := parseTwo(t)
	an, err := Prepare(cohorts, Config{Model: FixedInverseVariance, OverlapPolicy: OverlapFlag})
	if err != nil {
		t.Fatal(err)
	}
	res, err := an.PoolAll(Config{Model: FixedInverseVariance, OverlapPolicy: OverlapFlag})
	if err != nil {
		t.Fatal(err)
	}
	if len(res) == 0 {
		t.Fatal("no pooled variants")
	}
	for _, r := range res {
		t.Logf("POOL %s beta=%.4f se=%.4f p=%.3g q=%.3f i2=%.3f tau2=%.4f warnings=%d",
			r.Key, r.Beta, r.SE, r.PValue, r.Q, r.I2, r.Tau2, len(r.Warnings))
		for _, e := range r.Effects {
			t.Logf("   cohort=%s beta=%.3f sign=%v orient=%s",
				e.Cohort, e.Beta, e.Aligned.Sign, e.Aligned.Orient)
		}
	}
	// Variant 1:100 A/G in cohort A vs G/A in cohort B: B must be sign-flipped,
	// so pooled beta stays positive (both cohorts agree after alignment).
	var found bool
	for _, r := range res {
		if string(r.Key) == "1:100:A:G" {
			found = true
			var bBeta float64
			for _, e := range r.Effects {
				if e.Cohort == "B" {
					bBeta = e.Beta
				}
			}
			if bBeta <= 0 || math.Abs(bBeta-0.09) > 1e-9 {
				t.Fatalf("cohort B beta should flip from -0.09 to +0.09, got %v", bBeta)
			}
			// after flip it is +0.09; pooled should be positive
			if r.Beta <= 0 {
				t.Fatalf("pooled beta after alignment = %v, want > 0", r.Beta)
			}
		}
	}
	if !found {
		t.Fatal("variant 1:100 missing from pooled results")
	}
}

func TestRandomEffects(t *testing.T) {
	cohorts := parseTwo(t)
	an, err := Prepare(cohorts, Config{OverlapPolicy: OverlapFlag})
	if err != nil {
		t.Fatal(err)
	}
	res, err := an.PoolAll(Config{Model: RandomDerSimonianLaird, OverlapPolicy: OverlapFlag})
	if err != nil {
		t.Fatal(err)
	}
	for _, r := range res {
		if r.Tau2 < 0 {
			t.Fatalf("tau2 negative: %v", r.Tau2)
		}
	}
	t.Logf("random-effects pooled %d variants; first tau2=%.4f", len(res), res[0].Tau2)
}

func TestDuplicateCohortRejected(t *testing.T) {
	cohorts := parseTwo(t)
	cohorts[0].OverlapKey = "same-study"
	cohorts[1].OverlapKey = "same-study"
	_, err := Prepare(cohorts, Config{OverlapPolicy: OverlapReject})
	if !stgerr.Is(err, stgerr.KindDuplicateCohort) {
		t.Fatalf("expected duplicate_cohort, got %v", err)
	}
	t.Logf("JUDGE duplicate: %v", err)
}

func TestSampleOverlapFlaggedAndDedup(t *testing.T) {
	cohorts := parseTwo(t)
	cohorts[0].SampleIDs = []string{"s1", "s2", "s3"}
	cohorts[1].SampleIDs = []string{"s3", "s4", "s5"}
	an, err := Prepare(cohorts, Config{OverlapPolicy: OverlapFlag})
	if err != nil {
		t.Fatal(err)
	}
	if len(an.Overlaps) != 1 || an.Overlaps[0].Shared != 1 ||
		an.Overlaps[0].Kind != OverlapSamples {
		t.Fatalf("overlap finding wrong: %+v", an.Overlaps)
	}
	t.Logf("JUDGE overlap: %+v", an.Overlaps[0])
	r, err := an.PoolAll(Config{Model: FixedInverseVariance, OverlapPolicy: OverlapFlag})
	if err != nil {
		t.Fatal(err)
	}
	if len(r[0].Warnings) == 0 {
		t.Fatal("overlapping pooled result must carry a warning")
	}

	// Deduplicate: second cohort dropped entirely.
	an2, err := Prepare(cohorts, Config{OverlapPolicy: OverlapDeduplicate})
	if err != nil {
		t.Fatal(err)
	}
	if len(an2.Effective) != 1 || len(an2.Dropped) != 1 {
		t.Fatalf("dedup effective=%d dropped=%v", len(an2.Effective), an2.Dropped)
	}
}

func TestUnsupportedModel(t *testing.T) {
	cohorts := parseTwo(t)
	an, err := Prepare(cohorts, Config{OverlapPolicy: OverlapFlag})
	if err != nil {
		t.Fatal(err)
	}
	_, err = an.PoolAll(Config{Model: "bayes", OverlapPolicy: OverlapFlag})
	if !stgerr.Is(err, stgerr.KindUnsupportedOption) {
		t.Fatalf("unsupported model: %v", err)
	}
}

func TestParallelDeterminism(t *testing.T) {
	cohorts := parseTwo(t)
	an, _ := Prepare(cohorts, Config{OverlapPolicy: OverlapFlag})
	r1, err := an.PoolAll(Config{Model: RandomDerSimonianLaird, Workers: 1, OverlapPolicy: OverlapFlag})
	if err != nil {
		t.Fatal(err)
	}
	r4, err := an.PoolAll(Config{Model: RandomDerSimonianLaird, Workers: 4, OverlapPolicy: OverlapFlag})
	if err != nil {
		t.Fatal(err)
	}
	if len(r1) != len(r4) {
		t.Fatal("length drift")
	}
	for i := range r1 {
		if r1[i].Beta != r4[i].Beta || r1[i].SE != r4[i].SE {
			t.Fatalf("variant %d drift: serial=%v/%v parallel=%v/%v",
				i, r1[i].Beta, r1[i].SE, r4[i].Beta, r4[i].SE)
		}
	}
}
