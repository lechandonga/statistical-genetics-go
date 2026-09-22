package pipeline

import (
	"math"
	"testing"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/genotype"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/ldsc"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/meta"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/prs"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
)

func baseCfg() Config {
	return Config{
		SumstatFiles: []string{"../../../testdata/cohort_A.tsv", "../../../testdata/cohort_B.tsv"},
		CohortNames:  []string{"A", "B"},
		LDScoreFile:  "../../../testdata/ldscores.tsv",
		GenotypeFile: "../../../testdata/genotypes.tsv",
		LDConfig:     ldsc.Config{},
		MetaConfig: meta.Config{Model: meta.FixedInverseVariance,
			OverlapPolicy: meta.OverlapFlag, MinCohorts: 1},
		PRSConfig: prs.Config{
			MissingGenotype: genotype.DosageMeanImpute,
			MissingLocus:    prs.MissingLocusSkip,
		},
	}
}

func TestEndToEndReproducible(t *testing.T) {
	r1, err := Run(baseCfg())
	if err != nil {
		t.Fatal(err)
	}
	for _, rep := range r1.Reports {
		t.Logf("INPUT %s", rep.Summary())
	}
	for _, h := range r1.H2 {
		t.Logf("H2 %s h2=%.4f se=%.4f status=%s", h.Cohort, h.Result.H2, h.Result.SE, h.Result.Status)
	}
	for _, ov := range r1.Overlaps {
		t.Logf("OVERLAP %+v", ov)
	}
	if r1.Scores == nil || len(r1.Scores.Values) != 6 {
		t.Fatal("scores missing")
	}

	r2, err := Run(baseCfg())
	if err != nil {
		t.Fatal(err)
	}
	for i := range r1.Meta {
		if math.Float64bits(r1.Meta[i].Beta) != math.Float64bits(r2.Meta[i].Beta) {
			t.Fatal("meta result not reproducible")
		}
	}
	for i := range r1.Scores.Values {
		if math.Float64bits(r1.Scores.Values[i]) != math.Float64bits(r2.Scores.Values[i]) {
			t.Fatal("score not reproducible")
		}
	}
}

func TestEndToEndParallel(t *testing.T) {
	c1 := baseCfg()
	rSerial, err := Run(c1)
	if err != nil {
		t.Fatal(err)
	}
	c2 := baseCfg()
	c2.LDConfig.Workers = 4
	c2.MetaConfig.Workers = 4
	c2.PRSConfig.Workers = 4
	rPar, err := Run(c2)
	if err != nil {
		t.Fatal(err)
	}
	for i := range rSerial.Scores.Values {
		if math.Float64bits(rSerial.Scores.Values[i]) != math.Float64bits(rPar.Scores.Values[i]) {
			t.Fatalf("score drift sample %d: %v vs %v",
				i, rSerial.Scores.Values[i], rPar.Scores.Values[i])
		}
	}
	for i := range rSerial.Meta {
		if math.Float64bits(rSerial.Meta[i].Beta) != math.Float64bits(rPar.Meta[i].Beta) {
			t.Fatal("meta drift under parallelism")
		}
	}
}

func TestBudgetVariantLimit(t *testing.T) {
	cfg := baseCfg()
	cfg.Limits.MaxVariants = 3
	_, err := Run(cfg)
	if !stgerr.Is(err, stgerr.KindLimitExceeded) {
		t.Fatalf("expected limit exceeded, got %v", err)
	}
	t.Logf("JUDGE budget: %v", err)
}

func TestBudgetSampleLimit(t *testing.T) {
	cfg := baseCfg()
	cfg.Limits.MaxSamples = 2
	_, err := Run(cfg)
	if !stgerr.Is(err, stgerr.KindLimitExceeded) {
		t.Fatalf("expected sample limit, got %v", err)
	}
}

func TestBudgetCohortLimit(t *testing.T) {
	cfg := baseCfg()
	cfg.Limits.MaxCohorts = 1
	_, err := Run(cfg)
	if !stgerr.Is(err, stgerr.KindLimitExceeded) {
		t.Fatalf("expected cohort limit, got %v", err)
	}
}

func TestDuplicateCohortOverlapPolicy(t *testing.T) {
	cfg := baseCfg()
	cfg.OverlapKeys = []string{"study1", "study1"}
	cfg.MetaConfig.OverlapPolicy = meta.OverlapReject
	_, err := Run(cfg)
	if !stgerr.Is(err, stgerr.KindDuplicateCohort) {
		t.Fatalf("expected duplicate cohort, got %v", err)
	}
	t.Logf("JUDGE duplicate: %v", err)
}

func TestBudgetConsumeConcurrencySafe(t *testing.T) {
	b := NewBudget(Limits{MaxVariants: 1_000_000, MaxCohorts: 100,
		MaxSamples: 1_000_000, MaxFileBytes: 1 << 40})
	done := make(chan error, 100)
	for g := 0; g < 100; g++ {
		go func() {
			for i := 0; i < 100; i++ {
				if err := b.Consume(BudgetVariants, 1); err != nil {
					done <- err
					return
				}
			}
			done <- nil
		}()
	}
	for i := 0; i < 100; i++ {
		if err := <-done; err != nil {
			t.Fatal(err)
		}
	}
	if b.Remaining(BudgetVariants) != 990_000 {
		t.Fatalf("race/lost update: remaining=%d", b.Remaining(BudgetVariants))
	}
}
