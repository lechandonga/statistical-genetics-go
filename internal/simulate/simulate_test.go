package simulate

import (
	"bytes"
	"testing"
)

func demoSpec() Spec {
	return Spec{
		Seed: 42, NVariants: 120, NSamples: 80,
		FractionAmbiguous: 0.05, FractionFlip: 0.25, FractionSwap: 0.1,
		FractionMissingGT: 0.03,
		Cohorts: []CohortSpec{
			{Name: "c1", N: 60},
			{Name: "c2", N: 60, SharedWith: "c1", SharedN: 20},
		},
		DuplicateVariants: 1,
		InvalidRows:       2,
	}
}

func TestGenerateDeterministicBytes(t *testing.T) {
	d1, err := Generate(demoSpec())
	if err != nil {
		t.Fatal(err)
	}
	d2, err := Generate(demoSpec())
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(d1.GenotypeFile, d2.GenotypeFile) {
		t.Fatal("genotype file bytes differ across runs with same seed")
	}
	if !bytes.Equal(d1.LDScoreFile, d2.LDScoreFile) {
		t.Fatal("ld file bytes differ")
	}
	for name := range d1.SummaryFiles {
		if !bytes.Equal(d1.SummaryFiles[name], d2.SummaryFiles[name]) {
			t.Fatalf("summary %s bytes differ", name)
		}
	}
	t.Logf("[input] seed=42 fixed [judgement] genotype/ld/%d summaries byte-identical across runs",
		len(d1.SummaryFiles))
}

func TestGenerateTruthMarkers(t *testing.T) {
	d, err := Generate(demoSpec())
	if err != nil {
		t.Fatal(err)
	}
	if len(d.Truth.AmbiguousIDs) == 0 {
		t.Fatal("expected ambiguous sites")
	}
	if len(d.Truth.FlippedIDs) == 0 || len(d.Truth.SwappedIDs) == 0 {
		t.Fatalf("expected flip/swap markers: %+v", d.Truth)
	}
	if len(d.Truth.OverlapPairs) != 1 {
		t.Fatalf("overlap pairs=%v", d.Truth.OverlapPairs)
	}
	t.Logf("[input] spec flip=25%% swap=10%% overlap=20/60 [judgement] ambiguous=%d flipped=%d swapped=%d overlap=%v",
		len(d.Truth.AmbiguousIDs), len(d.Truth.FlippedIDs), len(d.Truth.SwappedIDs),
		d.Truth.OverlapPairs)
}

func TestGenerateRejectsBadSpec(t *testing.T) {
	if _, err := Generate(demoSpec()); err != nil {
		t.Fatal(err)
	}
	bad := demoSpec()
	bad.NVariants = 0
	if _, err := Generate(bad); err == nil {
		t.Fatal("zero variants must reject")
	}
	bad = demoSpec()
	bad.Cohorts[1].SharedN = 1000 // 超过基础队列
	if _, err := Generate(bad); err == nil {
		t.Fatal("SharedN > base cohort must reject")
	}
	t.Logf("[judgement] zero variants and oversized overlap rejected with explicit errors")
}
