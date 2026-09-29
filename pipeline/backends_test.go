package pipeline

import (
	"context"
	"path/filepath"
	"testing"

	"github.com/lechandonga/statistical-genetics-go/internal/simulate"
	"github.com/lechandonga/statistical-genetics-go/numwrap"
)

// TestPipelineCrossCheckRunsBothBackends：整条流程在 CrossCheck 下
// 必须让 h2/meta/prs 各自双后端一致（不一致即整体失败，不挑一边）。
func TestPipelineCrossCheckRunsBothBackends(t *testing.T) {
	dir := t.TempDir()
	spec := simulate.Spec{
		Seed: 99, NVariants: 200, NSamples: 120,
		Cohorts: []simulate.CohortSpec{
			{Name: "a", N: 100},
			{Name: "b", N: 100},
		},
	}
	ds, err := simulate.Generate(spec)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := ds.WriteToDir(dir); err != nil {
		t.Fatal(err)
	}
	in := Input{
		Summaries: []SummaryInput{
			{Cohort: "a", Path: filepath.Join(dir, "a_sumstats.tsv"), N: 100},
			{Cohort: "b", Path: filepath.Join(dir, "b_sumstats.tsv"), N: 100},
		},
		LDScores:  filepath.Join(dir, "ldscores.tsv"),
		Genotypes: filepath.Join(dir, "genotypes.tsv"),
	}
	for _, b := range []numwrap.Backend{numwrap.BackendStd, numwrap.BackendGonum} {
		cfg := DefaultConfig()
		cfg.Backend = b
		cfg.CrossCheck = true
		res, err := Run(context.Background(), in, cfg)
		if err != nil {
			t.Fatalf("backend=%s crosscheck pipeline: %v", b, err)
		}
		if !res.Meta.CrossChecked || !res.PRS.CrossChecked {
			t.Fatalf("backend=%s stages must report crosschecked", b)
		}
		if res.Meta.Backend != b {
			t.Fatalf("meta backend=%s want %s", res.Meta.Backend, b)
		}
	}
	t.Log("[judgement] end-to-end cross-check passes for both primary backends")

	// 未知后端显式拒绝，不回退
	cfg := DefaultConfig()
	cfg.Backend = numwrap.Backend("cuda")
	if _, err := Run(context.Background(), in, cfg); err == nil {
		t.Fatal("unknown backend must fail the pipeline without fallback")
	}
}
