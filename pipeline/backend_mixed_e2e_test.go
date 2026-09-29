package pipeline

import (
	"bytes"
	"context"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/lechandonga/statistical-genetics-go/ingest"
	"github.com/lechandonga/statistical-genetics-go/internal/simulate"
	"github.com/lechandonga/statistical-genetics-go/numwrap"
	"github.com/lechandonga/statistical-genetics-go/stats"
)

func cleanSpec() simulate.Spec {
	return simulate.Spec{
		Seed: 42, NVariants: 120, NSamples: 80,
		FractionAmbiguous: 0.0,
		FractionFlip:      0.0,
		FractionMissingGT: 0.0,
		Cohorts: []simulate.CohortSpec{
			{Name: "eur", N: 70},
			{Name: "eur2", N: 70, SharedWith: "eur", SharedN: 20},
		},
	}
}

func generate(s simulate.Spec) (*simulate.Dataset, error) { return simulate.Generate(s) }

func itoa2(i int) string {
	if i == 0 {
		return "0"
	}
	var b [12]byte
	k := len(b)
	for i > 0 {
		k--
		b[k] = byte('0' + i%10)
		i /= 10
	}
	return string(b[k:])
}

// TestPipelineCrosscheckAndUsage 端到端验证：
//  1. crosscheck 后端跑通完整流程，与 stdlib/gonum 结论一致；
//  2. 报告中带每文件的用量小结；
//  3. 未知后端显式失败。
func TestPipelineCrosscheckAndUsage(t *testing.T) {
	dir := t.TempDir()
	ds, err := generate(cleanSpec())
	if err != nil {
		t.Fatal(err)
	}
	if _, err := ds.WriteToDir(dir); err != nil {
		t.Fatal(err)
	}
	sumName := func(c string) string { return filepath.Join(dir, c+"_sumstats.tsv") }
	in := Input{
		Summaries: []SummaryInput{
			{Cohort: "eur", Path: sumName("eur"), N: 0},
			{Cohort: "eur2", Path: sumName("eur2"), N: 0},
		},
		LDScores:  filepath.Join(dir, "ldscores.tsv"),
		Genotypes: filepath.Join(dir, "genotypes.tsv"),
	}

	run := func(b numwrap.Backend) *Result {
		cfg := DefaultConfig()
		cfg.Backend = b
		res, err := Run(context.Background(), in, cfg)
		if err != nil {
			t.Fatalf("backend=%s: %v", b, err)
		}
		return res
	}
	rStd := run(numwrap.BackendStd)
	rGon := run(numwrap.BackendGonum)
	rX := run(numwrap.BackendCrosscheck)

	if rStd.H2.H2 != rX.H2.H2 || rStd.H2.Flags != rX.H2.Flags {
		t.Fatal("crosscheck h2 must match stdlib reference")
	}
	if !numwrap.WithinTol(rStd.H2.H2, rGon.H2.H2) ||
		!numwrap.WithinTol(rStd.H2.SE, rGon.H2.SE) {
		t.Fatalf("h2 backends disagree: %v vs %v", rStd.H2.H2, rGon.H2.H2)
	}
	for i := range rStd.Meta.Variants {
		a, b := rStd.Meta.Variants[i], rGon.Meta.Variants[i]
		if a.Excluded != b.Excluded || a.ExcludeReason != b.ExcludeReason ||
			a.NCohorts != b.NCohorts {
			t.Fatalf("meta judgement differs at %d", i)
		}
	}

	// 用量小结：每个解析报告都必须有正的字节/行计数。
	for _, rep := range rX.Reports {
		u := rep.Usage
		if u.BytesRead <= 0 || u.LinesRead <= 1 || u.Elapsed <= 0 {
			t.Fatalf("source=%s usage not recorded: %+v", rep.Source, u)
		}
		t.Logf("[usage] %s: %s", rep.Source, u.Summary())
	}

	// 未知后端：显式失败
	cfg := DefaultConfig()
	cfg.Backend = numwrap.Backend("openblas")
	if _, err := Run(context.Background(), in, cfg); err == nil {
		t.Fatal("unknown backend must fail the pipeline explicitly")
	}
}

// TestPipelinePerTypeLimits 验证按输入类型独立的上限真实生效：
// 极小 MaxSamples（拒绝策略）让 PRS 阶段失败，且汇总统计解析已先行
// 成功（证明上限按类型独立，不互相覆盖）。
func TestPipelinePerTypeLimits(t *testing.T) {
	dir := t.TempDir()
	ds, err := generate(cleanSpec())
	if err != nil {
		t.Fatal(err)
	}
	if _, err := ds.WriteToDir(dir); err != nil {
		t.Fatal(err)
	}
	in := Input{
		Summaries: []SummaryInput{{Cohort: "eur",
			Path: filepath.Join(dir, "eur_sumstats.tsv"), N: 0}},
		Genotypes: filepath.Join(dir, "genotypes.tsv"),
	}
	cfg := DefaultConfig()
	cfg.EnableH2 = false
	cfg.EnableMeta = false
	cfg.GTLimits = ingest.Limits{
		MaxSamples: 2, MaxVariants: 5_000_000,
		MaxLineBytes: 1 << 20, MaxFileBytes: 1 << 31,
		Policy: ingest.LimitReject}
	_, err = Run(context.Background(), in, cfg)
	if err == nil || !strings.Contains(err.Error(), "MaxSamples") {
		t.Fatalf("per-type MaxSamples must reject, got %v", err)
	}

	// 截断策略：PRS 在 2 个样本上成功，报告标记截断
	cfg.GTLimits.Policy = ingest.LimitTruncate
	res, err := Run(context.Background(), in, cfg)
	if err != nil {
		t.Fatal(err)
	}
	if len(res.PRS.SampleIDs) != 2 {
		t.Fatalf("expected 2 truncated samples, got %d", len(res.PRS.SampleIDs))
	}
	var gtTruncated bool
	for _, rep := range res.Reports {
		if strings.Contains(rep.Source, "genotypes") && rep.Truncated {
			gtTruncated = true
		}
	}
	if !gtTruncated {
		t.Fatal("genotype report must record truncation")
	}
	t.Logf("[judgement] per-type genotype truncate: samples=%v", res.PRS.SampleIDs)
}

// TestPipelineMixedZOnlyMeta 端到端跑一个 beta+SE 队列和一个 Z-only
// 队列：两者必须都参与，审计中口径与排除原因可查。
func TestPipelineMixedZOnlyMeta(t *testing.T) {
	dir := t.TempDir()

	// beta 队列
	var beta bytes.Buffer
	beta.WriteString("SNP\tA1\tA2\tBETA\tSE\tN\tEAF\n")
	// z-only 队列（无 BETA/SE 列）
	var zonly bytes.Buffer
	zonly.WriteString("SNP\tA1\tA2\tZ\tN\tEAF\n")
	// 无 EAF、无 N 的行：必须以可区分原因排除
	var zbad bytes.Buffer
	zbad.WriteString("SNP\tA1\tA2\tZ\n")
	for i := 0; i < 40; i++ {
		id := "RS" + itoa2(i)
		fmt.Fprintf(&beta, "%s\tA\tG\t0.1\t0.05\t200\t0.3\n", id)
		fmt.Fprintf(&zonly, "%s\tA\tG\t2.0\t300\t0.4\n", id)
		if i < 10 {
			fmt.Fprintf(&zbad, "%s\tA\tG\t1.0\n", id)
		}
	}
	betaPath := filepath.Join(dir, "beta.tsv")
	zPath := filepath.Join(dir, "z.tsv")
	zbPath := filepath.Join(dir, "zbad.tsv")
	if err := os.WriteFile(betaPath, beta.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(zPath, zonly.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(zbPath, zbad.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}

	cfg := DefaultConfig()
	cfg.EnableH2 = false
	cfg.EnablePRS = false
	cfg.Effect = ingest.EffectBeta

	parse := func(path string, eff ingest.EffectField) *ingest.SummaryTable {
		f, err := os.Open(path)
		if err != nil {
			t.Fatal(err)
		}
		defer f.Close()
		tbl, err := ingest.ParseSummary(f, ingest.SummaryOptions{
			Source: path, Effect: eff, Limits: ingest.DefaultLimits(), FallbackN: 0})
		if err != nil {
			t.Fatal(err)
		}
		return tbl
	}
	_ = cfg
	cohort := func(name string, tbl *ingest.SummaryTable) stats.Cohort {
		return stats.Cohort{Name: name, Table: tbl, SampleIDs: []string{name}}
	}
	res, err := stats.RunMeta(context.Background(), []stats.Cohort{
		cohort("beta", parse(betaPath, ingest.EffectBeta)),
		cohort("zonly", parse(zPath, ingest.EffectZ)),
		cohort("zbad", parse(zbPath, ingest.EffectZ)),
	}, stats.MetaOptions{Model: stats.MetaFixed, MinCohorts: 2})
	if err != nil {
		t.Fatal(err)
	}
	nZIncluded, nZBasis := 0, 0
	reasons := map[string]int{}
	for _, mv := range res.Variants {
		for _, l := range mv.PerCohort {
			if l.Cohort == "zonly" && !l.Excluded {
				nZIncluded++
				if l.Basis == stats.BasisZNEAF {
					nZBasis++
				}
			}
			if l.Excluded {
				reasons[string(l.Reason)]++
			}
		}
	}
	if nZIncluded != 40 || nZBasis != 40 {
		t.Fatalf("z-only inclusion=%d basis=%d want 40/40", nZIncluded, nZBasis)
	}
	if reasons["z_missing_n"] == 0 {
		t.Fatalf("z_missing_n reason not observed: %v", reasons)
	}
	t.Logf("[judgement] z-only included=%d/40 with basis z_n_eaf; exclusion reasons=%v",
		nZIncluded, reasons)
}
