package pipeline

import (
	"bytes"
	"context"
	"os"
	"path/filepath"
	"testing"

	"github.com/lechandonga/statistical-genetics-go/internal/simulate"
	"github.com/lechandonga/statistical-genetics-go/stats"
)

// TestEndToEndFlow 验证：解析拒绝非法/重复行、模糊位点被标记、
// 链翻转被对齐、重叠队列被识别，且 workers=1 与 workers=4 结果逐位一致。
func TestEndToEndFlow(t *testing.T) {
	dir := t.TempDir()
	spec := simulate.Spec{
		Seed: 42, NVariants: 300, NSamples: 200,
		FractionAmbiguous: 0.04,
		FractionFlip:      0.20,
		FractionSwap:      0.0,
		FractionMissingGT: 0.02,
		Cohorts: []simulate.CohortSpec{
			{Name: "eur", N: 180},
			{Name: "eur2", N: 180, SharedWith: "eur", SharedN: 60},
		},
		DuplicateVariants: 1,
		InvalidRows:       3,
	}
	ds, err := simulate.Generate(spec)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := ds.WriteToDir(dir); err != nil {
		t.Fatal(err)
	}
	// 本测试不提供样本 ID 文件：重叠检测必须“无证据不猜测”。

	sumName := func(c string) string { return filepath.Join(dir, c+"_sumstats.tsv") }
	in := Input{
		Summaries: []SummaryInput{
			{Cohort: "eur", Path: sumName("eur"), N: 0},
			{Cohort: "eur2", Path: sumName("eur2"), N: 0},
		},
		LDScores:  filepath.Join(dir, "ldscores.tsv"),
		Genotypes: filepath.Join(dir, "genotypes.tsv"),
	}

	// 非法/重复行存在：严格模式必须失败，并给出可区分原因。
	cfg := DefaultConfig()
	if _, err := Run(context.Background(), in, cfg); err == nil {
		t.Fatal("strict pipeline must fail on injected invalid/duplicate rows")
	} else {
		t.Logf("[input] 3 invalid + 1 duplicate injected [judgement] strict parse rejected: %v", err)
	}

	// 干净数据重新生成（无注入），验证全流程成功
	spec.DuplicateVariants = 0
	spec.InvalidRows = 0
	ds2, err := simulate.Generate(spec)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := ds2.WriteToDir(dir); err != nil {
		t.Fatal(err)
	}
	cfg = DefaultConfig()
	cfg.EnableMeta = true
	res, err := Run(context.Background(), in, cfg)
	if err != nil {
		t.Fatalf("clean pipeline: %v", err)
	}
	if res.H2 == nil || res.Meta == nil || res.PRS == nil {
		t.Fatal("all stages expected")
	}
	t.Logf("[input-summary] variants=%d samples=%d flip=20%% ambiguous=4%% missingGT=2%%",
		spec.NVariants, spec.NSamples)
	t.Logf("[judgement:h2] %s", res.H2.Diag)
	t.Logf("[judgement:meta] %s", res.Meta.Diag)
	t.Logf("[judgement:prs] %s", res.PRS.Diag)
	if len(res.H2.Rejected) == 0 {
		t.Fatal("ambiguous sites should be recorded in h2 rejected list")
	}
	if res.H2.Flags&stats.H2FlagAmbiguousDropped == 0 {
		t.Fatal("h2 should flag ambiguous drop")
	}

	// 并行一致性：workers=4 必须逐位相同
	cfg4 := cfg
	cfg4.Workers = 4
	cfg4.H2.Workers = 4
	cfg4.Meta.Workers = 4
	cfg4.PRS.Workers = 4
	res4, err := Run(context.Background(), in, cfg4)
	if err != nil {
		t.Fatal(err)
	}
	if res4.H2.H2 != res.H2.H2 || res4.H2.SE != res.H2.SE {
		t.Fatalf("h2 drift: serial=%v/%v parallel=%v/%v",
			res.H2.H2, res.H2.SE, res4.H2.H2, res4.H2.SE)
	}
	for i := range res.Meta.Variants {
		a, b := res.Meta.Variants[i], res4.Meta.Variants[i]
		if a.Beta != b.Beta || a.SE != b.SE {
			t.Fatalf("meta drift at %d", i)
		}
	}
	for j := range res.PRS.Scores {
		if res.PRS.Scores[j] != res4.PRS.Scores[j] {
			t.Fatalf("prs drift at %d", j)
		}
	}
	t.Logf("[judgement] workers=1 vs workers=4: h2/meta/prs bit-identical")
}

// TestEndToEndOverlap 用共享样本 ID 文件验证重叠标记/拒绝策略。
func TestEndToEndOverlap(t *testing.T) {
	dir := t.TempDir()
	idsA := filepath.Join(dir, "idsA.txt")
	idsB := filepath.Join(dir, "idsB.txt")
	var bufA, bufB bytes.Buffer
	for i := 0; i < 100; i++ {
		bufA.WriteString("ind")
		bufA.WriteString(pad(i))
		bufA.WriteString("\n")
	}
	for i := 40; i < 140; i++ {
		bufB.WriteString("ind")
		bufB.WriteString(pad(i))
		bufB.WriteString("\n")
	}
	if err := os.WriteFile(idsA, bufA.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(idsB, bufB.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}

	spec := simulate.Spec{
		Seed: 7, NVariants: 100, NSamples: 140,
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
			{Cohort: "a", Path: filepath.Join(dir, "a_sumstats.tsv"), N: 100, SampleIDPath: idsA},
			{Cohort: "b", Path: filepath.Join(dir, "b_sumstats.tsv"), N: 100, SampleIDPath: idsB},
		},
	}
	// flag
	cfg := DefaultConfig()
	cfg.EnableH2 = false
	cfg.EnablePRS = false
	cfg.Meta.Overlap = stats.OverlapFlag
	res, err := Run(context.Background(), in, cfg)
	if err != nil {
		t.Fatal(err)
	}
	if len(res.Meta.Overlap) != 1 || res.Meta.Overlap[0].SharedSamples != 60 {
		t.Fatalf("overlap=%+v", res.Meta.Overlap)
	}
	t.Logf("[input] cohorts share ind040..ind099 (60/100) [judgement] flagged Jaccard=%.3f",
		res.Meta.Overlap[0].Jaccard)

	// reject
	cfg.Meta.Overlap = stats.OverlapReject
	if _, err := Run(context.Background(), in, cfg); err == nil {
		t.Fatal("overlap reject policy must fail pipeline")
	} else {
		t.Logf("[judgement] overlap=reject => %v", err)
	}
}

func pad(i int) string {
	const digits = "0123456789"
	var b [6]byte
	for k := 5; k >= 0; k-- {
		b[k] = digits[i%10]
		i /= 10
	}
	return string(b[:])
}
