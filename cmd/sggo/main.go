// Command sggo 是本地遗传统计流程的命令行入口。
//
// 子命令：
//
//	sggo simulate -out <dir>         生成确定性演示数据集
//	sggo h2 -sum f -ld f [-n N]      估计 SNP 遗传力
//	sggo meta -in name=path ...      多队列元分析
//	sggo prs -sum f -gt f            构建个体风险评分（输出每样本 TSV）
//	sggo version                     打印版本与数值后端信息
package main

import (
	"context"
	"flag"
	"fmt"
	"os"
	"strings"

	"github.com/lechandonga/statistical-genetics-go/ingest"
	"github.com/lechandonga/statistical-genetics-go/internal/simulate"
	"github.com/lechandonga/statistical-genetics-go/numwrap"
	"github.com/lechandonga/statistical-genetics-go/pipeline"
	"github.com/lechandonga/statistical-genetics-go/stats"
	"github.com/lechandonga/statistical-genetics-go/variant"
)

func main() {
	if err := run(os.Args[1:]); err != nil {
		fmt.Fprintln(os.Stderr, "sggo:", err)
		os.Exit(1)
	}
}

func run(args []string) error {
	if len(args) == 0 {
		usage()
		return nil
	}
	switch args[0] {
	case "simulate":
		return cmdSimulate(args[1:])
	case "h2":
		return cmdH2(args[1:])
	case "meta":
		return cmdMeta(args[1:])
	case "prs":
		return cmdPRS(args[1:])
	case "version":
		v, _ := numwrap.Version(numwrap.DefaultBackend())
		fmt.Printf("sggo dev\ngo-module: github.com/lechandonga/statistical-genetics-go\nnumerical backend: %s (%s)\n",
			numwrap.DefaultBackend(), v)
		return nil
	default:
		usage()
		return fmt.Errorf("unknown subcommand %q", args[0])
	}
}

func usage() {
	fmt.Fprintln(os.Stderr, "usage: sggo {simulate|h2|meta|prs|version} [flags]")
}

func cmdSimulate(args []string) error {
	fs := flag.NewFlagSet("simulate", flag.ContinueOnError)
	out := fs.String("out", "testdata", "output directory")
	nv := fs.Int("variants", 300, "number of variants")
	ns := fs.Int("samples", 200, "number of genotype samples")
	seed := fs.Int64("seed", 42, "deterministic RNG seed")
	amb := fs.Float64("ambiguous", 0.04, "fraction of A/T,C/G sites")
	flip := fs.Float64("flip", 0.2, "fraction strand-flipped in sumstats")
	overlap := fs.Int("overlap", 60, "shared samples between the two demo cohorts")
	if err := fs.Parse(args); err != nil {
		return err
	}
	spec := simulate.Spec{
		Seed: *seed, NVariants: *nv, NSamples: *ns,
		FractionAmbiguous: *amb, FractionFlip: *flip, FractionMissingGT: 0.02,
		Cohorts: []simulate.CohortSpec{
			{Name: "cohort1", N: *ns - 20},
			{Name: "cohort2", N: *ns - 20, SharedWith: "cohort1", SharedN: *overlap},
		},
	}
	ds, err := simulate.Generate(spec)
	if err != nil {
		return err
	}
	paths, err := ds.WriteToDir(*out)
	if err != nil {
		return err
	}
	fmt.Println("generated deterministic dataset in", *out)
	for _, p := range paths {
		fmt.Println("  ", p)
	}
	fmt.Printf("truth: h2=%.2f ambiguous=%d flipped=%d overlapPairs=%d\n",
		ds.Truth.H2, len(ds.Truth.AmbiguousIDs), len(ds.Truth.FlippedIDs),
		len(ds.Truth.OverlapPairs))
	return nil
}

func openSummary(path string, n float64) (*ingest.SummaryTable, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return ingest.ParseSummary(f, ingest.SummaryOptions{
		Source: path, Effect: ingest.EffectBeta, FallbackN: n,
		Limits: ingest.DefaultLimits(),
	})
}

func openLD(path string) (*ingest.LDScoreTable, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return ingest.ParseLDScores(f, ingest.LDSOptions{Source: path, Limits: ingest.DefaultLimits()})
}

func cmdH2(args []string) error {
	fs := flag.NewFlagSet("h2", flag.ContinueOnError)
	sumPath := fs.String("sum", "", "summary statistics TSV")
	ldPath := fs.String("ld", "", "LD score TSV")
	n := fs.Float64("n", 0, "fallback sample size when N column absent")
	free := fs.Bool("free-intercept", true, "estimate intercept freely")
	boundary := fs.String("boundary", "report", "report|clamp|reject for h2 outside [0,h2max]")
	workers := fs.Int("workers", 1, "shard parallelism (results identical to serial)")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if *sumPath == "" || *ldPath == "" {
		return fmt.Errorf("-sum and -ld are required")
	}
	sum, err := openSummary(*sumPath, *n)
	if err != nil {
		return err
	}
	ld, err := openLD(*ldPath)
	if err != nil {
		return err
	}
	pol := stats.H2ReportRaw
	switch *boundary {
	case "report":
	case "clamp":
		pol = stats.H2Clamp
	case "reject":
		pol = stats.H2Reject
	default:
		return fmt.Errorf("unsupported boundary policy %q", *boundary)
	}
	res, err := stats.EstimateH2(context.Background(), sum, ld, stats.H2Options{
		NBar: *n, FreeIntercept: *free, Boundary: pol,
		Align: variant.DefaultAlignConfig(), Workers: *workers,
	})
	if err != nil {
		return err
	}
	fmt.Println(res.Diag)
	fmt.Printf("h2=%.6f se=%.6f ci95=[%.4f,%.4f] intercept=%.4f nsnps=%d flags=%s\n",
		res.H2, res.SE, res.CI95[0], res.CI95[1], res.Intercept, res.NSnps, res.Flags)
	return nil
}

// cohortFlag 支持重复的 -in name=path[:sampleIDFile]
type cohortFlag []pipeline.SummaryInput

func (c *cohortFlag) String() string {
	parts := make([]string, len(*c))
	for i, x := range *c {
		parts[i] = x.Cohort + "=" + x.Path
	}
	return strings.Join(parts, ",")
}

func (c *cohortFlag) Set(v string) error {
	eq := strings.Index(v, "=")
	if eq <= 0 {
		return fmt.Errorf("-in expects name=path[:idfile], got %q", v)
	}
	name, rest := v[:eq], v[eq+1:]
	si := pipeline.SummaryInput{Cohort: name}
	if colon := strings.Index(rest, ":"); colon >= 0 {
		si.Path, si.SampleIDPath = rest[:colon], rest[colon+1:]
	} else {
		si.Path = rest
	}
	*c = append(*c, si)
	return nil
}

func cmdMeta(args []string) error {
	fs := flag.NewFlagSet("meta", flag.ContinueOnError)
	var ins cohortFlag
	fs.Var(&ins, "in", "cohort input name=path[:sampleIDFile] (repeatable)")
	model := fs.String("model", "fixed", "fixed|random")
	overlap := fs.String("overlap", "flag", "flag|drop|reject")
	workers := fs.Int("workers", 1, "parallelism (identical results)")
	out := fs.String("out", "", "optional output TSV path")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if len(ins) < 2 {
		return fmt.Errorf("at least two -in cohorts required")
	}
	cohorts := make([]stats.Cohort, 0, len(ins))
	for _, si := range ins {
		tbl, err := openSummary(si.Path, 0)
		if err != nil {
			return fmt.Errorf("cohort %s: %w", si.Cohort, err)
		}
		var ids []string
		if si.SampleIDPath != "" {
			b, err := os.ReadFile(si.SampleIDPath)
			if err != nil {
				return err
			}
			for _, line := range strings.Fields(string(b)) {
				ids = append(ids, line)
			}
		}
		cohorts = append(cohorts, stats.Cohort{Name: si.Cohort, SampleIDs: ids, Table: tbl})
	}
	m := stats.MetaFixed
	if *model == "random" {
		m = stats.MetaRandomDL
	} else if *model != "fixed" {
		return fmt.Errorf("unsupported model %q", *model)
	}
	ov := stats.OverlapFlag
	switch *overlap {
	case "flag":
	case "drop":
		ov = stats.OverlapDrop
	case "reject":
		ov = stats.OverlapReject
	default:
		return fmt.Errorf("unsupported overlap policy %q", *overlap)
	}
	res, err := stats.RunMeta(context.Background(), cohorts, stats.MetaOptions{
		Model: m, Overlap: ov, Align: variant.DefaultAlignConfig(), Workers: *workers,
	})
	if err != nil {
		return err
	}
	fmt.Println(res.Diag)
	for _, e := range res.Overlap {
		fmt.Printf("overlap: %s/%s shared=%d jaccard=%.3f duplicate=%v\n",
			e.CohortA, e.CohortB, e.SharedSamples, e.Jaccard, e.Duplicate)
	}
	if *out != "" {
		return writeMetaTSV(*out, res)
	}
	nShow := 5
	if len(res.Variants) < nShow {
		nShow = len(res.Variants)
	}
	fmt.Println("variant\tbeta\tse\tp\tncohorts\tQ\tI2\texcluded")
	for i := 0; i < nShow; i++ {
		mv := res.Variants[i]
		fmt.Printf("%s\t%.5g\t%.5g\t%.3g\t%d\t%.3f\t%.3f\t%v\n",
			mv.V.ID, mv.Beta, mv.SE, mv.P, mv.NCohorts, mv.Q, mv.I2, mv.Excluded)
	}
	return nil
}

func writeMetaTSV(path string, res *stats.MetaResult) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	fmt.Fprintln(f, "variant\tbeta\tse\tci_low\tci_high\tp\tncohorts\tQ\tdf\tpQ\tI2\ttau2\texcluded")
	for _, mv := range res.Variants {
		fmt.Fprintf(f, "%s\t%.10g\t%.10g\t%.10g\t%.10g\t%.6g\t%d\t%.6g\t%d\t%.6g\t%.6g\t%.6g\t%v\n",
			mv.V.ID, mv.Beta, mv.SE, mv.CI95[0], mv.CI95[1], mv.P, mv.NCohorts,
			mv.Q, mv.Df, mv.PQ, mv.I2, mv.Tau2, mv.Excluded)
	}
	return nil
}

func cmdPRS(args []string) error {
	fs := flag.NewFlagSet("prs", flag.ContinueOnError)
	sumPath := fs.String("sum", "", "model summary statistics TSV")
	gtPath := fs.String("gt", "", "genotype matrix TSV")
	n := fs.Float64("n", 0, "fallback sample size")
	scale := fs.String("scale", "raw", "raw|std|perallele")
	missingModel := fs.String("missing-sites", "reject", "reject|ignore")
	missingGT := fs.String("missing-gt", "impute", "impute|exclude")
	pthr := fs.Float64("p", 0, "P-value threshold (0 = no filter)")
	enc := fs.String("encoding", "012", "012|acgt")
	out := fs.String("out", "", "output TSV path (required)")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if *sumPath == "" || *gtPath == "" || *out == "" {
		return fmt.Errorf("-sum, -gt and -out are required")
	}
	model, err := openSummary(*sumPath, *n)
	if err != nil {
		return err
	}
	f, err := os.Open(*gtPath)
	if err != nil {
		return err
	}
	encoding := ingest.GTDosage012
	if *enc == "acgt" {
		encoding = ingest.GTACGT
	} else if *enc != "012" {
		f.Close()
		return fmt.Errorf("unsupported encoding %q", *enc)
	}
	gt, err := ingest.ParseGenotype(f, ingest.GenotypeOptions{
		Source: *gtPath, Encoding: encoding, Limits: ingest.DefaultLimits()})
	f.Close()
	if err != nil {
		return err
	}
	mp := stats.ModelStrict
	if *missingModel == "ignore" {
		mp = stats.ModelIgnore
	} else if *missingModel != "reject" {
		return fmt.Errorf("unsupported missing-sites policy %q", *missingModel)
	}
	gp := stats.GTMeanImpute
	if *missingGT == "exclude" {
		gp = stats.GTExclude
	} else if *missingGT != "impute" {
		return fmt.Errorf("unsupported missing-gt policy %q", *missingGT)
	}
	res, err := stats.BuildPRS(context.Background(), model, gt, stats.PRSOptions{
		ScoreScale: *scale, Model: mp, GT: gp, PThreshold: *pthr,
		Align: variant.DefaultAlignConfig(),
	})
	if err != nil {
		return err
	}
	outf, err := os.Create(*out)
	if err != nil {
		return err
	}
	defer outf.Close()
	fmt.Fprintln(outf, "sample\tscore\tactive_sites")
	for i, id := range res.SampleIDs {
		fmt.Fprintf(outf, "%s\t%.10g\t%d\n", id, res.Scores[i], res.NActive[i])
	}
	fmt.Println(res.Diag)
	fmt.Println("scores written to", *out)
	return nil
}
