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
	"bufio"
	"bytes"
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
		fmt.Printf("sggo dev\ngo-module: github.com/lechandonga/statistical-genetics-go\n")
		fmt.Printf("default numerical backend: %s\n", numwrap.DefaultBackend())
		for _, b := range numwrap.AvailableBackends() {
			v, _ := numwrap.Version(b)
			fmt.Printf("backend: %-7s %s\n", b, v)
		}
		fmt.Printf("cross-backend tolerance: rel=%.0e abs=%.0e\n",
			numwrap.CrossCheckRelTol, numwrap.CrossCheckAbsTol)
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

// parseBackend 把 -backend 文本解析为已编译后端；
// 未知/未编译后端显式报错，绝不静默退回默认实现。
func parseBackend(name string) (numwrap.Backend, error) {
	if name == "" || name == string(numwrap.DefaultBackend()) {
		return numwrap.DefaultBackend(), nil
	}
	b := numwrap.Backend(name)
	if !numwrap.Available(b) {
		return "", fmt.Errorf("unsupported backend %q (available: stdlib, gonum); "+
			"no silent fallback to the default", name)
	}
	return b, nil
}

// limitFlags 注册三类输入共用的规模上限 flag，返回 *ingest.Limits 与 policy 文本。
func limitFlags(fs *flag.FlagSet) (*ingest.Limits, *string) {
	lim := ingest.DefaultLimits()
	fs.IntVar(&lim.MaxVariants, "max-variants", lim.MaxVariants, "max accepted variant rows")
	fs.IntVar(&lim.MaxSamples, "max-samples", lim.MaxSamples, "max genotype sample columns")
	fs.Int64Var(&lim.MaxLineBytes, "max-line-bytes", lim.MaxLineBytes, "max bytes per line")
	fs.Int64Var(&lim.MaxFileBytes, "max-file-bytes", lim.MaxFileBytes, "max bytes read per input file")
	fs.DurationVar(&lim.MaxElapsed, "max-elapsed", lim.MaxElapsed, "wall-clock budget per input file")
	policy := fs.String("limit-policy", "reject", "reject (fail whole input) or truncate (keep first N)")
	return &lim, policy
}

func applyPolicy(lim *ingest.Limits, policy string) error {
	switch policy {
	case "reject":
		lim.Policy = ingest.LimitReject
	case "truncate":
		lim.Policy = ingest.LimitTruncate
	default:
		return fmt.Errorf("unsupported limit-policy %q", policy)
	}
	return nil
}

// detectOrChooseEffect 选择汇总统计的效应字段类型：
// 显式给出（beta/or/z）时严格采用；否则按表头自动识别
// （有 SE 列 => beta；仅有 Z 列 => z-only）。
func detectOrChooseEffect(path, explicit string) (ingest.EffectField, error) {
	switch explicit {
	case "beta":
		return ingest.EffectBeta, nil
	case "or":
		return ingest.EffectOR, nil
	case "z":
		return ingest.EffectZ, nil
	}
	f, err := os.Open(path)
	if err != nil {
		return 0, err
	}
	defer f.Close()
	br := bufio.NewReader(f)
	line, _, err := br.ReadLine()
	if err != nil {
		return 0, fmt.Errorf("%s: cannot read header: %w", path, err)
	}
	has := func(names ...string) bool {
		for _, f := range bytes.FieldsFunc(line, func(r rune) bool { return r == '\t' }) {
			n := strings.ToLower(strings.TrimSpace(string(f)))
			n = strings.NewReplacer("_", "", ".", "", "-", "", " ", "").Replace(n)
			for _, want := range names {
				if n == want {
					return true
				}
			}
		}
		return false
	}
	switch {
	case has("or", "oddsratio"):
		return ingest.EffectOR, nil
	case has("z", "zscore") && !has("se", "sebeta", "standarderror", "stderr"):
		return ingest.EffectZ, nil
	default:
		return ingest.EffectBeta, nil
	}
}

func openSummary(path string, n float64, effect ingest.EffectField, lim ingest.Limits) (*ingest.SummaryTable, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return ingest.ParseSummary(f, ingest.SummaryOptions{
		Source: path, Effect: effect, FallbackN: n, Limits: lim,
	})
}

func openLD(path string, lim ingest.Limits) (*ingest.LDScoreTable, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return ingest.ParseLDScores(f, ingest.LDSOptions{Source: path, Limits: lim})
}

func cmdH2(args []string) error {
	fs := flag.NewFlagSet("h2", flag.ContinueOnError)
	sumPath := fs.String("sum", "", "summary statistics TSV")
	ldPath := fs.String("ld", "", "LD score TSV")
	n := fs.Float64("n", 0, "fallback sample size when N column absent")
	free := fs.Bool("free-intercept", true, "estimate intercept freely")
	boundary := fs.String("boundary", "report", "report|clamp|reject for h2 outside [0,h2max]")
	workers := fs.Int("workers", 1, "shard parallelism (results identical to serial)")
	backendName := fs.String("backend", "", "stdlib|gonum (per-invocation; default does not change)")
	crosscheck := fs.Bool("crosscheck", false, "run stdlib and gonum and require agreement within tolerance")
	lim, policy := limitFlags(fs)
	if err := fs.Parse(args); err != nil {
		return err
	}
	backend, err := parseBackend(*backendName)
	if err != nil {
		return err
	}
	if err := applyPolicy(lim, *policy); err != nil {
		return err
	}
	if *sumPath == "" || *ldPath == "" {
		return fmt.Errorf("-sum and -ld are required")
	}
	sum, err := openSummary(*sumPath, *n, ingest.EffectBeta, *lim)
	if err != nil {
		return err
	}
	ld, err := openLD(*ldPath, *lim)
	if err != nil {
		return err
	}
	fmt.Fprintf(os.Stderr, "usage summary: %s | %s\n", sum.Report.Usage, ld.Report.Usage)
	fmt.Fprintf(os.Stderr, "numerical backend: %s crosscheck=%v\n", backend, *crosscheck)
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
		Backend: backend, CrossCheck: *crosscheck,
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
	parts := strings.Split(rest, ":")
	si.Path = parts[0]
	if len(parts) >= 2 && parts[1] != "" {
		si.SampleIDPath = parts[1]
	}
	if len(parts) >= 3 && parts[2] != "" {
		// 显式效应来源：beta|or|z；不给则按表头自动识别
		// （有 SE 列即 beta/or，仅有 Z 列即 z-only）。
		switch parts[2] {
		case "beta", "or", "z":
			si.Effect = parts[2]
		default:
			return fmt.Errorf("cohort effect kind must be beta|or|z, got %q", parts[2])
		}
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
	backendName := fs.String("backend", "", "stdlib|gonum (per-invocation)")
	crosscheck := fs.Bool("crosscheck", false, "run stdlib and gonum and require agreement")
	audit := fs.Bool("audit", true, "include per-variant cohort audit columns in TSV")
	lim, policy := limitFlags(fs)
	if err := fs.Parse(args); err != nil {
		return err
	}
	backend, err := parseBackend(*backendName)
	if err != nil {
		return err
	}
	if err := applyPolicy(lim, *policy); err != nil {
		return err
	}
	if len(ins) < 2 {
		return fmt.Errorf("at least two -in cohorts required")
	}
	cohorts := make([]stats.Cohort, 0, len(ins))
	for _, si := range ins {
		effect, err := detectOrChooseEffect(si.Path, si.Effect)
		if err != nil {
			return err
		}
		tbl, err := openSummary(si.Path, 0, effect, *lim)
		if err != nil {
			return fmt.Errorf("cohort %s: %w", si.Cohort, err)
		}
		fmt.Fprintf(os.Stderr, "cohort %s effect=%s usage: %s\n", si.Cohort, effect, tbl.Report.Usage)
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
		Backend: backend, CrossCheck: *crosscheck,
	})
	if err != nil {
		return err
	}
	fmt.Println(res.Diag)
	fmt.Printf("backend=%s crosschecked=%v\n", res.Backend, res.CrossChecked)
	for _, e := range res.Overlap {
		fmt.Printf("overlap: %s/%s shared=%d jaccard=%.3f duplicate=%v\n",
			e.CohortA, e.CohortB, e.SharedSamples, e.Jaccard, e.Duplicate)
	}
	if *out != "" {
		return writeMetaTSV(*out, res, *audit)
	}
	nShow := 5
	if len(res.Variants) < nShow {
		nShow = len(res.Variants)
	}
	fmt.Println("variant\tbeta\tse\tp\tncohorts\tQ\tI2\texcluded")
	for i := 0; i < nShow; i++ {
		mv := res.Variants[i]
		id := variantID(mv)
		fmt.Printf("%s\t%.5g\t%.5g\t%.3g\t%d\t%.3f\t%.3f\t%v\n",
			id, mv.Beta, mv.SE, mv.P, mv.NCohorts, mv.Q, mv.I2, mv.Excluded)
	}
	return nil
}

func variantID(mv stats.MetaVariant) string {
	if mv.V != nil {
		return mv.V.ID
	}
	return "<excluded:" + string(mv.ExcludeReason) + ">"
}

func writeMetaTSV(path string, res *stats.MetaResult, auditCols bool) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	header := "variant\tbeta\tse\tci_low\tci_high\tp\tncohorts\tQ\tdf\tpQ\tI2\ttau2\texcluded"
	if auditCols {
		header += "\tincluded_cohorts\tper_cohort_detail\texcluded_detail"
	}
	fmt.Fprintln(f, header)
	for _, mv := range res.Variants {
		fmt.Fprintf(f, "%s\t%.10g\t%.10g\t%.10g\t%.10g\t%.6g\t%d\t%.6g\t%d\t%.6g\t%.6g\t%.6g\t%v",
			variantID(mv), mv.Beta, mv.SE, mv.CI95[0], mv.CI95[1], mv.P, mv.NCohorts,
			mv.Q, mv.Df, mv.PQ, mv.I2, mv.Tau2, mv.Excluded)
		if auditCols {
			fmt.Fprintf(f, "\t%s\t%s\t%s\n",
				strings.Join(mv.IncludedCohorts, ","),
				encodeCohortLines(mv.PerCohort, false),
				encodeCohortLines(mv.ExclusionLines, true))
		} else {
			fmt.Fprintln(f)
		}
	}
	return nil
}

// encodeCohortLines 把每队列审计行编码为 TSV 单元格内的分号分隔文本：
// 队列/来源/权重/方向依据（纳入）或队列/原因（排除）。
func encodeCohortLines(lines []stats.MetaCohortLine, excludedOnly bool) string {
	var parts []string
	for _, l := range lines {
		if excludedOnly {
			parts = append(parts, fmt.Sprintf("%s[reason=%s]", l.Cohort, l.Reason))
			continue
		}
		if l.Excluded {
			continue
		}
		parts = append(parts, fmt.Sprintf(
			"%s[kind=%s weight=%.6g basis=%q sign=%+.0f relation=%s]",
			l.Cohort, l.SourceKind, l.Weight, l.WeightBasis, l.Sign, l.Relation))
	}
	return strings.Join(parts, "; ")
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
	backendName := fs.String("backend", "", "stdlib|gonum (per-invocation)")
	crosscheck := fs.Bool("crosscheck", false, "run stdlib and gonum and require agreement")
	lim, policy := limitFlags(fs)
	if err := fs.Parse(args); err != nil {
		return err
	}
	backend, err := parseBackend(*backendName)
	if err != nil {
		return err
	}
	if err := applyPolicy(lim, *policy); err != nil {
		return err
	}
	if *sumPath == "" || *gtPath == "" || *out == "" {
		return fmt.Errorf("-sum, -gt and -out are required")
	}
	model, err := openSummary(*sumPath, *n, ingest.EffectBeta, *lim)
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
		Source: *gtPath, Encoding: encoding, Limits: *lim})
	f.Close()
	if err != nil {
		return err
	}
	fmt.Fprintf(os.Stderr, "usage summary: %s | %s\n", model.Report.Usage, gt.Report.Usage)
	fmt.Fprintf(os.Stderr, "numerical backend: %s crosscheck=%v\n", backend, *crosscheck)
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
		Align: variant.DefaultAlignConfig(), Backend: backend, CrossCheck: *crosscheck,
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
