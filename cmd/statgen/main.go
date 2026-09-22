// Command statgen runs the local statistical-genetics pipeline:
// parses summary statistics, LD scores and a genotype matrix, performs
// LD-score heritability estimation, multi-cohort meta-analysis and
// individual risk scoring, and emits deterministic JSON plus human-readable
// rejection/input summaries.
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"strings"
	"time"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/ldsc"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/meta"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/pipeline"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/prs"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/variant"
)

func main() {
	var (
		sumstats = flag.String("sumstats", "", "comma-separated cohort summary-statistics files")
		names    = flag.String("names", "", "comma-separated cohort names (aligned with -sumstats)")
		ldfile   = flag.String("ldscores", "", "LD-score file")
		geno     = flag.String("genotypes", "", "genotype dosage matrix file")
		model    = flag.String("meta-model", "fixed_iv", "fixed_iv|fixed_z|random_dl")
		overlap  = flag.String("overlap-policy", "flag", "reject|flag|deduplicate")
		amb      = flag.String("ambiguous", "reject", "reject|mark|force")
		neg      = flag.String("negative-h2", "report", "report|truncate_zero")
		scale    = flag.String("score-scale", "raw", "raw|mean_center|standardized")
		workers  = flag.Int("workers", 1, "parallel workers (output stays serial-identical)")
		maxVar   = flag.Int("max-variants", 0, "variant budget override")
		maxSamp  = flag.Int("max-samples", 0, "sample budget override")
		timeout  = flag.Duration("timeout", 0, "compute time budget (e.g. 5m)")
		out      = flag.String("out", "", "write JSON result to this path (default: stdout summary only)")
	)
	flag.Parse()

	if *sumstats == "" {
		fmt.Fprintln(os.Stderr, "statgen: -sumstats is required")
		flag.Usage()
		os.Exit(2)
	}

	cfg := pipeline.Config{
		SumstatFiles: splitCSV(*sumstats),
		CohortNames:  splitCSV(*names),
		LDScoreFile:  *ldfile,
		GenotypeFile: *geno,
		Aligner:      alignerFromFlag(*amb),
		LDConfig: ldsc.Config{
			NegPolicy: ldsc.NegativePolicy(*neg), Workers: *workers,
		},
		MetaConfig: meta.Config{
			Model:         meta.Model(*model),
			OverlapPolicy: meta.OverlapPolicy(*overlap),
			MinCohorts:    1, Workers: *workers,
		},
		PRSConfig: prs.Config{
			ScoreScale:   prs.Scale(*scale),
			MissingLocus: prs.MissingLocusSkip,
			Workers:      *workers,
		},
		Limits: pipeline.Limits{
			MaxVariants: *maxVar, MaxSamples: *maxSamp, Timeout: *timeout,
		},
	}

	res, err := pipeline.Run(cfg)
	if err != nil {
		se := stgerr.AsError(err)
		if se != nil {
			fmt.Fprintf(os.Stderr, "statgen: rejected: kind=%s source=%s reason=%s\n",
				se.Kind, se.Source, se.Reason)
		} else {
			fmt.Fprintf(os.Stderr, "statgen: %v\n", err)
		}
		os.Exit(1)
	}

	for _, r := range res.Reports {
		fmt.Fprintln(os.Stderr, "INPUT ", r.Summary())
	}
	for _, h := range res.H2 {
		r := h.Result
		fmt.Fprintf(os.Stderr, "H2[%s] h2=%.4f se=%.4f ci=[%.4f,%.4f] intercept=%.3f status=%s n=%d\n",
			h.Cohort, r.H2, r.SE, r.CI[0], r.CI[1], r.Intercept, r.Status, r.NVariants)
	}
	for _, ov := range res.Overlaps {
		fmt.Fprintf(os.Stderr, "OVERLAP %s~%s kind=%s shared=%d\n", ov.A, ov.B, ov.Kind, ov.Shared)
	}
	if res.Scores != nil {
		fmt.Fprintf(os.Stderr, "PRS %s\n", res.Scores.Report.String())
	}
	if *out != "" {
		f, ferr := os.Create(*out)
		if ferr != nil {
			fatal(ferr)
		}
		enc := json.NewEncoder(f)
		enc.SetIndent("", "  ")
		if jerr := enc.Encode(jsonResult(res)); jerr != nil {
			_ = f.Close()
			fatal(jerr)
		}
		_ = f.Close()
		fmt.Fprintf(os.Stderr, "wrote %s at %s\n", *out, time.Now().Format(time.RFC3339))
	}
}

func splitCSV(s string) []string {
	if s == "" {
		return nil
	}
	parts := strings.Split(s, ",")
	for i := range parts {
		parts[i] = strings.TrimSpace(parts[i])
	}
	return parts
}

func alignerFromFlag(s string) *variant.Aligner {
	a := variant.DefaultAligner()
	switch s {
	case "reject":
		a.Ambiguous = variant.AmbiguousReject
	case "mark":
		a.Ambiguous = variant.AmbiguousMark
	case "force":
		a.Ambiguous = variant.AmbiguousFlipForce
	default:
		fatalf("unsupported -ambiguous policy %q", s)
	}
	return a
}

type jsonOut struct {
	H2       []jsonH2       `json:"heritability"`
	Overlaps []meta.Overlap `json:"overlaps"`
	MetaN    int            `json:"meta_variants"`
	PRS      *jsonPRS       `json:"prs,omitempty"`
}
type jsonH2 struct {
	Cohort    string     `json:"cohort"`
	H2        float64    `json:"h2"`
	SE        float64    `json:"se"`
	CI        [2]float64 `json:"ci"`
	Intercept float64    `json:"intercept"`
	Status    string     `json:"status"`
	NVariants int        `json:"n_variants"`
}
type jsonPRS struct {
	N        int       `json:"n_samples"`
	UsedLoci int       `json:"used_loci"`
	Values   []float64 `json:"values"`
}

func jsonResult(r *pipeline.Result) jsonOut {
	o := jsonOut{MetaN: len(r.Meta), Overlaps: r.Overlaps}
	for _, h := range r.H2 {
		o.H2 = append(o.H2, jsonH2{
			Cohort: h.Cohort, H2: h.Result.H2, SE: h.Result.SE, CI: h.Result.CI,
			Intercept: h.Result.Intercept, Status: string(h.Result.Status),
			NVariants: h.Result.NVariants,
		})
	}
	if r.Scores != nil {
		o.PRS = &jsonPRS{N: len(r.Scores.Values), UsedLoci: r.Scores.UsedLoci,
			Values: r.Scores.Values}
	}
	return o
}

func fatal(err error) {
	fmt.Fprintln(os.Stderr, "statgen:", err)
	os.Exit(1)
}
func fatalf(format string, args ...any) {
	fmt.Fprintf(os.Stderr, "statgen: "+format+"\n", args...)
	os.Exit(2)
}
