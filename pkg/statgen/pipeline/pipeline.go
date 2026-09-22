package pipeline

import (
	"fmt"
	"os"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/genotype"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/ldsc"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/ldscore"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/meta"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/prs"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/sumstats"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/variant"
)

// Config is the full end-to-end run configuration.
type Config struct {
	Workdir      string
	SumstatFiles []string
	LDScoreFile  string
	GenotypeFile string
	CohortNames  []string
	SampleIDs    [][]string // optional de-identified overlap identifiers
	OverlapKeys  []string   // optional provenance keys
	Aligner      *variant.Aligner
	LDConfig     ldsc.Config
	MetaConfig   meta.Config
	PRSConfig    prs.Config
	GenotypeCfg  genotype.ParseConfig
	SumstatCfg   sumstats.ParseConfig
	LDScoreCfg   ldscore.ParseConfig
	Limits       Limits
}

// H2ByCohort is one cohort heritability result.
type H2ByCohort struct {
	Cohort string
	Result *ldsc.Result
}

// Result is the complete output of a successful pipeline run. It is only
// constructed after every stage succeeds, so a failed run leaves no
// partially populated Result with the caller.
type Result struct {
	H2       []H2ByCohort
	Meta     []meta.VariantResult
	Overlaps []meta.Overlap
	Scores   *prs.Scores
	Reports  []variant.Report
}

// Run executes the end-to-end pipeline.
func Run(cfg Config) (*Result, error) {
	budget := NewBudget(cfg.Limits)

	if err := budget.Consume(BudgetCohorts, int64(len(cfg.SumstatFiles))); err != nil {
		return nil, err
	}

	// Stage 1: parse all cohorts (local variables only until complete).
	cohorts := make([]meta.Cohort, 0, len(cfg.SumstatFiles))
	var reports []variant.Report
	for i, path := range cfg.SumstatFiles {
		if err := checkFileSize(path, budget); err != nil {
			return nil, err
		}
		scfg := cfg.SumstatCfg
		scfg.SkipReject = true // pipeline collects rejections in reports
		scfg.HasHeader = true
		if len(cfg.CohortNames) > i && cfg.CohortNames[i] != "" {
			scfg.CohortName = cfg.CohortNames[i]
		}
		if scfg.CohortName == "" {
			scfg.CohortName = cohortDefaultName(i)
		}
		cs, err := sumstats.ParseFile(os.ExpandEnv(path), scfg)
		if err != nil {
			return nil, err // file already closed by ParseFile's defer
		}
		if err := budget.Consume(BudgetVariants, int64(len(cs.Variants))); err != nil {
			return nil, err
		}
		reports = append(reports, cs.Report)
		mc := meta.Cohort{Summary: cs}
		if len(cfg.SampleIDs) > i {
			mc.SampleIDs = cfg.SampleIDs[i]
		}
		if len(cfg.OverlapKeys) > i {
			mc.OverlapKey = cfg.OverlapKeys[i]
		}
		cohorts = append(cohorts, mc)
	}

	// Stage 2: LD scores (optional for meta-only runs).
	var ld *ldscore.Table
	if cfg.LDScoreFile != "" {
		if err := checkFileSize(cfg.LDScoreFile, budget); err != nil {
			return nil, err
		}
		lcfg := cfg.LDScoreCfg
		lcfg.SkipReject = true
		lcfg.HasHeader = true
		if lcfg.Source == "" {
			lcfg.Source = cfg.LDScoreFile
		}
		t, err := ldscore.ParseFile(os.ExpandEnv(cfg.LDScoreFile), lcfg)
		if err != nil {
			return nil, err
		}
		if err := budget.Consume(BudgetVariants, int64(len(t.Scores))); err != nil {
			return nil, err
		}
		ld = t
		reports = append(reports, t.Report)
	}

	// Stage 3: heritability per cohort.
	var h2res []H2ByCohort
	if ld != nil {
		for _, mc := range cohorts {
			lcfg := cfg.LDConfig
			lcfg.Name = mc.Summary.Name
			r, err := ldsc.Estimate(mc.Summary, ld, lcfg)
			if err != nil && !stgerr.Is(err, stgerr.KindOutOfBounds) {
				return nil, err
			}
			// Negative/unidentified h2 is a valid interpretable outcome.
			h2res = append(h2res, H2ByCohort{Cohort: mc.Summary.Name, Result: r})
			reports = append(reports, r.Report)
		}
	}

	// Stage 4: meta-analysis (overlap explicitly handled).
	an, err := meta.Prepare(cohorts, cfg.MetaConfig)
	if err != nil {
		return nil, err
	}
	pooled, err := an.PoolAll(cfg.MetaConfig)
	if err != nil {
		return nil, err
	}

	res := &Result{
		H2:       h2res,
		Meta:     pooled,
		Overlaps: an.Overlaps,
	}

	// Stage 5 (optional): genotypes and PRS from pooled effects.
	if cfg.GenotypeFile != "" {
		if err := checkFileSize(cfg.GenotypeFile, budget); err != nil {
			return nil, err
		}
		gcfg := cfg.GenotypeCfg
		gcfg.SkipReject = true
		if gcfg.Source == "" {
			gcfg.Source = cfg.GenotypeFile
		}
		if !gcfg.HasHeader {
			gcfg.HasHeader = true
			gcfg.HasSampleCol = true
		}
		mat, err := genotype.ParseFile(os.ExpandEnv(cfg.GenotypeFile), gcfg)
		if err != nil {
			return nil, err
		}
		if err := budget.Consume(BudgetVariants, int64(mat.NumVariants())); err != nil {
			return nil, err
		}
		if err := budget.Consume(BudgetSamples, int64(mat.NumSamples())); err != nil {
			return nil, err
		}
		reports = append(reports, mat.Report)

		aligner := cfg.Aligner
		if aligner == nil {
			aligner = variant.DefaultAligner()
		}
		set, err := prs.FromPooled(pooled, mat.Variants, aligner)
		if err != nil {
			return nil, err
		}
		scores, err := set.Score(mat, cfg.PRSConfig)
		if err != nil {
			return nil, err
		}
		res.Scores = scores
	}
	res.Reports = reports
	return res, nil
}

func checkFileSize(path string, b *Budget) error {
	info, err := os.Stat(os.ExpandEnv(path))
	if err != nil {
		return stgerr.Wrap(stgerr.KindMalformedInput, "pipeline", "stat "+path, err)
	}
	if err := b.Consume(BudgetBytes, info.Size()); err != nil {
		return err
	}
	return nil
}

func cohortDefaultName(i int) string {
	return fmt.Sprintf("cohort_%d", i+1)
}
