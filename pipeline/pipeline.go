// Package pipeline 编排解析、对齐、遗传力估计、元分析与评分的完整本地流程。
//
// 全局保证：
//   - 相同输入与配置，任意次数、任意并行度下结果逐位一致；
//   - 共享状态并发安全，对外只暴露完整快照（无半更新结果）；
//   - 超限/取消时不残留部分统计状态；
//   - 所有错误均为 variant / numwrap / ingest / stats 的可区分错误类型。
package pipeline

import (
	"context"
	"fmt"
	"os"
	"time"

	"github.com/lechandonga/statistical-genetics-go/ingest"
	"github.com/lechandonga/statistical-genetics-go/numwrap"
	"github.com/lechandonga/statistical-genetics-go/stats"
	"github.com/lechandonga/statistical-genetics-go/variant"
)

// Config 是流程级配置。
type Config struct {
	Limits   ingest.Limits
	Align    variant.AlignConfig
	H2       stats.H2Options
	Meta     stats.MetaOptions
	PRS      stats.PRSOptions
	Workers  int           // 全局并行度上限（各阶段取 Min）
	Deadline time.Duration // 整体墙钟上限；0 表示不限
	// EffectField / FallbackN：汇总统计解析的共享设置。
	Effect    ingest.EffectField
	FallbackN float64
	// Backend 选择本次流程的数值后端；CrossCheck=true 时两侧各跑一遍
	// 并要求数值在声明容差内、离散判定完全一致，否则整体失败。
	Backend    numwrap.Backend
	CrossCheck bool
	// EnableH2/Meta/PRS 控制执行阶段（至少启用一个）。
	EnableH2, EnableMeta, EnablePRS bool
	// PRS 基因型编码（默认 0/1/2 剂量）。
	GTEncoding ingest.GTEncoding
}

// DefaultConfig 返回带保守上限的默认配置（三阶段全开，严格模式）。
func DefaultConfig() Config {
	align := variant.DefaultAlignConfig()
	return Config{
		Limits:     ingest.DefaultLimits(),
		Align:      align,
		H2:         stats.H2Options{Align: align, MinVariants: 20, Boundary: stats.H2ReportRaw, Workers: 1},
		Meta:       stats.MetaOptions{Align: align, Model: stats.MetaFixed, Overlap: stats.OverlapFlag, MinCohorts: 2, Workers: 1},
		PRS:        stats.PRSOptions{Align: align, ScoreScale: "raw", Workers: 1},
		Workers:    1,
		Effect:     ingest.EffectBeta,
		EnableH2:   true,
		EnableMeta: true,
		EnablePRS:  true,
		GTEncoding: ingest.GTDosage012,
	}
}

// Input 声明本地输入文件。
type Input struct {
	Summaries []SummaryInput
	LDScores  string // 路径（EnableH2 时必填）
	Genotypes string // 路径（EnablePRS 时必填）
}

// SummaryInput 声明单个队列的汇总统计文件。
type SummaryInput struct {
	Cohort string
	Path   string
	N      int
	// SampleIDPath 可选：每行一个样本 ID，用于队列重叠硬证据。
	SampleIDPath string
	// Effect 可选：beta|or|z；空 => 按 cfg.Effect。
	Effect string
}

// Result 是完整流程输出（成功时为完整快照，失败时为零值 + 错误）。
type Result struct {
	H2      *stats.H2Result
	Meta    *stats.MetaResult
	PRS     *stats.PRSResult
	Reports []*ingest.ParseReport
}

type pipelineError struct {
	stage string
	err   error
}

func (e *pipelineError) Error() string {
	return fmt.Sprintf("pipeline stage %s failed: %v", e.stage, e.err)
}
func (e *pipelineError) Unwrap() error { return e.err }

// Run 执行完整流程。
//
// 失败语义：任一阶段失败立即返回错误与零值 Result，已打开的文件句柄
// 通过 defer 全部关闭，不返回任何部分统计状态。
func Run(ctx context.Context, in Input, cfg Config) (*Result, error) {
	if cfg.Deadline > 0 {
		var cancel context.CancelFunc
		ctx, cancel = context.WithTimeout(ctx, cfg.Deadline)
		defer cancel()
	}
	if !cfg.EnableH2 && !cfg.EnableMeta && !cfg.EnablePRS {
		return nil, &pipelineError{stage: "config", err: fmt.Errorf("at least one stage must be enabled")}
	}
	if len(in.Summaries) == 0 {
		return nil, &pipelineError{stage: "config", err: fmt.Errorf("at least one summary file required")}
	}

	// 阶段 1：解析全部汇总统计。任何解析失败 => 整体失败（严格模式默认）。
	parsed := make([]*ingest.SummaryTable, len(in.Summaries))
	cohorts := make([]stats.Cohort, 0, len(in.Summaries))
	reports := []*ingest.ParseReport{}
	for i, si := range in.Summaries {
		effect := cfg.Effect
		switch si.Effect {
		case "":
		case "beta":
			effect = ingest.EffectBeta
		case "or":
			effect = ingest.EffectOR
		case "z":
			effect = ingest.EffectZ
		default:
			return nil, &pipelineError{stage: "config",
				err: fmt.Errorf("cohort %s: effect kind must be beta|or|z, got %q", si.Cohort, si.Effect)}
		}
		tbl, err := parseFile(si.Path, func(f *os.File) (*ingest.SummaryTable, error) {
			return ingest.ParseSummary(f, ingest.SummaryOptions{
				Source:    si.Path,
				Effect:    effect,
				FallbackN: float64(si.N),
				Limits:    cfg.Limits,
			})
		})
		if err != nil {
			return nil, &pipelineError{stage: "parse-summary", err: err}
		}
		parsed[i] = tbl
		rep := tbl.Report
		reports = append(reports, &rep)

		var ids []string
		if si.SampleIDPath != "" {
			idBytes, err := os.ReadFile(si.SampleIDPath)
			if err != nil {
				return nil, &pipelineError{stage: "read-sample-ids", err: err}
			}
			for _, line := range splitLines(string(idBytes)) {
				if line != "" {
					ids = append(ids, line)
				}
			}
		}
		cohorts = append(cohorts, stats.Cohort{
			Name:      si.Cohort,
			Samples:   si.N,
			SampleIDs: ids,
			Table:     tbl,
		})
	}

	// 阶段 2（可选）：LD 分数 + H2
	res := &Result{Reports: reports}
	if cfg.EnableH2 {
		if in.LDScores == "" {
			return nil, &pipelineError{stage: "config", err: fmt.Errorf("LDScores path required for h2 stage")}
		}
		ld, err := parseFile(in.LDScores, func(f *os.File) (*ingest.LDScoreTable, error) {
			return ingest.ParseLDScores(f, ingest.LDSOptions{
				Source: in.LDScores, Limits: cfg.Limits,
			})
		})
		if err != nil {
			return nil, &pipelineError{stage: "parse-ld", err: err}
		}
		rep := ld.Report
		reports = append(reports, &rep)
		h2opts := cfg.H2
		h2opts.Align = cfg.Align
		h2opts.Backend = cfg.Backend
		h2opts.CrossCheck = cfg.CrossCheck
		if cfg.Workers > 0 {
			h2opts.Workers = cfg.Workers
		}
		h2, err := stats.EstimateH2(ctx, parsed[0], ld, h2opts)
		if err != nil {
			return nil, &pipelineError{stage: "h2", err: err}
		}
		res.H2 = h2
	}

	// 阶段 3（可选）：多队列元分析
	if cfg.EnableMeta {
		if len(cohorts) < cfg.Meta.MinCohorts {
			return nil, &pipelineError{stage: "meta",
				err: fmt.Errorf("need >= %d cohorts, got %d", cfg.Meta.MinCohorts, len(cohorts))}
		}
		mopts := cfg.Meta
		mopts.Align = cfg.Align
		mopts.Backend = cfg.Backend
		mopts.CrossCheck = cfg.CrossCheck
		if cfg.Workers > 0 {
			mopts.Workers = cfg.Workers
		}
		meta, err := stats.RunMeta(ctx, cohorts, mopts)
		if err != nil {
			return nil, &pipelineError{stage: "meta", err: err}
		}
		res.Meta = meta
	}

	// 阶段 4（可选）：PRS
	if cfg.EnablePRS {
		if in.Genotypes == "" {
			return nil, &pipelineError{stage: "config", err: fmt.Errorf("Genotypes path required for prs stage")}
		}
		gt, err := parseFile(in.Genotypes, func(f *os.File) (*ingest.GenotypeMatrix, error) {
			return ingest.ParseGenotype(f, ingest.GenotypeOptions{
				Source: in.Genotypes, Encoding: cfg.GTEncoding, Limits: cfg.Limits,
			})
		})
		if err != nil {
			return nil, &pipelineError{stage: "parse-genotype", err: err}
		}
		rep := gt.Report
		reports = append(reports, &rep)
		popts := cfg.PRS
		popts.Align = cfg.Align
		popts.Backend = cfg.Backend
		popts.CrossCheck = cfg.CrossCheck
		if cfg.Workers > 0 {
			popts.Workers = cfg.Workers
		}
		prs, err := stats.BuildPRS(ctx, parsed[0], gt, popts)
		if err != nil {
			return nil, &pipelineError{stage: "prs", err: err}
		}
		res.PRS = prs
	}
	res.Reports = reports
	return res, nil
}

type parsedTable interface{}

func parseFile[T any](path string, parse func(*os.File) (T, error)) (T, error) {
	var zero T
	f, err := os.Open(path)
	if err != nil {
		return zero, err
	}
	defer f.Close()
	return parse(f)
}

func splitLines(s string) []string {
	var out []string
	cur := ""
	for _, r := range s {
		if r == '\n' || r == '\r' {
			if cur != "" {
				out = append(out, cur)
				cur = ""
			}
			continue
		}
		cur += string(r)
	}
	if cur != "" {
		out = append(out, cur)
	}
	return out
}
