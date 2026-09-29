package main

import (
	"flag"
	"time"

	"github.com/lechandonga/statistical-genetics-go/ingest"
	"github.com/lechandonga/statistical-genetics-go/numwrap"
)

// limitFlags 把 ingest.Limits 的各项开放到命令行；不同子命令用不同
// 前缀注册（-sum-max-variants / -ld-max-variants / -gt-max-variants），
// 实现“按输入类型分别配置”。
type limitFlags struct {
	maxVariants  int
	maxSamples   int
	maxLineBytes int64
	maxFileBytes int64
	maxElapsed   time.Duration
	policy       string
}

func registerLimitFlags(fs *flag.FlagSet, prefix string) *limitFlags {
	d := ingest.DefaultLimits()
	lf := &limitFlags{
		maxVariants:  d.MaxVariants,
		maxSamples:   d.MaxSamples,
		maxLineBytes: d.MaxLineBytes,
		maxFileBytes: d.MaxFileBytes,
		maxElapsed:   d.MaxElapsed,
		policy:       "reject",
	}
	fs.IntVar(&lf.maxVariants, prefix+"-max-variants", lf.maxVariants,
		"max data rows (variants) for this input type; <=0 uses default")
	fs.IntVar(&lf.maxSamples, prefix+"-max-samples", lf.maxSamples,
		"max sample columns (genotype inputs); <=0 uses default")
	fs.Int64Var(&lf.maxLineBytes, prefix+"-max-line-bytes", lf.maxLineBytes,
		"max bytes per single line; oversized lines are always rejected")
	fs.Int64Var(&lf.maxFileBytes, prefix+"-max-file-bytes", lf.maxFileBytes,
		"max total file bytes")
	fs.DurationVar(&lf.maxElapsed, prefix+"-max-elapsed", lf.maxElapsed,
		"wall-clock budget for parsing this input")
	fs.StringVar(&lf.policy, prefix+"-limit-policy", "reject",
		"reject (fail, no partial output) or truncate (cut at limit and report)")
	return lf
}

func (lf *limitFlags) limits() (ingest.Limits, error) {
	l := ingest.Limits{
		MaxVariants:  lf.maxVariants,
		MaxSamples:   lf.maxSamples,
		MaxLineBytes: lf.maxLineBytes,
		MaxFileBytes: lf.maxFileBytes,
		MaxElapsed:   lf.maxElapsed,
	}
	switch lf.policy {
	case "reject":
		l.Policy = ingest.LimitReject
	case "truncate":
		l.Policy = ingest.LimitTruncate
	default:
		return l, errUnknownLimitPolicy(lf.policy)
	}
	return l, nil
}

func errUnknownLimitPolicy(p string) error {
	return &userError{msg: "unknown limit policy " + p + " (want reject|truncate)"}
}

type userError struct{ msg string }

func (e *userError) Error() string { return e.msg }

// backendFlag 注册数值后端选择；显式拒绝未知后端，不退回默认实现。
func backendFlag(fs *flag.FlagSet) *string {
	return fs.String("backend", string(numwrap.BackendStd),
		"numerical backend: stdlib|gonum|crosscheck (unknown value is rejected)")
}

func resolveBackend(s string) (numwrap.Backend, error) {
	b, err := numwrap.ResolveBackend(s)
	if err != nil {
		return "", &userError{msg: err.Error()}
	}
	return b, nil
}
