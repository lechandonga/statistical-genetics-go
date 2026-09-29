package main

import (
	"fmt"
	"io"
	"os"

	"github.com/lechandonga/statistical-genetics-go/ingest"
	"github.com/lechandonga/statistical-genetics-go/stats"
)

// stringFlag 是可重复的字符串 flag（用于 -z 标记多个 Z-only 队列）。
type stringFlag []string

func (s *stringFlag) String() string {
	if s == nil {
		return ""
	}
	return fmt.Sprint([]string(*s))
}

func (s *stringFlag) Set(v string) error {
	*s = append(*s, v)
	return nil
}

func effectName(e ingest.EffectField) string {
	switch e {
	case ingest.EffectBeta:
		return "beta+se"
	case ingest.EffectOR:
		return "or"
	case ingest.EffectZ:
		return "z-only"
	default:
		return "unknown"
	}
}

// writeMetaAuditTSV 输出每个变体 × 每个队列的完整审计：
// 是否参与、权重、对齐方向、效应口径（beta_se / z_n_eaf）及换算依据、
// 排除原因。任何“参与/排除”都能在该文件中逐行核对。
func writeMetaAuditTSV(path string, res *stats.MetaResult) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	return WriteMetaAudit(f, res)
}

// WriteMetaAudit 便于测试直接写入任意 io.Writer。
func WriteMetaAudit(w io.Writer, res *stats.MetaResult) error {
	if _, err := fmt.Fprintln(w,
		"variant\tcohort\tincluded\tbasis\tbeta_aligned\tse\tweight\tsign\trelation\tN\tZ\tEAF\treason\tbasis_detail"); err != nil {
		return err
	}
	for _, mv := range res.Variants {
		vid := ""
		if mv.V != nil {
			vid = mv.V.ID
		}
		for _, l := range mv.PerCohort {
			reason := ""
			if l.Excluded {
				reason = string(l.Reason)
			}
			if _, err := fmt.Fprintf(w,
				"%s\t%s\t%v\t%s\t%.10g\t%.10g\t%.10g\t%g\t%s\t%.0f\t%.6g\t%.6g\t%s\t%s\n",
				vid, l.Cohort, !l.Excluded, l.Basis, l.BetaAligned, l.SE,
				l.Weight, l.Sign, l.Relation.String(), l.N, l.Z, l.EAF,
				reason, l.BasisDetail); err != nil {
				return err
			}
		}
	}
	return nil
}
