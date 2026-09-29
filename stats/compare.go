package stats

import (
	"fmt"

	"github.com/lechandonga/statistical-genetics-go/ingest"
	"github.com/lechandonga/statistical-genetics-go/numwrap"
)

// quantile975 按后端取标准正态 97.5% 分位数（错误已在 numwrap 层归类）。
func quantile975(b numwrap.Backend) float64 {
	z, _ := numwrap.NormalCDFInvBackend(0.975, b)
	return z
}

// compareH2Results 比较两个后端的遗传力结论：数值在容差内一致，
// 标志位、NSnps、拒绝明细的原因集合必须完全一致。
func compareH2Results(a, b *H2Result) error {
	pairs := []struct {
		name   string
		v1, v2 float64
	}{
		{"h2_raw", a.H2Raw, b.H2Raw},
		{"h2", a.H2, b.H2},
		{"se", a.SE, b.SE},
		{"ci_low", a.CI95[0], b.CI95[0]},
		{"ci_high", a.CI95[1], b.CI95[1]},
		{"intercept", a.Intercept, b.Intercept},
		{"intercept_se", a.InterceptSE, b.InterceptSE},
	}
	for _, p := range pairs {
		if !floatsAgree(p.v1, p.v2) {
			return fmt.Errorf("h2 %s differs: %.17g vs %.17g", p.name, p.v1, p.v2)
		}
	}
	if a.Flags != b.Flags {
		return fmt.Errorf("h2 flags differ: %s vs %s (judgement must match exactly)", a.Flags, b.Flags)
	}
	if a.NSnps != b.NSnps {
		return fmt.Errorf("h2 NSnps differ: %d vs %d", a.NSnps, b.NSnps)
	}
	if err := compareRejectRows(a.Rejected, b.Rejected); err != nil {
		return fmt.Errorf("h2 rejected rows differ: %w", err)
	}
	return nil
}

func compareRejectRows(a, b []ingest.RowRef) error {
	if len(a) != len(b) {
		return fmt.Errorf("rejected count %d != %d", len(a), len(b))
	}
	for i := range a {
		if a[i].Line != b[i].Line || a[i].ID != b[i].ID ||
			a[i].Reason != b[i].Reason {
			return fmt.Errorf("row %d differs: %+v vs %+v", i, a[i], b[i])
		}
	}
	return nil
}

// compareMetaResults 比较两个后端的元分析结论。
// 后端差异仅可能来自分位数（CI/P 值）；所有判定字段必须逐字段一致。
func compareMetaResults(a, b *MetaResult) error {
	if len(a.Variants) != len(b.Variants) {
		return fmt.Errorf("variant count %d != %d", len(a.Variants), len(b.Variants))
	}
	for i := range a.Variants {
		va, vb := a.Variants[i], b.Variants[i]
		nums := []struct {
			name   string
			v1, v2 float64
		}{
			{"beta", va.Beta, vb.Beta}, {"se", va.SE, vb.SE},
			{"ci_low", va.CI95[0], vb.CI95[0]}, {"ci_high", va.CI95[1], vb.CI95[1]},
			{"p", va.P, vb.P}, {"q", va.Q, vb.Q}, {"pq", va.PQ, vb.PQ},
			{"i2", va.I2, vb.I2}, {"tau2", va.Tau2, vb.Tau2},
			{"wfixed", va.WeightFixed, vb.WeightFixed},
			{"wrandom", va.WeightRandom, vb.WeightRandom},
		}
		for _, p := range nums {
			if !floatsAgree(p.v1, p.v2) {
				return fmt.Errorf("variant %d %s differs: %.17g vs %.17g", i, p.name, p.v1, p.v2)
			}
		}
		if va.Excluded != vb.Excluded || va.ExcludeReason != vb.ExcludeReason ||
			va.NCohorts != vb.NCohorts || va.Df != vb.Df {
			return fmt.Errorf("variant %d judgement differs: excluded=%v/%v reason=%s/%s n=%d/%d",
				i, va.Excluded, vb.Excluded, va.ExcludeReason, vb.ExcludeReason,
				va.NCohorts, vb.NCohorts)
		}
		if err := compareCohortLines(va.PerCohort, vb.PerCohort); err != nil {
			return fmt.Errorf("variant %d audit differs: %w", i, err)
		}
	}
	return nil
}

func compareCohortLines(a, b []MetaCohortLine) error {
	if len(a) != len(b) {
		return fmt.Errorf("per-cohort line count %d != %d", len(a), len(b))
	}
	for i := range a {
		x, y := a[i], b[i]
		if x.Cohort != y.Cohort || x.Excluded != y.Excluded ||
			x.Reason != y.Reason || x.Sign != y.Sign || x.Relation != y.Relation {
			return fmt.Errorf("line %d judgement differs: %+v vs %+v", i, x, y)
		}
		// 参与合并的队列：权重/对齐效应也必须在容差内一致。
		if !x.Excluded {
			if !floatsAgree(x.Weight, y.Weight) {
				return fmt.Errorf("line %d (%s) weight differs: %.17g vs %.17g",
					i, x.Cohort, x.Weight, y.Weight)
			}
			if !floatsAgree(x.BetaAligned, y.BetaAligned) {
				return fmt.Errorf("line %d (%s) aligned beta differs: %.17g vs %.17g",
					i, x.Cohort, x.BetaAligned, y.BetaAligned)
			}
		}
	}
	return nil
}

// comparePRSResults 比较两个后端的评分结论；PRS 中后端仅影响 std
// 尺度的均值/方差与评估回归，raw/perallele 分数应逐位相同。
func comparePRSResults(a, b *PRSResult) error {
	if len(a.Scores) != len(b.Scores) || len(a.NActive) != len(b.NActive) {
		return fmt.Errorf("sample count differs")
	}
	for i := range a.Scores {
		if !floatsAgree(a.Scores[i], b.Scores[i]) {
			return fmt.Errorf("score sample %d differs: %.17g vs %.17g",
				i, a.Scores[i], b.Scores[i])
		}
		if a.NActive[i] != b.NActive[i] {
			return fmt.Errorf("active count sample %d differs: %d vs %d",
				i, a.NActive[i], b.NActive[i])
		}
	}
	if a.Scale != b.Scale || len(a.Audit) != len(b.Audit) ||
		len(a.MissingModelSites) != len(b.MissingModelSites) {
		return fmt.Errorf("scale/audit/site-set judgement differs")
	}
	for i := range a.Audit {
		if a.Audit[i] != b.Audit[i] {
			return fmt.Errorf("audit line %d differs: %+v vs %+v", i, a.Audit[i], b.Audit[i])
		}
	}
	return nil
}
