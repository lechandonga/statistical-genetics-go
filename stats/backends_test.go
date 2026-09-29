package stats

import (
	"context"
	"math"
	"testing"

	"github.com/lechandonga/statistical-genetics-go/ingest"
	"github.com/lechandonga/statistical-genetics-go/numwrap"
	"github.com/lechandonga/statistical-genetics-go/variant"
)

// TestMetaCrossBackendsAgree：完整 (beta,se) 队列上，两个后端的
// 合并数值在容差内、纳入/排除判定完全一致；同后端逐位一致。
func TestMetaCrossBackendsAgree(t *testing.T) {
	var recs1, recs2 []ingest.SummaryRecord
	for i := 0; i < 120; i++ {
		recs1 = append(recs1, metaRec(idFor(i), "A", "G", 0.10+0.001*float64(i%3), 0.10, 100))
		recs2 = append(recs2, metaRec(idFor(i), "A", "G", 0.14-0.001*float64(i%4), 0.08, 120))
	}
	c1 := Cohort{Name: "c1", Table: &ingest.SummaryTable{Records: recs1}, SampleIDs: []string{"a"}}
	c2 := Cohort{Name: "c2", Table: &ingest.SummaryTable{Records: recs2}, SampleIDs: []string{"b"}}

	for _, model := range []MetaModel{MetaFixed, MetaRandomDL} {
		base := MetaOptions{Model: model, MinCohorts: 2, Backend: numwrap.BackendStd}
		rStd, err := RunMeta(context.Background(), []Cohort{c1, c2}, base)
		if err != nil {
			t.Fatal(err)
		}
		// 同后端重复跑逐位一致
		rStd2, _ := RunMeta(context.Background(), []Cohort{c1, c2}, base)
		for i := range rStd.Variants {
			if rStd.Variants[i].Beta != rStd2.Variants[i].Beta ||
				rStd.Variants[i].SE != rStd2.Variants[i].SE {
				t.Fatalf("model=%v stdlib not bit-reproducible at %d", model, i)
			}
		}
		base.Backend = numwrap.BackendGonum
		rGo, err := RunMeta(context.Background(), []Cohort{c1, c2}, base)
		if err != nil {
			t.Fatal(err)
		}
		// CrossCheck 本身也必须通过（它返回 mismatch 即硬失败）
		cc := base
		cc.CrossCheck = true
		if _, err := RunMeta(context.Background(), []Cohort{c1, c2}, cc); err != nil {
			t.Fatalf("model=%v crosscheck: %v", model, err)
		}
		if err := crossCheckMeta(rStd.Variants, rGo.Variants); err != nil {
			t.Fatalf("model=%v: %v", model, err)
		}
	}
	t.Log("[judgement] stdlib vs gonum meta: values within tolerance, decisions identical")
}

// TestMetaMixedZOnlyCohorts：一个完整队列 + 一个 Z-only 队列。
// Z-only 队列必须真正参与（有正权重、出现在纳入审计里），
// 缺 N / Z=0 的行按可区分原因排除。
func TestMetaMixedZOnlyCohorts(t *testing.T) {
	zv, _ := variant.New("RS1", "1", 100, "A", "G", math.NaN(), variant.StrandUnknown)
	full := Cohort{Name: "full", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{
		metaRec("RS1", "A", "G", 0.20, 0.10, 100),
	}}, SampleIDs: []string{"a"}}

	// 正常 Z-only：N=100, Z=2
	zok := ingest.SummaryRecord{V: zv, Beta: 0, SE: math.NaN(), Z: 2.0, N: 100, HasZ: true}
	zc := Cohort{Name: "zcohort", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{zok}}}
	res, err := RunMeta(context.Background(), []Cohort{full, zc},
		MetaOptions{Model: MetaFixed, MinCohorts: 2})
	if err != nil {
		t.Fatal(err)
	}
	mv := res.Variants[0]
	if mv.Excluded {
		t.Fatalf("mixed variant must participate: reason=%s", mv.ExcludeReason)
	}
	if mv.NCohorts != 2 {
		t.Fatalf("ncohorts=%d want 2", mv.NCohorts)
	}
	if len(mv.IncludedCohorts) != 2 || mv.IncludedCohorts[0] != "full" ||
		mv.IncludedCohorts[1] != "zcohort" {
		t.Fatalf("included audit wrong: %v", mv.IncludedCohorts)
	}
	var zline *MetaCohortLine
	for i := range mv.PerCohort {
		l := &mv.PerCohort[i]
		if l.Cohort == "zcohort" {
			zline = l
		}
	}
	if zline == nil || zline.SourceKind != "z_only" || !(zline.Weight > 0) ||
		zline.SE <= 0 || zline.WeightBasis == "" {
		t.Fatalf("z-only audit wrong: %+v", zline)
	}
	// 换算值必须有限
	if math.IsNaN(mv.Beta) || math.IsNaN(mv.SE) || math.IsNaN(mv.P) ||
		math.IsInf(mv.Beta, 0) {
		t.Fatal("combined values must be finite")
	}
	t.Logf("[input] full(beta/se) + z-only(Z=2,N=100) [judgement] pooled beta=%.4f se=%.4f; %s weight=%.3f basis=%s",
		mv.Beta, mv.SE, zline.Cohort, zline.Weight, zline.WeightBasis)

	// 缺 N：按 missing_sample_size 排除，变体不足两队 => 整体排除且原因可区分
	znN := zok
	znN.N = 0
	zcN := Cohort{Name: "zn", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{znN}}}
	res2, err := RunMeta(context.Background(), []Cohort{full, zcN},
		MetaOptions{Model: MetaFixed, MinCohorts: 2})
	if err != nil {
		t.Fatal(err)
	}
	mv2 := res2.Variants[0]
	if !mv2.Excluded || mv2.ExcludeReason != "insufficient_cohorts" {
		t.Fatalf("missing-N z cohort must drive exclusion, got excluded=%v reason=%s",
			mv2.Excluded, mv2.ExcludeReason)
	}
	if got := mv2.PerCohort[1]; !got.Excluded || got.Reason != "missing_sample_size" {
		t.Fatalf("missing-N reason wrong: %+v", got)
	}
	t.Logf("[input] z-only without N [judgement] line reason=%s variant reason=%s",
		mv2.PerCohort[1].Reason, mv2.ExcludeReason)

	// Z=0：方向不可定（undirected），可区分原因
	z0 := zok
	z0.Z = 0
	zc0 := Cohort{Name: "z0", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{z0}}}
	res3, _ := RunMeta(context.Background(), []Cohort{full, zc0},
		MetaOptions{Model: MetaFixed, MinCohorts: 2})
	if got := res3.Variants[0].PerCohort[1]; !got.Excluded ||
		got.Reason != "z_undirected" {
		t.Fatalf("z=0 reason wrong: %+v", got)
	}
	t.Logf("[input] z-only Z=0 [judgement] excluded reason=z_undirected")
}

// TestMetaMixedCrossCheckDecisionsIdentical：混合来源队列在两个后端下
// 参与/排除/权重口径的离散结论必须完全一致。
func TestMetaMixedCrossCheckDecisionsIdentical(t *testing.T) {
	zv, _ := variant.New("RS1", "1", 100, "A", "G", math.NaN(), variant.StrandUnknown)
	full := Cohort{Name: "full", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{
		metaRec("RS1", "A", "G", 0.20, 0.10, 100),
	}}}
	zc := Cohort{Name: "z", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{
		{V: zv, SE: math.NaN(), Z: -2.0, N: 100, HasZ: true},
	}}}
	op := MetaOptions{MinCohorts: 2, CrossCheck: true}
	res, err := RunMeta(context.Background(), []Cohort{full, zc}, op)
	if err != nil {
		t.Fatalf("mixed crosscheck must agree: %v", err)
	}
	l := res.Variants[0].PerCohort[1]
	if l.BetaAligned >= 0 || l.Weight <= 0 {
		t.Fatalf("negative z must yield negative aligned beta with positive weight: %+v", l)
	}
	t.Logf("[judgement] mixed-source crosscheck passed; z betaAligned=%.4f weight=%.4f",
		l.BetaAligned, l.Weight)
}

// TestH2CrossBackendsAgree：h2 在两个后端间一致，flags/拒绝原因一致。
func TestH2CrossBackendsAgree(t *testing.T) {
	sum, ld := synthH2Tables(t, 200, 4000, 0.3, 5, true)
	r1, err := EstimateH2(context.Background(), sum, ld, H2Options{
		NBar: 200, Boundary: H2ReportRaw, Backend: numwrap.BackendStd})
	if err != nil {
		t.Fatal(err)
	}
	r2, err := EstimateH2(context.Background(), sum, ld, H2Options{
		NBar: 200, Boundary: H2ReportRaw, Backend: numwrap.BackendGonum})
	if err != nil {
		t.Fatal(err)
	}
	if err := crossCheckH2(r1, r2); err != nil {
		t.Fatal(err)
	}
	if _, err := EstimateH2(context.Background(), sum, ld, H2Options{
		NBar: 200, Boundary: H2ReportRaw, Backend: numwrap.BackendStd, CrossCheck: true}); err != nil {
		t.Fatal(err)
	}
	t.Logf("[judgement] h2 std=%.6f gonum=%.6g flags agree", r1.H2Raw, r2.H2Raw)
}

// TestPRSCrossBackendsAgree：评分与离散审计两后端一致。
func TestPRSCrossBackendsAgree(t *testing.T) {
	samples := []string{"s1", "s2", "s3", "s4"}
	v1, _ := variant.New("RS1", "1", 1, "A", "G", math.NaN(), variant.StrandUnknown)
	v2, _ := variant.New("RS2", "1", 2, "C", "T", math.NaN(), variant.StrandUnknown)
	gt := buildGT(samples,
		[]*variant.Variant{v1, v2},
		[][]float64{{2, 0, 1, 2}, {1, 2, 0, 1}},
		[][]bool{{false, false, false, false}, {false, false, true, false}})
	model := &ingest.SummaryTable{Records: []ingest.SummaryRecord{
		prsModelRec("RS1", "A", "G", 0.5),
		prsModelRec("RS2", "C", "T", -0.25),
	}}
	r1, err := BuildPRS(context.Background(), model, gt, PRSOptions{
		ScoreScale: "std", Model: ModelIgnore, Backend: numwrap.BackendStd})
	if err != nil {
		t.Fatal(err)
	}
	r2, err := BuildPRS(context.Background(), model, gt, PRSOptions{
		ScoreScale: "std", Model: ModelIgnore, Backend: numwrap.BackendGonum})
	if err != nil {
		t.Fatal(err)
	}
	if err := crossCheckPRS(r1, r2); err != nil {
		t.Fatal(err)
	}
	if _, err := BuildPRS(context.Background(), model, gt, PRSOptions{
		ScoreScale: "std", Model: ModelIgnore, Backend: numwrap.BackendStd, CrossCheck: true}); err != nil {
		t.Fatal(err)
	}
	t.Log("[judgement] prs std vs gonum scores/audit agree")
}
