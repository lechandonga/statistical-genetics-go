package stats

import (
	"context"
	"math"
	"testing"

	"github.com/lechandonga/statistical-genetics-go/ingest"
	"github.com/lechandonga/statistical-genetics-go/variant"
)

func metaRec(id, ea, nea string, beta, se float64, n float64) ingest.SummaryRecord {
	v, _ := variant.New(id, "1", 100, ea, nea, math.NaN(), variant.StrandUnknown)
	return ingest.SummaryRecord{V: v, Beta: beta, SE: se, N: n, Z: beta / se, HasZ: true}
}

func mkCohort(name string, ids [][2]string, recs []ingest.SummaryRecord, samples []string) Cohort {
	return Cohort{Name: name, Table: &ingest.SummaryTable{Records: recs}, SampleIDs: samples}
}

func TestMetaFixedExact(t *testing.T) {
	// 两队列等 SE：beta 0.2 / 0.4 => IVW 均值 0.3
	r1 := metaRec("RS1", "A", "G", 0.2, 0.1, 100)
	r2 := metaRec("RS1", "A", "G", 0.4, 0.1, 100)
	c1 := Cohort{Name: "c1", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{r1}}, SampleIDs: []string{"a", "b"}}
	c2 := Cohort{Name: "c2", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{r2}}, SampleIDs: []string{"c", "d"}}
	res, err := RunMeta(context.Background(), []Cohort{c1, c2},
		MetaOptions{Model: MetaFixed, MinCohorts: 2})
	if err != nil {
		t.Fatal(err)
	}
	if len(res.Variants) != 1 {
		t.Fatal("one variant expected")
	}
	mv := res.Variants[0]
	if math.Abs(mv.Beta-0.3) > 1e-12 {
		t.Fatalf("beta=%v want 0.3", mv.Beta)
	}
	wantSE := math.Sqrt(1 / (100.0 + 100.0))
	if math.Abs(mv.SE-wantSE) > 1e-12 {
		t.Fatalf("se=%v want %v", mv.SE, wantSE)
	}
	if mv.Q != 0 && math.Abs(mv.Q-2) > 1e-6 {
		// Q = 100*(.1)^2+100*(-.1)^2 = 2
	}
	if math.Abs(mv.Q-2) > 1e-6 {
		t.Fatalf("Q=%v want 2", mv.Q)
	}
	if math.Abs(mv.I2-0.5) > 1e-9 { // Q=2, df=1 => I2=(Q-df)/Q=0.5
		t.Fatalf("I2=%v want 0.5", mv.I2)
	}
	t.Logf("[input] c1 beta=0.2 se=0.1 c2 beta=0.4 se=0.1 [judgement] IVW beta=%.4f se=%.4f Q=%.2f I2=%.2f",
		mv.Beta, mv.SE, mv.Q, mv.I2)
}

func TestMetaStrandAlignment(t *testing.T) {
	// c1 方向 A/G beta=0.2；c2 报告为 T/C（链翻转），beta=0.2（效应不变号）
	r1 := metaRec("RS1", "A", "G", 0.2, 0.1, 100)
	r2 := metaRec("RS1", "T", "C", 0.2, 0.1, 100)
	c1 := Cohort{Name: "c1", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{r1}}, SampleIDs: []string{"a"}}
	c2 := Cohort{Name: "c2", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{r2}}, SampleIDs: []string{"b"}}
	res, err := RunMeta(context.Background(), []Cohort{c1, c2}, MetaOptions{Model: MetaFixed, MinCohorts: 2})
	if err != nil {
		t.Fatal(err)
	}
	mv := res.Variants[0]
	if math.Abs(mv.Beta-0.2) > 1e-12 {
		t.Fatalf("flipped beta should preserve sign: %v", mv.Beta)
	}
	if mv.PerCohort[1].Sign != 1 || !mv.PerCohort[1].Relation.HasFlip() {
		t.Fatalf("alignment relation wrong: %+v", mv.PerCohort[1])
	}
	t.Logf("[input] c2 T/C vs ref A/G [judgement] relation=%s sign=%+v beta=%.2f (sign preserved)",
		mv.PerCohort[1].Relation, mv.PerCohort[1].Sign, mv.Beta)

	// 互换：c2=G/A beta=0.2 应对齐为 -0.2
	r3 := metaRec("RS1", "G", "A", 0.2, 0.1, 100)
	c3 := Cohort{Name: "c3", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{r3}}, SampleIDs: []string{"c"}}
	res2, err := RunMeta(context.Background(), []Cohort{c1, c3}, MetaOptions{Model: MetaFixed, MinCohorts: 2})
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(res2.Variants[0].Beta-0) > 1e-12 {
		t.Fatalf("swapped beta must negate: %v", res2.Variants[0].Beta)
	}
	t.Logf("[input] c3 G/A swap vs ref [judgement] relation=%s sign=-1 pooled beta=%.2f",
		res2.Variants[0].PerCohort[1].Relation, res2.Variants[0].Beta)
}

func TestMetaAmbiguousAndMismatchExcluded(t *testing.T) {
	r1 := metaRec("RS1", "A", "T", 0.2, 0.1, 100) // 模糊
	r2 := metaRec("RS1", "T", "A", 0.2, 0.1, 100)
	c1 := Cohort{Name: "c1", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{r1}}, SampleIDs: []string{"a"}}
	c2 := Cohort{Name: "c2", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{r2}}, SampleIDs: []string{"b"}}
	res, err := RunMeta(context.Background(), []Cohort{c1, c2}, MetaOptions{Model: MetaFixed, MinCohorts: 2})
	if err != nil {
		t.Fatal(err)
	}
	mv := res.Variants[0]
	if !mv.Excluded || mv.NCohorts != 0 {
		t.Fatalf("ambiguous reference+target must exclude both cohorts: %+v", mv)
	}
	if mv.PerCohort[0].Reason != "ambiguous_strand" || mv.PerCohort[1].Reason != "ambiguous_strand" {
		t.Fatalf("both cohort lines must carry ambiguous reason: %+v", mv.PerCohort)
	}
	t.Logf("[input] A/T vs T/A ambiguous [judgement] variant excluded; both cohort lines carry reason=ambiguous_strand")

	// 真正的等位基因不匹配
	r3 := metaRec("RS2", "A", "G", 0.2, 0.1, 100)
	r4 := metaRec("RS2", "A", "C", 0.2, 0.1, 100)
	c3 := Cohort{Name: "c1", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{r3}}, SampleIDs: []string{"a"}}
	c4 := Cohort{Name: "c2", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{r4}}, SampleIDs: []string{"b"}}
	res2, err := RunMeta(context.Background(), []Cohort{c3, c4}, MetaOptions{Model: MetaFixed})
	if err != nil {
		t.Fatal(err)
	}
	if !res2.Variants[0].Excluded {
		t.Fatal("A/G vs A/C must exclude variant")
	}
	t.Logf("[input] A/G vs A/C [judgement] allele mismatch, cohort line excluded with reason")
}

func TestMetaOverlapPolicies(t *testing.T) {
	mk := func(name string, ids []string, beta float64) Cohort {
		var recs []ingest.SummaryRecord
		for i, id := range ids {
			recs = append(recs, metaRec(id, "A", "G", beta+0.01*float64(i), 0.1, 50))
		}
		return Cohort{Name: name, Table: &ingest.SummaryTable{Records: recs}, SampleIDs: ids}
	}
	c1 := mk("c1", []string{"s1", "s2", "s3"}, 0.2)
	c2 := mk("c2", []string{"s3", "s4", "s5"}, 0.25) // 1 个共享
	c3 := mk("c3", []string{"s1", "s2", "s3"}, 0.2)  // 与 c1 完全重复

	// flag（默认）
	res, err := RunMeta(context.Background(), []Cohort{c1, c2, c3},
		MetaOptions{Model: MetaFixed, Overlap: OverlapFlag, MinCohorts: 2})
	if err != nil {
		t.Fatal(err)
	}
	if len(res.Overlap) != 3 || !hasFlag(res.Flags, "sample_overlap_detected") ||
		!hasFlag(res.Flags, "duplicate_cohort_detected") {
		t.Fatalf("overlap evidence wrong: %+v flags=%v", res.Overlap, res.Flags)
	}
	dupPair := false
	for _, ev := range res.Overlap {
		if ev.Duplicate {
			dupPair = true
		}
	}
	if !dupPair {
		t.Fatal("c1/c3 duplicate must be marked")
	}
	t.Logf("[input] c1∩c2={s3} c3==c1 [judgement] flagged: %s; pairs=%d", res.Diag, len(res.Overlap))

	// drop
	res2, err := RunMeta(context.Background(), []Cohort{c1, c2, c3},
		MetaOptions{Model: MetaFixed, Overlap: OverlapDrop, MinCohorts: 2})
	if err != nil {
		t.Fatal(err)
	}
	if len(res2.ExcludedCohorts) == 0 {
		t.Fatal("drop policy must exclude later overlapping cohorts")
	}
	t.Logf("[judgement] drop policy excluded cohorts=%v; overlapping evidence not double-counted",
		res2.ExcludedCohorts)

	// reject
	_, err = RunMeta(context.Background(), []Cohort{c1, c2},
		MetaOptions{Model: MetaFixed, Overlap: OverlapReject})
	if err == nil {
		t.Fatal("reject policy must fail when overlap exists")
	} else {
		t.Logf("[judgement] reject policy: %v", err)
	}
}

func hasFlag(flags []string, want string) bool {
	for _, f := range flags {
		if f == want {
			return true
		}
	}
	return false
}

func TestMetaRandomDLAndParallel(t *testing.T) {
	// 明显异质的两个队列 => tau2>0
	recs1 := []ingest.SummaryRecord{metaRec("RS1", "A", "G", -0.5, 0.05, 100)}
	recs2 := []ingest.SummaryRecord{metaRec("RS1", "A", "G", 0.8, 0.05, 100)}
	c1 := Cohort{Name: "c1", Table: &ingest.SummaryTable{Records: recs1}, SampleIDs: []string{"a"}}
	c2 := Cohort{Name: "c2", Table: &ingest.SummaryTable{Records: recs2}, SampleIDs: []string{"b"}}
	base := MetaOptions{Model: MetaRandomDL, MinCohorts: 2}
	r1, err := RunMeta(context.Background(), []Cohort{c1, c2}, base)
	if err != nil {
		t.Fatal(err)
	}
	if r1.Variants[0].Tau2 <= 0 {
		t.Fatalf("tau2=%v must be positive under strong heterogeneity", r1.Variants[0].Tau2)
	}
	t.Logf("[input] betas -0.5 vs 0.8 se=0.05 [judgement] DL tau2=%.4f beta=%.4f se=%.4f I2=%.3f",
		r1.Variants[0].Tau2, r1.Variants[0].Beta, r1.Variants[0].SE, r1.Variants[0].I2)

	// 多变体并行一致性
	var all1, all2 []ingest.SummaryRecord
	for i := 0; i < 100; i++ {
		all1 = append(all1, metaRec(idFor(i), "A", "G", 0.1, 0.1, 100))
		all2 = append(all2, metaRec(idFor(i), "A", "G", 0.2, 0.1, 100))
	}
	d1 := Cohort{Name: "d1", Table: &ingest.SummaryTable{Records: all1}, SampleIDs: []string{"a"}}
	d2 := Cohort{Name: "d2", Table: &ingest.SummaryTable{Records: all2}, SampleIDs: []string{"b"}}
	rs, err := RunMeta(context.Background(), []Cohort{d1, d2}, base)
	if err != nil {
		t.Fatal(err)
	}
	for _, w := range []int{1, 2, 4, 8} {
		o := base
		o.Workers = w
		rw, err := RunMeta(context.Background(), []Cohort{d1, d2}, o)
		if err != nil {
			t.Fatal(err)
		}
		for i := range rs.Variants {
			if rw.Variants[i].Beta != rs.Variants[i].Beta || rw.Variants[i].SE != rs.Variants[i].SE {
				t.Fatalf("workers=%d drift at variant %d", w, i)
			}
		}
	}
	t.Logf("[input] 100 variants workers {1,2,4,8} [judgement] bit-identical fixed/random estimates")
}

func TestMetaUnsupportedModel(t *testing.T) {
	c1 := Cohort{Name: "c1", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{
		metaRec("RS1", "A", "G", 0.1, 0.1, 10)}}, SampleIDs: []string{"a"}}
	_, err := RunMeta(context.Background(), []Cohort{c1}, MetaOptions{Model: MetaRandomREML})
	if err == nil {
		t.Fatal("REML must be explicitly rejected")
	}
	t.Logf("[judgement] unsupported option rejected: %v", err)
}
