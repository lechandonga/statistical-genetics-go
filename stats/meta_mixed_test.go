package stats

import (
	"context"
	"math"
	"testing"

	"github.com/lechandonga/statistical-genetics-go/ingest"
	"github.com/lechandonga/statistical-genetics-go/numwrap"
	"github.com/lechandonga/statistical-genetics-go/variant"
)

func zRec(id, ea, nea string, z, n float64, eaf float64) ingest.SummaryRecord {
	v, _ := variant.New(id, "1", 100, ea, nea, eaf, variant.StrandUnknown)
	return ingest.SummaryRecord{V: v, Beta: 0, SE: math.NaN(), N: n, Z: z, HasZ: true}
}

// TestMetaMixedZAndBetaCohorts 验证 Z+N(+EAF) 队列按统一口径真正参与
// 合并；beta/SE 队列照常参与。
func TestMetaMixedZAndBetaCohorts(t *testing.T) {
	// beta 队列：beta=0.2, se=0.1
	r1 := metaRec("RS1", "A", "G", 0.2, 0.1, 100)
	// Z-only 队列：Z=3, N=400, EAF=0.5 => SE=1/sqrt(2*400*0.25)=1/200^0.5≈0.0707
	// beta = 3*SE ≈ 0.21213
	r2 := zRec("RS1", "A", "G", 3.0, 400, 0.5)
	c1 := Cohort{Name: "beta", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{r1}}, SampleIDs: []string{"a"}}
	c2 := Cohort{Name: "zonly", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{r2}}, SampleIDs: []string{"b"}}
	res, err := RunMeta(context.Background(), []Cohort{c1, c2},
		MetaOptions{Model: MetaFixed, MinCohorts: 2, Backend: numwrap.BackendStd})
	if err != nil {
		t.Fatal(err)
	}
	mv := res.Variants[0]
	if mv.Excluded {
		t.Fatalf("mixed variant must be included: reason=%s", mv.ExcludeReason)
	}
	if mv.NCohorts != 2 {
		t.Fatalf("both cohorts must contribute, got %d", mv.NCohorts)
	}
	if l := mv.PerCohort[1]; l.Excluded || l.Basis != BasisZNEAF {
		t.Fatalf("z-only cohort line wrong: %+v", l)
	}
	seWant := 1 / math.Sqrt(2*400*0.5*0.5)
	if math.Abs(mv.PerCohort[1].SE-seWant) > 1e-12 {
		t.Fatalf("z-only SE=%v want %v", mv.PerCohort[1].SE, seWant)
	}
	if math.Abs(mv.PerCohort[1].BetaAligned-3*seWant) > 1e-12 {
		t.Fatalf("z-only beta=%v", mv.PerCohort[1].BetaAligned)
	}
	if l := mv.PerCohort[0]; l.Basis != BasisBetaSE {
		t.Fatalf("beta cohort basis=%s want beta_se", l.Basis)
	}
	t.Logf("[input] beta(0.2,0.1) + Z(3,N=400,EAF=0.5) [judgement] pooled beta=%.6f se=%.6f",
		mv.Beta, mv.SE)
	t.Logf("[audit] z-only: %s", mv.PerCohort[1].BasisDetail)

	// 结果必须全部有限（不允许“参与了却不可解释”）
	if math.IsNaN(mv.Beta) || math.IsNaN(mv.SE) || math.IsNaN(mv.P) {
		t.Fatal("pooled result must be finite")
	}
}

// TestMetaZOnlyExclusionReasons 验证三类排除原因可区分：
// 缺 N、EAF 缺失（无方向/尺度）、该变体在队列缺失。
func TestMetaZOnlyExclusionReasons(t *testing.T) {
	r1 := metaRec("RS1", "A", "G", 0.2, 0.1, 100)

	// N 缺失
	noN := zRec("RS1", "A", "G", 2.0, 0, 0.5)
	// EAF 缺失
	noEAF := zRec("RS1", "A", "G", 2.0, 100, math.NaN())
	// EAF=1（剂量方差 0）
	boundaryEAF := zRec("RS1", "A", "G", 2.0, 100, 1.0)

	mk := func(name string, r ingest.SummaryRecord) Cohort {
		return Cohort{Name: name,
			Table:     &ingest.SummaryTable{Records: []ingest.SummaryRecord{r}},
			SampleIDs: []string{name}}
	}
	c0 := Cohort{Name: "ref", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{r1}}, SampleIDs: []string{"x"}}

	for _, tc := range []struct {
		name   string
		rec    ingest.SummaryRecord
		reason variant.Reason
	}{
		{"noN", noN, ReasonZMissingN},
		{"noEAF", noEAF, ReasonZUnorientableMAF},
		{"eaf1", boundaryEAF, ReasonZUnorientableMAF},
	} {
		t.Run(tc.name, func(t *testing.T) {
			c1 := mk(tc.name, tc.rec)
			res, err := RunMeta(context.Background(), []Cohort{c0, c1},
				MetaOptions{Model: MetaFixed, MinCohorts: 2})
			if err != nil {
				t.Fatal(err)
			}
			mv := res.Variants[0]
			if !mv.Excluded || mv.ExcludeReason != "insufficient_cohorts" {
				t.Fatalf("variant should be excluded insufficient_cohorts: %+v reason=%s",
					mv.Excluded, mv.ExcludeReason)
			}
			var found bool
			for _, l := range mv.PerCohort {
				if l.Cohort == tc.name {
					if !l.Excluded || l.Reason != tc.reason {
						t.Fatalf("cohort %s: excluded=%v reason=%s want %s",
							l.Cohort, l.Excluded, l.Reason, tc.reason)
					}
					found = true
				}
			}
			if !found {
				t.Fatal("cohort missing from audit")
			}
			t.Logf("[judgement] %s -> variant excluded, cohort reason=%s (distinguishable)",
				tc.name, tc.reason)
		})
	}

	// 方向依据可审计：Z-only 队列在链互换位点上 Z 应取反
	r2 := zRec("RS1", "G", "A", 3.0, 400, 0.5)
	cSwap := Cohort{Name: "zswap",
		Table:     &ingest.SummaryTable{Records: []ingest.SummaryRecord{r2}},
		SampleIDs: []string{"s"}}
	res, err := RunMeta(context.Background(), []Cohort{c0, cSwap},
		MetaOptions{Model: MetaFixed, MinCohorts: 2})
	if err != nil {
		t.Fatal(err)
	}
	l := res.Variants[0].PerCohort[1]
	if l.Excluded || l.Sign != -1 || l.Z != -3 {
		t.Fatalf("z-only swap sign audit wrong: %+v", l)
	}
	t.Logf("[audit] swapped z-only cohort sign=-1 Z_aligned=%.1f basis=%s", l.Z, l.Basis)
}

// TestMetaCrosscheckJudgementConsistency 双后端必须给出一致数值与
// 完全一致的判定（纳入/排除、原因、权重审计）。
func TestMetaCrosscheckJudgementConsistency(t *testing.T) {
	r1 := metaRec("RS1", "A", "G", 0.2, 0.1, 100)
	r2 := zRec("RS1", "A", "G", 3.0, 400, 0.5)
	c1 := Cohort{Name: "a", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{r1}}, SampleIDs: []string{"a"}}
	c2 := Cohort{Name: "b", Table: &ingest.SummaryTable{Records: []ingest.SummaryRecord{r2}}, SampleIDs: []string{"b"}}

	rs, err := RunMeta(context.Background(), []Cohort{c1, c2},
		MetaOptions{Model: MetaRandomDL, MinCohorts: 2, Backend: numwrap.BackendStd})
	if err != nil {
		t.Fatal(err)
	}
	rg, err := RunMeta(context.Background(), []Cohort{c1, c2},
		MetaOptions{Model: MetaRandomDL, MinCohorts: 2, Backend: numwrap.BackendGonum})
	if err != nil {
		t.Fatal(err)
	}
	if err := compareMetaResults(rs, rg); err != nil {
		t.Fatal(err)
	}
	rc, err := RunMeta(context.Background(), []Cohort{c1, c2},
		MetaOptions{Model: MetaRandomDL, MinCohorts: 2, Backend: numwrap.BackendCrosscheck})
	if err != nil {
		t.Fatal(err)
	}
	if rc.Variants[0].Beta != rs.Variants[0].Beta {
		t.Fatal("crosscheck must return stdlib reference value")
	}
	t.Logf("[judgement] stdlib/gonum/crosscheck agree: beta=%.10f se=%.10f p=%.6g",
		rs.Variants[0].Beta, rs.Variants[0].SE, rs.Variants[0].P)
}

// TestMetaConcurrentBackendSwitchRaceFree 并发切换默认后端时，按调用
// 指定后端的任务互不干扰（-race 必须干净）。
func TestMetaConcurrentBackendSwitchRaceFree(t *testing.T) {
	var recs1, recs2 []ingest.SummaryRecord
	for i := 0; i < 80; i++ {
		recs1 = append(recs1, metaRec(idFor(i), "A", "G", 0.1, 0.1, 100))
		recs2 = append(recs2, zRec(idFor(i), "A", "G", 2.0, 200, 0.3))
	}
	c1 := Cohort{Name: "a", Table: &ingest.SummaryTable{Records: recs1}, SampleIDs: []string{"a"}}
	c2 := Cohort{Name: "b", Table: &ingest.SummaryTable{Records: recs2}, SampleIDs: []string{"b"}}

	done := make(chan struct{})
	go func() {
		defer close(done)
		for i := 0; i < 100; i++ {
			_ = numwrap.SetDefaultBackend(numwrap.BackendGonum)
			_ = numwrap.SetDefaultBackend(numwrap.BackendStd)
		}
	}()
	ref, err := RunMeta(context.Background(), []Cohort{c1, c2},
		MetaOptions{Model: MetaFixed, Backend: numwrap.BackendStd, Workers: 4})
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i < 100; i++ {
		b := numwrap.BackendStd
		if i%2 == 0 {
			b = numwrap.BackendGonum
		}
		got, err := RunMeta(context.Background(), []Cohort{c1, c2},
			MetaOptions{Model: MetaFixed, Backend: b, Workers: 2})
		if err != nil {
			t.Fatal(err)
		}
		for k := range got.Variants {
			if !numwrap.WithinTol(got.Variants[k].Beta, ref.Variants[k].Beta) {
				t.Fatalf("backend task contaminated by default switch at %d", k)
			}
		}
	}
	<-done
}
