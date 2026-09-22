package stats

import (
	"context"
	"math"
	"strings"
	"testing"

	"github.com/lechandonga/statistical-genetics-go/ingest"
	"github.com/lechandonga/statistical-genetics-go/variant"
)

// buildGT 直接由 (变体键 -> 样本剂量序列) 构造矩阵，避免重复走解析器。
func buildGT(samples []string, vars []*variant.Variant, doses [][]float64, miss [][]bool) *ingest.GenotypeMatrix {
	nS := len(samples)
	m := &ingest.GenotypeMatrix{Variants: vars, SampleIDs: samples}
	for i := 0; i < len(vars); i++ {
		for j := 0; j < nS; j++ {
			m.Dosage = append(m.Dosage, doses[i][j])
			m.Missing = append(m.Missing, miss[i][j])
		}
	}
	return m
}

func prsModelRec(id, ea, nea string, beta float64) ingest.SummaryRecord {
	v, _ := variant.New(id, "1", 100, ea, nea, math.NaN(), variant.StrandUnknown)
	return ingest.SummaryRecord{V: v, Beta: beta, SE: 0.1, P: 0.001, N: 100, Z: beta / 0.1, HasZ: true}
}

func TestPRSAlleleAlignment(t *testing.T) {
	samples := []string{"s1", "s2"}
	v1, _ := variant.New("RS1", "1", 1, "A", "G", math.NaN(), variant.StrandUnknown)
	v2, _ := variant.New("RS2", "1", 2, "T", "G", math.NaN(), variant.StrandUnknown)
	v3, _ := variant.New("RS3", "1", 3, "A", "G", math.NaN(), variant.StrandUnknown)
	gt := buildGT(samples,
		[]*variant.Variant{v1, v2, v3},
		[][]float64{{2, 0}, {2, 0}, {2, 0}},
		[][]bool{{false, false}, {false, false}, {false, false}})

	model := &ingest.SummaryTable{Records: []ingest.SummaryRecord{
		prsModelRec("RS1", "A", "G", 1), // 直接匹配：+2
		prsModelRec("RS2", "A", "C", 1), // 矩阵 T/G：翻转后直接匹配，符号不变：+2
		prsModelRec("RS3", "G", "A", 1), // 互换 -> sign -1：-2
	}}
	res, err := BuildPRS(context.Background(), model, gt, PRSOptions{ScoreScale: "raw"})
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("[input] model sites aligned: direct, flip+swap, swap [judgement] %s", res.Diag)
	// s1 剂量全为 2：+2(直接) +2(翻转,符号不变) -2(互换) = 2
	if res.Scores[0] != 2 {
		t.Fatalf("s1 score=%v want 2", res.Scores[0])
	}
	// s2 剂量全为 0：得分 0
	if res.Scores[1] != 0 {
		t.Fatalf("s2 score=%v want 0", res.Scores[1])
	}
	// 纯链翻转（不互换）应保留符号
	v4, _ := variant.New("RS4", "1", 4, "A", "C", math.NaN(), variant.StrandUnknown)
	gt2 := buildGT(samples, []*variant.Variant{v4},
		[][]float64{{2, 0}}, [][]bool{{false, false}})
	m2 := &ingest.SummaryTable{Records: []ingest.SummaryRecord{prsModelRec("RS4", "T", "G", 3)}}
	r2, err := BuildPRS(context.Background(), m2, gt2, PRSOptions{ScoreScale: "raw"})
	if err != nil {
		t.Fatal(err)
	}
	if r2.Scores[0] != 6 { // 2*3，符号保留
		t.Fatalf("flip-only must preserve sign: %v", r2.Scores[0])
	}
	t.Logf("[input] flip-only beta=3 dosage=2 [judgement] score=6, sign preserved")
}

func TestPRSAmbiguousRejected(t *testing.T) {
	samples := []string{"s1"}
	va, _ := variant.New("RS1", "1", 1, "A", "T", math.NaN(), variant.StrandUnknown)
	gt := buildGT(samples, []*variant.Variant{va}, [][]float64{{2}}, [][]bool{{false}})
	model := &ingest.SummaryTable{Records: []ingest.SummaryRecord{prsModelRec("RS1", "T", "A", 1)}}
	_, err := BuildPRS(context.Background(), model, gt, PRSOptions{ScoreScale: "raw"})
	if err == nil || !strings.Contains(err.Error(), "no_usable_sites") {
		t.Fatalf("ambiguous site unusable under reject policy: %v", err)
	}
	t.Logf("[input] A/T vs T/A in PRS [judgement] %v", err)
}

func TestPRSMissingGenotypePolicies(t *testing.T) {
	samples := []string{"s1", "s2"}
	v1, _ := variant.New("RS1", "1", 1, "A", "G", 0.3, variant.StrandUnknown)
	// EAF=0.3 -> 期望剂量 0.6
	gt := buildGT(samples, []*variant.Variant{v1},
		[][]float64{{math.NaN(), 2}}, [][]bool{{true, false}})
	model := &ingest.SummaryTable{Records: []ingest.SummaryRecord{prsModelRec("RS1", "A", "G", 1)}}

	r1, err := BuildPRS(context.Background(), model, gt, PRSOptions{ScoreScale: "raw", GT: GTMeanImpute})
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(r1.Scores[0]-0.6) > 1e-12 {
		t.Fatalf("mean impute score=%v want 0.6", r1.Scores[0])
	}
	if r1.NActive[0] != 1 {
		t.Fatalf("imputed site still counts active: %v", r1.NActive)
	}
	if r1.Audit[0].MissingGT != 1 {
		t.Fatalf("missing count=%d", r1.Audit[0].MissingGT)
	}
	t.Logf("[input] s1 missing GT, EAF=0.3 [judgement] impute 2*EAF=0.6; s2 observed score=2")

	r2, err := BuildPRS(context.Background(), model, gt, PRSOptions{ScoreScale: "raw", GT: GTExclude})
	if err != nil {
		t.Fatal(err)
	}
	if r2.Scores[0] != 0 || r2.NActive[0] != 0 {
		t.Fatalf("exclude policy: score=%v active=%d want 0/0", r2.Scores[0], r2.NActive[0])
	}
	t.Logf("[input] s1 missing GT, policy=exclude [judgement] score=0 active=0; s2 unaffected")
}

func TestPRSMissingModelSites(t *testing.T) {
	samples := []string{"s1"}
	v1, _ := variant.New("RS1", "1", 1, "A", "G", math.NaN(), variant.StrandUnknown)
	gt := buildGT(samples, []*variant.Variant{v1}, [][]float64{{2}}, [][]bool{{false}})
	model := &ingest.SummaryTable{Records: []ingest.SummaryRecord{
		prsModelRec("RS1", "A", "G", 1),
		prsModelRec("RS2", "A", "G", 2), // 矩阵中不存在
	}}
	// strict 默认拒绝
	_, err := BuildPRS(context.Background(), model, gt, PRSOptions{ScoreScale: "raw"})
	if err == nil || !strings.Contains(err.Error(), "effect_sites_missing") {
		t.Fatalf("strict policy must reject missing effect sites: %v", err)
	}
	t.Logf("[input] RS2 requested but absent in GT [judgement] %v", err)

	// ignore：仅用交集，显式列出
	r, err := BuildPRS(context.Background(), model, gt, PRSOptions{ScoreScale: "raw", Model: ModelIgnore})
	if err != nil {
		t.Fatal(err)
	}
	if r.Scores[0] != 2 || len(r.MissingModelSites) == 0 {
		t.Fatalf("intersection scoring wrong: %v missing=%v", r.Scores, r.MissingModelSites)
	}
	t.Logf("[input] same, policy=ignore [judgement] score on intersection=%v; missing=%v (caliber change explicit)",
		r.Scores, r.MissingModelSites)
}

func TestPRSScalesAndParallel(t *testing.T) {
	samples := []string{"s1", "s2", "s3", "s4"}
	v1, _ := variant.New("RS1", "1", 1, "A", "G", math.NaN(), variant.StrandUnknown)
	gt := buildGT(samples, []*variant.Variant{v1},
		[][]float64{{0, 1, 1, 2}}, [][]bool{{false, false, false, false}})
	model := &ingest.SummaryTable{Records: []ingest.SummaryRecord{prsModelRec("RS1", "A", "G", 4)}}
	// raw scores 0,4,4,8：mean 4, sample sd = sqrt(32/3)
	base := PRSOptions{ScoreScale: "std"}
	r1, err := BuildPRS(context.Background(), model, gt, base)
	if err != nil {
		t.Fatal(err)
	}
	sm := 0.0
	for _, x := range r1.Scores {
		sm += x
	}
	mean := sm / float64(len(r1.Scores))
	var ss float64
	for _, x := range r1.Scores {
		ss += (x - mean) * (x - mean)
	}
	sd := math.Sqrt(ss / float64(len(r1.Scores)-1))
	if math.Abs(mean) > 1e-12 || math.Abs(sd-1) > 1e-12 {
		t.Fatalf("std scale mean=%v sd=%v", mean, sd)
	}
	t.Logf("[input] raw 0/4/4/8 [judgement] standardized mean=0 sd=1")

	// 并行逐位一致
	for _, w := range []int{1, 2, 3, 4} {
		o := base
		o.Workers = w
		rw, err := BuildPRS(context.Background(), model, gt, o)
		if err != nil {
			t.Fatal(err)
		}
		for j := range r1.Scores {
			if rw.Scores[j] != r1.Scores[j] {
				t.Fatalf("workers=%d drift: %v vs %v", w, rw.Scores, r1.Scores)
			}
		}
	}
	t.Logf("[input] workers {1,2,3,4} [judgement] bit-identical standardized scores")

	// perallele
	rp, err := BuildPRS(context.Background(), model, gt, PRSOptions{ScoreScale: "perallele"})
	if err != nil {
		t.Fatal(err)
	}
	wantPA := []float64{0, 4, 4, 8} // 单位点：beta*dosage / 1
	for j, sc := range rp.Scores {
		if sc != wantPA[j] {
			t.Fatalf("perallele score[%d]=%v want %v", j, sc, wantPA[j])
		}
	}
	t.Logf("[input] perallele scale [judgement] scores 0/4/4/8 normalized per active site")
}

func TestPRSPThreshold(t *testing.T) {
	samples := []string{"s1"}
	v1, _ := variant.New("RS1", "1", 1, "A", "G", math.NaN(), variant.StrandUnknown)
	v2, _ := variant.New("RS2", "1", 2, "A", "G", math.NaN(), variant.StrandUnknown)
	gt := buildGT(samples, []*variant.Variant{v1, v2},
		[][]float64{{2}, {2}}, [][]bool{{false}, {false}})
	r1 := prsModelRec("RS1", "A", "G", 1)
	r1.P = 0.001
	r2 := prsModelRec("RS2", "A", "G", 1)
	r2.P = 0.5
	model := &ingest.SummaryTable{Records: []ingest.SummaryRecord{r1, r2}}
	// RS2 被阈值过滤，模型位点都在矩阵中，无缺失
	r, err := BuildPRS(context.Background(), model, gt, PRSOptions{
		ScoreScale: "raw", PThreshold: 0.05, Model: ModelStrict})
	if err != nil {
		t.Fatal(err)
	}
	if r.Scores[0] != 2 { // 仅 RS1: 1*2
		t.Fatalf("score=%v want 2", r.Scores[0])
	}
	t.Logf("[input] p threshold 0.05: RS1 p=.001 RS2 p=.5 [judgement] only RS1 used, score=2")

	// 不支持的尺度必须显式拒绝
	if _, err := BuildPRS(context.Background(), model, gt, PRSOptions{ScoreScale: "weird"}); err == nil {
		t.Fatal("unsupported scale must be rejected")
	}
}

func TestPRSEvaluate(t *testing.T) {
	prs := &PRSResult{Scores: []float64{-1.5, -0.5, 0.5, 1.5}}
	// y = 2*x + 1
	y := []float64{-2, 0, 2, 4}
	res, err := EvaluatePRS(context.Background(), prs, EvalOptions{Phenotype: y})
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(res.Beta-2) > 1e-9 || math.Abs(res.R2-1) > 1e-12 {
		t.Fatalf("beta=%v R2=%v", res.Beta, res.R2)
	}
	t.Logf("[input] y=2x+1 [judgement] eval beta=%.2f se=%.4f p=%.2g R2=%.3f", res.Beta, res.SE, res.P, res.R2)

	// 协变量调整当前不支持，须显式拒绝
	_, err = EvaluatePRS(context.Background(), prs, EvalOptions{
		Phenotype: y, Covariates: [][]float64{{1, 1, 1, 1}}})
	if err == nil {
		t.Fatal("covariates must be explicitly rejected")
	}
	t.Logf("[judgement] unsupported covariate option rejected: %v", err)

	// 长度不匹配
	_, err = EvaluatePRS(context.Background(), prs, EvalOptions{Phenotype: []float64{1, 2}})
	if err == nil {
		t.Fatal("dimension mismatch must error")
	}
}
