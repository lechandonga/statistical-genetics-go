package stats

import (
	"context"
	"math"
	"testing"

	"github.com/lechandonga/statistical-genetics-go/ingest"
	"github.com/lechandonga/statistical-genetics-go/variant"
)

// synthH2Tables 构造满足 chi2 ≈ 1 + slope*l2 关系的确定性数据，
// h2 真值 = slope*M/N。
func synthH2Tables(t *testing.T, nSNP int, nBar float64, h2 float64,
	flipEvery int, addAmbiguous bool) (*ingest.SummaryTable, *ingest.LDScoreTable) {
	t.Helper()
	M := float64(nSNP)
	slope := h2 * nBar / M
	sumRecs := make([]ingest.SummaryRecord, 0, nSNP+2)
	ldRecs := make([]ingest.LDScoreRecord, 0, nSNP)
	for i := 0; i < nSNP; i++ {
		// 非模糊等位基因循环 A/G, C/T, A/C, G/T
		pairs := [][2]string{{"A", "G"}, {"C", "T"}, {"A", "C"}, {"G", "T"}}
		p := pairs[i%4]
		id := idFor(i)
		ea, nea := p[0], p[1]
		l2 := 1.0 + float64(i%7) // 1..7
		// 无噪声确定性关系：chi2 = 1 + slope*l2；z=sqrt(chi2)
		chi2 := 1 + slope*l2
		z := math.Sqrt(chi2)
		se := 0.05
		beta := z * se

		// LD 参考保持 EA=p[0]
		lv, _ := variant.New(id, "1", int64(100+i), ea, nea, math.NaN(), variant.StrandUnknown)
		ldRecs = append(ldRecs, ingest.LDScoreRecord{V: lv, L2: l2, N: nBar})

		// 汇总侧按 flipEvery 注入链翻转（效应符号与 chi2 无关，但测对齐路径）
		sea, snea := ea, nea
		if flipEvery > 0 && i%flipEvery == 0 {
			c1, _ := variant.ComplementBase(ea[0])
			c2, _ := variant.ComplementBase(nea[0])
			sea, snea = string(c1), string(c2)
		}
		sv, _ := variant.New(id, "1", int64(100+i), sea, snea, math.NaN(), variant.StrandUnknown)
		sumRecs = append(sumRecs, ingest.SummaryRecord{
			V: sv, Beta: beta, SE: se, P: math.NaN(), N: nBar, Z: z, HasZ: true})
	}
	if addAmbiguous {
		id := "RSAMB"
		lv, _ := variant.New(id, "1", 999, "A", "T", math.NaN(), variant.StrandUnknown)
		ldRecs = append(ldRecs, ingest.LDScoreRecord{V: lv, L2: 3, N: nBar})
		sv, _ := variant.New(id, "1", 999, "T", "A", math.NaN(), variant.StrandUnknown)
		sumRecs = append(sumRecs, ingest.SummaryRecord{
			V: sv, Beta: 0.1, SE: 0.05, N: nBar, Z: 2, HasZ: true})
	}
	return &ingest.SummaryTable{Records: sumRecs},
		&ingest.LDScoreTable{Records: ldRecs}
}

func idFor(i int) string {
	// rs 编号大写归一化后保持唯一
	return "RS" + itoa(i)
}

func itoa(i int) string {
	if i == 0 {
		return "0"
	}
	var b [12]byte
	n := len(b)
	for i > 0 {
		n--
		b[n] = byte('0' + i%10)
		i /= 10
	}
	return string(b[n:])
}

func TestEstimateH2RecoversTruth(t *testing.T) {
	const nSNP = 200
	const nBar = 5000
	const truth = 0.4
	sum, ld := synthH2Tables(t, nSNP, nBar, truth, 0, false)
	opts := H2Options{FreeIntercept: true, MinVariants: 50, Align: variant.DefaultAlignConfig(), Workers: 1}
	res, err := EstimateH2(context.Background(), sum, ld, opts)
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("[input] synthetic M=%d N=%d truth h2=%.2f chi2=1+slope*l2 (no noise)", nSNP, nBar, truth)
	t.Logf("[judgement] %s", res.Diag)
	if math.Abs(res.H2-truth) > 0.02 {
		t.Fatalf("h2=%.4f want ~%.2f", res.H2, truth)
	}
	if math.Abs(res.Intercept-1) > 0.02 {
		t.Fatalf("intercept=%.4f want ~1", res.Intercept)
	}
}

func TestEstimateH2ParallelEqualsSerial(t *testing.T) {
	sum, ld := synthH2Tables(t, 300, 4000, 0.3, 5, false)
	base := H2Options{FreeIntercept: true, MinVariants: 50, Align: variant.DefaultAlignConfig()}

	r1, err := EstimateH2(context.Background(), sum, ld, base)
	if err != nil {
		t.Fatal(err)
	}
	for _, w := range []int{2, 3, 4, 8, 16, 100} {
		o := base
		o.Workers = w
		rw, err := EstimateH2(context.Background(), sum, ld, o)
		if err != nil {
			t.Fatal(err)
		}
		if rw.H2 != r1.H2 || rw.SE != r1.SE || rw.Intercept != r1.Intercept {
			t.Fatalf("workers=%d drift: h2 serial=%v parallel=%v", w, r1.H2, rw.H2)
		}
	}
	t.Logf("[input] M=300 with every-5th SNP strand-flipped; workers in {1..100}")
	t.Logf("[judgement] bit-identical h2/se/intercept across parallelism; flips aligned, h2=%.4f", r1.H2)
}

func TestEstimateH2NegativeAndBoundaryPolicies(t *testing.T) {
	// 构造负斜率：chi2 = 1 - slope*l2（小斜率，保持 chi2>0）
	nSNP, nBar := 200, 5000.0
	var sumRecs []ingest.SummaryRecord
	var ldRecs []ingest.LDScoreRecord
	for i := 0; i < nSNP; i++ {
		pairs := [][2]string{{"A", "G"}, {"C", "T"}}
		p := pairs[i%2]
		id := idFor(i)
		l2 := 1 + float64(i%5)
		chi2 := 1 - 0.03*l2
		z := math.Sqrt(chi2)
		lv, _ := variant.New(id, "1", int64(100+i), p[0], p[1], math.NaN(), variant.StrandUnknown)
		sv, _ := variant.New(id, "1", int64(100+i), p[0], p[1], math.NaN(), variant.StrandUnknown)
		ldRecs = append(ldRecs, ingest.LDScoreRecord{V: lv, L2: l2})
		sumRecs = append(sumRecs, ingest.SummaryRecord{V: sv, SE: 1, Z: z, HasZ: true, N: nBar})
	}
	sum := &ingest.SummaryTable{Records: sumRecs}
	ld := &ingest.LDScoreTable{Records: ldRecs}

	// report raw：负值保留并置标志
	o := H2Options{FreeIntercept: true, MinVariants: 50, Boundary: H2ReportRaw, NBar: nBar}
	r, err := EstimateH2(context.Background(), sum, ld, o)
	if err != nil {
		t.Fatal(err)
	}
	if r.H2Raw >= 0 || r.Flags&H2FlagNegativeRaw == 0 {
		t.Fatalf("negative raw not flagged: h2=%v flags=%s", r.H2Raw, r.Flags)
	}
	t.Logf("[input] chi2=1-0.03*l2 [judgement] raw h2=%.4f flag=%s, kept interpretable via flag", r.H2Raw, r.Flags)

	// clamp：截断为 0
	o.Boundary = H2Clamp
	r2, err := EstimateH2(context.Background(), sum, ld, o)
	if err != nil {
		t.Fatal(err)
	}
	if r2.H2 != 0 || r2.Flags&H2FlagClamped == 0 {
		t.Fatalf("clamp failed: h2=%v flags=%s", r2.H2, r2.Flags)
	}
	if r2.CI95[0] >= 0 {
		t.Fatal("CI must remain based on raw distribution (negative lower bound)")
	}
	t.Logf("[judgement] clamp policy: reported h2=0, CI stays [%.4f,%.4f], flag=%s",
		r2.CI95[0], r2.CI95[1], r2.Flags)

	// reject：显式错误
	o.Boundary = H2Reject
	if _, err := EstimateH2(context.Background(), sum, ld, o); err == nil {
		t.Fatal("reject policy must error on negative h2")
	} else {
		t.Logf("[judgement] reject policy: %v", err)
	}
}

func TestEstimateH2AmbiguousDropped(t *testing.T) {
	sum, ld := synthH2Tables(t, 100, 3000, 0.2, 0, true)
	o := H2Options{FreeIntercept: true, MinVariants: 50, Align: variant.DefaultAlignConfig()}
	r, err := EstimateH2(context.Background(), sum, ld, o)
	if err != nil {
		t.Fatal(err)
	}
	if r.Flags&H2FlagAmbiguousDropped == 0 || r.NSnps != 100 {
		t.Fatalf("ambiguous not dropped: n=%d flags=%s", r.NSnps, r.Flags)
	}
	if len(r.Rejected) == 0 {
		t.Fatal("rejected list must record ambiguous site")
	}
	t.Logf("[input] 100 SNPs + 1 A/T ambiguous [judgement] n=%d ambiguous flagged, reason recorded", r.NSnps)
}

func TestEstimateH2InsufficientAndMissingN(t *testing.T) {
	sum, ld := synthH2Tables(t, 10, 1000, 0.2, 0, false)
	o := H2Options{FreeIntercept: true, MinVariants: 50}
	if _, err := EstimateH2(context.Background(), sum, ld, o); err == nil {
		t.Fatal("fewer than MinVariants must reject")
	} else {
		t.Logf("[judgement] insufficient variants: %v", err)
	}

	// 清空汇总与 LD 两侧的 N
	for i := range sum.Records {
		sum.Records[i].N = 0
	}
	for i := range ld.Records {
		ld.Records[i].N = 0
	}
	o.MinVariants = 5
	if _, err := EstimateH2(context.Background(), sum, ld, o); err == nil {
		t.Fatal("missing N without NBar must reject")
	} else {
		t.Logf("[judgement] missing sample size: %v", err)
	}
}
