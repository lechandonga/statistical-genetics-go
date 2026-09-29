package stats

import (
	"context"
	"math"
	"testing"

	"github.com/lechandonga/statistical-genetics-go/ingest"
	"github.com/lechandonga/statistical-genetics-go/numwrap"
	"github.com/lechandonga/statistical-genetics-go/variant"
)

// noisyH2Tables 构造带确定性噪声的 chi2 = 1 + slope*l2 + eps 数据，
// 使两后端的残差方差都建立在真实 RSS 之上（避免完美拟合下的浮点伪差）。
func noisyH2Tables(t *testing.T, nSNP int, nBar, h2 float64) (*ingest.SummaryTable, *ingest.LDScoreTable) {
	t.Helper()
	M := float64(nSNP)
	slope := h2 * nBar / M
	sumRecs := make([]ingest.SummaryRecord, 0, nSNP)
	ldRecs := make([]ingest.LDScoreRecord, 0, nSNP)
	// 固定种子 xorshift，保证跨运行、跨后端一致。
	seed := uint64(0x9e3779b97f4a7c15)
	rnd := func() float64 {
		seed ^= seed << 13
		seed ^= seed >> 7
		seed ^= seed << 17
		// 映射到 [-1,1)
		return float64(int64(seed)>>11) / float64(1<<53) * 2
	}
	for i := 0; i < nSNP; i++ {
		pairs := [][2]string{{"A", "G"}, {"C", "T"}, {"A", "C"}, {"G", "T"}}
		p := pairs[i%4]
		id := idFor(i)
		l2 := 1.0 + float64(i%7)
		eps := rnd() * 0.15
		chi2 := 1 + slope*l2 + eps
		if chi2 <= 0 {
			chi2 = 0.5
		}
		z := math.Sqrt(chi2)
		se := 0.05
		lv, _ := variant.New(id, "1", int64(100+i), p[0], p[1], math.NaN(), variant.StrandUnknown)
		ldRecs = append(ldRecs, ingest.LDScoreRecord{V: lv, L2: l2, N: nBar})
		sv, _ := variant.New(id, "1", int64(100+i), p[0], p[1], math.NaN(), variant.StrandUnknown)
		sumRecs = append(sumRecs, ingest.SummaryRecord{
			V: sv, Beta: z * se, SE: se, N: nBar, Z: z, HasZ: true})
	}
	return &ingest.SummaryTable{Records: sumRecs},
		&ingest.LDScoreTable{Records: ldRecs}
}

// TestH2BackendsAgree 双后端在同一合成数据上的 h2/SE/CI/标志必须一致：
// 数值在声明容差内，Flags/拒绝明细完全一致。
func TestH2BackendsAgree(t *testing.T) {
	sum, ld := noisyH2Tables(t, 120, 300, 0.4)
	rs, err := EstimateH2(context.Background(), sum, ld, H2Options{
		NBar: 300, Boundary: H2ReportRaw, Workers: 1, Backend: numwrap.BackendStd})
	if err != nil {
		t.Fatal(err)
	}
	rg, err := EstimateH2(context.Background(), sum, ld, H2Options{
		NBar: 300, Boundary: H2ReportRaw, Workers: 1, Backend: numwrap.BackendGonum})
	if err != nil {
		t.Fatal(err)
	}
	if err := compareH2Results(rs, rg); err != nil {
		t.Fatal(err)
	}
	// 同一后端重复跑必须逐位一致。
	rs2, _ := EstimateH2(context.Background(), sum, ld, H2Options{
		NBar: 300, Boundary: H2ReportRaw, Workers: 4, Backend: numwrap.BackendStd})
	if rs.H2 != rs2.H2 || rs.SE != rs2.SE {
		t.Fatal("h2 rerun / parallel not bit-exact")
	}
	t.Logf("[judgement] h2 stdlib=%.6f(%.6f) gonum=%.6f(%.6f) flags=%s rejected=%d",
		rs.H2, rs.SE, rg.H2, rg.SE, rs.Flags, len(rs.Rejected))

	// crosscheck 后端必须成功且判定字段与 stdlib 相同。
	rc, err := EstimateH2(context.Background(), sum, ld, H2Options{
		NBar: 300, Boundary: H2ReportRaw, Workers: 3, Backend: numwrap.BackendCrosscheck})
	if err != nil {
		t.Fatal(err)
	}
	if rc.H2 != rs.H2 || rc.Flags != rs.Flags || rc.NSnps != rs.NSnps {
		t.Fatalf("crosscheck result must equal stdlib reference: %+v vs %+v", rc, rs)
	}

	// 错误路径的判定在两后端也必须一致：clamp 边界策略不改变结论结构。
	rcg, err := EstimateH2(context.Background(), sum, ld, H2Options{
		NBar: 300, Boundary: H2Clamp, Workers: 1, Backend: numwrap.BackendGonum})
	if err != nil {
		t.Fatal(err)
	}
	rcs, err := EstimateH2(context.Background(), sum, ld, H2Options{
		NBar: 300, Boundary: H2Clamp, Workers: 1, Backend: numwrap.BackendStd})
	if err != nil {
		t.Fatal(err)
	}
	if err := compareH2Results(rcs, rcg); err != nil {
		t.Fatal(err)
	}
}

// TestH2CrosscheckRejectsUnknownBackend 未知后端显式拒绝。
func TestH2CrosscheckRejectsUnknownBackend(t *testing.T) {
	sum, ld := noisyH2Tables(t, 60, 100, 0.2)
	_, err := EstimateH2(context.Background(), sum, ld, H2Options{
		NBar: 100, Backend: numwrap.Backend("lapack")})
	if err == nil {
		t.Fatal("unknown backend must be rejected, never fall back")
	}
}
