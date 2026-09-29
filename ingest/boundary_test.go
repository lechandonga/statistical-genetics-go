package ingest

import (
	"strings"
	"testing"
)

// TestMaxFileBytesBoundary：上限恰好等于文件大小时必须接受，
// 上限只大一个字节也必须接受；超限才拒绝（防止 off-by-one）。
func TestMaxFileBytesBoundary(t *testing.T) {
	in := "SNP\tA1\tA2\tBETA\tSE\tN\nrs1\tA\tG\t0.1\t0.1\t10\n"
	size := int64(len(in))
	for _, delta := range []int64{0, 1, 100} {
		op := optsBeta()
		op.FallbackN = 0
		op.Limits.MaxFileBytes = size + delta
		tab, err := ParseSummary(strings.NewReader(in), op)
		if err != nil {
			t.Fatalf("delta=%d size=%d: unexpected error: %v", delta, size, err)
		}
		if tab.Report.Usage.BytesRead != size {
			t.Fatalf("bytes=%d want %d", tab.Report.Usage.BytesRead, size)
		}
	}
	op := optsBeta()
	op.FallbackN = 0
	op.Limits.MaxFileBytes = size - 1
	if _, err := ParseSummary(strings.NewReader(in), op); err == nil {
		t.Fatal("one byte over the cap must reject")
	}
	t.Logf("[judgement] MaxFileBytes boundary exact/under/over correct (size=%d)", size)
}
