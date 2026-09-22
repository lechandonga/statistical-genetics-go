// Package ingest 负责把本地文件（TSV/简单文本矩阵）解析并归一化为统一变体模型。
//
// 任何字段缺失、非法数值、重复变体、无法识别的等位基因都会以 variant.RejectError
// 逐条记录在 ParseReport 中，绝不静默跳过。
package ingest

import (
	"errors"
	"fmt"
	"time"

	"github.com/lechandonga/statistical-genetics-go/variant"
)

// LimitPolicy 决定超过规模上限时的行为。
type LimitPolicy uint8

const (
	// LimitReject 超限即整体拒绝（默认）：不返回任何部分结果。
	LimitReject LimitPolicy = iota
	// LimitTruncate 截断到上限：ParseReport.Truncated 记录截断事实与位置。
	LimitTruncate
)

// errTruncated 是扫描器在截断策略下的正常结束信号（非错误）。
var errTruncated = errors.New("truncated at configured limit")

// Limits 是单次解析的规模与耗时上限。零值字段表示使用 DefaultLimits。
type Limits struct {
	MaxVariants  int           // 变体（行）数上限
	MaxSamples   int           // 基因型样本数上限
	MaxLineBytes int64         // 单行字节上限（防止畸形巨行）
	MaxFileBytes int64         // 单文件字节上限
	MaxElapsed   time.Duration // 解析墙钟耗时上限
	Policy       LimitPolicy
}

// DefaultLimits 返回推荐上限。
func DefaultLimits() Limits {
	return Limits{
		MaxVariants:  5_000_000,
		MaxSamples:   100_000,
		MaxLineBytes: 1 << 20, // 1 MiB
		MaxFileBytes: 1 << 31, // 2 GiB
		MaxElapsed:   5 * time.Minute,
		Policy:       LimitReject,
	}
}

// RowRef 定位一条被拒绝的输入行。
type RowRef struct {
	Line   int
	ID     string
	Reason variant.Reason
	Detail string
}

func (r RowRef) String() string {
	return fmt.Sprintf("line=%d id=%q reason=%s %s", r.Line, r.ID, r.Reason, r.Detail)
}

// ParseReport 汇总一次解析的拒绝、截断与输入摘要。
type ParseReport struct {
	Source     string
	TotalRows  int // 物理数据行数（不含表头；截断时为实际上游行数）
	Accepted   int
	Rejected   []RowRef
	Truncated  bool
	TruncateAt int      // Policy==LimitTruncate 时保留的行数上限
	Columns    []string // 识别到的列
	Summary    string   // 输入摘要（测试日志用）
}

// RejectError 将 RowRef 升级为 *variant.RejectError。
func (r RowRef) RejectError(source string) *variant.RejectError {
	return &variant.RejectError{Source: source, Line: r.Line, ID: r.ID, Reason: r.Reason, Detail: r.Detail}
}

// isLimitError / asLimitError 帮助各解析器区分超限与普通 IO 错误。
func asLimitError(err error) (*limitExceeded, bool) {
	var le *limitExceeded
	if errors.As(err, &le) {
		return le, true
	}
	return nil, false
}
