package ingest

import (
	"bufio"
	"bytes"
	"errors"
	"io"
	"strconv"
	"strings"
	"time"

	"github.com/lechandonga/statistical-genetics-go/variant"
)

// headerIndex 把列名（大小写不敏感、忽略下划线/点/空格差异）映射到列序号。
type headerIndex map[string]int

func normalizeHeaderName(s string) string {
	var b strings.Builder
	for _, r := range strings.ToLower(strings.TrimSpace(s)) {
		switch r {
		case '_', '.', ' ', '-':
			// 跳过分隔差异
		default:
			b.WriteRune(r)
		}
	}
	return b.String()
}

func buildHeader(fields [][]byte) (headerIndex, []string) {
	idx := headerIndex{}
	raw := make([]string, len(fields))
	for i, f := range fields {
		name := normalizeHeaderName(string(f))
		raw[i] = strings.TrimSpace(string(f))
		if _, exists := idx[name]; !exists {
			idx[name] = i
		}
	}
	return idx, raw
}

// lookup 在表头中查找第一个匹配的别名。
func (h headerIndex) lookup(names ...string) (int, bool) {
	for _, n := range names {
		if i, ok := h[normalizeHeaderName(n)]; ok {
			return i, true
		}
	}
	return -1, false
}

// limitExceeded 表示解析在到达上限后被拒绝/截断。
type limitExceeded struct{ what string }

func (e *limitExceeded) Error() string { return "limit exceeded: " + e.what }

// rowScanner 封装带字节/耗时/行数上限的流式扫描。
type rowScanner struct {
	r        *bufio.Reader
	lim      Limits
	report   *ParseReport
	deadline time.Time
	started  time.Time
	line     int // 物理行号（表头为 1）
	header   headerIndex
	rawCols  []string
	delim    byte
	cr       *countingReader
}

// countingReader 统计实际消费字节；MaxFileBytes 在这里强制执行：
// 用 LimitReader 把可见数据截到 limit+1 字节——
// 若恰好消费了 limit+1，说明真实文件严格大于上限（超限）；
// 只消费到 <=limit 则为“恰好在上限内/恰好等于上限”，正常接受。
// 这样上限恰好等于文件大小时不会误报（无 off-by-one）。
type countingReader struct {
	r       io.Reader
	limit   int64
	read    int64
	allowed int64
}

func newCountingReader(r io.Reader, limit int64) *countingReader {
	allowed := int64(-1)
	if limit > 0 {
		allowed = limit + 1
	}
	return &countingReader{r: io.LimitReader(r, allowed), limit: limit, allowed: allowed}
}

func (c *countingReader) Read(p []byte) (int, error) {
	n, err := c.r.Read(p)
	c.read += int64(n)
	if c.limit > 0 && c.read > c.limit {
		return n, &limitExceeded{what: "file larger than MaxFileBytes"}
	}
	return n, err
}

// overLimit 在读取结束后报告文件是否严格大于上限。
func (c *countingReader) overLimit() bool { return c.limit > 0 && c.read > c.limit }

func newRowScanner(r io.Reader, lim Limits, report *ParseReport, delim byte) *rowScanner {
	def := DefaultLimits()
	if lim.MaxLineBytes <= 0 {
		lim.MaxLineBytes = def.MaxLineBytes
	}
	if lim.MaxFileBytes <= 0 {
		lim.MaxFileBytes = def.MaxFileBytes
	}
	if delim == 0 {
		delim = '\t'
	}
	var dl time.Time
	if lim.MaxElapsed > 0 {
		dl = time.Now().Add(lim.MaxElapsed)
	}
	cr := newCountingReader(r, lim.MaxFileBytes)
	return &rowScanner{r: bufio.NewReaderSize(cr, 64*1024), lim: lim,
		report: report, deadline: dl, delim: delim, started: time.Now(), cr: cr}
}

// readLine 读取一行（不含行尾），强制执行 MaxLineBytes 与 MaxFileBytes。
func (s *rowScanner) readLine() ([]byte, error) {
	var buf []byte
	for {
		chunk, prefix, err := s.r.ReadLine()
		buf = append(buf, chunk...)
		if int64(len(buf)) > s.lim.MaxLineBytes {
			return nil, &limitExceeded{what: "line longer than MaxLineBytes at line " + strconv.Itoa(s.line+1)}
		}
		if err != nil {
			if le, ok := err.(*limitExceeded); ok {
				s.report.Usage.LimitKind = le.what
			}
			return buf, err
		}
		if !prefix {
			return buf, nil
		}
	}
}

// readHeader 读取并解析表头行。
func (s *rowScanner) readHeader() error {
	row, err := s.readLine()
	s.line = 1
	if err != nil {
		if errors.Is(err, io.EOF) {
			return errors.New("empty input: header row expected")
		}
		return err
	}
	s.header, s.rawCols = buildHeader(bytes.Split(row, []byte{s.delim}))
	s.report.Columns = s.rawCols
	return nil
}

// next 读取下一数据行，返回字段切片与物理行号；io.EOF 表示结束。
func (s *rowScanner) next() ([][]byte, int, error) {
	if !s.deadline.IsZero() && time.Now().After(s.deadline) {
		s.report.Usage.LimitKind = "MaxElapsed"
		return nil, s.line, &limitExceeded{what: "MaxElapsed"}
	}
	line, err := s.readLine()
	if err != nil {
		return nil, s.line, err
	}
	s.line++
	// 跳过纯空行（不计入数据行）
	if len(strings.TrimSpace(string(line))) == 0 {
		return s.next()
	}
	s.report.TotalRows++
	s.report.Usage.Lines = s.report.TotalRows
	if s.lim.MaxVariants > 0 && s.report.TotalRows > s.lim.MaxVariants {
		if s.lim.Policy == LimitTruncate {
			s.report.Truncated = true
			s.report.TruncateAt = s.lim.MaxVariants
			s.report.Usage.Truncated = true
			s.report.Usage.LimitKind = "MaxVariants"
			return nil, s.line, errTruncated
		}
		s.report.Usage.LimitKind = "MaxVariants"
		return nil, s.line, &limitExceeded{what: "MaxVariants"}
	}
	return bytes.Split(line, []byte{s.delim}), s.line, nil
}

// finalize 把扫描器累计的真实用量写回报告（在解析成功/失败两条路径都可调用）。
func (s *rowScanner) finalize(report *ParseReport, accepted, samples int) {
	report.Usage.BytesRead = s.cr.read
	report.Usage.Lines = report.TotalRows
	report.Usage.Variants = accepted
	report.Usage.Samples = samples
	report.Usage.Elapsed = time.Since(s.started)
	if s.report.Usage.LimitKind == "" && s.cr.limit > 0 && s.cr.read >= s.cr.limit {
		report.Usage.LimitKind = "MaxFileBytes"
	}
}

// field 取字段字符串；越界返回空串。
func field(fields [][]byte, i int) string {
	if i < 0 || i >= len(fields) {
		return ""
	}
	return strings.TrimSpace(string(fields[i]))
}

// parseFloatField 解析数值字段。
//
//	required 且为空 -> ReasonMissingField
//	非 required 且为空 -> (0,false,nil)
//	无法解析或为 NaN/Inf -> ReasonInvalidNumeric
func parseFloatField(fields [][]byte, i int, col string, line int, required bool) (float64, bool, error) {
	raw := field(fields, i)
	if raw == "" {
		if required {
			return 0, false, &rowError{line: line, reason: variant.ReasonMissingField,
				detail: "required numeric column empty: " + col}
		}
		return 0, false, nil
	}
	v, err := strconv.ParseFloat(raw, 64)
	if err != nil {
		return 0, false, &rowError{line: line, reason: variant.ReasonInvalidNumeric,
			detail: "column " + col + " not finite-numeric: " + raw}
	}
	if nanOrInf(v) {
		return 0, false, &rowError{line: line, reason: variant.ReasonInvalidNumeric,
			detail: "column " + col + " is NaN/Inf: " + raw}
	}
	return v, true, nil
}

func nanOrInf(v float64) bool { return v != v || v > 1e308 || v < -1e308 }

// rowError 是解析期的轻量行错误，由各 Parse 函数升级为 variant.RejectError。
type rowError struct {
	id     string
	line   int
	reason variant.Reason
	detail string
}

func (e *rowError) Error() string { return e.detail }
