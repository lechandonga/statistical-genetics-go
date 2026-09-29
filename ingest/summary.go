package ingest

import (
	"errors"
	"fmt"
	"io"
	"math"
	"strings"

	"github.com/lechandonga/statistical-genetics-go/variant"
)

// EffectField 声明输入汇总统计使用的效应量字段。
type EffectField uint8

const (
	EffectBeta EffectField = iota // beta（线性/对数尺度效应）
	EffectOR                      // 比值比 OR（内部转换为 logOR）
	EffectZ                       // 仅提供 Z 分数（退化模式，无 SE 反演）
)

// String 返回效应字段类型的 CLI 可读名称。
func (e EffectField) String() string {
	switch e {
	case EffectBeta:
		return "beta"
	case EffectOR:
		return "or"
	case EffectZ:
		return "z"
	default:
		return "unknown"
	}
}

// ColumnMap 定义列名覆盖；空字符串使用内置默认别名集合。
type ColumnMap struct {
	ID     string
	Chrom  string
	Pos    string
	EA     string
	NEA    string
	Effect string // beta / or / z 列名
	SE     string
	P      string
	N      string
	EAF    string
	Strand string
}

// SummaryOptions 控制汇总统计解析。
type SummaryOptions struct {
	Source    string
	Delimiter byte // 0 => tab
	Effect    EffectField
	Columns   ColumnMap
	FallbackN float64 // N 列缺失时使用的统一样本量；仍缺失则拒绝
	Strand    variant.Strand
	Limits    Limits
	// Lenient=true 时即使存在被拒绝行也返回已接受记录与报告；
	// 默认（false）为严格模式：任何拒绝都导致整体失败、不返回部分结果。
	Lenient bool
}

// SummaryRecord 是一行归一化后的汇总统计。
type SummaryRecord struct {
	V    *variant.Variant
	Beta float64 // 统一为对数尺度效应（OR 输入已取 log）
	SE   float64
	P    float64 // 缺失为 NaN
	N    float64
	Z    float64
	HasZ bool
}

// SummaryTable 是解析结果：记录顺序与文件行顺序严格一致（确定性要求）。
type SummaryTable struct {
	Records []SummaryRecord
	Report  ParseReport
}

// ParseSummary 从流式读取器解析全基因组汇总统计。
//
// 必填列：ID（或 CHROM+POS）、EA、NEA、效应列；
// beta/or 模式下 SE 必填；N 列或 FallbackN 至少有一个可用。
// 重复变体（相同 MatchKey）第二行起以 ReasonDuplicateVariant 拒绝。
// Strict 模式下任何拒绝都导致整体失败且不返回部分结果。
func ParseSummary(r io.Reader, opts SummaryOptions) (*SummaryTable, error) {
	if opts.Limits == (Limits{}) {
		opts.Limits = DefaultLimits()
	}
	delim := opts.Delimiter
	if delim == 0 {
		delim = '\t'
	}
	report := ParseReport{Source: opts.Source}
	sc := newRowScanner(r, opts.Limits, &report, delim)
	if err := sc.readHeader(); err != nil {
		return nil, fmt.Errorf("source=%s: %w", opts.Source, err)
	}

	col, err := resolveSummaryColumns(sc.header, opts)
	if err != nil {
		return nil, err
	}

	out := &SummaryTable{Report: report}
	seen := map[variant.MatchKey]int{} // key -> 已接受行号

	for {
		fields, line, rerr := sc.next()
		if rerr != nil {
			if errors.Is(rerr, io.EOF) || errors.Is(rerr, errTruncated) {
				break
			}
			if le, ok := asLimitError(rerr); ok {
				sc.finalize(&report, len(out.Records), 0)
				return nil, fmt.Errorf("source=%s: %w", opts.Source, le)
			}
			return nil, rerr
		}
		rec, perr := parseSummaryRow(fields, line, col, opts)
		if perr != nil {
			rr := RowRef{Line: perr.line, ID: perr.id, Reason: perr.reason, Detail: perr.detail}
			report.Rejected = append(report.Rejected, rr)
			continue
		}
		key := rec.V.Key()
		if first, dup := seen[key]; dup {
			report.Rejected = append(report.Rejected, RowRef{
				Line: line, ID: rec.V.ID, Reason: variant.ReasonDuplicateVariant,
				Detail: fmt.Sprintf("key=%s already accepted at line %d", key, first)})
			continue
		}
		seen[key] = line
		out.Records = append(out.Records, rec)
	}

	out.Report = report
	out.Report.Accepted = len(out.Records)
	sc.finalize(&out.Report, len(out.Records), 0)
	out.Report.Summary = summarize(out.Report)
	// 默认严格：存在拒绝即失败；opts.Strict==false 仅在显式构造时发生。
	// 注：零值 SummaryOptions{} 不合法（Limits 等必填），因此用
	// opts 来源是否设置 Lenient 标记；这里以字段语义“Strict=true 严格，
	// false 宽松”，ParseSummary 的公共约定见文档：默认严格。
	if !opts.Lenient {
		if err := out.Report.Err(); err != nil {
			return nil, err
		}
	}
	return out, nil
}

// summaryColumns 是解析后的列位置，-1 表示该列不存在。
type summaryColumns struct {
	id, chrom, pos, ea, nea, effect, se, p, n, eaf, strand int
}

func resolveSummaryColumns(h headerIndex, opts SummaryOptions) (summaryColumns, error) {
	c := summaryColumns{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
	pick := func(over string, aliases ...string) int {
		if over != "" {
			if i, ok := h.lookup(over); ok {
				return i
			}
			return -1
		}
		i, _ := h.lookup(aliases...)
		return i
	}
	c.id = pick(opts.Columns.ID, "snp", "snpid", "id", "variantid", "rsid", "rsids")
	c.chrom = pick(opts.Columns.Chrom, "chrom", "chr", "chromosome", "hg19chr", "hg38chr")
	c.pos = pick(opts.Columns.Pos, "pos", "bp", "position", "basepair")
	c.ea = pick(opts.Columns.EA, "a1", "ea", "effectallele", "alt", "minorallele")
	c.nea = pick(opts.Columns.NEA, "a2", "nea", "noneffectallele", "ref", "majorallele")
	c.eaf = pick(opts.Columns.EAF, "eaf", "freq", "frq", "a1freq", "maf", "freq1")
	c.strand = pick(opts.Columns.Strand, "strand")
	switch opts.Effect {
	case EffectBeta:
		c.effect = pick(opts.Columns.Effect, "beta", "effect", "effectsize", "b")
	case EffectOR:
		c.effect = pick(opts.Columns.Effect, "or", "oddsratio")
	case EffectZ:
		c.effect = pick(opts.Columns.Effect, "z", "zscore")
	default:
		return c, fmt.Errorf("source=%s: unsupported effect field kind %d", opts.Source, opts.Effect)
	}
	c.se = pick(opts.Columns.SE, "se", "sebeta", "standarderror", "stderr")
	c.p = pick(opts.Columns.P, "p", "pval", "pvalue", "p.value")
	c.n = pick(opts.Columns.N, "n", "n_samples", "samplesize", "neff", "weight")

	var missing []string
	if c.id < 0 && !(c.chrom >= 0 && c.pos >= 0) {
		missing = append(missing, "id(or chrom+pos)")
	}
	if c.ea < 0 {
		missing = append(missing, "effect allele (A1/EA/ALT)")
	}
	if c.nea < 0 {
		missing = append(missing, "non-effect allele (A2/NEA/REF)")
	}
	if c.effect < 0 {
		missing = append(missing, effectColName(opts.Effect))
	}
	if opts.Effect != EffectZ && c.se < 0 {
		missing = append(missing, "se")
	}
	if len(missing) > 0 {
		return c, &variant.RejectError{Source: opts.Source, Line: 1,
			Reason: variant.ReasonMissingField,
			Detail: "missing required columns: " + joinComma(missing)}
	}
	return c, nil
}

func effectColName(e EffectField) string {
	switch e {
	case EffectBeta:
		return "beta"
	case EffectOR:
		return "or"
	case EffectZ:
		return "z"
	default:
		return "?"
	}
}

func joinComma(xs []string) string {
	out := ""
	for i, x := range xs {
		if i > 0 {
			out += ", "
		}
		out += x
	}
	return out
}

func parseSummaryRow(fields [][]byte, line int, c summaryColumns, opts SummaryOptions) (SummaryRecord, *rowError) {
	id := field(fields, c.id)
	eaRaw := field(fields, c.ea)
	neaRaw := field(fields, c.nea)

	var chromRaw string
	var posV int64
	if c.chrom >= 0 || c.pos >= 0 {
		chromRaw = field(fields, c.chrom)
		if c.pos >= 0 {
			v, ok, err := parseFloatField(fields, c.pos, "POS", line, false)
			if err != nil {
				return SummaryRecord{}, asRow(err)
			}
			if ok {
				if v <= 0 || v != math.Trunc(v) {
					return SummaryRecord{}, &rowError{id: id, line: line,
						reason: variant.ReasonInvalidPosition,
						detail: fmt.Sprintf("position must be positive integer: %v", v)}
				}
				posV = int64(v)
			}
		}
	}

	var eaf = math.NaN()
	if c.eaf >= 0 {
		v, ok, err := parseFloatField(fields, c.eaf, "EAF", line, false)
		if err != nil {
			return SummaryRecord{}, asRow(err)
		}
		if ok {
			if v < 0 || v > 1 {
				return SummaryRecord{}, &rowError{id: id, line: line,
					reason: variant.ReasonInvalidFrequency,
					detail: fmt.Sprintf("EAF out of [0,1]: %v", v)}
			}
			eaf = v
		}
	}

	strand := opts.Strand
	if c.strand >= 0 {
		st := strings.ToUpper(strings.TrimSpace(field(fields, c.strand)))
		switch st {
		case "", "+":
			// 默认前向
		case "-":
			strand = variant.StrandReverse
		default:
			return SummaryRecord{}, &rowError{id: id, line: line,
				reason: variant.ReasonUnresolvableStrand,
				detail: "unrecognized strand marker: " + st}
		}
	}

	v, verr := variant.New(id, chromRaw, posV, eaRaw, neaRaw, eaf, strand)
	if verr != nil {
		if rr, ok := variant.AsReject(verr); ok {
			return SummaryRecord{}, &rowError{id: rr.ID, line: line, reason: rr.Reason, detail: rr.Detail}
		}
		return SummaryRecord{}, &rowError{id: id, line: line, reason: variant.ReasonInvalidAllele, detail: verr.Error()}
	}

	rec := SummaryRecord{V: v}

	// 效应量
	ev, _, err := parseFloatField(fields, c.effect, effectColName(opts.Effect), line, true)
	if err != nil {
		return SummaryRecord{}, asRow(err)
	}
	switch opts.Effect {
	case EffectBeta:
		rec.Beta = ev
	case EffectOR:
		if ev <= 0 {
			return SummaryRecord{}, &rowError{id: id, line: line,
				reason: variant.ReasonInvalidNumeric, detail: fmt.Sprintf("OR must be positive: %v", ev)}
		}
		rec.Beta = math.Log(ev)
	case EffectZ:
		rec.Z = ev
		rec.HasZ = true
	}

	if opts.Effect != EffectZ {
		sev, _, serr := parseFloatField(fields, c.se, "SE", line, true)
		if serr != nil {
			return SummaryRecord{}, asRow(serr)
		}
		if sev <= 0 {
			return SummaryRecord{}, &rowError{id: id, line: line,
				reason: variant.ReasonInvalidNumeric, detail: fmt.Sprintf("SE must be positive: %v", sev)}
		}
		rec.SE = sev
		rec.Z = rec.Beta / sev
		rec.HasZ = true
	} else {
		rec.SE = math.NaN()
	}

	if c.p >= 0 {
		pv, ok, perr := parseFloatField(fields, c.p, "P", line, false)
		if perr != nil {
			return SummaryRecord{}, asRow(perr)
		}
		if ok {
			if pv < 0 || pv > 1 {
				return SummaryRecord{}, &rowError{id: id, line: line,
					reason: variant.ReasonInvalidNumeric, detail: fmt.Sprintf("P out of [0,1]: %v", pv)}
			}
			rec.P = pv
		} else {
			rec.P = math.NaN()
		}
	} else {
		rec.P = math.NaN()
	}

	if c.n >= 0 {
		nv, ok, nerr := parseFloatField(fields, c.n, "N", line, false)
		if nerr != nil {
			return SummaryRecord{}, asRow(nerr)
		}
		if ok {
			if nv <= 0 {
				return SummaryRecord{}, &rowError{id: id, line: line,
					reason: variant.ReasonInvalidNumeric, detail: fmt.Sprintf("N must be positive: %v", nv)}
			}
			rec.N = nv
		}
	}
	if rec.N == 0 {
		if opts.FallbackN > 0 {
			rec.N = opts.FallbackN
		} else if opts.Effect != EffectZ {
			return SummaryRecord{}, &rowError{id: id, line: line,
				reason: variant.ReasonMissingField, detail: "sample size N missing and no FallbackN configured"}
		}
		// Z 模式允许缺 N：该行在元分析里以 missing_sample_size 原因按变体排除，
		// 不拖垮整个文件；h2/PRS 仍会在各自入口要求可用 N。
	}

	return rec, nil
}

func asRow(err error) *rowError {
	if re, ok := err.(*rowError); ok {
		return re
	}
	return &rowError{line: 0, reason: variant.ReasonInvalidNumeric, detail: err.Error()}
}
