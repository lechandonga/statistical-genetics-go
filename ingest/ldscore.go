package ingest

import (
	"errors"
	"fmt"
	"io"
	"math"

	"github.com/lechandonga/statistical-genetics-go/variant"
)

// LDScoreRecord 是一行 LD 分数注释（LDSC 风格 l2 列）。
type LDScoreRecord struct {
	V    *variant.Variant
	L2   float64 // LD score: sum over SNPs of r^2（必须 >=0）
	MAF  float64 // 可选
	Chi2 float64 // 可选
	N    float64 // 可选
}

// LDScoreTable 保持与文件一致的记录顺序。
type LDScoreTable struct {
	Records []LDScoreRecord
	Report  ParseReport
}

// LDSOptions 控制 LD 分数解析（SummaryOptions 的精简版）。
type LDSOptions struct {
	Source    string
	Delimiter byte
	Strand    variant.Strand
	Limits    Limits
	Lenient   bool
	// FallbackN 为可选的逐 SNP 样本量缺省值。
	FallbackN float64
}

// ParseLDScores 解析 LD 分数文件。
//
// 必填列：SNP（或 CHR+BP）、A1、A2、L2；
// 可选列：MAF、CHISQ、N。重复变体按 ReasonDuplicateVariant 拒绝。
// 当以 SummaryOptions 调用时（旧入口），仅采用其通用字段。
func ParseLDScores(r io.Reader, opts LDSOptions) (*LDScoreTable, error) {
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

	col, err := resolveLDColumns(sc.header)
	if err != nil {
		return nil, err
	}

	out := &LDScoreTable{}
	seen := map[variant.MatchKey]int{}

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
		rec, perr := parseLDRow(fields, line, col, opts)
		if perr != nil {
			report.Rejected = append(report.Rejected, RowRef{
				Line: perr.line, ID: perr.id, Reason: perr.reason, Detail: perr.detail})
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
	if !opts.Lenient {
		if err := out.Report.Err(); err != nil {
			return nil, err
		}
	}
	return out, nil
}

type ldColumns struct {
	id, chrom, pos, ea, nea, l2, maf, chi2, n int
}

func resolveLDColumns(h headerIndex) (ldColumns, error) {
	c := ldColumns{-1, -1, -1, -1, -1, -1, -1, -1, -1}
	c.id, _ = h.lookup("snp", "snpid", "id", "rsid")
	c.chrom, _ = h.lookup("chrom", "chr", "chromosome")
	c.pos, _ = h.lookup("pos", "bp", "position")
	c.ea, _ = h.lookup("a1", "ea", "alt")
	c.nea, _ = h.lookup("a2", "nea", "ref")
	c.l2, _ = h.lookup("l2", "ldscore", "ld", "ldscore_sum")
	c.maf, _ = h.lookup("maf")
	c.chi2, _ = h.lookup("chisq", "chi2", "chi_sq")
	c.n, _ = h.lookup("n")

	var missing []string
	if c.id < 0 && !(c.chrom >= 0 && c.pos >= 0) {
		missing = append(missing, "id(or chrom+pos)")
	}
	if c.ea < 0 {
		missing = append(missing, "A1")
	}
	if c.nea < 0 {
		missing = append(missing, "A2")
	}
	if c.l2 < 0 {
		missing = append(missing, "L2")
	}
	if len(missing) > 0 {
		return c, &variant.RejectError{Reason: variant.ReasonMissingField,
			Detail: "ld score missing required columns: " + joinComma(missing)}
	}
	return c, nil
}

func parseLDRow(fields [][]byte, line int, c ldColumns, opts LDSOptions) (LDScoreRecord, *rowError) {
	id := field(fields, c.id)
	var chromRaw string
	var posV int64
	if c.chrom >= 0 {
		chromRaw = field(fields, c.chrom)
	}
	if c.pos >= 0 {
		if v, ok, err := parseFloatField(fields, c.pos, "BP", line, false); err != nil {
			return LDScoreRecord{}, asRow(err)
		} else if ok {
			if v <= 0 || v != math.Trunc(v) {
				return LDScoreRecord{}, &rowError{id: id, line: line,
					reason: variant.ReasonInvalidPosition, detail: "bad position"}
			}
			posV = int64(v)
		}
	}
	v, verr := variant.New(id, chromRaw, posV,
		field(fields, c.ea), field(fields, c.nea), math.NaN(), opts.Strand)
	if verr != nil {
		if rr, ok := variant.AsReject(verr); ok {
			return LDScoreRecord{}, &rowError{id: rr.ID, line: line, reason: rr.Reason, detail: rr.Detail}
		}
		return LDScoreRecord{}, &rowError{id: id, line: line, reason: variant.ReasonInvalidAllele, detail: verr.Error()}
	}
	rec := LDScoreRecord{V: v, MAF: math.NaN(), Chi2: math.NaN()}
	l2, _, err := parseFloatField(fields, c.l2, "L2", line, true)
	if err != nil {
		return LDScoreRecord{}, asRow(err)
	}
	if l2 < 0 {
		return LDScoreRecord{}, &rowError{id: id, line: line,
			reason: variant.ReasonInvalidNumeric, detail: fmt.Sprintf("L2 must be >=0: %v", l2)}
	}
	rec.L2 = l2
	if c.maf >= 0 {
		if mv, ok, err := parseFloatField(fields, c.maf, "MAF", line, false); err != nil {
			return LDScoreRecord{}, asRow(err)
		} else if ok {
			if mv < 0 || mv > 1 {
				return LDScoreRecord{}, &rowError{id: id, line: line,
					reason: variant.ReasonInvalidFrequency, detail: fmt.Sprintf("MAF out of [0,1]: %v", mv)}
			}
			rec.MAF = mv
		}
	}
	if c.chi2 >= 0 {
		if cv, ok, err := parseFloatField(fields, c.chi2, "CHISQ", line, false); err != nil {
			return LDScoreRecord{}, asRow(err)
		} else if ok {
			if cv < 0 {
				return LDScoreRecord{}, &rowError{id: id, line: line,
					reason: variant.ReasonInvalidNumeric, detail: "CHISQ negative"}
			}
			rec.Chi2 = cv
		}
	}
	if c.n >= 0 {
		if nv, ok, err := parseFloatField(fields, c.n, "N", line, false); err != nil {
			return LDScoreRecord{}, asRow(err)
		} else if ok {
			if nv <= 0 {
				return LDScoreRecord{}, &rowError{id: id, line: line,
					reason: variant.ReasonInvalidNumeric, detail: "N must be positive"}
			}
			rec.N = nv
		}
	}
	if rec.N == 0 && opts.FallbackN > 0 {
		rec.N = opts.FallbackN
	}
	return rec, nil
}
