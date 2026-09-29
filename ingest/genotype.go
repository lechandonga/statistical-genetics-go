package ingest

import (
	"errors"
	"fmt"
	"io"
	"math"
	"strings"

	"github.com/lechandonga/statistical-genetics-go/variant"
)

// GTEncoding 声明基因型矩阵的编码方式。
type GTEncoding uint8

const (
	// GTDosage012：每格为 0/1/2 的效应等位基因剂量（EA 拷贝数）。
	GTDosage012 GTEncoding = iota
	// GTACGT：每格为 "AA"/"AT" 二联体，按该行 EA/NEA 判定 EA 拷贝数。
	GTACGT
	// GTSpaceAlleles：每格为 "A A" 形式的两等位基因。
	GTSpaceAlleles
)

// GenotypeOptions 控制基因型矩阵解析。
type GenotypeOptions struct {
	Source    string
	Delimiter byte
	Encoding  GTEncoding
	// AnnotationColumns 指定前导注释列名（按表头名识别样本列起点）。
	// 默认 {"SNP","CHR","BP","A1","A2"}。
	AnnotationColumns []string
	Limits            Limits
	// MissingTokens 标记缺失基因型的字符串集合，默认 NA/NaN/./?/空
	MissingTokens map[string]struct{}
	// Lenient=true 返回已接受行与报告；默认严格。
	Lenient bool
	// SampleColumns 显式声明样本列名集合；为空时取“注释列之后的全部列”。
	SampleColumns []string
}

// GenotypeMatrix 是个体 × 变体的本地基因型矩阵。
//
// 存储：Dosage[i*nSamples + j]（变体 i、样本 j），缺失为 NaN 且 Missing 同位置 true；
// Variants 顺序与文件行顺序一致。
type GenotypeMatrix struct {
	Variants  []*variant.Variant
	SampleIDs []string
	Dosage    []float64
	Missing   []bool
	Report    ParseReport
}

// NVariants / NSamples 返回维度。
func (m *GenotypeMatrix) NVariants() int {
	if m == nil {
		return 0
	}
	return len(m.Variants)
}

func (m *GenotypeMatrix) NSamples() int {
	if m == nil {
		return 0
	}
	return len(m.SampleIDs)
}

// At 返回 (variant i, sample j) 的剂量与缺失标记。
func (m *GenotypeMatrix) At(i, j int) (float64, bool) {
	idx := i*m.NSamples() + j
	return m.Dosage[idx], m.Missing[idx]
}

// ParseGenotype 解析“变体为行、样本为列”的本地基因型矩阵。
func ParseGenotype(r io.Reader, opts GenotypeOptions) (*GenotypeMatrix, error) {
	if opts.Limits == (Limits{}) {
		opts.Limits = DefaultLimits()
	}
	delim := opts.Delimiter
	if delim == 0 {
		delim = '\t'
	}
	missing := opts.MissingTokens
	if missing == nil {
		missing = map[string]struct{}{"NA": {}, "NAN": {}, ".": {}, "?": {}, "": {}, "-": {}}
	}
	ann := opts.AnnotationColumns
	if ann == nil {
		ann = []string{"SNP", "CHR", "CHROM", "BP", "POS", "A1", "A2", "EA", "NEA", "STRAND"}
	}
	annSet := map[string]struct{}{}
	for _, a := range ann {
		annSet[normalizeHeaderName(a)] = struct{}{}
	}

	report := ParseReport{Source: opts.Source}
	sc := newRowScanner(r, opts.Limits, &report, delim)
	defer sc.finish()
	if err := sc.readHeader(); err != nil {
		return nil, fmt.Errorf("source=%s: %w", opts.Source, err)
	}

	idC, idOk := sc.header.lookup("snp", "snpid", "id", "rsid")
	eaC, eaOk := sc.header.lookup("a1", "ea", "alt")
	neaC, neaOk := sc.header.lookup("a2", "nea", "ref")
	if !idOk || !eaOk || !neaOk {
		return nil, &variant.RejectError{Source: opts.Source, Reason: variant.ReasonMissingField,
			Detail: "genotype matrix requires SNP, A1, A2 annotation columns"}
	}
	chromC, _ := sc.header.lookup("chrom", "chr")
	posC, _ := sc.header.lookup("bp", "pos")

	// 样本列 = 非注释列（或显式集合）
	var sampleIdx []int
	if len(opts.SampleColumns) > 0 {
		for _, name := range opts.SampleColumns {
			i, ok := sc.header.lookup(name)
			if !ok {
				return nil, &variant.RejectError{Source: opts.Source,
					Reason: variant.ReasonMissingField, Detail: "sample column not found: " + name}
			}
			sampleIdx = append(sampleIdx, i)
		}
	} else {
		for i, name := range sc.rawCols {
			if _, isAnn := annSet[normalizeHeaderName(name)]; !isAnn {
				sampleIdx = append(sampleIdx, i)
			}
		}
	}
	if len(sampleIdx) == 0 {
		return nil, &variant.RejectError{Source: opts.Source, Reason: variant.ReasonMissingField,
			Detail: "no sample columns detected after annotations"}
	}
	report.Usage.SamplesOffered = len(sampleIdx)
	if opts.Limits.MaxSamples > 0 && len(sampleIdx) > opts.Limits.MaxSamples {
		if opts.Limits.Policy == LimitTruncate {
			// 按表头列顺序保留前 MaxSamples 个样本列，确定性截断。
			sampleIdx = sampleIdx[:opts.Limits.MaxSamples]
			report.Truncated = true
			report.TruncateAt = opts.Limits.MaxSamples
		} else {
			return nil, &limitExceeded{what: fmt.Sprintf("sample count %d > MaxSamples %d",
				len(sampleIdx), opts.Limits.MaxSamples)}
		}
	}

	out := &GenotypeMatrix{}
	for _, i := range sampleIdx {
		out.SampleIDs = append(out.SampleIDs, sc.rawCols[i])
	}
	nS := len(sampleIdx)
	report.Usage.Samples = nS
	seen := map[variant.MatchKey]int{}

	for {
		fields, line, rerr := sc.next()
		if rerr != nil {
			if errors.Is(rerr, io.EOF) || errors.Is(rerr, errTruncated) {
				break
			}
			if le, ok := asLimitError(rerr); ok {
				return nil, fmt.Errorf("source=%s: %w", opts.Source, le)
			}
			return nil, rerr
		}
		v, perr := parseGTVariant(fields, line, idC, chromC, posC, eaC, neaC)
		if perr != nil {
			report.Rejected = append(report.Rejected, RowRef{
				Line: line, ID: perr.id, Reason: perr.reason, Detail: perr.detail})
			continue
		}
		if first, dup := seen[v.Key()]; dup {
			report.Rejected = append(report.Rejected, RowRef{
				Line: line, ID: v.ID, Reason: variant.ReasonDuplicateVariant,
				Detail: fmt.Sprintf("key=%s already at line %d", v.Key(), first)})
			continue
		}
		seen[v.Key()] = line

		base := len(out.Variants) * nS
		rowDose := make([]float64, nS)
		rowMiss := make([]bool, nS)
		bad := false
		for j, ci := range sampleIdx {
			raw := strings.ToUpper(strings.TrimSpace(field(fields, ci)))
			if _, isMiss := missing[raw]; isMiss {
				rowDose[j] = math.NaN()
				rowMiss[j] = true
				continue
			}
			d, err := decodeDosage(raw, v, opts.Encoding, line)
			if err != nil {
				report.Rejected = append(report.Rejected, RowRef{
					Line: line, ID: v.ID, Reason: variant.ReasonInvalidNumeric,
					Detail: fmt.Sprintf("sample=%s genotype=%q: %s", out.SampleIDs[j], raw, err.Error())})
				bad = true
				break
			}
			rowDose[j] = d
			rowMiss[j] = false
		}
		if bad {
			continue
		}
		out.Variants = append(out.Variants, v)
		out.Dosage = append(out.Dosage, rowDose...)
		out.Missing = append(out.Missing, rowMiss...)
		_ = base
	}

	sc.finish()
	out.Report = report
	out.Report.Accepted = len(out.Variants)
	out.Report.Usage.Samples = nS
	out.Report.Usage.SamplesOffered = report.Usage.SamplesOffered
	out.Report.Summary = summarize(out.Report)
	if !opts.Lenient {
		if err := out.Report.Err(); err != nil {
			return nil, err
		}
	}
	return out, nil
}

func parseGTVariant(fields [][]byte, line, idC, chromC, posC, eaC, neaC int) (*variant.Variant, *rowError) {
	id := field(fields, idC)
	chrom := ""
	if chromC >= 0 {
		chrom = field(fields, chromC)
	}
	var pos int64
	if posC >= 0 {
		if pv, ok, err := parseFloatField(fields, posC, "BP", line, false); err != nil {
			return nil, asRow(err)
		} else if ok {
			pos = int64(pv)
		}
	}
	v, err := variant.New(id, chrom, pos, field(fields, eaC), field(fields, neaC),
		math.NaN(), variant.StrandUnknown)
	if err != nil {
		if rr, ok := variant.AsReject(err); ok {
			return nil, &rowError{id: rr.ID, line: line, reason: rr.Reason, detail: rr.Detail}
		}
		return nil, &rowError{id: id, line: line, reason: variant.ReasonInvalidAllele, detail: err.Error()}
	}
	return v, nil
}

// decodeDosage 把单格基因型解码为效应等位基因（v.EA）拷贝数 0/1/2。
func decodeDosage(raw string, v *variant.Variant, enc GTEncoding, line int) (float64, error) {
	switch enc {
	case GTDosage012:
		f, err := parseFloatStrict(raw)
		if err != nil {
			return 0, err
		}
		if f != 0 && f != 1 && f != 2 {
			return 0, fmt.Errorf("dosage must be 0/1/2, got %s", raw)
		}
		return f, nil
	case GTACGT:
		if len(raw) != 2 {
			return 0, fmt.Errorf("ACGT genotype must be 2 bases, got %q", raw)
		}
		return countEA(raw[0], raw[1], v)
	case GTSpaceAlleles:
		parts := strings.Fields(raw)
		if len(parts) != 2 || len(parts[0]) != 1 || len(parts[1]) != 1 {
			return 0, fmt.Errorf("space-allele genotype must be 'B B', got %q", raw)
		}
		return countEA(parts[0][0], parts[1][0], v)
	default:
		return 0, fmt.Errorf("unsupported genotype encoding %d", enc)
	}
}

func countEA(b1, b2 byte, v *variant.Variant) (float64, error) {
	ea, nea := byte(v.EA[0]), byte(v.NEA[0])
	count := 0.0
	for _, b := range []byte{b1, b2} {
		switch b {
		case ea:
			count++
		case nea:
			// ok
		default:
			return 0, fmt.Errorf("allele %c not in {%s,%s}", b, v.EA, v.NEA)
		}
	}
	return count, nil
}

func parseFloatStrict(raw string) (float64, error) {
	var f float64
	if _, err := fmt.Sscanf(raw, "%g", &f); err != nil {
		return 0, fmt.Errorf("not numeric: %q", raw)
	}
	if f != f {
		return 0, fmt.Errorf("NaN dosage")
	}
	return f, nil
}
