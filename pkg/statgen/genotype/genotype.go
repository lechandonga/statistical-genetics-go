// Package genotype parses local genotype matrices (0/1/2 effect-allele
// dosages with a documented missing code), provides missingness-aware
// per-sample genotype extraction and supports bounded parallel evaluation.
//
// Row layout: an optional header line names samples; the first column is a
// variant id "chr:pos:ref:effect". Variant ids are normalized through the
// canonical variant model so allele swaps/flips are visible, not implicit.
package genotype

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/internal/parseutil"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/variant"
)

// MissingCode marks a missing dosage call inside a Matrix.
const MissingCode = -1.0

// Matrix is a variant-major genotype matrix: Dosage[v][s].
type Matrix struct {
	Variants []variant.Variant
	Samples  []string
	Dosage   [][]float64 // rows == len(Variants), columns == len(Samples)
	Missing  [][]bool    // true where dosage was missing in the input
	Source   string
	Report   variant.Report
	index    map[variant.Key]int
}

// ParseConfig controls matrix parsing and limits.
type ParseConfig struct {
	Source       string
	Delim        string
	HasHeader    bool   // first non-comment line names the sample columns
	HasSampleCol bool   // first column of every row is a variant id
	MissingToken string // default "NA"
	MaxRows      int    // variant cap
	MaxSamples   int    // sample cap
	SkipReject   bool
}

const defaultMaxRows = 1_000_000

// Parse reads a delimited dosage matrix.
func Parse(r io.Reader, cfg ParseConfig) (*Matrix, error) {
	delim := cfg.Delim
	if delim == "" {
		delim = "\t"
	}
	miss := cfg.MissingToken
	if miss == "" {
		miss = "NA"
	}
	maxRows := cfg.MaxRows
	if maxRows == 0 {
		maxRows = defaultMaxRows
	}

	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 1024*1024), 32*1024*1024)

	lineNo := 0
	nextLine := func() (string, int, error) {
		for sc.Scan() {
			lineNo++
			t := sc.Text()
			if parseutil.IsBlank(t) || t[0] == '#' {
				continue
			}
			return t, lineNo, nil
		}
		if err := sc.Err(); err != nil {
			return "", lineNo, stgerr.Wrap(stgerr.KindMalformedInput, cfg.Source,
				"stream read failure", err)
		}
		return "", lineNo, io.EOF
	}

	m := &Matrix{Source: cfg.Source, index: map[variant.Key]int{},
		Report: variant.Report{Source: cfg.Source, RejectByKind: map[stgerr.Kind]int{}}}

	nSamples := -1
	if cfg.HasHeader {
		hdr, lno, err := nextLine()
		if err != nil {
			return nil, err
		}
		cols := parseutil.Split(hdr, delim)
		first := 0
		if cfg.HasSampleCol {
			first = 1
		}
		m.Samples = cols[first:]
		nSamples = len(m.Samples)
		if cfg.MaxSamples > 0 && nSamples > cfg.MaxSamples {
			return nil, stgerr.New(stgerr.KindLimitExceeded, cfg.Source,
				fmt.Sprintf("sample count %d exceeds limit %d (header line %d)",
					nSamples, cfg.MaxSamples, lno))
		}
	}

	for {
		text, lno, err := nextLine()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		if len(m.Variants) >= maxRows {
			if cfg.SkipReject {
				m.Report.Truncated = true
				break
			}
			return nil, stgerr.New(stgerr.KindLimitExceeded, cfg.Source,
				fmt.Sprintf("variant count exceeds limit %d", maxRows))
		}
		row := parseutil.Split(text, delim)
		idStart := 0
		if cfg.HasSampleCol {
			idStart = 1
		}
		var v variant.Variant
		if cfg.HasSampleCol {
			vv, verr := parseVariantID(row[0])
			if verr != nil {
				if !cfg.SkipReject {
					return nil, verr
				}
				m.Report.AddError(lno, row[0], verr)
				continue
			}
			v = vv
		}
		cells := row[idStart:]
		if nSamples < 0 {
			nSamples = len(cells)
			if cfg.MaxSamples > 0 && nSamples > cfg.MaxSamples {
				return nil, stgerr.New(stgerr.KindLimitExceeded, cfg.Source,
					fmt.Sprintf("sample count %d exceeds limit %d (line %d)",
						nSamples, cfg.MaxSamples, lno))
			}
			m.Samples = make([]string, nSamples)
			for i := range m.Samples {
				m.Samples[i] = fmt.Sprintf("sample_%d", i+1)
			}
		}
		if len(cells) != nSamples {
			verr := stgerr.New(stgerr.KindMalformedInput, cfg.Source,
				fmt.Sprintf("line %d has %d dosage fields, expected %d",
					lno, len(cells), nSamples))
			if !cfg.SkipReject {
				return nil, verr
			}
			m.Report.AddError(lno, row[0], verr)
			continue
		}

		dosage := make([]float64, nSamples)
		missing := make([]bool, nSamples)
		bad := false
		for i, cell := range cells {
			c := strings.TrimSpace(cell)
			if c == miss || strings.EqualFold(c, "NA") || c == "." || c == "-9" {
				dosage[i] = MissingCode
				missing[i] = true
				continue
			}
			x, ferr := strconv.ParseFloat(c, 64)
			if ferr != nil {
				verr := stgerr.New(stgerr.KindInvalidValue, cfg.Source,
					fmt.Sprintf("line %d column %d: unparsable dosage %q", lno, i+idStart+1, c))
				if !cfg.SkipReject {
					return nil, verr
				}
				m.Report.AddError(lno, row[0], verr)
				bad = true
				break
			}
			if verr := ValidateDosage(x, MissingCode); verr != nil {
				se := stgerr.AsError(verr)
				verr = stgerr.New(se.Kind, cfg.Source,
					fmt.Sprintf("line %d: %s", lno, se.Reason))
				if !cfg.SkipReject {
					return nil, verr
				}
				m.Report.AddError(lno, row[0], verr)
				bad = true
				break
			}
			dosage[i] = x
		}
		if bad {
			continue
		}

		k := v.Key()
		if _, dup := m.index[k]; dup {
			verr := stgerr.New(stgerr.KindDuplicateVariant, cfg.Source,
				"duplicate variant in matrix: "+string(k))
			if !cfg.SkipReject {
				return nil, verr
			}
			m.Report.AddError(lno, string(k), verr)
			continue
		}
		m.index[k] = len(m.Variants)
		m.Variants = append(m.Variants, v)
		m.Dosage = append(m.Dosage, dosage)
		m.Missing = append(m.Missing, missing)
		m.Report.Accepted++
	}

	return m, nil
}

// ParseFile reads and parses a named file, always closing the handle.
func ParseFile(path string, cfg ParseConfig) (*Matrix, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, stgerr.Wrap(stgerr.KindMalformedInput, cfg.Source,
			"cannot open "+path, err)
	}
	defer f.Close()
	return Parse(f, cfg)
}

// parseVariantID parses "chr:pos:ref:effect" into a canonical variant.
func parseVariantID(id string) (variant.Variant, error) {
	parts := strings.Split(strings.TrimSpace(id), ":")
	if len(parts) != 4 {
		return variant.Variant{}, stgerr.New(stgerr.KindMalformedInput, "genotype",
			"variant id must be chr:pos:ref:effect, got "+id)
	}
	chrom, err := variant.ParseChromosome(parts[0])
	if err != nil {
		return variant.Variant{}, err
	}
	pos, err := parseutil.Int(parts[1], "pos")
	if err != nil {
		return variant.Variant{}, err
	}
	return variant.New(chrom, pos, parts[2], parts[3])
}

// Row returns dosages and missing mask for one variant.
func (m *Matrix) Row(k variant.Key) ([]float64, []bool, bool) {
	i, ok := m.index[k]
	if !ok {
		return nil, nil, false
	}
	return m.Dosage[i], m.Missing[i], true
}

// VariantIndex returns the row index of k.
func (m *Matrix) VariantIndex(k variant.Key) (int, bool) {
	i, ok := m.index[k]
	return i, ok
}

// NumVariants returns the variant count.
func (m *Matrix) NumVariants() int { return len(m.Variants) }

// NumSamples returns the sample count.
func (m *Matrix) NumSamples() int { return len(m.Samples) }
