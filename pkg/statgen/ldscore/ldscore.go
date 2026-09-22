// Package ldscore parses per-variant linkage-disequilibrium score tables
// (the l2 / chi-square regression covariate) and normalizes them onto the
// canonical variant model.
package ldscore

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/internal/parseutil"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/variant"
)

// Scores are the LD-score annotations available per variant.
type Scores struct {
	Var       variant.Variant
	L2        float64 // sum of r^2 with all tagged variants
	MAF       float64 // optional minor allele frequency
	HasMAF    bool
	Weight    float64 // optional regression weight
	HasWeight bool
}

// Table is a normalized LD-score file.
type Table struct {
	Source    string
	Scores    map[variant.Key]Scores
	Report    variant.Report
	Truncated bool
}

// ColumnMap binds header names to canonical fields.
type ColumnMap struct {
	Chrom  string
	Pos    string
	Ref    string
	Effect string
	L2     string
	MAF    string
	Weight string
}

// ParseConfig controls parsing and limits.
type ParseConfig struct {
	Source     string
	Delim      string
	Columns    ColumnMap
	HasHeader  bool
	MaxRows    int
	MinL2      float64 // reject rows below this floor (default 0: non-negative)
	SkipReject bool
}

const defaultMaxRows = 2_000_000

// Parse reads a delimited LD-score stream.
func Parse(r io.Reader, cfg ParseConfig) (*Table, error) {
	delim := cfg.Delim
	if delim == "" {
		delim = "\t"
	}
	maxRows := cfg.MaxRows
	if maxRows == 0 {
		maxRows = defaultMaxRows
	}

	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 64*1024), 8*1024*1024)
	src := cfg.Source
	report := variant.Report{Source: src}
	out := make(map[variant.Key]Scores)
	cols := make(map[string]int)

	line := 0
	nextLine := func() (string, int, error) {
		for sc.Scan() {
			line++
			t := sc.Text()
			if parseutil.IsBlank(t) || t[0] == '#' {
				continue
			}
			return t, line, nil
		}
		if err := sc.Err(); err != nil {
			return "", line, stgerr.Wrap(stgerr.KindMalformedInput, src,
				"stream read failure", err)
		}
		return "", line, io.EOF
	}

	if cfg.HasHeader {
		hdr, _, err := nextLine()
		if err != nil {
			return nil, err
		}
		header := parseutil.Split(hdr, delim)
		aliases := map[string][]string{
			"chrom":  {"CHR", "CHROM"},
			"pos":    {"POS", "BP"},
			"ref":    {"REF", "A2", "NEA"},
			"effect": {"EFFECT", "ALT", "A1", "EA"},
			"l2":     {"L2", "SCORE", "LDSCORE", "LD_SCORE"},
			"maf":    {"MAF", "FREQ"},
			"weight": {"WEIGHT", "W", "L2_WEIGHT"},
		}
		custom := map[string]string{
			"chrom": cfg.Columns.Chrom, "pos": cfg.Columns.Pos,
			"ref": cfg.Columns.Ref, "effect": cfg.Columns.Effect,
			"l2": cfg.Columns.L2, "maf": cfg.Columns.MAF,
			"weight": cfg.Columns.Weight,
		}
		for f, names := range aliases {
			if c := custom[f]; c != "" {
				if i := parseutil.IndexOf(header, c); i >= 0 {
					cols[f] = i
				}
			} else if i := parseutil.IndexOf(header, names...); i >= 0 {
				cols[f] = i
			}
		}
	} else {
		for i, f := range []string{"chrom", "pos", "ref", "effect", "l2", "maf", "weight"} {
			cols[f] = i
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
		if lno-1 > maxRows {
			if cfg.SkipReject {
				report.Truncated = true
				break
			}
			return nil, stgerr.New(stgerr.KindLimitExceeded, src,
				fmt.Sprintf("more than %d data rows", maxRows))
		}
		row := parseutil.Split(text, delim)
		s, rerr := parseRow(row, cols, cfg)
		if rerr != nil {
			if !cfg.SkipReject {
				return nil, rerr
			}
			report.AddError(lno, keyHint(row, cols), rerr)
			continue
		}
		k := s.Var.Key()
		if _, dup := out[k]; dup {
			derr := stgerr.New(stgerr.KindDuplicateVariant, src,
				"duplicate variant "+string(k))
			if !cfg.SkipReject {
				return nil, derr
			}
			report.AddError(lno, string(k), derr)
			continue
		}
		out[k] = s
		report.Accepted++
	}

	return &Table{Source: src, Scores: out, Report: report, Truncated: report.Truncated}, nil
}

// ParseFile reads and parses a named file, always closing the handle.
func ParseFile(path string, cfg ParseConfig) (*Table, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, stgerr.Wrap(stgerr.KindMalformedInput, cfg.Source,
			"cannot open "+path, err)
	}
	defer f.Close()
	return Parse(f, cfg)
}

// Lookup returns scores for a variant key.
func (t *Table) Lookup(k variant.Key) (Scores, bool) {
	s, ok := t.Scores[k]
	return s, ok
}

func keyHint(row []string, cols map[string]int) string {
	get := func(f string) string {
		if i, ok := cols[f]; ok && i < len(row) {
			return strings.TrimSpace(row[i])
		}
		return ""
	}
	c, p, r, e := get("chrom"), get("pos"), get("ref"), get("effect")
	if c == "" && p == "" {
		return ""
	}
	return fmt.Sprintf("%s:%s:%s:%s", c, p, r, e)
}

func parseRow(row []string, cols map[string]int, cfg ParseConfig) (Scores, error) {
	var s Scores
	src := cfg.Source
	req := func(f string) (string, error) {
		i, ok := cols[f]
		if !ok {
			return "", stgerr.New(stgerr.KindMissingField, src, "no column bound to "+f)
		}
		return parseutil.Cell(row, i, f)
	}
	cs, err := req("chrom")
	if err != nil {
		return s, err
	}
	chrom, err := variant.ParseChromosome(cs)
	if err != nil {
		return s, err
	}
	ps, err := req("pos")
	if err != nil {
		return s, err
	}
	pos, err := parseutil.Int(ps, "pos")
	if err != nil {
		return s, err
	}
	rs, err := req("ref")
	if err != nil {
		return s, err
	}
	es, err := req("effect")
	if err != nil {
		return s, err
	}
	v, err := variant.New(chrom, pos, rs, es)
	if err != nil {
		return s, err
	}
	s.Var = v

	ls, err := req("l2")
	if err != nil {
		return s, err
	}
	l2, err := parseutil.Float(ls, "l2")
	if err != nil {
		return s, err
	}
	if l2 < cfg.MinL2 {
		return s, stgerr.New(stgerr.KindInvalidValue, src,
			fmt.Sprintf("L2=%.6g below floor %.6g", l2, cfg.MinL2))
	}
	s.L2 = l2

	if i, ok := cols["maf"]; ok && i < len(row) {
		if f := strings.TrimSpace(row[i]); f != "" && !strings.EqualFold(f, "NA") {
			maf, err := parseutil.Float(f, "maf")
			if err != nil {
				return s, err
			}
			if maf < 0 || maf > 1 {
				return s, stgerr.New(stgerr.KindInvalidValue, src,
					fmt.Sprintf("MAF must be in [0,1], got %v", maf))
			}
			s.MAF, s.HasMAF = maf, true
		}
	}
	if i, ok := cols["weight"]; ok && i < len(row) {
		if f := strings.TrimSpace(row[i]); f != "" && !strings.EqualFold(f, "NA") {
			w, err := parseutil.Float(f, "weight")
			if err != nil {
				return s, err
			}
			if w <= 0 {
				return s, stgerr.New(stgerr.KindInvalidValue, src,
					fmt.Sprintf("regression weight must be > 0, got %v", w))
			}
			s.Weight, s.HasWeight = w, true
		}
	}
	return s, nil
}
