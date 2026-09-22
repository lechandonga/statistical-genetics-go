// Package sumstats parses per-cohort genome-wide association summary
// statistics files and normalizes them onto the canonical variant model.
//
// Missing fields, illegal numerics, duplicate variants and unrecognizable
// alleles are never silently skipped: each rejected row appears in the
// returned variant.Report with a distinguishable stgerr.Kind. In strict
// mode (SkipReject=false) the first rejection aborts parsing; in tolerant
// mode all rejections are collected and parsing continues.
package sumstats

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

// VariantStats is one normalized association record.
type VariantStats struct {
	Var       variant.Variant
	Beta      float64 // effect size on log-odds / standardized scale
	SE        float64 // standard error of beta, strictly positive
	PValue    float64 // only valid when HasP
	HasP      bool
	EAF       float64 // effect-allele frequency; 0..1 when present
	HasEAF    bool
	N         float64 // per-variant sample size; 0 means "use cohort default"
	NCases    float64 // binary trait only
	NControls float64
}

// CohortSummary is one parsed cohort file.
type CohortSummary struct {
	Name     string
	Variants map[variant.Key]VariantStats
	DefaultN float64
	Report   variant.Report
}

// ColumnMap binds header names to canonical fields. Empty fields fall back
// to the built-in default header aliases.
type ColumnMap struct {
	Chrom     string
	Pos       string
	Ref       string
	Effect    string
	Beta      string
	SE        string
	PValue    string
	EAF       string
	N         string
	NCases    string
	NControls string
}

// ParseConfig controls parsing, validation and resource limits.
type ParseConfig struct {
	CohortName string
	Delim      string
	Columns    ColumnMap
	Comment    byte // lines beginning with this byte are skipped (default '#')
	HasHeader  bool
	DefaultN   float64 // used when rows carry no N column
	MaxRows    int     // hard row cap; 0 = unlimited unless budget applies
	SkipReject bool    // collect rejected rows instead of aborting
}

const defaultMaxRows = 2_000_000

// reader is the streaming line parser state.
type reader struct {
	sc     *bufio.Scanner
	cfg    ParseConfig
	report variant.Report
	out    map[variant.Key]VariantStats
	cols   map[string]int
	line   int
}

// Parse reads a delimited summary-statistics stream.
func Parse(r io.Reader, cfg ParseConfig) (*CohortSummary, error) {
	delim := cfg.Delim
	if delim == "" {
		delim = "\t"
	}
	if cfg.Comment == 0 {
		cfg.Comment = '#'
	}
	maxRows := cfg.MaxRows
	if maxRows == 0 {
		maxRows = defaultMaxRows
	}

	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 64*1024), 8*1024*1024)

	pr := &reader{
		sc:     sc,
		cfg:    cfg,
		report: variant.Report{Source: cfg.CohortName},
		out:    make(map[variant.Key]VariantStats),
		cols:   make(map[string]int),
	}

	if cfg.HasHeader {
		line, err := pr.nextDataLine()
		if err != nil {
			return nil, err
		}
		pr.resolveHeader(parseutil.Split(line, delim))
	} else {
		pr.positionalColumns()
	}

	for {
		line, err := pr.nextDataLine()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		if pr.line-1 > maxRows { // line 1 is the header
			if cfg.SkipReject {
				pr.report.Truncated = true
				break
			}
			return nil, stgerr.New(stgerr.KindLimitExceeded, cfg.CohortName,
				fmt.Sprintf("more than %d data rows", maxRows))
		}
		row := parseutil.Split(line, delim)
		st, rerr := pr.parseRow(row)
		if rerr != nil {
			if !cfg.SkipReject {
				return nil, rerr
			}
			pr.report.AddError(pr.line, pr.rowKeyHint(row), rerr)
			continue
		}
		k := st.Var.Key()
		if _, dup := pr.out[k]; dup {
			err := stgerr.New(stgerr.KindDuplicateVariant, cfg.CohortName,
				"duplicate variant "+string(k))
			if !cfg.SkipReject {
				return nil, err
			}
			pr.report.AddError(pr.line, string(k), err)
			continue
		}
		pr.out[k] = st
		pr.report.Accepted++
	}

	return &CohortSummary{
		Name:     cfg.CohortName,
		Variants: pr.out,
		DefaultN: cfg.DefaultN,
		Report:   pr.report,
	}, nil
}

// ParseFile reads and parses a named file. The handle is always closed,
// even when validation fails.
func ParseFile(path string, cfg ParseConfig) (*CohortSummary, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, stgerr.Wrap(stgerr.KindMalformedInput, cfg.CohortName,
			"cannot open "+path, err)
	}
	defer f.Close()
	return Parse(f, cfg)
}

// Lookup returns stats for a variant key.
func (c *CohortSummary) Lookup(k variant.Key) (VariantStats, bool) {
	s, ok := c.Variants[k]
	return s, ok
}

// nextDataLine advances past blank and comment lines, tracking physical lines.
func (pr *reader) nextDataLine() (string, error) {
	for pr.sc.Scan() {
		pr.line++
		line := pr.sc.Text()
		if parseutil.IsBlank(line) {
			continue
		}
		if pr.cfg.Comment != 0 && line[0] == pr.cfg.Comment {
			continue
		}
		return line, nil
	}
	if err := pr.sc.Err(); err != nil {
		return "", stgerr.Wrap(stgerr.KindMalformedInput, pr.cfg.CohortName,
			"stream read failure", err)
	}
	return "", io.EOF
}

// resolveHeader maps canonical field names onto column indices.
func (pr *reader) resolveHeader(header []string) {
	alias := map[string][]string{
		"chrom":     {"CHR", "CHROM", "CHROMOSOME"},
		"pos":       {"POS", "BP", "POSITION"},
		"ref":       {"REF", "A2", "A0", "NON_EFFECT_ALLELE", "NEA"},
		"effect":    {"EFFECT", "ALT", "A1", "EFFECT_ALLELE", "EA"},
		"beta":      {"BETA", "EFFECT_SIZE", "OR"},
		"se":        {"SE", "STDERR", "SE_BETA"},
		"p":         {"P", "PVAL", "PVALUE"},
		"eaf":       {"EAF", "MAF", "FREQ", "A1_FREQ"},
		"n":         {"N", "NEFF", "SAMPLESIZE"},
		"ncases":    {"N_CASE", "NCASES", "NCAS"},
		"ncontrols": {"N_CONTROL", "NCONTROLS", "NCON"},
	}
	override := map[string]string{
		"chrom": pr.cfg.Columns.Chrom, "pos": pr.cfg.Columns.Pos,
		"ref": pr.cfg.Columns.Ref, "effect": pr.cfg.Columns.Effect,
		"beta": pr.cfg.Columns.Beta, "se": pr.cfg.Columns.SE,
		"p": pr.cfg.Columns.PValue, "eaf": pr.cfg.Columns.EAF,
		"n": pr.cfg.Columns.N, "ncases": pr.cfg.Columns.NCases,
		"ncontrols": pr.cfg.Columns.NControls,
	}
	for field, names := range alias {
		if custom := override[field]; custom != "" {
			if i := parseutil.IndexOf(header, custom); i >= 0 {
				pr.cols[field] = i
			}
			continue
		}
		if i := parseutil.IndexOf(header, names...); i >= 0 {
			pr.cols[field] = i
		}
	}
}

// positionalColumns assigns the fixed headerless order:
// chrom pos ref effect beta se [p] [eaf] [n].
func (pr *reader) positionalColumns() {
	for i, f := range []string{"chrom", "pos", "ref", "effect", "beta", "se", "p", "eaf", "n"} {
		pr.cols[f] = i
	}
}

// rowKeyHint extracts a best-effort key string for rejection logs.
func (pr *reader) rowKeyHint(row []string) string {
	get := func(f string) string {
		if i, ok := pr.cols[f]; ok && i < len(row) {
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

// parseRow validates one data row and converts it to VariantStats.
func (pr *reader) parseRow(row []string) (VariantStats, error) {
	var st VariantStats
	get := func(f string) (string, error) {
		i, ok := pr.cols[f]
		if !ok {
			return "", stgerr.New(stgerr.KindMissingField, pr.cfg.CohortName,
				"no column bound to field "+f)
		}
		return parseutil.Cell(row, i, f)
	}
	opt := func(f string) (string, bool, error) {
		i, ok := pr.cols[f]
		if !ok || i >= len(row) {
			return "", false, nil
		}
		v := strings.TrimSpace(row[i])
		if v == "" || strings.EqualFold(v, "NA") {
			return "", false, nil
		}
		return v, true, nil
	}

	cs, err := get("chrom")
	if err != nil {
		return st, err
	}
	chrom, err := variant.ParseChromosome(cs)
	if err != nil {
		return st, err
	}
	ps, err := get("pos")
	if err != nil {
		return st, err
	}
	pos, err := parseutil.Int(ps, "pos")
	if err != nil {
		return st, err
	}
	rs, err := get("ref")
	if err != nil {
		return st, err
	}
	es, err := get("effect")
	if err != nil {
		return st, err
	}
	v, err := variant.New(chrom, pos, rs, es)
	if err != nil {
		return st, err
	}
	st.Var = v

	bs, err := get("beta")
	if err != nil {
		return st, err
	}
	beta, err := parseutil.Float(bs, "beta")
	if err != nil {
		return st, err
	}
	st.Beta = beta

	ss, err := get("se")
	if err != nil {
		return st, err
	}
	seV, err := parseutil.Float(ss, "se")
	if err != nil {
		return st, err
	}
	if seV <= 0 {
		return st, stgerr.New(stgerr.KindInvalidValue, pr.cfg.CohortName,
			fmt.Sprintf("SE must be > 0, got %v", seV))
	}
	st.SE = seV

	if pv, present, err := opt("p"); err != nil {
		return st, err
	} else if present {
		p, err := parseutil.Float(pv, "p")
		if err != nil {
			return st, err
		}
		if p < 0 || p > 1 {
			return st, stgerr.New(stgerr.KindInvalidValue, pr.cfg.CohortName,
				fmt.Sprintf("p-value must be in [0,1], got %v", p))
		}
		st.PValue, st.HasP = p, true
	}

	if ev, present, err := opt("eaf"); err != nil {
		return st, err
	} else if present {
		eaf, err := parseutil.Float(ev, "eaf")
		if err != nil {
			return st, err
		}
		if eaf < 0 || eaf > 1 {
			return st, stgerr.New(stgerr.KindInvalidValue, pr.cfg.CohortName,
				fmt.Sprintf("EAF must be in [0,1], got %v", eaf))
		}
		st.EAF, st.HasEAF = eaf, true
	}

	if nv, present, err := opt("n"); err != nil {
		return st, err
	} else if present {
		n, err := parseutil.Float(nv, "n")
		if err != nil {
			return st, err
		}
		if n <= 0 {
			return st, stgerr.New(stgerr.KindInvalidValue, pr.cfg.CohortName,
				fmt.Sprintf("N must be > 0, got %v", n))
		}
		st.N = n
	}
	for f, dst := range map[string]*float64{"ncases": &st.NCases, "ncontrols": &st.NControls} {
		if vs, present, err := opt(f); err != nil {
			return st, err
		} else if present {
			val, err := parseutil.Float(vs, f)
			if err != nil {
				return st, err
			}
			if val < 0 {
				return st, stgerr.New(stgerr.KindInvalidValue, pr.cfg.CohortName,
					fmt.Sprintf("%s must be >= 0, got %v", f, val))
			}
			*dst = val
		}
	}
	return st, nil
}
