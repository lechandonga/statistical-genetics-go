// Package variant defines the canonical variant model shared by every
// input source (summary statistics, LD scores, genotype matrices) and the
// allele/strand alignment rules.
//
// A variant is identified by (chromosome, position, reference allele,
// effect allele). All alleles are upper-cased. The canonical Key uses
// alphabetically ordered alleles so that allele swaps hash identically.
package variant

import (
	"fmt"
	"strings"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
)

// Chromosome is a validated chromosome label ("1".."22", "X", "Y", "MT").
type Chromosome string

// Variant is the normalized variant identifier.
type Variant struct {
	Chrom  Chromosome
	Pos    int
	Ref    string // reference allele, upper-case A/C/G/T sequence
	Effect string // effect/alternative allele, upper-case A/C/G/T sequence
}

// Key is the allele-order-independent canonical key.
type Key string

// validChromosomes is the closed set of accepted labels after normalization.
var validChromosomes = map[string]bool{
	"1": true, "2": true, "3": true, "4": true, "5": true,
	"6": true, "7": true, "8": true, "9": true, "10": true,
	"11": true, "12": true, "13": true, "14": true, "15": true,
	"16": true, "17": true, "18": true, "19": true, "20": true,
	"21": true, "22": true,
	"X": true, "Y": true, "MT": true,
}

// ParseChromosome validates and normalizes a chromosome label.
// Accepted inputs: "1".."22", "X", "Y", "M", "MT" with optional "chr"
// prefix, any capitalization. "M" is normalized to "MT".
func ParseChromosome(s string) (Chromosome, error) {
	t := strings.TrimSpace(s)
	t = strings.ToUpper(t)
	t = strings.TrimPrefix(t, "CHR")
	if t == "" {
		return "", stgerr.New(stgerr.KindMissingField, "variant", "empty chromosome")
	}
	if t == "M" {
		t = "MT"
	}
	if !validChromosomes[t] {
		return "", stgerr.New(stgerr.KindInvalidValue, "variant",
			fmt.Sprintf("unrecognized chromosome %q", s))
	}
	return Chromosome(t), nil
}

// complement maps each IUPAC DNA base to its Watson-Crick complement.
var complement = [256]byte{}

func init() {
	for i := range complement {
		complement[i] = 0
	}
	complement['A'] = 'T'
	complement['T'] = 'A'
	complement['C'] = 'G'
	complement['G'] = 'C'
}

// complementBase returns the Watson-Crick complement of a single base.
func complementBase(b byte) byte { return complement[b] }

// NormalizeAlleles upper-cases and validates two alleles. Both must be
// non-empty sequences over A/C/G/T and must differ from each other.
func NormalizeAlleles(a1, a2 string) (string, string, error) {
	r := strings.ToUpper(strings.TrimSpace(a1))
	e := strings.ToUpper(strings.TrimSpace(a2))
	if r == "" || e == "" {
		return "", "", stgerr.New(stgerr.KindMissingField, "variant",
			"empty allele")
	}
	if r == e {
		return "", "", stgerr.New(stgerr.KindInvalidValue, "variant",
			fmt.Sprintf("reference and effect alleles identical: %s", r))
	}
	for i := 0; i < len(r); i++ {
		if complement[r[i]] == 0 {
			return "", "", stgerr.New(stgerr.KindUnknownAllele, "variant",
				fmt.Sprintf("illegal base %q in allele %q", r[i], r))
		}
	}
	for i := 0; i < len(e); i++ {
		if complement[e[i]] == 0 {
			return "", "", stgerr.New(stgerr.KindUnknownAllele, "variant",
				fmt.Sprintf("illegal base %q in allele %q", e[i], e))
		}
	}
	return r, e, nil
}

// New constructs and validates a variant.
func New(chrom Chromosome, pos int, ref, effect string) (Variant, error) {
	if string(chrom) == "" {
		return Variant{}, stgerr.New(stgerr.KindMissingField, "variant", "empty chromosome")
	}
	if !validChromosomes[string(chrom)] {
		return Variant{}, stgerr.New(stgerr.KindInvalidValue, "variant",
			fmt.Sprintf("unrecognized chromosome %q", chrom))
	}
	if pos <= 0 {
		return Variant{}, stgerr.New(stgerr.KindInvalidValue, "variant",
			fmt.Sprintf("position must be positive, got %d", pos))
	}
	r, e, err := NormalizeAlleles(ref, effect)
	if err != nil {
		return Variant{}, err
	}
	return Variant{Chrom: chrom, Pos: pos, Ref: r, Effect: e}, nil
}

// orderedAlleles returns the two alleles in deterministic (alphabetical) order.
func orderedAlleles(a1, a2 string) (string, string) {
	if a1 <= a2 {
		return a1, a2
	}
	return a2, a1
}

// CanonicalKey returns the order-independent key for a variant.
func CanonicalKey(chrom Chromosome, pos int, a1, a2 string) Key {
	x, y := orderedAlleles(a1, a2)
	return Key(fmt.Sprintf("%s:%d:%s:%s", chrom, pos, x, y))
}

// Key returns the canonical key of v.
func (v Variant) Key() Key {
	return CanonicalKey(v.Chrom, v.Pos, v.Ref, v.Effect)
}

// IsSNP reports whether both alleles are single bases.
func (v Variant) IsSNP() bool { return len(v.Ref) == 1 && len(v.Effect) == 1 }

// IsAmbiguous reports whether v is an A/T or C/G SNP. Such loci are
// strand-ambiguous because complement maps them onto an allele swap.
func (v Variant) IsAmbiguous() bool {
	if !v.IsSNP() {
		return false
	}
	r, e := v.Ref[0], v.Effect[0]
	return (r == 'A' && e == 'T') || (r == 'T' && e == 'A') ||
		(r == 'C' && e == 'G') || (r == 'G' && e == 'C')
}

// IsSwapOf reports whether other equals v with ref/effect exchanged.
func (v Variant) IsSwapOf(other Variant) bool {
	return v.Chrom == other.Chrom && v.Pos == other.Pos &&
		v.Ref == other.Effect && v.Effect == other.Ref
}

// IsComplementOf reports whether other equals v with both alleles
// complemented without allele exchange: e.g. A/G vs T/C.
func (v Variant) IsComplementOf(other Variant) bool {
	cv, err := v.Complemented()
	if err != nil {
		return false
	}
	return cv.Chrom == other.Chrom && cv.Pos == other.Pos &&
		cv.Ref == other.Ref && cv.Effect == other.Effect
}

// Swapped returns v with ref and effect exchanged.
func (v Variant) Swapped() Variant {
	return Variant{Chrom: v.Chrom, Pos: v.Pos, Ref: v.Effect, Effect: v.Ref}
}

// Complemented returns v with each allele base-complemented.
// Indels return KindIncompatibleAlleles because a strand flip is only
// defined allele-wise for single bases and must not be guessed for indels.
func (v Variant) Complemented() (Variant, error) {
	if !v.IsSNP() {
		return Variant{}, stgerr.New(stgerr.KindIncompatibleAlleles, "variant",
			fmt.Sprintf("strand complement undefined for non-SNP %s", v.Key()))
	}
	r := complement[v.Ref[0]]
	e := complement[v.Effect[0]]
	if r == 0 || e == 0 {
		return Variant{}, stgerr.New(stgerr.KindUnknownAllele, "variant",
			fmt.Sprintf("cannot complement %s", v.Key()))
	}
	return Variant{Chrom: v.Chrom, Pos: v.Pos, Ref: string(r), Effect: string(e)}, nil
}
