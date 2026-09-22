package variant

import (
	"fmt"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
)

// Orientation describes how a source variant's alleles map to the target.
type Orientation int

const (
	// OrientMatch: alleles identical, same order.
	OrientMatch Orientation = iota
	// OrientSwap: alleles exchanged (effect sign must flip).
	OrientSwap
	// OrientFlip: both alleles strand-complemented, same order.
	OrientFlip
	// OrientSwapFlip: complement plus exchange (effect sign must flip).
	OrientSwapFlip
)

// Sign returns +1 for orientations preserving effect direction and -1 for
// orientations that reverse it (swap, swap+flip). A pure strand flip keeps
// the 5'->3' allele order and therefore does NOT reverse the effect sign:
// effects are attached to the named allele, which is tracked explicitly.
func (o Orientation) Sign() float64 {
	switch o {
	case OrientMatch, OrientFlip:
		return 1
	case OrientSwap, OrientSwapFlip:
		return -1
	default:
		return 0
	}
}

func (o Orientation) String() string {
	switch o {
	case OrientMatch:
		return "match"
	case OrientSwap:
		return "swap"
	case OrientFlip:
		return "flip"
	case OrientSwapFlip:
		return "swap+flip"
	default:
		return "unknown"
	}
}

// AmbiguousPolicy controls strand-ambiguous A/T and C/G SNPs.
type AmbiguousPolicy int

const (
	// AmbiguousReject refuses ambiguous loci (default, safest).
	AmbiguousReject AmbiguousPolicy = iota
	// AmbiguousMark keeps them but flags Aligned.Ambiguous=true.
	AmbiguousMark
	// AmbiguousFlipForce complements regardless; explicit, auditable.
	AmbiguousFlipForce
)

// IndelPolicy controls variants longer than one base.
type IndelPolicy int

const (
	// IndelSwapOnly aligns indels only by exact match or swap.
	IndelSwapOnly IndelPolicy = iota
	// IndelReject refuses all non-SNP variants.
	IndelReject
)

// Aligner holds allele-alignment configuration.
type Aligner struct {
	Ambiguous AmbiguousPolicy
	Indels    IndelPolicy
}

// DefaultAligner rejects ambiguous loci and aligns indels by swap only.
func DefaultAligner() *Aligner {
	return &Aligner{Ambiguous: AmbiguousReject, Indels: IndelSwapOnly}
}

// Aligned bundles the alignment outcome for a variant.
type Aligned struct {
	Target    Variant
	Orient    Orientation
	Sign      float64
	Ambiguous bool // true when the locus is A/T or C/G
	Marked    bool // kept despite ambiguity under AmbiguousMark
}

// Align maps source alleles onto target alleles. It never silently adds
// signed effects: every non-identity orientation is returned explicitly so
// callers multiply beta by Aligned.Sign themselves.
//
// Recognition order matters for ambiguous loci: a source A/T variant that
// was strand-flipped is indistinguishable from the same locus allele-swapped
// (both produce T/A). Such loci are handled by the AmbiguousPolicy and are
// never silently oriented.
func (a *Aligner) Align(target, source Variant) (Aligned, error) {
	if target.Chrom != source.Chrom || target.Pos != source.Pos {
		return Aligned{}, stgerr.New(stgerr.KindIncompatibleAlleles, string(target.Chrom),
			fmt.Sprintf("locus mismatch target=%s source=%s", target.Key(), source.Key()))
	}

	snp := target.IsSNP() && source.IsSNP()
	if !snp && a.Indels == IndelReject {
		return Aligned{}, stgerr.New(stgerr.KindIncompatibleAlleles, string(target.Chrom),
			fmt.Sprintf("non-SNP variant rejected by policy: %s", target.Key()))
	}

	amb := target.IsAmbiguous() || source.IsAmbiguous()

	// 1. Identity: exact same alleles and order.
	if target.Ref == source.Ref && target.Effect == source.Effect {
		if amb && a.Ambiguous == AmbiguousReject {
			return Aligned{}, ambiguousErr(target)
		}
		return Aligned{Target: target, Orient: OrientMatch, Sign: 1, Ambiguous: amb,
			Marked: amb && a.Ambiguous == AmbiguousMark}, nil
	}

	// 2. Pure allele swap (no complement): sign reverses.
	if source.Ref == target.Effect && source.Effect == target.Ref {
		if amb && a.Ambiguous == AmbiguousReject {
			return Aligned{}, ambiguousErr(target)
		}
		return Aligned{Target: target, Orient: OrientSwap, Sign: -1, Ambiguous: amb,
			Marked: amb && a.Ambiguous == AmbiguousMark}, nil
	}

	if !snp {
		// Indels cannot be strand-complemented; match/swap were the only legal maps.
		return Aligned{}, stgerr.New(stgerr.KindIncompatibleAlleles, string(target.Chrom),
			fmt.Sprintf("indel alleles align neither by match nor swap: target=%s source=%s",
				target.Key(), source.Key()))
	}

	// 3. Strand complement, same order: sign preserved, alleles tracked.
	if comp, err := source.Complemented(); err == nil {
		if comp.Ref == target.Ref && comp.Effect == target.Effect {
			// For ambiguous loci "complement, same order" coincides with swap
			// on the original alleles, so direction is unknowable.
			if amb {
				switch a.Ambiguous {
				case AmbiguousReject:
					return Aligned{}, ambiguousErr(target)
				case AmbiguousFlipForce:
					return Aligned{Target: target, Orient: OrientFlip, Sign: 1,
						Ambiguous: true}, nil
				case AmbiguousMark:
					return Aligned{Target: target, Orient: OrientFlip, Sign: 1,
						Ambiguous: true, Marked: true}, nil
				}
			}
			return Aligned{Target: target, Orient: OrientFlip, Sign: 1}, nil
		}
	}

	// 4. Strand complement plus allele swap: sign reverses.
	if c := source.Swapped(); true {
		if comp, err := c.Complemented(); err == nil {
			if comp.Ref == target.Ref && comp.Effect == target.Effect {
				if amb && a.Ambiguous == AmbiguousReject {
					return Aligned{}, ambiguousErr(target)
				}
				return Aligned{Target: target, Orient: OrientSwapFlip, Sign: -1,
					Ambiguous: amb, Marked: amb && a.Ambiguous == AmbiguousMark}, nil
			}
		}
	}

	return Aligned{}, stgerr.New(stgerr.KindIncompatibleAlleles, string(target.Chrom),
		fmt.Sprintf("alleles cannot be aligned: target=%s source=%s",
			target.Key(), source.Key()))
}

func ambiguousErr(v Variant) error {
	return stgerr.New(stgerr.KindStrandAmbiguous, string(v.Chrom),
		fmt.Sprintf("A/T or C/G locus %s is strand-ambiguous and rejected by policy", v.Key()))
}
