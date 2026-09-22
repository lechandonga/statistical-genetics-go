package variant

import (
	"testing"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
)

func mustV(t *testing.T, chrom string, pos int, ref, eff string) Variant {
	t.Helper()
	c, err := ParseChromosome(chrom)
	if err != nil {
		t.Fatalf("chrom: %v", err)
	}
	v, err := New(c, pos, ref, eff)
	if err != nil {
		t.Fatalf("new variant: %v", err)
	}
	return v
}

func TestParseChromosome(t *testing.T) {
	cases := map[string]string{
		"1": "1", "chr22": "22", "CHRx": "X", "y": "Y", "M": "MT", "chrm": "MT",
	}
	for in, want := range cases {
		c, err := ParseChromosome(in)
		if err != nil || string(c) != want {
			t.Fatalf("ParseChromosome(%q)=%q,%v want %q", in, c, err, want)
		}
	}
	for _, bad := range []string{"", "0", "23", "ZZ"} {
		if _, err := ParseChromosome(bad); !stgerr.Is(err, stgerr.KindInvalidValue) &&
			!stgerr.Is(err, stgerr.KindMissingField) {
			t.Fatalf("chrom %q: expected reject, got %v", bad, err)
		}
	}
}

func TestAlleleValidation(t *testing.T) {
	if _, _, err := NormalizeAlleles("a", "g"); err != nil {
		t.Fatalf("lowercase should normalize: %v", err)
	}
	if _, _, err := NormalizeAlleles("A", "N"); !stgerr.Is(err, stgerr.KindUnknownAllele) {
		t.Fatalf("N allele: %v", err)
	}
	if _, _, err := NormalizeAlleles("A", "A"); !stgerr.Is(err, stgerr.KindInvalidValue) {
		t.Fatalf("identical alleles: %v", err)
	}
}

func TestAmbiguous(t *testing.T) {
	at := mustV(t, "1", 1, "A", "T")
	cg := mustV(t, "1", 2, "C", "G")
	ag := mustV(t, "1", 3, "A", "G")
	if !at.IsAmbiguous() || !cg.IsAmbiguous() || ag.IsAmbiguous() {
		t.Fatal("ambiguity classification wrong")
	}
}

func TestAlignOrientations(t *testing.T) {
	a := DefaultAligner()
	tgt := mustV(t, "1", 10, "A", "G")
	cases := []struct {
		ref, eff string
		orient   Orientation
		sign     float64
	}{
		{"A", "G", OrientMatch, 1},     // identity
		{"G", "A", OrientSwap, -1},     // allele swap
		{"T", "C", OrientFlip, 1},      // pure strand complement
		{"C", "T", OrientSwapFlip, -1}, // complement + swap
	}
	for _, c := range cases {
		src := mustV(t, "1", 10, c.ref, c.eff)
		al, err := a.Align(tgt, src)
		if err != nil {
			t.Fatalf("%s/%s: %v", c.ref, c.eff, err)
		}
		if al.Orient != c.orient || al.Sign != c.sign {
			t.Fatalf("%s/%s: orient=%s sign=%v want %s/%v",
				c.ref, c.eff, al.Orient, al.Sign, c.orient, c.sign)
		}
	}
}

func TestAlignAmbiguousRejected(t *testing.T) {
	a := DefaultAligner()
	tgt := mustV(t, "1", 10, "A", "T")
	// T/A is simultaneously swap and strand flip => unknowable.
	src := mustV(t, "1", 10, "T", "A")
	if _, err := a.Align(tgt, src); !stgerr.Is(err, stgerr.KindStrandAmbiguous) {
		t.Fatalf("ambiguous: %v", err)
	}
}

func TestAlignAmbiguousMarkAndForce(t *testing.T) {
	tgt := mustV(t, "1", 10, "A", "T")
	src := mustV(t, "1", 10, "T", "A")
	mark := &Aligner{Ambiguous: AmbiguousMark}
	al, err := mark.Align(tgt, src)
	if err != nil || !al.Ambiguous || !al.Marked {
		t.Fatalf("mark policy: %+v err=%v", al, err)
	}
	force := &Aligner{Ambiguous: AmbiguousFlipForce}
	al2, err := force.Align(tgt, src)
	if err != nil || !al2.Ambiguous {
		t.Fatalf("force policy: %+v err=%v", al2, err)
	}
}

func TestAlignIncompatibleAndIndel(t *testing.T) {
	a := DefaultAligner()
	tgt := mustV(t, "1", 10, "A", "G")
	src := mustV(t, "1", 10, "A", "C") // different allele pair, no map
	if _, err := a.Align(tgt, src); !stgerr.Is(err, stgerr.KindIncompatibleAlleles) {
		t.Fatalf("incompatible: %v", err)
	}
	tgt2 := mustV(t, "1", 11, "AT", "G")
	src2 := mustV(t, "1", 11, "A", "G")
	if _, err := a.Align(tgt2, src2); !stgerr.Is(err, stgerr.KindIncompatibleAlleles) {
		t.Fatalf("indel mismatch: %v", err)
	}
	rejectIndel := &Aligner{Ambiguous: AmbiguousReject, Indels: IndelReject}
	ti := mustV(t, "1", 12, "AT", "G")
	si := mustV(t, "1", 12, "AT", "G")
	if _, err := rejectIndel.Align(ti, si); !stgerr.Is(err, stgerr.KindIncompatibleAlleles) {
		t.Fatalf("indel reject policy: %v", err)
	}
}

func TestCanonicalKeySwapInvariant(t *testing.T) {
	k1 := mustV(t, "1", 1, "A", "G").Key()
	k2 := mustV(t, "1", 1, "G", "A").Key()
	if k1 != k2 {
		t.Fatalf("keys differ across swap: %s vs %s", k1, k2)
	}
}

func TestReportSummaryDeterministic(t *testing.T) {
	r := Report{Source: "x", Accepted: 3}
	r.Add(Rejection{Kind: stgerr.KindDuplicateVariant, Reason: "d1"})
	r.Add(Rejection{Kind: stgerr.KindDuplicateVariant, Reason: "d2"})
	r.Add(Rejection{Kind: stgerr.KindUnknownAllele, Reason: "u1"})
	s1, s2 := r.Summary(), r.Summary()
	if s1 != s2 {
		t.Fatal("summary not deterministic")
	}
	want := "source=x accepted=3 rejected=3 rejects{duplicate_variant=2,unknown_allele=1}"
	if s1 != want {
		t.Fatalf("summary=%q want %q", s1, want)
	}
}
