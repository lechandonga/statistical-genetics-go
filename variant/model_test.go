package variant

import (
	"math"
	"testing"
)

func TestNormalizeAllele(t *testing.T) {
	cases := []struct {
		in      string
		want    Allele
		wantErr Reason
	}{
		{" a ", "A", ""},
		{"acgt", "ACGT", ""},
		{"", "", ReasonMissingAlleles},
		{"A I", "", ReasonInvalidAllele},
		{"D", "", ReasonInvalidAllele},
		{"A.B", "", ReasonInvalidAllele},
	}
	for _, c := range cases {
		got, err := NormalizeAllele(c.in)
		if c.wantErr != "" {
			r, ok := AsReject(err)
			if !ok || r.Reason != c.wantErr {
				t.Fatalf("input=%q want reason=%s got %v", c.in, c.wantErr, err)
			}
			t.Logf("[accept-reject] input=%q -> reason=%s detail=%s", c.in, r.Reason, r.Detail)
			continue
		}
		if err != nil || got != c.want {
			t.Fatalf("input=%q want %s got %s err=%v", c.in, c.want, got, err)
		}
		t.Logf("[accept] allele input=%q normalized=%s", c.in, got)
	}
}

func TestNewVariant(t *testing.T) {
	cases := []struct {
		name    string
		id      string
		chrom   string
		pos     int64
		ea, nea string
		eaf     float64
		strand  Strand
		wantR   Reason
		check   func(*testing.T, *Variant)
	}{
		{"id only", "rs1", "", 0, "A", "G", math.NaN(), StrandUnknown, "",
			func(t *testing.T, v *Variant) {
				if v.Key() != MatchKey("RS1") {
					t.Fatalf("key=%s", v.Key())
				}
			}},
		{"chr prefix + m", "rsM", "chrm", 10, "A", "G", 0.1, StrandUnknown, "",
			func(t *testing.T, v *Variant) {
				if v.Chrom != "MT" {
					t.Fatalf("chrom=%s", v.Chrom)
				}
			}},
		{"numeric chr 23", "rsX", "23", 1, "C", "T", 0.5, StrandUnknown, "",
			func(t *testing.T, v *Variant) {
				if v.Chrom != "X" {
					t.Fatalf("chrom=%s", v.Chrom)
				}
			}},
		{"pos without chrom", "rs2", "", 5, "A", "G", math.NaN(), StrandUnknown, ReasonInvalidChrom, nil},
		{"chrom without pos", "rs3", "1", 0, "A", "G", math.NaN(), StrandUnknown, ReasonInvalidPosition, nil},
		{"bad chrom", "rs4", "Z", 9, "A", "G", math.NaN(), StrandUnknown, ReasonInvalidChrom, nil},
		{"identical alleles", "rs5", "1", 9, "A", "A", math.NaN(), StrandUnknown, ReasonInvalidAllele, nil},
		{"eaf oob", "rs6", "1", 9, "A", "G", 1.5, StrandUnknown, ReasonInvalidFrequency, nil},
		{"eaf neg", "rs7", "1", 9, "A", "G", -0.001, StrandUnknown, ReasonInvalidFrequency, nil},
		{"no identity", "", "", 0, "A", "G", math.NaN(), StrandUnknown, ReasonMissingID, nil},
		{"reverse strand flipped", "rs8", "1", 9, "T", "C", math.NaN(), StrandReverse, "",
			func(t *testing.T, v *Variant) {
				if v.EA != "A" || v.NEA != "G" || v.Strand != StrandForward {
					t.Fatalf("reverse strand not normalized: %+v", v)
				}
			}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			v, err := New(c.id, c.chrom, c.pos, c.ea, c.nea, c.eaf, c.strand)
			t.Logf("[input] id=%s chrom=%s pos=%d ea=%s nea=%s strand=%d",
				c.id, c.chrom, c.pos, c.ea, c.nea, c.strand)
			if c.wantR != "" {
				r, ok := AsReject(err)
				if !ok || r.Reason != c.wantR {
					t.Fatalf("want reason=%s got err=%v", c.wantR, err)
				}
				t.Logf("[judgement] REJECT reason=%s detail=%s", r.Reason, r.Detail)
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			t.Logf("[judgement] ACCEPT key=%s chrom=%s pos=%d ea=%s nea=%s",
				v.Key(), v.Chrom, v.Pos, v.EA, v.NEA)
			if c.check != nil {
				c.check(t, v)
			}
		})
	}
}

func TestAmbiguous(t *testing.T) {
	mk := func(ea, nea string) *Variant {
		v, err := New("rs", "1", 1, ea, nea, math.NaN(), StrandUnknown)
		if err != nil {
			t.Fatal(err)
		}
		return v
	}
	if !mk("A", "T").IsAmbiguous() {
		t.Fatal("A/T must be ambiguous")
	}
	if !mk("g", "c").IsAmbiguous() {
		t.Fatal("G/C must be ambiguous")
	}
	if mk("A", "G").IsAmbiguous() {
		t.Fatal("A/G is not ambiguous")
	}
	if mk("AC", "GT").IsAmbiguous() {
		t.Fatal("indel is not SNV-ambiguous")
	}
	t.Logf("[input-summary] A/T, C/G marked ambiguous; A/G, indels unambiguous")
}
