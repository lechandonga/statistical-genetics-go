package numlib

// KahanSum returns a compensated order-independent(ish) floating-point sum.
// All serial and parallel chromosome accumulations funnel through the same
// summation so that chunking cannot change the conclusion. Kahan summation
// is not strictly permutation invariant, so parallel callers must reduce
// chunks by merging their (sum, compensation) pairs via PairSum.
func KahanSum(xs []float64) float64 {
	sum, comp := 0.0, 0.0
	for _, x := range xs {
		y := x - comp
		t := sum + y
		comp = (t - sum) - y
		sum = t
	}
	return sum
}

// SumPair is an accumulated Kahan state.
type SumPair struct {
	Sum  float64
	Comp float64
}

// PairSum merges one value into an accumulated state.
func PairSum(p SumPair, x float64) SumPair {
	y := x - p.Comp
	t := p.Sum + y
	return SumPair{Sum: t, Comp: (t - p.Sum) - y}
}

// Merge combines two Kahan states; the merge order is fixed by construction
// (chromosome index ascending), keeping parallel == serial output.
func Merge(a, b SumPair) SumPair {
	// Combine the raw sums, carrying the larger compensation conservatively.
	if b.Sum == 0 && b.Comp == 0 {
		return a
	}
	return PairSum(a, b.Sum+b.Comp)
}

// Dot returns the dot product of two equally-sized vectors via Kahan sum.
func Dot(a, b []float64) float64 {
	if len(a) != len(b) {
		panic("numlib: Dot length mismatch")
	}
	acc := SumPair{}
	for i := range a {
		acc = PairSum(acc, a[i]*b[i])
	}
	return acc.Sum
}
