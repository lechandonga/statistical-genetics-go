package prs

import (
	"math"
	"sort"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/numlib"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
)

// Phenotype holds case/control status (1 = case, 0 = control).
// Samples with unknown status are excluded from evaluation.
type Phenotype struct {
	SampleIDs []string
	Status    []float64
}

// EvalResult evaluates score discrimination.
type EvalResult struct {
	N         int
	NCases    int
	NControls int
	AUC       float64
	R2        float64 // linear R2 of status explained by score
	MeanCase  float64
	MeanCtrl  float64
}

// Evaluate scores against phenotypes, restricting to samples known in both.
func Evaluate(s *Scores, p Phenotype) (EvalResult, error) {
	if s == nil {
		return EvalResult{}, stgerr.New(stgerr.KindMissingField, "prs", "nil scores")
	}
	if len(p.SampleIDs) != len(p.Status) {
		return EvalResult{}, stgerr.New(stgerr.KindMalformedInput, "prs",
			"phenotype ids and status length differ")
	}
	idx := make(map[string]int, len(s.SampleIDs))
	for i, id := range s.SampleIDs {
		idx[id] = i
	}
	type pair struct {
		score  float64
		status float64
	}
	var pairs []pair
	for j, id := range p.SampleIDs {
		y := p.Status[j]
		if math.IsNaN(y) || (y != 0 && y != 1) {
			continue // unknown/illegal status: excluded, not assumed control
		}
		if i, ok := idx[id]; ok {
			pairs = append(pairs, pair{s.Values[i], y})
		}
	}
	if len(pairs) == 0 {
		return EvalResult{}, stgerr.New(stgerr.KindVariantNotFound, "prs",
			"no overlap between scored samples and phenotype samples")
	}

	var cases, controls int
	var sumCase, sumCtrl float64
	for _, pr := range pairs {
		if pr.status == 1 {
			cases++
			sumCase += pr.score
		} else {
			controls++
			sumCtrl += pr.score
		}
	}
	if cases == 0 || controls == 0 {
		return EvalResult{}, stgerr.New(stgerr.KindInvalidValue, "prs",
			"AUC/R2 require both cases and controls")
	}

	// AUC by rank with average handling of ties (Mann-Whitney U).
	ranked := make([]pair, len(pairs))
	copy(ranked, pairs)
	sort.SliceStable(ranked, func(i, j int) bool { return ranked[i].score < ranked[j].score })
	ranks := make([]float64, len(ranked))
	for i := 0; i < len(ranked); {
		j := i + 1
		for j < len(ranked) && ranked[j].score == ranked[i].score {
			j++
		}
		avg := (float64(i+1) + float64(j)) / 2
		for k := i; k < j; k++ {
			ranks[k] = avg
		}
		i = j
	}
	var rankSumCases float64
	for k, pr := range ranked {
		if pr.status == 1 {
			rankSumCases += ranks[k]
		}
	}
	u := rankSumCases - float64(cases)*(float64(cases)+1)/2
	auc := u / float64(cases*controls)

	// Linear R2 against the binary status on the scored sample set.
	mean := (sumCase + sumCtrl) / float64(len(pairs))
	tss := 0.0
	rss := 0.0
	accT, accR := numlib.SumPair{}, numlib.SumPair{}
	for _, pr := range pairs {
		accT = numlib.PairSum(accT, (pr.status-mean)*(pr.status-mean))
	}
	// Simple regression intercept+slope for R2.
	var sx, sy, sxx, sxy, sn float64
	for _, pr := range pairs {
		sx += pr.score
		sy += pr.status
		sxx += pr.score * pr.score
		sxy += pr.score * pr.status
		sn++
	}
	b1 := (sn*sxy - sx*sy) / (sn*sxx - sx*sx)
	b0 := (sy - b1*sx) / sn
	for _, pr := range pairs {
		pred := b0 + b1*pr.score
		accR = numlib.PairSum(accR, (pr.status-pred)*(pr.status-pred))
	}
	tss, rss = accT.Sum, accR.Sum
	r2 := 0.0
	if tss > 0 {
		r2 = 1 - rss/tss
	}

	return EvalResult{
		N: len(pairs), NCases: cases, NControls: controls,
		AUC: auc, R2: r2,
		MeanCase: sumCase / float64(cases),
		MeanCtrl: sumCtrl / float64(controls),
	}, nil
}
