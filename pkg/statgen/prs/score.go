package prs

import (
	"fmt"
	"math"
	"sync"

	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/genotype"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/numlib"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/stgerr"
	"github.com/lechandonga/statistical-genetics-go/pkg/statgen/variant"
)

// Score computes every sample's score deterministically. The variant loop
// is the fixed outer order (sorted by key in ScoreSet); each sample's value
// is accumulated with Kahan compensation. Parallel mode shards over samples
// with per-shard output slots, so it is bit-identical to serial.
func (set *ScoreSet) Score(m *genotype.Matrix, cfg Config) (*Scores, error) {
	if m == nil {
		return nil, stgerr.New(stgerr.KindMissingField, "prs", "nil genotype matrix")
	}
	scale := cfg.ScoreScale
	if scale == "" {
		scale = ScaleRaw
	}
	if scale != ScaleRaw && scale != ScaleMeanCenter && scale != ScaleStandardized {
		return nil, stgerr.New(stgerr.KindUnsupportedOption, "prs",
			"unsupported score scale: "+string(scale))
	}
	locusPolicy := cfg.MissingLocus
	if locusPolicy == "" {
		locusPolicy = MissingLocusError
	}
	if locusPolicy != MissingLocusError && locusPolicy != MissingLocusSkip {
		return nil, stgerr.New(stgerr.KindUnsupportedOption, "prs",
			"unsupported missing-locus policy: "+string(locusPolicy))
	}

	// Resolve effect loci against matrix rows exactly once.
	type resolved struct {
		e      Effect
		row    int
		center float64 // weighted mean dosage contribution
	}
	rs := make([]resolved, 0, len(set.Effects))
	report := ScoreReport{
		SignConvention: set.signConvention,
		MissingData:    map[variant.Key]int{},
	}
	for _, e := range set.Effects {
		idx, ok := m.VariantIndex(e.Target.Key())
		if !ok {
			if locusPolicy == MissingLocusError {
				return nil, stgerr.New(stgerr.KindVariantNotFound, "prs",
					"effect locus absent from genotype matrix: "+string(e.Target.Key()))
			}
			report.SkippedLoci = append(report.SkippedLoci, e.Target.Key())
			continue
		}
		rs = append(rs, resolved{e: e, row: idx})
	}
	if len(rs) == 0 {
		return nil, stgerr.New(stgerr.KindVariantNotFound, "prs",
			"no scoring loci present in genotype matrix")
	}

	// Prepare imputed dosage rows once (deterministic mean/zero policy),
	// recording per-locus imputation counts and weighted centering constants.
	type preparedRow struct {
		weight float64
		dosage []float64
		center float64
	}
	rows := make([]preparedRow, len(rs))
	for i, r := range rs {
		d, miss, _ := m.Row(r.e.Target.Key())
		imp, imputed, err := genotype.ImputedRow(d, miss, cfg.MissingGenotype)
		if err != nil {
			return nil, err
		}
		if imputed > 0 {
			report.MissingData[r.e.Target.Key()] = imputed
		}
		mean := numlib.KahanSum(imp) / float64(len(imp))
		rows[i] = preparedRow{
			weight: r.e.Weight,
			dosage: imp,
			center: r.e.Weight * mean,
		}
	}
	report.UsedLoci = len(rows)

	nSamp := m.NumSamples()
	values := make([]float64, nSamp)

	compute := func(lo, hi int) {
		for s := lo; s < hi; s++ {
			acc := numlib.SumPair{}
			center := numlib.SumPair{}
			for _, pr := range rows {
				acc = numlib.PairSum(acc, pr.weight*pr.dosage[s])
				center = numlib.PairSum(center, pr.center)
			}
			switch scale {
			case ScaleMeanCenter, ScaleStandardized:
				values[s] = acc.Sum - center.Sum
			default:
				values[s] = acc.Sum
			}
		}
	}

	workers := cfg.Workers
	if workers <= 1 || nSamp < 64 {
		compute(0, nSamp)
	} else {
		if workers > nSamp {
			workers = nSamp
		}
		var wg sync.WaitGroup
		for w := 0; w < workers; w++ {
			lo := w * nSamp / workers
			hi := (w + 1) * nSamp / workers
			wg.Add(1)
			go func() { defer wg.Done(); compute(lo, hi) }()
		}
		wg.Wait()
	}

	if scale == ScaleStandardized {
		mean := numlib.KahanSum(values) / float64(nSamp)
		var sq float64
		acc := numlib.SumPair{}
		for _, v := range values {
			d := v - mean
			acc = numlib.PairSum(acc, d*d)
		}
		sq = acc.Sum / float64(nSamp)
		sd := math.Sqrt(math.Max(sq, 0))
		if sd == 0 {
			return nil, stgerr.New(stgerr.KindOutOfBounds, "prs",
				"zero score variance: cannot standardize")
		}
		for i := range values {
			values[i] = (values[i] - mean) / sd
		}
	}

	return &Scores{
		SampleIDs: append([]string(nil), m.Samples...),
		Values:    values,
		UsedLoci:  len(rows),
		Report:    report,
	}, nil
}

// String gives a deterministic score-set summary for test logs.
func (rep ScoreReport) String() string {
	return fmt.Sprintf("used_loci=%d skipped_loci=%d imputed_loci=%d convention=%q",
		rep.UsedLoci, len(rep.SkippedLoci), len(rep.MissingData), rep.SignConvention)
}
