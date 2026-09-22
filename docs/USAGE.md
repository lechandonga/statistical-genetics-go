# Usage, compatibility and verification

## Building

Requires Go 1.25+ (developed with Go 1.26). The only external dependency is
pinned in `go.mod`:

```
gonum.org/v1/gonum v0.17.0
```

All other code is standard library. Build and test:

```sh
go build ./...
go test ./...
```

## Input formats

All inputs are delimited text files (default tab; `Delim: "whitespace"` is
also supported by the library parsers). Header names are matched
case-insensitively with aliases; use the per-package `ColumnMap` for
non-default headers.

### Summary statistics (`.tsv`)

Required columns: `CHR, POS, REF, ALT/BETA` columns: `BETA/OR/EFFECT_SIZE`,
`SE/STDERR`; optional: `P/PVAL`, `EAF/MAF/FREQ`, `N/NEFF`, `NCASES`,
`NCONTROLS`. Ref is treated as the non-effect allele and ALT as the effect
allele. `N` may be per-row or supplied once as `ParseConfig.DefaultN`.

Validation rules (each yields a distinguishable error kind):

- SE must be strictly > 0; p in [0,1]; EAF in [0,1]; N > 0.
- NaN/Inf numeric fields are rejected explicitly.
- Duplicate variant keys within one file are rejected.

### LD scores (`.tsv`)

Columns `CHR POS REF ALT L2` plus optional `MAF`, `WEIGHT`. L2 must be
non-negative (floor configurable via `MinL2`).

### Genotype matrix

- First line: optional header naming the samples.
- First column (with `HasSampleCol`): variant id `chr:pos:ref:effect`.
- Dosages in `{0,1,2}`; missing tokens `NA`, `.`, `-9` become the internal
  sentinel `-1` and a per-cell boolean mask (never a numeric 0 guess).

Rejected rows are collected in `variant.Report`, never silently skipped.
In strict mode parsing aborts at the first rejection; with
`SkipReject: true` parsing continues and `Report.Summary()` lists counts by
stable kind, e.g.:

```
source=A accepted=10 rejected=0 rejects{}
source=ld accepted=7 rejected=3 rejects{strand_ambiguous=3}
```

## CLI

```sh
go run ./cmd/statgen \
  -sumstats testdata/cohort_A.tsv,testdata/cohort_B.tsv \
  -names A,B \
  -ldscores testdata/ldscores.tsv \
  -genotypes testdata/genotypes.tsv \
  -meta-model fixed_iv|fixed_z|random_dl \
  -overlap-policy reject|flag|deduplicate \
  -ambiguous reject|mark|force \
  -negative-h2 report|truncate_zero \
  -score-scale raw|mean_center|standardized \
  -workers 4 \
  -max-variants 1000000 -max-samples 100000 -timeout 5m \
  -out result.json
```

Input summaries and judgments are printed to stderr; `-out` emits JSON.
Exit code 2 = bad invocation, 1 = run rejected (the stable error kind is
printed, e.g. `kind=strand_ambiguous` or `kind=limit_exceeded`).

## Library entry points

| Concern | API |
|---|---|
| Error taxonomy | `stgerr.Kind`, `stgerr.Is(err, kind)` |
| Variant/alignment | `variant.New`, `variant.ParseChromosome`, `Aligner.Align` |
| Summary stats | `sumstats.Parse`, `sumstats.ParseFile` |
| LD scores | `ldscore.Parse`, `ldscore.Match` |
| Genotypes | `genotype.Parse`, `genotype.ImputedRow` |
| Heritability | `ldsc.Estimate` (result `Status`, CI, weights, report) |
| Meta-analysis | `meta.Prepare`, `Analysis.Pool`, `Analysis.PoolAll` |
| PRS | `prs.NewScoreSet`, `prs.FromCohort`, `prs.FromPooled`, `ScoreSet.Score`, `prs.Evaluate` |
| Numerical seam | `numlib.WeightedOLS`, `numlib.Safe`, CDF/quantile wrappers |
| Orchestration | `pipeline.Run(pipeline.Config)` |

## Verification: what is tested

Tests concentrate on boundary and abnormal data. Every test logs an input
summary and the judgment basis.

- **Strand handling**: all four orientations, ambiguous A/T & C/G under
  reject/mark/force, indel refusal, canonical-key swap invariance.
- **Parsing rejections**: SE ≤ 0, NaN beta, p outside [0,1], N ≤ 0, bad
  chromosome, unknown allele, identical alleles, duplicates, ragged rows,
  malformed variant ids — each asserted by stable error kind.
- **Heritability**: negative h2 reported (`status=negative`) vs clamped,
  too-few-variants → out-of-bounds, positive SE/CI, repeat-run bit
  identity, serial == parallel bit identity.
- **Meta-analysis**: swap sign realignment (B cohort −0.09 → +0.09),
  fixed/random/Z models, Q/I2/tau2, duplicate cohort, one-sample overlap,
  reject/flag/deduplicate policies, unsupported model rejection.
- **PRS**: mean/zero/reject missing-genotype policies, missing-locus
  error/skip, ambiguous-locus refusal, set-change sign-convention
  invariance, standardized mean 0 / SD 1, perfect-separation AUC = 1.
- **Numerical seam**: exact known-line weighted OLS recovery, singular
  system rejection, unsupported method rejection, Kahan summation.
- **Pipeline**: full run reproducibility, serial/parallel bit equality for
  LD regression/meta/PRS, variant/sample/cohort budget rejection,
  duplicate-cohort policy, concurrent budget counter correctness.

```sh
go test ./... -v   # inspect per-case INPUT/JUDGE log lines
go vet ./...
```

## Determinism contract

Given identical file bytes and configuration:

1. rejected-row reports, h2, CIs, meta effects and PRS values are identical
   across runs;
2. `Workers=1` and `Workers=N` produce bit-identical floating-point results
   (asserted in `ldsc`, `meta`, `prs` and `pipeline` tests);
3. there is no random seed anywhere in the analysis path — the pipeline
   performs no sampling.
