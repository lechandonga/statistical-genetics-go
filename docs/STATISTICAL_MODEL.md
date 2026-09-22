# Statistical model and assumptions

This document states exactly what the pipeline computes and which
assumptions must hold for results to be interpretable. It is the reference
that error kinds and CLI flags are designed against.

## 1. Canonical variant model

A variant is `(chromosome, position, reference allele, effect allele)`.

- Chromosomes: `1..22, X, Y, MT`; an optional `chr` prefix and case
  differences are normalized; `M` is normalized to `MT`.
- Alleles are upper-cased sequences over `A/C/G/T`; empty alleles, unknown
  bases (`N`, IUPAC ambiguity codes, symbols) and ref==effect are rejected
  with `unknown_allele` / `invalid_value`, never skipped.
- The canonical key is `chr:pos:min(a1,a2):max(a1,a2)`, so allele swaps hash
  identically. The canonical key identifies a *locus*; allele *orientation*
  is tracked separately (see below).

## 2. Allele direction and strand alignment

Every cross-source effect alignment produces an explicit `Orientation`:

| Orientation | Source vs target            | Effect sign |
|-------------|-----------------------------|-------------|
| match       | identical alleles/order     | +1          |
| swap        | alleles exchanged           | -1          |
| flip        | both bases complemented     | +1          |
| swap+flip   | complement + exchange       | -1          |

- A pure strand flip (`A/G -> T/C`) keeps 5'→3' allele order, so the effect
  attached to the named allele keeps its sign. The aligned beta is
  `beta * Orientation.Sign()`, applied by the caller, never silently.
- `A/T` and `C/G` SNPs are **strand ambiguous**: complement maps them onto
  an allele swap, so orientation cannot be inferred from alleles alone.
  Policy `AmbiguousPolicy`:
  - `reject` (default): refuse with `strand_ambiguous`;
  - `mark`: keep but flag `Aligned.Ambiguous = Marked = true`; PRS refuses
    marked loci so a guessed direction never enters a score;
  - `force`: complement regardless, auditable in the result provenance.
- Indels (>1 bp allele) are aligned only by match/swap (`IndelSwapOnly`),
  or refused (`IndelReject`). Strand complementation is undefined for them.

## 3. LD-score heritability (`ldsc`)

Model (Bulik-Sullivan et al. 2015):

```
E[chi2] = intercept + (N / M) * h2 * l2
```

- `chi2 = (beta/SE)^2` after alignment (sign invariant by construction);
  aligned betas still share one allele convention for downstream use.
- Weights: two-step LDSC-style approximation
  `w = 1 / (2 (intercept + N*h2*l2/M)^2)`; a zero/non-finite weight
  contributes zero information.
- Outlier rule: variants with `chi2 > 80` (configurable) are rejected with
  a distinguishable reason and listed in the report.
- The regression is solved by QR through gonum with a Cholesky-based
  covariance check (`numlib`); singular systems return an
  `out_of_bounds`/library error instead of a numeric garbage estimate.
- Uncertainty: `SE(h2)` comes from the weighted-fit covariance; the CI is
  `h2 ± z(1-alpha/2) SE` with default alpha 0.05.

Interpretation status is never hidden by clipping:

| Status           | Condition                                  |
|------------------|--------------------------------------------|
| `valid`          | h2 in [0,1] and CI contained in [0,1]      |
| `negative`       | h2 < 0 (data consistent with zero)         |
| `above_one`      | h2 > 1 (model/data inconsistency)          |
| `ci_out_of_bounds` | point in [0,1] but CI leaves the range   |
| `unidentified`   | degenerate/too few variants                |

`NegativePolicy`: `report` keeps the negative estimate; `truncate_zero`
sets the point estimate to zero while preserving `status=negative`.

### Reproducibility

- Variants are consumed in canonical-key sorted order.
- Parallel mode shards contiguous key ranges and writes private slots,
  merged back in fixed order; all accumulations use Kahan compensation.
- Serial and parallel outputs are compared bit-for-bit in the test suite.

## 4. Multi-cohort meta-analysis (`meta`)

- `fixed_iv`: inverse-variance fixed effects,
  `b = sum w_i b_i / sum w_i`, `SE = 1/sqrt(sum w_i)`, `w_i = 1/SE_i^2`.
- `fixed_z`: Z-score pooling with weights `sqrt(N_i)`.
- `random_dl`: DerSimonian–Laird random effects,
  `tau2 = max(0, (Q - df)/C)`, weights `1/(SE_i^2 + tau2)`.
- Heterogeneity: Cochran's `Q`, `df = k-1`, Q p-value from the chi-square
  CDF, `I2 = max(0,(Q-df)/Q)`, `H = sqrt(Q/df)`.
- Per-cohort effects are aligned to the first effective cohort's allele
  convention; each contribution records cohort, orientation and sign.

### Overlap is never independent evidence

- Equal `OverlapKey` → `duplicate_cohort`; intersecting declared
  `SampleID` sets → `sample_overlap` with the shared count.
- Policies: `reject` (error kinds `duplicate_cohort` / `cohort_overlap`),
  `flag` (analysis proceeds, each affected variant carries a warning and
  the result includes the overlap table), `deduplicate` (the later cohort
  is dropped and recorded in `Analysis.Dropped`).
- Correlated standard errors under overlap are *not* mathematically
  corrected (no full covariance is available); the flag makes the
  assumption explicit. Do not interpret flagged SEs as independent.

## 5. Individual risk score (`prs`)

`score(s) = sum_i w_i dosage_i(s)`

- One immutable, validated `ScoreSet` per variant collection. The fixed
  sign convention is: **weight = per-copy effect of `Target.Effect`**.
  Changing the variant set builds a new set; removing a locus changes a
  score by exactly that locus's signed contribution (test-enforced).
- Missing genotype policies: mean imputation, zero imputation, or
  `reject_sample`. Counts of imputed calls per locus are reported.
- Missing effect loci: `error` or `skip` with an explicit skipped-key list.
- Scales: `raw`, `mean_center`, `standardized` (zero variance is an error,
  not a NaN score).
- Evaluation: AUC via the tie-averaged Mann–Whitney rank statistic, plus a
  linear R2 of status on score; samples with unknown phenotype status are
  excluded rather than assumed controls.

## 6. Third-party numerical code (`numlib`)

- The only third-party dependency is pinned `gonum v0.17.0`; no other
  package imports gonum directly.
- Unsupported methods/options are rejected before a backend call
  (`unsupported_option`). Backend panics are recovered and mapped to
  `library_failure`. NaN/Inf inputs and singular/non-positive-definite
  systems are rejected with `invalid_value` / `out_of_bounds`.

## 7. Resource limits and failure atomicity

`pipeline.Limits` caps variants, samples, cohorts, per-file bytes and wall
time. A thread-safe `Budget` refuses excess consumption (`limit_exceeded`;
`truncate` is an explicit, marked mode). A `Result` is assembled only after
every stage succeeds, so a failed run leaves no half-populated result;
parsers always close their file handles via `defer`, and parallel workers
are joined through `WaitGroup` before returning.

## 8. Known scope and non-goals

- Only bi-allelic SNP/indel identifiers from local delimited files; no VCF,
  BGEN or network fetchers.
- Overlap flagging is based on declared identifiers/keys, not genomic
  fingerprint reconstruction.
- LD-score regression intercept is unconstrained by default; h2 estimates
  with pathological intercepts (e.g. strongly negative on toy data) are
  reported numerically with full provenance, not declared causal truth.
