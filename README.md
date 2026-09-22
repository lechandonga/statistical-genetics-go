# statistical-genetics-go

A local, deterministic statistical-genetics analysis pipeline in Go:
LD-score heritability, multi-cohort meta-analysis and individual polygenic
risk scores from local files, with explicit allele-strand alignment and
distinguishable rejection reasons.

- Canonical variant model; no silently skipped rows (`pkg/statgen/variant`,
  `pkg/statgen/sumstats`, `pkg/statgen/ldscore`, `pkg/statgen/genotype`)
- Allele swap/strand-flip alignment with explicit A/T & C/G ambiguity
  policy (never silent sign addition)
- LD-score regression with uncertainty and explicit negative/out-of-range
  status (`pkg/statgen/ldsc`)
- Fixed/random-effect meta-analysis with Q, I2, tau2 and sample-overlap /
  duplicate-cohort detection (`pkg/statgen/meta`)
- PRS with documented sign convention, missing-data policies and AUC/R2
  evaluation (`pkg/statgen/prs`)
- Single pinned gonum seam (`pkg/statgen/numlib`); panics and unsupported
  options become structured errors (`pkg/statgen/stgerr`)
- Serial/parallel bit-identical results; configurable size/time budgets with
  failure atomicity (`pkg/statgen/pipeline`)

See [`docs/STATISTICAL_MODEL.md`](docs/STATISTICAL_MODEL.md) for model
assumptions and allele-alignment rules, and
[`docs/USAGE.md`](docs/USAGE.md) for formats, flags and verification.

Quick start:

```sh
go test ./...
go run ./cmd/statgen \
  -sumstats testdata/cohort_A.tsv,testdata/cohort_B.tsv -names A,B \
  -ldscores testdata/ldscores.tsv -genotypes testdata/genotypes.tsv \
  -workers 4 -out result.json
```
