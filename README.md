# highMLR 1.0.0

Machine learning feature selection for high dimensional survival data.

## What's new in 1.0.0

Complete rewrite, superseding the original CRAN release v0.1.1. v0.1.1
offered six hardcoded functions, each fitting a single Cox-based model in
a feature-by-feature loop. This version replaces that with a single,
flexible entry point that dispatches to modern ML backends:

| Method        | Backend              | Notes                                  |
|---------------|----------------------|----------------------------------------|
| `coxnet`      | `glmnet` (Cox EN)    | Fast; baseline for p >> n              |
| `rsf`         | `ranger`             | Permutation importance                 |
| `aorsf`       | `aorsf`              | Accelerated oblique RSF (Jaeger 2024)  |
| `xgboost`     | `xgboost`            | Gradient-boosted Cox                   |
| `stability`   | `stabs` + `glmnet`   | Stability selection with PFER control  |
| `univariate`  | `survival::coxph`    | v0.1.x-style screening baseline        |
| `pseudo`      | pseudo-observations  | Bridge to any regression learner       |
| `finegray`    | `cmprsk` / `survival`| Competing-risks subdistribution hazard |

Plus companions: `highmlr_compare()`, `highmlr_stability()`,
`highmlr_explain()` (time-dependent SHAP via `survex`),
`highmlr_screen()`, `highmlr_report()`, and the additional tools
`highmlr_causal()` (causal survival forest via `grf`, experimental) and
`highmlr_conformal()` (conformal survival prediction intervals).

## Install (development)

```r
# install.packages("remotes")
remotes::install_local("path/to/highMLR")
```

## Quick start

```r
library(highMLR)
data(hnscc)

# One-liner: Cox elastic net with 5-fold CV
fit <- highmlr(hnscc, time = "OS", status = "Death", method = "coxnet")
print(fit)
plot(fit)

# Compare three methods
cmp <- highmlr_compare(hnscc, "OS", "Death",
                       methods = c("coxnet", "rsf", "xgboost"))
cmp$summary
```

## Breaking changes from 0.1.1

All v0.1.1 functions (`mlhighCox`, `mlhighKap`, `mlhighFrail`,
`mlhighHet`, `mlclassCox`, `mlclassKap`) are removed. Use `highmlr()`
with the appropriate `method` argument. The `univariate` method
reproduces the v0.1.x screening behaviour.
