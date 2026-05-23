# highMLR 1.0.1

* Default feature set in `highmlr()` now excludes non-numeric columns
  (e.g. identifiers) automatically, preventing an error when `features`
  is left unspecified.
* The `coxnet` backend now validates that supplied features are numeric
  and reports a clear error otherwise.
* Learner backends (`rsf`, `aorsf`, `xgboost`) now default to at most two
  threads, configurable via `options()`, to respect multi-core policies
  in batch-checking environments.

# highMLR 1.0.0

Major rewrite. This release supersedes the original CRAN version 0.1.1
and is **not backward-compatible** with it. The previous six functions
(`mlhighCox`, `mlhighKap`, `mlhighFrail`, `mlhighHet`, `mlclassCox`,
`mlclassKap`) and their `mlr3`/`coxme`/`missForest` backend have been
replaced by a single, unified interface. Code written for 0.1.1 will
not run unchanged.

## New unified interface

* Single `highmlr()` entry point dispatching to eight feature-selection
  methods, all returning a common `highmlr_fit` S3 object:
  Cox elastic net (`coxnet`), random survival forest (`rsf`),
  accelerated oblique RSF (`aorsf`), gradient-boosted Cox (`xgboost`),
  stability selection (`stability`), univariate Cox screening
  (`univariate`), pseudo-observation bridging (`pseudo`), and Fine-Gray
  competing-risks selection (`finegray`).
* Common verbs for the fit object: `print()`, `summary()`, `plot()`,
  `coef()`, `predict()`.
* Companion functions: `highmlr_compare()`, `highmlr_stability()`,
  `highmlr_explain()`, `highmlr_screen()`, `highmlr_report()`.

## Feature-selection methods

* `method = "pseudo"` -- Andersen & Perme (2010) pseudo-observation
  bridge. Computes jackknife pseudo-values for S(t) at chosen evaluation
  times, then fits any regression learner (`ranger`, `xgboost`,
  `glmnet`) on the transformed outcome, allowing arbitrary regression ML
  on right-censored data without the proportional hazards assumption.
* `method = "finegray"` -- Fine-Gray subdistribution hazard feature
  selection for competing-risks data. Status coded 0/1/2+ (cause of
  interest configurable). Returns subdistribution hazard ratios with
  FDR-adjusted p-values.

## Causal, conformal, and explanation tools

* `highmlr_causal()` -- **EXPERIMENTAL.** Causal survival forest via
  `grf` for heterogeneous treatment effects on survival. Returns
  per-patient CATE estimates (RMST or survival-probability difference)
  with standard errors, the average treatment effect, and covariate
  importance.
* `highmlr_conformal()` -- conformal prediction intervals for survival
  time with finite-sample marginal coverage (Candes, Lei & Ren, 2023),
  using inverse-probability-of-censoring weighting on the calibration
  set.
* `highmlr_explain()` -- time-dependent SHAP via SurvSHAP(t) (Krzyzinski
  et al., 2023) by default, returning SHAP values per feature, per
  patient, across time, with aggregated mean |SHAP| importance.
  Permutation and break-down methods are also available.

## Dependencies

* Imports: `survival`, `glmnet`, `ranger`, `aorsf`, `xgboost`, `stabs`,
  `survex`, `grf`, `prodlim`, `cmprsk`, `future`, `future.apply`,
  `tibble`, `ggplot2`, `rlang`, `stats`, `utils`.
* Suggests: `knitr`, `rmarkdown`, `testthat`, `mice`, `riskRegression`.
* The `mlr3`, `mlr3learners`, `gtools`, `dplyr`, `coxme`, `missForest`,
  and `R6` dependencies from 0.1.1 are no longer used.

## Notes

* Neural network methods (DeepSurv, DeepHit) are deferred to a future
  release to keep this version free of torch/Python dependencies.
* Datasets `hnscc` and `srdata` are retained from 0.1.1, now with
  expanded and corrected documentation.
