#' Machine learning feature selection for high dimensional survival data
#'
#' Fits one of several survival ML methods and returns a unified
#' `highmlr_fit` object summarising the selected features, their
#' importance/coefficients, and (optionally) out-of-sample performance.
#'
#' @param data A data frame containing `time`, `status`, and the candidate
#'   features (or a superset). Rows with missing time/status are dropped.
#' @param time Character scalar: name of the survival time column.
#' @param status Character scalar: name of the event indicator column.
#'   For right-censored methods: 1 = event, 0 = censored.
#'   For Fine-Gray (method = "finegray"): 0 = censored, 1 = event of
#'   interest, 2+ = competing event(s).
#' @param features Character vector of candidate feature column names.
#'   If `NULL` (default), all columns except `time` and `status` are used.
#' @param method One of `"coxnet"`, `"rsf"`, `"aorsf"`, `"xgboost"`,
#'   `"stability"`, `"univariate"`, `"pseudo"`, `"finegray"`.
#' @param engine Optional engine override.
#' @param recipe Optional preprocessing recipe object (currently accepted
#'   for forward compatibility; not yet applied).
#' @param resampling One of `"cv"`, `"bootstrap"`, `"holdout"`, `"none"`.
#' @param folds Integer, number of CV folds (default 5).
#' @param tune Logical. Internal tuning (currently coxnet only).
#' @param top_n Integer. For ranking-based methods, keep this many top
#'   features (default 50).
#' @param parallel Logical. Use future-based parallelism for the
#'   embarrassingly parallel parts.
#' @param seed Optional integer for reproducibility.
#' @param ... Additional arguments passed to the method-specific fitter.
#'
#' @return An object of class `highmlr_fit`. See [new_highmlr_fit()].
#'
#' @examples
#' \donttest{
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   data(hnscc)
#'   fit <- highmlr(hnscc, time = "OS", status = "Death",
#'                  method = "coxnet", resampling = "cv", folds = 5)
#'   print(fit)
#' }
#' }
#'
#' @export
highmlr <- function(data,
                    time,
                    status,
                    features  = NULL,
                    method    = c("coxnet", "rsf", "aorsf", "xgboost",
                                  "stability", "univariate",
                                  "pseudo", "finegray"),
                    engine    = NULL,
                    recipe    = NULL,
                    resampling = c("cv", "bootstrap", "holdout", "none"),
                    folds     = 5L,
                    tune      = FALSE,
                    top_n     = 50L,
                    parallel  = FALSE,
                    seed      = NULL,
                    ...) {

  call       <- match.call()
  method     <- match.arg(method)
  resampling <- match.arg(resampling)

  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame.")
  }
  if (!is.character(time)   || length(time)   != 1L ||
      !is.character(status) || length(status) != 1L) {
    rlang::abort("`time` and `status` must be single column names.")
  }
  if (!all(c(time, status) %in% names(data))) {
    rlang::abort(sprintf("Columns '%s' and/or '%s' not found in data.",
                         time, status))
  }
  if (is.null(features)) {
    candidates <- setdiff(names(data), c(time, status))
    is_num     <- vapply(data[candidates], is.numeric, logical(1))
    features   <- candidates[is_num]
    if (!length(features)) {
      rlang::abort("No numeric candidate features found after excluding `time` and `status`.")
    }
  }
  missing_feats <- setdiff(features, names(data))
  if (length(missing_feats)) {
    rlang::abort(sprintf("Missing feature columns: %s",
                         paste(missing_feats, collapse = ", ")))
  }

  keep <- !is.na(data[[time]]) & !is.na(data[[status]])
  data <- data[keep, c(time, status, features), drop = FALSE]

  data[[status]] <- as.numeric(data[[status]])
  if (method == "finegray") {
    # Allow 0, 1, 2+ for competing risks
    if (any(data[[status]] < 0)) {
      rlang::abort("`status` for Fine-Gray must be non-negative integers.")
    }
  } else if (!all(data[[status]] %in% c(0, 1))) {
    rlang::abort("`status` must be coded 0 (censored) / 1 (event).")
  }

  data_summary <- list(
    n          = nrow(data),
    p          = length(features),
    events     = sum(data[[status]] == 1L),
    event_rate = mean(data[[status]] == 1L)
  )
  if (method == "finegray") {
    data_summary$competing_events <- sum(data[[status]] >= 2L)
  }

  if (!is.null(seed)) set.seed(seed)

  if (parallel) {
    old_plan <- future::plan(future::multisession)
    on.exit(future::plan(old_plan), add = TRUE)
  }

  fit <- switch(method,
    coxnet     = fit_coxnet(data, time, status, features,
                            recipe = recipe,
                            resampling = resampling, folds = folds,
                            tune = tune, ...),
    rsf        = fit_rsf(data, time, status, features,
                         engine = engine %||% "ranger",
                         top_n = top_n,
                         resampling = resampling, folds = folds, ...),
    aorsf      = fit_aorsf(data, time, status, features,
                           top_n = top_n,
                           resampling = resampling, folds = folds, ...),
    xgboost    = fit_xgboost(data, time, status, features,
                             top_n = top_n,
                             resampling = resampling, folds = folds, ...),
    stability  = fit_stability(data, time, status, features,
                               parallel = parallel, ...),
    univariate = fit_univariate(data, time, status, features,
                                top_n = top_n, parallel = parallel, ...),
    pseudo     = fit_pseudo(data, time, status, features,
                            engine = engine %||% "ranger",
                            top_n = top_n,
                            parallel = parallel, ...),
    finegray   = fit_finegray(data, time, status, features,
                              top_n = top_n,
                              parallel = parallel, ...)
  )

  fit$method       <- method
  fit$call         <- call
  fit$data_summary <- data_summary
  fit
}
