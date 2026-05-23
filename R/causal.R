#' Causal survival forest for heterogeneous treatment effects (experimental)
#'
#' Estimates patient-level conditional average treatment effects (CATEs)
#' on a survival outcome using `grf::causal_survival_forest`. Unlike the
#' rest of `highMLR`, this function answers a different question: not
#' "which features predict survival?" but "for which patients does
#' treatment T extend (or shorten) survival, and which features modify
#' that effect?".
#'
#' @section Experimental:
#' This function is marked experimental. The signature, defaults, and
#' return shape may change in a future release. Use with care in
#' published analyses, and report the package version.
#'
#' @param data A data frame.
#' @param time Character: name of the survival time column.
#' @param status Character: name of the event indicator (0/1).
#' @param treatment Character: name of the binary treatment column
#'   (0 = control, 1 = treated). Must be exactly two levels.
#' @param covariates Character vector of covariate column names. If
#'   `NULL`, all columns other than `time`, `status`, `treatment`.
#' @param horizon Numeric. The time horizon at which the treatment
#'   effect on the survival probability is estimated. Defaults to the
#'   median observed time.
#' @param num.trees Number of trees in the forest (default 2000).
#' @param target One of `"RMST"` (restricted mean survival time
#'   difference up to `horizon`) or `"survival.probability"`
#'   (difference in survival probability at `horizon`).
#' @param honesty Logical (default TRUE) -- honest splitting per `grf`.
#' @param seed Optional integer seed.
#' @param ... Passed to `grf::causal_survival_forest`.
#'
#' @return An object of class `highmlr_causal` containing the fitted
#'   forest, per-patient CATE estimates with standard errors, and
#'   covariate importance.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 500; p <- 10
#' X <- matrix(rnorm(n*p), n, p); colnames(X) <- paste0("V", 1:p)
#' W <- rbinom(n, 1, 0.5)
#' t <- rexp(n, rate = exp(0.3*W + 0.5*X[,1]*W))
#' c <- rexp(n, rate = 0.05)
#' d <- data.frame(OS = pmin(t,c), Death = as.integer(t<=c),
#'                 arm = W, X)
#' cf <- highmlr_causal(d, "OS", "Death", treatment = "arm",
#'                      covariates = paste0("V", 1:p))
#' print(cf); plot(cf)
#' }
#'
#' @export
highmlr_causal <- function(data,
                           time,
                           status,
                           treatment,
                           covariates = NULL,
                           horizon    = NULL,
                           num.trees  = 2000L,
                           target     = c("RMST", "survival.probability"),
                           honesty    = TRUE,
                           seed       = NULL,
                           ...) {

  target <- match.arg(target)

  if (!requireNamespace("grf", quietly = TRUE)) {
    rlang::abort("Package 'grf' required for highmlr_causal().")
  }
  if (!is.data.frame(data)) rlang::abort("`data` must be a data frame.")
  for (col in c(time, status, treatment)) {
    if (!col %in% names(data)) {
      rlang::abort(sprintf("Column '%s' not found in data.", col))
    }
  }
  if (is.null(covariates)) {
    covariates <- setdiff(names(data), c(time, status, treatment))
  }
  miss <- setdiff(covariates, names(data))
  if (length(miss)) {
    rlang::abort(sprintf("Missing covariates: %s",
                         paste(miss, collapse = ", ")))
  }

  keep <- !is.na(data[[time]]) & !is.na(data[[status]]) &
          !is.na(data[[treatment]])
  data <- data[keep, , drop = FALSE]

  W <- as.numeric(data[[treatment]])
  if (!all(W %in% c(0, 1)) || length(unique(W)) < 2L) {
    rlang::abort("`treatment` must be coded 0/1 with both levels present.")
  }

  Y <- data[[time]]
  D <- as.numeric(data[[status]])
  if (!all(D %in% c(0, 1))) {
    rlang::abort("`status` must be 0/1.")
  }

  X <- as.matrix(data[, covariates, drop = FALSE])
  X <- apply(X, 2, function(col) {
    if (anyNA(col)) col[is.na(col)] <- mean(col, na.rm = TRUE)
    col
  })

  if (is.null(horizon)) horizon <- stats::median(Y)

  if (!is.null(seed)) set.seed(seed)

  cf <- grf::causal_survival_forest(
    X        = X,
    Y        = Y,
    W        = W,
    D        = D,
    horizon  = horizon,
    target   = target,
    num.trees = num.trees,
    honesty  = honesty,
    ...
  )

  # Per-patient CATE with SE
  pred <- stats::predict(cf, estimate.variance = TRUE)
  cate <- as.numeric(pred$predictions)
  cate_se <- sqrt(pmax(0, as.numeric(pred$variance.estimates)))

  patient_effects <- tibble::tibble(
    cate    = cate,
    cate_se = cate_se,
    ci_lo   = cate - 1.96 * cate_se,
    ci_hi   = cate + 1.96 * cate_se
  )

  # Variable importance from grf
  vi <- grf::variable_importance(cf)
  imp <- tibble::tibble(
    covariate  = covariates,
    importance = as.numeric(vi)
  )
  imp <- imp[order(-imp$importance), ]

  # Average treatment effect (ATE)
  ate <- tryCatch(grf::average_treatment_effect(cf),
                  error = function(e) c(estimate = NA, std.err = NA))

  out <- list(
    forest          = cf,
    patient_effects = patient_effects,
    importance      = imp,
    ate             = ate,
    target          = target,
    horizon         = horizon,
    n               = length(Y),
    n_treated       = sum(W == 1L),
    n_events        = sum(D == 1L),
    covariates      = covariates,
    call            = match.call()
  )
  class(out) <- "highmlr_causal"
  out
}

#' @param x A `highmlr_causal` object.
#' @param n Number of top covariates to print (default 10).
#' @rdname highmlr_causal
#' @return `print()` invisibly returns `x`; `plot()` returns a
#'   `ggplot` object showing the distribution of estimated CATEs.
#' @export
print.highmlr_causal <- function(x, n = 10, ...) {
  cat("<highmlr_causal>  [EXPERIMENTAL]\n")
  cat("  Target:   ", x$target, " at horizon ", round(x$horizon, 3),
      "\n", sep = "")
  cat("  n:        ", x$n, "  (", x$n_treated, " treated, ",
      x$n - x$n_treated, " control, ", x$n_events, " events)\n", sep = "")
  if (!is.null(x$ate) && !any(is.na(x$ate))) {
    cat("  ATE:       ", signif(x$ate["estimate"], 3),
        "  (SE ", signif(x$ate["std.err"], 3), ")\n", sep = "")
  }
  cat("  CATE distribution:\n")
  print(round(stats::quantile(x$patient_effects$cate,
                              c(0.05, 0.25, 0.5, 0.75, 0.95)), 3))
  cat("\n  Top ", min(n, nrow(x$importance)),
      " covariates by importance:\n", sep = "")
  print(utils::head(x$importance, n))
  invisible(x)
}

#' @rdname highmlr_causal
#' @export
plot.highmlr_causal <- function(x, ...) {
  ggplot2::ggplot(x$patient_effects,
                  ggplot2::aes(x = .data$cate)) +
    ggplot2::geom_histogram(bins = 40, fill = "steelblue",
                            colour = "white") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey30") +
    ggplot2::labs(
      title = paste0("highMLR causal: distribution of CATEs (",
                     x$target, ")"),
      x = paste0("Estimated treatment effect (", x$target, ")"),
      y = "Patients"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}
