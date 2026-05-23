# Per-feature univariate Cox PH screening.
# Equivalent to the v0.1.x style approach; retained as a baseline.

fit_univariate <- function(data, time, status, features,
                           top_n = 50L,
                           rank_by = c("p_value", "logLik", "abs_coef"),
                           parallel = FALSE,
                           ...) {

  rank_by <- match.arg(rank_by)
  df <- data[, c(time, status, features), drop = FALSE]
  df <- impute_simple(df, features)

  one_feature <- function(f) {
    x <- df[[f]]
    if (stats::sd(x, na.rm = TRUE) == 0) {
      return(c(coef = NA, hr = NA, se = NA, z = NA, p = NA, logLik = NA))
    }
    fit <- tryCatch(
      survival::coxph(survival::Surv(df[[time]], df[[status]]) ~ x,
                      data = df),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(c(coef = NA, hr = NA, se = NA, z = NA, p = NA, logLik = NA))
    }
    s <- summary(fit)
    c(coef   = unname(s$coefficients[1, "coef"]),
      hr     = unname(s$coefficients[1, "exp(coef)"]),
      se     = unname(s$coefficients[1, "se(coef)"]),
      z      = unname(s$coefficients[1, "z"]),
      p      = unname(s$coefficients[1, "Pr(>|z|)"]),
      logLik = unname(stats::logLik(fit)))
  }

  results <- if (parallel) {
    do.call(rbind, future.apply::future_lapply(features, one_feature,
                                               future.seed = TRUE))
  } else {
    do.call(rbind, lapply(features, one_feature))
  }
  results <- as.data.frame(results)
  results$feature <- features

  results <- results[stats::complete.cases(results), ]

  # FDR-adjusted p-values (BH)
  results$p_adj <- stats::p.adjust(results$p, method = "BH")

  ord <- switch(rank_by,
    p_value  = order(results$p),
    logLik   = order(-results$logLik),
    abs_coef = order(-abs(results$coef))
  )
  results <- results[ord, ]
  results <- utils::head(results, top_n)

  selected <- tibble::tibble(
    feature       = results$feature,
    coef          = results$coef,
    hazard_ratio  = results$hr,
    se            = results$se,
    z             = results$z,
    p_value       = results$p,
    p_adjusted    = results$p_adj,
    importance    = -log10(results$p + 1e-300)
  )

  performance <- list(
    n_tested    = length(features),
    n_selected  = nrow(selected),
    rank_by     = rank_by,
    min_p       = min(results$p, na.rm = TRUE),
    min_p_adj   = min(results$p_adj, na.rm = TRUE)
  )

  new_highmlr_fit(
    selected    = selected,
    performance = performance,
    model       = list(results = results, features = features),
    meta        = list(rank_by = rank_by)
  )
}
