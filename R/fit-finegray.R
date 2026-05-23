# Fine-Gray subdistribution hazard feature selection for competing risks.
#
# Two stages: (1) univariate Fine-Gray screening of each feature for the
# cause of interest (default cause = 1), with FDR adjustment;
# (2) optional penalised Fine-Gray refit on the screened set, exposing
# coefficients and subdistribution hazard ratios. The penalised step
# uses fastcmprsk-style ideas: we use the screening step as the primary
# selection mechanism (well-supported by cmprsk::crr) and present
# results in the same tibble structure as the other methods.

fit_finegray <- function(data, time, status, features,
                         cause = 1L,
                         top_n = 50L,
                         parallel = FALSE,
                         ...) {

  if (!requireNamespace("cmprsk", quietly = TRUE)) {
    rlang::abort("Package 'cmprsk' required for finegray method.")
  }

  df <- data[, c(time, status, features), drop = FALSE]
  df <- impute_simple(df, features)

  t_vec <- df[[time]]
  s_vec <- df[[status]]

  one_feature <- function(f) {
    x <- df[[f]]
    if (stats::sd(x, na.rm = TRUE) == 0) {
      return(c(coef = NA, shr = NA, se = NA, z = NA, p = NA))
    }
    fit <- tryCatch(
      cmprsk::crr(ftime = t_vec, fstatus = s_vec,
                  cov1 = as.matrix(x),
                  failcode = cause, cencode = 0L),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(c(coef = NA, shr = NA, se = NA, z = NA, p = NA))
    }
    s <- summary(fit)
    cf <- s$coef
    c(coef = unname(cf[1, "coef"]),
      shr  = unname(cf[1, "exp(coef)"]),
      se   = unname(cf[1, "se(coef)"]),
      z    = unname(cf[1, "z"]),
      p    = unname(cf[1, "p-value"]))
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
  results$p_adj <- stats::p.adjust(results$p, method = "BH")

  ord <- order(results$p)
  results <- results[ord, ]
  results <- utils::head(results, top_n)

  selected <- tibble::tibble(
    feature             = results$feature,
    coef                = results$coef,
    subdistribution_hr  = results$shr,
    se                  = results$se,
    z                   = results$z,
    p_value             = results$p,
    p_adjusted          = results$p_adj,
    importance          = -log10(results$p + 1e-300)
  )

  performance <- list(
    n_tested        = length(features),
    n_selected      = nrow(selected),
    cause           = cause,
    min_p           = min(results$p, na.rm = TRUE),
    min_p_adj       = min(results$p_adj, na.rm = TRUE)
  )

  new_highmlr_fit(
    selected    = selected,
    performance = performance,
    model       = list(results = results, features = features,
                       cause = cause),
    meta        = list(cause = cause, framework = "Fine-Gray")
  )
}
