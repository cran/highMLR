# Accelerated Oblique Random Survival Forest backend.
# Jaeger et al. (2024), aorsf package on CRAN.

fit_aorsf <- function(data, time, status, features,
                      top_n = 50L,
                      resampling = "cv",
                      folds = 5L,
                      n_tree = 500L,
                      mtry = NULL,
                      importance = c("negate", "permute", "anova"),
                      ...) {

  importance <- match.arg(importance)
  df <- data[, c(time, status, features), drop = FALSE]
  df <- impute_simple(df, features)

  if (is.null(mtry)) {
    mtry <- max(1L, floor(sqrt(length(features))))
  }

  fit <- aorsf::orsf(
    data       = df,
    formula    = stats::as.formula(sprintf("%s + %s ~ .", time, status)),
    n_tree     = n_tree,
    mtry       = mtry,
    importance = importance,
    n_thread   = getOption("aorsf.n_thread", 2L),
    ...
  )

  imp <- aorsf::orsf_vi(fit)
  imp <- imp[order(-imp)]
  imp <- utils::head(imp, top_n)

  selected <- tibble::tibble(
    feature    = names(imp),
    importance = unname(imp)
  )

  performance <- list(
    oob_c_index = tryCatch(unname(fit$eval_oobag$stat_values[1]),
                           error = function(e) NA_real_),
    n_tree      = n_tree,
    mtry        = mtry
  )

  new_highmlr_fit(
    selected    = selected,
    performance = performance,
    model       = list(fit = fit, features = features,
                       imputation = attr(df, "imputation")),
    meta        = list(importance = importance)
  )
}
