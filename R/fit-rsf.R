# Random Survival Forest backend for highmlr().
# Uses ranger directly (fast, CRAN-clean, permutation importance).

fit_rsf <- function(data, time, status, features,
                    engine = "ranger",
                    top_n = 50L,
                    resampling = "cv",
                    folds = 5L,
                    num.trees = 500L,
                    mtry = NULL,
                    min.node.size = 15L,
                    importance = "permutation",
                    ...) {

  df <- data[, c(time, status, features), drop = FALSE]

  # ranger needs no NAs in predictors; mean-impute numerics, mode-impute factors
  df <- impute_simple(df, features)

  fml <- stats::as.formula(
    sprintf("survival::Surv(%s, %s) ~ .", time, status)
  )

  if (is.null(mtry)) {
    mtry <- max(1L, floor(sqrt(length(features))))
  }

  fit <- ranger::ranger(
    formula       = fml,
    data          = df,
    num.trees     = num.trees,
    mtry          = mtry,
    min.node.size = min.node.size,
    importance    = importance,
    splitrule     = "logrank",
    num.threads   = getOption("ranger.num.threads", 2L),
    ...
  )

  imp <- fit$variable.importance
  imp <- imp[order(-imp)]
  imp <- utils::head(imp, top_n)

  selected <- tibble::tibble(
    feature    = names(imp),
    importance = unname(imp)
  )

  performance <- list(
    oob_error = unname(fit$prediction.error),  # 1 - Harrell's C on OOB
    c_index   = 1 - unname(fit$prediction.error),
    num.trees = num.trees,
    mtry      = mtry
  )

  new_highmlr_fit(
    selected    = selected,
    performance = performance,
    model       = list(fit = fit, features = features,
                       imputation = attr(df, "imputation")),
    meta        = list(engine = engine, importance = importance)
  )
}

# Shared predict helper for parsnip-style survival models
predict_parsnip_surv <- function(model, new_data, type, ...) {
  fit <- model$fit
  feats <- model$features
  miss <- setdiff(feats, names(new_data))
  if (length(miss)) {
    rlang::abort(sprintf("new_data missing %d feature(s).", length(miss)))
  }
  new_data <- new_data[, feats, drop = FALSE]
  new_data <- apply_imputation(new_data, model$imputation)

  if (inherits(fit, "ranger")) {
    pred <- stats::predict(fit, data = new_data)
    return(switch(type,
      linear_pred = -rowMeans(pred$chf),  # higher chf = higher risk
      risk        = rowMeans(pred$chf),
      survival    = pred$survival
    ))
  }
  if (inherits(fit, "orsf_fit")) {
    return(switch(type,
      risk        = as.numeric(stats::predict(fit, new_data = new_data,
                                              pred_type = "risk")),
      linear_pred = as.numeric(stats::predict(fit, new_data = new_data,
                                              pred_type = "risk")),
      survival    = stats::predict(fit, new_data = new_data,
                                   pred_type = "surv")
    ))
  }
  if (inherits(fit, "xgb.Booster")) {
    Xn <- as.matrix(new_data)
    return(switch(type,
      risk        = stats::predict(fit, newdata = Xn),
      linear_pred = log(stats::predict(fit, newdata = Xn)),
      survival    = rlang::abort("type='survival' not supported for xgboost.")
    ))
  }
  rlang::abort("Unknown model class in predict_parsnip_surv().")
}

# Simple imputation utility: numeric -> mean, factor/character -> mode
impute_simple <- function(df, features) {
  info <- list()
  for (f in features) {
    x <- df[[f]]
    if (is.numeric(x)) {
      m <- mean(x, na.rm = TRUE)
      df[[f]][is.na(x)] <- m
      info[[f]] <- list(type = "numeric", value = m)
    } else {
      tab <- sort(table(x), decreasing = TRUE)
      m   <- names(tab)[1]
      df[[f]][is.na(x)] <- m
      info[[f]] <- list(type = "categorical", value = m)
    }
  }
  attr(df, "imputation") <- info
  df
}

apply_imputation <- function(df, info) {
  if (is.null(info)) return(df)
  for (f in names(info)) {
    if (f %in% names(df) && anyNA(df[[f]])) {
      df[[f]][is.na(df[[f]])] <- info[[f]]$value
    }
  }
  df
}
