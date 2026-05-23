# Gradient-boosted Cox backend via xgboost.

fit_xgboost <- function(data, time, status, features,
                        top_n = 50L,
                        resampling = "cv",
                        folds = 5L,
                        nrounds = 200L,
                        eta = 0.1,
                        max_depth = 4L,
                        early_stopping_rounds = 20L,
                        ...) {

  df <- data[, c(time, status, features), drop = FALSE]
  df <- impute_simple(df, features)

  # xgboost Cox encoding: time positive for events, negative for censored
  t <- df[[time]]
  s <- df[[status]]
  y <- ifelse(s == 1L, t, -t)
  X <- as.matrix(df[, features, drop = FALSE])

  dtrain <- xgboost::xgb.DMatrix(data = X, label = y)

  params <- list(
    objective   = "survival:cox",
    eval_metric = "cox-nloglik",
    eta         = eta,
    max_depth   = max_depth,
    nthread     = getOption("xgboost.nthread", 2L)
  )

  # CV to determine best nrounds, then refit
  cv <- xgboost::xgb.cv(
    params                = params,
    data                  = dtrain,
    nrounds               = nrounds,
    nfold                 = folds,
    early_stopping_rounds = early_stopping_rounds,
    verbose               = 0,
    ...
  )
  best <- cv$best_iteration %||% nrounds

  fit <- xgboost::xgb.train(
    params  = params,
    data    = dtrain,
    nrounds = best,
    verbose = 0
  )

  imp <- xgboost::xgb.importance(model = fit)
  imp <- imp[order(-imp$Gain), ]
  imp <- utils::head(imp, top_n)

  selected <- tibble::tibble(
    feature    = imp$Feature,
    importance = imp$Gain,
    cover      = imp$Cover,
    frequency  = imp$Frequency
  )

  performance <- list(
    best_iteration = best,
    cv_nloglik     = unname(cv$evaluation_log$test_cox_nloglik_mean[best]),
    nrounds        = nrounds
  )

  new_highmlr_fit(
    selected    = selected,
    performance = performance,
    model       = list(fit = fit, features = features,
                       imputation = attr(df, "imputation")),
    meta        = list(eta = eta, max_depth = max_depth, best_nrounds = best)
  )
}
