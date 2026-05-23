# Pseudo-observation bridge to arbitrary regression learners.
#
# Computes jackknife pseudo-observations for the survival probability at
# one or more chosen time points (Andersen & Perme, 2010), then fits a
# user-chosen regression learner (ranger, xgboost, glmnet) on the pseudo-
# values. This lets the user apply any standard regression ML algorithm
# to right-censored survival data and obtain feature importances /
# coefficients in the regression sense, with no proportional hazards
# assumption.

fit_pseudo <- function(data, time, status, features,
                       engine = c("ranger", "xgboost", "glmnet"),
                       eval_times = NULL,
                       top_n = 50L,
                       parallel = FALSE,
                       ...) {

  engine <- match.arg(engine)

  df <- data[, c(time, status, features), drop = FALSE]
  df <- impute_simple(df, features)

  t_vec <- df[[time]]
  s_vec <- df[[status]]

  # Default evaluation times: 25%, 50%, 75% quantiles of observed event times
  if (is.null(eval_times)) {
    event_times <- t_vec[s_vec == 1L]
    if (length(event_times) < 3L) {
      eval_times <- stats::quantile(t_vec, c(0.25, 0.5, 0.75))
    } else {
      eval_times <- stats::quantile(event_times, c(0.25, 0.5, 0.75))
    }
    eval_times <- unname(eval_times)
  }

  # Compute pseudo-observations via prodlim::jackknife
  if (!requireNamespace("prodlim", quietly = TRUE)) {
    rlang::abort("Package 'prodlim' required for pseudo method.")
  }
  fit_km <- prodlim::prodlim(
    survival::Surv(t_vec, s_vec) ~ 1,
    data = data.frame(t_vec = t_vec, s_vec = s_vec)
  )
  # pseudo: n x length(eval_times); each column is the pseudo-value of S(t)
  pseudo_mat <- prodlim::jackknife(fit_km, times = eval_times,
                                   cause = 1)
  if (is.null(dim(pseudo_mat))) pseudo_mat <- matrix(pseudo_mat, ncol = 1L)
  colnames(pseudo_mat) <- paste0("t", round(eval_times, 2))

  X <- as.matrix(df[, features, drop = FALSE])

  # Fit one regression learner per evaluation time, average importance
  one_time_fit <- function(j) {
    yj <- pseudo_mat[, j]
    fit_one_pseudo_regressor(X, yj, engine = engine, ...)
  }

  per_time <- if (parallel) {
    future.apply::future_lapply(seq_len(ncol(pseudo_mat)), one_time_fit,
                                future.seed = TRUE)
  } else {
    lapply(seq_len(ncol(pseudo_mat)), one_time_fit)
  }

  # Aggregate importances across time points (mean)
  imp_mat <- do.call(cbind, lapply(per_time, `[[`, "importance"))
  rownames(imp_mat) <- features
  mean_imp <- rowMeans(imp_mat, na.rm = TRUE)

  ord <- order(-abs(mean_imp))
  ord <- ord[seq_len(min(top_n, length(ord)))]

  selected <- tibble::tibble(
    feature    = features[ord],
    importance = unname(mean_imp[ord])
  )

  performance <- list(
    eval_times = eval_times,
    engine     = engine,
    n_times    = ncol(pseudo_mat)
  )

  new_highmlr_fit(
    selected    = selected,
    performance = performance,
    model       = list(per_time_fits = per_time,
                       eval_times    = eval_times,
                       engine        = engine,
                       features      = features,
                       imputation    = attr(df, "imputation")),
    meta        = list(engine = engine, eval_times = eval_times)
  )
}

# Helper: fit a single regression learner on pseudo-values
fit_one_pseudo_regressor <- function(X, y, engine, ...) {
  switch(engine,
    ranger = {
      df_local <- data.frame(y = y, X)
      g <- ranger::ranger(y ~ ., data = df_local,
                          num.trees = 300, importance = "permutation",
                          ...)
      list(fit = g, importance = g$variable.importance)
    },
    xgboost = {
      dtrain <- xgboost::xgb.DMatrix(X, label = y)
      g <- xgboost::xgb.train(
        params = list(objective = "reg:squarederror", eta = 0.1,
                      max_depth = 4),
        data = dtrain, nrounds = 100, verbose = 0
      )
      imp_df <- xgboost::xgb.importance(model = g)
      imp <- stats::setNames(rep(0, ncol(X)), colnames(X))
      imp[imp_df$Feature] <- imp_df$Gain
      list(fit = g, importance = imp)
    },
    glmnet = {
      g <- glmnet::cv.glmnet(X, y, alpha = 1, standardize = TRUE)
      beta <- as.numeric(stats::coef(g, s = "lambda.min")[-1L])
      names(beta) <- colnames(X)
      list(fit = g, importance = abs(beta))
    }
  )
}
