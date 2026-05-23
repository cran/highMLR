# Cox elastic net (glmnet) backend for highmlr().
# Internal -- not exported. See ?highmlr for user-facing docs.

fit_coxnet <- function(data, time, status, features,
                       recipe = NULL,
                       resampling = "cv",
                       folds = 5L,
                       tune = FALSE,
                       alpha = 1,
                       nlambda = 100L,
                       s = c("lambda.min", "lambda.1se"),
                       ...) {

  s <- match.arg(s)

  # Predictor matrix and survival response
  non_numeric <- features[!vapply(data[features], is.numeric, logical(1))]
  if (length(non_numeric)) {
    rlang::abort(sprintf(
      "method = 'coxnet' requires numeric features; non-numeric column(s): %s",
      paste(utils::head(non_numeric, 5), collapse = ", ")))
  }
  X <- as.matrix(data[, features, drop = FALSE])
  y <- survival::Surv(data[[time]], data[[status]])

  # Drop zero-variance columns (glmnet errors otherwise)
  sds  <- apply(X, 2, stats::sd, na.rm = TRUE)
  keep <- which(!is.na(sds) & sds > 0)
  if (!length(keep)) {
    rlang::abort("All features have zero variance after preprocessing.")
  }
  X <- X[, keep, drop = FALSE]
  feature_names <- colnames(X)

  # Mean-impute remaining NAs (glmnet does not accept NA)
  if (anyNA(X)) {
    col_means <- colMeans(X, na.rm = TRUE)
    na_idx <- which(is.na(X), arr.ind = TRUE)
    X[na_idx] <- col_means[na_idx[, 2]]
  }

  # Fit penalised Cox via cv.glmnet (handles lambda selection internally)
  cv_folds <- if (resampling == "cv") folds else 10L
  fit <- glmnet::cv.glmnet(
    x        = X,
    y        = y,
    family   = "cox",
    alpha    = alpha,
    nfolds   = cv_folds,
    nlambda  = nlambda,
    standardize = TRUE,
    ...
  )

  # Extract coefficients at chosen lambda
  beta <- as.numeric(stats::coef(fit, s = s))
  names(beta) <- feature_names

  nonzero <- beta[beta != 0]
  selected <- tibble::tibble(
    feature       = names(nonzero),
    coef          = unname(nonzero),
    hazard_ratio  = exp(unname(nonzero)),
    importance    = abs(unname(nonzero))
  )
  selected <- selected[order(-selected$importance), ]

  # Performance: partial-likelihood deviance at chosen lambda
  lambda_val <- fit[[s]]
  cvm_idx    <- which.min(abs(fit$lambda - lambda_val))
  performance <- list(
    partial_deviance = unname(fit$cvm[cvm_idx]),
    lambda           = lambda_val,
    n_nonzero        = nrow(selected)
  )

  # If user asked for resampling-based C-index, compute it
  if (resampling != "none") {
    c_idx <- tryCatch(
      coxnet_c_index(X, y, alpha, s, folds = folds, nlambda = nlambda),
      error = function(e) NA_real_
    )
    performance$c_index <- c_idx
  }

  # Compute Breslow baseline survival on the training data so we can
  # build survival curves later (used by highmlr_explain SurvSHAP path).
  baseline_haz <- tryCatch({
    lp_train <- as.numeric(stats::predict(fit, newx = X,
                                          s = lambda_val, type = "link"))
    bh <- compute_breslow_baseline(time = y[, 1], status = y[, 2],
                                   lp = lp_train)
    bh
  }, error = function(e) NULL)

  new_highmlr_fit(
    selected    = selected,
    performance = performance,
    model       = list(cv_fit = fit, s = s,
                       feature_names = feature_names,
                       col_means = colMeans(X),
                       baseline_haz = baseline_haz),
    meta        = list(alpha = alpha, s = s,
                       lambda_min = fit$lambda.min,
                       lambda_1se = fit$lambda.1se)
  )
}

# Breslow baseline cumulative hazard / survival.
# Returns a list with `time`, `cumhaz`, and `surv` (Breslow KM-like).
compute_breslow_baseline <- function(time, status, lp) {
  ord <- order(time)
  t_o <- time[ord]; d_o <- status[ord]; lp_o <- lp[ord]
  exp_lp <- exp(lp_o)
  # At each event time, increment = d_i / sum_{j: t_j >= t_i} exp(lp_j)
  unique_t <- sort(unique(t_o[d_o == 1L]))
  if (!length(unique_t)) return(NULL)
  cumhaz <- numeric(length(unique_t))
  for (k in seq_along(unique_t)) {
    tk <- unique_t[k]
    d_k <- sum(d_o == 1L & t_o == tk)
    risk_set <- sum(exp_lp[t_o >= tk])
    if (risk_set > 0) cumhaz[k] <- d_k / risk_set
  }
  cumhaz <- cumsum(cumhaz)
  list(time = unique_t, cumhaz = cumhaz, surv = exp(-cumhaz))
}

# Out-of-sample C-index by re-doing CV with held-out predictions
coxnet_c_index <- function(X, y, alpha, s, folds, nlambda) {

  n     <- nrow(X)
  fold  <- sample(rep_len(seq_len(folds), n))
  lp    <- numeric(n)

  for (k in seq_len(folds)) {
    tr <- which(fold != k)
    te <- which(fold == k)
    cvf <- glmnet::cv.glmnet(X[tr, , drop = FALSE], y[tr, ],
                             family = "cox", alpha = alpha,
                             nlambda = nlambda, standardize = TRUE)
    lp[te] <- as.numeric(
      stats::predict(cvf, newx = X[te, , drop = FALSE],
                     s = s, type = "link")
    )
  }
  # Harrell's C via survival::concordance on the linear predictor
  cc <- survival::concordance(y ~ lp, reverse = TRUE)
  unname(cc$concordance)
}

# predict.highmlr_fit dispatches here for method = "coxnet"
predict_coxnet <- function(model, new_data, type, ...) {
  cv_fit  <- model$cv_fit
  s       <- model$s
  feats   <- model$feature_names

  miss <- setdiff(feats, names(new_data))
  if (length(miss)) {
    rlang::abort(sprintf("new_data is missing %d feature column(s): %s",
                         length(miss),
                         paste(utils::head(miss, 5), collapse = ", ")))
  }
  Xn <- as.matrix(new_data[, feats, drop = FALSE])
  if (anyNA(Xn)) {
    cm <- model$col_means
    idx <- which(is.na(Xn), arr.ind = TRUE)
    Xn[idx] <- cm[idx[, 2]]
  }

  glmnet_type <- switch(type,
    linear_pred = "link",
    risk        = "response",
    survival    = rlang::abort("type = 'survival' not yet implemented for coxnet."),
  )
  out <- stats::predict(cv_fit, newx = Xn, s = s, type = glmnet_type)
  as.numeric(out)
}
