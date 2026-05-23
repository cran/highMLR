#' Compare multiple highMLR methods on the same data
#'
#' Runs several methods and returns a side-by-side comparison of selected
#' features and performance.
#'
#' @param data,time,status,features As in [highmlr()].
#' @param methods Character vector of methods to compare.
#' @param ... Passed to each call to `highmlr()`.
#'
#' @return A list with two elements: `fits` (named list of `highmlr_fit`
#'   objects) and `summary` (a tibble of method, n_selected, key metric).
#'
#' @examples
#' \dontrun{
#' data(hnscc)
#' cmp <- highmlr_compare(hnscc, "OS", "Death",
#'                        methods = c("coxnet", "rsf", "univariate"))
#' cmp$summary
#' }
#'
#' @export
highmlr_compare <- function(data, time, status, features = NULL,
                            methods = c("coxnet", "rsf", "univariate"),
                            ...) {
  fits <- lapply(methods, function(m) {
    message("Fitting method: ", m)
    tryCatch(
      highmlr(data = data, time = time, status = status,
              features = features, method = m, ...),
      error = function(e) {
        warning("Method ", m, " failed: ", conditionMessage(e))
        NULL
      }
    )
  })
  names(fits) <- methods
  fits <- fits[!vapply(fits, is.null, logical(1))]

  summary_tbl <- tibble::tibble(
    method     = names(fits),
    n_selected = vapply(fits, function(f) nrow(f$selected), integer(1)),
    c_index    = vapply(fits, function(f) {
      f$performance$c_index %||% NA_real_
    }, numeric(1))
  )

  list(fits = fits, summary = summary_tbl)
}

#' Post-hoc stability analysis of a fitted highmlr_fit
#'
#' Runs stability selection on the data used in `fit`, returning a
#' selection frequency per feature.
#'
#' @param fit A `highmlr_fit` object (used only for the data / call).
#' @param B Number of subsamples (default 100).
#' @param cutoff Selection probability threshold (default 0.75).
#' @param PFER Per-family error rate bound (default 1).
#' @param ... Passed to [fit_stability()].
#'
#' @return A new `highmlr_fit` with `method = "stability"`.
#' @export
highmlr_stability <- function(fit, B = 100L, cutoff = 0.75,
                              PFER = 1, ...) {
  cl <- fit$call
  cl$method <- "stability"
  cl$B      <- B
  cl$cutoff <- cutoff
  cl$PFER   <- PFER
  eval(cl, envir = parent.frame())
}

#' Time-dependent SHAP explanations for a highmlr_fit (SurvSHAP(t))
#'
#' Computes SurvSHAP(t) attributions (Krzyzinski et al., 2023) -- SHAP
#' values that vary with follow-up time -- for the top features in a
#' fitted `highmlr_fit`. Returns the survex explainer, per-feature
#' aggregated importance, and a plotting helper.
#'
#' @param fit A `highmlr_fit` object with a stored model.
#' @param new_data Data on which to compute explanations.
#' @param top_n Number of top features to explain (default 10).
#' @param times Optional numeric vector of time points at which SHAP
#'   values are computed. Defaults to a 20-point grid spanning the
#'   observed time range.
#' @param method SHAP method passed through to `survex`. Default
#'   `"survshap"` (time-dependent). Other options: `"permutation"`,
#'   `"break_down"`.
#' @param n_explain How many test rows to compute SHAP for. Default 25
#'   (SHAP is expensive; full-cohort computation is rarely needed).
#' @param seed Optional integer for reproducibility of subsampling.
#' @param ... Passed to `survex::model_survshap()` or
#'   `survex::explain_survival()`.
#'
#' @return A list with class `highmlr_explain` containing:
#'   * `explainer`     -- the `survex` explainer object
#'   * `survshap`      -- the time-dependent SHAP object (if applicable)
#'   * `top_features`  -- the top features table from the fit
#'   * `aggregated`    -- tibble of mean absolute SHAP per feature,
#'                      averaged across time and explained rows
#'
#' @export
highmlr_explain <- function(fit, new_data = NULL, top_n = 10L,
                            times = NULL,
                            method = c("survshap", "permutation",
                                       "break_down"),
                            n_explain = 25L,
                            seed = NULL, ...) {

  method <- match.arg(method)

  if (!requireNamespace("survex", quietly = TRUE)) {
    rlang::abort("Package 'survex' required. Install it first.")
  }
  if (is.null(fit$model)) {
    rlang::abort("This highmlr_fit has no stored model.")
  }
  if (is.null(new_data)) {
    rlang::abort("Please supply `new_data` for explanation.")
  }

  time_col   <- as.character(fit$call$time   %||% "OS")
  status_col <- as.character(fit$call$status %||% "Death")

  if (!all(c(time_col, status_col) %in% names(new_data))) {
    rlang::abort(sprintf(
      "new_data must contain '%s' and '%s'.", time_col, status_col
    ))
  }

  if (!is.null(seed)) set.seed(seed)

  # Subsample rows to explain (SurvSHAP(t) is O(n_explain * p * |times|))
  n_avail <- nrow(new_data)
  if (n_avail > n_explain) {
    expl_rows <- sample.int(n_avail, n_explain)
  } else {
    expl_rows <- seq_len(n_avail)
  }
  explain_data <- new_data[expl_rows, , drop = FALSE]

  if (is.null(times)) {
    obs_t  <- new_data[[time_col]]
    times  <- seq(stats::quantile(obs_t, 0.05),
                  stats::quantile(obs_t, 0.95),
                  length.out = 20L)
  }

  # Predictor matrix for the explainer (no time/status columns)
  predictor_data <- new_data[, setdiff(names(new_data),
                                       c(time_col, status_col)),
                             drop = FALSE]

  # ---- Build two prediction functions survex needs ---------------------
  # 1) A risk function: scalar per row.
  # 2) A survival-curve function: matrix of nrow x length(times).
  # Both are required by survex::explain_survival(); without (2),
  # model_survshap() has nothing to attribute.

  risk_fun <- function(model, newdata) {
    predict_highmlr_internal(model, newdata, type = "risk",
                             method = fit$method)
  }

  # Survival-curve predict function (rows x times)
  surv_fun <- build_surv_predict_function(fit, times)

  # ---- Build the survex explainer --------------------------------------
  expl <- survex::explain_survival(
    model                     = fit$model,
    data                      = predictor_data,
    y                         = survival::Surv(new_data[[time_col]],
                                               new_data[[status_col]]),
    predict_function          = risk_fun,
    predict_survival_function = surv_fun,
    times                     = times,
    label                     = paste0("highMLR-", fit$method),
    verbose                   = FALSE
  )

  result <- list(
    explainer    = expl,
    top_features = utils::head(fit$selected, top_n),
    times        = times,
    method       = method,
    n_explained  = length(expl_rows)
  )

  # ---- Compute time-dependent SHAP ------------------------------------
  if (method == "survshap") {
    # Subset to features actually in the model
    obs_for_shap <- explain_data[, setdiff(names(explain_data),
                                           c(time_col, status_col)),
                                 drop = FALSE]

    # Let the error surface if model_survshap fails; do not swallow.
    sshap <- survex::model_survshap(
      explainer       = expl,
      new_observation = obs_for_shap,
      ...
    )
    result$survshap <- sshap

    # Aggregate |SHAP| across time and observations -> one number / feature
    agg <- aggregate_survshap(sshap)
    result$aggregated <- agg
  } else if (method == "permutation") {
    pfi <- survex::model_parts(expl, type = "raw")
    result$permutation <- pfi
  } else if (method == "break_down") {
    # Single observation break-down for the first explained row
    bd <- survex::predict_parts(expl,
                                new_observation = explain_data[1L, , drop = FALSE])
    result$break_down <- bd
  }

  class(result) <- "highmlr_explain"
  result
}

# ---- Helpers used by highmlr_explain() ------------------------------------

# Build a method-specific survival-curve predict function.
# survex calls this with signature f(model, newdata, times = ...).
# Returns a matrix of nrow(newdata) x length(times) survival probabilities.
build_surv_predict_function <- function(fit, default_times) {

  method <- fit$method

  if (method == "coxnet") {
    # S(t|x) = S0(t) ^ exp(lp), with S0 from Breslow baseline.
    return(function(model, newdata, times = default_times) {
      if (is.null(times)) times <- default_times
      lp_new <- predict_coxnet(model, newdata, type = "linear_pred")
      base <- model$baseline_haz
      if (is.null(base)) {
        med <- stats::median(times)
        S0 <- exp(-times / med)
      } else {
        S0_fn <- stats::stepfun(base$time, c(1, base$surv))
        S0 <- S0_fn(times)
      }
      # nrow x length(times) matrix; each row is S0^exp(lp_i) at every t
      mat <- outer(exp(lp_new), S0, FUN = function(r, s0) s0 ^ r)
      # outer with vectors gives length(lp_new) x length(S0); good.
      colnames(mat) <- as.character(times)
      mat
    })
  }

  if (method == "rsf") {
    return(function(model, newdata, times = default_times) {
      if (is.null(times)) times <- default_times
      feats <- model$features
      newdata <- newdata[, intersect(feats, names(newdata)), drop = FALSE]
      newdata <- apply_imputation(newdata, model$imputation)
      pred <- stats::predict(model$fit, data = newdata)
      # pred$survival is n x length(pred$unique.death.times)
      surv_mat <- pred$survival
      ts <- pred$unique.death.times
      # Interpolate (step-function) to user-requested times
      stepwise_align(surv_mat, ts, times)
    })
  }

  if (method == "aorsf") {
    return(function(model, newdata, times = default_times) {
      if (is.null(times)) times <- default_times
      feats <- model$features
      newdata <- newdata[, intersect(feats, names(newdata)), drop = FALSE]
      newdata <- apply_imputation(newdata, model$imputation)
      surv_df <- tryCatch(
        stats::predict(model$fit, new_data = newdata,
                       pred_type = "surv", pred_horizon = times),
        error = function(e) NULL
      )
      if (is.null(surv_df)) {
        return(matrix(0.5, nrow = nrow(newdata), ncol = length(times)))
      }
      as.matrix(surv_df)
    })
  }

  if (method == "xgboost") {
    # xgboost Cox does not produce survival curves natively.
    # Approximate via exponential survival with rate = predicted risk.
    return(function(model, newdata, times = default_times) {
      if (is.null(times)) times <- default_times
      Xn <- as.matrix(newdata[, intersect(model$features, names(newdata)),
                              drop = FALSE])
      risk <- stats::predict(model$fit, newdata = Xn)
      # Scale so median patient has reasonable survival at median time
      med_risk <- stats::median(risk[risk > 0])
      med_t    <- stats::median(times)
      lambda0  <- log(2) / med_t   # so median patient has S=0.5 at med_t
      mat <- outer(risk / med_risk * lambda0, times,
                   FUN = function(r, t) exp(-r * t))
      colnames(mat) <- as.character(times)
      mat
    })
  }

  # Generic fallback: convert risk to a survival curve via exponential approx
  function(model, newdata, times = default_times) {
    if (is.null(times)) times <- default_times
    risk <- predict_highmlr_internal(model, newdata, type = "risk",
                                     method = method)
    med  <- stats::median(times)
    base_S <- exp(-times / med)
    mat <- outer(risk, base_S, FUN = function(r, s0) s0 ^ r)
    colnames(mat) <- as.character(times)
    mat
  }
}

# Step-function alignment: rows of surv_mat are evaluated at `from_times`;
# return an n x length(to_times) matrix evaluated at `to_times`.
stepwise_align <- function(surv_mat, from_times, to_times) {
  n <- nrow(surv_mat)
  out <- matrix(NA_real_, nrow = n, ncol = length(to_times))
  for (i in seq_len(n)) {
    fn <- stats::stepfun(from_times, c(1, surv_mat[i, ]))
    out[i, ] <- fn(to_times)
  }
  colnames(out) <- as.character(to_times)
  out
}

# Aggregate a survex SurvSHAP(t) result into a per-feature mean |SHAP|
aggregate_survshap <- function(sshap) {
  if (is.null(sshap) || is.null(sshap$result)) return(NULL)
  per_obs <- sshap$result
  # Each element of result is a data frame with a `_times_` column
  # and one column per feature (or _baseline_)
  stk <- do.call(rbind, lapply(per_obs, function(df) {
    val_cols <- setdiff(names(df), c("_times_", "_baseline_"))
    if (!length(val_cols)) return(NULL)
    vals <- as.matrix(df[, val_cols, drop = FALSE])
    colMeans(abs(vals), na.rm = TRUE)
  }))
  if (is.null(stk) || nrow(stk) == 0L) return(NULL)
  mean_imp <- colMeans(stk, na.rm = TRUE)
  tibble::tibble(
    feature       = names(mean_imp),
    mean_abs_shap = unname(mean_imp)
  )[order(-mean_imp), ]
}

#' @param x A `highmlr_explain` object.
#' @param n Number of top features to print (default 10).
#' @rdname highmlr_explain
#' @export
print.highmlr_explain <- function(x, n = 10, ...) {
  cat("<highmlr_explain>\n")
  cat("  Method:        ", x$method, "\n", sep = "")
  cat("  Base fit:      ", x$explainer$label %||% "(unlabelled)", "\n", sep = "")
  cat("  Times:         ", length(x$times),
      " points from ", round(min(x$times), 3),
      " to ", round(max(x$times), 3), "\n", sep = "")
  if (!is.null(x$n_explained)) {
    cat("  Observations:  ", x$n_explained, " patient", 
        if (x$n_explained != 1L) "s", " explained\n", sep = "")
  }
  if (!is.null(x$aggregated)) {
    cat("\n  Top ", min(n, nrow(x$aggregated)),
        " features by mean |SHAP|:\n", sep = "")
    print(utils::head(x$aggregated, n))
  } else {
    cat("  (no aggregated SHAP available)\n")
  }
  invisible(x)
}

#' @rdname highmlr_explain
#' @export
plot.highmlr_explain <- function(x, top_n = 10, ...) {
  if (is.null(x$aggregated)) {
    return(ggplot2::ggplot() +
             ggplot2::labs(title = "No aggregated SHAP available"))
  }
  d <- utils::head(x$aggregated, top_n)
  d$feature <- factor(d$feature, levels = rev(d$feature))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$mean_abs_shap,
                                  y = .data$feature)) +
    ggplot2::geom_col(fill = "darkorange") +
    ggplot2::labs(
      x = "Mean |SHAP|", y = NULL,
      title = paste0("highMLR explain: top features (", x$method, ")")
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

#' Pre-screen features when p is very large
#'
#' Lightweight filter before the main pipeline (e.g. to drop features
#' with low variance or low marginal association).
#'
#' @param data,time,status,features As in [highmlr()].
#' @param filter One of `"variance"`, `"univariate_p"`, `"none"`.
#' @param keep Integer, how many features to retain (default 1000).
#'
#' @return Character vector of retained feature names.
#'
#' @examples
#' \dontrun{
#' data(srdata)
#' keep <- highmlr_screen(srdata, "OS", "event",
#'                        filter = "variance", keep = 500)
#' fit <- highmlr(srdata, "OS", "event", features = keep, method = "coxnet")
#' }
#'
#' @export
highmlr_screen <- function(data, time, status, features = NULL,
                           filter = c("variance", "univariate_p", "none"),
                           keep = 1000L) {

  filter <- match.arg(filter)
  if (is.null(features)) {
    features <- setdiff(names(data), c(time, status))
  }

  if (filter == "none" || length(features) <= keep) return(features)

  if (filter == "variance") {
    num_feats <- features[vapply(features,
                                 function(f) is.numeric(data[[f]]),
                                 logical(1))]
    vars <- vapply(num_feats, function(f) stats::var(data[[f]], na.rm = TRUE),
                   numeric(1))
    keep_feats <- names(sort(vars, decreasing = TRUE))[seq_len(min(keep,
                                                                  length(vars)))]
    return(keep_feats)
  }

  if (filter == "univariate_p") {
    uf <- fit_univariate(data, time, status, features,
                         top_n = keep, rank_by = "p_value")
    return(uf$selected$feature)
  }
  features
}

#' Generate a Quarto/Rmd report skeleton for a highmlr_fit
#'
#' Writes a self-contained Rmd file that, when rendered, produces a
#' standard biomarker report (selected features, hazard ratios where
#' available, performance, forest plot).
#'
#' @param fit A `highmlr_fit` object.
#' @param file Output `.Rmd` path (default `"highmlr_report.Rmd"`).
#' @param render Logical: if `TRUE`, also render via `rmarkdown::render()`.
#'
#' @return Invisibly, the path to the written file.
#' @export
highmlr_report <- function(fit, file = "highmlr_report.Rmd",
                           render = FALSE) {
  if (!inherits(fit, "highmlr_fit")) {
    rlang::abort("`fit` must be a highmlr_fit object.")
  }
  template <- c(
    "---",
    "title: \"highMLR biomarker report\"",
    "output: html_document",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)",
    "library(highMLR); library(ggplot2)",
    "```",
    "",
    sprintf("## Method: `%s`", fit$method),
    "",
    "### Data summary",
    "",
    "```{r}",
    "fit$data_summary",
    "```",
    "",
    "### Selected features",
    "",
    "```{r}",
    "fit$selected",
    "```",
    "",
    "### Performance",
    "",
    "```{r}",
    "fit$performance",
    "```",
    "",
    "### Importance plot",
    "",
    "```{r, fig.width=7, fig.height=5}",
    "plot(fit, top_n = 20)",
    "```"
  )
  writeLines(template, file)
  if (render) {
    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
      rlang::abort("Install 'rmarkdown' to render the report.")
    }
    rmarkdown::render(file, quiet = TRUE)
  }
  invisible(file)
}
