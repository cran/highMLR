#' Constructor for highmlr_fit objects
#'
#' Internal constructor. Not exported. Use [highmlr()] to create fits.
#'
#' @param selected A `tibble` of selected features with columns
#'   `feature`, `importance`, and method-specific extras (e.g. `coef`,
#'   `hazard_ratio`, `selection_freq`).
#' @param performance A named list of out-of-sample performance metrics
#'   (e.g. `c_index`, `ibs`).
#' @param model The fitted underlying model object (parsnip fit,
#'   glmnet object, or list of these for stability selection).
#' @param method Character scalar, the method used.
#' @param call The matched call.
#' @param data_summary A small list summarising the input data
#'   (n, p, events, censoring rate).
#' @param meta Optional list for method-specific metadata.
#'
#' @return An object of class `highmlr_fit`.
#' @keywords internal
new_highmlr_fit <- function(selected,
                            performance = NULL,
                            model = NULL,
                            method = NA_character_,
                            call = NULL,
                            data_summary = list(),
                            meta = list()) {

  stopifnot(tibble::is_tibble(selected) || is.data.frame(selected))
  selected <- tibble::as_tibble(selected)

  structure(
    list(
      selected     = selected,
      performance  = performance,
      model        = model,
      method       = method,
      call         = call,
      data_summary = data_summary,
      meta         = meta
    ),
    class = "highmlr_fit"
  )
}

#' Print method for highmlr_fit
#'
#' @param x A `highmlr_fit` object.
#' @param n Number of top features to display (default 10).
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.highmlr_fit <- function(x, n = 10, ...) {
  cat("<highmlr_fit>\n")
  cat("  Method:    ", x$method, "\n", sep = "")
  ds <- x$data_summary
  if (length(ds)) {
    cat("  Data:      n = ", ds$n %||% NA,
        ", p = ", ds$p %||% NA,
        ", events = ", ds$events %||% NA,
        " (", round(100 * (ds$event_rate %||% NA), 1), "% event rate)\n",
        sep = "")
  }
  k <- nrow(x$selected)
  cat("  Selected:  ", k, " feature", if (k != 1L) "s", "\n", sep = "")
  if (!is.null(x$performance) && length(x$performance)) {
    perf_str <- paste(
      names(x$performance),
      vapply(x$performance, function(v) sprintf("%.4f", v), character(1)),
      sep = " = ", collapse = ", "
    )
    cat("  Performance: ", perf_str, "\n", sep = "")
  }
  if (k > 0) {
    cat("\n  Top ", min(n, k), " features:\n", sep = "")
    print(utils::head(x$selected, n))
  }
  invisible(x)
}

#' Summary method for highmlr_fit
#'
#' @param object A `highmlr_fit` object.
#' @param ... Unused.
#'
#' @return A list with the full selected feature table and performance.
#' @export
summary.highmlr_fit <- function(object, ...) {
  out <- list(
    method       = object$method,
    data_summary = object$data_summary,
    n_selected   = nrow(object$selected),
    selected     = object$selected,
    performance  = object$performance,
    call         = object$call
  )
  class(out) <- "summary.highmlr_fit"
  out
}

#' Coefficients from a highmlr_fit
#'
#' @param object A `highmlr_fit` object.
#' @param ... Unused.
#'
#' @return A named numeric vector of coefficients (where defined) or
#'   importance scores otherwise.
#' @export
coef.highmlr_fit <- function(object, ...) {
  s <- object$selected
  val_col <- if ("coef" %in% names(s)) "coef" else "importance"
  stats::setNames(s[[val_col]], s$feature)
}

#' Forest / importance plot for a highmlr_fit
#'
#' @param x A `highmlr_fit` object.
#' @param top_n Number of top features to plot (default 20).
#' @param ... Unused.
#'
#' @return A `ggplot` object.
#' @export
plot.highmlr_fit <- function(x, top_n = 20, ...) {
  s <- x$selected
  if (!nrow(s)) {
    return(ggplot2::ggplot() + ggplot2::labs(title = "No features selected"))
  }
  val_col <- if ("coef" %in% names(s)) "coef" else "importance"
  s <- utils::head(s[order(-abs(s[[val_col]])), ], top_n)
  s$feature <- factor(s$feature, levels = rev(s$feature))

  ggplot2::ggplot(s, ggplot2::aes(x = .data[[val_col]], y = .data$feature)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey40") +
    ggplot2::labs(
      x = val_col,
      y = NULL,
      title = paste0("highMLR: top ", nrow(s), " features (",
                     x$method, ")")
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

#' Predict from a highmlr_fit
#'
#' @param object A `highmlr_fit` object.
#' @param new_data A data frame containing the features used in fitting.
#' @param type One of `"linear_pred"`, `"survival"`, or `"risk"`.
#'   Availability depends on the underlying model.
#' @param ... Passed to the underlying model's predict method.
#'
#' @return Predicted values (vector or tibble depending on `type`).
#' @export
predict.highmlr_fit <- function(object, new_data,
                                type = c("linear_pred", "survival", "risk"),
                                ...) {
  type <- match.arg(type)
  if (is.null(object$model)) {
    rlang::abort("This highmlr_fit has no stored model; cannot predict.")
  }
  predict_highmlr_internal(object$model, new_data, type = type,
                           method = object$method, ...)
}

# Internal dispatcher; method-specific predict logic lives here.
predict_highmlr_internal <- function(model, new_data, type, method, ...) {
  switch(method,
    coxnet    = predict_coxnet(model, new_data, type = type, ...),
    rsf       = predict_parsnip_surv(model, new_data, type = type, ...),
    aorsf     = predict_parsnip_surv(model, new_data, type = type, ...),
    xgboost   = predict_parsnip_surv(model, new_data, type = type, ...),
    rlang::abort(paste0("predict() not implemented for method '", method, "'."))
  )
}

# Null-coalescing operator used internally
`%||%` <- function(a, b) if (is.null(a)) b else a
