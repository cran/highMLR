#' Conformal prediction intervals for survival times
#'
#' Computes calibrated lower bounds on survival time for each new
#' subject using a split-conformal procedure with inverse probability
#' of censoring weights (Candes, Lei and Ren, 2023). The returned
#' lower bound satisfies a marginal coverage guarantee approximately
#' equal to one minus alpha under standard conformal assumptions and
#' a consistent censoring model.
#'
#' @param fit A highmlr_fit object whose predict() method returns a
#'   linear predictor or risk score.
#' @param new_data Data frame on which to compute prediction intervals.
#' @param calibration_data Data frame on which to compute conformity
#'   scores. If NULL, a random calibration_split fraction of new_data
#'   is held out for calibration and the rest is used as the test set
#'   (split-conformal).
#' @param alpha Miscoverage level; default 0.1 (so 90 percent coverage).
#' @param calibration_split Fraction of new_data to use for calibration
#'   when calibration_data is NULL. Default 0.3.
#' @param time Name of the survival time column in calibration data.
#'   Defaults to the column used in fit.
#' @param status Name of the event column in calibration data.
#' @param seed Optional integer seed for the split.
#'
#' @return An object of class highmlr_conformal containing per-subject
#'   point predictions and lower confidence bounds for survival time.
#'
#' @examples
#' \dontrun{
#' fit  <- highmlr(d_train, "OS", "Death", method = "coxnet")
#' intv <- highmlr_conformal(fit, new_data = d_test, alpha = 0.1)
#' print(intv)
#' plot(intv)
#' }
#'
#' @export
highmlr_conformal <- function(fit,
                              new_data,
                              calibration_data  = NULL,
                              alpha             = 0.1,
                              calibration_split = 0.3,
                              time              = NULL,
                              status            = NULL,
                              seed              = NULL) {

  if (!inherits(fit, "highmlr_fit")) {
    rlang::abort("`fit` must be a highmlr_fit object.")
  }
  if (!is.data.frame(new_data)) {
    rlang::abort("`new_data` must be a data frame.")
  }

  if (is.null(time)) {
    time <- as.character(fit$call$time)
    if (!length(time)) time <- "OS"
  }
  if (is.null(status)) {
    status <- as.character(fit$call$status)
    if (!length(status)) status <- "Death"
  }

  if (!is.null(seed)) set.seed(seed)

  if (is.null(calibration_data)) {
    n <- nrow(new_data)
    cal_idx <- sample.int(n, ceiling(calibration_split * n))
    calibration_data <- new_data[cal_idx, , drop = FALSE]
    test_data        <- new_data[-cal_idx, , drop = FALSE]
  } else {
    test_data <- new_data
  }

  if (!all(c(time, status) %in% names(calibration_data))) {
    rlang::abort(sprintf(
      "Calibration data missing '%s' and/or '%s'.", time, status
    ))
  }

  lp_cal <- stats::predict(fit, new_data = calibration_data,
                           type = "linear_pred")

  t_cal <- calibration_data[[time]]
  d_cal <- as.numeric(calibration_data[[status]])

  km_cens <- survival::survfit(
    survival::Surv(t_cal, 1 - d_cal) ~ 1
  )
  G_at <- stats::stepfun(km_cens$time, c(1, km_cens$surv))
  w <- ifelse(d_cal == 1L, 1 / pmax(G_at(t_cal), 0.05), 0)

  log_t_cal <- log(pmax(t_cal, 1e-8))
  scores <- log_t_cal + lp_cal

  abs_scores <- abs(scores)
  q_alpha <- weighted_quantile(abs_scores[d_cal == 1L],
                               w[d_cal == 1L],
                               probs = 1 - alpha)

  lp_test <- stats::predict(fit, new_data = test_data,
                            type = "linear_pred")

  log_t_hat <- -lp_test
  log_t_lo  <- log_t_hat - q_alpha
  t_hat     <- exp(log_t_hat)
  t_lo      <- exp(log_t_lo)

  intervals <- tibble::tibble(
    row_id    = seq_len(nrow(test_data)),
    t_hat     = t_hat,
    t_lcb     = t_lo,
    log_t_hat = log_t_hat,
    log_t_lcb = log_t_lo
  )

  out <- list(
    intervals            = intervals,
    test_data            = test_data,
    calibration_data     = calibration_data,
    alpha                = alpha,
    coverage_target      = 1 - alpha,
    q_alpha              = q_alpha,
    fit_method           = fit$method,
    n_calibration        = nrow(calibration_data),
    n_test               = nrow(test_data),
    n_calibration_events = sum(d_cal == 1L),
    call                 = match.call()
  )
  class(out) <- "highmlr_conformal"
  out
}

# Weighted quantile helper used by highmlr_conformal()
weighted_quantile <- function(x, w, probs) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (!length(x)) return(NA_real_)
  ord <- order(x)
  x <- x[ord]; w <- w[ord]
  cw <- cumsum(w) / sum(w)
  idx <- which(cw >= probs)[1L]
  if (is.na(idx)) return(x[length(x)])
  x[idx]
}

#' Print method for highmlr_conformal objects
#'
#' @param x A highmlr_conformal object.
#' @param n Number of rows to display in the preview table (default 10).
#' @param ... Unused.
#'
#' @return Invisibly returns x.
#' @export
print.highmlr_conformal <- function(x, n = 10, ...) {
  cat("<highmlr_conformal>\n")
  cat("  Base method:    ", x$fit_method, "\n", sep = "")
  cat("  Coverage target:", round(100 * x$coverage_target, 1), "%\n",
      sep = " ")
  cat("  Calibration n:  ", x$n_calibration, " (",
      x$n_calibration_events, " events)\n", sep = "")
  cat("  Test n:         ", x$n_test, "\n", sep = "")
  cat("  Conformity q:   ", round(x$q_alpha, 3),
      " (log-time units)\n", sep = "")
  cat("\n  First ", min(n, nrow(x$intervals)),
      " predictions (point estimate and lower bound on T):\n", sep = "")
  print(utils::head(x$intervals[, c("row_id", "t_hat", "t_lcb")], n))
  invisible(x)
}

#' Plot method for highmlr_conformal objects
#'
#' @param x A highmlr_conformal object.
#' @param ... Unused.
#'
#' @return A ggplot object.
#' @export
plot.highmlr_conformal <- function(x, ...) {
  d <- x$intervals
  d$row_id <- factor(d$row_id, levels = d$row_id[order(d$t_hat)])
  ggplot2::ggplot(d, ggplot2::aes(x = .data$row_id)) +
    ggplot2::geom_point(ggplot2::aes(y = .data$t_hat),
                        colour = "steelblue") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$t_lcb, ymax = .data$t_hat),
      colour = "grey60", width = 0
    ) +
    ggplot2::labs(
      title = paste0("Conformal lower bounds (",
                     round(100 * x$coverage_target, 0), "% coverage)"),
      x = "Patient (sorted)",
      y = "Survival time"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
}
