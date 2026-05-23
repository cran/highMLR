# Tests for highmlr() across all methods on a small synthetic dataset.

make_synth <- function(n = 80, p = 30, seed = 1L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("g", seq_len(p))
  # First 5 features are truly informative
  beta <- c(rep(0.7, 5), rep(0, p - 5))
  lp <- as.numeric(X %*% beta)
  t  <- stats::rexp(n, rate = exp(lp - mean(lp)))
  c  <- stats::rexp(n, rate = 0.1)
  time   <- pmin(t, c)
  status <- as.integer(t <= c)
  data.frame(OS = time, Death = status, X)
}

test_that("coxnet runs and selects informative features", {
  d <- make_synth()
  skip_if_not_installed("glmnet")
  fit <- highmlr(d, time = "OS", status = "Death",
                 method = "coxnet", resampling = "none")
  expect_s3_class(fit, "highmlr_fit")
  expect_true(nrow(fit$selected) >= 1L)
  expect_true(all(c("feature", "coef", "hazard_ratio", "importance") %in%
                  names(fit$selected)))
})

test_that("rsf runs and produces importance", {
  skip_on_cran()
  d <- make_synth()
  skip_if_not_installed("ranger")
  fit <- highmlr(d, time = "OS", status = "Death",
                 method = "rsf", top_n = 10)
  expect_s3_class(fit, "highmlr_fit")
  expect_lte(nrow(fit$selected), 10L)
})

test_that("aorsf runs", {
  skip_on_cran()
  d <- make_synth()
  skip_if_not_installed("aorsf")
  fit <- highmlr(d, time = "OS", status = "Death",
                 method = "aorsf", top_n = 10)
  expect_s3_class(fit, "highmlr_fit")
})

test_that("xgboost runs", {
  skip_on_cran()
  d <- make_synth()
  skip_if_not_installed("xgboost")
  fit <- highmlr(d, time = "OS", status = "Death",
                 method = "xgboost", top_n = 10, folds = 3)
  expect_s3_class(fit, "highmlr_fit")
})

test_that("univariate runs and returns p-values", {
  d <- make_synth()
  fit <- highmlr(d, time = "OS", status = "Death",
                 method = "univariate", top_n = 10)
  expect_s3_class(fit, "highmlr_fit")
  expect_true("p_value" %in% names(fit$selected))
})

test_that("stability runs (small B for speed)", {
  d <- make_synth()
  skip_if_not_installed("stabs")
  fit <- highmlr(d, time = "OS", status = "Death",
                 method = "stability", B = 20, cutoff = 0.5, PFER = 5)
  expect_s3_class(fit, "highmlr_fit")
})

test_that("compare runs across methods", {
  d <- make_synth()
  cmp <- highmlr_compare(d, "OS", "Death",
                         methods = c("coxnet", "univariate"),
                         resampling = "none")
  expect_named(cmp, c("fits", "summary"))
  expect_s3_class(cmp$summary, "tbl_df")
})

test_that("input validation catches bad arguments", {
  d <- make_synth()
  expect_error(highmlr(d, time = "missing", status = "Death"))
  expect_error(highmlr(d, time = "OS", status = "missing"))
  expect_error(highmlr("not_a_df", time = "OS", status = "Death"))
})
