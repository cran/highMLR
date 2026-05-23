# Tests for v0.3.0 additions: pseudo, finegray, causal, conformal,
# upgraded explain.

make_synth <- function(n = 100, p = 20, seed = 1L,
                       competing = FALSE) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("g", seq_len(p))
  beta <- c(rep(0.7, 5), rep(0, p - 5))
  lp <- as.numeric(X %*% beta)
  t  <- stats::rexp(n, rate = exp(lp - mean(lp)))
  c  <- stats::rexp(n, rate = 0.1)
  time   <- pmin(t, c)
  if (competing) {
    # Three-state: 0 = censored, 1 = cause of interest, 2 = competing
    status <- ifelse(t <= c,
                     sample(c(1L, 2L), n, replace = TRUE, prob = c(0.7, 0.3)),
                     0L)
  } else {
    status <- as.integer(t <= c)
  }
  data.frame(OS = time, Death = status, X)
}

test_that("pseudo method with ranger runs", {
  skip_on_cran()
  skip_if_not_installed("prodlim")
  skip_if_not_installed("ranger")
  d <- make_synth()
  fit <- highmlr(d, time = "OS", status = "Death",
                 method = "pseudo", engine = "ranger", top_n = 10)
  expect_s3_class(fit, "highmlr_fit")
  expect_true(nrow(fit$selected) >= 1L)
  expect_true("importance" %in% names(fit$selected))
})

test_that("pseudo method with glmnet runs", {
  skip_if_not_installed("prodlim")
  skip_if_not_installed("glmnet")
  d <- make_synth()
  fit <- highmlr(d, time = "OS", status = "Death",
                 method = "pseudo", engine = "glmnet", top_n = 10)
  expect_s3_class(fit, "highmlr_fit")
})

test_that("finegray method on competing risks data runs", {
  skip_if_not_installed("cmprsk")
  d <- make_synth(competing = TRUE)
  fit <- highmlr(d, time = "OS", status = "Death",
                 method = "finegray", top_n = 10)
  expect_s3_class(fit, "highmlr_fit")
  expect_true("subdistribution_hr" %in% names(fit$selected))
})

test_that("highmlr_causal runs", {
  skip_on_cran()
  skip_if_not_installed("grf")
  set.seed(2)
  n <- 200; p <- 8
  X <- matrix(stats::rnorm(n*p), n, p); colnames(X) <- paste0("V", 1:p)
  W <- stats::rbinom(n, 1, 0.5)
  t <- stats::rexp(n, rate = exp(0.3*W + 0.5*X[,1]*W))
  c <- stats::rexp(n, rate = 0.05)
  d <- data.frame(OS = pmin(t,c), Death = as.integer(t<=c),
                  arm = W, X)
  cf <- highmlr_causal(d, "OS", "Death", treatment = "arm",
                       covariates = paste0("V", 1:p),
                       num.trees = 500L)
  expect_s3_class(cf, "highmlr_causal")
  expect_true("cate" %in% names(cf$patient_effects))
  expect_equal(nrow(cf$patient_effects), n)
})

test_that("highmlr_conformal returns valid intervals", {
  skip_if_not_installed("glmnet")
  d <- make_synth(n = 200)
  train_idx <- 1:120
  fit <- highmlr(d[train_idx, ], "OS", "Death", method = "coxnet",
                 resampling = "none")
  intv <- highmlr_conformal(fit, new_data = d[-train_idx, ],
                            alpha = 0.1)
  expect_s3_class(intv, "highmlr_conformal")
  expect_true(all(intv$intervals$t_lcb <= intv$intervals$t_hat))
})
