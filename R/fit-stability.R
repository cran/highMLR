# Stability selection backend (Meinshausen & Buehlmann 2010; Shah & Samworth 2013)
# Primary path: stabs::stabsel with a properly-shaped fit function.
# Fallback path: manual subsampling-based stability selection.
# The fallback ensures users always get a result even if the stabs interface
# rejects the fit function (which has happened in some versions).

fit_stability <- function(data, time, status, features,
                          B = 100L,
                          PFER = 1,
                          cutoff = 0.75,
                          assumption = "unimodal",
                          parallel = FALSE,
                          alpha = 1,
                          q = NULL,
                          ...) {

  df <- data[, c(time, status, features), drop = FALSE]
  df <- impute_simple(df, features)

  X <- as.matrix(df[, features, drop = FALSE])
  y <- survival::Surv(df[[time]], df[[status]])

  sds  <- apply(X, 2, stats::sd, na.rm = TRUE)
  keep <- which(sds > 0)
  X    <- X[, keep, drop = FALSE]
  feature_names <- colnames(X)
  p <- ncol(X)

  if (is.null(q)) {
    # Default q from PFER/cutoff per Meinshausen-Buhlmann formula
    q <- max(2L, ceiling(sqrt(0.8 * PFER * p * (2 * cutoff - 1))))
    q <- min(q, p)
  }

  # ---- Attempt 1: stabs::stabsel ----------------------------------------
  stab <- try_stabsel(X, y, alpha, cutoff, PFER, B, assumption,
                      parallel, ...)

  if (is.null(stab)) {
    # ---- Fallback: manual subsampling stability ------------------------
    message("Using manual stability selection (stabs unavailable for this fit).")
    probs <- manual_stability(X, y, alpha = alpha, B = B, q = q,
                              parallel = parallel)
    threshold_val <- cutoff
  } else {
    probs <- stab$max
    threshold_val <- stab$cutoff
  }

  sel_idx <- which(probs >= cutoff)
  sel_probs <- sort(probs[sel_idx], decreasing = TRUE)

  selected <- tibble::tibble(
    feature        = names(sel_probs),
    selection_freq = unname(sel_probs),
    importance     = unname(sel_probs)
  )

  performance <- list(
    n_selected = nrow(selected),
    B          = B,
    cutoff     = cutoff,
    PFER       = PFER,
    q          = q,
    threshold  = threshold_val
  )

  new_highmlr_fit(
    selected    = selected,
    performance = performance,
    model       = list(stab = stab, alpha = alpha,
                       probabilities = probs,
                       feature_names = feature_names),
    meta        = list(B = B, cutoff = cutoff, PFER = PFER,
                       assumption = assumption,
                       used_fallback = is.null(stab))
  )
}

# ---- stabs path -----------------------------------------------------------
try_stabsel <- function(X, y, alpha, cutoff, PFER, B, assumption,
                        parallel, ...) {

  fit_fun <- function(x, y, q, ...) {
    p_local <- ncol(x)
    feat_nms <- colnames(x); if (is.null(feat_nms)) feat_nms <- paste0("V", seq_len(p_local))
    g <- glmnet::glmnet(x, y, family = "cox", alpha = alpha,
                        pmax = q, standardize = TRUE)
    path_mat <- as.matrix(g$beta != 0)
    storage.mode(path_mat) <- "logical"
    rownames(path_mat) <- feat_nms
    colnames(path_mat) <- paste0("s", seq_len(ncol(path_mat)) - 1L)
    sel_vec <- if (ncol(path_mat) > 0L) path_mat[, ncol(path_mat)]
               else stats::setNames(rep(FALSE, p_local), feat_nms)
    list(selected = sel_vec, path = path_mat)
  }

  tryCatch(
    stabs::stabsel(
      x = X, y = y, fitfun = fit_fun,
      cutoff = cutoff, PFER = PFER, B = B,
      assumption = assumption, sampling.type = "SS",
      args.fitfun = list(),
      mc.cores = if (parallel) max(1L, future::availableCores() - 1L) else 1L,
      verbose = FALSE, ...
    ),
    error = function(e) NULL
  )
}

# ---- Manual fallback ------------------------------------------------------
manual_stability <- function(X, y, alpha = 1, B = 100L, q = 10L,
                             parallel = FALSE) {
  n <- nrow(X); p <- ncol(X)
  feat_nms <- colnames(X)
  half <- floor(n / 2)

  one_replicate <- function(b) {
    idx <- sample.int(n, half)
    Xs <- X[idx, , drop = FALSE]; ys <- y[idx, ]
    g <- tryCatch(
      glmnet::glmnet(Xs, ys, family = "cox", alpha = alpha,
                     pmax = q, standardize = TRUE),
      error = function(e) NULL
    )
    if (is.null(g) || ncol(g$beta) == 0L) return(rep(FALSE, p))
    sel_step <- as.matrix(g$beta[, ncol(g$beta)] != 0)
    as.logical(sel_step)
  }

  reps <- if (parallel) {
    do.call(cbind, future.apply::future_lapply(seq_len(B), one_replicate,
                                               future.seed = TRUE))
  } else {
    do.call(cbind, lapply(seq_len(B), one_replicate))
  }
  probs <- rowMeans(reps)
  names(probs) <- feat_nms
  probs
}
