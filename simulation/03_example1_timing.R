#!/usr/bin/env Rscript

# Example 1 timing. Default runs the requested grid once; --smoke is tiny.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", ".UTF-8"))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
script_dir <- dirname(normalizePath(script, winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "workflow_utils.R"), encoding = "UTF-8")
v2 <- v2_root(script_dir)
output <- file.path(v2, "results", "example1_timing", paste0(if ("--smoke" %in% commandArgs(TRUE)) "smoke_" else "run_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
writeLines("status=running", file.path(output, "run_status.txt"))
source(file.path(script_dir, "simulation_functions.R"), encoding = "UTF-8")
args <- commandArgs(trailingOnly = TRUE)
if (!all(args %in% "--smoke")) stop("Only --smoke is supported.")
smoke <- "--smoke" %in% args
grid <- if (smoke) {
  data.frame(n = 24L, p = 12L)
} else {
  expand.grid(n = c(200L, 500L), p = c(200L, 500L))
}
ratios <- c(0.3, 0.5, 0.7, 1)
kappa <- 0.5
tau <- 1e-10
base_seed <- 202611240L
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

check_simulation_packages()
stopifnot("CLARABEL" %in% CVXR::installed_solvers())
time_call <- function(expr) {
  start <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, seconds = unname(proc.time()[["elapsed"]] - start))
}

# Compute each unordered pair touching a landmark exactly once.
# Rows of W and the fitted scores follow the original variable order.
compute_landmark_columns <- function(X, landmarks) {
  p <- ncol(X)
  s <- length(landmarks)
  rest <- setdiff(seq_len(p), landmarks)
  A <- diag(s)
  B <- matrix(0, nrow = length(rest), ncol = s)
  for (i in seq_len(s)) {
    if (i < s) {
      j <- seq.int(i + 1L, s)
      values <- vapply(landmarks[j], function(k) {
        bcor_scalar(X[, landmarks[i]], X[, k])
      }, numeric(1L))
      A[j, i] <- A[i, j] <- values
    }
    B[, i] <- vapply(rest, function(k) {
      bcor_scalar(X[, landmarks[i]], X[, k])
    }, numeric(1L))
  }
  W <- matrix(0, p, s)
  W[landmarks, ] <- A
  W[rest, ] <- B
  list(W = W, A = A)
}

factorize_landmark_columns <- function(blocks) {
  spectral <- positive_spectral_inverse(blocks$A, tau)
  L <- sweep(blocks$W %*% spectral$vectors, 2L, sqrt(spectral$values), "/")
  diagonal <- rowSums(L^2)
  stopifnot(max(diagonal) <= 1 + 1e-8)
  L <- sweep(L, 1L, sqrt(pmax(diagonal, 1)), "/")
  d <- pmax(1 - rowSums(L^2), 0)
  list(L = L, diagonal_correction = d,
       qbar_all = (sum(colSums(L)^2) + sum(d)) / nrow(L)^2,
       rank = spectral$effective_rank)
}

# Load and initialize the common code paths outside the measured fits.
warm <- simulate_mixed_hermite_example(1L, 12L, 6L, kappa, base_seed)
F_warm <- compute_relevance(warm$X, warm$Y)
Q_warm <- compute_full_bcor_matrix(warm$X)
a_warm <- factorize_landmark_columns(compute_landmark_columns(warm$X, c(2L, 4L)))
invisible(solve_exact_bcorqpfs(F_warm, Q_warm, 0.5, solver = "CLARABEL"))
invisible(solve_nystrom_bcorqpfs(F_warm, a_warm, 0.5, solver = "CLARABEL"))

tag <- if (smoke) "smoke" else "full_grid"
csv_path <- file.path(output, paste0("example1_timing_", tag, ".csv"))
session_path <- file.path(output, paste0("example1_session_", tag, ".txt"))
rows <- list()
for (i in seq_len(nrow(grid))) {
  n <- grid$n[i]
  p <- grid$p[i]
  data_seed <- base_seed + i
  landmark_seed <- base_seed + 100L + i
  data <- simulate_mixed_hermite_example(1L, n, p, kappa, data_seed)
  relevance <- time_call(compute_relevance(data$X, data$Y))
  F <- relevance$value
  smoke_blocks <- list()
  for (r in ratios) {
    s <- as.integer(floor(r * p))
    method <- if (r < 1) "nystrom" else "exact"
    message(sprintf("Example 1, n=%d, p=%d, r=%.1f, %s", n, p, r, method))
    if (r < 1) {
      sampling <- time_call({
        set.seed(landmark_seed)
        # Common permutation gives nested uniform landmark sets within a grid cell.
        sample.int(p)[seq_len(s)]
      })
      landmarks <- sampling$value
      blocks <- time_call(compute_landmark_columns(data$X, landmarks))
      factor <- time_call(factorize_landmark_columns(blocks$value))
      approximation <- factor$value
      alpha <- choose_balanced_alpha(approximation$qbar_all, F)
      qp <- time_call(solve_nystrom_bcorqpfs(F, approximation, alpha, solver = "CLARABEL"))
      matrix_seconds <- blocks$seconds
      factor_seconds <- factor$seconds
      sampling_seconds <- sampling$seconds
      retained_rank <- approximation$rank
      if (smoke) smoke_blocks[[as.character(r)]] <- list(landmarks, blocks$value)
    } else {
      full <- time_call(compute_full_bcor_matrix(data$X))
      alpha <- choose_balanced_alpha(mean(full$value), F)
      qp <- time_call(solve_exact_bcorqpfs(F, full$value, alpha, solver = "CLARABEL"))
      matrix_seconds <- full$seconds
      factor_seconds <- sampling_seconds <- 0
      retained_rank <- NA_integer_
      if (smoke) {
        for (check in smoke_blocks) {
          stopifnot(max(abs(check[[2L]]$W - full$value[, check[[1L]], drop = FALSE])) < 1e-12)
        }
      }
    }
    post_relevance_seconds <- sampling_seconds + matrix_seconds + factor_seconds + qp$seconds
    rows[[length(rows) + 1L]] <- data.frame(
      example = 1L, n, p, kappa, r, s, method, retained_rank, tau, alpha,
      data_seed, landmark_seed = if (r < 1) landmark_seed else NA_integer_,
      pairs_evaluated = choose(s, 2) + s * (p - s),
      relevance_seconds = relevance$seconds, sampling_seconds, matrix_seconds,
      factor_seconds, qp_seconds = qp$seconds, post_relevance_seconds,
      total_seconds = relevance$seconds + post_relevance_seconds,
      solver = qp$value$solver, status = qp$value$status
    )
    write.csv(do.call(rbind, rows), csv_path, row.names = FALSE)
  }
}
capture.output(sessionInfo(), file = session_path)
writeLines("status=completed", file.path(output, "run_status.txt"))
print(do.call(rbind, rows)[, c("n", "p", "r", "pairs_evaluated", "matrix_seconds",
                              "factor_seconds", "qp_seconds", "total_seconds", "status")])
if (smoke) message("Smoke check passed: all four fits ran and on-demand columns match the full matrix.")
