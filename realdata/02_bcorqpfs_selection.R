#!/usr/bin/env Rscript

# Run bcorQPFS on saved BCor-screened ADNI predictors.

suppressPackageStartupMessages({
  library(Ball)
  library(CVXR)
})

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_PATH <- if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else
  tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
SCRIPT_DIR <- if (!is.null(SCRIPT_PATH) && nzchar(SCRIPT_PATH)) {
  dirname(SCRIPT_PATH)
} else {
  getwd()
}
PROJECT_ROOT <- file.path(SCRIPT_DIR, "..", "..")
DATA_DIR <- file.path(PROJECT_ROOT, "code", "data")
RESULT_DIR <- file.path(DATA_DIR, "result", "BcorQPFS_0729")
CACHE_DIR <- file.path(RESULT_DIR, "cache")

# Set the selection, Nystrom, numerical, and parallel parameters.
M_SELECT <- 30L
LANDMARK_PROPORTION <- 0.30
SPECTRAL_REL_TOL <- 1e-10
SEED <- 20260729L
N_CORES <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
CVXR_SOLVER <- NULL
SIDES <- c("left", "right")

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

log_message <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  paste0(..., collapse = "")))
}

load_rdata_object <- function(path, preferred = NULL) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = env)
  if (!is.null(preferred) && preferred %in% object_names) return(env[[preferred]])
  if ("selected_X" %in% object_names) return(env[["selected_X"]])
  if (length(object_names) == 1L) return(env[[object_names]])
  sizes <- vapply(object_names, function(nm) object.size(env[[nm]]), numeric(1))
  env[[object_names[[which.max(sizes)]]]]
}

as_numeric_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x
}

find_screened_file <- function(side) {
  candidates <- c(
    file.path(DATA_DIR, sprintf("ADNI_%s_screening_selected.RData", side)),
    file.path(DATA_DIR, sprintf("ADNI_COMBINE_%s_selected.RData", side))
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("No saved screening result found for ", side, ". Expected one of: ",
         paste(basename(candidates), collapse = ", "), call. = FALSE)
  }
  hit[[1L]]
}

align_X_Y <- function(X, Y) {
  X <- as_numeric_matrix(X)
  Y <- as_numeric_matrix(Y)
  x_ids <- rownames(X)
  y_ids <- rownames(Y)

  if (!is.null(x_ids) && !is.null(y_ids)) {
    common <- intersect(x_ids, y_ids)
    if (!length(common)) stop("X and Y have no common row names.", call. = FALSE)
    X <- X[common, , drop = FALSE]
    Y <- Y[common, , drop = FALSE]
  } else if (nrow(X) != nrow(Y)) {
    stop("X and Y cannot be aligned: row names are absent and row counts differ.",
         call. = FALSE)
  }

  keep <- apply(X, 2L, function(z) {
    all(is.finite(z)) && length(unique(z)) > 1L
  })
  if (!all(keep)) {
    log_message("Removed ", sum(!keep), " non-finite or constant screened predictors.")
    X <- X[, keep, drop = FALSE]
  }
  if (!ncol(X)) stop("No usable predictors remain.", call. = FALSE)
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))
  list(X = X, Y = Y)
}

make_cluster <- function(n_cores) {
  if (n_cores <= 1L) return(NULL)
  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(Ball)))
  cl
}

parallel_map <- function(X, FUN, cl) {
  if (is.null(cl)) lapply(X, FUN) else parallel::parLapplyLB(cl, X, FUN)
}

compute_relevance <- function(X, Y, cl = NULL) {
  DY <- stats::dist(Y)
  if (!is.null(cl)) {
    parallel::clusterExport(cl, c("X", "DY"), envir = environment())
  }
  values <- parallel_map(seq_len(ncol(X)), function(j) {
    as.numeric(Ball::bcor(stats::dist(X[, j, drop = FALSE]), DY,
                          distance = TRUE))
  }, cl)
  F <- unlist(values, use.names = FALSE)
  names(F) <- colnames(X)
  F
}

# Compute landmark and cross blocks of the predictor dependence matrix.
compute_landmark_blocks <- function(X, landmark_idx, cl = NULL) {
  p <- ncol(X)
  s <- length(landmark_idx)
  other_idx <- setdiff(seq_len(p), landmark_idx)
  perm <- c(landmark_idx, other_idx)
  Xp <- X[, perm, drop = FALSE]

  if (!is.null(cl)) {
    parallel::clusterExport(cl, c("Xp", "s", "p"), envir = environment())
  }
  rows <- parallel_map(seq_len(s), function(i) {
    vals <- numeric(p)
    vals[i] <- 1
    if (i < p) {
      for (j in seq.int(i + 1L, p)) {
        vals[j] <- as.numeric(Ball::bcor(Xp[, i], Xp[, j]))
      }
    }
    vals
  }, cl)

  A <- diag(1, s)
  B <- matrix(0, nrow = s, ncol = p - s)
  for (i in seq_len(s)) {
    if (i < s) {
      A[i, (i + 1L):s] <- rows[[i]][(i + 1L):s]
      A[(i + 1L):s, i] <- rows[[i]][(i + 1L):s]
    }
    if (p > s) B[i, ] <- rows[[i]][(s + 1L):p]
  }
  list(A = A, B = B, perm = perm)
}

# Build the thresholded Nystrom factor and diagonal correction.
spectral_nystrom_factor <- function(A, B, rel_tol = 1e-10) {
  eig <- eigen((A + t(A)) / 2, symmetric = TRUE)
  cutoff <- rel_tol * max(eig$values)
  keep <- which(eig$values > cutoff)
  if (!length(keep)) stop("A has no eigenvalue above the spectral threshold.",
                          call. = FALSE)

  V <- eig$vectors[, keep, drop = FALSE]
  sigma <- eig$values[keep]
  W <- rbind(A, t(B))

  A_spectral_ginv <- V %*% (t(V) / sigma)
  L <- sweep(W %*% V, 2L, sqrt(sigma), "/")
  d <- 1 - rowSums(L * L)
  numerical_tol <- max(1e-8, 100 * .Machine$double.eps)
  if (min(d) < -numerical_tol) {
    stop("The Nyström diagonal residual is materially negative: min(d) = ",
         signif(min(d), 6), call. = FALSE)
  }
  d[d < 0] <- 0

  list(L = L, d = d, A_spectral_ginv = A_spectral_ginv,
       eigenvalues = eig$values, retained = keep, cutoff = cutoff)
}

# Balance predictor redundancy and response relevance.
choose_alpha <- function(L, d, F) {
  p <- length(F)
  qbar <- (sum(colSums(L)^2) + sum(d)) / (p * p)
  fbar <- mean(F)
  if (!is.finite(fbar) || fbar <= 0) {
    stop("mean(F) must be positive to compute alpha.", call. = FALSE)
  }
  alpha <- qbar / (qbar + fbar)
  list(alpha = alpha, qbar = qbar, fbar = fbar)
}

# Solve the factorized simplex-constrained quadratic program.
solve_factorized_simplex_qp <- function(L, d, F, alpha, solver = NULL) {
  p <- length(F)
  x <- CVXR::Variable(p)
  quadratic <- CVXR::sum_squares(t(L) %*% x) +
    sum(d * CVXR::square(x))
  objective <- CVXR::Minimize(
    0.5 * (1 - alpha) * quadratic - alpha * t(F) %*% x
  )
  problem <- CVXR::Problem(objective, list(x >= 0, sum(x) == 1))
  opt_value <- if (is.null(solver)) CVXR::psolve(problem) else
    CVXR::psolve(problem, solver = solver)
  solve_status <- CVXR::status(problem)
  if (!(solve_status %in% c("optimal", "optimal_inaccurate"))) {
    stop("CVXR failed with status: ", solve_status, call. = FALSE)
  }
  score <- as.numeric(CVXR::value(x))
  score[score < 0 & score > -1e-8] <- 0
  list(score = score, status = solve_status, objective = opt_value)
}

# Fit one hippocampal side and save its ranking and diagnostics.
fit_side <- function(side) {
  screened_path <- find_screened_file(side)
  y_path <- file.path(DATA_DIR, sprintf("%s_res_ori.Rdata", side))
  preferred_x <- sprintf("ADNI_COMBINE_%s_selected", side)

  log_message(side, ": loading saved screening result ", basename(screened_path))
  aligned <- align_X_Y(
    load_rdata_object(screened_path, preferred_x),
    load_rdata_object(y_path, sprintf("%s_res_ori", side))
  )
  X <- aligned$X
  Y <- aligned$Y
  p <- ncol(X)
  if (!is.finite(LANDMARK_PROPORTION) ||
      LANDMARK_PROPORTION <= 0 || LANDMARK_PROPORTION > 1) {
    stop("LANDMARK_PROPORTION must lie in (0, 1].", call. = FALSE)
  }
  s <- min(p, max(2L, floor(LANDMARK_PROPORTION * p)))
  if (s < 2L) stop("At least two screened predictors are required.", call. = FALSE)

  cache_path <- file.path(
    CACHE_DIR,
    sprintf("%s_factor_inputs_s%d.rds", side, s)
  )
  if (file.exists(cache_path)) {
    log_message(side, ": reusing cached F/A/B from ", basename(cache_path))
    cached <- readRDS(cache_path)
    if (!identical(cached$variables, colnames(X)) ||
        !identical(cached$n, nrow(X)) || !identical(cached$s, s)) {
      stop("Cache does not match the current screened data; remove: ", cache_path,
           call. = FALSE)
    }
    F <- cached$F
    blocks <- cached$blocks
  } else {
    set.seed(SEED + match(side, SIDES))
    landmark_idx <- sort(sample.int(p, s, replace = FALSE))
    cl <- make_cluster(N_CORES)
    on.exit(if (!is.null(cl)) parallel::stopCluster(cl), add = TRUE)
    log_message(side, ": computing relevance F on screened predictors")
    F <- compute_relevance(X, Y, cl)
    log_message(side, ": computing landmark blocks A and B")
    blocks <- compute_landmark_blocks(X, landmark_idx, cl)
    saveRDS(list(F = F, blocks = blocks, variables = colnames(X),
                 n = nrow(X), s = s), cache_path)
  }

  fac <- spectral_nystrom_factor(blocks$A, blocks$B, SPECTRAL_REL_TOL)
  F_perm <- F[blocks$perm]
  alpha_info <- choose_alpha(fac$L, fac$d, F_perm)
  log_message(side, ": alpha = ", signif(alpha_info$alpha, 8),
              ", retained rank = ", length(fac$retained))

  opt <- solve_factorized_simplex_qp(
    fac$L, fac$d, F_perm, alpha_info$alpha, CVXR_SOLVER
  )
  score_original <- numeric(p)
  score_original[blocks$perm] <- opt$score
  names(score_original) <- colnames(X)
  ord <- order(score_original, decreasing = TRUE)
  top <- head(ord, min(M_SELECT, p))

  ranking <- data.frame(
    rank = seq_len(p),
    variable = names(score_original)[ord],
    score = score_original[ord],
    relevance = F[ord],
    selected_top_M = seq_len(p) <= min(M_SELECT, p),
    stringsAsFactors = FALSE
  )
  csv_path <- file.path(RESULT_DIR, sprintf("BcorQPFS_%s_top%d.csv", side, M_SELECT))
  rds_path <- file.path(RESULT_DIR, sprintf("BcorQPFS_%s_full_result.rds", side))
  utils::write.csv(ranking, csv_path, row.names = FALSE)
  saveRDS(list(
    side = side,
    screened_input = screened_path,
    selected_variables = names(score_original)[top],
    ranking = ranking,
    alpha = alpha_info,
    spectral = list(
      landmark_proportion = LANDMARK_PROPORTION,
      landmark_count = s,
      relative_tolerance = SPECTRAL_REL_TOL,
      cutoff = fac$cutoff,
      retained_rank = length(fac$retained),
      eigenvalues_A = fac$eigenvalues
    ),
    landmark_indices = blocks$perm[seq_len(s)],
    diagonal_residual = fac$d,
    optimization = opt[c("status", "objective")],
    constraints = c("x >= 0", "sum(x) = 1")
  ), rds_path)
  log_message(side, ": saved ", basename(csv_path), " and ", basename(rds_path))
  invisible(ranking)
}

main <- function() {
  for (side in SIDES) fit_side(side)
  log_message("All requested ADNI BcorQPFS analyses completed.")
}

if (!identical(Sys.getenv("BCORQPFS_SKIP_MAIN"), "1")) main()
