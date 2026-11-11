#!/usr/bin/env Rscript

# Fixed-screening ADNI experiment. Run with Rscript from any working directory.
suppressPackageStartupMessages(library(Matrix))

if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", ".UTF-8"))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
script_dir <- dirname(normalizePath(script, winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "workflow_utils.R"), encoding = "UTF-8")
root <- v2_root(script_dir)
data_dir <- v2_external_data(root)
archive <- v2_bcor_archive(root)
args <- commandArgs(TRUE)
if (!all(args %in% "--smoke")) stop("Supported option: --smoke")
smoke <- "--smoke" %in% args
results <- file.path(root, "results", "adni_stability", paste0(if (smoke) "smoke_" else "run_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(results, recursive = TRUE, showWarnings = FALSE)
writeLines("status=running", file.path(results, "run_status.txt"))
# Reuse the paper's factorization only, without running the archived script.
source_path <- file.path(root, "realdata", "02_bcorqpfs_selection.R")
factor_definition <- Filter(function(expr) {
  is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
    identical(expr[[2L]], as.name("spectral_nystrom_factor"))
}, as.list(parse(source_path)))
stopifnot(length(factor_definition) == 1L)
eval(factor_definition[[1L]])

reps <- if (smoke) 2L else 100L
top_m <- 20L
proportion <- 0.30
spectral_tol <- 1e-10
near_tie_tol <- 1e-7
sides <- c("left", "right")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

log_step <- function(...) {
  message(format(Sys.time(), "%F %T %Z"), " ", paste0(..., collapse = ""))
}
write_result <- function(x, name) {
  write.csv(x, file.path(results, paste0(name, ".csv")), row.names = FALSE)
}
elapsed <- function(start) unname(proc.time()[["elapsed"]] - start)

# Group the empirical V-statistic by the nine joint genotype cells. Each
# ordered pair of anchor cells has weight p(anchor1) * p(anchor2).
# The masks use the paper's closed metric balls, including tied distances.
full_bcor_dosage <- function(X) {
  stopifnot(all(X %in% 0:2))
  q <- ncol(X)
  pairs <- which(upper.tri(matrix(FALSE, q, q), diag = TRUE), arr.ind = TRUE)
  cells <- expand.grid(u = 0:2, v = 0:2)
  indicator <- lapply(0:2, function(g) 1.0 * (X == g))
  mass <- vapply(seq_len(nrow(cells)), function(a) {
    u <- cells$u[a] + 1L
    v <- cells$v[a] + 1L
    crossprod(indicator[[u]], indicator[[v]])[pairs] / nrow(X)
  }, numeric(nrow(pairs)))

  bcov <- numeric(nrow(pairs))
  for (a in seq_len(nrow(cells))) {
    for (b in seq_len(nrow(cells))) {
      in_u <- abs(cells$u - cells$u[a]) <= abs(cells$u[b] - cells$u[a])
      in_v <- abs(cells$v - cells$v[a]) <= abs(cells$v[b] - cells$v[a])
      if (all(in_u) || all(in_v)) next
      delta <- rowSums(mass[, in_u & in_v, drop = FALSE]) -
        rowSums(mass[, in_u, drop = FALSE]) * rowSums(mass[, in_v, drop = FALSE])
      bcov <- bcov + mass[, a] * mass[, b] * delta^2
    }
  }

  self <- bcov[pairs[, 1L] == pairs[, 2L]]
  stopifnot(all(self > 0))
  upper <- bcov / sqrt(self[pairs[, 1L]] * self[pairs[, 2L]])
  Q <- matrix(0, q, q, dimnames = list(colnames(X), colnames(X)))
  Q[pairs] <- upper
  Q <- Q + t(Q)
  diag(Q) <- 1
  Q
}

solve_qp <- function(Q, F) {
  q <- length(F)
  alpha <- mean(Q) / (mean(Q) + mean(F))
  # Positive scaling improves numerical accuracy while preserving minimizers.
  H <- Q * (1 - alpha) * q
  linear <- -alpha * F * q
  start <- proc.time()[["elapsed"]]
  fit <- clarabel::clarabel(
    A = rbind(Matrix(rep(1, q), nrow = 1L, sparse = TRUE), -Diagonal(q)),
    b = c(1, rep(0, q)), q = linear,
    P = forceSymmetric(Matrix(H, sparse = TRUE)),
    cones = list(z = 1L, l = q),
    control = clarabel::clarabel_control(
      verbose = FALSE, tol_gap_abs = 1e-11, tol_gap_rel = 1e-11,
      tol_feas = 1e-11, max_threads = 1L
    )
  )
  status <- names(clarabel::solver_status_descriptions())[fit$status]
  x <- fit$x
  g <- drop(H %*% x) + linear
  simplex_error <- max(abs(sum(x) - 1), -min(x))
  # For a feasible simplex vector this first-order gap bounds suboptimality.
  gap <- (sum(g * x) - min(g)) / q
  stopifnot(status %in% c("Solved", "AlmostSolved"),
            all(is.finite(x)), simplex_error < 1e-8,
            abs(gap) < 1e-8)
  list(x = x, alpha = alpha, status = status, simplex_error = simplex_error,
       gap = gap, seconds = elapsed(start))
}

get_order <- function(x) order(-x, seq_along(x))
boundary_gap <- function(x) {
  ord <- get_order(x)
  x[ord[top_m]] - x[ord[top_m + 1L]]
}
objective <- function(x, Q, F, alpha) {
  (1 - alpha) * sum(x * drop(Q %*% x)) / 2 - alpha * sum(F * x)
}
pair_stats <- function(a, b) {
  n <- length(intersect(a, b))
  c(overlap = n, jaccard = n / length(union(a, b)))
}

run_side <- function(side, side_index) {
  start <- proc.time()[["elapsed"]]
  log_step(side, ": loading frozen screening and archival relevance")
  old_path <- file.path(archive, sprintf("BcorQPFS_%s_full_result.rds", side))
  old <- readRDS(old_path)
  screened <- old$screened_input
  if (is.null(screened) || !file.exists(screened)) {
    screened <- tryCatch(v2_data_file(root, sprintf("ADNI_%s_screening_selected.RData", side)),
                         error = function(e) v2_data_file(root, sprintf("ADNI_COMBINE_%s_selected.RData", side)))
  }
  env <- new.env(parent = emptyenv()); objects <- load(screened, envir = env)
  xname <- if ("selected_X" %in% objects) "selected_X" else sprintf("ADNI_COMBINE_%s_selected", side)
  X <- as.matrix(env[[xname]])
  variables <- colnames(X)
  # A newly generated 01 result has companion metadata; archived COMBINE files carry these fields inside RData.
  meta_path <- file.path(dirname(screened), sprintf("ADNI_%s_screening_metadata.rds", side))
  if (file.exists(meta_path)) {
    metadata <- readRDS(meta_path)
    indices <- metadata$original_index[match(variables, metadata$variables)]
    threshold <- metadata$threshold
  } else {
    indices <- env[[paste0(side, "_keep_idx")]]
    threshold <- env[[paste0(side, "_thr")]]
  }
  q <- ncol(X); s <- as.integer(floor(proportion * q))
  stopifnot(q >= top_m, nrow(X) >= 2L, length(indices) == q,
            !anyNA(indices), !anyDuplicated(indices), !anyDuplicated(variables),
            length(threshold)==1L, is.finite(threshold))
  cached_path <- file.path(archive, "cache", sprintf("%s_factor_inputs_s%d.rds", side, s))
  dy_path <- file.path(archive, "cache", sprintf("%s_DY_dist.rds", side))
  if (!file.exists(dy_path)) dy_path <- v2_data_file(root, sprintf("%s_DY_dist.rds", side))
  cached <- readRDS(cached_path)
  F <- cached$F
  original_rank <- match(variables, old$ranking$variable)
  stopifnot(identical(cached$variables, variables), identical(names(F), variables),
            cached$n == nrow(X), cached$s == s, !anyNA(original_rank),
            all(is.finite(F)), all(F >= threshold),
            max(abs(F - old$ranking$relevance[original_rank])) < 1e-12)

  # Check the archived F against saved phenotype distances at spaced columns.
  # The frozen SNP matrices inherit their archival row order.
  DY <- readRDS(dy_path)
  stopifnot(attr(DY, "Size") == nrow(X))
  probes <- unique(as.integer(round(seq(1, q, length.out = 8L))))
  f_check <- vapply(probes, function(j) {
    as.numeric(Ball::bcor(dist(X[, j, drop = FALSE]), DY, distance = TRUE))
  }, numeric(1L))
  f_error <- max(abs(f_check - F[probes]))
  stopifnot(f_error < 1e-12)
  write_result(data.frame(side, candidate_index = seq_len(q), original_index = indices,
                          variable = variables, relevance = F,
                          screening_threshold = threshold),
               paste0(side, "_frozen_candidates"))

  q_start <- proc.time()[["elapsed"]]
  Q <- full_bcor_dosage(X)
  q_seconds <- elapsed(q_start)
  a <- cached$blocks$perm[seq_len(s)]
  b <- cached$blocks$perm[-seq_len(s)]
  historical_error <- max(abs(Q[a, a] - cached$blocks$A),
                          abs(Q[a, b] - cached$blocks$B))
  set.seed(716L + side_index)
  probe_pairs <- replicate(24L, sample.int(q, 2L))
  direct <- apply(probe_pairs, 2L, function(ij) {
    as.numeric(Ball::bcor(X[, ij[1L]], X[, ij[2L]]))
  })
  q_error <- max(abs(direct - Q[t(probe_pairs)]))
  spectrum <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
  log_step(side, ": Q validation errors, archival=", historical_error, ", direct=", q_error, ", min eigenvalue=", min(spectrum))
  stopifnot(historical_error < 1e-9, q_error < 1e-9, min(spectrum) > -1e-10)
  saveRDS(Q, file.path(results, paste0(side, "_full_Q.rds")))
  log_step(side, ": full Q checked in ", q_seconds, " s (", choose(q, 2), " pairs)")

  ld <- cor(X)^2
  count_ld <- function(top) sum(ld[top, top][upper.tri(matrix(FALSE, top_m, top_m))] >= 0.8)
  exact <- solve_qp(Q, F)
  exact_rank <- match(seq_len(q), get_order(exact$x))
  exact_top <- which(exact_rank <= top_m)
  old_top <- which(original_rank <= top_m)
  exact_obj <- objective(exact$x, Q, F, exact$alpha)

  # Keep the archived landmark draw as a separate reference, outside the repeats.
  seeds <- 202611140L + (side_index - 1L) * 1000L + seq_len(reps)
  samples <- c(list(a), lapply(seeds, function(seed) {
    set.seed(seed)
    sort(sample.int(q, s))
  }))
  stopifnot(!anyDuplicated(vapply(samples[-1L], paste, "", collapse = ",")))

  scores <- ranks <- matrix(NA_real_, q, length(samples))
  metadata <- vector("list", length(samples))
  for (k in seq_along(samples)) {
    landmarks <- samples[[k]]
    rest <- setdiff(seq_len(q), landmarks)
    perm <- c(landmarks, rest)
    factor_start <- proc.time()[["elapsed"]]
    fac <- spectral_nystrom_factor(Q[landmarks, landmarks], Q[landmarks, rest],
                                  spectral_tol)
    approx <- tcrossprod(fac$L) + diag(fac$d)
    factor_seconds <- elapsed(factor_start)
    fit <- solve_qp(approx, F[perm])
    x <- numeric(q)
    x[perm] <- fit$x
    scores[, k] <- x
    ranks[, k] <- match(seq_len(q), get_order(x))
    top <- which(ranks[, k] <= top_m)
    err <- Q[perm, perm] - approx
    d <- sort(fac$d, decreasing = TRUE)
    stats <- pair_stats(top, exact_top)
    metadata[[k]] <- data.frame(
      side, run = k - 1L, seed = c(NA_integer_, seeds)[k],
      q, s, retained_rank = length(fac$retained), alpha = fit$alpha,
      entrywise_error = max(abs(err)),
      entrywise_bound = sqrt(d[1L] * d[2L]),
      spectral_bound = sqrt(max(0, sum(d)^2 - sum(d^2))),
      relative_frobenius_error = norm(err, "F") / norm(Q, "F"),
      overlap_exact = unname(stats["overlap"]),
      jaccard_exact = unname(stats["jaccard"]),
      overlap_archived = length(intersect(top, old_top)),
      rank_spearman_exact = cor(ranks[, k], exact_rank),
      top20_gap = boundary_gap(x), high_ld_pairs = count_ld(top),
      exact_objective_excess = objective(x, Q, F, exact$alpha) - exact_obj,
      factor_seconds, qp_seconds = fit$seconds, status = fit$status,
      simplex_error = fit$simplex_error, first_order_gap = fit$gap
    )
    stopifnot(metadata[[k]]$entrywise_error <= metadata[[k]]$entrywise_bound + 1e-8,
              metadata[[k]]$exact_objective_excess >= -1e-8)
    if ((k - 1L) %% 10L == 0L || k == length(samples)) {
      log_step(side, ": finished ", k - 1L, "/", reps, " new draws")
    }
  }
  draw_summary <- do.call(rbind, metadata)
  repeat_scores <- scores[, -1L, drop = FALSE]
  repeat_ranks <- ranks[, -1L, drop = FALSE]
  mean_score <- rowMeans(repeat_scores)
  mean_rank <- match(seq_len(q), get_order(mean_score))
  top_sets <- lapply(seq_len(reps), function(j) which(repeat_ranks[, j] <= top_m))
  pair_grid <- combn(reps, 2L)
  pairwise <- t(apply(pair_grid, 2L, function(ij) {
    pair_stats(top_sets[[ij[1L]]], top_sets[[ij[2L]]])
  }))
  pairwise <- data.frame(side, first = pair_grid[1L, ], second = pair_grid[2L, ], pairwise)
  snps <- data.frame(
    side, candidate_index = seq_len(q), original_index = indices, variable = variables,
    relevance = F, exact_rank, exact_weight = exact$x, archived_rank = original_rank,
    archived_refit_rank = ranks[, 1L], mean_weight_rank = mean_rank, mean_weight = mean_score,
    frequency = rowMeans(repeat_ranks <= top_m),
    mean_rank = rowMeans(repeat_ranks),
    min_rank = apply(repeat_ranks, 1L, min), max_rank = apply(repeat_ranks, 1L, max)
  )
  exact_boundary <- sort(exact$x, decreasing = TRUE)[top_m]
  snps$near_exact_boundary <- abs(exact$x - exact_boundary) <= near_tie_tol
  snps <- snps[order(snps$exact_rank), ]
  snps$top20_count <- as.integer(round(snps$frequency * reps))
  snps$top20_fraction <- snps$top20_count / reps
  stopifnot(sum(snps$top20_count) == top_m * reps)
  frequency_order <- order(-snps$top20_count, snps$mean_rank, snps$original_index)
  snps$frequency_rank <- match(seq_len(q), frequency_order)
  write_result(snps[seq_len(top_m),
                    c("side", "exact_rank", "variable", "original_index", "exact_weight",
                      "relevance", "archived_rank", "top20_count", "top20_fraction",
                      "frequency_rank", "mean_rank")],
               paste0(side, "_exact_qp_top20"))
  write_result(snps[frequency_order,
                    c("side", "frequency_rank", "variable", "original_index", "top20_count",
                      "top20_fraction", "mean_rank", "min_rank", "max_rank", "mean_weight",
                      "mean_weight_rank", "exact_rank", "archived_rank")],
               paste0(side, "_landmark_frequency_ranking"))
  ld_idx <- which(upper.tri(ld) & ld >= 0.8, arr.ind = TRUE)
  ld_pairs <- data.frame(
    side = rep(side, nrow(ld_idx)), variable1 = variables[ld_idx[, 1L]], variable2 = variables[ld_idx[, 2L]],
    r2 = ld[ld_idx], squared_bcor = Q[ld_idx],
    exact_both = exact_rank[ld_idx[, 1L]] <= top_m & exact_rank[ld_idx[, 2L]] <= top_m,
    archived_both = original_rank[ld_idx[, 1L]] <= top_m &
      original_rank[ld_idx[, 2L]] <= top_m,
    repeat_both_frequency = rowMeans(repeat_ranks[ld_idx[, 1L], , drop = FALSE] <= top_m &
      repeat_ranks[ld_idx[, 2L], , drop = FALSE] <= top_m)
  )
  write_result(snps, paste0(side, "_snp_stability"))
  write_result(draw_summary, paste0(side, "_draw_diagnostics"))
  write_result(pairwise, paste0(side, "_pairwise_stability"))
  write_result(ld_pairs, paste0(side, "_high_ld_pairs"))
  write_result(data.frame(
    side, run = rep(seq_along(samples) - 1L, each = q),
    variable = rep(variables, length(samples)), rank = as.vector(ranks),
    weight = as.vector(scores), landmark = as.vector(vapply(
      samples, function(i) seq_len(q) %in% i, logical(q)
    ))
  ), paste0(side, "_draw_rankings"))
  write_result(data.frame(side, n = nrow(X), q, screening_threshold = threshold,
                          min_relevance = min(F), s, draws = reps, spectral_tol,
                          near_tie_tol, f_error, q_error, historical_error,
                          min_full_eigenvalue = min(spectrum), full_Q_seconds = q_seconds,
                          exact_QP_seconds = exact$seconds, exact_alpha = exact$alpha,
                          exact_status = exact$status, exact_gap = exact$gap,
                          exact_top20_gap = boundary_gap(exact$x),
                          exact_high_ld_pairs = count_ld(exact_top),
                          archived_high_ld_pairs = count_ld(old_top),
                          exact_overlap_archived = length(intersect(exact_top, old_top)),
                          mean_weight_overlap_exact = sum(exact_rank <= top_m & mean_rank <= top_m),
                          mean_weight_top20_gap = boundary_gap(mean_score),
                          elapsed_seconds = elapsed(start)), paste0(side, "_side_summary"))
  provenance_paths <- c(screened, cached_path, old_path, dy_path, source_path)
  write_result(data.frame(file = normalizePath(provenance_paths, winslash = "/", mustWork = TRUE),
                          md5 = unname(tools::md5sum(provenance_paths))),
               paste0(side, "_input_provenance"))
  log_step(side, ": results saved")
}

for (i in seq_along(sides)) run_side(sides[i], i)
capture.output(sessionInfo(), file = file.path(results, "session_info.txt"))
writeLines("status=completed", file.path(results, "run_status.txt"))
log_step("All requested fits completed")
