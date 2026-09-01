# Shared functions for the simulation studies.

check_simulation_packages <- function(include_hsic = FALSE) {
  packages <- c("Ball", "CVXR")
  if (include_hsic) packages <- c(packages, "glmnet")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing R package(s): ", paste(missing, collapse = ", "),
      ". Install them before running the simulation.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

standardize_columns <- function(X) {
  X <- as.matrix(X)
  centers <- colMeans(X)
  scales <- apply(X, 2L, stats::sd)
  scales[!is.finite(scales) | scales <= 0] <- 1
  sweep(sweep(X, 2L, centers, "-"), 2L, scales, "/")
}

bcor_scalar <- function(x, y) {
  value <- suppressWarnings(as.numeric(Ball::bcor(x, y))[1L])
  if (!is.finite(value)) 0 else min(1, max(0, value))
}

# Generate predictors from the mixed Hermite design.
simulate_mixed_hermite_example <- function(
    example, n, p, kappa, seed, rho = 0.8) {
  if (!example %in% 1:3) stop("example must be 1, 2, or 3.", call. = FALSE)
  if (n < 10L || p < 5L) stop("n and p are too small.", call. = FALSE)
  if (!is.finite(kappa) || kappa < 0 || kappa > 1) {
    stop("kappa must lie in [0, 1].", call. = FALSE)
  }
  if (!is.finite(rho) || abs(rho) >= 1) {
    stop("rho must lie in (-1, 1).", call. = FALSE)
  }

  set.seed(seed)
  innovations <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  Z <- innovations
  if (p > 1L) {
    innovation_scale <- sqrt(1 - rho^2)
    for (j in 2:p) {
      Z[, j] <- rho * Z[, j - 1L] + innovation_scale * innovations[, j]
    }
  }

  H2 <- (Z^2 - 1) / sqrt(2)
  X_raw <- sqrt(1 - kappa) * Z + sqrt(kappa) * H2
  epsilon <- stats::rnorm(n)

  if (example == 1L) {
    omega <- 3 * X_raw[, 2L]^2 - 3 * X_raw[, 3L]^2 +
      5 * X_raw[, 5L] + epsilon
    Y <- 1 + 2 * as.numeric(omega > 0) * omega
    active <- c("X2", "X3", "X5")
  } else if (example == 2L) {
    omega <- 3 * X_raw[, 1L] * X_raw[, 2L] + 5 * X_raw[, 5L] + epsilon
    Y <- 1 + 2 * as.numeric(omega < 0) * omega
    active <- c("X1", "X2", "X5")
  } else {
    exponent <- X_raw[, 2L]^2
    exponent_shift <- max(exponent)
    response_scale <- exp(-exponent_shift)
    Y <- (2 * sin(3 * X_raw[, 1L] * X_raw[, 3L]) + epsilon) *
      response_scale - exp(exponent - exponent_shift)
    active <- c("X1", "X2", "X3")
  }

  response_scale <- max(abs(Y))
  if (is.finite(response_scale) && response_scale > 0) Y <- Y / response_scale
  X <- standardize_columns(X_raw)
  colnames(X) <- paste0("X", seq_len(p))
  list(X = X, Y = as.numeric(Y), active = active)
}

# Compute predictor relevance scores.
compute_relevance <- function(X, Y) {
  X <- as.matrix(X)
  scores <- vapply(
    seq_len(ncol(X)),
    function(j) bcor_scalar(X[, j], Y),
    numeric(1)
  )
  names(scores) <- colnames(X)
  scores
}

# Compute the predictor Ball-correlation matrix.
compute_full_bcor_matrix <- function(X) {
  X <- as.matrix(X)
  p <- ncol(X)
  Q <- diag(1, p)
  if (p > 1L) {
    for (j in seq_len(p - 1L)) {
      values <- vapply(
        (j + 1L):p,
        function(ell) bcor_scalar(X[, j], X[, ell]),
        numeric(1)
      )
      Q[j, (j + 1L):p] <- values
      Q[(j + 1L):p, j] <- values
    }
  }
  dimnames(Q) <- list(colnames(X), colnames(X))
  Q
}

summarize_bcor_matrix <- function(Q) {
  Q <- as.matrix(Q)
  p <- nrow(Q)
  list(
    qbar_all = min(1, max(0, sum(Q) / p^2)),
    qbar_offdiag = if (p > 1L) {
      min(1, max(0, (sum(Q) - sum(diag(Q))) / (p * (p - 1L))))
    } else {
      0
    }
  )
}

# Extract uniformly sampled landmark columns.
sample_landmark_columns <- function(Q, landmark_size, seed) {
  Q <- as.matrix(Q)
  p <- nrow(Q)
  landmark_size <- as.integer(landmark_size)
  if (ncol(Q) != p || landmark_size < 1L || landmark_size > p) {
    stop("Invalid landmark matrix dimensions.", call. = FALSE)
  }
  set.seed(seed)
  landmarks <- sample.int(p, landmark_size, replace = FALSE)
  list(
    W = Q[, landmarks, drop = FALSE],
    A = Q[landmarks, landmarks, drop = FALSE],
    landmarks = landmarks
  )
}

# Compute the thresholded positive spectral inverse.
positive_spectral_inverse <- function(A, tau = 1e-10) {
  if (!is.finite(tau) || tau <= 0 || tau >= 1) {
    stop("tau must lie in (0, 1).", call. = FALSE)
  }
  eig <- eigen(0.5 * (A + t(A)), symmetric = TRUE)
  sigma <- as.numeric(eig$values)
  sigma_max <- max(sigma, 0)
  if (!is.finite(sigma_max) || sigma_max <= 0) {
    stop("The landmark block has no positive eigenvalue.", call. = FALSE)
  }
  threshold <- tau * sigma_max
  keep <- which(sigma > threshold)
  if (!length(keep)) stop("The retained rank is zero.", call. = FALSE)
  list(
    vectors = eig$vectors[, keep, drop = FALSE],
    values = sigma[keep],
    threshold = threshold,
    effective_rank = length(keep),
    minimum_eigenvalue = min(sigma)
  )
}

# Build the diagonally corrected Nystrom factorization.
build_nystrom_approximation <- function(
    Q, landmark_size, seed, tau = 1e-10, diagonal_tolerance = 1e-8) {
  sampled <- sample_landmark_columns(Q, landmark_size, seed)
  spectral <- positive_spectral_inverse(sampled$A, tau)
  L <- sampled$W %*% spectral$vectors
  L <- sweep(L, 2L, sqrt(spectral$values), "/")

  raw_diagonal <- rowSums(L^2)
  raw_correction <- 1 - raw_diagonal
  if (any(raw_correction < -diagonal_tolerance)) {
    stop("The Nystrom diagonal exceeds one beyond tolerance.", call. = FALSE)
  }
  row_scale <- 1 / sqrt(pmax(raw_diagonal, 1))
  L <- sweep(L, 1L, row_scale, "*")
  diagonal_correction <- pmax(1 - rowSums(L^2), 0)

  p <- nrow(L)
  one_L <- as.numeric(crossprod(L, rep(1, p)))
  total_sum <- sum(one_L^2) + sum(diagonal_correction)
  diagonal_sum <- sum(rowSums(L^2) + diagonal_correction)
  qbar_all <- total_sum / p^2
  qbar_offdiag <- if (p > 1L) {
    (total_sum - diagonal_sum) / (p * (p - 1L))
  } else {
    0
  }

  list(
    L = L,
    diagonal_correction = diagonal_correction,
    landmarks = sampled$landmarks,
    qbar_all = min(1, max(0, qbar_all)),
    qbar_offdiag = min(1, max(0, qbar_offdiag)),
    effective_rank = spectral$effective_rank,
    minimum_eigenvalue = spectral$minimum_eigenvalue,
    minimum_raw_diagonal_correction = min(raw_correction)
  )
}

# Select the data-adaptive trade-off parameter.
choose_balanced_alpha <- function(qbar_all, relevance, endpoint = 1e-6) {
  denominator <- qbar_all + mean(relevance)
  raw_alpha <- if (is.finite(denominator) && denominator > 0) {
    qbar_all / denominator
  } else {
    0.5
  }
  min(1 - endpoint, max(endpoint, raw_alpha))
}

solve_cvxr_problem <- function(problem, variable, solver = "CLARABEL") {
  installed <- CVXR::installed_solvers()
  candidates <- unique(c(toupper(solver), "CLARABEL", "OSQP", "SCS"))
  candidates <- candidates[candidates %in% installed]
  if (!length(candidates)) stop("No supported CVXR solver is installed.", call. = FALSE)

  failures <- character()
  for (candidate in candidates) {
    result <- tryCatch({
      raw_result <- CVXR::psolve(problem, solver = candidate)
      legacy_api <- is.list(raw_result) && !is.null(raw_result$status) &&
        is.function(raw_result$getValue)
      if (legacy_api) {
        status <- raw_result$status
        value <- raw_result$getValue(variable)
      } else {
        status <- CVXR::status(problem)
        value <- CVXR::value(variable)
      }
      status <- tolower(as.character(status)[1L])
      if (!status %in% c("optimal", "optimal_inaccurate")) {
        stop("Solver status: ", status, call. = FALSE)
      }
      value <- as.numeric(value)
      if (!length(value) || any(!is.finite(value))) {
        stop("The solver returned invalid values.", call. = FALSE)
      }
      list(value = value, status = status, solver = candidate)
    }, error = function(e) {
      failures <<- c(failures, paste(candidate, conditionMessage(e), sep = ": "))
      NULL
    })
    if (!is.null(result)) return(result)
  }
  stop("CVXR failed: ", paste(failures, collapse = " | "), call. = FALSE)
}

check_cvxr_solver <- function(solver = "CLARABEL") {
  x <- CVXR::Variable(2L)
  problem <- CVXR::Problem(
    CVXR::Minimize(CVXR::sum_squares(x)),
    constraints = list(x >= 0, sum(x) == 1)
  )
  result <- solve_cvxr_problem(problem, x, solver)
  if (abs(sum(result$value) - 1) > 1e-6 || min(result$value) < -1e-6) {
    stop("The CVXR simplex check failed.", call. = FALSE)
  }
  invisible(result)
}

# Solve the exact bcorQPFS problem.
solve_exact_bcorqpfs <- function(
    relevance, Q, alpha, solver = "CLARABEL", tolerance = 1e-6) {
  variable_names <- names(relevance)
  relevance <- as.numeric(relevance)
  Q <- 0.5 * (as.matrix(Q) + t(as.matrix(Q)))
  p <- length(relevance)
  x <- CVXR::Variable(p)
  objective <- CVXR::Minimize(
    0.5 * (1 - alpha) * CVXR::quad_form(x, Q, assume_PSD = TRUE) -
      alpha * sum(relevance * x)
  )
  problem <- CVXR::Problem(objective, constraints = list(x >= 0, sum(x) == 1))
  solution <- solve_cvxr_problem(problem, x, solver)
  score <- solution$value
  names(score) <- variable_names
  if (abs(sum(score) - 1) > tolerance || min(score) < -tolerance) {
    stop("The exact solution violates the simplex tolerance.", call. = FALSE)
  }
  list(
    score = score,
    ranking = names(sort(score, decreasing = TRUE)),
    status = solution$status,
    solver = solution$solver
  )
}

# Solve the diagonally corrected Nystrom problem.
solve_nystrom_bcorqpfs <- function(
    relevance, approximation, alpha, solver = "CLARABEL", tolerance = 1e-6) {
  variable_names <- names(relevance)
  relevance <- as.numeric(relevance)
  L <- approximation$L
  diagonal_correction <- approximation$diagonal_correction
  p <- length(relevance)
  x <- CVXR::Variable(p)
  objective <- CVXR::Minimize(
    0.5 * (1 - alpha) * (
      CVXR::sum_squares(t(L) %*% x) +
        CVXR::sum_squares(sqrt(diagonal_correction) * x)
    ) - alpha * sum(relevance * x)
  )
  problem <- CVXR::Problem(objective, constraints = list(x >= 0, sum(x) == 1))
  solution <- solve_cvxr_problem(problem, x, solver)
  score <- solution$value
  names(score) <- variable_names
  if (abs(sum(score) - 1) > tolerance || min(score) < -tolerance) {
    stop("The Nystrom solution violates the simplex tolerance.", call. = FALSE)
  }
  list(
    score = score,
    ranking = names(sort(score, decreasing = TRUE)),
    status = solution$status,
    solver = solution$solver
  )
}

# Rank predictors with Gaussian mutual-information mRMR.
rank_mrmr <- function(X, Y, M) {
  X <- as.matrix(X)
  M <- min(as.integer(M), ncol(X))
  gaussian_mi <- function(correlation) {
    correlation <- pmin(1 - 1e-12, pmax(-1 + 1e-12, correlation))
    -0.5 * log1p(-correlation^2)
  }
  relevance_cor <- vapply(
    seq_len(ncol(X)),
    function(j) stats::cor(X[, j], Y),
    numeric(1)
  )
  relevance_cor[!is.finite(relevance_cor)] <- 0
  relevance <- gaussian_mi(relevance_cor)
  selected <- integer(0)
  remaining <- seq_len(ncol(X))

  while (length(selected) < M) {
    if (!length(selected)) {
      criterion <- relevance[remaining]
    } else {
      redundancy_cor <- stats::cor(
        X[, remaining, drop = FALSE],
        X[, selected, drop = FALSE]
      )
      redundancy_cor[!is.finite(redundancy_cor)] <- 0
      redundancy <- gaussian_mi(redundancy_cor)
      if (is.null(dim(redundancy))) {
        redundancy <- matrix(redundancy, nrow = length(remaining))
      }
      criterion <- relevance[remaining] - rowMeans(redundancy)
    }
    chosen <- remaining[which.max(criterion)]
    selected <- c(selected, chosen)
    remaining <- setdiff(remaining, chosen)
  }
  colnames(X)[selected]
}

standardize_hsic_input <- function(x) {
  x <- as.matrix(x)
  centers <- colMeans(x)
  scales <- apply(x, 2L, stats::sd)
  if (any(!is.finite(scales) | scales <= 0)) {
    stop("HSIC Lasso received a constant variable.", call. = FALSE)
  }
  sweep(sweep(x, 2L, centers, "-"), 2L, scales, "/")
}

matlab_median_sigma <- function(values) {
  values <- as.numeric(values)
  if (length(values) > 500L) values <- values[seq_len(500L)]
  distances <- outer(values, values, function(a, b) (a - b)^2)
  positive <- distances[upper.tri(distances) & distances > 0]
  if (!length(positive)) stop("The kernel bandwidth is undefined.", call. = FALSE)
  sqrt(0.5 * stats::median(positive))
}

center_kernel <- function(K) {
  sweep(sweep(K, 1L, rowMeans(K), "-"), 2L, colMeans(K), "-") + mean(K)
}

gaussian_kernel <- function(values) {
  sigma <- matlab_median_sigma(values)
  distances <- outer(values, values, function(a, b) (a - b)^2)
  exp(-distances / (2 * sigma^2))
}

# Rank predictors with full-kernel nonnegative HSIC Lasso.
rank_hsic_lasso <- function(X, Y, M, nlambda = 120L) {
  X <- as.matrix(X)
  M <- min(as.integer(M), ncol(X))
  Xs <- standardize_hsic_input(X)
  Ys <- as.numeric(standardize_hsic_input(matrix(Y, ncol = 1L)))
  n <- nrow(Xs)
  p <- ncol(Xs)
  design <- matrix(0, nrow = n * n, ncol = p)
  colnames(design) <- colnames(Xs)
  for (j in seq_len(p)) {
    design[, j] <- as.vector(center_kernel(gaussian_kernel(Xs[, j])))
  }
  response <- as.vector(center_kernel(gaussian_kernel(Ys)))

  fit <- glmnet::glmnet(
    x = design,
    y = response,
    family = "gaussian",
    alpha = 1,
    intercept = FALSE,
    standardize = FALSE,
    lower.limits = 0,
    nlambda = max(120L, as.integer(nlambda)),
    lambda.min.ratio = 1e-6,
    control = list(thresh = 1e-10, maxit = 1000000L)
  )
  beta_path <- as.matrix(fit$beta)
  support_size <- colSums(beta_path > 1e-10)
  eligible <- which(support_size >= M)
  path_index <- if (length(eligible)) {
    smallest_support <- min(support_size[eligible])
    eligible[which(support_size[eligible] == smallest_support)[1L]]
  } else {
    ncol(beta_path)
  }
  coefficients <- beta_path[, path_index]
  names(coefficients) <- rownames(beta_path)
  ranking <- names(sort(coefficients[coefficients > 1e-10], decreasing = TRUE))
  if (length(ranking) < M) {
    marginal_hsic <- as.numeric(crossprod(design, response))
    names(marginal_hsic) <- colnames(design)
    ranking <- c(
      ranking,
      head(setdiff(names(sort(marginal_hsic, decreasing = TRUE)), ranking), M - length(ranking))
    )
  }
  list(
    ranking = head(ranking, M),
    lambda = fit$lambda[path_index] * nrow(design),
    nonzero = support_size[path_index]
  )
}

# Record individual and joint active-variable recovery.
selection_record <- function(ranking, active, M) {
  selected <- head(ranking, M)
  row <- data.frame(M = as.integer(M), stringsAsFactors = FALSE)
  for (variable in c("X1", "X2", "X3", "X5")) {
    row[[variable]] <- if (variable %in% active) {
      as.integer(variable %in% selected)
    } else {
      NA_integer_
    }
  }
  row$All <- as.integer(all(active %in% selected))
  row$selected <- paste(selected, collapse = ";")
  row
}

measure_call <- function(expression) {
  start <- proc.time()[["elapsed"]]
  value <- force(expression)
  list(value = value, seconds = unname(proc.time()[["elapsed"]] - start))
}

bind_result_component <- function(results, name) {
  pieces <- lapply(results, `[[`, name)
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (!length(pieces)) NULL else do.call(rbind, pieces)
}

safe_mean <- function(x) {
  if (!length(x) || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

# Run simulation jobs sequentially or with a local cluster.
run_job_grid <- function(
    jobs, worker_function, config, workers, functions_file,
    checkpoint_function = NULL) {
  workers <- max(1L, min(as.integer(workers), nrow(jobs)))
  results <- vector("list", nrow(jobs))
  if (workers == 1L) {
    for (i in seq_len(nrow(jobs))) {
      message(sprintf("Job %d of %d", i, nrow(jobs)))
      results[[i]] <- worker_function(jobs[i, , drop = FALSE], config)
      if (!is.null(checkpoint_function)) checkpoint_function(results, i)
    }
    return(results)
  }

  cluster <- parallel::makeCluster(workers)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterExport(cluster, "functions_file", envir = environment())
  parallel::clusterEvalQ(cluster, source(functions_file, encoding = "UTF-8"))
  batches <- split(seq_len(nrow(jobs)), ceiling(seq_len(nrow(jobs)) / workers))
  completed <- 0L
  for (indices in batches) {
    batch_jobs <- lapply(indices, function(i) jobs[i, , drop = FALSE])
    results[indices] <- parallel::parLapplyLB(
      cluster, batch_jobs, worker_function, config = config
    )
    completed <- completed + length(indices)
    message(sprintf("Completed %d of %d jobs", completed, nrow(jobs)))
    if (!is.null(checkpoint_function)) checkpoint_function(results, completed)
  }
  results
}
