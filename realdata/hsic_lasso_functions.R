# Utilities for HSIC Lasso variable selection.

options(stringsAsFactors = FALSE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Manage packages, paths, input objects, and parallel settings.
install_if_missing <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    install.packages(miss, dependencies = TRUE)
  }
}

install_if_missing(c("openxlsx", "pbapply"))

get_script_dir <- function() {

  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(sub("^--file=", "", file_arg[1])))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    ok <- tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE)
    if (isTRUE(ok)) {
      p <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                    error = function(e) "")
      if (nzchar(p)) {
        return(dirname(p))
      }
    }
  }

  getwd()
}

find_file_case_insensitive <- function(dir_path, target_basename) {
  if (!dir.exists(dir_path)) {
    stop("目录不存在：", dir_path)
  }

  all_files <- list.files(dir_path, full.names = TRUE)
  hit <- all_files[tolower(basename(all_files)) == tolower(target_basename)]

  if (length(hit) == 0) {
    stop("在目录中找不到文件：", target_basename, "\n目录：", dir_path)
  }

  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

load_rdata_auto <- function(path, preferred_name = NULL) {
  e <- new.env(parent = emptyenv())
  objs <- load(path, envir = e)

  if (!is.null(preferred_name) && preferred_name %in% objs) {
    return(e[[preferred_name]])
  }

  is_good <- vapply(objs, function(nm) {
    obj <- e[[nm]]
    is.matrix(obj) || is.data.frame(obj)
  }, logical(1))

  if (!any(is_good)) {
    stop("文件中没有 matrix/data.frame 对象：", path,
         "\n对象有：", paste(objs, collapse = ", "))
  }

  e[[objs[which(is_good)[1]]]]
}

as_numeric_matrix <- function(x, name = "object") {
  if (is.data.frame(x)) {
    bad <- !vapply(x, is.numeric, logical(1))
    if (any(bad)) {
      stop(name, " 中存在非数值列：",
           paste(names(x)[bad], collapse = ", "))
    }
    x <- as.matrix(x)
  } else {
    x <- as.matrix(x)
  }

  old_na <- is.na(x)
  suppressWarnings(storage.mode(x) <- "double")
  new_na <- is.na(x)

  if (any(new_na & !old_na)) {
    stop(name, " 转成 numeric matrix 时引入了新的 NA。")
  }

  if (!is.null(colnames(x))) colnames(x) <- trimws(colnames(x))
  if (!is.null(rownames(x))) rownames(x) <- trimws(rownames(x))

  x
}

check_unique_colnames <- function(x, name = "object") {
  cn <- colnames(x)
  if (is.null(cn)) {
    stop(name, " 没有列名，无法做变量名匹配。")
  }
  if (anyDuplicated(cn)) {
    dupn <- unique(cn[duplicated(cn)])
    stop(name, " 存在重复列名，无法唯一匹配。重复示例：",
         paste(head(dupn, 10), collapse = ", "))
  }
}

choose_n_cores <- function(n_tasks,
                           reserve_cores = 1L,
                           use_physical_cores = TRUE) {

  dc <- tryCatch(
    parallel::detectCores(logical = !use_physical_cores),
    error = function(e) NA_integer_
  )

  if (is.na(dc) || dc <= 1L) {
    return(1L)
  }

  n_cores <- max(1L, dc - reserve_cores)
  n_cores <- min(n_cores, n_tasks)
  as.integer(n_cores)
}

# Construct normalized block-kernel representations.
kernel_gaussian_py <- function(X1, X2, sigma) {
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)

  n1 <- ncol(X1)
  n2 <- ncol(X2)

  x12 <- colSums(X1^2)
  x22 <- colSums(X2^2)

  dist2 <- matrix(x22, nrow = n1, ncol = n2, byrow = TRUE) +
    matrix(x12, nrow = n1, ncol = n2, byrow = FALSE) -
    2 * crossprod(X1, X2)

  exp(-dist2 / (2 * sigma^2))
}

kernel_delta_norm_py <- function(X1, X2) {
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)

  if (nrow(X1) != 1L || nrow(X2) != 1L) {
    stop("Delta kernel 这里只用于单行离散标签。")
  }

  v1 <- as.vector(X1)
  v2 <- as.vector(X2)

  n1 <- length(v1)
  n2 <- length(v2)
  K <- matrix(0, nrow = n1, ncol = n2)

  u_list <- unique(v1)
  for (ind in u_list) {
    ind1 <- which(v1 == ind)
    ind2 <- which(v2 == ind)
    c1 <- sqrt(length(ind1))
    c2 <- sqrt(length(ind2))
    if (length(ind1) > 0 && length(ind2) > 0) {
      K[ind1, ind2] <- 1 / c1 / c2
    }
  }
  K
}

compute_kernel_py <- function(x, kernel = c("Gaussian", "Delta"),
                              B = 0L, M = 1L, discarded = NULL) {
  kernel <- match.arg(kernel)
  x <- as.matrix(x)

  d <- nrow(x)
  n <- ncol(x)

  if (B <= 0L || B > n) B <- n
  if (is.null(discarded)) discarded <- n %% B

  H <- diag(B) - matrix(1 / B, nrow = B, ncol = B)

  numblocks <- n %/% B

  M_eff <- if (numblocks > 1L) M else 1L

  K <- numeric(n * B * M_eff)

  if (kernel == "Gaussian") {
    x <- x / (stats::sd(as.vector(x)) + 1e-19)
  }

  st <- 1L
  ed <- B * B
  index <- seq_len(n)

  for (m in 0:(M_eff - 1L)) {

    set.seed(m)
    index <- sample(index, length(index), replace = FALSE)

    if (numblocks > 0L) {
      for (i in seq.int(1L, n - discarded, by = B)) {
        j <- min(n, i + B - 1L)
        idx <- index[i:j]

        if (kernel == "Gaussian") {
          k <- kernel_gaussian_py(
            x[, idx, drop = FALSE],
            x[, idx, drop = FALSE],
            sigma = sqrt(d)
          )
        } else {
          k <- kernel_delta_norm_py(
            x[, idx, drop = FALSE],
            x[, idx, drop = FALSE]
          )
        }

        k <- H %*% k %*% H

        k <- k / (norm(k, type = "F") + 1e-9)

        K[st:ed] <- as.vector(t(k))

        st <- st + B * B
        ed <- ed + B * B
      }
    }
  }

  list(
    K = K,
    B = B,
    M_eff = M_eff,
    numblocks = numblocks,
    discarded = discarded
  )
}

# Build the HSIC Lasso design matrix in parallel.
hsic_lasso_prepare_py_parallel <- function(X_dn, Y_qn,
                                           y_kernel = c("Gaussian", "Delta"),
                                           x_kernel = c("Gaussian", "Delta"),
                                           B = 20L, M = 3L,
                                           n_cores = 1L,
                                           show_progress = TRUE) {
  y_kernel <- match.arg(y_kernel)
  x_kernel <- match.arg(x_kernel)

  X_dn <- as.matrix(X_dn)
  Y_qn <- as.matrix(Y_qn)

  d <- nrow(X_dn)
  n <- ncol(X_dn)

  if (ncol(Y_qn) != n) {
    stop("X_dn 和 Y_qn 的样本数不一致。")
  }

  if (B <= 0L || B > n) B <- n
  discarded <- n %% B
  numblocks <- n %/% B
  M_eff <- if (numblocks > 1L) M else 1L

  if (discarded > 0L) {
    warning(sprintf(
      "B = %d 不能整除样本数 %d；将有 %d 个样本处于 block 余数部分。为完全规整，建议改用能整除样本数的 B（例如 10）。",
      B, n, discarded
    ))
  }

  message(sprintf(
    "Block HSIC Lasso 参数：B = %d, M = %d, numblocks = %d, discarded = %d",
    B, M_eff, numblocks, discarded
  ))

  L_obj <- compute_kernel_py(
    Y_qn,
    kernel = y_kernel,
    B = B,
    M = M_eff,
    discarded = discarded
  )
  L_vec <- matrix(L_obj$K, ncol = 1)

  total_len <- length(L_obj$K)
  feat_index <- seq_len(d)

  compute_one_feature <- function(k) {
    compute_kernel_py(
      matrix(X_dn[k, ], nrow = 1),
      kernel = x_kernel,
      B = B,
      M = M_eff,
      discarded = discarded
    )$K
  }

  if (n_cores <= 1L) {
    if (show_progress) {
      pb <- txtProgressBar(min = 0, max = d, style = 3)
      on.exit(close(pb), add = TRUE)
    }

    K_list <- vector("list", length = d)
    for (k in feat_index) {
      K_list[[k]] <- compute_one_feature(k)
      if (show_progress) setTxtProgressBar(pb, k)
    }

  } else {

    cl <- parallel::makeCluster(n_cores, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(
      cl,
      varlist = c(
        "X_dn",
        "x_kernel",
        "B",
        "M_eff",
        "discarded",
        "compute_kernel_py",
        "kernel_gaussian_py",
        "kernel_delta_norm_py"
      ),
      envir = environment()
    )

    if (show_progress) {
      K_list <- pbapply::pblapply(feat_index, compute_one_feature, cl = cl)
    } else {
      K_list <- parallel::parLapply(cl, feat_index, compute_one_feature)
    }
  }

  K_mat <- do.call(cbind, K_list)

  KtL <- crossprod(K_mat, L_vec)

  list(
    K = K_mat,
    KtL = KtL,
    L = L_vec,
    B = B,
    M_eff = M_eff,
    numblocks = numblocks,
    discarded = discarded,
    n_cores = n_cores
  )
}

# Fit the nonnegative LARS path used by HSIC Lasso.
nnlars_py_clone <- function(X, Xty, num_feat = 30L, max_neighbors = 10L) {
  X <- as.matrix(X)
  Xty <- as.matrix(Xty)

  n <- nrow(X)
  d <- ncol(X)

  A <- integer(0)

  A_neighbors <- list()
  A_neighbors_score <- list()

  beta <- matrix(0, nrow = d, ncol = 1)

  path <- matrix(0, nrow = d, ncol = max(4L * d, 1L))

  lam <- numeric(max(4L * d, 1L))

  Iset <- seq_len(d)

  XtXbeta <- crossprod(X, X %*% beta)
  cvec <- as.numeric(Xty - XtXbeta)

  j_local <- which.max(cvec)
  C <- cvec[j_local]

  A <- c(A, Iset[j_local])
  Iset <- Iset[-j_local]

  lam[1] <- if (length(C) == 0) 0 else C
  k <- 0L

  repeat {

    if (length(A) == 0L) break

    if (!is.finite(mean(cvec[A])) || mean(cvec[A]) < 1e-9) break

    if (length(A) >= (num_feat + 1L)) break

    s <- matrix(1, nrow = length(A), ncol = 1)

    XA <- X[, A, drop = FALSE]
    GAA <- crossprod(XA)

    w <- tryCatch(
      solve(GAA, s),
      error = function(e) {
        XA2 <- XA + matrix(
          rnorm(length(XA), mean = 0, sd = 1e-9),
          nrow = nrow(XA), ncol = ncol(XA)
        )
        solve(crossprod(XA2), s)
      }
    )

    XtXw <- crossprod(X, XA %*% w)

    gamma1 <- if (length(Iset) > 0L) {
      (C - cvec[Iset]) / (XtXw[A[1], 1] - XtXw[Iset, 1])
    } else {
      numeric(0)
    }

    gamma2 <- -beta[A, 1] / as.numeric(w)

    gamma3 <- cvec[A[1]] / XtXw[A[1], 1]

    gamma <- c(gamma1, gamma2, gamma3)

    gamma[gamma <= 1e-9] <- Inf

    t_idx <- which.min(gamma)
    mu <- gamma[t_idx]

    if (!is.finite(mu)) break

    beta[A, 1] <- beta[A, 1] + mu * as.numeric(w)

    if (t_idx > length(gamma1) && t_idx < (length(gamma1) + length(gamma2) + 1L)) {
      lasso_cond <- 1L
      j_remove <- t_idx - length(gamma1)
      Iset <- c(Iset, A[j_remove])
      A <- A[-j_remove]
    } else {
      lasso_cond <- 0L
    }

    XtXbeta <- crossprod(X, X %*% beta)
    cvec <- as.numeric(Xty - XtXbeta)

    if (length(Iset) > 0L) {
      j_local <- which.max(cvec[Iset])
      C <- max(cvec[Iset])
    } else {
      j_local <- integer(0)
      C <- numeric(0)
    }

    k <- k + 1L
    path[, k + 1L] <- beta[, 1]
    lam[k + 1L] <- if (length(C) == 0) 0 else C

    if (lasso_cond == 0L && length(Iset) > 0L) {
      A <- c(A, Iset[j_local])
      Iset <- Iset[-j_local]
    }
  }

  if (length(A) > num_feat) {
    A <- A[-length(A)]
  }

  if (length(A) == 0L) {
    return(list(
      path = path[, seq_len(k + 1L), drop = FALSE],
      beta = beta,
      A = integer(0),
      lam = lam[seq_len(k + 1L)],
      A_neighbors = list(),
      A_neighbors_score = list()
    ))
  }

  svals <- beta[A, 1]
  sort_index <- order(svals, decreasing = TRUE)
  A_sorted <- A[sort_index]

  XtXA <- crossprod(X, X[, A_sorted, drop = FALSE])
  num_neighbors <- max_neighbors + 1L

  for (i in seq_along(A_sorted)) {
    tmp <- XtXA[, i]
    ord <- order(tmp, decreasing = TRUE)
    keep <- ord[seq_len(min(num_neighbors, length(ord)))]
    A_neighbors[[i]] <- keep
    A_neighbors_score[[i]] <- tmp[keep]
  }

  list(
    path = path[, seq_len(k + 1L), drop = FALSE],
    beta = beta,
    A = A_sorted,
    lam = lam[seq_len(k + 1L)],
    A_neighbors = A_neighbors,
    A_neighbors_score = A_neighbors_score
  )
}

# Run HSIC Lasso and return ranked features and coefficients.
run_hsiclasso_regression_py_clone_parallel <- function(X_nxp, Y_nxq,
                                                       featname = NULL,
                                                       num_feat = 30L,
                                                       B = 10L,
                                                       M = 3L,
                                                       max_neighbors = 10L,
                                                       n_cores = 1L,
                                                       show_progress = TRUE) {
  X_nxp <- as_numeric_matrix(X_nxp, "X")
  Y_nxq <- as_numeric_matrix(Y_nxq, "Y")

  if (nrow(X_nxp) != nrow(Y_nxq)) {
    stop("X 和 Y 的样本数不一致。")
  }

  if (is.null(featname)) {
    featname <- colnames(X_nxp)
  }
  if (is.null(featname)) {
    featname <- paste0("V", seq_len(ncol(X_nxp)))
  }

  prep <- hsic_lasso_prepare_py_parallel(
    X_dn = t(X_nxp),
    Y_qn = t(Y_nxq),
    y_kernel = "Gaussian",
    x_kernel = "Gaussian",
    B = B,
    M = M,
    n_cores = n_cores,
    show_progress = show_progress
  )

  scale_fac <- prep$numblocks * prep$M_eff
  if (scale_fac <= 0) {
    stop("无有效 block，请检查 B。")
  }

  X_scaled <- prep$K * sqrt(1 / scale_fac)
  Xty_scaled <- prep$KtL / scale_fac

  fit <- nnlars_py_clone(
    X = X_scaled,
    Xty = Xty_scaled,
    num_feat = num_feat,
    max_neighbors = max_neighbors
  )

  selected_idx <- fit$A
  selected_features <- featname[selected_idx]
  selected_scores <- as.numeric(fit$beta[selected_idx, 1])

  list(
    selected_indices = selected_idx,
    selected_features = selected_features,
    selected_scores = selected_scores,
    beta = fit$beta,
    path = fit$path,
    lam = fit$lam,
    neighbors = fit$A_neighbors,
    neighbors_score = fit$A_neighbors_score,
    prep = prep
  )
}

# Format selected variables for export.
make_export_df <- function(selected_features, adni_noid_colnames) {
  data.frame(
    important_variable = selected_features,
    column_in_ADNI_1_3_COMBINE1215_after_drop_id = match(selected_features, adni_noid_colnames),
    stringsAsFactors = FALSE
  )
}
