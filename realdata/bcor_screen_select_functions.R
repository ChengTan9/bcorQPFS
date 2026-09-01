# Utilities for BCor screening and bcorQPFS variable selection.

# Manage packages, files, progress reporting, matrices, and parallel workers.
get_default_cores <- function() {
  k <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (!is.finite(k) || is.na(k) || k < 1L) return(1L)
  max(1L, as.integer(k) - 3L)
}

check_required_packages <- function() {
  pkgs <- c("Ball", "hdf5r", "openxlsx", "CVXR", "progress")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    stop(
      sprintf("缺少以下 R 包，请先安装：%s", paste(miss, collapse = ", ")),
      call. = FALSE
    )
  }
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

save_rdata_named <- function(object, file, object_name = "obj") {
  ensure_dir(dirname(file))
  e <- new.env(parent = emptyenv())
  assign(object_name, object, envir = e)
  save(list = object_name, file = file, envir = e)
  invisible(file)
}

load_first_object_from_rdata <- function(file, object_name = NULL) {
  if (!file.exists(file)) {
    stop(sprintf("文件不存在：%s", file), call. = FALSE)
  }
  e <- new.env(parent = emptyenv())
  nm <- load(file, envir = e)
  if (length(nm) == 0L) {
    stop(sprintf("RData 文件为空：%s", file), call. = FALSE)
  }
  if (!is.null(object_name)) {
    if (!object_name %in% nm) {
      stop(sprintf("对象 '%s' 不在文件 %s 中。", object_name, file), call. = FALSE)
    }
    return(get(object_name, envir = e))
  }
  get(nm[1], envir = e)
}

make_eta_pb <- function(total, label, show_progress = TRUE) {
  if (!show_progress || total <= 0L) return(NULL)
  pb <- progress::progress_bar$new(
    total = total,
    clear = FALSE,
    width = 100,
    show_after = 0,
    format = paste0(label, " [:bar] :current/:total (:percent) elapsed=:elapsed eta=:eta")
  )
  pb$tick(0)
  pb
}

tick_pb <- function(pb, n = 1L) {
  if (!is.null(pb)) pb$tick(n)
  invisible(NULL)
}

split_into_blocks <- function(x, block_size) {
  if (length(x) == 0L) return(list())
  split(x, ceiling(seq_along(x) / block_size))
}

safe_matrix_inverse <- function(A, tol = 1e-8) {
  A <- as.matrix(A)
  ed <- eigen(A, symmetric = TRUE)
  val <- Re(ed$values)
  vec <- Re(ed$vectors)
  inv_val <- ifelse(val > tol, 1 / val, 0)
  vec %*% diag(inv_val, nrow = length(inv_val), ncol = length(inv_val)) %*% t(vec)
}

open_parallel_cluster <- function(n_cores, packages = character(0)) {
  if (n_cores <= 1L) return(NULL)

  cl <- parallel::makeCluster(n_cores)

  if (length(packages) > 0L) {
    parallel::clusterCall(
      cl,
      function(pkgs) {
        invisible(lapply(pkgs, require, character.only = TRUE))
        NULL
      },
      pkgs = packages
    )
  }

  cl
}

# Convert input data to finite numeric objects.
sanitize_vector <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) x <- as.integer(x)
  suppressWarnings(as.numeric(x))
}

sanitize_df_like <- function(x) {
  if (is.vector(x) && !is.list(x)) {
    out <- data.frame(V1 = sanitize_vector(x), check.names = FALSE)
    return(out)
  }
  if (is.matrix(x)) x <- as.data.frame(x, check.names = FALSE)
  if (!is.data.frame(x)) {
    stop("输入必须是 vector / matrix / data.frame。", call. = FALSE)
  }
  out <- as.data.frame(lapply(x, sanitize_vector), check.names = FALSE, stringsAsFactors = FALSE)
  rownames(out) <- rownames(x)
  colnames(out) <- colnames(x)
  out
}

as_numeric_df <- function(x) {
  out <- sanitize_df_like(x)
  if (anyNA(out)) {
    stop("数据中存在无法转换为数值的元素，转换后出现 NA。", call. = FALSE)
  }
  out
}

is_constant_vector <- function(x) {
  x <- sanitize_vector(x)
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(TRUE)
  length(unique(x)) <= 1L
}

# Compute Ball-correlation quantities.
compute_distance_matrix <- function(Y) {
  if (is.vector(Y) && !is.list(Y)) {
    Y <- matrix(sanitize_vector(Y), ncol = 1L)
  } else {
    Y <- as.matrix(as_numeric_df(Y))
  }

  if (anyNA(Y)) {
    stop("Y 转为数值后存在 NA，无法计算距离矩阵。", call. = FALSE)
  }

  as.matrix(stats::dist(Y))
}

compute_single_bcor <- function(x, DY) {
  x <- sanitize_vector(x)
  if (anyNA(x)) stop("x 转为数值后存在 NA。", call. = FALSE)
  if (is_constant_vector(x)) return(0)

  DX <- as.matrix(stats::dist(matrix(x, ncol = 1L)))

  out <- tryCatch(
    suppressWarnings(Ball::bcor(DX, DY, distance = TRUE)),
    error = function(e) NA_real_
  )
  as.numeric(out)
}

compute_pair_bcor_xx <- function(x1, x2) {
  x1 <- sanitize_vector(x1)
  x2 <- sanitize_vector(x2)

  if (anyNA(x1) || anyNA(x2)) {
    stop("x1 或 x2 转为数值后存在 NA。", call. = FALSE)
  }
  if (is_constant_vector(x1) || is_constant_vector(x2)) return(0)

  DX1 <- as.matrix(stats::dist(matrix(x1, ncol = 1L)))
  DX2 <- as.matrix(stats::dist(matrix(x2, ncol = 1L)))

  out <- tryCatch(
    suppressWarnings(Ball::bcor(DX1, DX2, distance = TRUE)),
    error = function(e) NA_real_
  )
  as.numeric(out)
}

# Estimate memory use and manage optional HDF5 caching.
estimate_screening_memory <- function(
    X,
    n_cores = get_default_cores(),
    memory_limit_gb = 8,
    reserve_gb = 1.5) {
  X <- as_numeric_df(X)
  n <- nrow(X)
  p <- ncol(X)

  bytes_X <- n * p * 8
  bytes_DY <- n * n * 8
  bytes_parallel_overhead <- n_cores * (bytes_DY + min(bytes_X, 5e7))
  bytes_misc <- 5e8

  total_bytes <- bytes_X + bytes_DY + bytes_parallel_overhead + bytes_misc
  total_gb <- total_bytes / 1024^3
  threshold_gb <- max(0.5, memory_limit_gb - reserve_gb)

  list(
    n = n,
    p = p,
    n_cores = n_cores,
    estimated_total_bytes = total_bytes,
    estimated_total_gb = total_gb,
    threshold_gb = threshold_gb,
    use_h5 = isTRUE(total_gb > threshold_gb)
  )
}

cache_matrix_to_h5 <- function(
    X,
    h5_path,
    dataset_name = "X",
    show_progress = TRUE) {
  check_required_packages()
  X <- as_numeric_df(X)

  n <- nrow(X)
  p <- ncol(X)

  ensure_dir(dirname(h5_path))
  if (file.exists(h5_path)) file.remove(h5_path)

  rn <- rownames(X)
  cn <- colnames(X)
  if (is.null(rn)) rn <- paste0("row_", seq_len(n))
  if (is.null(cn)) cn <- paste0("V", seq_len(p))

  h5f <- hdf5r::H5File$new(h5_path, mode = "w")
  on.exit(try(h5f$close_all(), silent = TRUE), add = TRUE)

  dset <- h5f$create_dataset(
    name = dataset_name,
    dtype = hdf5r::h5types$H5T_NATIVE_DOUBLE,
    dims = c(n, p),
    chunk_dims = c(min(n, 256L), 1L)
  )

  h5f$create_dataset(name = paste0(dataset_name, "_rownames"), robj = rn)
  h5f$create_dataset(name = paste0(dataset_name, "_colnames"), robj = cn)

  pb <- make_eta_pb(total = p, label = "screening-h5-write", show_progress = show_progress)

  for (j in seq_len(p)) {
    dset[, j] <- as.numeric(X[[j]])
    tick_pb(pb, 1L)
  }

  invisible(h5_path)
}

read_h5_column <- function(h5_path, j, dataset_name = "X") {
  if (!file.exists(h5_path)) {
    stop(sprintf("h5 文件不存在：%s", h5_path), call. = FALSE)
  }
  h5f <- hdf5r::H5File$new(h5_path, mode = "r")
  on.exit(try(h5f$close_all(), silent = TRUE), add = TRUE)
  as.numeric(h5f[[dataset_name]][, j])
}

read_h5_dimnames <- function(h5_path, dataset_name = "X") {
  if (!file.exists(h5_path)) {
    stop(sprintf("h5 文件不存在：%s", h5_path), call. = FALSE)
  }
  h5f <- hdf5r::H5File$new(h5_path, mode = "r")
  on.exit(try(h5f$close_all(), silent = TRUE), add = TRUE)

  rn <- tryCatch(h5f[[paste0(dataset_name, "_rownames")]][], error = function(e) NULL)
  cn <- tryCatch(h5f[[paste0(dataset_name, "_colnames")]][], error = function(e) NULL)

  list(rownames = rn, colnames = cn)
}

can_reuse_h5_cache <- function(
    X,
    h5_path,
    dataset_name = "X") {

  if (is.null(h5_path) || !nzchar(h5_path) || !file.exists(h5_path)) {
    return(FALSE)
  }

  X <- as_numeric_df(X)
  n <- nrow(X)
  p <- ncol(X)
  x_colnames <- colnames(X)

  ok <- tryCatch({
    h5f <- hdf5r::H5File$new(h5_path, mode = "r")
    on.exit(try(h5f$close_all(), silent = TRUE), add = TRUE)

    if (!(dataset_name %in% names(h5f))) {
      return(FALSE)
    }

    dset <- h5f[[dataset_name]]
    dims <- dset$dims

    if (length(dims) != 2L) return(FALSE)
    if (!identical(as.integer(dims), as.integer(c(n, p)))) return(FALSE)

    cn_name <- paste0(dataset_name, "_colnames")
    if (cn_name %in% names(h5f)) {
      h5_colnames <- tryCatch(h5f[[cn_name]][], error = function(e) NULL)
      if (!is.null(h5_colnames) && !is.null(x_colnames)) {
        if (!identical(as.character(h5_colnames), as.character(x_colnames))) {
          return(FALSE)
        }
      }
    }

    TRUE
  }, error = function(e) FALSE)

  ok
}

# Evaluate BCor scores in parallel from memory or HDF5.
parallel_bcor_from_matrix <- function(
    X,
    DY,
    indices,
    cl = NULL,
    show_progress = TRUE,
    label = "bcor-matrix",
    block_size = 16L) {

  indices <- as.integer(indices)
  n_task <- length(indices)
  out <- numeric(n_task)

  if (n_task == 0L) return(out)

  pb <- make_eta_pb(total = n_task, label = label, show_progress = show_progress)

  if (is.null(cl)) {
    for (ii in seq_along(indices)) {
      out[ii] <- compute_single_bcor(X[[indices[ii]]], DY)
      tick_pb(pb, 1L)
    }
    return(out)
  }

  parallel::clusterExport(
    cl,
    varlist = c("X", "DY", "compute_single_bcor", "sanitize_vector", "is_constant_vector"),
    envir = environment()
  )

  blocks <- split_into_blocks(seq_along(indices), block_size = block_size)

  pos <- 1L
  for (bb in seq(1L, length(blocks), by = length(cl))) {
    block_ids <- bb:min(length(blocks), bb + length(cl) - 1L)
    this_blocks <- blocks[block_ids]
    this_idx <- lapply(this_blocks, function(id_sub) indices[id_sub])

    res <- parallel::parLapplyLB(
      cl,
      this_idx,
      function(idx) {
        vals <- numeric(length(idx))
        for (u in seq_along(idx)) {
          vals[u] <- compute_single_bcor(X[[idx[u]]], DY)
        }
        vals
      }
    )

    for (rr in seq_along(res)) {
      len <- length(res[[rr]])
      out[pos:(pos + len - 1L)] <- res[[rr]]
      pos <- pos + len
      tick_pb(pb, len)
    }
  }

  out
}

parallel_bcor_from_h5 <- function(
    h5_path,
    dataset_name,
    DY,
    indices,
    cl = NULL,
    show_progress = TRUE,
    label = "bcor-h5",
    block_size = 256L) {

  indices <- as.integer(indices)
  n_task <- length(indices)
  out <- numeric(n_task)

  if (n_task == 0L) return(out)

  pb <- make_eta_pb(total = n_task, label = label, show_progress = show_progress)

  if (is.null(cl)) {
    h5f <- hdf5r::H5File$new(h5_path, mode = "r")
    on.exit(try(h5f$close_all(), silent = TRUE), add = TRUE)

    dset <- h5f[[dataset_name]]

    for (ii in seq_along(indices)) {
      xj <- as.numeric(dset[, indices[ii]])
      out[ii] <- compute_single_bcor(xj, DY)
      tick_pb(pb, 1L)
    }
    return(out)
  }

  parallel::clusterExport(
    cl,
    varlist = c("h5_path", "dataset_name", "DY", "compute_single_bcor",
                "sanitize_vector", "is_constant_vector"),
    envir = environment()
  )
  parallel::clusterEvalQ(cl, {
    library(hdf5r)
    NULL
  })

  index_blocks <- split_into_blocks(indices, block_size = block_size)

  pos <- 1L
  for (bb in seq(1L, length(index_blocks), by = length(cl))) {
    block_ids <- bb:min(length(index_blocks), bb + length(cl) - 1L)
    task_list <- index_blocks[block_ids]

    res <- parallel::parLapplyLB(
      cl,
      task_list,
      function(idx_block) {
        h5f <- hdf5r::H5File$new(h5_path, mode = "r")
        on.exit(try(h5f$close_all(), silent = TRUE), add = TRUE)

        dset <- h5f[[dataset_name]]
        vals <- numeric(length(idx_block))

        for (u in seq_along(idx_block)) {
          xj <- as.numeric(dset[, idx_block[u]])
          vals[u] <- compute_single_bcor(xj, DY)
        }

        vals
      }
    )

    for (rr in seq_along(res)) {
      len <- length(res[[rr]])
      out[pos:(pos + len - 1L)] <- res[[rr]]
      pos <- pos + len
      tick_pb(pb, len)
    }
  }

  out
}

# Apply permutation-calibrated BCor screening.
screen_bcor <- function(
    X,
    Y,
    n_cores = get_default_cores(),
    save_dir = NULL,
    output_prefix = NULL,
    save_options = c("screened_only"),
    memory_limit_gb = 8,
    reserve_gb = 1.5,
    force_h5 = FALSE,
    keep_h5_cache = TRUE,
    show_progress = TRUE,
    seed = 123,
    dataset_name = "X") {

  check_required_packages()
  set.seed(seed)

  X <- as_numeric_df(X)
  n <- nrow(X)
  p <- ncol(X)

  Y_n <- if (is.vector(Y) && !is.list(Y)) length(Y) else nrow(as.data.frame(Y))
  if (n != Y_n) stop("X 与 Y 的样本量不一致。", call. = FALSE)

  if (!is.null(save_dir)) ensure_dir(save_dir)

  message(sprintf("[screen_bcor] 开始：n = %d, p = %d, n_cores = %d", n, p, n_cores))
  message("[screen_bcor] 第一步：计算 Y 的距离矩阵 DY（这一步本身仍为串行）...")
  DY <- compute_distance_matrix(Y)
  message("[screen_bcor] DY 计算完成。")

  mem_plan <- estimate_screening_memory(
    X = X,
    n_cores = n_cores,
    memory_limit_gb = memory_limit_gb,
    reserve_gb = reserve_gb
  )
  use_h5 <- isTRUE(force_h5) || isTRUE(mem_plan$use_h5)

  message(sprintf(
    "[screen_bcor] 预计总内存约 %.3f GB；是否启用 h5：%s",
    mem_plan$estimated_total_gb,
    ifelse(use_h5, "TRUE", "FALSE")
  ))

  sampled_idx <- sample.int(p, size = p, replace = TRUE)

  threshold_path <- NULL
  bcor_path <- NULL
  screened_path <- NULL
  h5_path <- NULL

  if (!is.null(output_prefix) && !is.null(save_dir)) {
    threshold_path <- file.path(save_dir, paste0(output_prefix, "_threshold.csv"))
    bcor_path <- file.path(save_dir, paste0(output_prefix, "_bcor.csv"))
    screened_path <- file.path(save_dir, paste0(output_prefix, "_selected.RData"))
    h5_path <- file.path(save_dir, paste0(output_prefix, "_cache.h5"))
  }

  if (use_h5 && is.null(h5_path)) {
    h5_path <- tempfile(pattern = "screen_cache_", fileext = ".h5")
  }

  cl <- open_parallel_cluster(n_cores = n_cores, packages = c("Ball", "hdf5r"))
  on.exit(if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  if (use_h5) {
    if (can_reuse_h5_cache(X = X, h5_path = h5_path, dataset_name = dataset_name)) {
      message(sprintf(
        "[screen_bcor] 第二步：检测到已有可复用 h5 文件，直接读取：%s",
        h5_path
      ))
    } else {
      message("[screen_bcor] 第二步：生成 h5 文件（ETA 进度条）...")
      cache_matrix_to_h5(
        X = X,
        h5_path = h5_path,
        dataset_name = dataset_name,
        show_progress = show_progress
      )
    }

    message("[screen_bcor] 第三步：并行计算阈值所需 sampled_bcor（ETA 进度条）...")
    sampled_bcor <- parallel_bcor_from_h5(
      h5_path = h5_path,
      dataset_name = dataset_name,
      DY = DY,
      indices = sampled_idx,
      cl = cl,
      show_progress = show_progress,
      label = "screening-threshold-h5"
    )

    tau_tilde <- max(sampled_bcor, na.rm = TRUE)
    tau_n <- (sqrt(p) / (2 * sqrt(p) + 10)) * tau_tilde

    message(sprintf("[screen_bcor] 阈值完成：tau_tilde = %.6f, tau_n = %.6f", tau_tilde, tau_n))

    message("[screen_bcor] 第四步：并行扫描全部变量 BCor（ETA 进度条）...")
    scan_bcor <- parallel_bcor_from_h5(
      h5_path = h5_path,
      dataset_name = dataset_name,
      DY = DY,
      indices = seq_len(p),
      cl = cl,
      show_progress = show_progress,
      label = "screening-scan-h5"
    )

  } else {
    message("[screen_bcor] 第二步：并行计算阈值所需 sampled_bcor（ETA 进度条）...")
    sampled_bcor <- parallel_bcor_from_matrix(
      X = X,
      DY = DY,
      indices = sampled_idx,
      cl = cl,
      show_progress = show_progress,
      label = "screening-threshold"
    )

    tau_tilde <- max(sampled_bcor, na.rm = TRUE)
    tau_n <- (sqrt(p) / (2 * sqrt(p) + 10)) * tau_tilde

    message(sprintf("[screen_bcor] 阈值完成：tau_tilde = %.6f, tau_n = %.6f", tau_tilde, tau_n))

    message("[screen_bcor] 第三步：并行扫描全部变量 BCor（ETA 进度条）...")
    scan_bcor <- parallel_bcor_from_matrix(
      X = X,
      DY = DY,
      indices = seq_len(p),
      cl = cl,
      show_progress = show_progress,
      label = "screening-scan"
    )
  }

  keep_flag <- is.finite(scan_bcor) & (scan_bcor > tau_n)
  selected_names <- colnames(X)[keep_flag]
  X_screened <- X[, keep_flag, drop = FALSE]

  message(sprintf("[screen_bcor] 完成：保留变量 %d / %d", ncol(X_screened), p))

  if (!is.null(save_dir) && !is.null(output_prefix)) {
    if ("threshold" %in% save_options) {
      utils::write.csv(
        data.frame(tau_tilde = tau_tilde, tau_n = tau_n, p = p, stringsAsFactors = FALSE),
        file = threshold_path,
        row.names = FALSE
      )
    }

    if ("bcor" %in% save_options) {
      utils::write.csv(
        data.frame(
          variable = colnames(X),
          bcor = scan_bcor,
          selected = keep_flag,
          stringsAsFactors = FALSE
        ),
        file = bcor_path,
        row.names = FALSE
      )
    }

    if ("screened_only" %in% save_options) {
      save_rdata_named(X_screened, file = screened_path, object_name = "selected_X")
    }
  }

  if (use_h5 && !keep_h5_cache && !is.null(h5_path) && file.exists(h5_path)) {
    file.remove(h5_path)
    h5_path <- NULL
  }

  list(
    X_screened = X_screened,
    selected_names = selected_names,
    keep_flag = keep_flag,
    bcor = scan_bcor,
    tau_tilde = tau_tilde,
    tau_n = tau_n,
    DY = DY,
    used_h5 = use_h5,
    h5_path = h5_path,
    memory_plan = mem_plan,
    saved_files = list(
      threshold_csv = threshold_path,
      bcor_csv = bcor_path,
      screened_rdata = screened_path
    )
  )
}

# Construct relevance and redundancy quantities for variable selection.
compute_F_vector <- function(
    X,
    DY = NULL,
    Y = NULL,
    cl = NULL,
    show_progress = TRUE) {

  X <- as_numeric_df(X)
  if (is.null(DY)) {
    if (is.null(Y)) stop("DY 与 Y 不能同时为空。", call. = FALSE)
    DY <- compute_distance_matrix(Y)
  }

  F_vec <- parallel_bcor_from_matrix(
    X = X,
    DY = DY,
    indices = seq_len(ncol(X)),
    cl = cl,
    show_progress = show_progress,
    label = "selection-F"
  )
  names(F_vec) <- colnames(X)
  F_vec
}

parallel_pair_bcor_matrix <- function(
    X,
    pair_mat,
    cl = NULL,
    show_progress = TRUE,
    label = "pair-bcor",
    block_size = 16L) {

  n_task <- nrow(pair_mat)
  if (n_task == 0L) {
    return(data.frame(i = integer(0), j = integer(0), val = numeric(0)))
  }

  pb <- make_eta_pb(total = n_task, label = label, show_progress = show_progress)

  if (is.null(cl)) {
    out_i <- integer(n_task)
    out_j <- integer(n_task)
    out_v <- numeric(n_task)

    for (tt in seq_len(n_task)) {
      i <- pair_mat[tt, 1]
      j <- pair_mat[tt, 2]
      out_i[tt] <- i
      out_j[tt] <- j
      out_v[tt] <- compute_pair_bcor_xx(X[[i]], X[[j]])
      tick_pb(pb, 1L)
    }

    return(data.frame(i = out_i, j = out_j, val = out_v))
  }

  parallel::clusterExport(
    cl,
    varlist = c("X", "compute_pair_bcor_xx", "sanitize_vector", "is_constant_vector"),
    envir = environment()
  )

  row_blocks <- split_into_blocks(seq_len(n_task), block_size = block_size)
  out_list <- vector("list", length(row_blocks))

  for (bb in seq(1L, length(row_blocks), by = length(cl))) {
    block_ids <- bb:min(length(row_blocks), bb + length(cl) - 1L)

    task_list <- lapply(block_ids, function(idb) pair_mat[row_blocks[[idb]], , drop = FALSE])

    res <- parallel::parLapplyLB(
      cl,
      task_list,
      function(pair_block) {
        ii <- pair_block[, 1]
        jj <- pair_block[, 2]
        vv <- numeric(length(ii))

        for (u in seq_along(ii)) {
          vv[u] <- compute_pair_bcor_xx(X[[ii[u]]], X[[jj[u]]])
        }

        data.frame(i = ii, j = jj, val = vv)
      }
    )

    for (kk in seq_along(res)) {
      out_list[[block_ids[kk]]] <- res[[kk]]
      tick_pb(pb, nrow(res[[kk]]))
    }
  }

  do.call(rbind, out_list)
}

compute_Q_matrix_exact <- function(
    X,
    cl = NULL,
    show_progress = TRUE) {

  X <- as_numeric_df(X)
  p <- ncol(X)
  Q <- diag(1, p)
  rownames(Q) <- colnames(X)
  colnames(Q) <- colnames(X)

  if (p <= 1L) return(Q)

  pair_mat <- utils::combn(p, 2L)
  pair_mat <- t(pair_mat)

  pair_res <- parallel_pair_bcor_matrix(
    X = X,
    pair_mat = pair_mat,
    cl = cl,
    show_progress = show_progress,
    label = "selection-Q-exact"
  )

  for (r in seq_len(nrow(pair_res))) {
    i <- pair_res$i[r]
    j <- pair_res$j[r]
    v <- pair_res$val[r]
    if (!is.finite(v)) v <- 0
    Q[i, j] <- v
    Q[j, i] <- v
  }

  diag(Q) <- 1
  Q
}

compute_C_from_AB_parallel <- function(
    A_inv,
    B,
    cl = NULL,
    show_progress = TRUE,
    label = "selection-nystrom-C",
    block_size = 16L) {

  q <- ncol(B)
  if (q == 0L) return(matrix(0, 0, 0))

  pb <- make_eta_pb(total = q, label = label, show_progress = show_progress)

  if (is.null(cl)) {
    C <- matrix(0, q, q)
    for (st in seq(1L, q, by = block_size)) {
      ed <- min(q, st + block_size - 1L)
      rows <- st:ed
      C[rows, ] <- t(B[, rows, drop = FALSE]) %*% A_inv %*% B
      tick_pb(pb, length(rows))
    }
    return(C)
  }

  parallel::clusterExport(
    cl,
    varlist = c("A_inv", "B"),
    envir = environment()
  )

  row_blocks <- split_into_blocks(seq_len(q), block_size = block_size)
  res_list <- vector("list", length(row_blocks))

  for (bb in seq(1L, length(row_blocks), by = length(cl))) {
    block_ids <- bb:min(length(row_blocks), bb + length(cl) - 1L)
    this_blocks <- row_blocks[block_ids]

    res <- parallel::parLapplyLB(
      cl,
      this_blocks,
      function(rows) {
        subC <- t(B[, rows, drop = FALSE]) %*% A_inv %*% B
        list(rows = rows, subC = subC)
      }
    )

    for (kk in seq_along(res)) {
      res_list[[block_ids[kk]]] <- res[[kk]]
      tick_pb(pb, length(res[[kk]]$rows))
    }
  }

  C <- matrix(0, q, q)
  for (res in res_list) {
    C[res$rows, ] <- res$subC
  }
  C
}

compute_Q_matrix_nystrom <- function(
    X,
    k = NULL,
    seed = 123,
    cl = NULL,
    show_progress = TRUE) {

  X <- as_numeric_df(X)
  p <- ncol(X)
  cn <- colnames(X)

  if (is.null(k)) k <- max(2L, floor(0.3 * p))
  k <- max(2L, min(as.integer(k), p))

  set.seed(seed)
  anchor_idx <- sort(sample.int(p, size = k, replace = FALSE))
  nonanchor_idx <- setdiff(seq_len(p), anchor_idx)

  A <- diag(1, k)
  rownames(A) <- cn[anchor_idx]
  colnames(A) <- cn[anchor_idx]

  if (k >= 2L) {
    A_pairs <- utils::combn(anchor_idx, 2L)
    A_pairs <- t(A_pairs)

    A_res <- parallel_pair_bcor_matrix(
      X = X,
      pair_mat = A_pairs,
      cl = cl,
      show_progress = show_progress,
      label = "selection-nystrom-A"
    )

    pos_map_a <- match(anchor_idx, anchor_idx)
    names(pos_map_a) <- as.character(anchor_idx)

    for (r in seq_len(nrow(A_res))) {
      i0 <- as.character(A_res$i[r])
      j0 <- as.character(A_res$j[r])
      ii <- pos_map_a[[i0]]
      jj <- pos_map_a[[j0]]
      v <- A_res$val[r]
      if (!is.finite(v)) v <- 0
      A[ii, jj] <- v
      A[jj, ii] <- v
    }
  }

  q <- length(nonanchor_idx)
  if (q > 0L) {
    B_pairs <- expand.grid(i = anchor_idx, j = nonanchor_idx, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    B_pairs <- as.matrix(B_pairs)

    B_res <- parallel_pair_bcor_matrix(
      X = X,
      pair_mat = B_pairs,
      cl = cl,
      show_progress = show_progress,
      label = "selection-nystrom-B"
    )

    B <- matrix(0, nrow = k, ncol = q)
    rownames(B) <- cn[anchor_idx]
    colnames(B) <- cn[nonanchor_idx]

    pos_map_a <- seq_along(anchor_idx); names(pos_map_a) <- as.character(anchor_idx)
    pos_map_b <- seq_along(nonanchor_idx); names(pos_map_b) <- as.character(nonanchor_idx)

    for (r in seq_len(nrow(B_res))) {
      ii <- pos_map_a[[as.character(B_res$i[r])]]
      jj <- pos_map_b[[as.character(B_res$j[r])]]
      v <- B_res$val[r]
      if (!is.finite(v)) v <- 0
      B[ii, jj] <- v
    }
  } else {
    B <- matrix(0, nrow = k, ncol = 0L)
  }

  ridge <- 1e-8
  A_inv <- safe_matrix_inverse(A + diag(ridge, nrow(A)))
  C_approx <- compute_C_from_AB_parallel(
    A_inv = A_inv,
    B = B,
    cl = cl,
    show_progress = show_progress,
    label = "selection-nystrom-C"
  )

  Q <- matrix(0, nrow = p, ncol = p)
  rownames(Q) <- cn
  colnames(Q) <- cn

  Q[anchor_idx, anchor_idx] <- A
  if (length(nonanchor_idx) > 0L) {
    Q[anchor_idx, nonanchor_idx] <- B
    Q[nonanchor_idx, anchor_idx] <- t(B)
    Q[nonanchor_idx, nonanchor_idx] <- C_approx
  }
  diag(Q) <- 1

  list(
    Q = Q,
    anchor_idx = anchor_idx,
    nonanchor_idx = nonanchor_idx,
    A = A,
    B = B
  )
}

# Solve the quadratic program and return selected variables.
eigen_decompose_positive <- function(Q, tol = 1e-8) {
  ed <- eigen(Q, symmetric = TRUE)
  lambda <- Re(ed$values)
  U <- Re(ed$vectors)
  keep <- which(lambda > tol)

  if (length(keep) == 0L) {
    return(list(
      lambda = numeric(0),
      U = matrix(0, nrow = nrow(Q), ncol = 0L),
      keep_idx = integer(0),
      full_eigen = ed
    ))
  }

  list(
    lambda = lambda[keep],
    U = U[, keep, drop = FALSE],
    keep_idx = keep,
    full_eigen = ed
  )
}

solve_selection_qp <- function(
    F,
    U,
    lambda,
    alpha = 0.5,
    solver = NULL) {

  check_required_packages()

  F <- as.numeric(F)
  lambda <- as.numeric(lambda)

  if (length(lambda) != ncol(U)) stop("lambda 长度必须等于 U 的列数。", call. = FALSE)
  if (length(F) != nrow(U)) stop("F 长度必须等于 U 的行数。", call. = FALSE)
  if (length(lambda) == 0L) stop("没有正特征值，无法继续优化。", call. = FALSE)

  linear_coef <- alpha * as.numeric(crossprod(U, F))

  y <- CVXR::Variable(length(lambda))
  Lambda_mat <- diag(lambda, nrow = length(lambda), ncol = length(lambda))

  objective <- CVXR::Minimize(
    0.5 * (1 - alpha) * CVXR::quad_form(y, Lambda_mat) - t(linear_coef) %*% y
  )
  constraints <- list(U %*% y >= 0)
  problem <- CVXR::Problem(objective, constraints)

  opt_value <- if (is.null(solver)) CVXR::psolve(problem) else CVXR::psolve(problem, solver = solver)
  solve_status <- CVXR::status(problem)

  if (!(solve_status %in% c("optimal", "optimal_inaccurate"))) {
    stop(sprintf("CVXR 求解失败，状态为: %s", solve_status), call. = FALSE)
  }

  y_hat <- as.numeric(CVXR::value(y))
  x_hat <- as.numeric(U %*% y_hat)

  list(
    y = y_hat,
    x = x_hat,
    status = solve_status,
    value = opt_value
  )
}

select_variables_bcor <- function(
    X,
    Y = NULL,
    DY = NULL,
    M = 30,
    alpha = 0.5,
    use_nystrom = FALSE,
    k = NULL,
    n_cores = get_default_cores(),
    seed = 123,
    solver = NULL,
    show_progress = TRUE) {

  X <- as_numeric_df(X)
  p <- ncol(X)
  M <- min(as.integer(M), p)

  if (is.null(DY)) {
    if (is.null(Y)) stop("DY 与 Y 不能同时为空。", call. = FALSE)
    DY <- compute_distance_matrix(Y)
  }

  cl <- open_parallel_cluster(n_cores = n_cores, packages = c("Ball"))
  on.exit(if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  message("[select_variables_bcor] 第一步：并行计算 F（ETA 进度条）...")
  F_vec <- compute_F_vector(
    X = X,
    DY = DY,
    cl = cl,
    show_progress = show_progress
  )

  if (use_nystrom) {
    message("[select_variables_bcor] 第二步：并行计算 Nystrom 的 A / B / C（ETA 进度条）...")
    Q_obj <- compute_Q_matrix_nystrom(
      X = X,
      k = k,
      seed = seed,
      cl = cl,
      show_progress = show_progress
    )
    Q <- Q_obj$Q
  } else {
    message("[select_variables_bcor] 第二步：并行计算精确 Q（ETA 进度条）...")
    Q_obj <- NULL
    Q <- compute_Q_matrix_exact(
      X = X,
      cl = cl,
      show_progress = show_progress
    )
  }

  eig_obj <- eigen_decompose_positive(Q)
  lambda_pos <- eig_obj$lambda
  U_pos <- eig_obj$U

  if (length(lambda_pos) == 0L || ncol(U_pos) == 0L) {
    score <- as.numeric(F_vec)
    names(score) <- colnames(X)
    ord <- order(score, decreasing = TRUE)
    selected_variables <- colnames(X)[ord][seq_len(M)]

    return(list(
      selected_variables = selected_variables,
      score = score,
      rank_order = names(score)[ord],
      F = F_vec,
      Q = Q,
      eigen = eig_obj,
      qp_result = NULL,
      use_nystrom = use_nystrom,
      nystrom_detail = Q_obj
    ))
  }

  message("[select_variables_bcor] 第三步：CVXR 求解选择优化问题...")
  qp_res <- solve_selection_qp(
    F = F_vec,
    U = U_pos,
    lambda = lambda_pos,
    alpha = alpha,
    solver = solver
  )

  score <- as.numeric(qp_res$x)
  names(score) <- colnames(X)

  ord <- order(score, decreasing = TRUE)
  selected_variables <- colnames(X)[ord][seq_len(M)]

  list(
    selected_variables = selected_variables,
    score = score,
    rank_order = names(score)[ord],
    F = F_vec,
    Q = Q,
    eigen = eig_obj,
    qp_result = qp_res,
    use_nystrom = use_nystrom,
    nystrom_detail = Q_obj
  )
}

# Generate simulation data and summarize selection results.
simulate_example_data <- function(example_id, n, p, seed = NULL, rho = 0.8) {

  if (!is.null(seed)) set.seed(seed)

  n <- as.integer(n)
  p <- as.integer(p)
  rho <- as.numeric(rho)

  if (p < 5L) {
    stop("为了匹配 example 构造，p 至少需要 >= 5。", call. = FALSE)
  }
  if (!is.finite(rho) || abs(rho) >= 1) {
    stop("rho 需要满足 |rho| < 1。", call. = FALSE)
  }

  idx <- seq_len(p)
  Sigma <- rho ^ abs(outer(idx, idx, "-"))

  Z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  R <- chol(Sigma)
  X <- Z %*% R
  colnames(X) <- paste0("X", seq_len(p))

  eps <- stats::rnorm(n)

  if (example_id == 1L) {
    omega <- 3 * X[, 2]^2 - 3 * X[, 3]^2 + 5 * X[, 5] + eps
    Y <- 1 + 2 * as.numeric(omega > 0) * omega
    active_set <- c("X2", "X3", "X5")
  } else if (example_id == 2L) {
    omega <- 3 * X[, 1] * X[, 2] + 5 * X[, 5] + eps
    Y <- 1 + 2 * as.numeric(omega < 0) * omega
    active_set <- c("X1", "X2", "X5")
  } else if (example_id == 3L) {
    Y <- 2 * sin(3 * X[, 1] * X[, 3]) - exp(X[, 2]^2) + eps
    active_set <- c("X1", "X2", "X3")
  } else {
    stop("example_id 只能是 1, 2, 3。", call. = FALSE)
  }

  list(
    X = as.data.frame(X, check.names = FALSE),
    Y = Y,
    active_set = active_set
  )
}

compute_selection_metrics <- function(selected_list, active_set) {
  active_set <- as.character(active_set)

  per_active <- vapply(
    active_set,
    function(v) mean(vapply(selected_list, function(z) v %in% z, logical(1))),
    numeric(1)
  )

  all_active <- mean(vapply(
    selected_list,
    function(z) all(active_set %in% z),
    logical(1)
  ))

  list(
    per_active = per_active,
    all_active = all_active,
    P_s = per_active,
    P_a = all_active
  )
}

run_one_simulation_setting <- function(
    example_id,
    n,
    p,
    M_values = c(3, 10, 20),
    n_repeat = 100,
    alpha = 0.5,
    use_nystrom = FALSE,
    k = NULL,
    selection_cores = get_default_cores(),
    seed = 123,
    show_progress = TRUE,
    solver = NULL) {

  M_values <- as.integer(M_values)
  result_rows <- list()
  selected_detail <- list()

  all_selected_by_M <- lapply(M_values, function(m) vector("list", n_repeat))
  names(all_selected_by_M) <- paste0("M", M_values)
  active_set_global <- NULL

  if (show_progress) {
    message("------------------------------------------------------------")
    message(sprintf(
      "开始 simulation：example = %d, n = %d, p = %d, use_nystrom = %s, n_repeat = %d",
      example_id, n, p, use_nystrom, n_repeat
    ))
    message("------------------------------------------------------------")
  }

  for (b in seq_len(n_repeat)) {
    if (show_progress) {
      message(sprintf(
        "[simulation] 第 %d / %d 次实验开始",
        b, n_repeat
      ))
    }

    dat <- simulate_example_data(
      example_id = example_id,
      n = n,
      p = p,
      seed = seed + b - 1L
    )

    X <- dat$X
    Y <- dat$Y
    active_set_global <- dat$active_set

    DY_now <- compute_distance_matrix(Y)

    sel_res <- select_variables_bcor(
      X = X,
      DY = DY_now,
      M = max(M_values),
      alpha = alpha,
      use_nystrom = use_nystrom,
      k = k,
      n_cores = selection_cores,
      seed = seed + 10000L + b - 1L,
      solver = solver,
      show_progress = FALSE
    )

    ranked_vars <- sel_res$rank_order

    for (ii in seq_along(M_values)) {
      m <- M_values[ii]
      all_selected_by_M[[ii]][[b]] <- ranked_vars[seq_len(min(m, length(ranked_vars)))]
    }

    if (show_progress) {
      message(sprintf(
        "[simulation] 第 %d / %d 次实验完成",
        b, n_repeat
      ))
    }
  }

  for (ii in seq_along(M_values)) {
    m <- M_values[ii]
    metric <- compute_selection_metrics(all_selected_by_M[[ii]], active_set = active_set_global)

    row_now <- data.frame(
      example = example_id,
      n = n,
      p = p,
      M = m,
      use_nystrom = use_nystrom,
      stringsAsFactors = FALSE
    )

    for (v in active_set_global) {
      row_now[[v]] <- unname(metric$P_s[v])
    }
    row_now[["All"]] <- metric$P_a

    result_rows[[ii]] <- row_now
    selected_detail[[ii]] <- list(
      M = m,
      selected_list = all_selected_by_M[[ii]],
      per_active = metric$per_active,
      all_active = metric$all_active
    )
  }

  list(
    summary = do.call(rbind, result_rows),
    detail = selected_detail,
    active_set = active_set_global
  )
}

# Align real-data matrices and format result tables.
strip_first_id_column <- function(df) {
  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(df) <= 1L) return(df)
  df[, -1, drop = FALSE]
}

align_by_rownames <- function(X, Y) {

  rx <- rownames(X)
  ry <- rownames(Y)

  if (is.null(rx) || is.null(ry)) {
    return(list(X = X, Y = Y))
  }

  common_id <- intersect(rx, ry)

  if (length(common_id) > 0L) {
    return(list(
      X = X[common_id, , drop = FALSE],
      Y = Y[common_id, , drop = FALSE]
    ))
  }

  warning("X 与 Y 的行名没有交集，已按当前原始顺序直接对应。", call. = FALSE)

  if (nrow(X) != nrow(Y)) {
    stop("X 与 Y 的行名没有交集，且样本量也不一致，无法按原顺序对应。", call. = FALSE)
  }

  list(X = X, Y = Y)
}

build_selection_frequency_table <- function(selected_list) {
  vars <- sort(unique(unlist(selected_list)))
  freq <- vapply(vars, function(v) mean(vapply(selected_list, function(z) v %in% z, logical(1))), numeric(1))
  out <- data.frame(variable = vars, proportion = freq, stringsAsFactors = FALSE)
  out[order(out$proportion, decreasing = TRUE), , drop = FALSE]
}

build_repeat_selection_sheet1 <- function(selected_list) {
  freq_tab <- build_selection_frequency_table(selected_list)
  vars_union <- freq_tab$variable

  M_each <- max(vapply(selected_list, length, integer(1)))
  n_rep <- length(selected_list)

  rep_df <- as.data.frame(matrix(NA_character_, nrow = M_each, ncol = n_rep), stringsAsFactors = FALSE)
  colnames(rep_df) <- sprintf("rep_%03d", seq_len(n_rep))
  for (b in seq_len(n_rep)) {
    z <- selected_list[[b]]
    rep_df[seq_along(z), b] <- z
  }

  max_row <- max(length(vars_union), nrow(rep_df))
  out <- data.frame(
    variable = c(vars_union, rep(NA_character_, max_row - length(vars_union))),
    proportion = c(freq_tab$proportion, rep(NA_real_, max_row - nrow(freq_tab))),
    stringsAsFactors = FALSE
  )

  rep_df_pad <- as.data.frame(matrix(NA_character_, nrow = max_row, ncol = ncol(rep_df)), stringsAsFactors = FALSE)
  colnames(rep_df_pad) <- colnames(rep_df)
  rep_df_pad[seq_len(nrow(rep_df)), ] <- rep_df

  cbind(out, rep_df_pad, stringsAsFactors = FALSE)
}

map_variables_to_reference_index <- function(variables, reference_df) {
  ref <- as.data.frame(reference_df, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(ref) >= 2L) ref <- ref[, -1, drop = FALSE]
  ref_names <- colnames(ref)

  idx <- match(variables, ref_names)

  data.frame(
    variable = variables,
    reference_col_index = idx,
    stringsAsFactors = FALSE
  )
}
