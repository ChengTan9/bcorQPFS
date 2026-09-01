# Screen ADNI predictors with dHSIC for both hippocampal outcomes.

rm(list = ls())

options(stringsAsFactors = FALSE)

# Load the required packages.
install_if_missing <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    install.packages(miss, dependencies = TRUE)
  }
}

install_if_missing(c("dHSIC", "pbapply"))

library(dHSIC)
library(pbapply)
library(parallel)

# Locate files and validate loaded objects.
`%||%` <- function(a, b) if (!is.null(a)) a else b

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
    stop(name, " 没有列名，无法输出变量名。")
  }
  if (anyDuplicated(cn)) {
    dupn <- unique(cn[duplicated(cn)])
    stop(name, " 存在重复列名，无法唯一匹配。重复示例：",
         paste(head(dupn, 10), collapse = ", "))
  }
}

choose_n_cores <- function(reserve_cores = 3L) {
  dc <- tryCatch(parallel::detectCores(logical = TRUE),
                 error = function(e) NA_integer_)
  if (is.na(dc) || dc <= 1L) return(1L)
  as.integer(max(1L, dc - reserve_cores))
}

# Construct Gaussian Gram matrices for outcomes and predictors.
gaussian_gram_from_distance <- function(D, sigma = NULL) {
  D <- as.matrix(D)

  if (!is.numeric(D)) {
    stop("D 必须是数值型距离矩阵。")
  }
  if (nrow(D) != ncol(D)) {
    stop("D 必须是方阵。")
  }

  if (is.null(sigma)) {
    upper_vals <- D[upper.tri(D)]
    upper_vals <- upper_vals[is.finite(upper_vals) & upper_vals > 0]

    if (length(upper_vals) == 0) {
      sigma <- 1
    } else {
      sigma <- stats::median(upper_vals)
      if (!is.finite(sigma) || sigma <= 0) sigma <- 1
    }
  }

  K <- exp(-(D^2) / (2 * sigma^2))
  diag(K) <- 1

  list(K = K, sigma = sigma)
}

gaussian_gram_for_1d_feature <- function(x) {
  x <- as.numeric(x)

  if (anyNA(x)) {
    stop("单个特征中存在 NA。")
  }

  if (all(x == x[1])) {
    n <- length(x)
    return(list(K = matrix(1, n, n), sigma = 0))
  }

  Dx <- abs(outer(x, x, "-"))
  gaussian_gram_from_distance(Dx, sigma = NULL)
}

# Cache the outcome Gram matrix for repeated screening runs.
prepare_y_kernel_cache <- function(Y,
                                   cache_rds_path = NULL,
                                   screen_name = "Y") {
  Y <- as_numeric_matrix(Y, screen_name)

  if (!is.null(cache_rds_path) && file.exists(cache_rds_path)) {
    cat("[Y-kernel] 检测到已有缓存，直接读取：", cache_rds_path, "\n")
    obj <- readRDS(cache_rds_path)

    if (!is.null(obj$n) && obj$n == nrow(Y)) {
      return(obj)
    } else {
      cat("[Y-kernel] 缓存维度不匹配，重新计算。\n")
    }
  }

  cat("[Y-kernel] 开始计算 ", screen_name, " 的 DY ...\n", sep = "")
  DY <- as.matrix(stats::dist(Y))

  cat("[Y-kernel] 开始由 DY 构造 Ky ...\n")
  Ky_obj <- gaussian_gram_from_distance(DY, sigma = NULL)

  out <- list(
    DY = DY,
    Ky = Ky_obj$K,
    sigma_y = Ky_obj$sigma,
    n = nrow(Y)
  )

  if (!is.null(cache_rds_path)) {
    saveRDS(out, cache_rds_path)
    cat("[Y-kernel] 缓存已写入：", cache_rds_path, "\n")
  }

  out
}

# Compute dHSIC for one predictor and the fixed outcome kernel.
compute_dhsic_one_feature <- function(x, Ky) {
  x <- as.numeric(x)

  if (all(x == x[1])) {
    return(0)
  }

  Kx_obj <- gaussian_gram_for_1d_feature(x)
  Kx <- Kx_obj$K

  val <- tryCatch(
    dHSIC::dhsic(K = list(Kx, Ky))$dHSIC,
    error = function(e) NA_real_
  )

  as.numeric(val)
}

# Screen predictors in column chunks and save scores and selected variables.
screen_dhsic_direct_chunk <- function(adni_rdata_path,
                                      y_rdata_path,
                                      score_csv_path,
                                      selected_rdata_path,
                                      y_cache_rds_path,
                                      chunk_cols = 500L,
                                      n_cores = NULL,
                                      screen_name = "screening") {
  cat("\n====================================================\n")
  cat("[", screen_name, "] 开始 screening\n", sep = "")
  cat("====================================================\n")

  Y_raw <- load_rdata_auto(y_rdata_path)
  Y <- as_numeric_matrix(Y_raw, basename(y_rdata_path))

  n <- nrow(Y)
  top_k <- floor(n / log(n))
  top_k <- max(1L, top_k)

  cat("[", screen_name, "] Y 维度：", nrow(Y), " x ", ncol(Y), "\n", sep = "")
  cat("[", screen_name, "] 计划选出变量数 floor(n/log(n)) = ", top_k, "\n", sep = "")

  y_cache <- prepare_y_kernel_cache(
    Y = Y,
    cache_rds_path = y_cache_rds_path,
    screen_name = screen_name
  )
  Ky <- y_cache$Ky

  cat("[", screen_name, "] 开始读取 ADNI 原始数据 ...\n", sep = "")

  adni_raw <- load_rdata_auto(adni_rdata_path, preferred_name = "ADNI_1_3_COMBINE1215")

  adni_df <- as.data.frame(adni_raw, check.names = FALSE)

  if (ncol(adni_df) < 2) {
    stop("[", screen_name, "] ADNI 数据至少需要两列。")
  }

  cat("[", screen_name, "] ADNI 原始维度：", nrow(adni_df), " x ", ncol(adni_df), "\n", sep = "")
  cat("[", screen_name, "] 原始第一列列名：", colnames(adni_df)[1], "\n", sep = "")
  cat("[", screen_name, "] 现在先删除第一列，再做数值化转换。\n", sep = "")

  adni_noid_df <- adni_df[, -1, drop = FALSE]

  X <- as_numeric_matrix(adni_noid_df, "ADNI_1_3_COMBINE1215_after_drop_first")

  rm(adni_raw, adni_df, adni_noid_df)
  invisible(gc())

  check_unique_colnames(X, "ADNI_1_3_COMBINE1215_after_drop_first")

  if (nrow(X) != n) {
    stop("[", screen_name, "] X 和 Y 的样本数不一致：X 有 ", nrow(X), " 行，Y 有 ", n, " 行。")
  }

  p <- ncol(X)
  x_colnames <- colnames(X)

  cat("[", screen_name, "] X 维度（删除第一列后）：", nrow(X), " x ", p, "\n", sep = "")

  if (is.null(n_cores)) {
    n_cores <- choose_n_cores(reserve_cores = 3L)
  }
  n_cores <- max(1L, min(as.integer(n_cores), p))

  cat("[", screen_name, "] 并行核心数：", n_cores, "\n", sep = "")
  cat("[", screen_name, "] 每块列数 chunk_cols：", chunk_cols, "\n", sep = "")

  total_chunks <- ceiling(p / chunk_cols)

  dhsic_scores <- rep(NA_real_, p)

  pb_chunk <- txtProgressBar(min = 0, max = total_chunks, style = 3)
  on.exit(close(pb_chunk), add = TRUE)

  cl <- NULL
  if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterEvalQ(cl, {
      library(dHSIC)
      NULL
    })

    parallel::clusterExport(
      cl,
      varlist = c(
        "Ky",
        "compute_dhsic_one_feature",
        "gaussian_gram_for_1d_feature",
        "gaussian_gram_from_distance"
      ),
      envir = environment()
    )
  }

  for (chunk_id in seq_len(total_chunks)) {
    start_col <- (chunk_id - 1L) * chunk_cols + 1L
    end_col   <- min(chunk_id * chunk_cols, p)

    cat("\n[", screen_name, "] 正在处理第 ", chunk_id, "/", total_chunks,
        " 个块：列 ", start_col, " - ", end_col, "\n", sep = "")

    X_block <- X[, start_col:end_col, drop = FALSE]
    p_block <- ncol(X_block)

    if (!is.null(cl)) {

      parallel::clusterExport(cl, varlist = c("X_block"), envir = environment())

      block_scores <- unlist(
        pbapply::pblapply(
          X = seq_len(p_block),
          FUN = function(j) compute_dhsic_one_feature(X_block[, j], Ky),
          cl = cl
        )
      )
    } else {

      pb_inner <- txtProgressBar(min = 0, max = p_block, style = 3)
      on.exit(close(pb_inner), add = TRUE)

      block_scores <- numeric(p_block)
      for (j in seq_len(p_block)) {
        block_scores[j] <- compute_dhsic_one_feature(X_block[, j], Ky)
        setTxtProgressBar(pb_inner, j)
      }
      close(pb_inner)
    }

    dhsic_scores[start_col:end_col] <- block_scores

    setTxtProgressBar(pb_chunk, chunk_id)

    rm(X_block, block_scores)
    invisible(gc())
  }

  score_df <- data.frame(
    column_in_ADNI_after_drop_first = seq_len(p),
    variable_name = x_colnames,
    dhsic = dhsic_scores,
    stringsAsFactors = FALSE
  )

  score_df <- score_df[order(score_df$dhsic, decreasing = TRUE, na.last = TRUE), , drop = FALSE]

  write.csv(score_df, file = score_csv_path, row.names = FALSE)
  cat("\n[", screen_name, "] dHSIC 全部结果已保存到 CSV：", score_csv_path, "\n", sep = "")

  top_k_use <- min(top_k, nrow(score_df))

  selected_idx   <- score_df$column_in_ADNI_after_drop_first[seq_len(top_k_use)]
  selected_names <- score_df$variable_name[seq_len(top_k_use)]

  X_selected <- X[, selected_idx, drop = FALSE]
  colnames(X_selected) <- selected_names

  screening_result <- list(
    X_selected = X_selected,
    selected_index = selected_idx,
    selected_names = selected_names,
    top_k = top_k_use,
    n = n,
    rule = "floor(n / log(n))",
    score_table_top = score_df[seq_len(top_k_use), , drop = FALSE],
    y_kernel_cache = list(
      cache_rds_path = y_cache_rds_path,
      sigma_y = y_cache$sigma_y
    )
  )

  save(screening_result, file = selected_rdata_path)
  cat("[", screen_name, "] 选中变量已保存到 RData：", selected_rdata_path, "\n", sep = "")

  invisible(list(
    score_df = score_df,
    screening_result = screening_result
  ))
}

# Define input, cache, and output paths.
script_dir <- get_script_dir()

project_root <- file.path(script_dir, "..", "..")
data_dir <- file.path(project_root, "code", "data")

if (!dir.exists(data_dir)) {
  stop("未找到 data 文件夹：", data_dir)
}

cat("script_dir = ", script_dir, "\n", sep = "")
cat("data_dir   = ", data_dir, "\n", sep = "")

path_adni  <- find_file_case_insensitive(data_dir, "ADNI_1_3_COMBINE1215.Rdata")

path_right <- find_file_case_insensitive(data_dir, "right_res_ori.Rdata")
path_left  <- find_file_case_insensitive(data_dir, "left_res_ori.Rdata")

cat("\n找到输入文件：\n")
cat("ADNI  : ", path_adni,  "\n", sep = "")
cat("RIGHT : ", path_right, "\n", sep = "")
cat("LEFT  : ", path_left,  "\n", sep = "")

path_right_y_cache <- file.path(data_dir, "right_res_ori_DY_Ky_cache.rds")
path_right_csv     <- file.path(data_dir, "ADNI_dHSIC_screening_vs_right_res_ori.csv")
path_right_rdata   <- file.path(data_dir, "ADNI_dHSIC_selected_vs_right_res_ori.RData")

path_left_y_cache <- file.path(data_dir, "left_res_ori_DY_Ky_cache.rds")
path_left_csv     <- file.path(data_dir, "ADNI_dHSIC_screening_vs_left_res_ori.csv")
path_left_rdata   <- file.path(data_dir, "ADNI_dHSIC_selected_vs_left_res_ori.RData")

cat("\n输出文件将写入：\n")
cat("RIGHT Y缓存 : ", path_right_y_cache, "\n", sep = "")
cat("RIGHT CSV   : ", path_right_csv, "\n", sep = "")
cat("RIGHT RData : ", path_right_rdata, "\n", sep = "")
cat("LEFT Y缓存  : ", path_left_y_cache, "\n", sep = "")
cat("LEFT CSV    : ", path_left_csv, "\n", sep = "")
cat("LEFT RData  : ", path_left_rdata, "\n", sep = "")

chunk_cols <- 500L

n_cores <- choose_n_cores(reserve_cores = 3L)

cat("\n运行参数：\n")
cat("chunk_cols = ", chunk_cols, "\n", sep = "")
cat("n_cores    = ", n_cores, "\n", sep = "")

# Run screening for the right and left outcomes.
res_right <- screen_dhsic_direct_chunk(
  adni_rdata_path = path_adni,
  y_rdata_path = path_right,
  score_csv_path = path_right_csv,
  selected_rdata_path = path_right_rdata,
  y_cache_rds_path = path_right_y_cache,
  chunk_cols = chunk_cols,
  n_cores = n_cores,
  screen_name = "RIGHT"
)

res_left <- screen_dhsic_direct_chunk(
  adni_rdata_path = path_adni,
  y_rdata_path = path_left,
  score_csv_path = path_left_csv,
  selected_rdata_path = path_left_rdata,
  y_cache_rds_path = path_left_y_cache,
  chunk_cols = chunk_cols,
  n_cores = n_cores,
  screen_name = "LEFT"
)

cat("\n====================================================\n")
cat("全部完成。\n")
cat("====================================================\n")

cat("\nRIGHT 前 10 个变量：\n")
print(utils::head(res_right$score_df, 10), row.names = FALSE)

cat("\nLEFT 前 10 个变量：\n")
print(utils::head(res_left$score_df, 10), row.names = FALSE)

cat("\n生成的文件：\n")
cat(path_right_y_cache, "\n")
cat(path_right_csv, "\n")
cat(path_right_rdata, "\n")
cat(path_left_y_cache, "\n")
cat(path_left_csv, "\n")
cat(path_left_rdata, "\n")
