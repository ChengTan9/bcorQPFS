if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", ".UTF-8"))
# Run BCor screening for the left and right hippocampal outcomes.

rm(list = ls())
options(stringsAsFactors = FALSE)

# Resolve project and data paths from the script location.
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_PATH <- if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else
  tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
SCRIPT_DIR <- if (!is.null(SCRIPT_PATH) && nzchar(SCRIPT_PATH)) {
  dirname(SCRIPT_PATH)
} else {
  getwd()
}
PROJECT_ROOT <- file.path(SCRIPT_DIR, "..", "..")
source(file.path(SCRIPT_DIR, "workflow_utils.R"), encoding = "UTF-8")
V2_ROOT <- v2_root(SCRIPT_DIR)
DATA_DIR <- v2_processed_data(V2_ROOT)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# Load the screening functions and verify package dependencies.
source(file.path(SCRIPT_DIR, "bcor_screen_select_functions.R"))
check_required_packages()

# Set screening, cache, and reproducibility parameters.
N_CORES <- 6L

FORCE_H5 <- TRUE

KEEP_H5_CACHE <- TRUE

SHOW_PROGRESS <- TRUE

SAVE_OPTIONS <- c("screened_only", "threshold", "bcor")

SEED_RIGHT <- 20260318
SEED_LEFT  <- 20260319

# Define the genotype matrix and hippocampal outcome files.
X_FILE <- v2_data_file(V2_ROOT, "ADNI_1_3_COMBINE1215.Rdata")
Y_RIGHT_FILE <- v2_data_file(V2_ROOT, "right_res_ori.Rdata")
Y_LEFT_FILE  <- v2_data_file(V2_ROOT, "left_res_ori.Rdata")

# Load the genotype matrix and move subject IDs to row names.
message("开始读取 X 文件：ADNI_1_3_COMBINE1215.Rdata")
X_raw <- load_first_object_from_rdata(X_FILE)
X_raw <- as.data.frame(X_raw, check.names = FALSE, stringsAsFactors = FALSE)

if (ncol(X_raw) < 2L) {
  stop("ADNI_1_3_COMBINE1215.Rdata 至少应包含 id 列和一个变量列。", call. = FALSE)
}

id_vec <- as.character(X_raw[[1]])
if (any(is.na(id_vec)) || any(id_vec == "")) {
  stop("ADNI_1_3_COMBINE1215.Rdata 第一列 id 存在 NA 或空字符串。", call. = FALSE)
}
if (anyDuplicated(id_vec)) {
  stop("ADNI_1_3_COMBINE1215.Rdata 第一列 id 存在重复，无法作为行名对齐。", call. = FALSE)
}

rownames(X_raw) <- id_vec

X_base <- strip_first_id_column(X_raw)
X_base <- as_numeric_df(X_base)

message(sprintf("X 读取完成：样本量 n = %d，变量数 p = %d", nrow(X_base), ncol(X_base)))

# Load and align the right hippocampal outcome.
message("开始读取右侧 Y 文件：right_res_ori.Rdata")
Y_right_raw <- load_first_object_from_rdata(Y_RIGHT_FILE)

if (is.matrix(Y_right_raw)) {
  Y_right_raw <- as.data.frame(Y_right_raw, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  Y_right_raw <- as.data.frame(Y_right_raw, check.names = FALSE, stringsAsFactors = FALSE)
}

if (is.null(rownames(Y_right_raw))) {
  stop("right_res_ori.Rdata 没有行名，无法按 id 与 X 对齐。", call. = FALSE)
}

aligned_right <- align_by_rownames(X_base, Y_right_raw)
X_right <- as_numeric_df(aligned_right$X)
Y_right <- as_numeric_df(aligned_right$Y)

message(sprintf(
  "右侧数据对齐完成：X 维度 = %d x %d，Y 维度 = %d x %d",
  nrow(X_right), ncol(X_right), nrow(Y_right), ncol(Y_right)
))

# Load and align the left hippocampal outcome.
message("开始读取左侧 Y 文件：left_res_ori.Rdata")
Y_left_raw <- load_first_object_from_rdata(Y_LEFT_FILE)

if (is.matrix(Y_left_raw)) {
  Y_left_raw <- as.data.frame(Y_left_raw, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  Y_left_raw <- as.data.frame(Y_left_raw, check.names = FALSE, stringsAsFactors = FALSE)
}

if (is.null(rownames(Y_left_raw))) {
  stop("left_res_ori.Rdata 没有行名，无法按 id 与 X 对齐。", call. = FALSE)
}

aligned_left <- align_by_rownames(X_base, Y_left_raw)
X_left <- as_numeric_df(aligned_left$X)
Y_left <- as_numeric_df(aligned_left$Y)

message(sprintf(
  "左侧数据对齐完成：X 维度 = %d x %d，Y 维度 = %d x %d",
  nrow(X_left), ncol(X_left), nrow(Y_left), ncol(Y_left)
))

# Screen predictors for the right hippocampal outcome.
message("============================================================")
message("开始进行右侧 screening")
message(sprintf("当前并行核数 N_CORES = %d", N_CORES))
message(sprintf("FORCE_H5 = %s，KEEP_H5_CACHE = %s", FORCE_H5, KEEP_H5_CACHE))
message("注意：screening 一开始会先计算 Y 的距离矩阵 DY，这一步本身没有进度条。")
message("如果 FORCE_H5 = TRUE，则后续会依次出现三段进度条：生成 h5、计算阈值、扫描全部变量 BCor。")
message("============================================================")

right_screen_res <- screen_bcor(
  X = X_right,
  Y = Y_right,
  n_cores = N_CORES,
  save_dir = DATA_DIR,
  output_prefix = "ADNI_right_screening",
  save_options = SAVE_OPTIONS,
  memory_limit_gb = 8,
  reserve_gb = 1.5,
  force_h5 = FORCE_H5,
  keep_h5_cache = KEEP_H5_CACHE,
  show_progress = SHOW_PROGRESS,
  seed = SEED_RIGHT,
  dataset_name = "X"
)

message("右侧 screening 完成。")
message(sprintf("右侧保留变量数：%d", ncol(right_screen_res$X_screened)))
message(sprintf("右侧阈值 tau_n：%.6f", right_screen_res$tau_n))
message("右侧输出文件：")
print(right_screen_res$saved_files)
if (!is.null(right_screen_res$h5_path)) {
  message(sprintf("右侧 h5 文件：%s", right_screen_res$h5_path))
}

# Screen predictors for the left hippocampal outcome.
message("============================================================")
message("开始进行左侧 screening")
message(sprintf("当前并行核数 N_CORES = %d", N_CORES))
message(sprintf("FORCE_H5 = %s，KEEP_H5_CACHE = %s", FORCE_H5, KEEP_H5_CACHE))
message("注意：screening 一开始会先计算 Y 的距离矩阵 DY，这一步本身没有进度条。")
message("如果 FORCE_H5 = TRUE，则后续会依次出现三段进度条：生成 h5、计算阈值、扫描全部变量 BCor。")
message("============================================================")

left_screen_res <- screen_bcor(
  X = X_left,
  Y = Y_left,
  n_cores = N_CORES,
  save_dir = DATA_DIR,
  output_prefix = "ADNI_left_screening",
  save_options = SAVE_OPTIONS,
  memory_limit_gb = 8,
  reserve_gb = 1.5,
  force_h5 = FORCE_H5,
  keep_h5_cache = KEEP_H5_CACHE,
  show_progress = SHOW_PROGRESS,
  seed = SEED_LEFT,
  dataset_name = "X"
)

message("左侧 screening 完成。")
message(sprintf("左侧保留变量数：%d", ncol(left_screen_res$X_screened)))
message(sprintf("左侧阈值 tau_n：%.6f", left_screen_res$tau_n))
message("左侧输出文件：")
print(left_screen_res$saved_files)
if (!is.null(left_screen_res$h5_path)) {
  message(sprintf("左侧 h5 文件：%s", left_screen_res$h5_path))
}

# Save complete screening objects for later inspection.
save_rdata_named(
  right_screen_res,
  file.path(DATA_DIR, "ADNI_right_screening_full_result.RData"),
  object_name = "right_screen_res"
)

save_rdata_named(
  left_screen_res,
  file.path(DATA_DIR, "ADNI_left_screening_full_result.RData"),
  object_name = "left_screen_res"
)

message("============================================================")
message("全部 screening 已完成。")
message("已保存的核心结果包括：")
message("1. ADNI_right_screening_selected.RData")
message("2. ADNI_right_screening_threshold.csv")
message("3. ADNI_right_screening_bcor.csv")
message("4. ADNI_left_screening_selected.RData")
message("5. ADNI_left_screening_threshold.csv")
message("6. ADNI_left_screening_bcor.csv")
message("另外还保存了完整结果对象：")
message("7. ADNI_right_screening_full_result.RData")
message("8. ADNI_left_screening_full_result.RData")
message("============================================================")

# Metadata consumed by 02/05: preserve global SNP indices and screening cutoff.
for (side in c("left", "right")) {
  result <- if (side == "left") left_screen_res else right_screen_res
  variables <- colnames(result$X_screened)
  original_index <- match(variables, colnames(X_base))
  stopifnot(!anyNA(original_index))
  saveRDS(list(variables = variables, original_index = original_index,
               threshold = result$tau_n, n = nrow(result$X_screened)),
          file.path(DATA_DIR, paste0("ADNI_", side, "_screening_metadata.rds")))
}
