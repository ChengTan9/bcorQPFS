if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", ".UTF-8"))
# Run HSIC Lasso on the dHSIC-screened ADNI predictors.

rm(list = ls())

# Load output and progress packages.
library(openxlsx)
library(pbapply)

# Resolve project, data, and output paths.
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else
  tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_dir <- if (!is.null(script_path) && nzchar(script_path)) {
  dirname(script_path)
} else {
  getwd()
}
source(file.path(script_dir, "hsic_lasso_functions.R"))
project_root <- file.path(script_dir, "..", "..")
source(file.path(script_dir, "workflow_utils.R"), encoding = "UTF-8")
V2_ROOT <- v2_root(script_dir)
data_dir <- v2_processed_data(V2_ROOT)
out_dir <- file.path(V2_ROOT, "results", "hsic_lasso_selection", format(Sys.time(), "%Y%m%d_%H%M%S"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("script_dir =", script_dir, "\n")
cat("data_dir   =", data_dir, "\n")
cat("out_dir    =", out_dir, "\n")

# Locate screening results, outcomes, and the reference genotype matrix.
path_right_x <- v2_data_file(V2_ROOT, "ADNI_dHSIC_selected_vs_right_res_ori.RData")
path_right_y <- v2_data_file(V2_ROOT, "right_res_ori.Rdata")
path_left_x  <- v2_data_file(V2_ROOT, "ADNI_dHSIC_selected_vs_left_res_ori.RData")
path_left_y  <- v2_data_file(V2_ROOT, "left_res_ori.Rdata")
path_adni    <- v2_data_file(V2_ROOT, "ADNI_1_3_COMBINE1215.Rdata")

cat("\n找到文件：\n")
cat("right_x:", path_right_x, "\n")
cat("right_y:", path_right_y, "\n")
cat("left_x :", path_left_x,  "\n")
cat("left_y :", path_left_y,  "\n")
cat("adni   :", path_adni,    "\n")

# Extract the selected predictor matrix from a screening result.
load_hsic_screened_x <- function(path) {
  e <- new.env(parent = emptyenv())
  objs <- load(path, envir = e)
  if ("screening_result" %in% objs &&
      is.list(e$screening_result) &&
      !is.null(e$screening_result$X_selected)) {
    return(e$screening_result$X_selected)
  }
  stop("HSIC screening 文件缺少 screening_result$X_selected：", path)
}

# Load screened predictors, outcomes, and reference variable names.
right_x <- as_numeric_matrix(load_hsic_screened_x(path_right_x), "right_dHSIC_screened_X")
right_y <- as_numeric_matrix(load_rdata_auto(path_right_y), "right_res_ori")
left_x  <- as_numeric_matrix(load_hsic_screened_x(path_left_x),  "left_dHSIC_screened_X")
left_y  <- as_numeric_matrix(load_rdata_auto(path_left_y),  "left_res_ori")

adni_raw <- load_rdata_auto(path_adni)
adni_df  <- as.data.frame(adni_raw, check.names = FALSE)

if (ncol(adni_df) < 2) {
  stop("ADNI_1_3_COMBINE1215.Rdata 至少需要两列。")
}

adni_noid <- adni_df[, -1, drop = FALSE]
adni_noid_colnames <- trimws(colnames(adni_noid))

check_unique_colnames(right_x, "right_x")
check_unique_colnames(left_x,  "left_x")

if (anyDuplicated(adni_noid_colnames)) {
  dupn <- unique(adni_noid_colnames[duplicated(adni_noid_colnames)])
  stop("ADNI_1_3_COMBINE1215 删除第一列后仍存在重复列名，无法唯一匹配。重复示例：",
       paste(head(dupn, 10), collapse = ", "))
}

colnames(right_x) <- trimws(colnames(right_x))
colnames(left_x)  <- trimws(colnames(left_x))

cat("\n数据维度检查：\n")
cat("right_x :", paste(dim(right_x), collapse = " x "), "\n")
cat("right_y :", paste(dim(right_y), collapse = " x "), "\n")
cat("left_x  :", paste(dim(left_x),  collapse = " x "), "\n")
cat("left_y  :", paste(dim(left_y),  collapse = " x "), "\n")
cat("adni_noid :", paste(dim(adni_noid), collapse = " x "), "\n")

if (nrow(right_x) != nrow(right_y)) stop("右侧 X/Y 行数不一致。")
if (nrow(left_x)  != nrow(left_y))  stop("左侧 X/Y 行数不一致。")

# Set HSIC Lasso selection and parallel parameters.
TOP_K <- 30L

B_use <- 0L

M_use <- 1L

n_cores_right <- choose_n_cores(n_tasks = ncol(right_x), reserve_cores = 5L, use_physical_cores = TRUE)
n_cores_left  <- choose_n_cores(n_tasks = ncol(left_x),  reserve_cores = 5L, use_physical_cores = TRUE)

cat("\n并行设置：\n")
cat("右侧使用核心数 =", n_cores_right, "\n")
cat("左侧使用核心数 =", n_cores_left,  "\n")
cat("说明：这里只在“单个特征 kernel 计算”层面并行；左右两侧仍顺序执行，避免内存过高。\n")

# Fit HSIC Lasso for the right outcome.
cat("\n================ 开始右侧 HSIC-Lasso ================\n")

res_right <- run_hsiclasso_regression_py_clone_parallel(
  X_nxp = right_x,
  Y_nxq = right_y,
  featname = colnames(right_x),
  num_feat = TOP_K,
  B = B_use,
  M = M_use,
  max_neighbors = 10L,
  n_cores = n_cores_right,
  show_progress = TRUE
)

export_right <- make_export_df(
  selected_features = res_right$selected_features,
  adni_noid_colnames = adni_noid_colnames
)

if (anyNA(export_right$column_in_ADNI_1_3_COMBINE1215_after_drop_id)) {
  bad_vars <- export_right$important_variable[
    is.na(export_right$column_in_ADNI_1_3_COMBINE1215_after_drop_id)
  ]
  warning("右侧以下变量未能在 ADNI_1_3_COMBINE1215 删除第一列后匹配到：",
          paste(bad_vars, collapse = ", "))
}

# Fit HSIC Lasso for the left outcome.
cat("\n================ 开始左侧 HSIC-Lasso ================\n")

res_left <- run_hsiclasso_regression_py_clone_parallel(
  X_nxp = left_x,
  Y_nxq = left_y,
  featname = colnames(left_x),
  num_feat = TOP_K,
  B = B_use,
  M = M_use,
  max_neighbors = 10L,
  n_cores = n_cores_left,
  show_progress = TRUE
)

export_left <- make_export_df(
  selected_features = res_left$selected_features,
  adni_noid_colnames = adni_noid_colnames
)

if (anyNA(export_left$column_in_ADNI_1_3_COMBINE1215_after_drop_id)) {
  bad_vars <- export_left$important_variable[
    is.na(export_left$column_in_ADNI_1_3_COMBINE1215_after_drop_id)
  ]
  warning("左侧以下变量未能在 ADNI_1_3_COMBINE1215 删除第一列后匹配到：",
          paste(bad_vars, collapse = ", "))
}

# Export the selected variables and reference column indices.
right_xlsx <- file.path(out_dir, "HSICLasso_right_top30.xlsx")
left_xlsx  <- file.path(out_dir, "HSICLasso_left_top30.xlsx")

write.xlsx(export_right, right_xlsx, overwrite = TRUE)
write.xlsx(export_left,  left_xlsx,  overwrite = TRUE)

cat("\n================ 全部完成 ================\n")
cat("右侧 Excel：", right_xlsx, "\n")
cat("左侧 Excel：", left_xlsx,  "\n")

cat("\n右侧前 10 个结果：\n")
print(utils::head(export_right, 10), row.names = FALSE)

cat("\n左侧前 10 个结果：\n")
print(utils::head(export_left, 10), row.names = FALSE)

