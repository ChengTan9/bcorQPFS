# Code V2 paths. Controlled-access ADNI inputs are kept outside this release.
v2_root <- function(script_dir) {
  candidate <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
  stopifnot(dir.exists(file.path(candidate, "realdata")), dir.exists(file.path(candidate, "simulation")))
  candidate
}
v2_external_data <- function(root) {
  Sys.getenv("BCORQPFS_DATA_DIR", unset = v2_processed_data(root))
}
v2_processed_data <- function(root) file.path(root, "data", "processed")
v2_data_file <- function(root, name) {
  for (folder in c(v2_processed_data(root), v2_external_data(root))) {
    if (!dir.exists(folder)) next
    files <- list.files(folder, full.names = TRUE)
    hits <- files[tolower(basename(files)) == tolower(name)]
    if (length(hits)) return(normalizePath(hits[[1L]], winslash = "/", mustWork = TRUE))
  }
  stop("Missing input ", name, ". Set BCORQPFS_DATA_DIR to your authorized processed ADNI directory.", call. = FALSE)
}
v2_bcor_archive <- function(root) {
  explicit <- Sys.getenv("BCORQPFS_BCOR_ARCHIVE_DIR", unset = "")
  if (nzchar(explicit)) return(explicit)
  folder <- file.path(root, "results", "adni_single_run")
  runs <- if (dir.exists(folder)) list.dirs(folder, recursive = FALSE, full.names = TRUE) else character(0)
  required <- c("BcorQPFS_left_full_result.rds", "BcorQPFS_right_full_result.rds")
  runs <- runs[vapply(runs, function(p) all(file.exists(file.path(p, required))), logical(1))]
  runs <- runs[vapply(runs, function(p) {
    status <- file.path(p, "run_status.txt")
    file.exists(status) && any(readLines(status, warn=FALSE)=="status=completed")
  }, logical(1))]
  if (length(runs)) return(runs[order(basename(runs), decreasing = TRUE)][[1L]])
  stop("No complete code_V2 single-fit run found. Run realdata/02_bcorqpfs_selection.R first, or explicitly set BCORQPFS_BCOR_ARCHIVE_DIR for historical reproduction.", call. = FALSE)
}
v2_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required R packages before running: ", paste(missing, collapse = ", "), call. = FALSE)
}

v2_pick_run <- function(root, family, required) {
  folder <- file.path(root, "results", family)
  dirs <- list.dirs(folder, full.names = TRUE, recursive = FALSE)
  dirs <- dirs[!grepl("^smoke_", basename(dirs))]
  dirs <- dirs[vapply(dirs, function(p) all(file.exists(file.path(p, required))), logical(1))]
  dirs <- dirs[vapply(dirs, function(p) {
    status <- file.path(p, "run_status.txt")
    if (file.exists(status)) return(any(readLines(status, warn=FALSE)=="status=completed"))
    startsWith(basename(p), "archived_")
  }, logical(1))]
  if (!length(dirs)) stop("No completed ", family, " run found. Required: ", paste(required, collapse = ", "))
  dirs[order(file.info(file.path(dirs, required[[1L]]))$mtime, decreasing = TRUE)][[1L]]
}
v2_summary_input <- function(root, family, required) {
  args <- commandArgs(TRUE)
  if (!length(args)) return(v2_pick_run(root, family, required))
  if (length(args)!=1L || !startsWith(args, "--input=")) stop("Supported argument: --input=RUN_DIRECTORY")
  normalizePath(sub("^--input=", "", args), winslash = "/", mustWork = TRUE)
}
