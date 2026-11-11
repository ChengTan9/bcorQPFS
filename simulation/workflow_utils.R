# Code V2 paths. Controlled-access ADNI inputs are kept outside this release.
v2_root <- function(script_dir) {
  candidate <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
  stopifnot(dir.exists(file.path(candidate, "realdata")), dir.exists(file.path(candidate, "simulation")))
  candidate
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
