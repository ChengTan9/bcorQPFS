#!/usr/bin/env Rscript
# Summarize the one-fit Example 1 timing grid; no rerun of the simulations.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", ".UTF-8"))
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]), winslash = "/"))
source(file.path(script_dir, "workflow_utils.R"), encoding = "UTF-8")
root <- v2_root(script_dir)
input <- v2_summary_input(root, "example1_timing", "example1_timing_full_grid.csv")
csv_path <- file.path(input, "example1_timing_full_grid.csv")
if (!file.exists(csv_path)) csv_path <- file.path(input, "example1_timing_smoke.csv")
if (!file.exists(csv_path)) stop("Input directory contains neither example1_timing_full_grid.csv nor example1_timing_smoke.csv")
x <- read.csv(csv_path, check.names = FALSE)
baseline <- x[x$r==1, c("n", "p", "total_seconds")]
names(baseline)[[3L]] <- "complete_qp_seconds"
stopifnot(!anyDuplicated(paste(baseline$n, baseline$p)))
tab <- merge(x, baseline, by = c("n", "p"), sort = FALSE)
tab$relative_percentage <- 100*tab$total_seconds/tab$complete_qp_seconds
tab$reduction_percentage <- 100-tab$relative_percentage
mean_tab <- aggregate(cbind(relative_percentage, reduction_percentage) ~ r, tab, mean)
out <- file.path(root, "results", "summaries", "example1_timing", basename(input))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(tab, file.path(out, "timing_with_relative_percentages.csv"), row.names = FALSE)
write.csv(mean_tab, file.path(out, "mean_relative_percentages.csv"), row.names = FALSE)
wide <- reshape(tab[,c("n","p","r","total_seconds","relative_percentage")], idvar=c("n","p"), timevar="r", direction="wide")
write.csv(wide, file.path(out, "timing_table_wide.csv"), row.names = FALSE)
writeLines(c(paste0("source_run=", normalizePath(input, winslash = "/")),paste0("source_csv=",basename(csv_path)),
 "Timing includes relevance, required dependence entries, factorization, and QP formulation/solution.",
 "Data generation, warm-up, preprocessing and screening are excluded. Each n,p grid point shares one generated data set across methods.",
 "Mean percentages give each n,p setting equal weight; these are recorded single-fit measurements."), file.path(out, "timing_contract.txt"))
print(mean_tab, row.names = FALSE)
cat("Saved:", out, "\n")
