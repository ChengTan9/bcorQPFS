#!/usr/bin/env Rscript
# Recreate the complete-QP and frequency tables from recorded results, without loading participant data.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", ".UTF-8"))
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]), winslash = "/"))
source(file.path(script_dir, "workflow_utils.R"), encoding = "UTF-8")
root <- v2_root(script_dir)
input <- v2_summary_input(root, "adni_stability", c("left_side_summary.csv", "right_side_summary.csv"))
out <- file.path(root, "results", "summaries", "adni", basename(input))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
annotations <- read.csv(file.path(script_dir, "variant_annotations_grch37.csv"), fileEncoding = "UTF-8-BOM", check.names = FALSE)
annotation_key <- paste(annotations$chromosome, annotations$position, sep = ":")
annotate <- function(x) {
  parts <- strsplit(x$variable, ".", fixed = TRUE)
  key <- vapply(parts, function(p) paste(sub("^X", "", p[[1L]]), p[[2L]], sep = ":"), character(1))
  matched <- match(key, annotation_key)
  labels <- annotations[matched, , drop = FALSE]
  rownames(labels) <- NULL
  cbind(x, labels, main_text_selected = x$archived_rank <= 20L,
        annotation_status = ifelse(is.na(matched), "not_in_reviewed_annotation_set", "reviewed_GRCh37"))
}
full <- frequent <- summaries <- list()
for (side in c("left", "right")) {
  exact <- read.csv(file.path(input, paste0(side, "_exact_qp_top20.csv")), check.names = FALSE)
  freq <- read.csv(file.path(input, paste0(side, "_landmark_frequency_ranking.csv")), check.names = FALSE)
  diag <- read.csv(file.path(input, paste0(side, "_draw_diagnostics.csv")), check.names = FALSE)
  info <- read.csv(file.path(input, paste0(side, "_side_summary.csv")), check.names = FALSE)
  stopifnot(nrow(exact)==20L, sum(freq$archived_rank <= 20L)==20L)
  full[[side]] <- annotate(exact)
  # Retain the published frequency Top-20 and every reported single-draw SNP.
  keep <- freq$frequency_rank <= 20L | freq$archived_rank <= 20L
  frequent[[side]] <- annotate(freq[keep, , drop = FALSE])
  frequent[[side]]$selection_percentage <- 100 * frequent[[side]]$top20_fraction
  summaries[[side]] <- data.frame(side, n = info$n, q = info$q, landmarks = info$s, draws = info$draws,
    reported_exact_overlap = info$exact_overlap_archived,
    frequency_top20_exact_overlap = sum(freq$frequency_rank<=20L & freq$exact_rank<=20L),
    always_selected = sum(freq$top20_count==info$draws),
    unique_ever_selected = sum(freq$top20_count>0L),
    successful_additional_fits = sum(diag$run>0 & diag$status %in% c("Solved", "AlmostSolved")),
    stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, full), file.path(out, "complete_qp_top20_both_sides.csv"), row.names = FALSE)
write.csv(do.call(rbind, frequent), file.path(out, "frequency_top20_plus_reported_both_sides.csv"), row.names = FALSE)
write.csv(do.call(rbind, summaries), file.path(out, "stability_summary.csv"), row.names = FALSE)
writeLines(c(paste0("source_run=", normalizePath(input, winslash = "/")),
 "The main-text list is the archived single landmark draw; repeated draws are an auxiliary sensitivity analysis.",
 "Gene/allele labels are copied from the reviewed GRCh37 manuscript tables, not inferred from genotype column coding."),
 file.path(out, "summary_provenance.txt"))
print(do.call(rbind, summaries), row.names = FALSE)
cat("Saved:", out, "\n")
