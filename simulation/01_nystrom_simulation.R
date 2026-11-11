if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", ".UTF-8"))
# Run the Nystrom simulation reported in Tables 1--3.

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
candidate_dirs <- c(
  if (length(file_arg)) dirname(sub("^--file=", "", file_arg[[1L]])) else character(0),
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    dirname(rstudioapi::getSourceEditorContext()$path)
  } else {
    character(0)
  },
  getwd(),
  file.path(getwd(), "code_V2", "simulation")
)
candidate_dirs <- unique(candidate_dirs[nzchar(candidate_dirs)])
script_dir <- candidate_dirs[file.exists(file.path(candidate_dirs, "simulation_functions.R"))][1L]
if (is.na(script_dir)) stop("Cannot locate simulation_functions.R.", call. = FALSE)
functions_file <- file.path(script_dir, "simulation_functions.R")
source(functions_file, encoding = "UTF-8")

# Set the manuscript simulation parameters.
config <- list(
  examples = 1:3,
  kappas = c(0, 0.3, 0.5, 0.7, 1),
  n = 200L,
  p_values = c(500L, 1000L, 2000L, 3000L),
  landmark_proportions = c(0.30, 0.40, 0.50, 0.75, 0.90),
  top_m = c(3L, 10L, 20L),
  repetitions = 100L,
  rho = 0.8,
  tau = 1e-10,
  diagonal_tolerance = 1e-8,
  solver = "CLARABEL",
  workers = 1L,
  seed = 20260727L
)

check_simulation_packages()
check_cvxr_solver(config$solver)

result_dir <- file.path(
  script_dir,
  "results",
  paste0("nystrom_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

run_nystrom_job <- function(job, config) {
  example <- as.integer(job$example)
  kappa <- config$kappas[[as.integer(job$kappa_index)]]
  p <- config$p_values[[as.integer(job$p_index)]]
  replication <- as.integer(job$rep)
  data_seed <- config$seed + 100000000L * example +
    1000000L * as.integer(job$kappa_index) +
    10000L * as.integer(job$p_index) + replication

  tryCatch({
    generated <- measure_call(simulate_mixed_hermite_example(
      example = example,
      n = config$n,
      p = p,
      kappa = kappa,
      seed = data_seed,
      rho = config$rho
    ))
    data <- generated$value
    relevance_result <- measure_call(compute_relevance(data$X, data$Y))
    relevance <- relevance_result$value
    matrix_result <- measure_call(compute_full_bcor_matrix(data$X))
    Q <- matrix_result$value

    detail <- list()
    diagnostics <- list()
    errors <- list()
    detail_id <- diagnostic_id <- error_id <- 1L

    for (r_index in seq_along(config$landmark_proportions)) {
      r <- config$landmark_proportions[[r_index]]
      landmark_size <- max(2L, min(p, as.integer(floor(r * p))))
      landmark_seed <- data_seed + 50000000L + 1000L * r_index
      fit_result <- tryCatch({
        approximation_result <- measure_call(build_nystrom_approximation(
          Q = Q,
          landmark_size = landmark_size,
          seed = landmark_seed,
          tau = config$tau,
          diagonal_tolerance = config$diagonal_tolerance
        ))
        approximation <- approximation_result$value
        alpha <- choose_balanced_alpha(approximation$qbar_all, relevance)
        optimization_result <- measure_call(solve_nystrom_bcorqpfs(
          relevance = relevance,
          approximation = approximation,
          alpha = alpha,
          solver = config$solver
        ))
        fit <- optimization_result$value

        for (M in config$top_m) {
          record <- selection_record(fit$ranking, data$active, M)
          record$rep <- replication
          record$example <- example
          record$kappa <- kappa
          record$n <- config$n
          record$p <- p
          record$r <- r
          record$landmark_size <- landmark_size
          record$alpha <- alpha
          detail[[detail_id]] <- record[, c(
            "rep", "example", "kappa", "n", "p", "r", "landmark_size",
            "alpha", "M", "X1", "X2", "X3", "X5", "All", "selected"
          )]
          detail_id <- detail_id + 1L
        }

        diagnostics[[diagnostic_id]] <- data.frame(
          rep = replication,
          example = example,
          kappa = kappa,
          n = config$n,
          p = p,
          r = r,
          landmark_size = landmark_size,
          alpha = alpha,
          effective_rank = approximation$effective_rank,
          minimum_landmark_eigenvalue = approximation$minimum_eigenvalue,
          minimum_raw_diagonal_correction =
            approximation$minimum_raw_diagonal_correction,
          data_seconds = generated$seconds,
          relevance_seconds = relevance_result$seconds,
          matrix_seconds = matrix_result$seconds,
          approximation_seconds = approximation_result$seconds,
          optimization_seconds = optimization_result$seconds,
          solver = fit$solver,
          status = fit$status,
          stringsAsFactors = FALSE
        )
        diagnostic_id <- diagnostic_id + 1L
        TRUE
      }, error = function(e) e)

      if (inherits(fit_result, "error")) {
        errors[[error_id]] <- data.frame(
          rep = replication,
          example = example,
          kappa = kappa,
          n = config$n,
          p = p,
          r = r,
          stage = "nystrom_fit",
          message = conditionMessage(fit_result),
          stringsAsFactors = FALSE
        )
        error_id <- error_id + 1L
      }
    }

    list(
      detail = if (length(detail)) do.call(rbind, detail) else NULL,
      diagnostics = if (length(diagnostics)) do.call(rbind, diagnostics) else NULL,
      errors = if (length(errors)) do.call(rbind, errors) else NULL
    )
  }, error = function(e) {
    list(
      detail = NULL,
      diagnostics = NULL,
      errors = data.frame(
        rep = replication,
        example = example,
        kappa = kappa,
        n = config$n,
        p = p,
        r = NA_real_,
        stage = "data_or_bcor_matrix",
        message = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    )
  })
}

jobs <- expand.grid(
  rep = seq_len(config$repetitions),
  p_index = seq_along(config$p_values),
  kappa_index = seq_along(config$kappas),
  example = config$examples,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

write_checkpoint <- function(results, completed) {
  detail <- bind_result_component(results, "detail")
  diagnostics <- bind_result_component(results, "diagnostics")
  errors <- bind_result_component(results, "errors")
  if (!is.null(detail)) {
    utils::write.csv(detail, file.path(result_dir, "checkpoint_selection_detail.csv"), row.names = FALSE)
  }
  if (!is.null(diagnostics)) {
    utils::write.csv(diagnostics, file.path(result_dir, "checkpoint_diagnostics.csv"), row.names = FALSE)
  }
  if (!is.null(errors)) {
    utils::write.csv(errors, file.path(result_dir, "checkpoint_errors.csv"), row.names = FALSE)
  }
  utils::write.csv(
    data.frame(completed_jobs = completed, total_jobs = nrow(jobs)),
    file.path(result_dir, "checkpoint_status.csv"),
    row.names = FALSE
  )
}

started_at <- Sys.time()
results <- run_job_grid(
  jobs = jobs,
  worker_function = run_nystrom_job,
  config = config,
  workers = config$workers,
  functions_file = functions_file,
  checkpoint_function = write_checkpoint
)

detail <- bind_result_component(results, "detail")
diagnostics <- bind_result_component(results, "diagnostics")
errors <- bind_result_component(results, "errors")
if (is.null(detail) || !nrow(detail)) stop("No simulation fit completed.", call. = FALSE)

group_columns <- c("example", "kappa", "n", "p", "r", "landmark_size", "M")
group_key <- interaction(detail[group_columns], drop = TRUE, lex.order = TRUE)
summary_rows <- lapply(split(detail, group_key), function(part) {
  row <- part[1L, group_columns, drop = FALSE]
  row$successful_repetitions <- length(unique(part$rep))
  row$mean_alpha <- safe_mean(part$alpha)
  row$P_M_s_X1 <- safe_mean(part$X1)
  row$P_M_s_X2 <- safe_mean(part$X2)
  row$P_M_s_X3 <- safe_mean(part$X3)
  row$P_M_s_X5 <- safe_mean(part$X5)
  row$P_M_a <- safe_mean(part$All)
  row
})
summary_table <- do.call(rbind, summary_rows)
summary_table <- summary_table[order(
  summary_table$example,
  summary_table$kappa,
  summary_table$p,
  summary_table$r,
  summary_table$M
), , drop = FALSE]
rownames(summary_table) <- NULL

utils::write.csv(detail, file.path(result_dir, "selection_detail.csv"), row.names = FALSE)
utils::write.csv(diagnostics, file.path(result_dir, "diagnostics.csv"), row.names = FALSE)
utils::write.csv(summary_table, file.path(result_dir, "tables1_3_summary.csv"), row.names = FALSE)
if (!is.null(errors)) {
  utils::write.csv(errors, file.path(result_dir, "errors.csv"), row.names = FALSE)
}
saveRDS(
  list(config = config, jobs = jobs, detail = detail, diagnostics = diagnostics,
       summary = summary_table, errors = errors),
  file.path(result_dir, "all_results.rds")
)

elapsed <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
message(sprintf("Completed in %.2f hours.", elapsed / 3600))
message("Results: ", normalizePath(result_dir, winslash = "/"))
