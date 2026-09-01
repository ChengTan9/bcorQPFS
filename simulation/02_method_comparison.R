# Run the exact bcorQPFS, mRMR, and HSIC Lasso comparison in Tables 4--5.

candidate_dirs <- c(
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    dirname(rstudioapi::getSourceEditorContext()$path)
  } else {
    character(0)
  },
  getwd(),
  file.path(getwd(), "code_V1", "simulation")
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
  n_values = c(200L, 500L),
  p_values = c(200L, 500L, 1000L),
  top_m = 10L,
  repetitions = 100L,
  rho = 0.8,
  solver = "CLARABEL",
  hsic_nlambda = 120L,
  workers = 1L,
  seed = 20260728L,
  method_labels = c(
    bcorqpfs = "bcorQPFS",
    mrmr = "mRMR",
    hsiclasso = "HSICLasso"
  )
)

check_simulation_packages(include_hsic = TRUE)
check_cvxr_solver(config$solver)

result_dir <- file.path(
  script_dir,
  "results",
  paste0("method_comparison_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

method_labels <- config$method_labels

run_comparison_job <- function(job, config) {
  example <- as.integer(job$example)
  kappa <- config$kappas[[as.integer(job$kappa_index)]]
  n <- config$n_values[[as.integer(job$n_index)]]
  p <- config$p_values[[as.integer(job$p_index)]]
  replication <- as.integer(job$rep)
  data_seed <- config$seed + 100000000L * example +
    1000000L * as.integer(job$kappa_index) +
    100000L * as.integer(job$n_index) +
    10000L * as.integer(job$p_index) + replication

  tryCatch({
    generated <- measure_call(simulate_mixed_hermite_example(
      example = example,
      n = n,
      p = p,
      kappa = kappa,
      seed = data_seed,
      rho = config$rho
    ))
    data <- generated$value
    detail <- list()
    diagnostics <- list()
    errors <- list()
    detail_id <- diagnostic_id <- error_id <- 1L

    add_result <- function(method_id, ranking, elapsed, alpha = NA_real_,
                           solver = NA_character_, status = "completed",
                           lambda = NA_real_, nonzero = NA_integer_) {
      record <- selection_record(ranking, data$active, config$top_m)
      record$rep <- replication
      record$example <- example
      record$kappa <- kappa
      record$n <- n
      record$p <- p
      record$method_id <- method_id
      record$method <- unname(config$method_labels[[method_id]])
      detail[[detail_id]] <<- record[, c(
        "rep", "example", "kappa", "n", "p", "method_id", "method",
        "M", "X1", "X2", "X3", "X5", "All", "selected"
      )]
      detail_id <<- detail_id + 1L

      diagnostics[[diagnostic_id]] <<- data.frame(
        rep = replication,
        example = example,
        kappa = kappa,
        n = n,
        p = p,
        method_id = method_id,
        method = unname(config$method_labels[[method_id]]),
        alpha = alpha,
        lambda = lambda,
        nonzero = nonzero,
        data_seconds = generated$seconds,
        method_seconds = elapsed,
        solver = solver,
        status = status,
        stringsAsFactors = FALSE
      )
      diagnostic_id <<- diagnostic_id + 1L
    }

    bcor_result <- tryCatch({
      timed <- measure_call({
        relevance <- compute_relevance(data$X, data$Y)
        Q <- compute_full_bcor_matrix(data$X)
        matrix_summary <- summarize_bcor_matrix(Q)
        alpha <- choose_balanced_alpha(matrix_summary$qbar_all, relevance)
        fit <- solve_exact_bcorqpfs(relevance, Q, alpha, config$solver)
        list(fit = fit, alpha = alpha)
      })
      add_result(
        method_id = "bcorqpfs",
        ranking = timed$value$fit$ranking,
        elapsed = timed$seconds,
        alpha = timed$value$alpha,
        solver = timed$value$fit$solver,
        status = timed$value$fit$status
      )
      TRUE
    }, error = function(e) e)
    if (inherits(bcor_result, "error")) {
      errors[[error_id]] <- data.frame(
        rep = replication, example = example, kappa = kappa, n = n, p = p,
        method = "bcorQPFS", message = conditionMessage(bcor_result),
        stringsAsFactors = FALSE
      )
      error_id <- error_id + 1L
    }

    mrmr_result <- tryCatch({
      timed <- measure_call(rank_mrmr(data$X, data$Y, config$top_m))
      add_result("mrmr", timed$value, timed$seconds)
      TRUE
    }, error = function(e) e)
    if (inherits(mrmr_result, "error")) {
      errors[[error_id]] <- data.frame(
        rep = replication, example = example, kappa = kappa, n = n, p = p,
        method = "mRMR", message = conditionMessage(mrmr_result),
        stringsAsFactors = FALSE
      )
      error_id <- error_id + 1L
    }

    hsic_result <- tryCatch({
      timed <- measure_call(rank_hsic_lasso(
        data$X, data$Y, config$top_m, config$hsic_nlambda
      ))
      add_result(
        method_id = "hsiclasso",
        ranking = timed$value$ranking,
        elapsed = timed$seconds,
        lambda = timed$value$lambda,
        nonzero = timed$value$nonzero
      )
      TRUE
    }, error = function(e) e)
    if (inherits(hsic_result, "error")) {
      errors[[error_id]] <- data.frame(
        rep = replication, example = example, kappa = kappa, n = n, p = p,
        method = "HSICLasso", message = conditionMessage(hsic_result),
        stringsAsFactors = FALSE
      )
      error_id <- error_id + 1L
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
        n = n,
        p = p,
        method = "data_generation",
        message = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    )
  })
}

jobs <- expand.grid(
  rep = seq_len(config$repetitions),
  p_index = seq_along(config$p_values),
  n_index = seq_along(config$n_values),
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
  worker_function = run_comparison_job,
  config = config,
  workers = config$workers,
  functions_file = functions_file,
  checkpoint_function = write_checkpoint
)

detail <- bind_result_component(results, "detail")
diagnostics <- bind_result_component(results, "diagnostics")
errors <- bind_result_component(results, "errors")
if (is.null(detail) || !nrow(detail)) stop("No comparison fit completed.", call. = FALSE)

group_columns <- c("example", "kappa", "n", "p", "method_id", "method", "M")
group_key <- interaction(detail[group_columns], drop = TRUE, lex.order = TRUE)
summary_rows <- lapply(split(detail, group_key), function(part) {
  row <- part[1L, group_columns, drop = FALSE]
  row$successful_repetitions <- length(unique(part$rep))
  row$P_M_s_X1 <- safe_mean(part$X1)
  row$P_M_s_X2 <- safe_mean(part$X2)
  row$P_M_s_X3 <- safe_mean(part$X3)
  row$P_M_s_X5 <- safe_mean(part$X5)
  row$P_M_a <- safe_mean(part$All)
  row
})
summary_table <- do.call(rbind, summary_rows)
summary_table <- summary_table[order(
  summary_table$n,
  summary_table$p,
  summary_table$example,
  summary_table$kappa,
  match(summary_table$method_id, names(method_labels))
), , drop = FALSE]
rownames(summary_table) <- NULL

utils::write.csv(detail, file.path(result_dir, "selection_detail.csv"), row.names = FALSE)
utils::write.csv(diagnostics, file.path(result_dir, "diagnostics.csv"), row.names = FALSE)
utils::write.csv(summary_table, file.path(result_dir, "tables4_5_summary.csv"), row.names = FALSE)
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
