# Plot selection frequencies across Nystrom landmark proportions.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

find_script_directory <- function() {
  local_candidates <- c(
    ".",
    file.path("code_V1", "simulation", "response")
  )
  local_candidates <- local_candidates[file.exists(
    file.path(local_candidates, "03_plot_selection_frequency_by_r.R")
  )]
  if (length(local_candidates)) return(local_candidates[[1L]])

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    path <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(path)) return(dirname(path))
  }
  for (frame in rev(sys.frames())) {
    if (!is.null(frame$ofile)) {
      return(dirname(frame$ofile))
    }
  }
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_argument)) {
    path <- sub("^--file=", "", file_argument[[1L]])
    return(dirname(path))
  }
  stop("Cannot locate the script directory.", call. = FALSE)
}

# Set the figure parameters.
input_csv <- NULL
target_M <- 10L
target_kappa <- 0.5
target_n <- 200L
target_p <- 1000L
target_examples <- 1:3
target_r <- c(0.30, 0.40, 0.50, 0.75, 0.90)
variables_per_panel <- 7L

script_dir <- find_script_directory()
output_dir <- file.path(script_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Locate a compatible selection-detail file.
find_selection_file <- function(
    input_csv, script_dir, target_M, target_kappa, target_n, target_p,
    target_examples, target_r) {
  if (!is.null(input_csv)) {
    if (!file.exists(input_csv)) stop("The specified input file does not exist.", call. = FALSE)
    return(input_csv)
  }

  search_roots <- unique(c(
    file.path(script_dir, "..", "results"),
    file.path(script_dir, "..", "..", "simulation_V0", "results"),
    file.path("code_V1", "simulation", "results"),
    file.path("code_V1", "simulation_V0", "results"),
    "results"
  ))
  search_roots <- search_roots[dir.exists(search_roots)]
  candidates <- unique(unlist(lapply(
    search_roots,
    function(root) list.files(
      root,
      pattern = "^selection_detail\\.csv$",
      recursive = TRUE,
      full.names = TRUE
    )
  )))
  if (!length(candidates)) stop("Cannot find selection_detail.csv.", call. = FALSE)

  candidates <- candidates[order(file.info(candidates)$mtime, decreasing = TRUE)]
  required <- c("rep", "example", "kappa", "n", "p", "r", "M", "selected")
  for (candidate in candidates) {
    columns <- tryCatch(
      names(utils::read.csv(candidate, nrows = 0L, check.names = FALSE)),
      error = function(e) character(0)
    )
    if (!all(required %in% columns)) next
    candidate_data <- tryCatch(
      utils::read.csv(candidate, check.names = FALSE)[, required, drop = FALSE],
      error = function(e) NULL
    )
    if (is.null(candidate_data)) next
    candidate_data <- candidate_data %>%
      filter(
        M == target_M,
        near(kappa, target_kappa),
        n == target_n,
        p == target_p,
        example %in% target_examples
      )
    available <- expand_grid(example = target_examples, r = target_r) %>%
      anti_join(
        distinct(candidate_data, example, r),
        by = c("example", "r")
    )
    if (!nrow(available)) {
      return(candidate)
    }
  }
  stop("No compatible selection-detail file was found.", call. = FALSE)
}

csv_file <- find_selection_file(
  input_csv = input_csv,
  script_dir = script_dir,
  target_M = target_M,
  target_kappa = target_kappa,
  target_n = target_n,
  target_p = target_p,
  target_examples = target_examples,
  target_r = target_r
)
message("Input: ", csv_file)

raw <- utils::read.csv(csv_file, check.names = FALSE)
required_columns <- c("rep", "example", "kappa", "n", "p", "r", "M", "selected")
missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns)) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

target_r <- sort(unique(target_r))
fixed_rows <- raw %>%
  filter(
    M == target_M,
    near(kappa, target_kappa),
    n == target_n,
    p == target_p,
    example %in% target_examples
  )

# Match stored landmark proportions to the requested grid.
r_index <- vapply(fixed_rows$r, function(value) {
  match_index <- which(near(value, target_r))
  if (length(match_index)) match_index[[1L]] else NA_integer_
}, integer(1))
selected_rows <- fixed_rows %>%
  mutate(r = target_r[r_index]) %>%
  filter(!is.na(r))

expected_panels <- expand_grid(example = target_examples, r = target_r)
missing_panels <- expected_panels %>%
  anti_join(distinct(selected_rows, example, r), by = c("example", "r"))
if (nrow(missing_panels)) {
  stop("The requested example and r grid is incomplete.", call. = FALSE)
}

duplicate_runs <- selected_rows %>%
  count(example, r, rep) %>%
  filter(n > 1L)
if (nrow(duplicate_runs)) stop("Duplicate simulation rows were found.", call. = FALSE)

# Calculate panel-specific selection frequencies.
denominators <- selected_rows %>%
  count(example, r, name = "n_rep")
selection_long <- selected_rows %>%
  select(example, r, rep, selected) %>%
  separate_longer_delim(selected, delim = ";") %>%
  transmute(example, r, rep, variable = str_trim(selected)) %>%
  distinct(example, r, rep, variable)

if (any(!str_detect(selection_long$variable, "^X[0-9]+$"))) {
  stop("Invalid variable names were found.", call. = FALSE)
}
selection_size <- selection_long %>% count(example, r, rep)
if (any(selection_size$n != target_M)) {
  stop("A selected set has an unexpected size.", call. = FALSE)
}

selection_counts <- selection_long %>%
  count(example, r, variable, name = "selected_count")
selection_frequency <- expand_grid(
  example = target_examples,
  r = target_r,
  variable = paste0("X", seq_len(target_p))
) %>%
  left_join(selection_counts, by = c("example", "r", "variable")) %>%
  mutate(selected_count = replace_na(selected_count, 0L)) %>%
  left_join(denominators, by = c("example", "r")) %>%
  mutate(frequency = selected_count / n_rep)

active_variables <- tibble(
  example = rep(1:3, each = 3L),
  variable = c("X2", "X3", "X5", "X1", "X2", "X5", "X1", "X2", "X3")
)

# Select and order variables within each panel.
plot_data <- selection_frequency %>%
  left_join(
    active_variables %>% mutate(is_active = TRUE),
    by = c("example", "variable")
  ) %>%
  mutate(
    is_active = replace_na(is_active, FALSE),
    role = factor(
      if_else(is_active, "Active variable", "Other variable"),
      levels = c("Active variable", "Other variable")
    ),
    variable_number = as.integer(str_remove(variable, "^X"))
  ) %>%
  arrange(r, example, desc(frequency), desc(is_active), variable_number) %>%
  group_by(example, r) %>%
  slice_head(n = variables_per_panel) %>%
  mutate(rank = row_number(), x_position = rank, is_top_three = rank <= 3L) %>%
  ungroup() %>%
  mutate(
    example_label = factor(
      example,
      levels = target_examples,
      labels = paste("Example", target_examples)
    ),
    r_label = factor(r, levels = target_r, labels = paste("r =", target_r))
  )

label_y <- -0.038
highlight_bottom <- -0.072
frequency_plot <- ggplot(plot_data, aes(x = x_position)) +
  geom_tile(
    aes(y = frequency / 2, height = frequency, fill = role),
    width = 0.62,
    colour = NA
  ) +
  geom_hline(yintercept = 0, colour = "#7F8C8D", linewidth = 0.45) +
  geom_text(
    aes(y = label_y, label = variable),
    size = 3.8,
    colour = "#39434A",
    vjust = 0.5
  ) +
  geom_tile(
    data = plot_data %>% filter(is_top_three),
    aes(
      x = x_position,
      y = (frequency + highlight_bottom) / 2,
      width = 0.74,
      height = frequency - highlight_bottom,
      colour = "Three highest-frequency variables"
    ),
    inherit.aes = FALSE,
    fill = NA,
    linewidth = 1.15
  ) +
  facet_grid(rows = vars(r_label), cols = vars(example_label), drop = FALSE) +
  scale_x_continuous(
    limits = c(0.5, variables_per_panel + 0.5),
    breaks = NULL,
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    labels = label_percent(accuracy = 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0.015))
  ) +
  coord_cartesian(ylim = c(0, 1.02), clip = "off") +
  scale_fill_manual(
    values = c("Active variable" = "#1B9E77", "Other variable" = "#A7B6C2"),
    drop = FALSE
  ) +
  scale_colour_manual(
    values = c("Three highest-frequency variables" = "#D95F02")
  ) +
  labs(x = NULL, y = "Selection frequency", fill = NULL, colour = "Colored border") +
  theme_minimal(base_size = 11.5, base_family = "sans") +
  theme(
    strip.text = element_text(face = "bold", colour = "#22313F"),
    strip.text.y = element_text(angle = 0),
    strip.background = element_rect(fill = "#EDF2F4", colour = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = grid::unit(1.2, "lines"),
    panel.spacing.y = grid::unit(2.8, "lines"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.justification = "left",
    plot.margin = margin(12, 16, 34, 12)
  )

figure_height <- 5.8 + 3.2 * (length(target_r) - 1L)
ggsave(
  file.path(output_dir, "selection_frequency.pdf"),
  plot = frequency_plot,
  width = 14,
  height = figure_height,
  units = "in",
  device = "pdf"
)
grDevices::png(
  filename = file.path(output_dir, "selection_frequency.png"),
  width = 14,
  height = figure_height,
  units = "in",
  res = 320,
  bg = "white"
)
print(frequency_plot)
grDevices::dev.off()
utils::write.csv(
  plot_data %>% select(
    example, r, variable, rank, role, is_top_three,
    selected_count, n_rep, frequency
  ),
  file.path(output_dir, "selection_frequency_data.csv"),
  row.names = FALSE
)

message("Figure directory: ", output_dir)
