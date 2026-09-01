# Simulation Code

`simulation_functions.R` contains data generation, bcorQPFS, mRMR, HSIC Lasso, evaluation, and execution functions.

`01_nystrom_simulation.R` runs the Nyström experiments for Tables 1--3 with the manuscript parameter grid.

`02_method_comparison.R` compares exact bcorQPFS, mRMR, and HSIC Lasso for Tables 4--5.

`response/03_plot_selection_frequency_by_r.R` creates the selection-frequency figure used in the response.

Required R packages are `Ball`, `CVXR`, and `glmnet`.

The response figure also uses `dplyr`, `tidyr`, `stringr`, `ggplot2`, and `scales`.

Open either numbered script in RStudio and run the complete file. Default parameters are defined in the `config` list. Results are saved in a timestamped folder under `results`.
