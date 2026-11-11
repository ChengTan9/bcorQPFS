# bcorQPFS

R code for **Quadratic Programming Variable Selection Based on Ball Correlation**.

## Repository contents

```text
realdata/     ADNI preprocessing, screening, selection and landmark sensitivity analysis
simulation/   Simulation studies, Example 1 timing and selection-frequency plots
README.md     Setup and execution instructions
```

Only source code and the public SNP annotation lookup are included. Scripts create their output directories when run. Historical code, local run results and internal verification records are not part of this repository.

## Requirements

The code was checked with R 4.5.1, Ball 1.3.13, CVXR 1.8.1 and clarabel 0.11.2. Install the required packages before running:

```r
install.packages(c("Ball", "CVXR", "clarabel", "glmnet", "Matrix",
                   "openxlsx", "pbapply", "dHSIC", "hdf5r", "progress"))
# Additional packages for simulation/05_plot_selection_frequency_by_r.R:
install.packages(c("dplyr", "tidyr", "stringr", "ggplot2", "scales"))
```

The plotting script uses `tidyr::separate_longer_delim` (tidyr >= 1.3.0). The simulation scripts use the CLARABEL solver through CVXR. `parallel`, `stats`, `utils` and the graphics packages supplied with R are also used. No script installs or updates packages automatically.

Commands below are run from the repository root. Entry scripts locate their supporting files relative to their own locations and can also be called using an absolute script path.

## Simulation studies

```sh
Rscript --vanilla simulation/01_nystrom_simulation.R
Rscript --vanilla simulation/02_method_comparison.R
Rscript --vanilla simulation/03_example1_timing.R
Rscript --vanilla simulation/04_summarize_timing.R
Rscript --vanilla simulation/05_plot_selection_frequency_by_r.R
```

| Script | Experiment and output |
|---|---|
| `01_nystrom_simulation.R` | Examples 1--3; n=200; p=500,1000,2000,3000; kappa=0,0.3,0.5,0.7,1; r=0.30,0.40,0.50,0.75,0.90; M=3,10,20; 100 replications. Saves `simulation/results/nystrom_TIMESTAMP/selection_detail.csv`, summaries and diagnostics. |
| `02_method_comparison.R` | Complete QP, mRMR and HSIC Lasso; n=200,500; p=200,500,1000; the same kappa grid; fixed Top-10; 100 replications. Saves `simulation/results/method_comparison_TIMESTAMP/`. |
| `03_example1_timing.R` | Example 1, kappa=0.5; n=200,500; p=200,500; r=0.3,0.5,0.7,1; one fit per setting. Saves `results/example1_timing/run_TIMESTAMP/example1_timing_full_grid.csv` and session information. |
| `04_summarize_timing.R` | Reads the latest completed full timing grid and saves times, percentages relative to complete QP and summary tables in `results/summaries/example1_timing/RUN_NAME/`. |
| `05_plot_selection_frequency_by_r.R` | Reads `selection_detail.csv` from `simulation/results/` matching the target parameters at the top of the script. Saves PDF, PNG and plotting data to `simulation/figures/`. |

The main simulation defaults and numerical routines are retained from the previous code version. Simulations do not use permutation screening. Script `01` reuses a complete dependence matrix to study selection behavior; computational savings from evaluating only the required landmark entries are assessed separately by `03`. Within each n,p timing setting, all fits use the same generated data. Timing includes response relevance, dependence evaluation, factor construction, and QP formulation and solution; data generation is excluded.

For a small timing check:

```sh
Rscript --vanilla simulation/03_example1_timing.R --smoke
Rscript --vanilla simulation/04_summarize_timing.R --input="results/example1_timing/smoke_TIMESTAMP"
```

Replace `smoke_TIMESTAMP` with the directory actually created. Smoke results are excluded from the default full-grid summaries. `04` can also read an explicitly specified complete run via `--input=RUN_DIRECTORY`. All full experiments use the seeds recorded in their scripts; the smoke option changes only the size of the small check.

## ADNI analysis

Participant-level ADNI data are controlled-access and are not distributed here. Authorized users must supply aligned genotype and residual hippocampal phenotype files. Set `BCORQPFS_DATA_DIR` to their directory, for example:

```r
Sys.setenv(BCORQPFS_DATA_DIR = "/path/to/authorized/processed_data")
```

For separate `Rscript` commands, set the variable in the calling shell. In PowerShell:

```powershell
$env:BCORQPFS_DATA_DIR = "D:/authorized/processed_data"
```

Required processed inputs are:

- `ADNI_1_3_COMBINE1215.Rdata`: object `ADNI_COMBINE`, with participant ID in the first column and SNP dosage columns thereafter.
- `left_res_ori.Rdata`: left residual phenotype matrix `left_res_ori`.
- `right_res_ori.Rdata`: right residual phenotype matrix `right_res_ori`.

The subject identities and order must agree across inputs. New screening outputs are written to `data/processed/` and take precedence over external files of the same name. Without the environment variable, inputs are read from `data/processed/`. To reproduce preprocessing from authorized raw data, set `BCORQPFS_RAW_DATA_DIR` and run `realdata/lm_adni1-3.R`; its source specifies the required cohort-specific raw files and directory layout.

### Ball-correlation screening and bcorQPFS

```sh
Rscript --vanilla realdata/01_bcor_screening.R
Rscript --vanilla realdata/02_bcorqpfs_selection.R
Rscript --vanilla realdata/05_adni_exact_and_landmarks.R
Rscript --vanilla realdata/07_summarize_adni.R
```

| Producer | Output consumed by the next script |
|---|---|
| `01_bcor_screening.R` | `data/processed/ADNI_left_screening_selected.RData` and the right-side equivalent, containing `selected_X`; companion `ADNI_SIDE_screening_metadata.rds` files record variable names, original indices and screening information. |
| `02_bcorqpfs_selection.R` | `results/adni_single_run/TIMESTAMP/BcorQPFS_SIDE_full_result.rds`; the same run's `cache/SIDE_factor_inputs_sS.rds` and `cache/SIDE_DY_dist.rds`. |
| `05_adni_exact_and_landmarks.R` | `results/adni_stability/run_TIMESTAMP/SIDE_exact_qp_top20.csv`, `SIDE_landmark_frequency_ranking.csv`, `SIDE_side_summary.csv`, per-draw rankings and diagnostics. |
| `07_summarize_adni.R` | `results/summaries/adni/RUN_NAME/complete_qp_top20_both_sides.csv`, `frequency_top20_plus_reported_both_sides.csv` and `stability_summary.csv`. |

Here `SIDE` is `left` or `right`, and `S` is the number of landmarks. Script `02` retains the original r=0.3, relative spectral threshold 1e-10 and Top-30 export; the manuscript displays the first 20. Script `05` reads the most recent completed bilateral `02` run and uses the recorded screened inputs and caches. It fits the complete empirical QP and performs 100 additional r=0.3 landmark draws with fixed seeds. The saved single draw is a separate reference (run 0) and is excluded from repeat-selection proportions. `05 --smoke` performs only two additional draws; its output directory is excluded from the default summaries. Use `07 --input=RUN_DIRECTORY` to summarize an explicitly selected run, including a smoke run.

`realdata/variant_annotations_grch37.csv` contains the reviewed public SNP annotations used by the manuscript. The summary script joins annotations by chromosome and position; newly selected SNPs outside this lookup are marked `not_in_reviewed_annotation_set` rather than assigned an inferred gene.

For reproducing an already saved analysis, `BCORQPFS_BCOR_ARCHIVE_DIR` may explicitly point to a complete prior bilateral fit directory with its caches; the corresponding authorized input files must also be available. Leave this variable unset for the normal new-data workflow. Without a completed `02` run or an explicit prior fit, `05` stops with an input error.

### HSIC screening and HSIC Lasso

```sh
Rscript --vanilla realdata/03_hsic_screening.R
Rscript --vanilla realdata/04_hsic_lasso_selection.R
```

`03` saves `data/processed/ADNI_dHSIC_selected_vs_SIDE_res_ori.RData`, containing `screening_result` with `X_selected`, `selected_index`, `selected_names`, `top_k` and `n`. `04` reads these exact files and saves `results/hsic_lasso_selection/TIMESTAMP/HSICLasso_SIDE_top30.xlsx`. The original full-kernel settings B=0, M=1 and 30-variable path budget are retained. This workflow uses its own HSIC-screened candidate set.

## Outputs and reproducibility

Generated `data/`, `results/`, `simulation/results/` and `simulation/figures/` directories are local runtime outputs and should be excluded from a source-only upload, especially participant-level files and caches. They are absent from this distribution. Default seeds, selection budgets and numerical settings are declared in the scripts. Reading and writing paths have been checked with a separate synthetic-data workflow; the full simulation grid is not rerun by a smoke check.
