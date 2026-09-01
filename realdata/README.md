# ADNI Real-Data Analysis

This folder contains the preprocessing and feature-selection code for the ADNI analysis.

## Preprocessing

- `lm_adni1-3.R`: Combines ADNI1, ADNI2, and ADNI3, adjusts the hippocampal measurements for age, sex, APOE-4, and five genetic principal components, and saves `left_res_ori.Rdata` and `right_res_ori.Rdata`.

## Data Layout

Place the downloaded ADNI files and prepared genotype and hippocampal files under `code_V1/realdata/data/raw`:

```text
data/raw/
├── Hippocampus_all/
│   ├── sample information/
│   │   ├── ADNI1_hip_clinical_information.csv
│   │   ├── ADNI1_SNP_ID.csv
│   │   ├── ADNI1_hip_ID.csv
│   │   └── ADNI1Top20Pcs.txt
│   └── hipp_combine_information/
│       ├── ADNI1_MATCH_id23.Rdata
│       ├── left_hipp.csv
│       └── right_hipp.csv
├── ADNI2-3/
│   ├── ADNIMERGE_25Jan2024.csv
│   ├── hipp_information_radial/
│   │   ├── left_hipp_adni2.Rdata
│   │   ├── right_hipp_adni2.Rdata
│   │   ├── left_hipp_adni3.Rdata
│   │   └── right_hipp_adni3.Rdata
│   ├── ADNI2_SNP/
│   │   ├── ADNI2_MATCH_id13.Rdata
│   │   └── ADNI2_hipMADNI1.csv
│   └── ADNI3_SNP/
│       ├── ADNI3_MATCH_id12.Rdata
│       └── ADNI3_SamplehippMADNI1.csv
└── ADNI1-3/
    ├── left_G_21.Rdata
    └── right_G_21.Rdata
```

The preprocessing outputs are saved in `code/data`. Subject ID files are saved in `code/data/hipp_id`. The analysis scripts `01`--`04` read their inputs from `code/data`.

## Function Files

- `bcor_screen_select_functions.R`: Provides the functions used for BCor screening.
- `hsic_lasso_functions.R`: Provides the functions used for HSIC Lasso variable selection.

## Analysis Scripts

- `01_bcor_screening.R`: Performs BCor screening for the left and right hippocampal outcomes.
- `02_bcorqpfs_selection.R`: Runs bcorQPFS on the screened variables and selects the top SNPs.
- `03_hsic_screening.R`: Performs dHSIC screening for the left and right hippocampal outcomes.
- `04_hsic_lasso_selection.R`: Runs HSIC Lasso on the screened variables and selects the top SNPs.
