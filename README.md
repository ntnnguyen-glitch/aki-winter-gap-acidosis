# Winter gap and in-hospital mortality in acute kidney injury (MIMIC-IV and eICU-CRD)

R code for the study "Inadequate respiratory compensation for metabolic acidosis and in-hospital mortality in critically ill patients with acute kidney injury: a two-database retrospective cohort study".

The repository contains code only. It contains no patient-level data, no intermediate data files and no results derived from individual records.

## Data access

Both databases are available from PhysioNet to credentialed users who have completed human-subjects research training and signed the data use agreement:

- MIMIC-IV version 3.1: https://doi.org/10.13026/kpb9-mt58
- eICU Collaborative Research Database version 2.0: https://doi.org/10.13026/C2WM1R

Download the compressed CSV files and set `DATA_ROOT` in `config.R` to the folder that contains `mimic-iv-3.1/` (with `hosp/` and `icu/`) and `data eICU/`.

## Software

R 4.6.1 with: DBI, duckdb 2.0.0, data.table 1.18.6.1, survival 3.8.6, splines, sandwich, mice 3.19.0, WeightIt 2.1.0, MatchIt 4.8.2, EValue 4.1.4, cmprsk 2.2.12, ggplot2 4.0.4, survminer 0.5.2, patchwork, cowplot, jsonlite. The full session information is in `sessionInfo.txt`.

## Run order

Run the scripts from the repository root, in this order:

| Script | Purpose | Main output |
|---|---|---|
| `R/01_extract_base.R` | Reads the PhysioNet CSV files into a local DuckDB database (stays, laboratory values, charted vital signs, procedures, infusions, diagnoses; eICU-CRD equivalents) | `w.duckdb` |
| `R/02_extract_windows.R` | Extracts charted ventilation signals, procedures and eICU-CRD vital signs within 30 h of ICU admission for all stays with a PaCO2 measurement | tables added to `w.duckdb` |
| `R/03_cohort.R` | Applies the eligibility criteria (KDIGO creatinine AKI within 48 h, arterial specimen, metabolic acidosis, no prior kidney replacement therapy), classifies ventilatory support at the index gas and derives covariates and outcomes | `D2.rds`, `fl2.rds` |
| `R/04_analysis.R` | Multiple imputation (40 datasets), logistic, Cox, Fine-Gray, propensity score and E-value analyses, sensitivity and subgroup analyses, with hospital-clustered standard errors | `R2.rds`, `P0.rds` |
| `R/05_spline_tests.R` | Natural cubic spline curve and likelihood ratio tests (overall, non-linearity, gap categories) combined across imputations with the D2 method | updates `R2.rds` |
| `R/06_tables_figures.R` | Baseline tables and Figures 1-4 | `output/` |

Intermediate `.rds` files and `w.duckdb` contain patient-level data and must not be shared or committed; `.gitignore` excludes them.

## Notes

- eICU-CRD files are read with `ignore_errors=true`; rows that cannot be parsed are skipped by DuckDB.
- The random seed is fixed (`set.seed(2026)`, mice seed 2026).
- Variable definitions, item identifiers and time windows are listed in Additional file 2, Table S1 of the article.

## Citation

If you use this code, please cite the article and the archived release: https://doi.org/[ZENODO DOI]

## Licence

MIT (see `LICENSE`).
