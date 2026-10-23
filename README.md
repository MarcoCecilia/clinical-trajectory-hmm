[README.md](https://github.com/user-attachments/files/33881020/README.md)
# Latent-state modelling of inpatient biomarker trajectories with Hidden Markov Models

Code accompanying the Master's thesis in Mathematical Engineering (Politecnico di Milano)
by **Marco Cecilia**.

The repository implements a four-stage workflow that turns sparse, irregular in-hospital
blood draws into regular daily trajectories, fits a multivariate-normal Hidden Markov Model
(HMM) on them, and uses the decoded latent states for phenotype-based outcome analyses.

| Stage | File | Language | Role |
|---|---|---|---|
| 1 | `Final_pipeline.ipynb` | Python | Extraction, filtering and merging of the raw sources; long-format longitudinal dataset |
| 2 | `Imputation.Rmd` | R | Per-patient GAM reconstruction on a fixed daily grid; variance-stabilising transform |
| 3 | `hmm_training_tesi.Rmd` | R | Response-set selection, selection of the number of states, final HMM fit, diagnostics |
| 4 | `hmm_posthoc.Rmd` | R | Stayers/switchers, depth-based phenotypes, classification anchors, outcome models, alert rule |

The stages must be run **in this order**: each one consumes files written by the previous one.

---

## Data availability

The analyses use real-world clinical data from IRCCS Ospedale San Raffaele (Milan).
**No patient-level data are included in this repository**, and none of the intermediate
CSV/RDS files listed below may be committed. The code is published so that the analysis
can be inspected and re-run by anyone with authorised access to the source extracts.

> _[Add here the data-access statement / ethics approval reference required by the hospital.]_

Recommended `.gitignore` entries:

```gitignore
*.parquet
*.csv
*.rds
figures/
figures_ch4/
tables_ch4/
*_cache/
*_files/
```

---

## Repository structure

```
.
├── Final_pipeline.ipynb        # Stage 1 (Python)
├── Imputation.Rmd              # Stage 2 (R)
├── hmm_training_tesi.Rmd       # Stage 3 (R)
├── hmm_posthoc.Rmd             # Stage 4 (R)
├── README.md
└── .gitignore
```

All notebooks read and write files in their **working directory**. The simplest set-up is to
run all four from the same folder.

---

## Software requirements

### Python (Stage 1)

- Python 3
- `pandas`, `numpy`
- `scikit-learn` (`IterativeImputer`, `RandomForestRegressor`)
- a Parquet engine for `pandas.read_parquet` (`pyarrow` or `fastparquet`)
- Jupyter (to execute the notebook)

```bash
pip install pandas numpy scikit-learn pyarrow jupyter
```

### R (Stages 2–4)

The reference run used **R 4.5.3** (x86_64-apple-darwin20, macOS Sonoma). Package versions
of the reference run:

| Package | Version | Used in |
|---|---|---|
| `tidyverse` | 2.0.0 | 2, 3, 4 |
| `mgcv` | 1.9-4 | 2 (also a dependency of `hmmTMB`) |
| `splines` | base | 2 |
| `patchwork` | 1.3.2 | 2, 3, 4 |
| `hmmTMB` | 1.1.2 | 3, 4 |
| `TMB` | 1.9.21 | 3, 4 (dependency of `hmmTMB`) |
| `mvtnorm` | 1.3-7 | 3 |
| `ggalluvial` | 0.12.6 | 3, 4 |
| `DepthProc` | 2.1.6 | 4 |
| `aplpack` | 1.3.5 | 4 |
| `MASS` | 7.3-65 | 4 |
| `broom` | 1.0.12 | 4 |
| `survival` | 3.8-6 | 4 |
| `survminer` | 0.5.2 | 4 |
| `rmarkdown`, `knitr` | 2.31, 1.51 | rendering |

```r
install.packages(c(
  "tidyverse", "mgcv", "patchwork", "hmmTMB", "mvtnorm", "ggalluvial",
  "DepthProc", "aplpack", "MASS", "broom", "survival", "survminer",
  "rmarkdown", "knitr", "tinytex"
))
```

PDF rendering uses `xelatex` (e.g. via `tinytex::install_tinytex()`). To pin the exact
versions above, `renv` is recommended (`renv::init()` followed by `renv::snapshot()`).

---

## Pipeline overview

```
labs_clean_source_copy.parquet ┐
01_clinical.parquet            ├─► [1] Final_pipeline.ipynb ─► df_longitudinal.csv
new_df_episodekey_deathdate.csv┘                            └► biomarker_reference_ranges.csv
                                                                     │
df_longitudinal.csv ─► [2] Imputation.Rmd  (run twice, see below)    │
   ├─ run A: HGB, CREA, ALT            ─► df_imputed_wide_log1p_HgbCreaAlt.csv
   │                                      bm_log1p_decisions_HgbCreaAlt.csv
   └─ run B: HGB, CREA, ALT, WBC, CRP  ─► df_imputed_wide_log1p_HgbCreaAltWbcCrp.csv
                                          bm_log1p_decisions_HgbCreaAltWbcCrp.csv
                                                                     │
both cohorts + biomarker_reference_ranges.csv ─► [3] hmm_training_tesi.Rmd
                                                  └► hmm_fit_objects_tesi4.rds, figures/
                                                                     │
hmm_fit_objects_tesi4.rds ─► [4] hmm_posthoc.Rmd ─► figures_ch4/, tables_ch4/
```

---

## Running the pipeline

### Stage 1 — `Final_pipeline.ipynb`

**Inputs** (same folder as the notebook):

| File | Content | Required columns |
|---|---|---|
| `labs_clean_source_copy.parquet` | laboratory results, one row per draw | `Subject_identifier, Date, Method_text, Measure, Numeric_value, Range_string` |
| `01_clinical.parquet` | admission registry | `CloudPatientId, startDate, endDate` + admission covariates |
| `new_df_episodekey_deathdate.csv` | mortality registry | `CloudPatientId, dateOfDeath` |

**Configuration** (first cells): `DATE_START = '2024-02-01'`, `DATE_END = '2026-01-31'`
(only admissions whose entire stay falls within the window are kept), the biomarker
whitelist `BLOOD_TESTS_TO_KEEP`, and `MICE_RANDOM_STATE = 42`.

> **Before running:** `OUTPUT_DIR` is hard-coded to an absolute path of the original
> computing environment. Set it to the working directory, e.g. `OUTPUT_DIR = BASE_DIR`.

Main steps: biomarker whitelist and unit harmonisation; per-patient selection of the last
admission entirely within the window; merge of the trusted date of death (the
`dateOfDeath` column of the clinical parquet is dropped before the merge); pivot to one
column per analysis; `day` = days since admission (day 0 = admission day); patient-level
imputation of admission covariates with `IterativeImputer` + random forests (dates and
outcomes are never imputed); parsing of laboratory reference ranges; ITA→ENG renaming.

**Outputs:** `df_longitudinal.csv`, `biomarker_reference_ranges.csv`.

```bash
jupyter nbconvert --to notebook --execute Final_pipeline.ipynb --output Final_pipeline_executed.ipynb
```

### Stage 2 — `Imputation.Rmd`

**Input:** `df_longitudinal.csv`.

Each (patient × biomarker) series is reconstructed on the fixed grid of days 0–14 with a
Gamma GAM with log link (`mgcv::gam`, cubic regression spline, REML smoothing), with a
day-0 anchor (and a day-14 anchor for short stays) and LOCF as fallback for series with
fewer than four valid observations. Skewed biomarkers (|skewness| > 1) are then
`log1p`-transformed; the decision table is exported for the back-transform used
downstream. Survival outcomes (`discharge_dead`, `discharge_mortality_90d`,
`days_death_from_admission`, `days_death_from_discharge`, `is_dead`) are derived here from
the date of death.

The notebook is run **twice**, changing only the `BIOMARKERS_TO_USE` whitelist, because
the patient-exclusion rule depends on the biomarkers in scope:

| Run | `BIOMARKERS_TO_USE` | Rename the outputs to |
|---|---|---|
| A (analysis cohort) | `c("HGB (g/dL)", "CREA (mg/dL)", "ALT (U/L)")` | `df_imputed_wide_log1p_HgbCreaAlt.csv`, `bm_log1p_decisions_HgbCreaAlt.csv` |
| B (CRP/WBC subcohort, response-set selection only) | `c("HGB (g/dL)", "CREA (mg/dL)", "ALT (U/L)", "WBC (10^9/L)", "CRP (mg/L)")` | `df_imputed_wide_log1p_HgbCreaAltWbcCrp.csv`, `bm_log1p_decisions_HgbCreaAltWbcCrp.csv` |

The notebook always writes `df_imputed_wide.csv`, `df_imputed_wide_log1p.csv`,
`bm_log1p_decisions.csv` and `df_spline_long.csv`; rename the two files required by
Stage 3 after each run, before starting the next one.

```bash
Rscript -e 'rmarkdown::render("Imputation.Rmd")'
mv df_imputed_wide_log1p.csv df_imputed_wide_log1p_HgbCreaAlt.csv
mv bm_log1p_decisions.csv    bm_log1p_decisions_HgbCreaAlt.csv
# edit BIOMARKERS_TO_USE, re-render, then rename with the _HgbCreaAltWbcCrp suffix
```

### Stage 3 — `hmm_training_tesi.Rmd`

**Inputs:** the four files from Stage 2 and `biomarker_reference_ranges.csv`.

| Section | Content |
|---|---|
| 0–2 | Configuration, libraries, data preparation (z-score, `log1p` metadata, reference intervals) |
| 3A | Stage 1 of model selection: does the inflammatory axis (CRP) separate states? Fitted on the CRP/WBC subcohort at K = 4 |
| 3B | Stage 2 of model selection: K ∈ {3, 4, 5} on the analysis cohort |
| 4–5 | Final model: K = 4, multivariate-normal emissions on HGB, CREA, ALT; 1 k-means start + 25 perturbed starts |
| 6–10 | Decoding, state-dependent distributions, transition matrix, pseudo-residuals, stayers/switchers |
| 11 | Uncertainty quantification |
| 12 | Export of the fitted-model bundle and artefact manifest |

Key configuration: `SEED = 10`, `TEST_FRAC = 0` (no train/test split: every model is fitted
on the full cohort of its stage), `APPLY_ZSCORE = TRUE`, `IS_a_MVN = TRUE`,
`OBS_COVARIATES = character(0)` (covariates are not HMM inputs),
`K_SELECTION_COHORT = FINAL_COHORT = "full"`, `SMOKE_TEST_N = 0`.

**Outputs:** `hmm_fit_objects_tesi4.rds` and the figures/tables in `figures/`.

```bash
Rscript -e 'rmarkdown::render("hmm_training_tesi.Rmd")'
```

### Stage 4 — `hmm_posthoc.Rmd`

**Input:** `hmm_fit_objects_tesi4.rds`.

Sections: stayer/switcher composition; depth-based admission phenotypes (Mahalanobis depth,
central subpopulations, deepest points); classification anchors (admission depth, day-7
online decoding on the panel truncated at day 7); in-hospital and 90-day mortality
(contingency tables, logistic regression); Kaplan–Meier survival; switcher trajectories and
the state-change alert rule.

Key configuration: `STAYER_WINDOW = 14`, `CENTRAL_FRAC = 0.50`, `DAY_TARGET = 7`,
`REF_STATE = "Stable"`, `SURVIVAL_FOLLOWUP = 90`, `P_ADJUST_METHOD = "BH"`,
`PERM_SEED = 41703192`.

**Outputs:** figures (PNG + PDF pairs) in `figures_ch4/`, CSV tables in `tables_ch4/`.

```bash
Rscript -e 'rmarkdown::render("hmm_posthoc.Rmd")'
```

---

## Reproducibility notes

**Random seeds.** Stage 1: `MICE_RANDOM_STATE = 42` (imputer and random forest). Stage 2:
`set.seed(42)`. Stage 3: `SEED = 10`; the r-th perturbed start uses `SEED + 7919·r`;
k-means initialisation uses `nstart = 25`. Stage 4: `PERM_SEED = 41703192`.

**Sanity checks against the reference run.** A correct execution reproduces:

| Check | Expected value | Where printed |
|---|---|---|
| Analysis cohort (run A) | 1039 patients, 15 585 rows, days 0–14 | Stage 3, Section 2.3 |
| CRP/WBC subcohort (run B) | 1007 patients, 15 105 rows | Stage 3, Section 2.3 |
| `log1p`-transformed axes (run A) | ALT, CREA | Stage 3, Section 2.3 |
| Best log-likelihood, final model (26 starts) | −45 052.92 | Stage 3, Section 5 |
| Stayers / switchers | 887 / 152 | Stage 4, Section 7 |
| Day-7 anchor | 1036 patients | Stage 4, Section 7 |

Minor numerical differences may arise on other platforms or BLAS/LAPACK builds, because the
HMM likelihood is multimodal and is maximised numerically.

**Runtime.** On the reference machine a single fit of the final specification takes
roughly 4–5 minutes; Stage 3 performs the model-selection fits plus 26 starts of the final
model, so it is by far the slowest stage. Set `SMOKE_TEST_N > 0` in Stage 3 to check that
the notebook runs end to end on a few patients — **every number produced under a smoke
test is meaningless; reset it to 0 for real runs.**

**Known pitfalls.**
- `MASS::select` masks `dplyr::select`; all dplyr verbs in Stage 4 are namespace-qualified.
- `hmmTMB`'s `viterbi()` takes no `newdata`; decoding of truncated panels is done by
  rebuilding an HMM with the trained parameters.
- `biomarker_reference_ranges.csv` must be in the Stage 3 working directory; if missing,
  the notebook only warns and silently drops the reference bands.

---

## Citation

If you use this code, please cite the thesis:

```bibtex
@mastersthesis{cecilia2026hmm,
  author = {Cecilia, Marco},
  title  = {[Thesis title]},
  school = {Politecnico di Milano},
  year   = {2026},
  type   = {Master's thesis in Mathematical Engineering}
}
```

Core software: `hmmTMB` (Michelot, 2025), `TMB` (Kristensen et al., 2016),
`mgcv` (Wood, 2017).

## License

_[Choose a license, e.g. MIT.]_
