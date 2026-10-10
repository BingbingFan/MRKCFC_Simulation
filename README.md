# MRKCFC_Simulation

R code for the simulation study and the real-data analysis of the manuscript

> Fan B, et al. Multivariate Robust K-Centers Functional Clustering Model for
> Sparsely Observed Non-Gaussian Multivariate Longitudinal Data. Under review.

The method itself is the R package
[MRKCFC](https://github.com/BingbingFan/MRKCFC). The results in the manuscript
were produced with MRKCFC 0.2.0, fdapace 0.6.0, fda 6.3.0 and R 4.5.1.

## Requirements

```r
install.packages(c("remotes", "fdapace", "fda", "mclust", "MASS", "lcmm", "teigen",
                   "survival", "mgcv", "ggplot2", "patchwork"))
remotes::install_github("BingbingFan/MRKCFC@v0.2.0")
# optional, for the time-dependent AUC in the survival benchmark
install.packages("timeROC")
```

Run every script from the root of this repository, because the scripts find each
other and the configuration through relative paths. Each script states its
arguments in its first lines.

## Layout

| Folder | Contents |
|---|---|
| `config/` | The design as tables: global parameters with the master seed (`params_global.csv`), the mean functions and eigenfunctions (`basis.csv`), the score dispersion matrix (`Z.csv`), the eight methods (`methods.csv`), the twelve scenarios of the method comparison (`scenarios.csv`, six score distributions at two separations), the thirteen one-factor configurations of the design study (`scenarios_b.csv`), and the two cells of the sensitivity studies (`scenarios_c.csv`). `validate_config.R` checks them. |
| `R/` | Data generation (`generate.R`), the eight methods (`methods.R`), the metrics (`metrics.R`), the robust separation measure (`separation.R`), one driver and one combiner per study, the table generator (`make_tables.R`), and the real-data scripts (`realdata_*.R`, `make_rd_tables.R`, `plot_fig1_traj.R`). |

## How a study is run

Each driver divides its grid of cells into `n_tasks` slices and runs slice
`task_id`. Every cell takes its seed from the master seed in
`config/params_global.csv` and its own identity, so a result does not depend on
how the grid is divided, and any slice or cell can be rerun alone. In an array
job the two numbers come from the scheduler; on one machine they are given by
hand. For example, the method comparison in one process:

```sh
Rscript R/run_cell.R 1 1 results/studyA
Rscript R/combine.R results/studyA results/studyA/studyA
```

or in 250 slices, `Rscript R/run_cell.R <i> 250 results/studyA` for `i` = 1,
..., 250, in any order and on any number of cores, then the same combine.

The checks to run before a full study:

```sh
Rscript config/validate_config.R
Rscript R/test_generate.R
Rscript R/test_metrics.R
Rscript R/test_separation.R
Rscript R/test_methods.R
Rscript R/test_study_b.R
Rscript R/schema_audit.R
```

When the studies have run and been combined, the simulation tables of the
manuscript are written by

```sh
Rscript R/make_tables.R results tables
```

where `results` holds the combined summary files named in the table below.

A fit of MRKCFC-S at `n = 500` with `M = 4` took a median of 20 to 70 minutes
on one core depending on the scenario, and about 70 minutes when the loop ran
to its iteration cap; a fit of MKCFC took about one minute. The method
comparison alone is 9,600 fits, so the studies are best run in slices on many
cores.

## The simulation studies

All studies share the generator in `R/generate.R`, the methods in
`R/methods.R`, the metrics in `R/metrics.R` and the configuration tables.
`n = 500`, `M = 4` retained components, `K = 2` throughout, except where a
study varies them.

| Study | What it varies | Driver | Combiner | Summary file | Manuscript |
|---|---|---|---|---|---|
| Method comparison (Study A) | Six score distributions at two separations, eight methods, 100 replicates | `Rscript R/run_cell.R <task_id> <n_tasks> <outdir>` | `Rscript R/combine.R <outdir> <outstem>` | `studyA_summary.csv` | Table 2, Tables S4 to S6 |
| Design study (Study B) | One factor at a time from the contaminated-normal baseline: sample size, sparsity, separation, class proportions, number of variables; 50 replicates | `Rscript R/run_cell_b.R <task_id> <n_tasks> <outdir>` | `Rscript R/combine.R <outdir> <outstem>` | `studyB_summary.csv` | Table 4, Table S8 |
| Estimation on the true partition | The three second-stage constructions fitted under the true labels; principal angles and canonical correlations | `Rscript R/study_estim.R <task_id> <n_tasks> <outdir> [reps] [n]` | `Rscript R/combine_estim.R <outdir> <outstem>` | `estim_summary.csv` | Table 3, Table S7 |
| Component selection | The variance-explained rule on a global fit, at `n` = 500 and 1000 | `Rscript R/study_fve.R <task_id> <n_tasks> <outdir> [reps] [n] [obs_min]` | `Rscript R/combine_fve.R <outdir> <outstem> [sweep_summary.csv]` | `fve_summary.csv` | Table S1 |
| Component sweep | `M` from 1 to 6 on every scenario, 25 replicates | `Rscript R/study_sweep.R <task_id> <n_tasks> <outdir> [reps] [M_list] [n_list] [obs_min]` | `Rscript R/combine_sweep.R <outdir> <outstem> [fve.csv]` | `sweep_summary.csv` | Tables S2 and S3 |
| Initialization (C1) | K-means, hierarchical, oracle and random starts, 10 replicates | `Rscript R/study_c1.R <scenario_ids> <methods> <reps> <n_restarts> [outdir]` | `Rscript R/combine_c1.R <outdir> <outstem>` | `c1_summary.csv` | Table S9 |
| Class-specific against class-common eigenfunctions (C2) | The basis shared across classes or estimated per class, 50 replicates | `Rscript R/study_c2.R <task_id> <n_tasks> <outdir>` | `Rscript R/combine_c.R <outdir> <outstem>` | `c2_summary.csv` | Table S10 |
| Selecting K (C3) | Silhouette and Calinski-Harabasz over `K` = 2 to 6 on the multivariate scores, 50 replicates | `Rscript R/study_c3.R <task_id> <n_tasks> <outdir>` | `Rscript R/combine_c.R <outdir> <outstem>` | `c3_summary.csv` | Table S11 |
| Convergence thresholds (C4) | `conv.prop` of 0, 0.005, 0.01 and 0.02 on four scenarios, 25 replicates | `Rscript R/study_c4b.R <task_id> <n_tasks> <outdir>` | `Rscript R/combine_c4b.R <outdir> <outstem>` | `c4b_summary.csv` | Table S12 |
| Cost (C5) | Median wall time as `n` varies | `Rscript R/study_c5.R <task_id> <n_tasks> <outdir>` | `Rscript R/combine_c.R <outdir> <outstem>` | `c5_summary.csv` | Section 2.4 |

The multivariate silhouette and the Calinski-Harabasz index of the K-selection
study, and the robust Mahalanobis separation of `separation.R` that the design
tables quote, are computed on the multivariate scores of one global fit.

`R/make_tables.R` writes one file per table, named by content.

| File | Table | File | Table |
|---|---|---|---|
| `tab_studyA_ari.tex` | 2, adjusted Rand index by scenario | `tab_studyB_convergence.tex` | S8, convergence under the perturbations |
| `tab_estim_angle.tex` | 3, principal angles | `tab_c1_init.tex` | S9, starting partition |
| `tab_studyB_ari.tex` | 4, adjusted Rand index under perturbations | `tab_c2_basis.tex` | S10, class-specific against class-common |
| `tab_fve_selection.tex` | S1, components selected | `tab_c3_kselect.tex` | S11, recovery of K |
| `tab_sweep_M.tex` | S2, adjusted Rand index by M | `tab_c4b_agree.tex` | S12, convergence thresholds |
| `tab_sweep_convergence.tex` | S3, convergence by M | | |
| `tab_studyA_ccr.tex` | S4, correct classification rate | | |
| `tab_studyA_rre.tex` | S5, relative reconstruction error | | |
| `tab_studyA_convergence.tex` | S6, convergence by scenario | | |
| `tab_estim_cc.tex` | S7, canonical correlation | | |

Table 1 of the manuscript (the six score distributions) is typed from
`config/scenarios.csv`.

## The real-data analysis

The cohort data cannot be shared publicly because of patient privacy, and they
are available from the corresponding author on reasonable request. The loader
`R/realdata_load.R` expects the file `data/data_yn_use_digit1_CEA_199_242.rds`
and documents its contents in its first lines; it returns the 36-month analysis
set with an integer patient identifier. No result file of the real-data
analysis is included here, because those files hold a class label for every
patient.

The analysis runs on the original scale with the time axis transformed through
the empirical distribution function of the visit times and the Kendall-surface
bandwidth chosen by cross-validation (the primary arm), and again on the log
scale at `M` = 3 and `M` = 4, and with the bandwidth fixed at 0.20 (the
sensitivity arms of Table S21). Each arm writes to its own directory. The
steps, in order:

| Step | Scripts | Output |
|---|---|---|
| 1. Number of components | `Rscript R/realdata_fve.R <outstem> [scales] [time]` | The variance-explained rule on the cohort, `M` = 4 on the original scale and 3 on the log scale (Table S16) |
| 2. Fits | `Rscript R/realdata_fit_task.R <task_id> [outdir] [nK] [log] [n_restart] [bwK] [kmax] [time]`, one task per (method, `K`, restart): 3 methods, `K` from 2 to 4, 10 random restarts, 90 tasks per arm; a task resumes from its own checkpoint if interrupted | One checkpoint per task with the best iterate and the objective trace |
| 3. Selection | `Rscript R/realdata_combine.R <outdir> <outstem> 10 4` | The restart with the smallest objective per cell, `K` by the silhouette on the multivariate scores of a global fit, the labels (Table S17) |
| 4. Stability | `Rscript R/realdata_subsample_task.R <draw> <results-dir> <combined-stem> <outdir>` for 100 draws of 80 per cent, then `Rscript R/realdata_stability.R <results-dir> <combined-stem>` | Agreement among restarts, clusterwise Jaccard similarity over restarts and subsamples (Table S18) |
| 5. Agreement between methods | `Rscript R/realdata_agreement.R <results-dir> <combined-stem>` | Adjusted Rand index between the partitions, including the comparators fitted on the global scores (Table S19) |
| 6. Survival | `Rscript R/realdata_survival.R <results-dir> <combined-stem> OS 36 60 <covariates>` and the same with `RFS` | Cox models by class with and without the covariates (Tables 5 and S20) |
| 7. Diagnostics | `Rscript R/realdata_bwprobe.R`, `Rscript R/realdata_classcheck.R <outstem> <dir>`, `Rscript R/realdata_recon_trim.R <rd-dir> <out-dir>` | Where the robust smoother can be fitted on this cohort, which classes the loop formed, and the share of the reconstruction error carried by the most extreme observations (Section 3.2) |
| 8. Tables and figures | `Rscript R/make_rd_tables.R <rd-dir> <tables-dir> <figures-dir>` and `Rscript R/plot_fig1_traj.R <rd-dir> <out-stem>` | Table 5, Tables S14 to S21, Figures 1, 2 and S1 |

## License

GPL-3, as for the MRKCFC package.
