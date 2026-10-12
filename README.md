# MRKCFC_Simulation

R code for the simulation studies of the manuscript

> Fan B, et al. Multivariate Robust K-Centers Functional Clustering Model for
> Sparsely Observed Non-Gaussian Multivariate Longitudinal Data. Under review.

The method itself is the R package
[MRKCFC](https://github.com/BingbingFan/MRKCFC). The results in the manuscript
were produced with MRKCFC 0.2.0, fdapace 0.6.0, fda 6.3.0 and R 4.5.1.

This repository holds the code that produces the raw results of the simulation
studies.

## Requirements

```r
install.packages(c("remotes", "fdapace", "fda", "mclust", "MASS", "lcmm", "teigen"))
remotes::install_github("BingbingFan/MRKCFC@v0.2.0")
```

Run every script from the root of this repository, because the scripts find each
other and the configuration through relative paths. Each script states its
arguments in its first lines.

## Layout

| Folder | Contents |
|---|---|
| `config/` | The design as tables: global parameters with the master seed (`params_global.csv`), the mean functions and eigenfunctions (`basis.csv`), the score dispersion matrix (`Z.csv`), the eight methods (`methods.csv`), the twelve scenarios of the method comparison (`scenarios.csv`, six score distributions at two separations), the thirteen one-factor configurations of the design study (`scenarios_b.csv`), and the two cells of the sensitivity studies (`scenarios_c.csv`). `validate_config.R` checks them. |
| `R/` | Data generation (`generate.R`), the eight methods (`methods.R`), the metrics (`metrics.R`), the robust separation measure (`separation.R`), the estimation measures (`estimation.R`), and one driver and one combiner per study. |

## How a study is run

Each driver divides its grid of cells into `n_tasks` slices and runs slice
`task_id`, both given on the command line. Every cell takes its seed from the
master seed in `config/params_global.csv` and its own identity, so a result does
not depend on how the grid is divided, and any slice or cell can be rerun alone.
For example, the method comparison in one process:

```sh
Rscript R/run_cell.R 1 1 results/studyA
Rscript R/combine.R results/studyA results/studyA/studyA
```

or in 250 slices, `Rscript R/run_cell.R <i> 250 results/studyA` for `i` = 1,
..., 250, in any order and on any number of cores, then the same combine.

Before a full study, check the configuration:

```sh
Rscript config/validate_config.R
```

A fit of MRKCFC-S at `n = 500` with `M = 4` took a median of 20 to 70 minutes
on one core depending on the scenario, and about 70 minutes when the loop ran
to its iteration cap; a fit of MKCFC took about one minute. The method
comparison alone is 9,600 fits, so the studies are best run in slices on many
cores.

## The simulation studies

All studies share the generator in `R/generate.R`, the methods in
`R/methods.R`, the metrics in `R/metrics.R` and the configuration tables.
`n = 500`, `M = 4` retained components, `K = 2` throughout, except where a
study varies them. Each combiner merges the per-task result files of its study
and writes the summary file named below.

| Study | What it varies | Driver | Combiner | Summary file | Manuscript |
|---|---|---|---|---|---|
| Method comparison (Study A; manuscript: method comparison) | Six score distributions at two separations, eight methods, 100 replicates | `Rscript R/run_cell.R <task_id> <n_tasks> <outdir>` | `Rscript R/combine.R <outdir> <outstem>` | `studyA_summary.csv` | Table 2, Tables S5 to S7 |
| Design study (Study B; manuscript: design perturbations) | One factor at a time from the contaminated-normal baseline: sample size, sparsity, separation, class proportions, number of variables; 50 replicates | `Rscript R/run_cell_b.R <task_id> <n_tasks> <outdir>` | `Rscript R/combine.R <outdir> <outstem>` | `studyB_summary.csv` | Table 4, Table S9 |
| Estimation on the true partition | The three second-stage constructions fitted under the true labels; principal angles and canonical correlations | `Rscript R/study_estim.R <task_id> <n_tasks> <outdir> [reps] [n]` | `Rscript R/combine_estim.R <outdir> <outstem>` | `estim_summary.csv` | Table 3, Table S8 |
| Component selection | The variance-explained rule on a global fit, at `n` = 500 and 1000 | `Rscript R/study_fve.R <task_id> <n_tasks> <outdir> [reps] [n] [obs_min]` | `Rscript R/combine_fve.R <outdir> <outstem> [sweep_summary.csv]` | `fve_summary.csv` | Table S2 |
| Component sweep | `M` from 1 to 6 on every scenario, 25 replicates | `Rscript R/study_sweep.R <task_id> <n_tasks> <outdir> [reps] [M_list] [n_list] [obs_min]` | `Rscript R/combine_sweep.R <outdir> <outstem> [fve.csv]` | `sweep_summary.csv` | Tables S3 and S4 |
| Initialization (C1; manuscript: sensitivity studies, initialization) | K-means, hierarchical, oracle and random starts, 10 replicates | `Rscript R/study_c1.R <scenario_ids> <methods> <reps> <n_restarts> [outdir]` | `Rscript R/combine_c1.R <outdir> <outstem>` | `c1_summary.csv` | Table S10 |
| Class-specific against class-common eigenfunctions (C2; manuscript: sensitivity studies, class-specific against class-common eigenfunctions) | The basis shared across classes or estimated per class, 50 replicates | `Rscript R/study_c2.R <task_id> <n_tasks> <outdir>` | `Rscript R/combine_c.R <outdir> <outstem>` | `c2_summary.csv` | Table S11 |
| Selecting K (C3; manuscript: sensitivity studies, selection of K) | Silhouette and Calinski-Harabasz over `K` = 2 to 6 on the multivariate scores, 50 replicates | `Rscript R/study_c3.R <task_id> <n_tasks> <outdir>` | `Rscript R/combine_c.R <outdir> <outstem>` | `c3_summary.csv` | Table S12 |
| Convergence thresholds (C4; manuscript: sensitivity studies, convergence threshold) | `conv.prop` of 0, 0.005, 0.01 and 0.02 on four scenarios, 25 replicates | `Rscript R/study_c4b.R <task_id> <n_tasks> <outdir>` | `Rscript R/combine_c4b.R <outdir> <outstem>` | `c4b_summary.csv` | Table S13 |
| Cost (C5; manuscript: Section 2.4, cost) | Median wall time as `n` varies | `Rscript R/study_c5.R <task_id> <n_tasks> <outdir>` | `Rscript R/combine_c.R <outdir> <outstem>` | `c5_summary.csv` | Section 2.4 |

The multivariate silhouette and the Calinski-Harabasz index of the K-selection
study, and the robust Mahalanobis separation of `separation.R`, are computed on
the multivariate scores of one global fit.

## License

GPL-3, as for the MRKCFC package.
