# Study C2: how sensitive is the method to how the basis is chosen?
#
#   Rscript R/study_c2.R <task_id> [n_tasks] [outdir] [n_reps_override]
#
# TWO AXES, crossed with the two C cells:
#
#   M in {1, 2, 3, 4}          components retained per variable (optns$nK)
#   basis in {specific, common} whether the eigenfunctions are re-estimated
#                               per cluster or shared across them
#
# M_true = 2, so the M axis brackets the truth: M = 1 is under-specified, M = 2
# correct, M = 3 and 4 over-specified. What C2 reports is how fast performance
# decays away from the truth, which is the question a reader has when choosing M
# on data where the truth is unknown.
#
# THE PVE AXIS WAS DROPPED, decided 2026-08-28, and it was not a capability gap.
# KFPCA() takes nK, a component COUNT, and has no variance-explained argument,
# but it does compute eigen(Kend) internally and merely discards the values, so
# returning them plus a share-to-count helper would have been about fifteen
# lines. It was dropped because with only two true components the four
# thresholds would very likely all map to the same M and produce a degenerate
# axis. If it is ever reinstated: a smoothed PACE-type covariance surface
# routinely has negative eigenvalues, which must be truncated before any
# cumulative share is taken.
#
# THE BASIS AXIS IS MISSPECIFIED BY CONSTRUCTION HERE, and that is the point.
# The design gives variables 1 and 2 sines for class 1 and cosines for class 2,
# so a common basis cannot represent both classes well. The DIRECTION of the
# result is therefore predictable; the SIZE is not, and the size is what turns
# the paper's design choice from an assertion into a measurement.

source(file.path("R", "study_c.R"))

args     <- commandArgs(trailingOnly = TRUE)
task_id  <- if (length(args) >= 1) as.integer(args[1]) else 1L
n_tasks  <- if (length(args) >= 2) as.integer(args[2]) else 1L
outdir   <- if (length(args) >= 3) args[3] else "results_c2"
rep_over <- if (length(args) >= 4) as.integer(args[4]) else NA_integer_
# Subject count override, 5th argument. For SMOKE TESTS only: it makes a cell
# cheap enough to exercise the driver without paying for a full n = 500 fit.
# Like the replicate override it is an argument rather than a config edit, so
# the reportable value in params_global.csv is never mutated to run a check.
n_over   <- if (length(args) >= 5) as.integer(args[5]) else NA_integer_

cfg     <- load_config()
if (!is.na(n_over)) {
  cfg$global$n <- n_over
  cat(sprintf("SMOKE TEST: n overridden to %d (reportable value is %s)
", n_over, load_config()$global$n))
}
configs <- load_c_configs()
n_reps  <- c_reps(cfg, rep_over)

# The MKCFC family only. The basis axis is a property of the K-centers loop:
# the score-then-cluster comparators fit one global MFPCA by construction and
# have no per-class basis to vary, so they would contribute a flat line.
use_methods <- c("M1", "M2", "M3")

# THE M AXIS WAS CUT ON 2026-09-04, leaving 200 cells rather than 800.
#
# This study is two studies. The BASIS ablation is the strongest positive result
# in Study C -- class-specific estimation buys about +0.887 aRI for MRKCFC-R2 --
# and it was only ever measured at `nK_fit = 2`, so it needs re-running at the
# current operating point. The M axis does NOT: stage 1 sweeps M = 1..6 over all
# twelve scenarios and all eight methods, against C2's two configs and three
# methods, and fully dominates it. Re-running M = 1..3 here would spend about
# 1,200 core-hours reproducing a weaker version of a result already in hand.
#
# `M_GRID` is the operating point from the config, NOT a literal 4, so this
# tracks `nK_fit` if the operating point moves again. Restore `seq_len(4)` only
# if C2 is deliberately being made an independent replication of stage 1 -- and
# then say so in the write-up, because the table's framing changes with it.
M_GRID <- as.integer(cfg$global$nK_fit)

grid <- expand.grid(rep = seq_len(n_reps), M = M_GRID,
                    basis_common = c(FALSE, TRUE),
                    config_id = configs$config_id, stringsAsFactors = FALSE)
grid <- grid[order(grid$config_id, grid$basis_common, grid$M, grid$rep), ]
rownames(grid) <- NULL
mine <- c_slice(grid, task_id, n_tasks)
if (!nrow(mine)) { cat("nothing to do\n"); quit(save = "no", status = 0) }
c_header("C2 components and basis", nrow(mine), nrow(grid), task_id, n_tasks)

started <- Sys.time()
out <- vector("list", nrow(mine))
for (a in seq_len(nrow(mine))) {
  cell <- mine[a, ]
  row  <- configs[configs$config_id == cell$config_id, ]
  t0 <- Sys.time()
  res <- tryCatch({
    c2 <- apply_config_c(cfg, row)
    c2$global$nK_fit <- cell$M          # the component count, per variable
    dat <- generate_replicate(c2, as_scenario_c(row), rep = cell$rep)
    rows <- do.call(rbind, lapply(use_methods, function(mid) {
      spec <- cfg$methods[cfg$methods$method_id == mid, ]
      # Common random numbers across methods, M and basis within this replicate.
      c_seed_fit(c2, cell$config_id, cell$rep)
      tt0 <- Sys.time()
      fit <- tryCatch(MRKCFC::MKCFC(
        Ly = dat$Ly, Lt = dat$Lt, K = as.integer(c2$global$K_true),
        robust = as.logical(spec$robust),
        optns = build_optns(c2, as.logical(spec$robust)),
        mfpca_scatter = spec$mfpca_scatter,
        init.method = c2$global$init_method,
        maxIter = as.integer(c2$global$maxIter),
        conv.prop = as.numeric(c2$global$conv_prop),
        LOO = isTRUE(as.logical(c2$global$LOO)),
        basis_common = cell$basis_common,
        parallel = FALSE, trace = 0), error = function(e) e)
      el <- as.numeric(difftime(Sys.time(), tt0, units = "secs"))
      if (inherits(fit, "error"))
        return(data.frame(method = mid, CCR = NA_real_, aRI = NA_real_,
                          converged = NA, iter = NA_integer_, seconds = el,
                          error = conditionMessage(fit), stringsAsFactors = FALSE))
      g <- as.integer(fit$res_clust$g)
      data.frame(method = mid, CCR = ccr(dat$labels, g), aRI = ari(dat$labels, g),
                 converged = isTRUE(fit$conv), iter = fit$iter, seconds = el,
                 error = NA_character_, stringsAsFactors = FALSE)
    }))
    rows
  }, error = function(e) data.frame(
    method = NA_character_, CCR = NA_real_, aRI = NA_real_, converged = NA,
    iter = NA_integer_, seconds = NA_real_,
    error = paste("cell failed:", conditionMessage(e)), stringsAsFactors = FALSE))

  res$config_id <- cell$config_id; res$role <- row$role
  res$M <- cell$M; res$basis <- if (cell$basis_common) "common" else "specific"
  res$rep <- cell$rep
  out[[a]] <- res
  cat(sprintf("  [%3d/%3d] %-7s M=%d %-8s rep %3d  %6.1fs\n", a, nrow(mine),
              cell$config_id, cell$M, if (cell$basis_common) "common" else "specific",
              cell$rep, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  utils::flush.console()
}

results <- do.call(rbind, out)
results$method_id <- results$method
results$method <- cfg$methods$label[match(results$method_id, cfg$methods$method_id)]
c_write(results, outdir, "studyC2", task_id, started)
