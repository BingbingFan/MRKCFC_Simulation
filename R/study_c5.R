# Study C5: cost. Median wall clock as n varies.
#
#   Rscript R/study_c5.R <task_id> [n_tasks] [outdir] [n_reps_override]
#
# SCOPE, decided 2026-08-28: n ONLY. The plan asks for n, P and K, but P in
# {5,10} and K in {3,4,5} need config work that belongs to Study B -- and K also
# needs the centroid geometry settled and runs on code paths never executed
# above K = 2. Making a Tier 1 sub-study wait on Tier 2 work risks it being
# dropped when Tier 2 is cut, so C5 measures the axis that runs today.
#
# THE OTHER AXES ARE NOT LOST. Study A's `seconds` column already gives median
# cost per method per scenario across all twelve cells at no extra compute,
# covering the arm and separation axes; Study B's results cover P once it runs.
# C5 supplies the n axis, which neither of them varies. Say plainly in the paper
# which axes were measured and where each came from.
#
# WHAT TO REPORT. Median rather than mean wall clock: the distribution is
# strongly right-skewed because a non-converging fit runs to the iteration cap
# while a converging one stops at 6 to 13, so a mean is dominated by the
# failures and describes no typical run. Report the median with the
# non-convergence rate beside it, since the two together are what a reader needs
# to predict their own cost.

source(file.path("R", "study_c.R"))

args     <- commandArgs(trailingOnly = TRUE)
task_id  <- if (length(args) >= 1) as.integer(args[1]) else 1L
n_tasks  <- if (length(args) >= 2) as.integer(args[2]) else 1L
outdir   <- if (length(args) >= 3) args[3] else "results_c5"
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
# Fewer replicates than the rest of C: a median wall clock is far less variable
# than a convergence rate, and these are the most expensive cells in the project
# at n = 2000. Overridable with the fourth argument.
n_reps  <- c_reps(cfg, if (is.na(rep_over)) 10L else rep_over)
N_GRID  <- c(250L, 500L, 1000L, 2000L)
# C5 sets n per CELL from N_GRID, so the n override cannot work the way it does
# in the other drivers -- cfg$global$n is overwritten below. Collapse the grid
# instead, or a "smoke test" would run the n = 2000 cells, which are the most
# expensive in the project.
if (!is.na(n_over)) N_GRID <- as.integer(n_over)
use_methods <- c("M1", "M2", "M3", "M8")

grid <- expand.grid(rep = seq_len(n_reps), n = N_GRID,
                    config_id = configs$config_id, stringsAsFactors = FALSE)
grid <- grid[order(grid$config_id, grid$n, grid$rep), ]
rownames(grid) <- NULL
mine <- c_slice(grid, task_id, n_tasks)
if (!nrow(mine)) { cat("nothing to do\n"); quit(save = "no", status = 0) }
c_header("C5 cost against n", nrow(mine), nrow(grid), task_id, n_tasks)

started <- Sys.time()
out <- vector("list", nrow(mine))
for (a in seq_len(nrow(mine))) {
  cell <- mine[a, ]
  row  <- configs[configs$config_id == cell$config_id, ]
  t0 <- Sys.time()
  res <- tryCatch({
    c2 <- apply_config_c(cfg, row)
    c2$global$n <- cell$n
    dat <- generate_replicate(c2, as_scenario_c(row), rep = cell$rep)
    rows <- run_all_methods(dat, c2, methods = use_methods)
    rows[, c("method", "CCR", "aRI", "converged", "iter", "seconds", "error")]
  }, error = function(e) data.frame(
    method = NA_character_, CCR = NA_real_, aRI = NA_real_, converged = NA,
    iter = NA_integer_, seconds = NA_real_,
    error = paste("cell failed:", conditionMessage(e)), stringsAsFactors = FALSE))

  res$config_id <- cell$config_id; res$role <- row$role
  res$n <- cell$n; res$rep <- cell$rep
  out[[a]] <- res
  cat(sprintf("  [%3d/%3d] %-7s n=%-5d rep %3d  %6.1fs\n", a, nrow(mine),
              cell$config_id, cell$n, cell$rep,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  utils::flush.console()
}

results <- do.call(rbind, out)
results$method_id <- results$method
results$method <- cfg$methods$label[match(results$method_id, cfg$methods$method_id)]

cat("\n=== median seconds per fit, by method and n ===\n")
tab <- stats::aggregate(seconds ~ method + n, results, stats::median,
                        na.action = stats::na.omit)
print(stats::reshape(tab, idvar = "method", timevar = "n", direction = "wide"),
      row.names = FALSE)
cat("\nRead beside the non-convergence rate: a failing fit runs to the cap, so\n")
cat("cost and convergence are not independent.\n")

c_write(results, outdir, "studyC5", task_id, started)
