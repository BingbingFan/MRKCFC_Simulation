# Estimation accuracy on the TRUE partition, across the twelve Study A scenarios.
#
#   Rscript R/study_estim.R <task> [n_tasks] [outdir] [reps] [n]
#
# The measurement itself lives in R/estimation.R as oracle_estimation(), which
# the sweep will also call once the two are merged. This file is the driver:
# it walks the grid, writes one file, and summarises NOTHING.
#
# WHY IT NO LONGER SUMMARISES. The previous version printed three tables from
# `results`, which is only THIS TASK'S SLICE -- correct at n_tasks = 1 and
# silently wrong at any array size. It is run as 100 tasks, so
# every one of those hundred tasks printed a table that looked complete and was
# 99 percent empty. That is the same failure as the `tail -25` that hid the
# C_HARD rows for three days. Read the study with R/combine_estim.R.
#
# WHY THE RESULTS FROM BEFORE 2026-09-02 ARE VOID. The measurement used
# `M <- cfg$global$M_true` to index the FITTED rotation. That was harmless while
# nK_fit == M_true == 2 and is silently wrong at nK_fit = 4, where the rotation
# is 12 x 12 with blocks of size 4 and the old code read rows 1:2, 3:4, 5:6. It
# does not error; it returns meaningless numbers. See R/estimation.R for the fix.
#
# `results_estim_M2/` remains valid AS A nK_fit = 2 RESULT and is not comparable
# to anything produced after the config change of 2026-09-02.

source(file.path("R", "study_c.R"))        # load_config, build_optns, c_seed_fit, c_*
source(file.path("R", "estimation.R"))     # oracle_estimation, and the argument for it

args    <- commandArgs(trailingOnly = TRUE)
task_id <- if (length(args) >= 1) as.integer(args[1]) else 1L
n_tasks <- if (length(args) >= 2) as.integer(args[2]) else 1L
outdir  <- if (length(args) >= 3) args[3] else "results_estim"
n_reps  <- if (length(args) >= 4) as.integer(args[4]) else 100L
n_over  <- if (length(args) >= 5) as.integer(args[5]) else NA_integer_

cfg <- load_config()
if (!is.na(n_over)) cfg$global$n <- n_over

grid <- expand.grid(rep = seq_len(n_reps), scenario_id = cfg$scenarios$scenario_id,
                    stringsAsFactors = FALSE)
grid <- grid[order(grid$scenario_id, grid$rep), ]
rownames(grid) <- NULL
mine <- c_slice(grid, task_id, n_tasks)
if (!nrow(mine)) { cat("nothing to do\n"); quit(save = "no", status = 0) }
c_header("estimation accuracy, function space, true partition",
         nrow(mine), nrow(grid), task_id, n_tasks)
cat(sprintf("  nK_fit = %s (the FIT) | M_true = %s (the TRUTH) | n = %s | obs_min = %s\n",
            cfg$global$nK_fit, cfg$global$M_true, cfg$global$n, cfg$global$obs_min))

started <- Sys.time()
out <- list()
for (a in seq_len(nrow(mine))) {
  cell <- mine[a, ]
  sc <- cfg$scenarios[cfg$scenarios$scenario_id == cell$scenario_id, ]
  t0 <- Sys.time()

  rows <- tryCatch({
    dat <- generate_replicate(cfg, sc, rep = cell$rep)
    oracle_estimation(dat, cfg, d = sc$d)
  }, error = function(e) data.frame(
    variant = NA_character_, class = NA_integer_, angle = NA_real_,
    angle_lead3 = NA_real_, angle_lead4 = NA_real_,
    score_cc = NA_real_, nK_used = NA_integer_,
    error = paste("cell failed:", conditionMessage(e)), stringsAsFactors = FALSE))

  rows$scenario_id <- cell$scenario_id
  rows$arm         <- sc$arm
  rows$arm_label   <- sc$arm_label
  rows$d           <- sc$d
  rows$rep         <- cell$rep
  rows$nK_fit      <- as.integer(cfg$global$nK_fit)
  rows$M_true      <- as.integer(cfg$global$M_true)
  rows$n           <- as.integer(cfg$global$n)
  rows$obs_min     <- as.integer(cfg$global$obs_min)
  out[[length(out) + 1L]] <- rows
  cat(sprintf("  [%3d/%3d] %-4s rep %3d  %6.1fs%s\n", a, nrow(mine),
              cell$scenario_id, cell$rep,
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              if (any(!is.na(rows$error))) sprintf("  (%d err)", sum(!is.na(rows$error))) else ""))
  utils::flush.console()
}

results <- do.call(rbind, out)

# A completeness count, not a summary: every figure below is a property of THIS
# task and needs no sight of any other.
cat(sprintf("\ncells %d | rows %d | errors %d | NA angle %d | %.1f min\n",
            nrow(mine), nrow(results), sum(!is.na(results$error)),
            sum(is.na(results$angle)),
            as.numeric(difftime(Sys.time(), started, units = "mins"))))
c_write(results, outdir, "estim", task_id, started)
cat("Read the study with R/combine_estim.R -- this file does not summarise.\n")
