# Study C4: convergence. Fixes conv.prop, and reports the non-convergence rate.
#
#   Rscript R/study_c4.R <task_id> [n_tasks] [outdir] [n_reps_override]
#
# conv.prop in {0, 0.005, 0.01, 0.02} at maxIter = 100, across ALL TWELVE Study
# A scenarios, M1/M2/M3, 50 replicates. Decided 2026-08-28. All twelve rather
# than the four t arms so that non-convergence being confined to the heavy tails
# is VERIFIED rather than assumed.
#
# READ `iter`, NOT JUST `converged`. converged = FALSE has two causes needing
# opposite remedies, and MKCFC() only updates `iteration` AFTER the <10-curves
# break, so they are distinguishable:
#
#   iter == maxIter  the loop exhausted the cap
#   iter <  maxIter  it broke early because a group fell below ten curves,
#                    which NO convergence threshold can fix, because the break
#                    happens before the threshold is ever tested
#
# The 10-replicate round did not record iter, which is why its 22.5 percent
# non-convergence figure was never actionable.
#
# -----------------------------------------------------------------------------
# WARNING, READ BEFORE SUBMITTING. maxIter = 100 is expensive and the evidence
# for it has weakened since it was chosen.
#
# Budget is roughly 1,000 core-hours, about 2.5x what maxIter = 20 costs and
# more than twice Study A. Meanwhile c4b.rds (2026-08-28, replicates 4-10 of
# t nu=1.5 d=3, M2 alone) found:
#
#   conv.prop   converged   median aRI   cap-bound at iter=20
#   0           4/7         +0.002       1
#   0.005       5/7         +0.005       0
#   0.02        5/7         +0.001       0
#
# Only ONE replicate in seven was cap-bound at 20; replicates 5 and 9 failed at
# iter 8 and 9, identically at every threshold, which is the <10-curves break
# that maxIter cannot touch. And relaxing the threshold bought convergence
# without buying accuracy: median aRI sits at about zero throughout, and the
# one run rescued into convergence landed at aRI 0.005. Results are bimodal,
# roughly 0.7-0.8 or roughly 0, and conv.prop never moves a run between those
# basins.
#
# The implication is that INITIALISATION, not convergence, is what drives arm 6
# -- consistent with C1's finding that the two-staged fallback draws a SINGLE
# random start.
#
# RECONSIDERED AND CONFIRMED 2026-08-28, with the evidence above in view:
# maxIter STAYS AT 100. The argument for keeping it is that it removes the cap
# as a confound entirely, so any non-convergence C4 reports is unambiguously
# attributable to the threshold rather than to the iteration budget -- which is
# what makes the reported rate defensible. The cost is accepted.
# `iter` still records whether the cap bound anything, so the question stays
# answerable from the output either way. Setting MAXITER <- 20L remains a
# one-line change if the budget is later needed elsewhere.
# -----------------------------------------------------------------------------

source(file.path("R", "study_c.R"))

args     <- commandArgs(trailingOnly = TRUE)
task_id  <- if (length(args) >= 1) as.integer(args[1]) else 1L
n_tasks  <- if (length(args) >= 2) as.integer(args[2]) else 1L
outdir   <- if (length(args) >= 3) args[3] else "results_c4"
rep_over <- if (length(args) >= 4) as.integer(args[4]) else NA_integer_
# Subject count override, 5th argument. For SMOKE TESTS only: it makes a cell
# cheap enough to exercise the driver without paying for a full n = 500 fit.
# Like the replicate override it is an argument rather than a config edit, so
# the reportable value in params_global.csv is never mutated to run a check.
n_over   <- if (length(args) >= 5) as.integer(args[5]) else NA_integer_

cfg    <- load_config()
if (!is.na(n_over)) {
  cfg$global$n <- n_over
  cat(sprintf("SMOKE TEST: n overridden to %d (reportable value is %s)
", n_over, load_config()$global$n))
}
n_reps <- c_reps(cfg, rep_over)
THRESHOLDS <- c(0, 0.005, 0.01, 0.02)
MAXITER    <- 100L
use_methods <- c("M1", "M2", "M3")   # conv.prop is a property of the K-centers loop

grid <- expand.grid(rep = seq_len(n_reps), conv_prop = THRESHOLDS,
                    scenario_id = cfg$scenarios$scenario_id,
                    stringsAsFactors = FALSE)
grid <- grid[order(grid$scenario_id, grid$conv_prop, grid$rep), ]
rownames(grid) <- NULL
mine <- c_slice(grid, task_id, n_tasks)
if (!nrow(mine)) { cat("nothing to do\n"); quit(save = "no", status = 0) }
c_header(sprintf("C4 convergence (maxIter = %d)", MAXITER),
         nrow(mine), nrow(grid), task_id, n_tasks)

started <- Sys.time()
out <- vector("list", nrow(mine))
for (a in seq_len(nrow(mine))) {
  cell <- mine[a, ]
  sc <- cfg$scenarios[cfg$scenarios$scenario_id == cell$scenario_id, ]
  t0 <- Sys.time()
  res <- tryCatch({
    c2 <- cfg
    c2$global$maxIter   <- MAXITER
    c2$global$conv_prop <- cell$conv_prop
    dat <- generate_replicate(c2, sc, rep = cell$rep)
    rows <- run_all_methods(dat, c2, methods = use_methods)
    rows[, c("method", "CCR", "aRI", "converged", "iter", "seconds", "error")]
  }, error = function(e) data.frame(
    method = NA_character_, CCR = NA_real_, aRI = NA_real_, converged = NA,
    iter = NA_integer_, seconds = NA_real_,
    error = paste("cell failed:", conditionMessage(e)), stringsAsFactors = FALSE))

  res$scenario_id <- cell$scenario_id
  res$conv_prop <- cell$conv_prop
  res$maxIter <- MAXITER
  res$rep <- cell$rep
  # The diagnosis, computed once here so no downstream reader has to remember
  # the rule. NA when the fit converged or failed outright.
  res$failure <- ifelse(is.na(res$converged) | res$converged, NA_character_,
                 ifelse(res$iter >= MAXITER, "cap", "few_curves"))
  out[[a]] <- res
  cat(sprintf("  [%3d/%3d] %-4s cp=%-6s rep %3d  %6.1fs\n", a, nrow(mine),
              cell$scenario_id, format(cell$conv_prop), cell$rep,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  utils::flush.console()
}

results <- do.call(rbind, out)
results$method_id <- results$method
results$method <- cfg$methods$label[match(results$method_id, cfg$methods$method_id)]

cat("\n=== non-convergence, and WHY, by threshold ===\n")
cat("  cap        = exhausted maxIter; a higher cap might help\n")
cat("  few_curves = broke on the <10-curves rule; no threshold can help\n\n")
for (cp in THRESHOLDS) {
  s <- results[results$conv_prop == cp & !is.na(results$converged), ]
  if (!nrow(s)) next
  cat(sprintf("  conv.prop=%-6s converged %3d/%3d  median aRI %+.3f  cap %d  few_curves %d\n",
              format(cp), sum(s$converged), nrow(s),
              stats::median(s$aRI, na.rm = TRUE),
              sum(s$failure == "cap", na.rm = TRUE),
              sum(s$failure == "few_curves", na.rm = TRUE)))
}

c_write(results, outdir, "studyC4", task_id, started)
