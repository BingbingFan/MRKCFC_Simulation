# C4 at M = 4: is the non-convergence COSMETIC or SUBSTANTIVE?
#
#   Rscript R/study_c4b.R <task> [n_tasks] [outdir] [reps] [scenarios] [conv_props]
#
# Defaults: reps 25, scenarios "A02,A07,A11,A12", conv_props "0,0.005,0.01,0.02".
#
# WHY THIS IS NOT A RE-RUN OF study_c4.R
#
# C4 at nK_fit = 2 concluded that conv.prop "buys convergence, not accuracy":
# relaxing 0 -> 0.02 moved convergence 0.948 -> 0.960 and median aRI 0.559 ->
# 0.535, i.e. slightly WORSE. That conclusion was correct and is now irrelevant,
# because the failure MODE has changed.
#
# `converged = FALSE` has two causes needing opposite remedies:
#
#   iter == maxIter  the loop exhausted the cap        -> a threshold can help
#   iter <  maxIter  a group fell below ten curves     -> NO threshold can help,
#                    because the break happens before the threshold is tested
#
# At nK_fit = 2, 293 of 350 non-converged fits stopped at iter <= 10. Those are
# early breaks, so conv.prop was aimed at a mode that barely existed.
#
# At nK_fit = 4 the median iter for MRKCFC-R2 runs 47 to 97 against a cap of 100
# -- A11 sits at 97 with convergence 0.16, A07 at 87 with 0.56 -- where at
# nK_fit = 3 the same fits settle in about ten. The failures are now CAP-BOUND,
# which is exactly the mode conv.prop addresses and exactly the mode the earlier
# study did not have.
#
# THE STATISTIC THAT DECIDES IT: `agree_strict`
#
# Reporting aRI against the truth at each threshold is NOT enough. Two different
# partitions can score the same aRI against the truth, so matching columns would
# not show that the threshold left the answer alone.
#
# So every fit is also scored against the conv.prop = 0 partition FROM THE SAME
# REPLICATE, in `agree_strict`. That is the whole question:
#
#   agree_strict ~ 1.0   the partition stabilised long ago and a handful of
#                        boundary subjects oscillate forever. The
#                        non-convergence is COSMETIC and can be reported as such.
#   agree_strict << 1    the partition genuinely has not settled, and M = 4 has
#                        a real problem that no threshold disguises.
#
# THE TRAP THIS MUST NOT FALL INTO. Do not relax conv.prop until the method
# converges and then report that rate -- it would be loosening the definition
# of success until it was met. The
# defensible structure is conv.prop = 0 as the PRIMARY result, with the relaxed
# thresholds as a sensitivity analysis showing the partitions are unchanged.
# That makes this evidence about oscillation, not a redefinition.
#
# WHY FOUR SCENARIOS. At nK_fit = 4 the MKCFC family costs about 7,400 s per
# replicate against 1,000 s at nK_fit = 2, so the original twelve-scenario
# design is now 4,934 core-hours -- more than Study A. These four span the full
# range of M = 4 convergence and include both infinite-variance arms:
#
#   A02  Gaussian d = 3           conv 0.96, iter 17   healthy control
#   A07  contaminated d = 1.5     conv 0.56, iter 87
#   A12  t(1.5) d = 3             conv 0.36, iter 61
#   A11  t(1.5) d = 1.5           conv 0.16, iter 97   the flagship failure
#
# A02 is not padding: if the threshold moves partitions where the loop ALREADY
# converges, the instrument is suspect and nothing else here can be trusted.
#
# SLICING is on (scenario, replicate, method), so all four thresholds for one
# fit land in ONE task and `agree_strict` never crosses a task boundary.
#
# BUDGET about 822 core-hours: 4 x 25 x 3 = 300 cells, 4 fits each. M1 cells run
# roughly 5.4 h, M2 about 2.8 h, M3 minutes.

source(file.path("R", "study_c.R"))

args     <- commandArgs(trailingOnly = TRUE)
task_id  <- if (length(args) >= 1) as.integer(args[1]) else 1L
n_tasks  <- if (length(args) >= 2) as.integer(args[2]) else 1L
outdir   <- if (length(args) >= 3) args[3] else "results_c4b"
n_reps   <- if (length(args) >= 4) as.integer(args[4]) else 25L
SCEN     <- if (length(args) >= 5) strsplit(args[5], ",")[[1]] else
              c("A02", "A07", "A11", "A12")
CP       <- if (length(args) >= 6) as.numeric(strsplit(args[6], ",")[[1]]) else
              c(0, 0.005, 0.01, 0.02)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cfg <- load_config()
K   <- as.integer(cfg$global$K_true)
stopifnot(CP[1] == 0)   # the reference partition must be the strict fit

# The MKCFC family only. The score-then-cluster comparators have no K-centers
# loop, so conv.prop does not exist for them.
VAR <- list(M1 = list(robust = TRUE,  scatter = "moment"),
            M2 = list(robust = TRUE,  scatter = "sign"),
            M3 = list(robust = FALSE, scatter = "moment"))

grid <- expand.grid(rep = seq_len(n_reps), method_id = names(VAR),
                    scenario_id = SCEN, stringsAsFactors = FALSE)
grid <- grid[order(grid$scenario_id, grid$method_id, grid$rep), ]
rownames(grid) <- NULL
mine <- c_slice(grid, task_id, n_tasks)
if (!nrow(mine)) { cat("nothing to do\n"); quit(save = "no", status = 0) }
c_header(sprintf("C4 at nK_fit = %s: cosmetic or substantive", cfg$global$nK_fit),
         nrow(mine), nrow(grid), task_id, n_tasks)
cat(sprintf("  scenarios %s | conv.prop %s | %d reps | n = %s, obs_min = %s\n",
            paste(SCEN, collapse = ","), paste(CP, collapse = ","), n_reps,
            cfg$global$n, cfg$global$obs_min))

started <- Sys.time()
out <- list()
for (a in seq_len(nrow(mine))) {
  cell <- mine[a, ]
  sc <- cfg$scenarios[cfg$scenarios$scenario_id == cell$scenario_id, ]
  v  <- VAR[[cell$method_id]]
  t0 <- Sys.time()

  rows <- tryCatch({
    dat <- generate_replicate(cfg, sc, rep = cell$rep)
    ref <- NULL                       # the conv.prop = 0 partition

    do.call(rbind, lapply(CP, function(cp) {
      # Identical RNG state for every threshold, so a partition difference is
      # the THRESHOLD and never the draw. `two-staged` initialisation draws
      # twice -- k-means centres, then the random fallback -- so this matters.
      c_seed_fit(cfg, cell$scenario_id, cell$rep)
      s0 <- Sys.time()
      fit <- tryCatch(MRKCFC::MKCFC(
        Ly = dat$Ly, Lt = dat$Lt, K = K, robust = v$robust,
        optns = build_optns(cfg, v$robust), mfpca_scatter = v$scatter,
        init.method = cfg$global$init_method,
        maxIter = as.integer(cfg$global$maxIter),
        parallel = FALSE, trace = 0, conv.prop = cp,
        LOO = isTRUE(as.logical(cfg$global$LOO))), error = function(e) e)
      el <- as.numeric(difftime(Sys.time(), s0, units = "secs"))

      if (inherits(fit, "error"))
        return(data.frame(conv_prop = cp, aRI = NA_real_, CCR = NA_real_,
                          agree_strict = NA_real_, converged = NA,
                          iter = NA_integer_, seconds = el,
                          error = conditionMessage(fit), stringsAsFactors = FALSE))

      lab <- as.integer(fit$res_clust$g)
      if (cp == 0) ref <<- lab        # first, by the stopifnot above

      data.frame(conv_prop = cp, aRI = ari(dat$labels, lab),
                 CCR = ccr(dat$labels, lab),
                 # THE STATISTIC THIS STUDY EXISTS FOR: agreement with the
                 # strict fit on the SAME replicate. 1.0 means the threshold
                 # changed nothing and the non-convergence is cosmetic.
                 agree_strict = if (is.null(ref)) NA_real_ else ari(ref, lab),
                 converged = isTRUE(fit$conv),
                 iter = if (is.null(fit$iter)) NA_integer_ else as.integer(fit$iter),
                 seconds = el, error = NA_character_, stringsAsFactors = FALSE)
    }))
  }, error = function(e) data.frame(
    conv_prop = NA_real_, aRI = NA_real_, CCR = NA_real_, agree_strict = NA_real_,
    converged = NA, iter = NA_integer_, seconds = NA_real_,
    error = paste("cell failed:", conditionMessage(e)), stringsAsFactors = FALSE))

  rows$scenario_id <- cell$scenario_id
  rows$arm <- sc$arm
  rows$arm_label <- sc$arm_label
  rows$d <- sc$d
  rows$method_id <- cell$method_id
  rows$rep <- cell$rep
  rows$nK_fit <- as.integer(cfg$global$nK_fit)
  rows$maxIter <- as.integer(cfg$global$maxIter)
  out[[length(out) + 1L]] <- rows
  cat(sprintf("  [%3d/%3d] %-4s %s rep %3d  %7.1fs  conv %s\n", a, nrow(mine),
              cell$scenario_id, cell$method_id, cell$rep,
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              paste(ifelse(is.na(rows$converged), "?",
                           ifelse(rows$converged, "T", "F")), collapse = "")))
  utils::flush.console()
}

results <- do.call(rbind, out)

# A completeness count, not a summary. Everything a summary would need spans
# tasks; read it with R/combine_c4b.R.
cat(sprintf("\ncells %d | rows %d | errors %d | %.1f min\n",
            nrow(mine), nrow(results), sum(!is.na(results$error)),
            as.numeric(difftime(Sys.time(), started, units = "mins"))))
c_write(results, outdir, "c4b", task_id, started)
cat("Read the study with R/combine_c4b.R -- this file does not summarise.\n")
