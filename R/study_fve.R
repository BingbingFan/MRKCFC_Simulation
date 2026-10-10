# Fraction of variance explained by the first M components, per variable, per
# scenario, for both the robust and the non-robust stage-one FPCA.
#
#   Rscript R/study_fve.R <task> [n_tasks] [outdir] [reps] [n] [obs_min]
#
# Defaults: reps 20, n from the config, obs_min 4.
#
#   Rscript R/study_fve.R 1 60 results_fve 20 500 4
#
# WHY THIS EXISTS
#
# The stage-1 sweep settles which M performs best. It cannot settle which M a
# READER could have chosen, because aRI needs the labels a practitioner does not
# have. FVE uses only the observed curves, so "retain the smallest M with
# cumulative FVE >= 90 percent" can be stated in advance, applied identically to
# every method, and reproduced from the data alone. That converts M from a tuned
# constant into an estimated quantity.
#
# It matters more after stage 1 than before it. The sweep found that EVERY
# method gains from raising M, so the component count is a property of the
# design; and that MRKCFC-R2's convergence falls to 0.05 at M = 5 and 0.00 at
# M = 6, so M = 4 is the largest usable value. If the FVE rule independently
# lands on 4, the choice is made by a label-free criterion that happens to
# coincide with the convergence boundary -- much stronger than either argument
# alone. If it lands elsewhere, report the FVE choice and demote the sweep to
# sensitivity. DO NOT move the threshold to make them agree.
#
# obs_min DEFAULTS TO 4, matching the revised design and the sweep. Running at
# the config's obs_min = 2 would compute the rule on DIFFERENT DATA from the
# results it is meant to justify -- the schedules differ (mean 6.85 against
# 7.20 observations per subject), so the spectra differ.
#
# TWO THINGS TO KNOW ABOUT THE NUMBERS
#
# 1. KFPCA discarded its eigenvalues until 2026-09-01, so the robust branch had
#    no FVE at all and M could not have been chosen from the data on the arms
#    the paper is about. `KFPCA()` now returns `$values`.
#
# 2. NEGATIVE EIGENVALUES ARE COMMON AND MUST BE TRUNCATED. A smoothed
#    PACE-type covariance surface is not guaranteed positive semi-definite; on a
#    typical fit here 25 of 51 eigenvalues are negative. Cumulating without
#    truncation gives shares that overshoot and then fall back below one, which
#    would make any FVE threshold behave unpredictably. Everything below uses
#    pmax(values, 0).
#
# THE FIT IS GLOBAL, over all subjects, not per cluster. That is deliberate: M
# has to be chosen BEFORE the clustering that would define the clusters, so a
# global fit is what a practitioner can actually do. The method itself refits
# per cluster inside the loop, so the within-cluster FVE will differ somewhat --
# state that rather than implying the two are the same.
#
# THIS DRIVER WRITES; IT DOES NOT SUMMARISE. An earlier version printed median
# FVE tables computed from `results`, which is only THIS TASK'S SLICE -- correct
# at n_tasks = 1 and silently wrong at any array size. The smoke test at task 1
# of 12 printed all twelve scenarios with A01 populated and A02 to A12 showing
# `>6` at every threshold; those were empty cells formatted like results. Read
# the output with R/combine_fve.R, which sees every task.

suppressMessages(library(MRKCFC))
source(file.path("R", "generate.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "methods.R"))
for (f in c("gen_scores", "vm_intermediate", "rmvn", "FLEISHMAN", "IQR_SCALE"))
  assign(f, get(f, envir = asNamespace("MRKCFC")), envir = globalenv())
KFPCA <- get("KFPCA", envir = asNamespace("MRKCFC"))

args    <- commandArgs(trailingOnly = TRUE)
task_id <- if (length(args) >= 1) as.integer(args[1]) else 1L
n_tasks <- if (length(args) >= 2) as.integer(args[2]) else 1L
outdir  <- if (length(args) >= 3) args[3] else "results_fve"
n_reps  <- if (length(args) >= 4) as.integer(args[4]) else 20L
n_over  <- if (length(args) >= 5) as.integer(args[5]) else NA_integer_
obs_min <- if (length(args) >= 6) as.integer(args[6]) else 4L
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cfg <- load_config()
if (!is.na(n_over)) cfg$global$n <- n_over
cfg$global$obs_min <- obs_min
P      <- as.integer(cfg$global$P)
M_MAX  <- 6L

# Cumulative FVE at 1..M_MAX from a raw eigenvalue spectrum, negatives dropped.
cum_fve <- function(v) {
  v <- pmax(v, 0); s <- sum(v)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, M_MAX))
  cs <- cumsum(v) / s
  c(cs, rep(1, M_MAX))[seq_len(M_MAX)]
}

grid <- expand.grid(rep = seq_len(n_reps), scenario_id = cfg$scenarios$scenario_id,
                    stringsAsFactors = FALSE)
grid <- grid[order(grid$scenario_id, grid$rep), ]
rownames(grid) <- NULL
mine <- which((seq_len(nrow(grid)) - 1L) %% n_tasks + 1L == task_id)
if (!length(mine)) { cat("nothing to do\n"); quit(save = "no", status = 0) }
cat(sprintf("FVE | task %d of %d: %d cells of %d | n = %s, obs_min = %d, %d reps\n",
            task_id, n_tasks, length(mine), nrow(grid), cfg$global$n, obs_min, n_reps))
cat(sprintf("R %s | MRKCFC %s | host %s\n", getRversion(),
            utils::packageVersion("MRKCFC"), Sys.info()[["nodename"]]))

started <- Sys.time()
out <- list()
for (a in seq_along(mine)) {
  cell <- grid[mine[a], ]
  sc <- cfg$scenarios[cfg$scenarios$scenario_id == cell$scenario_id, ]
  t0 <- Sys.time()
  rows <- tryCatch({
    dat <- generate_replicate(cfg, sc, rep = cell$rep)
    ro <- build_optns(cfg, TRUE); no <- build_optns(cfg, FALSE)
    do.call(rbind, lapply(seq_len(P), function(j) {
      # robust: KFPCA, the branch MRKCFC and MRKCFC-R2 both use
      kf <- tryCatch(KFPCA(Lt = dat$Lt[[j]], Ly = dat$Ly[[j]],
                           interval = ro$interval, nK = ro$nK, bw = ro$bw,
                           bwK = ro$bwK, nRegGrid = ro$nRegGrid,
                           fdParobj = ro$basis, more = FALSE),
                     error = function(e) NULL)
      # non-robust: fdapace, the branch MKCFC uses. cumFVE comes ready made.
      fp <- tryCatch(fdapace::FPCA(dat$Ly[[j]], dat$Lt[[j]],
                                   optns = utils::modifyList(no, list(methodSelectK = M_MAX))),
                     error = function(e) NULL)
      rbind(
        data.frame(branch = "robust (KFPCA)", variable = j, M = seq_len(M_MAX),
                   FVE = if (is.null(kf)) rep(NA_real_, M_MAX) else cum_fve(kf$values),
                   stringsAsFactors = FALSE),
        data.frame(branch = "non-robust (fdapace)", variable = j, M = seq_len(M_MAX),
                   FVE = if (is.null(fp)) rep(NA_real_, M_MAX)
                         else c(fp$cumFVE, rep(1, M_MAX))[seq_len(M_MAX)],
                   stringsAsFactors = FALSE))
    }))
  }, error = function(e) data.frame(
    branch = NA_character_, variable = NA_integer_, M = NA_integer_,
    FVE = NA_real_, stringsAsFactors = FALSE))

  rows$scenario_id <- cell$scenario_id
  rows$arm <- sc$arm
  rows$rep <- cell$rep
  rows$n <- as.integer(cfg$global$n)
  rows$obs_min <- obs_min
  out[[length(out) + 1L]] <- rows
  cat(sprintf("  [%3d/%3d] %-4s rep %3d  %6.1fs\n", a, length(mine),
              cell$scenario_id, cell$rep,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  utils::flush.console()
}

results <- do.call(rbind, out)
results$task_id <- task_id
attr(results, "session") <- list(
  R = as.character(getRversion()), MRKCFC = as.character(utils::packageVersion("MRKCFC")),
  host = Sys.info()[["nodename"]], started = started, finished = Sys.time(),
  n = cfg$global$n, obs_min = obs_min, n_reps = n_reps)

f <- file.path(outdir, sprintf("fve_task_%04d.rds", task_id))
saveRDS(results, f)
utils::write.csv(results, sub("[.]rds$", ".csv", f), row.names = FALSE)

# A completeness count, not a summary: everything below is a property of THIS
# task and needs no sight of any other.
cat(sprintf("\ncells %d | rows %d | NA FVE %d | %.1f min\nwrote %s\n",
            length(mine), nrow(results), sum(is.na(results$FVE)),
            as.numeric(difftime(Sys.time(), started, units = "mins")), f))
cat("Read the study with R/combine_fve.R -- this file does not summarise.\n")
