# The component-count sweep: every scenario, M = 1..6, n in {500, 1000},
# obs_min = 4, all eight methods.
#
#   Rscript R/study_sweep.R <task> <n_tasks> <outdir> [reps] [M_list] [n_list] [obs_min]
#
# Defaults: reps 25, M "1,2,3,4,5,6", n "500", obs_min 4.
#
#   Rscript R/study_sweep.R 1 960 results_sweep 25 1,2,3,4,5,6 500 4
#   Rscript R/study_sweep.R 1 960 results_sweep 50 3,4        1000 4
#
# WHAT THIS ANSWERS
#
# Study C2 swept M for the MKCFC family on two cells and found MRKCFC-R2 rising
# from 0.003 at M = 2 to 0.702 at M = 4 on the hard cell. Three things were left
# open, and all three have to close before M can be written down in the paper:
#
#   1. Does it hold on all twelve scenarios, or only the two C cells?
#   2. Do the COMPARATORS gain as much? If they do, the M effect belongs to the
#      design -- the count was set too low for everyone -- and the ranking is
#      unchanged. mclust, k-means and teigen consume the same MFPCA scores, so
#      at M = 4 they receive 12 coordinates instead of 6.
#   3. Does the FVE >= 0.90 rule land on the same M that maximises aRI? If it
#      does, M is an ESTIMATED quantity rather than a tuned one, and the
#      tuning-on-the-truth objection goes away. `R/study_fve.R` supplies the
#      other half of that comparison.
#
# THREE THINGS THAT ARE HELD FIXED, and why
#
# M_true STAYS AT 2. The generating model is unchanged; only `nK_fit`, the count
# RETAINED when fitting, moves. M > 2 is deliberate over-retention, which is the
# whole question -- the second stage rotates P * M stacked scores, so M = 4 hands
# it twelve coordinates where the truth needs six.
#
# THE DATA IS IDENTICAL ACROSS M within a (scenario, n, replicate).
# generate_replicate() keys its seed on (seed_master, scenario_id, rep) and reads
# only M_true, so nothing it produces depends on nK_fit. The M comparison is
# therefore PAIRED: a difference between M = 2 and M = 4 is the fit, never the
# draw.
#
# PAIRING WITH STUDY A -- REWRITTEN 2026-09-04, THE OLD CLAIM IS DEAD. This
# said the (M = 2, n = 500) rows were a clean one-factor obs_min comparison
# against Study A, which was true only while Study A ran at obs_min = 2 and
# nK_fit = 2. Study A was re-run at obs_min = 4 and nK_fit = 4, so the sweep and
# Study A now AGREE on obs_min and there is no obs_min contrast left anywhere in
# this design.
#
# What survives is a REPLICATION CHECK at (M = 4, n = 500): same scenario, same
# seed, same obs_min, same M, differing only in replicate count (25 here against
# Study A's 100). Those rows should agree to sampling error. A systematic gap
# means one of the two runs is not the configuration it claims to be. Do NOT
# read any other cell as a comparison against Study A.
#
# THE lcmm SHORTCUT. m_multlcmm() reshapes the raw curves and never touches
# nK_fit, so its result is invariant to M -- flagged as `m_invariant` in
# methods.csv. It runs at the anchor M only, once per (scenario, n, replicate),
# and the other M levels carry no lcmm row. That is roughly 250 core-hours not
# spent re-deriving a known non-convergence six times. Read its rows as applying
# to every M.
#
#   This is only SAFE because run_all_methods() re-seeds before every method
#   (fixed 2026-08-28). Before that fix, dropping a method from the list moved
#   the RNG for every method after it, and gating lcmm would have silently
#   changed the other seven.
#
# WHY IT CALLS run_all_methods() UNCHANGED. Sharing the global MFPCA across the
# four score-then-cluster comparators would save about 2 percent. It would also
# create a second definition of what each method is, which is precisely how this
# project acquired two divergent copies of silhouette_score(). One code path is
# worth more than 2 percent.
#
# BUDGET, from measured numbers rather than guesses. Per replicate at n = 500,
# M = 2, all eight methods: 1483 s. Cost rises steeply with M because the loop
# stops converging -- on C_HARD, MRKCFC-R2 goes 1334 s at M = 2 to 5818 s at
# M = 4 with median iter pinned at the maxIter = 100 cap. Summed over M = 1..6
# that is about 6.4 core-hours per (scenario, replicate) at n = 500 and about
# 17.2 at n = 1000.
#
#   n = 500,  25 reps, M = 1..6   ~1,900 core-h    ~8 h on 250 CPUs
#   n = 500,  50 reps, M = 1..6   ~3,800 core-h   ~15 h
#   both n,   50 reps, M = 1..6  ~14,200 core-h   ~57 h
#
# RUN STAGE 1 FIRST. n = 500 settles which M the rules pick; n = 1000 at two M
# levels then costs a quarter of the blind factorial and answers a question that
# has been narrowed.

suppressMessages(library(MRKCFC))
source(file.path("R", "generate.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "methods.R"))
source(file.path("R", "separation.R"))
for (f in c("gen_scores", "vm_intermediate", "rmvn", "FLEISHMAN", "IQR_SCALE"))
  assign(f, get(f, envir = asNamespace("MRKCFC")), envir = globalenv())

args    <- commandArgs(trailingOnly = TRUE)
task_id <- if (length(args) >= 1) as.integer(args[1]) else 1L
n_tasks <- if (length(args) >= 2) as.integer(args[2]) else 1L
outdir  <- if (length(args) >= 3) args[3] else "results_sweep"
n_reps  <- if (length(args) >= 4) as.integer(args[4]) else 25L
M_list  <- if (length(args) >= 5) as.integer(strsplit(args[5], ",")[[1]]) else 1:6
n_list  <- if (length(args) >= 6) as.integer(strsplit(args[6], ",")[[1]]) else 500L
obs_min <- if (length(args) >= 7) as.integer(args[7]) else 4L
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cfg <- load_config()

# Methods come from the config, never from a list here: an exclusion has to be
# visible in methods.csv with its reason attached.
use_methods <- names(METHODS)
if ("include" %in% names(cfg$methods)) {
  keep <- cfg$methods$method_id[toupper(as.character(cfg$methods$include)) == "TRUE"]
  use_methods <- intersect(names(METHODS), keep)
}
invariant <- if ("m_invariant" %in% names(cfg$methods))
  intersect(use_methods, cfg$methods$method_id[
    toupper(as.character(cfg$methods$m_invariant)) == "TRUE"]) else character(0)
M_anchor <- min(M_list)          # the one M level that also runs lcmm

grid <- expand.grid(rep = seq_len(n_reps), M = M_list, n = n_list,
                    scenario_id = cfg$scenarios$scenario_id,
                    stringsAsFactors = FALSE)
grid <- grid[order(grid$scenario_id, grid$n, grid$M, grid$rep), ]
rownames(grid) <- NULL

# Strided, not contiguous. Cost varies by almost an order of magnitude between
# M = 1 and M = 4, so a contiguous slice would hand some tasks only cheap cells
# and others only expensive ones; striding mixes them and the array finishes
# together.
mine <- which((seq_len(nrow(grid)) - 1L) %% n_tasks + 1L == task_id)
if (!length(mine)) { cat("nothing to do\n"); quit(save = "no", status = 0) }

cat(sprintf("M sweep | task %d of %d: %d cells of %d\n",
            task_id, n_tasks, length(mine), nrow(grid)))
cat(sprintf("  M = %s | n = %s | obs_min = %d | %d reps | %d scenarios\n",
            paste(M_list, collapse = ","), paste(n_list, collapse = ","),
            obs_min, n_reps, nrow(cfg$scenarios)))
cat(sprintf("  methods: %s%s\n", paste(use_methods, collapse = ","),
            if (length(invariant))
              sprintf("  (%s at M = %d only, M-invariant)",
                      paste(invariant, collapse = ","), M_anchor) else ""))
cat(sprintf("R %s | MRKCFC %s | host %s\n", getRversion(),
            utils::packageVersion("MRKCFC"), Sys.info()[["nodename"]]))

started <- Sys.time()
out <- vector("list", length(mine))
for (a in seq_along(mine)) {
  cell <- grid[mine[a], ]
  sc <- cfg$scenarios[cfg$scenarios$scenario_id == cell$scenario_id, ]
  t0 <- Sys.time()

  # The three overrides that define this cell. nK_fit is the swept factor; n and
  # obs_min are the revised design. Everything else stays at the config value,
  # so a parameter that is not named here is the reportable one.
  cf <- cfg
  cf$global$nK_fit  <- cell$M
  cf$global$n       <- cell$n
  cf$global$obs_min <- obs_min

  run <- if (cell$M == M_anchor) use_methods else setdiff(use_methods, invariant)

  res <- tryCatch({
    dat <- generate_replicate(cf, sc, rep = cell$rep)
    rows <- run_all_methods(dat, cf, methods = run)
    sep <- separation_of(dat, cf)
    rows$delta <- sep$delta; rows$bayes_ccr <- sep$bayes
    rows
  }, error = function(e) data.frame(
    scenario_id = cell$scenario_id, rep = cell$rep, method = NA_character_,
    CCR = NA_real_, aRI = NA_real_, NMI = NA_real_, AMI = NA_real_,
    MPE = NA_real_, RRE = NA_real_, RRE_trim = NA_real_, converged = NA,
    iter = NA_integer_, seconds = NA_real_,
    error = paste("cell failed:", conditionMessage(e)),
    delta = NA_real_, bayes_ccr = NA_real_, stringsAsFactors = FALSE))

  res$M <- cell$M; res$n <- cell$n; res$obs_min <- obs_min
  res$arm <- sc$arm; res$arm_label <- sc$arm_label; res$d <- sc$d
  out[[a]] <- res
  cat(sprintf("  [%4d/%4d] %-4s n=%4d M=%d rep %3d  %7.1fs%s\n", a, length(mine),
              cell$scenario_id, cell$n, cell$M, cell$rep,
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              if (any(!is.na(res$error))) sprintf("  (%d err)", sum(!is.na(res$error))) else ""))
  utils::flush.console()
}

results <- do.call(rbind, out)
results$task_id <- task_id
attr(results, "session") <- list(
  R = as.character(getRversion()), MRKCFC = as.character(utils::packageVersion("MRKCFC")),
  host = Sys.info()[["nodename"]], started = started, finished = Sys.time(),
  M_list = M_list, n_list = n_list, obs_min = obs_min, n_reps = n_reps)

f <- file.path(outdir, sprintf("sweep_task_%04d.rds", task_id))
saveRDS(results, f)
utils::write.csv(results, sub("[.]rds$", ".csv", f), row.names = FALSE)
cat(sprintf("\nwrote %s  (%d rows, %.1f min)\n", f, nrow(results),
            as.numeric(difftime(Sys.time(), started, units = "mins"))))
