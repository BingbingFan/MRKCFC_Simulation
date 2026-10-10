# Study D: ONE (variant, K, restart) fit on the 36-month cohort. One array task.
#
#   Rscript R/realdata_fit_task.R <task_id> [outdir] [nK] [log] [n_restart] \
#                                 [bwK] [kmax] [time]
#
# PRIMARY ARM     90 tasks = 3 variants x K in {2,3,4} x 10 restarts,
#                 bwK GCV, warped time, raw marker values
# SENSITIVITY ARM 90 tasks, the same grid at a fixed bwK of 0.20, warped time
# LOG ARM         90 tasks, the primary settings with every marker log-
#                 transformed before normalisation (arg 4 = TRUE), 2026-09-11
#
# The two arms write to SEPARATE DIRECTORIES and are combined separately; the
# comparison is made on the shared K in {2,3} sub-grid. Read by
# realdata_combine.R, which selects best-of-N and then selects K.
#
# ===========================================================================
# REWRITTEN 2026-09-09, after a first run of 154 tasks produced not one usable
# fit. Five things it measured, and what each changed here.
#
# 1. LOO WAS THE COST BUG, AND IT WAS MINE.
#
#    params_global.csv sets LOO = FALSE and says why: "package default is TRUE,
#    manuscript never mentions it, and IT IS THE DOMINANT COST DRIVER". The
#    simulation reads it (methods.R:110). This driver hard-coded LOO = TRUE, so
#    the application was running a different and far more expensive estimator
#    than the simulation it exists to match. With LOO = TRUE the reassignment
#    step re-estimates class components once per subject per class -- n * K
#    component fits per iteration instead of K.
#
#    MEASURED CONSEQUENCE: a median of 153 minutes per iteration, range 58 to
#    350. Against maxIter = 100 that is 10.6 days per fit and a 48 h wall
#    permits about 19 iterations. No fit in the grid could ever have terminated.
#
#    maxIter, conv.prop, LOO and init_method now all come from the config, so
#    the application cannot drift from the simulation again. This is the same
#    correction that deleted rd_optns(): where the simulation has an answer, do
#    not hand-write a second one.
#
# 2. THE K GRID IS 2:4, AND THE BANDWIDTH IS THE REASON IT IS NOT 2:6.
#
#    Fits that had not degenerated after 38 h, out of 20 restarts per cell:
#
#      K = 2   59 of 59   100%      K = 5    8 of 60   13%
#      K = 3   51 of 60    85%      K = 6    2 of 60    3%
#      K = 4   23 of 60    38%
#
#    THAT TABLE DOES NOT MEAN WHAT IT LOOKS LIKE IT MEANS. An earlier bandwidth
#    check measured the smallest bandwidth leaving ten raw covariance pairs in every
#    window at each class size: 0.10 at K=2, 0.15 at K=3, 0.20 at K=4, 0.25-0.30
#    at K=5 and K=6. The run above used 0.10 throughout, so from K=3 upward the
#    ROBUST branch was being asked to smooth a covariance surface it could not
#    estimate -- and at K=5 and K=6 windows were literally EMPTY. All 112
#    bandwidth failures were MRKCFC or MRKCFC-R2 and none were MKCFC, because
#    build_optns() fixes bw only on the robust branch.
#
#    So "K >= 4 does not work on this cohort" IS WITHDRAWN for the robust
#    variants; they were never fairly tested there. What survives is narrower:
#    MKCFC degenerates at K >= 4 through the ten-curve guard, which is a
#    different mechanism and is unaffected by bw.
#
#    K in 2..4 at RD_BWK = 0.20 is the grid where every cell is estimable at one
#    common bandwidth. K = 5 and K = 6 would need 0.25 to 0.30, and raising the
#    bandwidth that far to chase them would smooth the class covariances toward
#    each other -- biasing K selection toward the small K it was raised to
#    reach. That is a circularity worth refusing.


#
# 3. THE DRIVER WROTE ONCE AT THE END, so a wall hit lost the whole fit.
#
#    The loop now runs in chunks of CHUNK accepted reassignments, carrying the
#    partition forward through init.labels, and writes a self-contained
#    checkpoint after each one. A wall hit now costs at most CHUNK iterations,
#    and re-submitting the same task RESUMES from the checkpoint rather than
#    restarting. That is what makes maxIter = 100 survivable without knowing
#    the per-iteration cost in advance.
#
#    WHAT CHUNKING COSTS: MKCFC declares convergence by comparing a candidate
#    against every earlier labelling in the SAME call, so a restart narrows
#    cycle detection to within-chunk plus the chunk boundaries checked below.
#    Fixed points are still caught exactly, since the chunk's own starting
#    labels are in its history. Longer cycles straddling a boundary out of
#    phase can be missed. At conv.prop = 0 the dominant outcome is exhausting
#    maxIter anyway, so this trades an unlikely miss for a recoverable run.
#
# 4. iter == 1 MEANS THE RANDOM START, UNTOUCHED.
#
#    MKCFC breaks BEFORE `final_labels <- candidate` (MKCFC.R:273), so a fit
#    whose first reassignment degenerates returns iteration 1 and the INITIAL
#    labelling. Fourteen fits did exactly this, and their cluster sizes are
#    uniform random splits of 1405: 271/271/289/276/298 and the like. Those are
#    not partitions, they are noise, and a balanced random split can score well
#    on squared reconstruction error precisely because it is balanced.
#
#    Such a fit is now recorded with status "init-only" AND with labels NULL
#    and objective NA, so a combiner that ignores status still cannot select
#    it. Defence in depth, because realdata_combine.R selects on the objective.
#
# 5. trace = 0 HID ALL OF THE ABOVE. It is 1 now, and every chunk prints its
#    own elapsed time, so per-iteration cost is visible from the first log
#    rather than reconstructed from wall times after the fact.
#
# 6. THE TIME AXIS IS WARPED, 2026-09-10, AND THE BANDWIDTH IS GCV AGAIN.
#
#    Two bandwidth fixes failed after this file was rewritten (a misattributed
#    `bw`, then a fixed `bwK` that failed 60 of 60 robust fits). The robust
#    smoother was then MEASURED directly (R/realdata_bwprobe.R): on the
#    rescaled axis it fails on the small, non-random classes the loop forms
#    because eighty percent of visits sit in the first quarter of the domain;
#    with time mapped through the empirical distribution of visit times it
#    succeeds at every class size and bandwidth tested, with GCV selecting
#    0.10-0.19. So `time` and `bwK` are arm arguments, the primary arm runs
#    warped/GCV, the sensitivity arm warped/0.20, and every checkpoint records
#    both. The table is in realdata_load.R.
# ===========================================================================
#
# THE DEFAULT INITIALISATION DOES NOT WORK ON THIS COHORT. MEASURED.
# `init.method = "kmeans"` is MKCFC()'s default and it ABORTS on the full 1405
# subjects with "Each cluster must contain at least 10 curves", and at n = 600
# too. Ward the same. Study C1 predicted it: over 460 simulation fits k-means
# aborted 14 times and Ward 19, while random starts aborted 0 of 400.
#
# BUT THE APPLICATION QUALIFIES C1 RATHER THAN CONFIRMING IT. Random starts do
# not abort at K = 2 or K = 3 here, and they fail almost always at K >= 4. The
# defensible claim is narrower than the one C1 alone supports: the k-means start
# fails outright, and random restarts work where the data supports the K being
# asked for. Section 3.2 must say the narrow version.
#
# EVERY RESTART IS SEEDED from (variant, K, restart), so a fit is reproducible
# wherever it runs and independent of how the array is sliced.

source(file.path("R", "study_c.R"))       # load_config, build_optns,
                                          # reconstruct_given_partition
source(file.path("R", "realdata_load.R"))

args      <- commandArgs(trailingOnly = TRUE)
task_id   <- if (length(args) >= 1) as.integer(args[1]) else 1L
outdir    <- if (length(args) >= 2) args[2] else "results_rd"
nK        <- if (length(args) >= 3) as.integer(args[3]) else 4L
uselog    <- if (length(args) >= 4) as.logical(args[4]) else FALSE
n_restart <- if (length(args) >= 5) as.integer(args[5]) else 10L
# THE ARMS. Defaults are the primary run: bwK by GCV, as the simulation has it,
# in WARPED time, K in 2..4. The sensitivity arm passes a fixed bwK of 0.20 so
# that K selection can be shown not to depend on the bandwidth. `time` exists
# so the rescaled axis stays runnable for comparison; the probe table in
# realdata_load.R is why warped is the default.
arm_bwK  <- if (length(args) >= 6) args[6] else "GCV"
if (!identical(arm_bwK, "GCV")) {
  arm_bwK <- as.numeric(arm_bwK)
  if (!is.finite(arm_bwK) || arm_bwK <= 0)
    stop("bwK must be 'GCV' or a positive number")
}
kmax     <- if (length(args) >= 7) as.integer(args[7]) else 4L
if (kmax < 2L) stop("kmax must be at least 2")
arm_time <- if (length(args) >= 8) args[8] else "warped"
if (!arm_time %in% c("warped", "rescaled"))
  stop("time must be 'warped' or 'rescaled'")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

WINDOW <- 36L
K_GRID <- 2:kmax             # see note 2; every cell estimable at this bandwidth
CHUNK  <- 5L                 # accepted reassignments between checkpoints
VARIANTS <- list(
  MKCFC       = list(robust = FALSE, scatter = "moment"),
  MRKCFC      = list(robust = TRUE,  scatter = "moment"),
  `MRKCFC-R2` = list(robust = TRUE,  scatter = "sign"))

grid <- expand.grid(restart = seq_len(n_restart), K = K_GRID,
                    variant = names(VARIANTS), stringsAsFactors = FALSE)
grid <- grid[order(match(grid$variant, names(VARIANTS)), grid$K, grid$restart), ]
rownames(grid) <- NULL
if (task_id < 1L || task_id > nrow(grid))
  stop("task_id must be 1..", nrow(grid))
cell <- grid[task_id, ]

ncore <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
if (is.na(ncore) || ncore < 1L) ncore <- 1L

cat(sprintf("task %d of %d: %s K=%d restart %02d | %d core(s) | nK=%d | %s\n",
            task_id, nrow(grid), cell$variant, cell$K, cell$restart, ncore, nK,
            if (uselog) "log" else "raw"))

ckpt <- file.path(outdir, sprintf("rd_task_%03d.rds", task_id))

# EARLY EXIT, BEFORE THE DATA IS TOUCHED. Re-submitting the whole array is the
# intended way to resume after a wall hit, so the common case is 60 tasks of
# which most are already finished. Reading the checkpoint before rd_load()
# keeps that from re-reading and re-preparing the cohort sixty times over to
# learn there was nothing to do.
# A FILE AT THIS PATH IS NOT NECESSARILY THIS TASK'S CHECKPOINT. An earlier
# 300-task run wrote the SAME rd_task_%03d.rds names into the SAME directory
# under a different grid, so its file 002 is a K=2 restart 02 fit of one design
# sitting where task 2 of another design will look. Those files carry no
# `status` field, and `!identical(NULL, "incomplete")` is TRUE, so an unguarded
# check would take them as terminal and exit -- printing nothing, because
# sprintf on NULL yields character(0), and reporting COMPLETED. Eighty-nine
# tasks would silently do nothing.
#
# So a checkpoint is only adopted when it agrees with this task on every field
# that changes what a fit MEANS. Anything else is announced and overwritten,
# because a stale file that is merely ignored comes back the next time someone
# forgets which directory they are in.
ckpt_mismatch <- function(p) {
  if (!is.character(p$status) || length(p$status) != 1L)
    return("no status field -- written by an earlier driver")
  want <- list(variant = cell$variant, K = as.integer(cell$K),
               restart = as.integer(cell$restart), nK = as.integer(nK),
               window = as.integer(WINDOW), log = as.logical(uselog))
  for (f in names(want)) {
    got <- p[[f]]
    if (is.null(got) || !identical(as.vector(got), want[[f]]))
      return(sprintf("%s differs: checkpoint %s, this task %s", f,
                     if (is.null(got)) "missing" else paste(got, collapse = ","),
                     paste(want[[f]], collapse = ",")))
  }
  # The time axis changes what every fit means, on both branches.
  if (!identical(as.character(p$time), arm_time))
    return(sprintf("time axis differs: checkpoint %s, this task %s",
                   if (is.null(p$time)) "missing" else as.character(p$time),
                   arm_time))
  # bwK binds on the ROBUST branch only, so a non-robust checkpoint carrying a
  # different bandwidth is still this task's and still valid. Compared as text
  # so "GCV" and a number are handled alike.
  if (isTRUE(VARIANTS[[cell$variant]]$robust) &&
      !identical(as.character(p$bwK), as.character(arm_bwK)))
    return(sprintf("bwK differs: checkpoint %s, this task %s",
                   if (is.null(p$bwK)) "missing" else as.character(p$bwK),
                   as.character(arm_bwK)))
  ""
}

prev <- NULL
if (file.exists(ckpt)) prev <- tryCatch(readRDS(ckpt), error = function(e) NULL)
if (!is.null(prev)) {
  why <- ckpt_mismatch(prev)
  if (nzchar(why)) {
    cat(sprintf("IGNORING the file at %s\n  %s\n  It will be overwritten; move it aside now if it is wanted.\n",
                ckpt, why))
    prev <- NULL
  }
}
if (!is.null(prev) && !identical(prev$status, "incomplete")) {
  cat(sprintf("checkpoint present and terminal (%s, %s iterations); nothing to do\n",
              prev$status,
              if (is.null(prev$iter)) "?" else format(prev$iter)))
  quit(save = "no", status = 0L)
}

dat  <- rd_load(window = WINDOW)
if (task_id == 1L) rd_floor_screen(dat, log = uselog)
prep <- rd_prepare(dat, log = uselog, verbose = task_id == 1L, time = arm_time)
cfg  <- load_config(); cfg$global$nK_fit <- nK
g    <- cfg$global

# FROM THE CONFIG, NOT FROM HERE. See note 1. The previous version hard-coded
# maxIter = 100L, conv.prop = 0 and LOO = TRUE; only the first two happened to
# match, and the third was the whole cost problem.
MAXITER   <- as.integer(g$maxIter)
CONV_PROP <- as.numeric(g$conv_prop)
USE_LOO   <- isTRUE(as.logical(g$LOO))
cat(sprintf("config: maxIter=%d conv.prop=%g LOO=%s  (all from params_global.csv)\n",
            MAXITER, CONV_PROP, USE_LOO))

v    <- VARIANTS[[cell$variant]]
# bwK: the config's GCV unless the arm passes a number. It binds on the ROBUST
# branch only. Both it and the time axis are recorded in the checkpoint so every
# result carries the operating point it was produced at.
optn <- rd_opt(nK, v$robust, cfg,
               bwK = if (identical(arm_bwK, "GCV")) NULL else arm_bwK)
cat(sprintf("bwK: %s | time axis: %s%s\n", as.character(arm_bwK), arm_time,
            if (v$robust) "" else " -- non-robust branch, fdapace selects its own bandwidth"))
seed <- 20260906L + 7919L * match(cell$variant, names(VARIANTS)) +
  131L * cell$K + cell$restart

# ---------------------------------------------------------------- resume ---
# `prev` was read above, before the data was loaded. Only an "incomplete"
# checkpoint reaches here, and it carries the partition the wall interrupted.
labels      <- NULL     # the CURRENT iterate: what the next chunk continues from
total_iter  <- 0L
chunks      <- 0L
secs        <- 0
history     <- list()
best_obj    <- Inf      # and the BEST iterate seen; see write_ckpt() below
best_labels <- NULL
best_iter   <- NA_integer_
trace_iter  <- integer(0)
trace_obj   <- numeric(0)
if (!is.null(prev) && (!is.null(prev$final_labels) || !is.null(prev$labels))) {
  # RESUME FROM final_labels, NOT labels. Since 2026-09-09 `labels` holds the
  # best iterate seen, which is not where the trajectory was; continuing from it
  # would silently restart the search from an earlier point.
  #
  # A checkpoint written between the chunking change and the best-iterate change
  # -- the 2026-09-09 probe wrote several -- has `labels` and no `final_labels`.
  # Falling back keeps those resumable; requiring `final_labels` would have
  # restarted them from zero without saying so.
  if (is.null(prev$final_labels))
    cat("NOTE: checkpoint predates best-iterate tracking; resuming from its `labels`.\n")
  labels      <- as.integer(if (!is.null(prev$final_labels)) prev$final_labels
                            else prev$labels)
  total_iter  <- as.integer(prev$iter)
  chunks      <- as.integer(prev$chunks)
  secs        <- as.numeric(prev$seconds)
  history     <- if (is.null(prev$history)) list(labels) else prev$history
  best_obj    <- if (is.null(prev$objective) || !is.finite(prev$objective)) Inf
                 else as.numeric(prev$objective)
  best_labels <- prev$labels
  best_iter   <- prev$best_iter
  if (!is.null(prev$obj_trace)) {
    trace_iter <- as.integer(prev$obj_trace$iter)
    trace_obj  <- as.numeric(prev$obj_trace$objective)
  }
  cat(sprintf("RESUMING from checkpoint at iteration %d (%d chunk(s), %.1f h, best %.1f at iter %s)\n",
              total_iter, chunks, secs / 3600, best_obj, format(best_iter)))
}

set.seed(seed + 1000L * chunks)   # distinct stream per resume, still deterministic

# Self-contained checkpoint. `objective` is recomputed every time so that any
# checkpoint, including one written just before a wall hit, is directly usable
# by the combiner. It costs about one LOO-free reassignment per chunk.
write_ckpt <- function(status, labels, err = NA_character_) {
  obj <- NA_real_
  if (!is.null(labels) && !identical(status, "init-only")) {
    obj <- tryCatch({
      r <- reconstruct_given_partition(prep$Ly, prep$Lt, labels, cfg)
      sum((r$y - r$yhat)^2)
    }, error = function(e) NA_real_)
  }
  # KEEP THE BEST ITERATE, NOT THE LAST, 2026-09-09.
  #
  # THE LOOP DOES NOT CONVERGE ON THIS COHORT -- IT OSCILLATES. Measured on the
  # probe, MKCFC K=2, objective at successive checkpoints:
  #
  #   iter   5  735869 | 10  756688 | 15  741132 | 20  748199
  #   iter  25  667667 | 30  747302 | 35  760056
  #
  # No trend, a 13.8 percent spread, and the smaller class swinging
  # 143/71/93/77/121/112/114. Monotone descent is not guaranteed here: the
  # refit step estimates TRUNCATED FPCA components per class rather than
  # minimising reconstruction error directly, so a Lloyd-style argument does
  # not close. At conv.prop = 0 nothing will satisfy the criterion, so the fit
  # runs to maxIter and stops wherever it happens to be.
  #
  # RETURNING THAT LAST ITERATE WOULD BREAK BEST-OF-N. The combiner selects the
  # minimum objective across restarts; if each restart returns an arbitrary
  # point in a +/-14 percent oscillation, that selects the luckiest STOPPING
  # POINT rather than the best PARTITION.
  #
  # So `labels` and `objective` are the BEST SEEN, and the final iterate is
  # carried separately. The objective is already computed at every checkpoint,
  # so this costs nothing beyond the comparison. `obj_trace` keeps the whole
  # sequence: it is the evidence for the non-convergence statement Section 3.2
  # now has to make, and it is cheap to carry.
  if (is.finite(obj)) {
    trace_iter <<- c(trace_iter, total_iter)
    trace_obj  <<- c(trace_obj, obj)
    if (obj < best_obj) {
      best_obj    <<- obj
      best_labels <<- labels
      best_iter   <<- total_iter
    }
  }
  out <- list(variant = cell$variant, K = cell$K, restart = cell$restart,
              seed = seed, ncore = ncore, nK = nK, log = uselog,
              window = WINDOW, ids = prep$ids, n = prep$n,
              t_range = prep$t_range,
              maxIter = MAXITER, conv.prop = CONV_PROP, LOO = USE_LOO,
              bwK = arm_bwK, bwK_binds = v$robust, time = arm_time,
              time_map = prep$time_map,
              # `labels` and `objective` are the BEST ITERATE SEEN -- this is
              # what realdata_combine.R reads. The final iterate is carried
              # beside it, and resume continues from `final_labels`.
              labels = best_labels,
              objective = if (is.finite(best_obj)) best_obj else NA_real_,
              best_iter = best_iter,
              final_labels = labels, final_objective = obj,
              obj_trace = data.frame(iter = trace_iter, objective = trace_obj),
              converged = identical(status, "converged"),
              iter = total_iter, chunks = chunks, seconds = secs,
              status = status, error = err, history = history,
              when = Sys.time())
  # PID IN THE STAGING NAME. A fixed ".tmp" is one path per TASK, not per
  # process, so two jobs running the same task interleave into it and rename a
  # corrupt file into place. That is not hypothetical -- two jobs were once
  # submitted for the same tasks on 2026-09-09, and re-submitting the array is
  # the intended way to resume, so the collision has a standing route back.
  tmp <- sprintf("%s.tmp%d", ckpt, Sys.getpid())
  saveRDS(out, tmp); file.rename(tmp, ckpt)
  cat(sprintf("checkpoint: status=%-10s iter=%-4d obj=%s best=%s@%s -> %s\n",
              status, total_iter,
              if (is.na(obj)) "NA" else sprintf("%.1f", obj),
              if (is.finite(best_obj)) sprintf("%.1f", best_obj) else "NA",
              format(best_iter), ckpt))
  invisible(out)
}

# ------------------------------------------------------------------ loop ---
#
# STATUS VOCABULARY, and which of these realdata_combine.R may select from:
#
#   converged    reached a fixed point or a repeat            USE
#   maxiter      exhausted the config's maxIter without one   USE, flag the rate
#   degenerate   guard fired after >= 1 accepted iteration    USE with care: a
#                real partition, but the loop could not go on
#   error        a fit error after >= 1 accepted iteration    USE with care
#   init-only    the guard fired on the FIRST reassignment;   NEVER
#                labels are NULL. MKCFC reports this exactly.
#   error-early  a fit error before any chunk was banked;     NEVER
#                labels are NULL. NOT the same as init-only: a
#                chunk is five reassignments inside one call,
#                and an error on its fourth loses the three
#                accepted ones with it. Until 2026-09-10 this
#                case was written as init-only, which made the
#                robust branch's mid-trajectory failures read
#                as failures of the random start.
#   incomplete   the wall interrupted it; resume, do not use  NEVER
#
# "maxiter" is a legitimate outcome, not a failure: the simulation runs the same
# conv.prop = 0 and maxIter = 100 and reports its non-convergence rate as a
# result. Excluding those fits here would make the application and the
# simulation report different things under the same name.
status <- "incomplete"
repeat {
  remaining <- MAXITER - total_iter
  if (remaining < 1L) { status <- "maxiter"; break }
  m <- min(CHUNK, remaining)

  t0 <- Sys.time()
  fit <- tryCatch(
    MRKCFC::MKCFC(Ly = prep$Ly, Lt = prep$Lt, K = cell$K, robust = v$robust,
                 optns = optn, mfpca_scatter = v$scatter,
                 init.method = "random", init.labels = labels,
                 maxIter = m + 1L,          # m+1 gives m accepted reassignments
                 parallel = ncore > 1L,
                 nmcores = if (ncore > 1L) ncore else NULL,
                 trace = 1, conv.prop = CONV_PROP, LOO = USE_LOO),
    error = function(e) e)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  secs <- secs + el
  chunks <- chunks + 1L

  if (inherits(fit, "error")) {
    cat(sprintf("ERROR in chunk %d after %.1f min: %s\n",
                chunks, el / 60, conditionMessage(fit)))
    # An error with iterations already banked still leaves a real partition.
    status <- if (total_iter == 0L) "error-early" else "error"
    write_ckpt(status, if (total_iter == 0L) NULL else labels,
               conditionMessage(fit))
    break
  }

  gained <- as.integer(fit$iter) - 1L      # accepted reassignments this chunk
  new    <- as.integer(fit$res_clust$g)
  cat(sprintf("chunk %d: +%d iter (total %d) %.1f min = %.1f min/iter | sizes %s\n",
              chunks, gained, total_iter + gained, el / 60,
              if (gained > 0) el / 60 / gained else NA_real_,
              paste(table(new), collapse = "/")))

  if (gained == 0L) {
    # The first reassignment of this chunk degenerated. With nothing banked,
    # `new` IS the random start -- see note 4.
    status <- if (total_iter == 0L) "init-only" else "degenerate"
    write_ckpt(status, if (total_iter == 0L) NULL else labels)
    break
  }

  total_iter <- total_iter + gained
  labels     <- new

  # Boundary-level cycle detection, recovering part of what chunking costs.
  repeated <- any(vapply(history, function(h) identical(h, labels), logical(1)))
  history[[length(history) + 1L]] <- labels

  if (isTRUE(fit$conv) || repeated) {
    status <- "converged"
    write_ckpt(status, labels); break
  }
  if (gained < m) {                        # guard fired mid-chunk
    status <- "degenerate"
    write_ckpt(status, labels); break
  }
  write_ckpt("incomplete", labels)         # chunk exhausted; keep going
}

if (status == "maxiter") write_ckpt("maxiter", labels)

cat(sprintf("DONE %s K=%d restart %02d: %s after %d iterations, %.2f h\n",
            cell$variant, cell$K, cell$restart, status, total_iter,
            secs / 3600))
# The spread is the reportable quantity, not a diagnostic aside: it is how far
# the returned partition could have been from the one the loop stopped on.
if (length(trace_obj) >= 2L) {
  cat(sprintf("objective over %d checkpoints: best %.1f at iter %s, last %.1f, spread %.1f%%\n",
              length(trace_obj), best_obj, format(best_iter),
              trace_obj[length(trace_obj)],
              100 * (max(trace_obj) - min(trace_obj)) / min(trace_obj)))
  if (!identical(best_labels, labels))
    cat("  the BEST iterate is not the last one; `labels` holds the best.\n")
}
