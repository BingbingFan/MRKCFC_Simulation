# Study D1b: ONE subsample draw. One array task.
#
#   Rscript R/realdata_subsample_task.R <draw> <results-dir> <combined-stem> \
#           <outdir> [variant] [frac] [init]
#
# init   "reference" (default) or "random" -- see THE INITIALISATION CHOICE below
#
# ===========================================================================
# WHY THIS IS A SEPARATE QUESTION FROM realdata_stability.R
#
#   INITIALISATION stability (D1)  refit the SAME data from many random starts.
#                                  Instability is an ALGORITHMIC statement: the
#                                  loop finds different local optima on
#                                  identical input.
#   SUBSAMPLE stability (this)     refit 80 percent draws. Instability is a
#                                  statement about the DATA: the structure is
#                                  not supported by the cohort robustly enough
#                                  to survive dropping a fifth of it.
#
# A partition can be perfectly initialisation-stable and subsample-unstable.
# Reporting only the first would mislead, which is why both are reported.
#
# IT IS AN ARRAY BECAUSE IT HAS TO BE. A refit costs 40 min to 3.3 h at 6 cores,
# so the pre-data version's 500 draws inside a for-loop was roughly 600
# core-hours pretending to be a loop. One draw per task at %30 is affordable.
#
# THE INITIALISATION CHOICE, and it is a real one.
#
#   "reference"  start each refit from the reported partition restricted to the
#                drawn subjects. Initialisation noise is removed, so what is
#                measured is the DATA's support for the structure -- which is
#                the question D1 does not answer. It is CONDITIONAL stability
#                and it is optimistic relative to a full-procedure bootstrap:
#                the refit begins at the answer it is being asked about.
#   "random"     the full procedure, as Hennig's convention has it. Honest but
#                confounded: instability mixes the data question with the
#                initialisation question that D1 already answers separately.
#
# DEFAULT IS "reference", because D1 supplies the other half and the pair is
# more informative than either alone. THE REPORT MUST NAME WHICH WAS USED and
# state the optimism. Do not report a conditional stability as though it were a
# bootstrap.
#
# A BANDWIDTH HAZARD SPECIFIC TO THIS ANALYSIS. An 80 percent draw shrinks every
# class by a fifth. The refit inherits the combined object's time axis and
# bandwidth (warped time and GCV since 2026-09-10; see realdata_load.R), under
# which the robust smoother succeeded at every class size the probe tested --
# but a real class formed by the loop is less favourable than a random subset,
# and a draw of one may still fail on "No enough points in local window".
# THAT IS A RESULT, NOT A CRASH: a partition whose smallest class cannot survive
# a subsample is a partition to qualify. The failure is recorded with its
# denominator rather than dropped.
# ===========================================================================

source(file.path("R", "study_c.R"))       # load_config, build_optns,
                                          # reconstruct_given_partition
source(file.path("R", "realdata_load.R"))
source(file.path("R", "realdata_fitloop.R"))

args    <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L)
  stop("usage: realdata_subsample_task.R <draw> <results-dir> <combined-stem> <outdir> [variant] [frac] [init]",
       call. = FALSE)
draw    <- as.integer(args[1])
inp     <- args[2]
stem    <- args[3]
outdir  <- args[4]
variant <- if (length(args) >= 5) args[5] else "MRKCFC-R2"
frac    <- if (length(args) >= 6) as.numeric(args[6]) else 0.80
initmode<- if (length(args) >= 7) args[7] else "reference"
if (!initmode %in% c("reference", "random"))
  stop("init must be 'reference' or 'random'", call. = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ckpt <- file.path(outdir, sprintf("rd_sub_%04d.rds", draw))
if (file.exists(ckpt)) {
  old <- tryCatch(readRDS(ckpt), error = function(e) NULL)
  if (!is.null(old) && is.character(old$status) &&
      identical(old$variant, variant) && isTRUE(all.equal(old$frac, frac)) &&
      identical(old$init, initmode)) {
    cat(sprintf("draw %d already done (%s); nothing to do\n", draw, old$status))
    quit(save = "no", status = 0L)
  }
  cat(sprintf("IGNORING %s: it does not match this draw's variant/frac/init.\n", ckpt))
}

C <- readRDS(paste0(stem, "_combined.rds"))
if (is.na(C$K_star)) stop("the combine did not determine K.", call. = FALSE)
K <- as.integer(C$K_star)
if (!variant %in% names(C$labels))
  stop(variant, " has no reported partition at K = ", K, call. = FALSE)
ref_all <- as.integer(C$labels[[variant]])

VARIANTS <- list(MKCFC       = list(robust = FALSE, scatter = "moment"),
                 MRKCFC      = list(robust = TRUE,  scatter = "moment"),
                 `MRKCFC-R2` = list(robust = TRUE,  scatter = "sign"))
v <- VARIANTS[[variant]]

ncore <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
if (is.na(ncore) || ncore < 1L) ncore <- 1L

dat  <- rd_load(window = C$window)
prep <- rd_prepare(dat, log = C$log, verbose = FALSE, time = C$time)
cfg  <- load_config(); cfg$global$nK_fit <- C$nK
g    <- cfg$global
optn <- rd_opt(C$nK, v$robust, cfg,
               bwK = if (is.na(C$bwK) || identical(C$bwK, "GCV")) NULL else as.numeric(C$bwK))

stopifnot(length(ref_all) == prep$n)

# One draw, seeded on the draw index alone so a task is reproducible wherever it
# runs and independent of how the array is sliced.
set.seed(20260909L + 977L * draw)
m   <- max(1L, round(frac * prep$n))
idx <- sort(sample.int(prep$n, m))
ref <- ref_all[idx]

cat(sprintf("draw %04d: %s K=%d, %d of %d subjects (%.0f%%), init=%s, bwK=%s, %d core(s)\n",
            draw, variant, K, m, prep$n, 100 * frac, initmode,
            format(C$bwK), ncore))
cat(sprintf("  reference sizes in the draw: %s\n",
            paste(table(ref), collapse = "/")))

# Subsetting a P-LIST: element [[p]] is the per-variable list of n subjects, so
# the draw indexes INSIDE each element. `Ly[idx]` would take the first m
# VARIABLES, which is the bug that was in the pre-data version of
# realdata_stability.R and would have refitted on a subset of the markers.
sub_Ly <- lapply(prep$Ly, function(el) el[idx])
sub_Lt <- lapply(prep$Lt, function(el) el[idx])
sub_cfg <- cfg

objective <- function(lab) {
  r <- reconstruct_given_partition(sub_Ly, sub_Lt, lab, sub_cfg)
  sum((r$y - r$yhat)^2)
}

init <- if (initmode == "reference") {
  # Relabel to 1..K over the classes actually present in the draw; MKCFC
  # requires init.labels in 1..K with every class at ten or more members.
  as.integer(factor(ref, levels = sort(unique(ref))))
} else NULL
if (!is.null(init) && (length(unique(init)) < K || min(table(init)) < 10L)) {
  cat("  reference start is degenerate within this draw; falling back to random\n")
  init <- NULL
}

fit <- fit_best_iterate(
  Ly = sub_Ly, Lt = sub_Lt, K = K, robust = v$robust, scatter = v$scatter,
  optns = optn, objective = objective, init.labels = init,
  maxIter = as.integer(g$maxIter), conv.prop = as.numeric(g$conv_prop),
  LOO = isTRUE(as.logical(g$LOO)), ncore = ncore, chunk = 5L, trace = 1,
  label = sprintf("  draw %04d ", draw))

# CLUSTER-WISE JACCARD, computed on the drawn subjects only. For each reference
# cluster, the best Jaccard against any cluster of the refit. A global index
# would let one small cluster dissolve unnoticed.
jac <- rep(NA_real_, K)
if (!is.null(fit$labels)) {
  for (k in seq_len(K)) {
    A <- which(ref == k)
    if (!length(A)) next
    jac[k] <- max(vapply(sort(unique(fit$labels)), function(l) {
      B <- which(fit$labels == l)
      length(intersect(A, B)) / length(union(A, B))
    }, numeric(1)), na.rm = TRUE)
  }
  cat(sprintf("  Jaccard by reference cluster: %s\n",
              paste(sprintf("%.3f", jac), collapse = " ")))
} else {
  cat(sprintf("  NO PARTITION (%s): %s\n", fit$status,
              if (is.na(fit$error)) "" else fit$error))
}

out <- list(draw = draw, variant = variant, K = K, frac = frac, init = initmode,
            n_draw = m, n_total = prep$n, idx = idx, bwK = C$bwK, time = C$time,
            nK = C$nK,
            window = C$window, log = C$log,
            ref_sizes = as.integer(table(ref)),
            labels = fit$labels, jaccard = jac,
            status = fit$status, error = fit$error,
            iter = fit$iter, best_iter = fit$best_iter,
            objective = fit$objective, obj_trace = fit$obj_trace,
            seconds = fit$seconds, ncore = ncore, when = Sys.time())
tmp <- sprintf("%s.tmp%d", ckpt, Sys.getpid())
saveRDS(out, tmp); file.rename(tmp, ckpt)
cat(sprintf("wrote %s  [%s after %d iterations, %.2f h]\n",
            ckpt, fit$status, fit$iter, fit$seconds / 3600))
