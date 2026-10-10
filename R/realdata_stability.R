# Study D1: is the real-data partition stable to initialisation, and to
# subsampling the cohort?
#
#   Rscript R/realdata_stability.R <results-dir> <combined-stem> [variant]
#
# e.g.  Rscript R/realdata_stability.R .../results_rd .../rd  MRKCFC-R2
#
# ===========================================================================
# REWRITTEN 2026-09-09. IT REFITS NOTHING.
#
# An earlier version ran 100 fresh random starts plus 500 subsample refits; at
# the cohort's measured cost that is on the order of 600 hours, and it could
# not have run.
#
# THE ARRAY ALREADY RAN THE EXPERIMENT. Every cell of Study D is ten independent
# random starts on identical data, which is exactly what initialisation
# stability asks about, and each returns the best iterate of its own trajectory
# -- the same quantity the paper reports. So the answer is a read, not a run.
#
# TWO QUESTIONS, AND THEY FAIL DIFFERENTLY. This file reports both, but it
# COMPUTES only the first:
#
#   INITIALISATION   read from the array's ten random starts per cell. No
#                    refitting. Instability is an ALGORITHMIC statement: the
#                    loop finds different local optima on identical input.
#   SUBSAMPLE        read from this arm's results_sub<suffix>/ (or RD_SUB_OUT),
#                    written by realdata_subsample_task.R, which must be
#                    run separately because a refit costs 40 min to 3.3 h;
#                    draws from another arm's design are detected and ignored.
#                    Instability is a statement about the DATA: the structure
#                    does not survive dropping a fifth of the cohort.
#
# A partition can be perfectly initialisation-stable and subsample-unstable.
# If the draws are absent this file says so, and the second must not then be
# claimed.
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L)
  stop("usage: realdata_stability.R <results-dir> <combined-stem> [variant]",
       call. = FALSE)
inp     <- args[1]
stem    <- args[2]
primary <- if (length(args) >= 3) args[3] else "MRKCFC-R2"

source(file.path("R", "metrics.R"))       # ari()

cf <- paste0(stem, "_combined.rds")
if (!file.exists(cf))
  stop("run realdata_combine.R first; no ", cf, call. = FALSE)
C <- readRDS(cf)
if (is.na(C$K_star))
  stop("the combine did not determine K; nothing to assess stability of.",
       call. = FALSE)
K <- as.integer(C$K_star)
cat(sprintf("Study D1: initialisation stability at K = %d, bwK = %s, %d restarts per cell\n",
            K, format(C$bwK), C$n_restart))

# ------------------------------------------------------- the restarts ------
# Canonicalise by first appearance so two partitions differing only by a label
# permutation compare as identical. aRI is already permutation-invariant; this
# is for the modal-partition count.
canon <- function(l) as.integer(factor(l, levels = unique(l)))

# Cluster-wise Jaccard, not a global index. A GLOBAL aRI CAN LOOK HEALTHY WHILE
# ONE SMALL CLUSTER DISSOLVES ENTIRELY, and a small unstable cluster is exactly
# the thing a reader must not build an interpretation on. Hennig's convention:
# below 0.5 dissolved, 0.6-0.75 a pattern rather than a stable cluster, above
# 0.85 highly stable.
jac_row <- function(ref, other, K) {
  vapply(seq_len(K), function(k) {
    A <- which(ref == k)
    if (!length(A)) return(NA_real_)
    max(vapply(sort(unique(other)), function(l) {
      B <- which(other == l)
      length(intersect(A, B)) / length(union(A, B))
    }, numeric(1)), na.rm = TRUE)
  }, numeric(1))
}
verdict <- function(j) if (!is.finite(j)) "---" else
  if (j < 0.50) "DISSOLVED" else if (j < 0.60) "unstable" else
  if (j < 0.75) "pattern, not a stable cluster" else
  if (j < 0.85) "stable" else "highly stable"

assess <- function(variant) {
  d <- C$info[C$info$variant == variant & C$info$K == K, ]
  n_attempt <- nrow(d)
  d <- d[d$usable, ]
  if (nrow(d) < 2L) {
    cat(sprintf("\n=== %s ===\n  %d of %d restarts usable -- too few to assess.\n",
                variant, nrow(d), n_attempt))
    return(NULL)
  }
  L <- lapply(file.path(inp, d$file), function(f) as.integer(readRDS(f)$labels))
  # THE REFERENCE IS THE REPORTED PARTITION -- the best-of-N winner, not an
  # arbitrary member -- because the question is whether THAT partition is
  # reproducible, not whether the restarts resemble each other in general.
  b   <- C$best[C$best$variant == variant & C$best$K == K, ]
  ref <- as.integer(readRDS(file.path(inp, b$file))$labels)

  cb <- utils::combn(length(L), 2)
  pw <- vapply(seq_len(ncol(cb)),
               function(i) ari(L[[cb[1, i]]], L[[cb[2, i]]]), numeric(1))
  to_ref <- vapply(L, function(l) ari(ref, l), numeric(1))
  tb <- sort(table(vapply(L, function(l) paste(canon(l), collapse = ""),
                          character(1))), decreasing = TRUE)
  J <- do.call(rbind, lapply(L, function(l) jac_row(ref, l, K)))

  cat(sprintf("\n=== %s ===\n", variant))
  cat(sprintf("  %d of %d restarts usable (%d pairs)\n",
              length(L), n_attempt, ncol(cb)))
  cat(sprintf("  pairwise aRI among restarts:      median %.3f  min %.3f  IQR %.3f\n",
              stats::median(pw), min(pw), stats::IQR(pw)))
  cat(sprintf("  aRI against the reported partition: median %.3f  min %.3f\n",
              stats::median(to_ref), min(to_ref)))
  cat(sprintf("  modal partition: %d of %d starts (%.0f%%), %d distinct partitions\n",
              tb[1], length(L), 100 * tb[1] / length(L), length(tb)))
  cat(sprintf("\n  %-8s %6s %8s %8s   %s\n", "cluster", "n", "mean J", "sd",
              "verdict (Hennig)"))
  jm <- colMeans(J, na.rm = TRUE)
  for (k in seq_len(K))
    cat(sprintf("  %-8d %6d %8.3f %8.3f   %s\n", k, sum(ref == k), jm[k],
                stats::sd(J[, k], na.rm = TRUE), verdict(jm[k])))

  if (stats::median(to_ref) < 0.9)
    cat("\n  WARNING: the reported partition is NOT reliably recovered from\n",
        "  random starts. Say so and qualify the interpretation; do not present\n",
        "  it as if it were the unique solution.\n", sep = "")
  if (any(jm < 0.6, na.rm = TRUE))
    cat("\n  ACTION: at least one cluster falls below Hennig's 0.6. Report the\n",
        "  value, NAME the cluster, and qualify every claim resting on it.\n",
        "  The manuscript reports the values and names the clusters.\n", sep = "")

  list(variant = variant, K = K, n_ok = length(L), n_attempt = n_attempt,
       pairwise_ari = pw, ari_to_reported = to_ref,
       modal_share = as.numeric(tb[1]) / length(L), n_distinct = length(tb),
       jaccard = J, jaccard_mean = jm, ref = ref)
}

res <- Filter(Negate(is.null), lapply(unique(C$info$variant), assess))
names(res) <- vapply(res, `[[`, "", "variant")

# The three variants share a stage-one MFPCA, so agreement between them is not
# independent evidence -- but a variant that is markedly LESS stable than the
# others is a real difference, and it is the kind a reader will look for given
# the paper's robustness claim.
if (length(res) >= 2L) {
  cat("\n=== across variants, at K =", K, "===\n")
  cat(sprintf("  %-12s %10s %12s %12s\n", "variant", "median aRI",
              "modal share", "min mean J"))
  for (r in res)
    cat(sprintf("  %-12s %10.3f %11.0f%% %12.3f\n", r$variant,
                stats::median(r$pairwise_ari), 100 * r$modal_share,
                min(r$jaccard_mean, na.rm = TRUE)))
}

# ------------------------------------------------- subsample stability -----
#
# THE OTHER HALF, and it means something different. Everything above is
# stability to INITIALISATION -- an algorithmic statement. This is stability to
# the SAMPLE: whether the structure survives dropping a fifth of the cohort,
# which is a statement about the DATA. A partition can be perfectly
# initialisation-stable and subsample-unstable.
#
# Produced by realdata_subsample_task.R, one draw per task, because a
# refit costs 40 min to 3.3 h. Absent unless that has run.
#
# WHICH DIRECTORY. The default follows the task script's convention: the
# primary arm's draws sit in results_sub beside results_rd, and every other
# arm's in results_sub<suffix> beside results_rd<suffix> (results_rd_log ->
# results_sub_log, results_rd_log_nK3 -> results_sub_log_nK3). RD_SUB_OUT
# overrides. UNTIL 2026-09-14 THE DEFAULT WAS results_sub FOR EVERY ARM, so
# the log arm's M = 3 run once read the primary arm's hundred draws
# and reported them as its own: rd_log_nK3_subsample.csv and the `subsample`
# element of rd_log_nK3_stability.rds from that run are the PRIMARY arm's
# numbers (0.77 / 0.64). No table or sentence used them.
#
# AND A GUARD. Every draw records the design it was refitted under (nK, log,
# window, time axis, bwK, K; realdata_subsample_task.R). Draws whose design
# differs from this combine's are announced and IGNORED, so a wrong directory
# yields "initialisation only" rather than another arm's answer. Draws written
# before the design fields existed cannot be checked and are accepted with a
# note naming the directory.
sub_default <- {
  b <- basename(inp)
  file.path(dirname(inp), if (startsWith(b, "results_rd"))
    paste0("results_sub", substring(b, nchar("results_rd") + 1L)) else "results_sub")
}
subdir <- Sys.getenv("RD_SUB_OUT", sub_default)
sfs <- if (dir.exists(subdir))
  list.files(subdir, pattern = "^rd_sub_[0-9]+[.]rds$", full.names = TRUE) else
  character(0)
sub_summary <- NULL
D <- if (length(sfs)) lapply(sfs, readRDS) else list()
ignored <- FALSE
if (length(D)) {
  design <- list(nK = C$nK, log = C$log, window = C$window, time = C$time,
                 bwK = C$bwK, K = K)
  cmp <- vapply(D, function(z) vapply(names(design), function(f)
    if (is.null(z[[f]])) NA else identical(format(z[[f]]), format(design[[f]])),
    NA), logical(length(design)))
  cmp <- matrix(cmp, nrow = length(design), dimnames = list(names(design), NULL))
  bad <- !is.na(cmp) & !cmp
  if (any(bad)) {
    i <- which(colSums(bad) > 0)[1]
    cat(sprintf("\n=== subsample stability ===\n  the %d draws in %s were refitted under a DIFFERENT design:\n",
                length(D), subdir))
    for (f in rownames(bad)[bad[, i]])
      cat(sprintf("    %-7s draws %-8s this combine %s\n", f,
                  format(D[[i]][[f]]), format(design[[f]])))
    cat("  They are another arm's and are IGNORED. Point RD_SUB_OUT at this arm's\n",
        "  draws, or run realdata_subsample_task.R with RD_SUB_RES,\n",
        "  RD_SUB_STEM and RD_SUB_OUT set for this arm. Until then the stability\n",
        "  answer covers INITIALISATION ONLY.\n", sep = "")
    D <- list(); ignored <- TRUE
  } else if (all(is.na(cmp))) {
    cat(sprintf("\n  NOTE: the draws in %s record no design fields and cannot be checked\n  against this combine; they are taken to be this arm's.\n", subdir))
  }
}
if (!length(D)) {
  if (!ignored) {
    cat(sprintf("\n=== subsample stability ===\n  no draws in %s\n", subdir))
    cat("  Run realdata_subsample_task.R for each draw. Until then the stability\n",
        "  answer covers INITIALISATION ONLY.\n", sep = "")
  }
} else {
  mode <- unique(vapply(D, function(z) z$init, ""))
  frac <- unique(vapply(D, function(z) z$frac, numeric(1)))
  vsub <- unique(vapply(D, function(z) z$variant, ""))
  if (length(mode) > 1L || length(frac) > 1L || length(vsub) > 1L)
    stop("the draws in ", subdir, " mix variant/frac/init; separate them.",
         call. = FALSE)
  okd <- vapply(D, function(z) !is.null(z$labels), TRUE)
  J <- do.call(rbind, lapply(D[okd], `[[`, "jaccard"))
  cat(sprintf("\n=== subsample stability: %s, %d draws at %.0f%%, %s start ===\n  from %s\n",
              vsub, length(D), 100 * frac, mode, subdir))
  # THE DENOMINATOR, beside the rate. A refit that fails is not a missing
  # observation: an 80%% draw shrinks every class by a fifth, and a class too
  # thin for the bandwidth is itself the finding.
  cat(sprintf("  %d of %d draws produced a partition\n", sum(okd), length(D)))
  if (any(!okd)) {
    tb <- table(vapply(D[!okd], function(z) z$status, ""))
    cat(sprintf("  failures: %s\n",
                paste(sprintf("%s:%d", names(tb), as.integer(tb)),
                      collapse = " ")))
    msg <- unique(stats::na.omit(vapply(D[!okd], function(z) z$error, "")))
    if (length(msg)) cat(sprintf("    %s\n", msg[1]))
    cat("  A high failure rate at this K IS A RESULT: the smallest class does\n",
        "  not survive a subsample at this bandwidth. Report it.\n", sep = "")
  }
  if (sum(okd) >= 2L) {
    jm <- colMeans(J, na.rm = TRUE)
    cat(sprintf("\n  %-8s %8s %8s   %s\n", "cluster", "mean J", "sd",
                "verdict (Hennig)"))
    for (k in seq_len(ncol(J)))
      cat(sprintf("  %-8d %8.3f %8.3f   %s\n", k, jm[k],
                  stats::sd(J[, k], na.rm = TRUE), verdict(jm[k])))
    if (mode == "reference")
      cat("\n  These are CONDITIONAL on the reported partition: each refit began\n",
          "  from it, so initialisation noise is removed and the numbers are\n",
          "  OPTIMISTIC relative to a full-procedure bootstrap. Say so.\n", sep = "")
    sub_summary <- data.frame(variant = vsub, cluster = seq_len(ncol(J)),
                              jaccard_mean = jm,
                              jaccard_sd = apply(J, 2, stats::sd, na.rm = TRUE),
                              n_ok = sum(okd), n_draw = length(D),
                              frac = frac, init = mode,
                              verdict = vapply(jm, verdict, ""),
                              stringsAsFactors = FALSE)
    utils::write.csv(sub_summary, paste0(stem, "_subsample.csv"),
                     row.names = FALSE)
  }
}

out <- file.path(dirname(stem), paste0(basename(stem), "_stability.rds"))
saveRDS(list(K = K, bwK = C$bwK, n_restart = C$n_restart, primary = primary,
             subsample = sub_summary,
             per_variant = res, when = Sys.time()), out)
tab <- do.call(rbind, lapply(res, function(r) data.frame(
  variant = r$variant, K = r$K, n_ok = r$n_ok, n_attempt = r$n_attempt,
  median_pairwise_ari = stats::median(r$pairwise_ari),
  min_pairwise_ari = min(r$pairwise_ari),
  median_ari_to_reported = stats::median(r$ari_to_reported),
  modal_share = r$modal_share, n_distinct = r$n_distinct,
  cluster = seq_len(r$K), cluster_n = as.integer(table(r$ref)),
  jaccard_mean = r$jaccard_mean,
  verdict = vapply(r$jaccard_mean, verdict, ""),
  stringsAsFactors = FALSE)))
utils::write.csv(tab, paste0(stem, "_stability.csv"), row.names = FALSE)
cat(sprintf("\nwrote %s and %s_stability.csv\n", out, stem))

# ===========================================================================
# THE SUBSAMPLE QUESTION, WHICH THIS FILE DOES NOT ANSWER
#
# Everything above is stability to INITIALISATION. Whether the structure
# survives dropping a fifth of the cohort is a different question with a
# different failure meaning -- instability there is a statement about the DATA,
# not the algorithm -- and a partition can be perfectly initialisation-stable
# and subsample-unstable. Reporting only the first would mislead.
#
# BUILT 2026-09-09 as R/realdata_subsample_task.R, one draw per task, and read
# back by the section above.
#
# It cannot be a loop in this script. A refit costs 40 min to 3.3 h at 6 cores,
# so the previous version's 500 draws was roughly 600 core-hours pretending to
# be a for-loop. Done properly it is an array on the pattern of
# realdata_fit_task.R -- one draw per task, checkpointed -- and at 100 draws x
# 6 cores it is affordable now in a way it was not this morning.
#
# If it is run, size it from `rd_restarts.csv`: the per-fit `hours` column
# prices it exactly, and the `sizes` column says whether the smallest cluster is
# large enough for an 80 percent draw to leave the bandwidth viable. That second
# check matters -- a class of 82 at K = 4 becomes about 66 in a draw, and the
# bandwidth of 0.20 is sized for classes far larger.
# ===========================================================================
