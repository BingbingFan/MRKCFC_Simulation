# Study D2: do the alternative methods find the same partition as the proposed
# method, and do the identified classes differ from the ones identified by the
# alternative clustering methods?
#
#   Rscript R/realdata_agreement.R <results-dir> <combined-stem>
#
# e.g.  Rscript R/realdata_agreement.R .../results_rd .../rd
#
# ===========================================================================
# REWRITTEN 2026-09-09. THE K-CENTRES VARIANTS ARE NO LONGER REFITTED HERE.
#
# The previous version refitted MKCFC and MRKCFC inside this script so that
# "every row of the matrix comes from the same code path". The array now fits
# all three variants at every K with ten restarts each, so refitting would give
# a DIFFERENT answer from the one the paper reports -- one start instead of ten,
# and no best-of-N. Reading the reported partitions is both cheaper and more
# honest.
#
# That removes every pre-data defect this file carried:
#
#   - `rd_optns(F$nK)`            deleted; wrong for BOTH branches
#   - `LOO = TRUE`                the comparability break with the simulation
#   - `init.method = "kmeans"`    ABORTS on this cohort -- measured
#   - `length(prep$Lt)`           THE SILENT P-LIST BUG, and here it was the
#                                 subject count AND the validity check, so
#                                 every correct comparator would have been
#                                 rejected for "returning 1405 labels for 3
#                                 subjects" -- a plausible-looking wrong answer
#                                 rather than an error
#
# It also means this script loads no data at all: the scores and the reported
# partitions are both in the combined object.
#
# WHAT AGREEMENT DOES AND DOES NOT SHOW. A high cross-method aRI means the
# methods find the same structure; it does NOT mean the structure is real, since
# methods sharing a stage-one MFPCA agree partly BECAUSE they share it. A low
# aRI means they disagree and says nothing about which is right. Neither
# direction is evidence of correctness -- that is what D3's external survival
# criterion is for, and the manuscript says so.
#
# THE COMPARATOR SET IS THE SIMULATION'S, deliberately. Using a different set
# here would invite the question of why. Methods that cannot be fitted on this
# cohort are reported as not fitted, with the reason, rather than omitted.
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L)
  stop("usage: realdata_agreement.R <results-dir> <combined-stem>", call. = FALSE)
inp  <- args[1]
stem <- args[2]

source(file.path("R", "metrics.R"))      # ari()
# ATTACHED, as methods.R does: mclust::Mclust() looks up mclustBIC() on the
# search path and fails with "could not find function" when the package is
# only installed. It cost the mclust comparator on all three arms (2026-09-12).
suppressMessages(library(mclust))

cf <- paste0(stem, "_combined.rds")
if (!file.exists(cf))
  stop("run realdata_combine.R first; no ", cf, call. = FALSE)
C <- readRDS(cf)
if (is.na(C$K_star)) stop("the combine did not determine K.", call. = FALSE)
K      <- as.integer(C$K_star)
n      <- length(C$ids)
# THE COMPARATORS' SCORE SPACE. In the simulation, mclust, k-means and teigen
# run on the NON-ROBUST global MFPCA scores (methods.R: M4, M5, M8) and only
# M6 is k-means on the robust scores. Until 2026-09-12 this script handed all
# of them the robust scores, whose least-squares score step explodes for the
# cohort's four-visit subjects (the four-visit score problem, Discussion): on every
# arm k-means isolated ONE subject and teigen found no converging model. The
# combine now saves both spaces; an older combine falls back with a note.
scores_rb <- C$scores
scores_nr <- C$scores_nonrobust
if (is.null(scores_nr)) {
  cat("NOTE: this combine predates scores_nonrobust; the comparators run on the\n",
      "      ROBUST scores, which is NOT the simulation's definition of them.\n", sep = "")
  scores_nr <- scores_rb
}
cat(sprintf("Study D2: K = %d, n = %d subjects, bwK = %s\n", K, n, format(C$bwK)))

SEED <- 20260904L
have <- function(p) requireNamespace(p, quietly = TRUE)

# The three K-centres variants come from the array: the best-of-N partition at
# K*, which is what the paper reports. `C$labels` carries one column per variant
# that had a usable fit at K*.
kc <- setdiff(names(C$labels), "ID")
methods <- setNames(lapply(kc, function(v) function() as.integer(C$labels[[v]])), kc)

methods <- c(methods, list(
  `mclust`  = function() {
    if (!have("mclust")) stop("mclust not installed")
    as.integer(mclust::Mclust(scores_nr, G = K, verbose = FALSE)$classification) },
  `k-means` = function() {
    set.seed(SEED); as.integer(stats::kmeans(scores_nr, K, nstart = 25)$cluster) },
  `k-means/rob` = function() {                 # the simulation's M6
    set.seed(SEED); as.integer(stats::kmeans(scores_rb, K, nstart = 25)$cluster) },
  `teigen`  = function() {
    if (!have("teigen")) stop("teigen not installed")
    set.seed(SEED)
    as.integer(teigen::teigen(scores_nr, Gs = K, scale = FALSE,
                              verbose = FALSE)$classification) },
  `lcmm`    = function() {
    # EXCLUDED BY DECISION, 2026-09-09. Kept in the list so the exclusion is
    # visible in every run's output rather than silent. multlcmm on this cohort
    # is a modelling exercise, not a wrapper call -- it needs a link function
    # and a fixed/random specification per marker -- and params_global.csv
    # records that it does not converge on the simulation design at 100, 300
    # or 1000 iterations. No claim is made about it for the real-data
    # agreement; it remains a comparator in the simulation study.
    stop("excluded by decision 2026-09-09: needs a considered per-marker ",
         "specification; not used for the real-data comparison")
  }))

labs <- list(); notes <- character(0)
for (nm in names(methods)) {
  r <- tryCatch(methods[[nm]](), error = function(e) e)
  if (inherits(r, "condition")) {
    notes[nm] <- conditionMessage(r)
    cat(sprintf("  %-12s NOT FITTED: %s\n", nm, conditionMessage(r)))
  } else if (length(r) != n) {
    notes[nm] <- sprintf("returned %d labels for %d subjects", length(r), n)
    cat(sprintf("  %-12s NOT USED: %s\n", nm, notes[nm]))
  } else {
    labs[[nm]] <- as.integer(r)
    cat(sprintf("  %-12s sizes %s\n", nm, paste(table(r), collapse = "/")))
  }
}
if (length(labs) < 2) stop("fewer than two methods produced a partition",
                           call. = FALSE)

# ------------------------------------------------ the agreement matrix -------
nms <- names(labs)
A <- matrix(NA_real_, length(nms), length(nms), dimnames = list(nms, nms))
for (i in seq_along(nms)) for (j in seq_along(nms))
  A[i, j] <- if (i == j) 1 else ari(labs[[i]], labs[[j]])

cat("\n=== cross-method agreement, adjusted Rand index ===\n")
cat(sprintf("%-12s %s\n", "", paste(sprintf("%11s", substr(nms, 1, 11)), collapse = "")))
for (i in seq_along(nms))
  cat(sprintf("%-12s %s\n", substr(nms[i], 1, 12),
              paste(sprintf("%11.3f", A[i, ]), collapse = "")))

off <- A[upper.tri(A)]
cat(sprintf("\noff-diagonal: median %.3f, min %.3f, max %.3f\n",
            stats::median(off), min(off), max(off)))

# ------------------------------- AGAINST WITHIN-METHOD STABILITY -------------
#
# THIS IS THE COMPARISON THAT MAKES THE MATRIX READABLE, and without it the
# matrix can be badly over-read. A method that agrees with ITSELF at 0.60 across
# its own random starts cannot meaningfully be said to disagree with another
# method at 0.55. The cross-method number only carries information above the
# method's own reproducibility floor.
sf <- paste0(stem, "_stability.rds")
if (file.exists(sf)) {
  S <- readRDS(sf)
  cat("\n=== cross-method agreement against within-method stability ===\n")
  cat(sprintf("  %-12s %14s %16s   %s\n", "variant", "own restarts",
              "best other method", "reading"))
  for (v in names(S$per_variant)) {
    if (!v %in% nms) next
    own <- stats::median(S$per_variant[[v]]$pairwise_ari)
    oth <- max(A[v, setdiff(nms, v)], na.rm = TRUE)
    cat(sprintf("  %-12s %14.3f %16.3f   %s\n", v, own, oth,
                if (oth >= own) "cross-method agreement EXCEEDS self-agreement"
                else "self-agreement is higher, as expected"))
  }
  cat("\n  Where cross-method agreement exceeds a method's agreement with its own\n",
      "  restarts, the partition is not distinguishing the methods -- it is\n",
      "  being driven by the shared stage-one MFPCA and the score geometry.\n",
      "  Report that rather than the matrix alone.\n", sep = "")
} else {
  cat("\nRun realdata_stability.R first to read this matrix against each\n",
      "method's own reproducibility floor; a cross-method aRI is not\n",
      "interpretable without it.\n", sep = "")
}

cat("\nRead the matrix as DESCRIPTIVE. Methods sharing a stage-one MFPCA agree\n")
cat("partly BECAUSE they share it, and no value here is evidence that any\n")
cat("partition is correct. The external criterion is Study D3.\n")

if (length(notes))
  cat(sprintf("\nnot fitted: %s\n", paste(names(notes), collapse = ", ")))

out <- file.path(dirname(stem), paste0(basename(stem), "_agreement.rds"))
saveRDS(list(K = K, n = n, bwK = C$bwK, labels = labs, agreement = A,
             notes = notes, when = Sys.time()), out)
utils::write.csv(as.data.frame(A), paste0(stem, "_agreement.csv"))
cat(sprintf("\nwrote %s and %s_agreement.csv\n", out, stem))
