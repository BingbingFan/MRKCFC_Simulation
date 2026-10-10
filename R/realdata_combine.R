# Study D: read the 90 (variant, K, restart) fits, select the restart, select K,
# and report.
#
#   Rscript R/realdata_combine.R <dir-of-rd_task_*.rds> [outstem] [n_restart] [kmax]
#
# ONE ARM PER CALL. The primary run (warped time, bwK GCV, raw values), the
# sensitivity arm (fixed bwK 0.20) and the log arm (every marker log-transformed
# before normalisation, 2026-09-11) write to separate directories and are
# combined separately; pooling them would compare partitions fitted under
# different smoothing or on different scales. The design check below refuses
# it anyway, on the `time`, `bwK` and `log` fields.
#
# The outstem is a PATH, like combine_estim.R's and unlike combine_c.R's.
#
# THE SELECTION AND THE AGREEMENT LIVE HERE, not in the task driver: no task
# holds more than one fit, and the criteria compare partitions of a fixed score
# space across K. The task driver writes only.
#
# ===========================================================================
# REWRITTEN 2026-09-09 for the re-designed grid. The previous version read 15
# files -- 3 variants x K in 2..6, ONE fit each, no restarts -- and knew nothing
# about `status`. It would have pooled an earlier run's files with the new
# ones and selected among partitions that are, in fourteen recorded cases, an
# UNTOUCHED RANDOM START.
#
# TWO SELECTIONS, TWO CRITERIA, AND THEY MUST NOT BE SWAPPED:
#
#   restart, within a (variant, K) cell -> the OBJECTIVE
#       Total squared reconstruction error under the partition. It uses no
#       outcome and no truth, so it cannot be circular. Valid here because
#       every restart in a cell fits the same K to the same data.
#
#   K, across cells -> the SILHOUETTE
#       The objective FALLS MONOTONICALLY IN K -- more classes always
#       reconstruct better -- so using it to choose K would always return the
#       largest K offered. That is the same defect that made the uncorrected
#       CH_mv monotone (see R/kselect.R and mvclass_enumeration.R), and it is
#       worth naming twice because it looks like a reasonable thing to do.
# ===========================================================================

source(file.path("R", "study_c.R"))       # load_config, build_optns via methods.R
source(file.path("R", "realdata_load.R"))

args <- commandArgs(trailingOnly = TRUE)
inp  <- if (length(args) >= 1) args[1] else "results_rd"
stem <- if (length(args) >= 2) args[2] else file.path(inp, "rd")
n_restart <- if (length(args) >= 3) as.integer(args[3]) else 10L
kmax      <- if (length(args) >= 4) as.integer(args[4]) else 4L

VARIANT_ORDER <- c("MKCFC", "MRKCFC", "MRKCFC-R2")
K_GRID        <- 2:kmax
N_EXPECT      <- length(VARIANT_ORDER) * length(K_GRID) * n_restart

# Statuses that carry a real partition. `init-only` never does -- its labels are
# the random start and the driver writes them as NULL -- and `incomplete` is a
# fit the wall interrupted, which is to be RESUMED by re-submitting the array,
# not reported. `maxiter` IS usable: the simulation runs the same conv.prop = 0
# and maxIter and reports its non-convergence rate as a result, so excluding
# those here would have the two sections reporting different things under the
# same name.
STATUS_USABLE  <- c("converged", "maxiter", "degenerate", "error")
STATUS_STRICT  <- c("converged", "maxiter")   # the sensitivity set, see below
STATUS_NEVER   <- c("init-only", "error-early", "incomplete")

fs <- list.files(inp, pattern = "^rd_task_[0-9]+[.]rds$", full.names = TRUE)
if (!length(fs)) stop("no rd_task_*.rds under ", inp, call. = FALSE)
R <- lapply(fs, readRDS)

# ------------------------------------------------ DESIGN CONSISTENCY -------
# THE FILENAMES ARE NOT UNIQUE ACROSS DESIGNS. rd_task_002.rds means task 2 of
# whatever grid wrote it, and an earlier 300-task run used the same names in
# the same directory. Anything without a `status` field predates this driver.
# Refusing here rather than warning: a pooled combine produces a plausible
# number, and a plausible wrong number is the failure mode this project keeps
# hitting.
no_status <- vapply(R, function(z) !is.character(z$status), TRUE)
if (any(no_status)) {
  stop(sum(no_status), " of ", length(R), " files carry no `status` field and ",
       "were written by an earlier driver.\n  First: ",
       basename(fs[which(no_status)[1]]),
       "\n  Move them aside; do not combine two designs.", call. = FALSE)
}
one <- function(f, subset = rep(TRUE, length(R))) {
  # A file without the field says "<missing>" rather than tripping vapply, so
  # the mixed-design message names the field instead of a type error.
  v <- unique(vapply(R[subset], function(z)
    if (is.null(z[[f]])) "<missing>" else as.character(z[[f]]), ""))
  if (length(v) != 1L)
    stop("files disagree on `", f, "`: ", paste(v, collapse = ", "),
         "\n  Two designs are mixed in ", inp, ". Separate them.", call. = FALSE)
  v
}
win    <- as.integer(one("window"))
nK     <- as.integer(one("nK"))
uselog <- as.logical(one("log"))
tm     <- one("time")            # "warped" or "rescaled"; every file must agree
# bw binds on the ROBUST branch only, so it is checked among those files alone;
# a non-robust fit carrying a different bw is not evidence of a mixed design.
is_rob <- vapply(R, function(z) isTRUE(z$bwK_binds), TRUE)
bwK    <- if (any(is_rob)) one("bwK", is_rob) else NA_character_   # "GCV" or a number, as text

ids <- R[[1]]$ids
cat(sprintf("read %d of %d expected task files\n", length(R), N_EXPECT))
cat(sprintf("cohort: %d-month, %d subjects, nK = %d, %s scale, %s time, bwK = %s (robust branch)\n",
            win, length(ids), nK, if (uselog) "log" else "raw", tm,
            if (is.na(bwK)) "n/a" else bwK))
if (length(R) < N_EXPECT)
  cat(sprintf("  INCOMPLETE: %d file(s) missing. K selection over a partial grid\n",
              N_EXPECT - length(R)),
      "  is not the same selection; the per-cell denominators below say where.\n",
      sep = "")

# ------------------------------------------------------- the fit table -----
# One row per fit. A named function rather than an inline lapply body so the
# classification can be exercised directly; it is the part of this file that a
# wrong answer would flow through into the paper.
fit_row <- function(z, fname) {
  lab <- z$labels
  # A partition that lost a class has fewer than K distinct labels. select_K()
  # derives K as length(unique(g)), so such a fit would be silently compared as
  # though it were a smaller-K fit. Caught here instead.
  kobs <- if (is.null(lab)) NA_integer_ else length(unique(lab))
  intact <- !is.null(lab) && isTRUE(kobs == as.integer(z$K))
  data.frame(
    file = fname, variant = z$variant, K = as.integer(z$K),
    restart = as.integer(z$restart), status = z$status,
    usable = z$status %in% STATUS_USABLE && intact,
    strict = z$status %in% STATUS_STRICT && intact,
    K_obs = kobs, collapsed = !is.null(lab) && isTRUE(kobs < as.integer(z$K)),
    converged = isTRUE(z$converged), iter = as.integer(z$iter),
    chunks = as.integer(z$chunks), objective = as.numeric(z$objective),
    # `objective` above is the BEST iterate's; these say how far the loop
    # wandered to reach it. The spread is reportable: on the probe the
    # objective moved 13.8% with no trend, so "the fit" is the best of a
    # trajectory, not a converged point. See write_ckpt() in the driver.
    best_iter = if (is.null(z$best_iter)) NA_integer_ else as.integer(z$best_iter),
    final_objective = if (is.null(z$final_objective)) NA_real_
                      else as.numeric(z$final_objective),
    spread_pct = if (is.null(z$obj_trace) || nrow(z$obj_trace) < 2L) NA_real_
                 else 100 * (max(z$obj_trace$objective) -
                             min(z$obj_trace$objective)) /
                      min(z$obj_trace$objective),
    hours = as.numeric(z$seconds) / 3600, ncore = as.integer(z$ncore),
    sizes = if (is.null(lab)) NA_character_ else paste(table(lab), collapse = "/"),
    error = if (is.null(z$error)) NA_character_ else z$error,
    stringsAsFactors = FALSE)
}
info <- do.call(rbind, lapply(seq_along(R),
                              function(i) fit_row(R[[i]], basename(fs[i]))))
info <- info[order(match(info$variant, VARIANT_ORDER), info$K, info$restart), ]
rownames(info) <- NULL

cat("\n=== outcomes by cell (denominators shown; see the survivorship rule) ===\n")
cat(sprintf("%-10s %3s %9s %9s %8s %8s   %s\n",
            "variant", "K", "files", "usable", "converg", "collapsed", "statuses"))
for (v in VARIANT_ORDER) for (k in K_GRID) {
  d <- info[info$variant == v & info$K == k, ]
  if (!nrow(d)) { cat(sprintf("%-10s %3d %9s\n", v, k, "-- ABSENT --")); next }
  tb <- table(d$status)
  cat(sprintf("%-10s %3d %4d/%-4d %4d/%-4d %4d/%-4d %8d   %s\n",
              v, k, nrow(d), n_restart, sum(d$usable), nrow(d),
              sum(d$converged), nrow(d), sum(d$collapsed),
              paste(sprintf("%s:%d", names(tb), as.integer(tb)), collapse = " ")))
}
if (any(info$status %in% STATUS_NEVER))
  cat(sprintf("\n%d fit(s) are init-only or incomplete and are excluded by rule.\n",
              sum(info$status %in% STATUS_NEVER)))
if (any(info$collapsed))
  cat(sprintf("%d fit(s) lost a class (K_obs < K) and are excluded: a partition\n",
              sum(info$collapsed)),
      "  with fewer classes than it was fitted at is not a fit at that K.\n", sep = "")

# ------------------------------------------------- best-of-N per cell ------
# BY THE OBJECTIVE, WHICH IS VALID WITHIN A CELL AND NOT ACROSS THEM. Every
# restart here fits the same K to the same data, so the comparison is of
# partitions and nothing else.
pick <- function(mask, tag) {
  out <- list()
  for (v in VARIANT_ORDER) for (k in K_GRID) {
    d <- which(info$variant == v & info$K == k & mask & is.finite(info$objective))
    if (!length(d)) next
    b <- d[which.min(info$objective[d])]
    out[[length(out) + 1L]] <- data.frame(
      variant = v, K = k, set = tag, file = info$file[b],
      restart = info$restart[b], status = info$status[b],
      objective = info$objective[b], n_ok = length(d),
      n_attempt = sum(info$variant == v & info$K == k),
      stringsAsFactors = FALSE)
  }
  if (!length(out)) NULL else do.call(rbind, out)
}
best        <- pick(info$usable, "usable")
best_strict <- pick(info$strict, "strict")
if (is.null(best)) stop("no usable fit in any cell; nothing to combine.", call. = FALSE)

cat("\n=== best restart per cell, by the objective ===\n")
cat(sprintf("%-10s %3s %8s %9s %14s   %s\n",
            "variant", "K", "restart", "of", "objective", "status"))
for (i in seq_len(nrow(best)))
  cat(sprintf("%-10s %3d %8d %4d/%-4d %14.1f   %s\n",
              best$variant[i], best$K[i], best$restart[i],
              best$n_ok[i], best$n_attempt[i], best$objective[i], best$status[i]))

# HOW FAR THE LOOP WANDERED, which decides how much the selection above means.
# If the within-fit spread is comparable to the between-restart spread, the
# best-of-N is choosing among trajectory noise and should be reported as such.
sp <- info$spread_pct[info$usable & is.finite(info$spread_pct)]
if (length(sp)) {
  cat(sprintf("\nwithin-fit objective spread across checkpoints: median %.1f%%, max %.1f%%\n",
              stats::median(sp), max(sp)))
  for (v in VARIANT_ORDER) for (k in K_GRID) {
    d <- info[info$variant == v & info$K == k & info$usable, ]
    if (nrow(d) < 2L) next
    btw <- 100 * (max(d$objective) - min(d$objective)) / min(d$objective)
    wth <- stats::median(d$spread_pct[is.finite(d$spread_pct)])
    cat(sprintf("  %-10s K=%d  between restarts %6.1f%%   within fits %6.1f%%%s\n",
                v, k, btw, wth,
                if (is.finite(wth) && wth >= btw)
                  "   <- selection is inside the noise" else ""))
  }
}
nconv <- sum(info$usable & !info$converged)
if (nconv)
  cat(sprintf("\n%d of %d usable fits did NOT converge; `labels` is the best iterate\n",
              nconv, sum(info$usable)),
      "  seen, not a fixed point. Report the rate, as Study A does.\n", sep = "")

# SENSITIVITY, NOT DECORATION. `degenerate` and `error` fits carry a real
# partition but stopped early, so whether they are allowed to win a cell is a
# judgement. Reporting both makes it auditable rather than silent.
if (!is.null(best_strict)) {
  m <- merge(best[, c("variant", "K", "restart")],
             best_strict[, c("variant", "K", "restart")],
             by = c("variant", "K"), suffixes = c("_usable", "_strict"))
  diff <- m[m$restart_usable != m$restart_strict, ]
  cat(sprintf("\nrestricting to converged/maxiter changes the winner in %d of %d cell(s)",
              nrow(diff), nrow(m)))
  cat(if (nrow(diff)) ":\n" else ".\n")
  if (nrow(diff)) print(diff, row.names = FALSE)
}

# --------------------------------------------------------- selecting K -----
#
# ONE COMMON SCORE SPACE, built once over all subjects, so the criteria compare
# partitions of a FIXED geometry. Rebuilding it per K would move the space and
# the partition together and the comparison would no longer be about K.
#
# THE SEARCH CANNOT INCLUDE K = 1 and that is a limitation, not a setting.
# Neither criterion is defined at a single class, so this answers "how many
# classes, given at least two" and CANNOT report the absence of clustering. Say
# so in the Discussion, and rest the interpretation on the survival benchmark
# rather than on internal fit.
#
# THE SILHOUETTE GOVERNS. Study C3 measured both: on the hard simulation cell
# Calinski-Harabasz fails completely paired with one comparator, recovering the
# true K in NONE of fifty replicates while the silhouette recovers it in
# forty-nine on the same fits. CH is reported beside it, not as an equal.
dat  <- rd_load(window = win)
prep <- rd_prepare(dat, log = uselog, verbose = FALSE, time = tm)
gl <- MRKCFC::MFPCA_k(1L, prep$Ly, prep$Lt, rep(1L, prep$n),
                     rd_opt(nK, TRUE, bwK = if (is.na(bwK) || identical(bwK, "GCV")) NULL
                                            else as.numeric(bwK)),
                     TRUE, mfpca_scatter = "sign")
scores <- gl$multi_score
# THE SAME COHORT THROUGH THE NON-ROBUST BRANCH, for two uses. The comparators
# in realdata_agreement.R are defined on these scores in the simulation
# (methods.R M4, M5, M8), and the silhouette is printed on both spaces below,
# because the robust branch estimates scores by least squares on each
# subject's own observations and that solve explodes for the cohort's
# four-visit subjects (the four-visit score problem, 2026-09-12). K
# SELECTION ITSELF STAYS ON THE ROBUST SPACE, as Study C3 had it; the second
# column is there so a reader can see how much the space moves the number.
gl_nr <- tryCatch(MRKCFC::MFPCA_k(1L, prep$Ly, prep$Lt, rep(1L, prep$n),
                                 rd_opt(nK, FALSE), FALSE, mfpca_scatter = "moment"),
                  error = function(e) e)
scores_nr <- if (inherits(gl_nr, "error")) {
  cat(sprintf("non-robust global MFPCA FAILED (%s); no second silhouette, and the\n  comparators will fall back to the robust scores\n",
              conditionMessage(gl_nr)))
  NULL
} else gl_nr$multi_score

by_file <- setNames(R, basename(fs))
sel_all <- list()
cat("\n=== K selection, per variant (silhouette governs) ===\n")
for (v in VARIANT_ORDER) {
  b <- best[best$variant == v, ]
  if (nrow(b) < 2L) {
    cat(sprintf("\n  %s: %d cell(s) with a usable fit -- K cannot be selected.\n",
                v, nrow(b)))
    next
  }
  parts <- lapply(b$file, function(f) by_file[[f]]$labels)
  names(parts) <- as.character(b$K)
  s <- MRKCFC:::select_K(scores, parts, distance = "euclidean", robust = TRUE)
  s$silhouette_nonrobust <- if (is.null(scores_nr)) NA_real_ else
    MRKCFC:::select_K(scores_nr, parts, distance = "euclidean", robust = TRUE)$silhouette
  s$variant <- v; s$K <- b$K; s$restart <- b$restart; s$n_ok <- b$n_ok
  s$n_attempt <- b$n_attempt
  sel_all[[v]] <- s
  ks <- s$K[which.max(s$silhouette)]; kc <- s$K[which.max(s$CH)]
  cat(sprintf("\n  %s  (K grid covered: %s of %s)\n", v,
              paste(b$K, collapse = ","), paste(K_GRID, collapse = ",")))
  cat(sprintf("    %-4s %12s %12s %14s %10s\n", "K", "CH", "silhouette",
              "sil(non-robust)", "best of"))
  for (i in seq_len(nrow(s)))
    cat(sprintf("    %-4d %12.3f %12.3f %14.3f %6d/%-4d\n", s$K[i], s$CH[i],
                s$silhouette[i], s$silhouette_nonrobust[i], s$n_ok[i], s$n_attempt[i]))
  cat(sprintf("    silhouette -> K = %d ; CH -> K = %d%s\n", ks, kc,
              if (ks != kc) "   DISAGREE" else ""))
  if (length(b$K) < length(K_GRID))
    cat("    NOTE: selected over a PARTIAL K grid; the missing cells had no\n",
        "          usable fit, which is itself a result -- report it.\n", sep = "")
}
sel <- if (length(sel_all)) do.call(rbind, sel_all) else NULL
if (!is.null(sel)) rownames(sel) <- NULL

K_star <- NA_integer_
if ("MRKCFC-R2" %in% names(sel_all)) {
  s <- sel_all[["MRKCFC-R2"]]
  K_star <- s$K[which.max(s$silhouette)]
  cat(sprintf("\nPRIMARY K = %d, from MRKCFC-R2 on the silhouette.\n", K_star))
} else {
  cat("\nPRIMARY K UNDETERMINED: MRKCFC-R2 has no selectable K.\n")
}
kk <- vapply(sel_all, function(s) s$K[which.max(s$silhouette)], integer(1))
if (length(kk) > 1L && length(unique(kk)) > 1L)
  cat("  NOTE: the variants do NOT agree on K (",
      paste(sprintf("%s:%d", names(kk), kk), collapse = ", "),
      "). Report that rather than presenting one as the answer.\n", sep = "")

# ------------------------------------------------------ agreement ----------
#
# DESCRIPTIVE ONLY. A high value means the variants find the same structure; it
# is NOT evidence that the structure is right, since all three share a stage-one
# MFPCA and agree partly BECAUSE they share it. The external criterion is the
# survival benchmark in realdata_survival.R.
lab <- data.frame(ID = ids, stringsAsFactors = FALSE)
if (!is.na(K_star)) {
  b <- best[best$K == K_star, ]
  if (nrow(b) >= 2L) {
    P <- lapply(b$file, function(f) by_file[[f]]$labels)
    nm <- b$variant
    cat(sprintf("\n=== agreement at K = %d (adjusted Rand index) ===\n", K_star))
    cat(sprintf("  %-12s %s\n", "", paste(sprintf("%12s", nm), collapse = "")))
    for (i in seq_along(P))
      cat(sprintf("  %-12s %s\n", nm[i], paste(vapply(seq_along(P), function(j)
        sprintf("%12.3f", mclust::adjustedRandIndex(P[[i]], P[[j]])),
        character(1)), collapse = "")))
  }
  for (i in seq_len(nrow(b))) lab[[b$variant[i]]] <- by_file[[b$file[i]]]$labels
}

utils::write.csv(info, paste0(stem, "_restarts.csv"), row.names = FALSE)
utils::write.csv(rbind(best, best_strict), paste0(stem, "_fits.csv"),
                 row.names = FALSE)
if (!is.null(sel)) utils::write.csv(sel, paste0(stem, "_selection.csv"),
                                    row.names = FALSE)
utils::write.csv(lab, paste0(stem, "_labels.csv"), row.names = FALSE)
saveRDS(list(info = info, best = best, best_strict = best_strict,
             selection = sel, K_star = K_star, labels = lab, scores = scores,
             scores_nonrobust = scores_nr,
             ids = ids, window = win, nK = nK, log = uselog, bwK = bwK,
             time = tm, time_map = prep$time_map,
             n_restart = n_restart, when = Sys.time()),
        paste0(stem, "_combined.rds"))
cat(sprintf("\nwrote %s_restarts.csv, %s_fits.csv, %s_selection.csv, %s_labels.csv and %s_combined.rds\n",
            stem, stem, stem, stem, stem))
