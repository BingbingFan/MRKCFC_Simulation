# Combine and summarise a Study C sub-study.
#
#   Rscript R/combine_c.R <resultsdir> [outstem]
#
# WHY THIS IS SEPARATE FROM combine.R, which is not a duplication
#
# combine.R summarises one fixed shape: scenario x method, with CCR, aRI, MPE
# and the rest. Study A and Study B both produce that shape, which is why they
# share it. C's four sub-studies do not share it with each other, let alone with
# A:
#
#   C2  adds M and basis; the question is decay away from M_true
#   C3  has NO CCR or aRI at all -- it reports criterion curves over K
#   C4  adds conv.prop and a failure mode; the question is a RATE
#   C5  adds n; the question is wall clock, not accuracy
#
# A single summary over those would be meaningless, so each gets its own. The
# sub-study is detected from the file names rather than passed as an argument,
# so it cannot be mislabelled.

args   <- commandArgs(trailingOnly = TRUE)
resdir <- if (length(args) >= 1) args[1] else "results_c2"
stem   <- if (length(args) >= 2) args[2] else NULL

files <- list.files(resdir, pattern = "^studyC[2-5]_task_.*[.]rds$", full.names = TRUE)
if (!length(files)) stop("no studyC*_task_*.rds files in ", resdir, call. = FALSE)
sub <- unique(sub("^studyC([2-5])_task_.*$", "C\\1", basename(files)))
if (length(sub) != 1L)
  stop("directory mixes sub-studies: ", paste(sub, collapse = ", "),
       ". Keep each in its own results directory.", call. = FALSE)
if (is.null(stem)) stem <- file.path(resdir, paste0("study", sub))

all <- do.call(rbind, lapply(files, readRDS))
cat(sprintf("%s: combining %d task files, %d rows\n", sub, length(files), nrow(all)))

med_iqr <- function(x) sprintf("%.3f (%.3f)", stats::median(x, na.rm = TRUE),
                               stats::IQR(x, na.rm = TRUE))

# --- completeness and failures, before anything is summarised ----------------
if ("error" %in% names(all)) {
  e <- all[!is.na(all$error), ]
  cat(sprintf("\nfailed fits: %d of %d\n", nrow(e), nrow(all)))
  if (nrow(e)) print(table(e$method, substr(e$error, 1, 45)))
}
if ("converged" %in% names(all)) {
  cat("\nnon-convergence by method:\n")
  print(round(tapply(all$converged, all$method,
                     function(z) 1 - mean(z, na.rm = TRUE)), 4))
}

grp <- function(d, by, metric, fun = med_iqr) {
  s <- split(d, d[by], drop = TRUE)
  out <- do.call(rbind, lapply(s, function(x) {
    r <- x[1, by, drop = FALSE]
    for (m in metric) r[[m]] <- fun(x[[m]])
    r
  }))
  rownames(out) <- NULL
  # unname() is load-bearing: `by` contains a column called "method", and
  # do.call(order, <named list>) would pass it as order()'s own `method`
  # argument, failing with "'arg' must be of length 1".
  out[do.call(order, unname(as.list(out[by]))), ]
}

if (sub == "C2") {
  # M_true = 2, so the M axis brackets the truth. What C2 reports is how fast
  # performance decays away from it -- the question a reader has when choosing
  # M on data where the truth is unknown.
  cat("\n=== median aRI (IQR) by configuration, M and basis ===\n")
  summ <- grp(all, c("config_id", "method", "M", "basis"), c("aRI", "CCR"))
  print(summ, row.names = FALSE)
  cat("\n=== the basis ablation: class-specific minus class-common, median aRI ===\n")
  cat("The design gives variables 1 and 2 class-specific eigenfunctions, so a\n")
  cat("common basis is misspecified BY CONSTRUCTION. The direction is expected;\n")
  cat("the SIZE is the result -- it is what class-specific estimation buys.\n\n")

  # --- C2 SUMMARY FILE, added 2026-09-04 -------------------------------------
  #
  # WRITTEN, not just printed, because make_tables.R reads summaries and must
  # not reach past them into the raw file -- that is how this project acquired
  # two divergent copies of silhouette_score(). Before this, C2 produced only
  # `_raw.csv` and its table could not be generated at all.
  #
  # IT CARRIES A RATE AND A CONDITIONAL MEDIAN, NOT ONLY A POOLED MEDIAN.
  # C_HARD is cleanly bimodal: over 50 replicates MRKCFC-R2 puts 17 at about
  # zero and 26 above 0.75, with an almost empty middle. The pooled median of
  # 0.774 sits inside the UPPER mode only because 58 percent of replicates land
  # there, and a few replicates fewer would report 0.00 -- so a pooled median
  # comparison against another method whose median falls in the LOWER mode
  # (MKCFC, 0.106) reads as a 7x accuracy effect when it is nothing of the kind.
  # The honest statistic is how OFTEN each method succeeds, and how well it does
  # when it does. This is the same treatment arm 6 gets everywhere else.
  SUCCESS <- 0.5                      # the empty middle of the split
  c2s <- do.call(rbind, lapply(
    split(all, list(all$config_id, all$method_id, all$M, all$basis), drop = TRUE),
    function(d) {
      v <- d$aRI[is.na(d$error)]
      s <- v > SUCCESS
      data.frame(
        config_id = d$config_id[1], method_id = d$method_id[1],
        method = d$method[1], M = d$M[1], basis = d$basis[1],
        n_fits = length(v),
        aRI_median = stats::median(v, na.rm = TRUE),
        aRI_iqr = stats::IQR(v, na.rm = TRUE),
        success_rate = mean(s, na.rm = TRUE),
        aRI_med_success = if (any(s)) stats::median(v[s]) else NA_real_,
        aRI_med_failure = if (any(!s)) stats::median(v[!s]) else NA_real_,
        CCR_median = stats::median(d$CCR[is.na(d$error)], na.rm = TRUE),
        stringsAsFactors = FALSE)
    }))
  rownames(c2s) <- NULL
  utils::write.csv(c2s, paste0(stem, "_summary.csv"), row.names = FALSE)
  cat(sprintf("wrote %s_summary.csv (%d rows)\n\n", stem, nrow(c2s)))
  for (cid in sort(unique(all$config_id))) for (m in sort(unique(all$method))) {
    v <- sapply(c("specific", "common"), function(b)
      stats::median(all$aRI[all$config_id == cid & all$method == m & all$basis == b],
                    na.rm = TRUE))
    cat(sprintf("  %-7s %-18s specific %6.3f  common %6.3f  gain %+6.3f\n",
                cid, m, v[1], v[2], v[1] - v[2]))
  }

} else if (sub == "C3") {
  # Grid is K = 2..6. Neither criterion is defined at K = 1, so these cannot
  # report the absence of clustering -- state that limitation in the paper.
  crits <- c("CH_euclid", "sil_euclid", "CH_mahal", "sil_mahal")
  Ktrue <- unique(all$K_true)[1]
  cat(sprintf("\n=== K selected by argmax over K = 2..6 (K_true = %d) ===\n", Ktrue))
  # `n_rep` IS NOT DECORATION. A recovery rate is computed over the replicates on
  # which that method returned a partition at some K, and a method that fails on
  # the hard cell contributes far fewer. The t-mixture recovers K = 2 in 0.80 of
  # its C_HARD replicates -- over FIVE of fifty, while every other method-config
  # cell uses all fifty. Printed without the count, that 0.80 reads like the
  # others. Same survivorship trap as the Study A tables and the C1 summary; it
  # has now appeared three times, so do not drop this column.
  n_all_reps <- length(unique(all$rep))
  cat(sprintf("%-8s %-20s %10s %9s %6s  %s\n", "config", "method", "criterion",
              "recovery", "n_rep", "distribution of selected K"))
  for (cid in sort(unique(all$config_id))) for (m in sort(unique(all$method)))
    for (cr in crits) {
      d <- all[all$config_id == cid & all$method == m, ]
      sp <- split(d, d$rep, drop = TRUE)
      kk <- vapply(sp, function(x) if (all(is.na(x[[cr]]))) NA_integer_
                   else x$K[which.max(x[[cr]])], integer(1))
      if (!length(kk)) next
      tb <- table(kk)
      nk <- sum(!is.na(kk))
      cat(sprintf("%-8s %-20s %10s %9.2f %6d%s  %s\n", cid, m, cr,
                  mean(kk == Ktrue, na.rm = TRUE), nk,
                  if (nk < n_all_reps) "*" else " ",
                  paste(sprintf("K%s:%d", names(tb), as.integer(tb)), collapse = " ")))
    }
  cat(sprintf("  * fewer than the %d replicates in the design: the rate is\n",
              n_all_reps),
      "    conditional on that method returning a partition at all.\n", sep = "")
  cat("\n=== do the two DISTANCES agree? ===\n")
  cat("C is a sensitivity study: a selected K that depends on the metric is\n")
  cat("itself the finding, not a nuisance.\n")
  for (pair in list(c("CH_euclid","CH_mahal"), c("sil_euclid","sil_mahal"))) {
    sp <- split(all, list(all$config_id, all$method, all$rep), drop = TRUE)
    ag <- vapply(sp, function(x) {
      a <- if (all(is.na(x[[pair[1]]]))) NA_integer_ else x$K[which.max(x[[pair[1]]])]
      b <- if (all(is.na(x[[pair[2]]]))) NA_integer_ else x$K[which.max(x[[pair[2]]])]
      if (is.na(a) || is.na(b)) NA else a == b
    }, logical(1))
    cat(sprintf("  %-11s vs %-11s agree on %.1f%% of fits\n",
                pair[1], pair[2], 100*mean(ag, na.rm = TRUE)))
  }

  # --- C3 SUMMARY FILE, added 2026-09-06 -------------------------------------
  #
  # WRITTEN, not just printed. `make_tables.R` reads summaries and must not
  # reach past them into a raw file; recomputing the recovery rate there would
  # create a second definition of it that can drift from the one printed above.
  # That is the failure which gave this project two divergent copies of
  # `silhouette_score()`. C2 got a writer for the same reason on 2026-09-04.
  #
  # `n_rep` TRAVELS WITH THE RATE, always. The t-mixture returns a partition on
  # five of fifty C_HARD replicates, so its recovery there is conditional on
  # being fittable at all. A rate without its denominator has now misled three
  # times in this project.
  c3s <- do.call(rbind, lapply(sort(unique(all$config_id)), function(cid)
    do.call(rbind, lapply(sort(unique(all$method)), function(m)
      do.call(rbind, lapply(crits, function(cr) {
        d  <- all[all$config_id == cid & all$method == m, ]
        sp <- split(d, d$rep, drop = TRUE)
        kk <- vapply(sp, function(x) if (all(is.na(x[[cr]]))) NA_integer_
                     else x$K[which.max(x[[cr]])], integer(1))
        if (!length(kk)) return(NULL)
        tb <- table(kk)
        data.frame(config_id = cid, method = m, criterion = cr, K_true = Ktrue,
                   recovery = mean(kk == Ktrue, na.rm = TRUE),
                   n_rep = sum(!is.na(kk)), n_design = n_all_reps,
                   K_dist = paste(sprintf("K%s:%d", names(tb), as.integer(tb)),
                                  collapse = " "),
                   stringsAsFactors = FALSE)
      }))))))
  rownames(c3s) <- NULL
  utils::write.csv(c3s, paste0(stem, "_summary.csv"), row.names = FALSE)
  cat(sprintf("\nwrote %s_summary.csv (%d rows)\n", stem, nrow(c3s)))

} else if (sub == "C4") {
  # The headline is a RATE, and `failure` is what makes it actionable:
  # cap means the iteration budget bound the fit, few_curves means it broke on
  # the cluster-size rule BEFORE the convergence test and no threshold reaches
  # it. Those need opposite remedies.
  cat("\n=== convergence by threshold, and WHY fits failed ===\n")
  cat(sprintf("%-6s %-20s %8s %9s %8s %11s %11s\n", "cp", "method", "conv",
              "median aRI", "cap", "few_curves", "n"))
  for (cp in sort(unique(all$conv_prop))) for (m in sort(unique(all$method))) {
    d <- all[all$conv_prop == cp & all$method == m & !is.na(all$converged), ]
    if (!nrow(d)) next
    cat(sprintf("%-6s %-20s %7.3f %9.3f %8d %11d %11d\n", format(cp), m,
                mean(d$converged), stats::median(d$aRI, na.rm = TRUE),
                sum(d$failure == "cap", na.rm = TRUE),
                sum(d$failure == "few_curves", na.rm = TRUE), nrow(d)))
  }
  cat("\n=== by scenario: is non-convergence confined to the heavy tails? ===\n")
  tb <- tapply(all$converged, list(all$scenario_id, all$method),
               function(z) round(1 - mean(z, na.rm = TRUE), 3))
  print(tb)

} else if (sub == "C5") {
  # Median rather than mean: the distribution is strongly right-skewed because a
  # non-converging fit runs to the cap while a converging one stops early, so a
  # mean is dominated by failures and describes no typical run.
  cat("\n=== median seconds per fit, by method and n ===\n")
  tb <- tapply(all$seconds, list(all$method, all$n), stats::median, na.rm = TRUE)
  print(round(tb, 1))
  cat("\n=== scaling: seconds relative to the smallest n ===\n")
  print(round(tb / tb[, 1], 2))
  cat("\n=== non-convergence by method and n (cost and convergence are linked) ===\n")
  print(round(tapply(all$converged, list(all$method, all$n),
                     function(z) 1 - mean(z, na.rm = TRUE)), 3))
  cat("\n=== median aRI, so cost is read beside what it buys ===\n")
  print(round(tapply(all$aRI, list(all$method, all$n), stats::median, na.rm = TRUE), 3))
}

utils::write.csv(all, paste0(stem, "_raw.csv"), row.names = FALSE)
cat(sprintf("\nwrote %s_raw.csv\n", stem))
