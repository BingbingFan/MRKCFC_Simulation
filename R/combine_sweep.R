# Read the M sweep and answer the three questions it was run to settle.
#
#   Rscript R/combine_sweep.R <dir-or-csv> [outstem] [fve-csv]
#
#   Rscript R/combine_sweep.R $ROOT/results_sweep $ROOT/results_sweep/sweep
#   Rscript R/combine_sweep.R $ROOT/results_sweep $ROOT/results_sweep/sweep \
#           $ROOT/results_fve/fve_task_0001.csv
#
# Reads every sweep_task_*.rds under the directory, writes a raw CSV and a
# summary CSV, and prints five tables. Pass the FVE output as a third argument
# to get the last one, which is the one the paper needs.

args   <- commandArgs(trailingOnly = TRUE)
inp    <- if (length(args) >= 1) args[1] else "results_sweep"
stem   <- if (length(args) >= 2) args[2] else file.path(inp, "sweep")
fvecsv <- if (length(args) >= 3) args[3] else NA_character_

if (dir.exists(inp)) {
  fs <- list.files(inp, pattern = "^sweep_task_[0-9]+[.]rds$", full.names = TRUE)
  if (!length(fs)) stop("no sweep_task_*.rds under ", inp)
  d <- do.call(rbind, lapply(fs, readRDS))
  cat(sprintf("read %d task files, %d rows\n", length(fs), nrow(d)))
} else {
  d <- utils::read.csv(inp, stringsAsFactors = FALSE)
  cat(sprintf("read %s, %d rows\n", inp, nrow(d)))
}

d <- d[!is.na(d$method), ]
lab <- c(M1 = "MRKCFC", M2 = "MRKCFC-R2", M3 = "MKCFC", M4 = "mclust",
         M5 = "k-means", M6 = "Rk-means", M7 = "lcmm", M8 = "teigen")
# run_all_methods() writes the method_id into `method`; the labels are for display.
d$mid <- d$method
Ms <- sort(unique(d$M)); ns <- sort(unique(d$n))
mids <- intersect(names(lab), unique(d$mid))

# Completeness first. A missing cell is a silent bias -- the expensive cells are
# exactly the ones that time out, and they are the high-M heavy-tail cells the
# whole sweep is about.
cat("\n=== completeness: replicates present per (scenario, n, M) ===\n")
for (nn in ns) {
  cat(sprintf("-- n = %d --\n", nn))
  tb <- table(d$scenario_id[d$n == nn & d$mid == "M2"], d$M[d$n == nn & d$mid == "M2"])
  print(tb)
}
if (any(!is.na(d$error)))
  cat(sprintf("\nerrors: %d rows. top messages:\n", sum(!is.na(d$error))))
if (any(!is.na(d$error)))
  print(utils::head(sort(table(d$error[!is.na(d$error)]), decreasing = TRUE), 5))


# --- Q1: does the M effect hold on all twelve scenarios? ---------------------
cat("\n\n=== Q1. median aRI by method x M, per scenario ===\n")
for (nn in ns) for (s in sort(unique(d$scenario_id))) {
  x <- d[d$n == nn & d$scenario_id == s, ]
  if (!nrow(x)) next
  cat(sprintf("\n-- %s, n = %d, arm %s (%s), d = %s --\n", s, nn,
              x$arm[1], x$arm_label[1], x$d[1]))
  cat(sprintf("%-12s %s %9s\n", "method",
              paste(sprintf("%8s", paste0("M=", Ms)), collapse = ""), "best M"))
  for (m in mids) {
    v <- vapply(Ms, function(MM)
      stats::median(x$aRI[x$mid == m & x$M == MM], na.rm = TRUE), numeric(1))
    bm <- if (all(is.na(v))) NA else Ms[which.max(v)]
    cat(sprintf("%-12s %s %9s\n", lab[[m]],
                paste(sprintf("%8.3f", v), collapse = ""),
                if (is.na(bm)) "--" else bm))
  }
}

# --- Q2: do the comparators gain as much? ------------------------------------
# The M4 - M2 style contrast, pooled over scenarios. If every method gains
# comparably the finding is about the DESIGN and the ranking is unchanged; if
# only the robust second stage gains, it is a claim about the method.
cat("\n\n=== Q2. gain from raising M, pooled over scenarios (median aRI) ===\n")
for (nn in ns) {
  x <- d[d$n == nn, ]
  cat(sprintf("\n-- n = %d --\n", nn))
  cat(sprintf("%-12s %s %10s\n", "method",
              paste(sprintf("%8s", paste0("M=", Ms)), collapse = ""),
              sprintf("M%d-M2", max(Ms))))
  for (m in mids) {
    v <- vapply(Ms, function(MM)
      stats::median(x$aRI[x$mid == m & x$M == MM], na.rm = TRUE), numeric(1))
    g <- if (2 %in% Ms) v[length(v)] - v[Ms == 2] else NA_real_
    cat(sprintf("%-12s %s %+10.3f\n", lab[[m]],
                paste(sprintf("%8.3f", v), collapse = ""), g))
  }
  cat("  lcmm is M-invariant by construction and runs at the anchor M only.\n")
}

# --- the cost and convergence the accuracy is bought with --------------------
cat("\n\n=== convergence rate and median seconds, MKCFC family ===\n")
for (nn in ns) {
  x <- d[d$n == nn & d$mid %in% c("M1", "M2", "M3"), ]
  if (!nrow(x)) next
  cat(sprintf("\n-- n = %d --\n", nn))
  cat(sprintf("%-12s %-6s %s\n", "method", "stat",
              paste(sprintf("%9s", paste0("M=", Ms)), collapse = "")))
  for (m in intersect(c("M1","M2","M3"), mids)) {
    for (st in c("conv", "secs", "iter")) {
      v <- vapply(Ms, function(MM) {
        y <- x[x$mid == m & x$M == MM, ]
        if (!nrow(y)) return(NA_real_)
        switch(st, conv = mean(y$converged, na.rm = TRUE),
                   secs = stats::median(y$seconds, na.rm = TRUE),
                   iter = stats::median(y$iter, na.rm = TRUE))
      }, numeric(1))
      cat(sprintf("%-12s %-6s %s\n", if (st == "conv") lab[[m]] else "", st,
                  paste(sprintf("%9.2f", v), collapse = "")))
    }
  }
  cat("  READ THIS ALONGSIDE Q1. An aRI that improves while convergence falls and\n")
  cat("  iter pins at maxIter is bought from fits that were still moving when the\n")
  cat("  loop stopped. That is the defect maxIter 20 -> 100 was raised to remove.\n")
}

# --- obs_min: the one-factor comparison against Study A ----------------------
cat("\n\n=== the (M = 4, n = 500) rows against Study A ===\n")
cat("EARLIER PREMISE, corrected 2026-09-04. This block was written when Study\n")
cat("A ran at obs_min = 2 and nK_fit = 2, and it printed the (M = 2, n = 500)\n")
cat("rows as a clean one-factor obs_min comparison. Study A was re-run at\n")
cat("obs_min = 4 and nK_fit = 4, so THAT COMPARISON NO LONGER EXISTS: the sweep\n")
cat("and Study A now agree on obs_min, and the (M = 2) rows differ from Study A\n")
cat("in M instead.\n\n")
cat("What the (M = 4, n = 500) rows now are is a REPLICATION CHECK -- same\n")
cat("scenario, same seed, same obs_min, same M, differing only in replicate\n")
cat("count (25 here against Study A's 100). They should agree to sampling error;\n")
cat("a systematic gap means one of the two runs is not what its config says.\n")
cat("Study A medians are in the round-two tables under the results root; no\n")
cat("automatic join is attempted because the two runs are separate files and a\n")
cat("silent mismatch would be worse than a manual read.\n")
if (4 %in% Ms && 500 %in% ns) {
  x <- d[d$M == 4 & d$n == 500, ]
  cat(sprintf("\n%-12s %s\n", "method",
              paste(sprintf("%6s", sort(unique(x$scenario_id))), collapse = "")))
  for (m in mids) {
    v <- vapply(sort(unique(x$scenario_id)), function(s)
      stats::median(x$aRI[x$mid == m & x$scenario_id == s], na.rm = TRUE), numeric(1))
    cat(sprintf("%-12s %s\n", lab[[m]], paste(sprintf("%6.3f", v), collapse = "")))
  }
}

# --- Q3: does the FVE rule pick the same M as the aRI argmax? ----------------
if (!is.na(fvecsv) && file.exists(fvecsv)) {
  fv <- utils::read.csv(fvecsv, stringsAsFactors = FALSE)
  cat("\n\n=== Q3. THE TABLE THE PAPER NEEDS ===\n")
  cat("M chosen by a LABEL-FREE rule against M that maximises aRI. If they agree,\n")
  cat("M is an estimated quantity and the tuning-on-the-truth objection is gone.\n")
  cat("If they disagree, report the FVE choice and let the aRI sweep be the\n")
  cat("sensitivity analysis -- do NOT move the threshold to make them agree.\n\n")
  cat(sprintf("%-5s %10s %10s %10s %10s\n", "scen", "FVE>=.90", "FVE>=.95",
              "argmax aRI", "argmax aRI"))
  cat(sprintf("%-5s %10s %10s %10s %10s\n", "", "(all mth)", "(all mth)",
              "MRKCFC-R2", "MKCFC"))
  for (s in sort(unique(d$scenario_id))) {
    need <- vapply(c(0.90, 0.95), function(thr) {
      # worst variable, worst branch: one M has to serve every method
      k <- vapply(unique(fv$branch), function(br)
        max(vapply(sort(unique(fv$variable)), function(j) {
          v <- vapply(Ms, function(MM) stats::median(
            fv$FVE[fv$scenario_id == s & fv$branch == br & fv$variable == j &
                   fv$M == MM], na.rm = TRUE), numeric(1))
          w <- which(v >= thr)[1]
          if (is.na(w)) max(Ms) + 1L else Ms[w]
        }, integer(1))), integer(1))
      max(k)
    }, integer(1))
    am <- vapply(c("M2", "M3"), function(m) {
      v <- vapply(Ms, function(MM)
        stats::median(d$aRI[d$scenario_id == s & d$mid == m & d$M == MM],
                      na.rm = TRUE), numeric(1))
      if (all(is.na(v))) NA_integer_ else Ms[which.max(v)]
    }, integer(1))
    f <- function(k) if (is.na(k)) "--" else if (k > max(Ms)) paste0(">", max(Ms)) else k
    cat(sprintf("%-5s %10s %10s %10s %10s\n", s, f(need[1]), f(need[2]),
                f(am[1]), f(am[2])))
  }
} else {
  cat("\n\n=== Q3 skipped: no FVE csv passed ===\n")
  cat("Run R/study_fve.R and pass its .csv as the third argument. Without it the\n")
  cat("sweep can say which M performs best but not which M a reader could have\n")
  cat("chosen without the labels, which is the objection it exists to answer.\n")
}

summ <- do.call(rbind, lapply(split(d, list(d$scenario_id, d$n, d$M, d$mid), drop = TRUE),
  function(x) data.frame(
    scenario_id = x$scenario_id[1], arm = x$arm[1], d_sep = x$d[1], n = x$n[1],
    M = x$M[1], obs_min = x$obs_min[1], method_id = x$mid[1], label = lab[[x$mid[1]]],
    reps = nrow(x), n_error = sum(!is.na(x$error)),
    aRI = stats::median(x$aRI, na.rm = TRUE), CCR = stats::median(x$CCR, na.rm = TRUE),
    NMI = stats::median(x$NMI, na.rm = TRUE), MPE = stats::median(x$MPE, na.rm = TRUE),
    RRE = stats::median(x$RRE, na.rm = TRUE),
    conv_rate = mean(x$converged, na.rm = TRUE),
    iter = stats::median(x$iter, na.rm = TRUE),
    seconds = stats::median(x$seconds, na.rm = TRUE),
    stringsAsFactors = FALSE)))
rownames(summ) <- NULL
utils::write.csv(d,    paste0(stem, "_raw.csv"),     row.names = FALSE)
utils::write.csv(summ, paste0(stem, "_summary.csv"), row.names = FALSE)
cat(sprintf("\nwrote %s_raw.csv (%d rows) and %s_summary.csv (%d rows)\n",
            stem, nrow(d), stem, nrow(summ)))
