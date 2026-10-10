# Read the FVE study and report the M that a LABEL-FREE rule selects.
#
#   Rscript R/combine_fve.R <dir-or-csv> [outstem] [sweep-summary-csv]
#
#   Rscript R/combine_fve.R $ROOT/results_fve $ROOT/results_fve/fve
#   Rscript R/combine_fve.R $ROOT/results_fve $ROOT/results_fve/fve \
#           $ROOT/results_sweep/n500/sweep_summary.csv
#
# Pass the sweep summary as a third argument to get the last table, which is the
# one the paper needs: the M chosen without labels against the M that maximises
# aRI. Agreement makes M an estimated quantity and removes the
# tuning-on-the-truth objection. Disagreement means reporting the FVE choice and
# demoting the sweep to sensitivity -- NOT moving the threshold until they
# agree.

args  <- commandArgs(trailingOnly = TRUE)
inp   <- if (length(args) >= 1) args[1] else "results_fve"
stem  <- if (length(args) >= 2) args[2] else file.path(inp, "fve")
swcsv <- if (length(args) >= 3) args[3] else NA_character_

if (dir.exists(inp)) {
  fs <- list.files(inp, pattern = "^fve_task_[0-9]+[.]rds$", full.names = TRUE)
  if (!length(fs)) stop("no fve_task_*.rds under ", inp)
  d <- do.call(rbind, lapply(fs, readRDS))
  cat(sprintf("read %d task files, %d rows\n", length(fs), nrow(d)))
} else {
  d <- utils::read.csv(inp, stringsAsFactors = FALSE)
  cat(sprintf("read %s, %d rows\n", inp, nrow(d)))
}

d   <- d[!is.na(d$branch), ]
Ms  <- sort(unique(d$M))
brs <- sort(unique(d$branch))
scs <- sort(unique(d$scenario_id))
vs  <- sort(unique(d$variable))
cat(sprintf("n = %s | obs_min = %s | scenarios %d | variables %d | branches %d\n",
            paste(unique(d$n), collapse = ","), paste(unique(d$obs_min), collapse = ","),
            length(scs), length(vs), length(brs)))

# Completeness. A short cell is a missing task, not a result.
cat("\n=== completeness: replicates per scenario (should be equal throughout) ===\n")
print(table(d$scenario_id[d$M == 1 & d$variable == vs[1]],
            d$branch[d$M == 1 & d$variable == vs[1]]))
if (any(is.na(d$FVE)))
  cat(sprintf("\nNA FVE: %d rows (%.1f%%) -- a failed stage-one fit\n",
              sum(is.na(d$FVE)), 100 * mean(is.na(d$FVE))))

med <- function(s, br, j, MM)
  stats::median(d$FVE[d$scenario_id == s & d$branch == br &
                      d$variable == j & d$M == MM], na.rm = TRUE)

# --- the curves ---------------------------------------------------------------
for (br in brs) {
  cat(sprintf("\n\n=== median cumulative FVE, %s ===\n", br))
  cat(sprintf("%-5s %-4s %s\n", "scen", "var",
              paste(sprintf("%8s", paste0("M=", Ms)), collapse = "")))
  for (s in scs) for (j in vs)
    cat(sprintf("%-5s %-4d %s\n", s, j,
                paste(sprintf("%8.3f", vapply(Ms, function(MM) med(s, br, j, MM),
                                              numeric(1))), collapse = "")))
}

# Smallest M reaching `thr`. The aggregation is the one real choice in this
# file: WORST VARIABLE governs, because the method takes ONE nK_fit and a single
# M has to serve all P variables. Reported per branch, then across branches,
# because a single M also has to serve every method being compared.
need_M <- function(s, br, thr) {
  k <- vapply(vs, function(j) {
    v <- vapply(Ms, function(MM) med(s, br, j, MM), numeric(1))
    w <- which(v >= thr)[1]
    if (is.na(w)) max(Ms) + 1L else as.integer(Ms[w])
  }, integer(1))
  max(k)
}
fmt <- function(k) if (is.na(k)) "--" else if (k > max(Ms)) paste0(">", max(Ms)) else as.character(k)

cat("\n\n=== smallest M reaching each threshold (worst variable governs) ===\n")
cat(sprintf("%-5s %-22s %10s %10s %10s\n", "scen", "branch",
            "FVE>=0.90", "FVE>=0.95", "FVE>=0.99"))
for (s in scs) for (br in brs)
  cat(sprintf("%-5s %-22s %10s %10s %10s\n", s, br,
              fmt(need_M(s, br, 0.90)), fmt(need_M(s, br, 0.95)),
              fmt(need_M(s, br, 0.99))))

cat("\n\n=== ONE M FOR EVERY METHOD: worst variable AND worst branch ===\n")
cat("The robust and non-robust branches disagree, and the methods are being\n")
cat("COMPARED, so they must all receive the same M. The maximum is the choice\n")
cat("that satisfies the rule for every method at once.\n\n")
cat(sprintf("%-5s %10s %10s %10s\n", "scen", "FVE>=0.90", "FVE>=0.95", "FVE>=0.99"))
allm <- t(vapply(scs, function(s)
  vapply(c(0.90, 0.95, 0.99), function(thr)
    max(vapply(brs, function(br) need_M(s, br, thr), integer(1))), integer(1)),
  integer(3)))
for (i in seq_along(scs))
  cat(sprintf("%-5s %10s %10s %10s\n", scs[i], fmt(allm[i,1]), fmt(allm[i,2]), fmt(allm[i,3])))
cat(sprintf("\n%-5s %10s %10s %10s   <- the rule's answer, if one M is fixed\n", "MODE",
            fmt(as.integer(names(which.max(table(allm[,1]))))),
            fmt(as.integer(names(which.max(table(allm[,2]))))),
            fmt(as.integer(names(which.max(table(allm[,3])))))))
cat(sprintf("%-5s %10s %10s %10s   <- if the rule must hold on EVERY scenario\n", "MAX",
            fmt(max(allm[,1])), fmt(max(allm[,2])), fmt(max(allm[,3]))))

# --- the table the paper needs ------------------------------------------------
if (!is.na(swcsv) && file.exists(swcsv)) {
  sw <- utils::read.csv(swcsv, stringsAsFactors = FALSE)
  cat("\n\n=== LABEL-FREE M AGAINST THE aRI ARGMAX ===\n")
  cat("Stage 1 established that M >= 5 is unreportable: MRKCFC-R2 convergence is\n")
  cat("0.05 at M = 5 and 0.00 at M = 6, with median iter pinned on the maxIter\n")
  cat("cap. An aRI argmax of 5 or 6 below is an ARTEFACT, not a recommendation.\n")
  cat("The usable range is M <= 4.\n\n")
  cat(sprintf("%-5s %9s %9s %11s %11s %11s\n", "scen", "FVE>=.90", "FVE>=.95",
              "argmax R2", "argmax MKCFC", "R2 conv@4"))
  for (i in seq_along(scs)) {
    s <- scs[i]
    am <- vapply(c("M2", "M3"), function(m) {
      v <- vapply(Ms, function(MM)
        sw$aRI[sw$scenario_id == s & sw$method_id == m & sw$M == MM][1], numeric(1))
      if (all(is.na(v))) NA_integer_ else as.integer(Ms[which.max(v)])
    }, integer(1))
    cv <- sw$conv_rate[sw$scenario_id == s & sw$method_id == "M2" & sw$M == 4][1]
    cat(sprintf("%-5s %9s %9s %11s %11s %11.2f\n", s,
                fmt(allm[i,1]), fmt(allm[i,2]), fmt(am[1]), fmt(am[2]),
                if (length(cv)) cv else NA_real_))
  }
  cat("\nIf the FVE column says 4 and the convergence boundary says 4, M = 4 is\n")
  cat("chosen by a rule a reader can apply, not by the outcome it produces.\n")
} else {
  cat("\n\n=== sweep comparison skipped: no sweep summary passed ===\n")
  cat("Pass results_sweep/n500/sweep_summary.csv as the third argument.\n")
}

summ <- do.call(rbind, lapply(split(d, list(d$scenario_id, d$branch, d$variable, d$M),
                                   drop = TRUE), function(x) data.frame(
  scenario_id = x$scenario_id[1], arm = x$arm[1], branch = x$branch[1],
  variable = x$variable[1], M = x$M[1], n = x$n[1], obs_min = x$obs_min[1],
  reps = nrow(x), FVE = stats::median(x$FVE, na.rm = TRUE),
  FVE_q25 = stats::quantile(x$FVE, 0.25, na.rm = TRUE),
  FVE_q75 = stats::quantile(x$FVE, 0.75, na.rm = TRUE),
  stringsAsFactors = FALSE)))
rownames(summ) <- NULL
utils::write.csv(d,    paste0(stem, "_raw.csv"),     row.names = FALSE)
utils::write.csv(summ, paste0(stem, "_summary.csv"), row.names = FALSE)
cat(sprintf("\nwrote %s_raw.csv (%d rows) and %s_summary.csv (%d rows)\n",
            stem, nrow(d), stem, nrow(summ)))
