# Study D: the paper's own component-count rule applied to the real cohort.
#
#   Rscript R/realdata_fve.R <outstem> [scales] [time] [n_sub]
#
#   scales  raw | log | raw+log                   (default raw+log)
#           PLUS-separated, NOT comma-separated, so that the list survives a
#           job wrapper that splits its arguments on commas.
#   time    warped | rescaled                     (default warped)
#   n_sub   SMOKE TEST ONLY: fit a random subset of this many subjects
#           (default 0 = the whole cohort; anything else is not a result)
#
# ===========================================================================
# WHY, 2026-09-12. The application fits nK = 4 components because the FVE rule
# ("retain the smallest M with cumulative FVE >= 0.90 on the worst variable
# and worst branch", study_fve.R) selected 4 on every simulation scenario. The
# rule was never applied to the cohort itself. It matters here for a reason
# the simulation cannot show: KFPCA estimates a subject's scores by least
# squares on that subject's own observations (FPCscoreLSE, predict.KFPCA), and
# 229 of the 1405 subjects have EXACTLY four visits, crowded into the first
# months. Four components from four nearly-collinear observations is an
# ill-conditioned solve; those subjects' scores reach 10^4 IQR units on both
# scales, they are the least stable subjects across restarts under MRKCFC and
# MRKCFC-R2 and not under MKCFC (whose scores are conditional expectations),
# and the K-selection indices computed in that score space are dominated by
# them (the four-visit score problem, Discussion).
#
# If the rule selects M < 4 on this cohort, the application runs at that M by
# the paper's own label-free criterion and every subject has more observations
# than components. If it selects 4 or more, the defect stands and is reported.
# The rule is applied exactly as in study_fve.R: one global fit per marker per
# branch on ALL subjects, spectra with negative eigenvalues truncated at zero,
# cumulative share at M = 1..6.
#
# THIS SCRIPT DOES NOT MOVE THE THRESHOLD. 0.90 is the paper's number; 0.95 is
# printed beside it for the reader, as the simulation table does.
# ===========================================================================

source(file.path("R", "study_c.R"))       # load_config, build_optns via methods.R
source(file.path("R", "realdata_load.R"))

args   <- commandArgs(trailingOnly = TRUE)
stem   <- if (length(args) >= 1) args[1] else "rd_fve"
scales <- if (length(args) >= 2) strsplit(args[2], "[+,]")[[1]] else c("raw", "log")
tm     <- if (length(args) >= 3) args[3] else "warped"
n_sub  <- if (length(args) >= 4) as.integer(args[4]) else 0L
if (!all(scales %in% c("raw", "log"))) stop("scales must be from raw,log")
if (!tm %in% c("warped", "rescaled")) stop("time must be warped or rescaled")
dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)

M_MAX <- 6L
cfg   <- load_config()
ro    <- build_optns(cfg, TRUE)
no    <- build_optns(cfg, FALSE)

# study_fve.R's helper, verbatim: negatives dropped, then the cumulative share.
cum_fve <- function(v) {
  v <- pmax(v, 0); s <- sum(v)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, M_MAX))
  cs <- cumsum(v) / s
  c(cs, rep(1, M_MAX))[seq_len(M_MAX)]
}

dat <- rd_load(window = 36L)
rows <- list()
for (sc in scales) {
  prep <- rd_prepare(dat, log = sc == "log", verbose = TRUE, time = tm)
  idx  <- seq_len(prep$n)
  if (n_sub > 0L && n_sub < prep$n) {
    set.seed(20260912L); idx <- sort(sample.int(prep$n, n_sub))
    cat(sprintf("SMOKE TEST on %d of %d subjects; not a result\n", n_sub, prep$n))
  }
  nv <- lengths(prep$Lt[[1]][idx])
  for (j in seq_along(RD_MARKERS)) {
    m  <- RD_MARKERS[j]
    Lt <- prep$Lt[[j]][idx]; Ly <- prep$Ly[[j]][idx]
    t0 <- Sys.time()
    kf <- tryCatch(MRKCFC:::KFPCA(Lt = Lt, Ly = Ly, interval = ro$interval,
                                 nK = ro$nK, bw = ro$bw, bwK = ro$bwK,
                                 nRegGrid = ro$nRegGrid, fdParobj = ro$basis,
                                 more = FALSE),
                   error = function(e) e)
    t1 <- Sys.time()
    fp <- tryCatch(fdapace::FPCA(Ly, Lt, optns = utils::modifyList(
                     no, list(methodSelectK = M_MAX))),
                   error = function(e) e)
    t2 <- Sys.time()
    kf_ok <- !inherits(kf, "error"); fp_ok <- !inherits(fp, "error")
    cat(sprintf("%s %-6s robust %s (%.0fs, bwK %s) | non-robust %s (%.0fs)\n", sc, m,
                if (kf_ok) "ok" else paste("FAILED:", substr(conditionMessage(kf), 1, 60)),
                as.numeric(difftime(t1, t0, units = "secs")),
                if (kf_ok) sprintf("%.3f", kf$bwK) else "NA",
                if (fp_ok) "ok" else paste("FAILED:", substr(conditionMessage(fp), 1, 60)),
                as.numeric(difftime(t2, t1, units = "secs"))))
    rows[[length(rows) + 1L]] <- rbind(
      data.frame(scale = sc, time = tm, marker = m, branch = "robust (KFPCA)",
                 M = seq_len(M_MAX),
                 FVE = if (kf_ok) cum_fve(kf$values) else rep(NA_real_, M_MAX),
                 bwK = if (kf_ok) as.numeric(kf$bwK) else NA_real_,
                 n = length(idx), stringsAsFactors = FALSE),
      data.frame(scale = sc, time = tm, marker = m, branch = "non-robust (fdapace)",
                 M = seq_len(M_MAX),
                 FVE = if (fp_ok) c(fp$cumFVE, rep(1, M_MAX))[seq_len(M_MAX)]
                       else rep(NA_real_, M_MAX),
                 bwK = NA_real_, n = length(idx), stringsAsFactors = FALSE))
  }
  # The score-step consideration, stated beside the rule so the two are read
  # together: how many subjects a candidate M leaves with no more observations
  # than components.
  cat(sprintf("\n%s scale: subjects with <= M visits, M = 1..%d: %s\n", sc, M_MAX,
              paste(sprintf("M=%d:%d", seq_len(M_MAX),
                            vapply(seq_len(M_MAX), function(M) sum(nv <= M), 0L)),
                    collapse = "  ")))
}
R <- do.call(rbind, rows)
utils::write.csv(R, paste0(stem, ".csv"), row.names = FALSE)

# ----------------------------------------------------------------- report --
rule <- function(worst, thr) { k <- which(worst >= thr); if (length(k)) min(k) else NA_integer_ }
for (sc in scales) {
  cat(sprintf("\n=== cumulative FVE, %s scale, %s time (n = %d) ===\n", sc, tm, R$n[1]))
  cat(sprintf("%-22s %-6s %s\n", "branch", "marker", paste(sprintf("  M=%d", seq_len(M_MAX)), collapse = "")))
  for (br in unique(R$branch)) {
    for (m in RD_MARKERS) {
      v <- R$FVE[R$scale == sc & R$branch == br & R$marker == m]
      cat(sprintf("%-22s %-6s %s\n", br, m, paste(sprintf("%6.3f", v), collapse = "")))
    }
    worst <- vapply(seq_len(M_MAX), function(M) min(R$FVE[R$scale == sc & R$branch == br & R$M == M]), 0)
    cat(sprintf("%-22s %-6s %s   <- worst marker; rule at 0.90 -> M = %s, at 0.95 -> M = %s\n",
                br, "WORST", paste(sprintf("%6.3f", worst), collapse = ""),
                format(rule(worst, 0.90)), format(rule(worst, 0.95))))
  }
  worst_all <- vapply(seq_len(M_MAX), function(M) min(R$FVE[R$scale == sc & R$M == M]), 0)
  cat(sprintf("THE RULE (worst marker, worst branch) on the %s scale: M = %s at 0.90, M = %s at 0.95\n",
              sc, format(rule(worst_all, 0.90)), format(rule(worst_all, 0.95))))
}
cat(sprintf("\nThe simulation's rule selected M = 4 on all twelve scenarios (HPC/fve_summary.csv).\n"))
cat(sprintf("wrote %s\n", paste0(stem, ".csv")))
