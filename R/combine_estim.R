# Read the estimation study. The one reader; the driver summarises nothing.
#
#   Rscript R/combine_estim.R <dir-or-csv> [outstem]

args <- commandArgs(trailingOnly = TRUE)
inp  <- if (length(args) >= 1) args[1] else "results_estim"
stem <- if (length(args) >= 2) args[2] else file.path(inp, "estim")

if (dir.exists(inp)) {
  fs <- list.files(inp, pattern = "^estim_task_[0-9]+[.]rds$", full.names = TRUE)
  if (!length(fs)) stop("no estim_task_*.rds under ", inp)
  d <- do.call(rbind, lapply(fs, readRDS))
  cat(sprintf("read %d task files, %d rows\n", length(fs), nrow(d)))
} else {
  d <- utils::read.csv(inp, stringsAsFactors = FALSE)
}

d <- d[!is.na(d$variant), ]
lab <- c(M1 = "MRKCFC (moment)", M2 = "MRKCFC-R2 (sign)", M3 = "MKCFC (moment)")
scs <- sort(unique(d$scenario_id))
vids <- intersect(names(lab), unique(d$variant))
cat(sprintf("nK_fit = %s (the FIT) | M_true = %s (the TRUTH) | n = %s | obs_min = %s\n",
            paste(unique(d$nK_fit), collapse = ","),
            paste(unique(d$M_true), collapse = ","),
            paste(unique(d$n), collapse = ","),
            paste(unique(d$obs_min), collapse = ",")))
if (length(unique(d$nK_fit)) > 1)
  cat("WARNING: mixed nK_fit in one directory. These are not comparable.\n")
if (length(unique(d$nK_used)) > 1)
  cat(sprintf("note: per-variable component counts differ (%s) -- fdapace selectK\n",
              paste(sort(unique(na.omit(d$nK_used))), collapse = ",")))

cat("\n=== completeness: replicates per scenario x variant (2 classes each) ===\n")
print(table(d$scenario_id, d$variant) / 2)
if (any(!is.na(d$error))) {
  cat(sprintf("\nerrors: %d rows\n", sum(!is.na(d$error))))
  print(utils::head(sort(table(d$error[!is.na(d$error)]), decreasing = TRUE), 5))
}

# Ellipticity governs where the theory applies at all. Arms 2 and 3 are
# Fleishman via Vale-Maurelli and are NOT elliptical, so the SSCM has no
# guarantee there and a loss on those arms is not evidence against the method.
ell <- c(`1` = "yes", `2` = "NO", `3` = "NO", `4` = "yes", `5` = "yes", `6` = "yes")

show <- function(metric, title, better) {
  cat(sprintf("\n\n=== %s ===\n(%s)\n", title, better))
  cat(sprintf("%-4s %-6s %-24s %s %10s\n", "arm", "ellip?", "scenario",
              paste(sprintf("%18s", lab[vids]), collapse = ""), "R2 - moment"))
  for (s in scs) {
    v <- vapply(vids, function(x) stats::median(
      d[[metric]][d$scenario_id == s & d$variant == x], na.rm = TRUE), numeric(1))
    arm <- as.character(d$arm[d$scenario_id == s][1])
    gain <- if (all(c("M1","M2") %in% vids)) v[["M1"]] - v[["M2"]] else NA_real_
    cat(sprintf("%-4s %-6s %-24s %s %10.2f\n", arm, ell[[arm]], s,
                paste(sprintf("%18.2f", v), collapse = ""), gain))
  }
}
show("angle", "eigenfunction recovery, leading 2 (ILL-CONDITIONED: 6.4% eigen-gap)",
     "degrees, LOWER is better; last column positive = sign beats moment")
show("angle_lead3", "eigenfunction recovery, leading 3 (well conditioned: 48% gap)",
     "degrees, LOWER is better; last column positive = sign beats moment")
show("angle_lead4", "eigenfunction recovery, leading 4 (15% gap: 3x worse conditioned)",
     "degrees, LOWER is better; sign beats moment 12/12 here too, but the three-way ranking does NOT survive: it shifts on A07, A09, A10")
show("score_cc", "score recovery, mean canonical correlation",
     "HIGHER is better; last column NEGATIVE = sign beats moment")

cat("\n\n=== how to read this ===\n")
cat("The three truncations answer different questions. Z's eigenvalues, with the\n")
cat("relative gap BELOW each:\n\n")
cat("    9.154   8.470   7.927  |  4.139   3.519   3.292\n")
cat("         7.5%    6.4%   47.8%     15.0%    6.5%\n\n")
cat("By Davis-Kahan the error in an estimated eigen-subspace scales as\n")
cat("perturbation over eigen-gap, so the 47.8 percent drop after the third is the\n")
cat("only clean cut: leading-2 (6.4%, where Step 3(d) truncates) is barely\n")
cat("identifiable, leading-3 is well separated, and leading-4 (15.0%) is about 3x\n")
cat("worse conditioned than leading-3.\n\n")
cat("WHAT LEAD4 ACTUALLY SHOWED, once measured (2026-09-04). It was added to\n")
cat("answer 'the study fits M=4, so why score 3?', on the expectation that it\n")
cat("would reproduce the lead3 ordering. IT DOES NOT:\n\n")
cat("  sign vs moment (M2 vs M1)  12/12 at lead3 AND 12/12 at lead4 -- unchanged\n")
cat("  three-way winner           flips on A07, A09, A10, all elliptical, where\n")
cat("                             M2's lead3 win over M3 reverses\n")
cat("  every angle degrades       sign 14-16 deg -> 35-47; M3 10-47 -> 27-68\n\n")
cat("The degradation is what a 15.0 percent gap predicts against 47.8. So lead4\n")
cat("supports the SIGN-AGAINST-MOMENT claim, which is what this study exists to\n")
cat("test, and does NOT support a truncation-invariant three-way ranking. Do not\n")
cat("promote it to the headline, and do not state stability it does not have.\n\n")
cat("  winning at lead3 but not lead2 -> the estimator is sound and the TARGET\n")
cat("                                    is unidentifiable at that truncation\n")
cat("  losing at both                 -> the estimator is not helping\n\n")
cat("ELLIPTICITY GATES THE CLAIM. The SSCM's eigenvectors coincide with the\n")
cat("covariance's only for elliptically symmetric distributions. Arms 2 and 3 are\n")
cat("NOT elliptical, so no guarantee applies and a loss there is not evidence\n")
cat("against the method. Report those arms anyway: a robustness claim that holds\n")
cat("only under ellipticity should say so.\n")

summ <- do.call(rbind, lapply(split(d, list(d$scenario_id, d$variant), drop = TRUE),
  function(x) data.frame(
    scenario_id = x$scenario_id[1], arm = x$arm[1], arm_label = x$arm_label[1],
    d_sep = x$d[1], elliptical = ell[[as.character(x$arm[1])]],
    variant = x$variant[1], label = lab[[x$variant[1]]],
    nK_fit = x$nK_fit[1], M_true = x$M_true[1], n = x$n[1], obs_min = x$obs_min[1],
    # The per-variable component count the fit ACTUALLY returned. Carried
    # through so a silent reversion of nK_fit, or a ragged fdapace selectK,
    # is visible in the summary rather than only in the raw files.
    nK_used = paste(sort(unique(stats::na.omit(x$nK_used))), collapse = "/"),
    fits = nrow(x), n_error = sum(!is.na(x$error)),
    angle = stats::median(x$angle, na.rm = TRUE),
    angle_lead3 = stats::median(x$angle_lead3, na.rm = TRUE),
    angle_lead4 = stats::median(x$angle_lead4, na.rm = TRUE),
    score_cc = stats::median(x$score_cc, na.rm = TRUE),
    stringsAsFactors = FALSE)))
rownames(summ) <- NULL
utils::write.csv(d,    paste0(stem, "_raw.csv"),     row.names = FALSE)
utils::write.csv(summ, paste0(stem, "_summary.csv"), row.names = FALSE)
cat(sprintf("\nwrote %s_raw.csv (%d rows) and %s_summary.csv (%d rows)\n",
            stem, nrow(d), stem, nrow(summ)))
