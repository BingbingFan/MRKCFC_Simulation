# Read C4-at-M=4 and answer the one question: cosmetic or substantive?
#
#   Rscript R/combine_c4b.R <dir-or-csv> [outstem]

args <- commandArgs(trailingOnly = TRUE)
inp  <- if (length(args) >= 1) args[1] else "results_c4b"
stem <- if (length(args) >= 2) args[2] else file.path(inp, "c4b")

if (dir.exists(inp)) {
  fs <- list.files(inp, pattern = "^c4b_task_[0-9]+[.]rds$", full.names = TRUE)
  if (!length(fs)) stop("no c4b_task_*.rds under ", inp)
  d <- do.call(rbind, lapply(fs, readRDS))
  cat(sprintf("read %d task files, %d rows\n", length(fs), nrow(d)))
} else {
  d <- utils::read.csv(inp, stringsAsFactors = FALSE)
}

d   <- d[!is.na(d$conv_prop), ]
lab <- c(M1 = "MRKCFC", M2 = "MRKCFC-R2", M3 = "MKCFC")
CP  <- sort(unique(d$conv_prop)); scs <- sort(unique(d$scenario_id))
mids <- intersect(names(lab), unique(d$method_id))
cat(sprintf("nK_fit = %s | maxIter = %s | conv.prop %s\n",
            paste(unique(d$nK_fit), collapse = ","),
            paste(unique(d$maxIter), collapse = ","),
            paste(CP, collapse = ", ")))

cat("\n=== completeness: replicates per scenario x method (all conv.prop) ===\n")
print(table(d$scenario_id[d$conv_prop == 0], d$method_id[d$conv_prop == 0]))
if (any(!is.na(d$error)))
  cat(sprintf("\nerrors: %d rows\n", sum(!is.na(d$error))))

g <- function(sc, m, cp, col, f) f(d[[col]][d$scenario_id == sc & d$method_id == m &
                                           d$conv_prop == cp], na.rm = TRUE)

cat("\n\n=== THE ANSWER: agreement with the conv.prop = 0 partition ===\n")
cat("Median aRI between each relaxed fit and the STRICT fit on the SAME\n")
cat("replicate. Near 1.0 means the threshold changed nothing, so the\n")
cat("non-convergence is a handful of boundary subjects oscillating -- COSMETIC.\n")
cat("Well below 1.0 means the partition genuinely had not settled.\n\n")
for (m in mids) {
  cat(sprintf("-- %s --\n", lab[[m]]))
  cat(sprintf("%-5s %s\n", "scen",
              paste(sprintf("%12s", paste0("cp=", CP[-1])), collapse = "")))
  for (sc in scs)
    cat(sprintf("%-5s %s\n", sc, paste(sprintf("%12.3f",
        vapply(CP[-1], function(cp) g(sc, m, cp, "agree_strict", stats::median),
               numeric(1))), collapse = "")))
}

cat("\n\n=== convergence rate, and what it costs in accuracy ===\n")
for (m in mids) {
  cat(sprintf("\n-- %s --\n", lab[[m]]))
  cat(sprintf("%-5s %-6s %s\n", "scen", "stat",
              paste(sprintf("%10s", paste0("cp=", CP)), collapse = "")))
  for (sc in scs) for (st in c("conv", "aRI", "iter")) {
    v <- vapply(CP, function(cp) switch(st,
          conv = g(sc, m, cp, "converged", mean),
          aRI  = g(sc, m, cp, "aRI", stats::median),
          iter = g(sc, m, cp, "iter", stats::median)), numeric(1))
    cat(sprintf("%-5s %-6s %s\n", if (st == "conv") sc else "", st,
                paste(sprintf("%10.3f", v), collapse = "")))
  }
}

cat("\n\n=== how to read this, and the trap ===\n")
cat("PRIMARY RESULT IS conv.prop = 0. The relaxed columns are a SENSITIVITY\n")
cat("ANALYSIS, never the headline. Reporting a convergence rate obtained by\n")
cat("loosening the threshold until the method converged is loosening the\n")
cat("definition of success until it was met.\n\n")
cat("What may be claimed, IF agree_strict is near 1.0: the partition is stable\n")
cat("and the loop's failure to declare convergence reflects a small number of\n")
cat("boundary subjects reassigning indefinitely, not an unsettled answer.\n\n")
cat("Check the control scenario. If a threshold moves partitions where the loop\n")
cat("ALREADY converges at cp = 0, the instrument is suspect and nothing else\n")
cat("here can be trusted.\n")

summ <- do.call(rbind, lapply(split(d, list(d$scenario_id, d$method_id, d$conv_prop),
                                    drop = TRUE), function(x) data.frame(
  scenario_id = x$scenario_id[1], arm = x$arm[1], arm_label = x$arm_label[1],
  d_sep = x$d[1], method_id = x$method_id[1], label = lab[[x$method_id[1]]],
  conv_prop = x$conv_prop[1], nK_fit = x$nK_fit[1], reps = nrow(x),
  n_error = sum(!is.na(x$error)),
  aRI = stats::median(x$aRI, na.rm = TRUE), CCR = stats::median(x$CCR, na.rm = TRUE),
  agree_strict = stats::median(x$agree_strict, na.rm = TRUE),
  conv_rate = mean(x$converged, na.rm = TRUE),
  iter = stats::median(x$iter, na.rm = TRUE),
  seconds = stats::median(x$seconds, na.rm = TRUE), stringsAsFactors = FALSE)))
rownames(summ) <- NULL
utils::write.csv(d,    paste0(stem, "_raw.csv"),     row.names = FALSE)
utils::write.csv(summ, paste0(stem, "_summary.csv"), row.names = FALSE)
cat(sprintf("\nwrote %s_raw.csv (%d rows) and %s_summary.csv (%d rows)\n",
            stem, nrow(d), stem, nrow(summ)))
