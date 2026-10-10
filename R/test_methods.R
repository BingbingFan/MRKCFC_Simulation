# Smoke test for the eight method wrappers, at reduced n so API errors surface
# in minutes rather than hours. Correctness of the metrics is tested separately
# in test_metrics.R; this asks only whether each method runs and returns a
# usable partition.
#
# Run from study/:   Rscript R/test_methods.R [n]

suppressMessages(library(MRKCFC))
source(file.path("R", "generate.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "methods.R"))
for (f in c("gen_scores", "vm_intermediate", "rmvn", "FLEISHMAN", "IQR_SCALE"))
  assign(f, get(f, envir = asNamespace("MRKCFC")), envir = globalenv())

args <- commandArgs(trailingOnly = TRUE)
n_small <- if (length(args)) as.integer(args[1]) else 80L

cfg <- load_config()
cfg$global$n <- n_small
sc <- cfg$scenarios[cfg$scenarios$scenario_id == "A02", ]   # Gaussian, d = 8, easiest
dat <- generate_replicate(cfg, sc, rep = 1L)
cat(sprintf("smoke test at n = %d, scenario %s (arm %d, d = %g)\n\n",
            n_small, sc$scenario_id, sc$arm, sc$d))

K <- as.integer(cfg$global$K_true)
lab <- cfg$methods$label; names(lab) <- cfg$methods$method_id

fail <- 0L
for (mid in names(METHODS)) {
  t0 <- Sys.time()
  res <- tryCatch(METHODS[[mid]](dat$Ly, dat$Lt, K, cfg),
                  error = function(e) conditionMessage(e))
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.character(res)) {
    fail <- fail + 1L
    cat(sprintf("%-3s %-30s **FAIL** %5.1fs  %s\n", mid, lab[[mid]], el, res))
  } else {
    good <- length(res$labels) == n_small && all(res$labels %in% seq_len(K))
    if (!good) fail <- fail + 1L
    cat(sprintf("%-3s %-30s %s %5.1fs  CCR %.3f  aRI %+.3f  sizes %s\n",
                mid, lab[[mid]], if (good) "  ok   " else "**BAD**", el,
                ccr(dat$labels, res$labels), ari(dat$labels, res$labels),
                paste(table(res$labels), collapse = "/")))
  }
}

cat("\n=== common reconstruction, used for MPE and RRE ===\n")
t0 <- Sys.time()
rec <- tryCatch(reconstruct_given_partition(dat$Ly, dat$Lt, dat$labels, cfg),
                error = function(e) conditionMessage(e))
if (is.character(rec)) {
  fail <- fail + 1L
  cat("**FAIL**", rec, "\n")
} else {
  cat(sprintf("  ok    %5.1fs  MPE %.4f  RRE %.4f  RRE(trim) %.4f\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              mpe(rec$y, rec$yhat),
              rre_overall(rec$y, rec$yhat, rec$variable),
              rre_overall(rec$y, rec$yhat, rec$variable, trim = 0.01)))
  cat(sprintf("  scored on the TRUE partition, so RRE well below 1 is expected\n"))
}

cat("\n=== every method starts its fit from the same RNG state ===\n")
# run_all_methods() re-seeds to a common per-cell state before each fit. Without
# it, methods inherit whatever state the previous method left, so toggling
# methods.csv's `include` or reordering METHODS silently moves every method
# downstream, and a recorded seed no longer reproduces a fit.
#
# Tested with PROBES rather than with real methods, deliberately. The obvious
# test -- run M5 first, then second, and compare -- is VACUOUS: k-means at
# nstart = 25 on an easy cell is effectively deterministic, so it passes whether
# or not the fix is present. Verified: with the re-seed removed, M5 was still
# bit-identical run alone, after M4 and after M6. The leak only shows through
# `two-staged` initialisation on the heavy-tailed cells, where a fit costs
# 200-400s, which is far too slow for a smoke test. The probes test the property
# itself and cost nothing.
seen <- new.env(parent = emptyenv())
mk_probe <- function(tag, burn) function(Ly, Lt, K, cfg) {
  assign(tag, .Random.seed, envir = seen)
  if (burn) stats::runif(10000L)          # consume RNG, as a real method does
  list(labels = rep_len(seq_len(K), length(Ly[[1]])), converged = TRUE)
}
METHODS$ZP1 <- mk_probe("p1", burn = TRUE)
METHODS$ZP2 <- mk_probe("p2", burn = FALSE)
invisible(run_all_methods(dat, cfg, methods = c("ZP1", "ZP2")))
METHODS$ZP1 <- METHODS$ZP2 <- NULL
if (!identical(seen$p1, seen$p2)) {
  fail <- fail + 1L
  cat("**FAIL** the second method inherited RNG state from the first\n")
} else {
  cat("  ok    identical .Random.seed at both fits, despite the first burning 10,000 draws\n")
}

cat("\n")
if (fail > 0L) stop(sprintf("method smoke test: %d FAILED", fail), call. = FALSE)
cat("method smoke test: PASS\n")
