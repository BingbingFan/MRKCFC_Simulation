# Study B smoke test. Run from study/:   Rscript R/test_study_b.R [n]
#
# Three things, in increasing cost:
#   1. every configuration generates data of the shape it claims
#   2. Delta does not drift with P  <- the assumption the whole P axis rests on
#   3. the driver's own code path runs a cell end to end
#
# Check 2 is the one that matters. If added variables are not uninformative,
# Delta rises with P, accuracy improves for reasons unrelated to the method, and
# the P axis measures the opposite of what it claims to.

suppressMessages(library(MRKCFC))
source(file.path("R", "generate.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "methods.R"))
source(file.path("R", "separation.R"))
for (f in c("gen_scores", "vm_intermediate", "rmvn", "FLEISHMAN", "IQR_SCALE"))
  assign(f, get(f, envir = asNamespace("MRKCFC")), envir = globalenv())

args <- commandArgs(trailingOnly = TRUE)
n_small <- if (length(args)) as.integer(args[1]) else 200L

cfg <- load_config()
configs <- utils::read.csv(file.path("config", "scenarios_b.csv"),
                           stringsAsFactors = FALSE)
fail <- 0L
ok <- function(lab, cond, extra = "") {
  good <- isTRUE(cond); if (!good) fail <<- fail + 1L
  cat(sprintf("  %-46s %s %s\n", lab, if (good) "ok  " else "FAIL", extra))
}

# apply_config_b() and as_scenario_b() come from R/generate.R, so this test
# exercises the SAME helpers run_cell_b.R uses rather than a copy of them.

cat("=== 1. every configuration generates what it claims ===\n")
for (i in seq_len(nrow(configs))) {
  row <- configs[i, ]
  c2 <- apply_config_b(cfg, row); c2$global$n <- n_small   # reduced n for speed
  dat <- try(generate_replicate(c2, as_scenario_b(row), rep = 1L), silent = TRUE)
  if (inherits(dat, "try-error")) {
    ok(sprintf("%s %s", row$config_id, row$label), FALSE,
       sub("\n.*", "", as.character(dat)))
    next
  }
  pi_ <- as.numeric(strsplit(row$pi, ",")[[1]])
  shape_ok <- length(dat$Ly) == row$P &&
              length(dat$Ly[[1]]) == n_small &&
              length(dat$labels) == n_small
  # observation counts should track the configured Poisson mean
  lam_hat <- mean(dat$n_i)
  lam_ok <- abs(lam_hat - row$obs_poisson_lambda) < 1.5
  # class balance should track pi
  p_hat <- mean(dat$labels == 1L)
  pi_ok <- abs(p_hat - pi_[1]) < 0.12
  ok(sprintf("%s %-14s P=%2d n=%d", row$config_id, row$label, row$P, n_small),
     shape_ok && lam_ok && pi_ok,
     sprintf("lambda_hat %.1f, pi_hat %.2f", lam_hat, p_hat))
}

cat("\n=== 2. Delta must NOT drift with P (added variables are noise) ===\n")
cat("  If Delta rises with P, the added variables carry signal and the P axis\n")
cat("  measures added information rather than added dimension.\n")
cat("  Run at the CONFIGURED n, not the reduced one, for the reason below.\n\n")
cat(sprintf("  %-6s %-4s %6s %8s %10s\n", "config", "P", "n", "delta", "ref_ccr"))
dl <- c()
for (cid in c("B01", "B11", "B12")) {
  row <- configs[configs$config_id == cid, ]
  c2 <- apply_config_b(cfg, row)          # configured n, deliberately not n_small
  dat <- generate_replicate(c2, as_scenario_b(row), rep = 1L)
  s <- separation_of(dat, c2)
  dl <- c(dl, s$delta)
  cat(sprintf("  %-6s %-4d %6d %8.3f %10.3f\n", cid, row$P, row$n, s$delta, s$bayes))
  utils::flush.console()
}
names(dl) <- c("P3", "P5", "P10")
# WHY THIS RUNS AT FULL n. The sample Mahalanobis distance between two class
# centroids is upward biased by roughly p/n, and the score space is P * M
# dimensional, so the bias grows with P. Measured 2026-08-28 on B01 against B12:
#
#     n = 200   delta(P=10)/delta(P=3) = 1.496
#     n = 500   delta(P=10)/delta(P=3) = 0.957
#
# At the reduced n this check fails on the ESTIMATOR while the data is perfectly
# uninformative; at the configured n the drift is gone. The first version of
# this test ran at n_small and reported a false positive.
#
# The same caution applies when reading results: `delta` is not comparable
# across cells that differ in both n and P. No Study B cell combines large P
# with small n, so it does not bite here.
ok("Delta at P=10 within 25% of Delta at P=3",
   is.finite(dl["P10"]) && is.finite(dl["P3"]) &&
     abs(dl["P10"] / dl["P3"] - 1) < 0.25,
   sprintf("ratio %.2f (sqrt(10/3) = 1.83 would mean informative)",
           dl["P10"] / dl["P3"]))

cat("\n=== 3. the driver's own path, one cell end to end ===\n")
row <- configs[configs$config_id == "B01", ]
c2 <- apply_config_b(cfg, row); c2$global$n <- n_small
dat <- generate_replicate(c2, as_scenario_b(row), rep = 1L)
t0 <- Sys.time()
# M8 only: the MKCFC family costs 200-400s a fit and this is a smoke test.
rows <- try(run_all_methods(dat, c2, methods = "M8"), silent = TRUE)
if (inherits(rows, "try-error")) {
  ok("run_all_methods on a B config", FALSE, sub("\n.*", "", as.character(rows)))
} else {
  ok("run_all_methods on a B config", nrow(rows) == 1L && !is.na(rows$aRI),
     sprintf("aRI %+.3f in %.1fs", rows$aRI,
             as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  ok("iter column present for the combine step", "iter" %in% names(rows))
}

cat("\n")
if (fail > 0L) stop(sprintf("Study B smoke test: %d FAILED", fail), call. = FALSE)
cat("Study B smoke test: PASS\n")
