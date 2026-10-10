# Checks the metrics against cases with known answers.
# Run from study/:   Rscript R/test_metrics.R

suppressMessages(library(MRKCFC))
source(file.path("R", "metrics.R"))

fail <- 0L
ok <- function(lab, cond, extra = "") {
  good <- isTRUE(cond); if (!good) fail <<- fail + 1L
  cat(sprintf("%-56s %s %s\n", lab, if (good) "PASS" else "**FAIL**", extra))
}
near <- function(a, b, tol = 1e-8) isTRUE(abs(a - b) < tol)

set.seed(1)
truth <- rep(1:2, each = 50)

cat("=== perfect and relabelled agreement ===\n")
ok("CCR = 1 on an exact match",   near(ccr(truth, truth), 1))
ok("CCR = 1 under relabelling",   near(ccr(truth, 3L - truth), 1),
   "(label switching must not be penalised)")
ok("aRI = 1 under relabelling",   near(ari(truth, 3L - truth), 1))
ok("NMI = 1 under relabelling",   near(nmi(truth, 3L - truth), 1))
ok("AMI = 1 under relabelling",   near(ami(truth, 3L - truth), 1))

cat("\n=== independence ===\n")
rand <- sample(rep(1:2, each = 50))
cat(sprintf("   random partition: aRI %+.4f  NMI %.4f  AMI %+.4f\n",
            ari(truth, rand), nmi(truth, rand), ami(truth, rand)))
ok("aRI near 0 under independence", abs(ari(truth, rand)) < 0.1)
ok("AMI near 0 under independence", abs(ami(truth, rand)) < 0.1)

cat("\n=== the reason AMI was added ===\n")
# A method that over-splits: correct structure, but 4 clusters where truth has 2.
oversplit <- rep(1:4, each = 25)
cat(sprintf("   over-split into 4: aRI %+.4f  NMI %.4f  AMI %+.4f\n",
            ari(truth, oversplit), nmi(truth, oversplit), ami(truth, oversplit)))
ok("plain NMI rewards over-splitting more than AMI does",
   nmi(truth, oversplit) > ami(truth, oversplit),
   "(this is why NMI must be read next to a corrected measure)")

cat("\n=== chance correction is doing work ===\n")
# many small random clusters: uncorrected NMI stays high, AMI collapses
many <- sample(rep(1:20, each = 5))
cat(sprintf("   20 random clusters: NMI %.4f  AMI %+.4f\n",
            nmi(truth, many), ami(truth, many)))
ok("AMI penalises random over-clustering, NMI does not",
   ami(truth, many) < 0.1 && nmi(truth, many) > ami(truth, many))

cat("\n=== reconstruction error ===\n")
v <- rep(c("A", "B", "C"), each = 200)
y <- c(rnorm(200, 5, 2), rnorm(200, 30, 12), rnorm(200, 20, 6))
ok("MPE = 0 for a perfect fit", near(mpe(y, y), 0))
ok("RRE = 0 for a perfect fit", near(rre_overall(y, y, v), 0))
mu <- ave(y, v, FUN = mean)
ok("RRE = 1 when predicting each variable's mean", near(rre_overall(y, mu, v), 1),
   "(the conventional reference point)")
ok("MPE is NOT 1 there, being unit-dependent", !near(mpe(y, mu), 1),
   sprintf("MPE = %.2f", mpe(y, mu)))
noisy <- mu + rnorm(600)
ok("trimmed RRE runs and is finite", is.finite(rre_overall(y, noisy, v, trim = 0.01)))

cat("\n=== assembly ===\n")
r <- score_fit(truth, 3L - truth, y = y, yhat = mu, variable = v)
ok("all seven columns present",
   identical(sort(names(r)), sort(c("CCR","aRI","NMI","AMI","MPE","RRE","RRE_trim"))))
r2 <- score_fit(truth, 3L - truth)
ok("reconstruction columns are NA when no fit is supplied",
   all(is.na(c(r2$MPE, r2$RRE, r2$RRE_trim))) && !is.na(r2$CCR),
   "(never a fabricated zero)")

cat("\n")
if (fail > 0L) stop(sprintf("metric tests: %d FAILED", fail), call. = FALSE)
cat("metric tests: PASS\n")
