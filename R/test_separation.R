# Checks the robust separation measure against cases with known answers.
# Run from study/:   Rscript R/test_separation.R

suppressMessages(library(MRKCFC))
source(file.path("R", "generate.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "methods.R"))
source(file.path("R", "separation.R"))
for (f in c("gen_scores", "vm_intermediate", "rmvn", "FLEISHMAN", "IQR_SCALE"))
  assign(f, get(f, envir = asNamespace("MRKCFC")), envir = globalenv())

fail <- 0L
ok <- function(lab, cond, extra = "") {
  good <- isTRUE(cond); if (!good) fail <<- fail + 1L
  cat(sprintf("%-58s %s %s\n", lab, if (good) "PASS" else "**FAIL**", extra))
}

set.seed(11)
n <- 600L; p <- 4L
mk <- function(shift) {
  lab <- rep(1:2, each = n / 2)
  X <- matrix(rnorm(n * p), n, p)
  X[lab == 2, 1] <- X[lab == 2, 1] + shift
  list(X = X, lab = lab)
}

cat("=== recovers a known separation ===\n")
# identity scatter and a shift of `s` on one coordinate gives Delta = s
for (s in c(1, 2, 4)) {
  d <- mk(s)
  got <- robust_separation(d$X, d$lab)$min
  ok(sprintf("shift %g gives Delta near %g", s, s), abs(got - s) < 0.25,
     sprintf("got %.3f", got))
}

cat("\n=== behaves as a distance ===\n")
d <- mk(2)
r <- robust_separation(d$X, d$lab)
ok("symmetric", isTRUE(all.equal(r$matrix, t(r$matrix))))
ok("zero on the diagonal", all(abs(diag(r$matrix)) < 1e-8))
ok("no separation gives Delta near 0", robust_separation(mk(0)$X, mk(0)$lab)$min < 0.35)

cat("\n=== scale invariance, which is why Mahalanobis and not Euclidean ===\n")
d <- mk(2)
big <- d$X; big[, 2] <- big[, 2] * 50   # inflate an irrelevant coordinate
ok("unchanged when an unrelated coordinate is rescaled",
   abs(robust_separation(big, d$lab)$min - robust_separation(d$X, d$lab)$min) < 0.25,
   "(a Euclidean distance would be swamped)")

cat("\n=== defined without moments, which Euclidean-on-covariance is not ===\n")
lab <- rep(1:2, each = n / 2)
Xh <- matrix(rt(n * p, df = 1.2), n, p)     # infinite variance
Xh[lab == 2, 1] <- Xh[lab == 2, 1] + 3
r <- robust_separation(Xh, lab)
ok("finite Delta under an infinite-variance law", is.finite(r$min),
   sprintf("Delta = %.2f", r$min))
ok("a covariance-based scatter would be unusable here",
   !is.finite(stats::var(Xh[, 1])) || stats::var(Xh[, 1]) > 50,
   sprintf("sample var = %.0f, unstable by construction", stats::var(Xh[, 1])))

cat("\n=== Bayes reference ===\n")
ok("reference CCR at Delta = 0 is chance", abs(reference_ccr(0) - 0.5) < 1e-9)
ok("Delta = 4 implies about 0.977", abs(reference_ccr(4) - 0.9772) < 1e-3)
ok("monotone in Delta", reference_ccr(5) > reference_ccr(3))
ok("NA for K > 2 rather than a wrong number", is.na(reference_ccr(3, K = 3L)))

cat("\n=== more classes: the minimum pair is what binds ===\n")
lab3 <- rep(1:3, each = 200)
X3 <- matrix(rnorm(600 * 2), 600, 2)
X3[lab3 == 2, 1] <- X3[lab3 == 2, 1] + 6      # far
X3[lab3 == 3, 1] <- X3[lab3 == 3, 1] + 1.5    # close to class 1
r3 <- robust_separation(X3, lab3)
ok("min is the hardest pair, not the average", r3$min < r3$mean,
   sprintf("min %.2f, mean %.2f", r3$min, r3$mean))

cat("\n")
if (fail > 0L) stop(sprintf("separation tests: %d FAILED", fail), call. = FALSE)
cat("separation tests: PASS\n")
