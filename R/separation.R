# Robust Mahalanobis separation between classes.
#
# `d` is a poor measure of how hard a scenario is. It sits on each variable's
# raw scale, ignores the within-class score dispersion Z, ignores how many
# variables carry the signal, and — decisively here — presumes a variance that
# arm 6 does not have.
#
# This measures separation where the algorithm actually works: in the estimated
# multivariate FPC score space, AFTER sparse observation and score estimation.
# With Poisson(7) points per subject the separation you can exploit is not the
# separation of the noiseless mean functions.
#
#   Delta^2 = t(rho_1 - rho_2) %*% solve(S) %*% (rho_1 - rho_2)
#
# Location is the spatial median per class and scatter is the pooled
# within-class robust estimate, both defined without moment assumptions, so the
# measure survives arms 5 and 6 where a covariance-based one does not.
#
# The scatter is built as C diag(nu) t(C) from robust_stage2() rather than from
# sscm() directly. The SSCM is the covariance of unit vectors and has trace 1,
# so inverting it would put Delta on an arbitrary scale; robust_stage2()'s `nu`
# recovers the spread along each direction on the score scale, which is what
# makes Delta comparable to a conventional Mahalanobis distance.
#
# Two uses:
#
#   1. A reference value. For two classes, Phi(Delta/2) is what the optimal
#      LINEAR rule achieves from centroid separation alone. It supplies a
#      conventional reference value, so a CCR of 0.98 can be read against a
#      yardstick rather than in a vacuum. It is NOT a bound here: see reference_ccr() below.
#
#   2. A check that the IQR-matching worked. The arms are matched on IQR so that
#      a difference cannot be attributed to signal-to-noise rather than tails.
#      Delta should therefore be near-constant ACROSS arms at fixed
#      d. If it is not, the arms are not exchangeable and tail effects are
#      confounded with separation, which is the objection the matching exists to
#      foreclose. This turns that assertion into evidence.

# Pairwise robust Mahalanobis distances between class centroids.
# Returns the full matrix plus the minimum, which is the binding pair: the
# hardest two classes to separate set the difficulty, not the average.
robust_separation <- function(scores, labels) {
  sm  <- MRKCFC:::spatial_median
  rs2 <- MRKCFC:::robust_stage2
  cls <- sort(unique(labels))
  if (length(cls) < 2L) return(list(matrix = NULL, min = NA_real_, mean = NA_real_))

  cent <- lapply(cls, function(k) sm(scores[labels == k, , drop = FALSE]))

  # pooled within-class scatter: centre each class on its own location first, so
  # the between-class shift does not inflate it
  resid <- do.call(rbind, lapply(seq_along(cls), function(i)
    sweep(scores[labels == cls[i], , drop = FALSE], 2L, cent[[i]], "-")))
  e <- rs2(resid, robust = TRUE)
  S <- e$vectors %*% diag(e$values, nrow = length(e$values)) %*% t(e$vectors)
  Sinv <- tryCatch(solve(S), error = function(err) MASS::ginv(S))

  M <- matrix(NA_real_, length(cls), length(cls), dimnames = list(cls, cls))
  for (i in seq_along(cls)) for (j in seq_along(cls)) {
    dv <- cent[[i]] - cent[[j]]
    M[i, j] <- sqrt(max(0, as.numeric(t(dv) %*% Sinv %*% dv)))
  }
  off <- M[upper.tri(M)]
  list(matrix = M, min = min(off), mean = mean(off))
}

# Reference CCR implied by a separation Delta, for two equal-prior classes.
#
# NOT a Bayes bound, and it was mislabelled as one initially. Phi(Delta/2) is
# the error rate of the optimal LINEAR rule for two Gaussians with a COMMON
# covariance, i.e. the best achievable from centroid separation alone.
#
# This design violates that assumption deliberately: variables 1 and 2 have
# class-specific eigenfunctions, sines for class 1 and cosines for class 2, so
# the classes differ in covariance structure as well as in location. The true
# Bayes rule there is quadratic and does better, and methods exploiting the
# shape difference legitimately EXCEED this reference. Observed at d = 1.5:
# realised CCR reached 107 to 111 percent of it on arms 1 and 4.
#
# So read it as a yardstick, not a ceiling. Exceeding it is informative — it
# says a method is using more than the centroid gap. Falling well short of it
# says the cell is a floor. Both readings are useful; treating it as a bound is
# not.
#
# Undefined for K > 2, so NA is returned rather than a number quietly wrong.
reference_ccr <- function(delta, K = 2L) {
  if (K != 2L) return(NA_real_)
  stats::pnorm(delta / 2)
}

# Retained under the old name so nothing silently breaks, but the name is
# misleading; prefer reference_ccr().
bayes_ccr <- reference_ccr

# Convenience: separation for one generated replicate, measured on the same
# global multivariate scores the score-then-cluster comparators consume.
separation_of <- function(dat, cfg, robust = TRUE) {
  sc <- tryCatch(global_mfpca_scores(dat$Ly, dat$Lt, cfg, robust),
                 error = function(e) NULL)
  if (is.null(sc)) return(list(delta = NA_real_, bayes = NA_real_))
  r <- robust_separation(sc, dat$labels)
  list(delta = r$min, bayes = reference_ccr(r$min, length(unique(dat$labels))))
}
