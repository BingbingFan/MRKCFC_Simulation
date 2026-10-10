# Metrics for Study A: CCR, aRI, NMI, AMI, MPE, RRE.
#
# Two families, and they fail in different places, which is why all of them are
# carried:
#
#   CCR   accuracy after optimal label matching. Needs a bijection, so it is
#         only defined when K-hat equals K. Study A fixes K at the truth, so it
#         is always defined here; it is not in Study C3.
#   aRI   pair-counting agreement, chance-corrected.
#   NMI   information-theoretic. NOT chance-corrected: it rewards over-splitting,
#         so a method returning too many clusters is flattered. Read next to aRI.
#   AMI   NMI with aRI's chance correction applied. Added because Study C3 lets
#         K-hat vary, where plain NMI would mislead.
#
#   MPE   manuscript Equation (21), mean squared reconstruction error.
#   RRE   the revision's replacement, per-variable relative error, unit-free.
#
# MPE and RRE are EVALUATION metrics. Equation (15), the reassignment rule, is
# a different object living in clust_update()'s seCost. They are both sums of
# squared residuals, which is exactly why they are easy to conflate. Nothing
# here touches the assignment rule.

# --- clustering agreement ----------------------------------------------------

# All permutations of seq_len(k). K is small here (2 in Study A, up to 5 in
# Study B), so exhaustive matching is cheaper and clearer than Hungarian.
.perms <- function(k) {
  if (k == 1L) return(matrix(1L, 1L, 1L))
  sub <- .perms(k - 1L)
  do.call(rbind, lapply(seq_len(k), function(i)
    cbind(i, matrix(c(seq_len(k)[-i])[sub], nrow(sub)))))
}

# Correct classification rate under the best label permutation.
ccr <- function(truth, est) {
  truth <- as.integer(factor(truth)); est <- as.integer(factor(est))
  k <- max(length(unique(truth)), length(unique(est)))
  if (k > 7L) stop("ccr: too many clusters for exhaustive matching", call. = FALSE)
  P <- .perms(k)
  max(apply(P, 1L, function(p) mean(p[est] == truth, na.rm = TRUE)))
}

ari <- function(truth, est) mclust::adjustedRandIndex(truth, est)

# Shared machinery for NMI and AMI.
.mi_parts <- function(truth, est) {
  tab <- table(truth, est)
  N <- sum(tab)
  a <- rowSums(tab); b <- colSums(tab)
  pij <- tab / N; pi_ <- a / N; pj_ <- b / N
  nz <- tab > 0
  mi <- sum(pij[nz] * log(pij[nz] / outer(pi_, pj_)[nz]))
  hu <- -sum(pi_[pi_ > 0] * log(pi_[pi_ > 0]))
  hv <- -sum(pj_[pj_ > 0] * log(pj_[pj_ > 0]))
  list(tab = tab, N = N, a = a, b = b, mi = mi, hu = hu, hv = hv)
}

# Normalized mutual information, arithmetic normalizer (scikit-learn's default).
# Recorded in study/config; state it in Methods, since geometric and max
# normalizers give values a few percent apart.
nmi <- function(truth, est) {
  p <- .mi_parts(truth, est)
  den <- (p$hu + p$hv) / 2
  if (den <= 0) return(if (p$mi <= 0) 1 else 0)
  p$mi / den
}

# Adjusted mutual information (Vinh, Epps & Bailey 2010). The expectation is
# taken over the hypergeometric model of random partitions with the observed
# margins fixed, the same null aRI uses.
ami <- function(truth, est) {
  p <- .mi_parts(truth, est)
  N <- p$N; a <- p$a; b <- p$b
  emi <- 0
  for (i in seq_along(a)) for (j in seq_along(b)) {
    lo <- max(1L, a[i] + b[j] - N); hi <- min(a[i], b[j])
    if (hi < lo) next
    nij <- lo:hi
    lterm <- lgamma(a[i] + 1) + lgamma(b[j] + 1) +
             lgamma(N - a[i] + 1) + lgamma(N - b[j] + 1) -
             lgamma(N + 1) - lgamma(nij + 1) - lgamma(a[i] - nij + 1) -
             lgamma(b[j] - nij + 1) - lgamma(N - a[i] - b[j] + nij + 1)
    emi <- emi + sum((nij / N) * log((N * nij) / (a[i] * b[j])) * exp(lterm))
  }
  den <- (p$hu + p$hv) / 2 - emi
  if (abs(den) < .Machine$double.eps) return(0)
  (p$mi - emi) / den
}

# --- reconstruction error ----------------------------------------------------

# Manuscript Equation (21): the mean squared error over every observation of
# every variable of every subject.
mpe <- function(y, yhat) mean((y - yhat)^2)

# The revision's replacement. Per-variable relativization makes it unit-free
# and comparable between simulation and application; 0 is perfect and 1 is no
# better than predicting the variable's own mean. Trimming is on how extreme
# the OBSERVATION is, never on how large the error is, so competing methods are
# scored on the same retained set.
rre_overall <- function(y, yhat, variable, trim = 0) {
  unname(MRKCFC:::rre(y, yhat, variable, trim = trim)["overall"])
}

# --- assembly ----------------------------------------------------------------

# One row of results. yhat may be NULL for a method that produces no trajectory
# fit, in which case the reconstruction columns are NA rather than a fabricated
# zero.
score_fit <- function(truth, est, y = NULL, yhat = NULL, variable = NULL) {
  out <- list(
    CCR = ccr(truth, est),
    aRI = ari(truth, est),
    NMI = nmi(truth, est),
    AMI = ami(truth, est),
    MPE = NA_real_, RRE = NA_real_, RRE_trim = NA_real_
  )
  if (!is.null(yhat)) {
    out$MPE      <- mpe(y, yhat)
    out$RRE      <- rre_overall(y, yhat, variable, trim = 0)
    out$RRE_trim <- rre_overall(y, yhat, variable, trim = 0.01)
  }
  out
}
