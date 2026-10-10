# Estimation accuracy in FUNCTION space, scored on the TRUE partition.
#
# ONE definition, sourced by both R/study_estim.R and (once merged) the sweep.
# Two copies of a scoring function is how this project acquired a divergent
# silhouette_score() that went unnoticed for weeks.
#
# WHY THIS STUDY EXISTS
#
# Every other study scores CLUSTERING, and clustering on the heavy-tailed arms is
# dominated by initialisation failure. Fitting under the TRUE partition removes
# the initialiser, the K-centers loop and convergence, leaving the question the
# method is actually about: does the spatial-sign scatter estimate the second
# stage better than the moment scatter?
#
# WHY FUNCTION SPACE AND NOT SCORE SPACE. Comparing the estimated rotation C
# against eigen(Z)$vectors is WRONG and a first version did it: C operates on
# ESTIMATED stacked scores whose components carry arbitrary sign, arbitrary
# scale, and an order set by each variable's own estimated eigenvalues, while the
# true rho blocks are ordered by the generating model. Unaligned bases gave
# 50-85 degrees on every arm with no pattern. EIGENFUNCTIONS have no such
# ambiguity: estimated and true both live on the time grid, and principal angles
# already handle sign and within-subspace ordering.
#
# THE SCORE COMPARISON IS NOT A SECOND OPINION ON STAGE TWO. Canonical
# correlation was adopted here because it is invariant to any nonsingular linear
# transform of either set, which was believed to sidestep the alignment problem.
# It does -- and it sidesteps the second stage along with it, because the
# multivariate scores ARE the stage-one scores post-multiplied by an invertible
# rotation. Verified 2026-09-03; see the note at the score_cc computation below.
# It is retained as a check on STAGE ONE and must not be reported as evidence
# about the scatter estimator.
#
# ============================================================================
# THE M_true / nK_fit SEPARATION, added 2026-09-02. READ THIS BEFORE EDITING.
# ============================================================================
#
# The previous version used ONE variable, `M <- cfg$global$M_true`, for three
# different jobs: the size of the true rotation, the block indexing of the FITTED
# eigenfunctions, and the score truncation. That was harmless only while
# nK_fit == M_true == 2.
#
# With nK_fit = 4 it is silently WRONG. `fit$multi_rotation` is
# (P*nK_fit) x (P*nK_fit) = 12 x 12 with blocks of size 4 -- rows 1:4, 5:8, 9:12.
# Indexing it with ((j-1)*2+1):(j*2) takes rows 1:2, 3:4, 5:6, all inside the
# first two variables' blocks. It does NOT error. It returns numbers, and they
# are meaningless.
#
# So the two sides are now sized independently:
#
#   FITTED  block j has ncol(phi_hat_j) columns -- read from the fit, NOT
#           assumed to be nK_fit, because on the non-robust branch fdapace
#           returns `selectK` which can differ between variables
#   TRUE    block j has M_true columns, from the generating design
#
# `keep` is fixed at 2, 3 and 4 on BOTH sides, and is NOT tied to nK_fit. The
# true subspace does not change when more components are retained; only the
# quality of its estimate does. That is precisely what makes angles comparable
# across nK_fit. Tying `keep` to nK_fit would make the TARGET move with the
# setting under test -- the same class of error as the M_true/nK_fit conflation
# fixed on 2026-09-02, which reversed the conclusion on the test cell.
#
# THREE TRUNCATIONS, because the leading-M target is ILL-CONDITIONED and that
# must be separated from the estimator's performance. Z's eigenvalues, and the
# relative gap BELOW each:
#
#   9.154   8.470   7.927  |  4.139   3.519   3.292
#        7.5%    6.4%   47.8%     15.0%    6.5%
#
# By Davis-Kahan the error in an estimated eigen-subspace scales as perturbation
# over eigen-gap, so the 47.8 percent drop after the third is the ONLY clean cut
# in this spectrum:
#
#   lead2  6.4% gap  -- what Step 3(d) actually uses; barely identifiable
#   lead3 47.8% gap  -- the well-posed cut, and the headline
#   lead4 15.0% gap  -- about 3x worse conditioned than lead3
#
# Winning at lead3 but not lead2 means the estimator is sound and the TARGET is
# unidentifiable; losing at both means the estimator is not helping.
#
# lead4 was added 2026-09-04 to answer the obvious question -- the study
# fits M = 4, so why is the headline scored at 3?
#
# MEASURED, AND IT DID NOT COME OUT AS PREDICTED. The expectation written here
# first was that lead4 would reproduce the lead3 ordering. It does not:
#
#   M2 (sign) vs M1 (moment)   12/12 at lead3 AND 12/12 at lead4  -- unchanged
#   three-way winner           changes on A07, A09, A10 (all elliptical), where
#                              M2's lead3 win over M3 reverses
#   every angle degrades       sign 14-16 deg -> 35-47, M3 10-47 -> 27-68
#
# The degradation is what a 15.0 percent gap predicts against 47.8. So the claim
# lead4 supports is the SIGN-AGAINST-MOMENT one, which is what this study exists
# to test; it does NOT support a truncation-invariant three-way ranking. Do not
# promote lead4 to the headline, and do not state stability it does not have.
#
# The full rotation is NOT a usable third comparison: C is invertible, so it
# spans the same subspace as Phi and the angle collapses to measuring stage one.
#
# WHAT THEORY PREDICTS, AND WHERE IT DOES NOT APPLY. The SSCM's eigenvectors
# coincide with the covariance's for ELLIPTICALLY symmetric distributions, while
# depending only on the directions of centred observations, so they stay defined
# when moments do not. Only four arms are elliptical:
#
#   arm 1 Gaussian            elliptical      expect a small LOSS (efficiency)
#   arm 2 Fleishman(0,4)      NOT elliptical  no guarantee
#   arm 3 Fleishman(1.5,4)    NOT elliptical  skewed; no guarantee
#   arm 4 contaminated normal elliptical      expect a WIN
#   arm 5 t(3)                elliptical      expect a WIN (no 4th moment)
#   arm 6 t(1.5)              elliptical      expect the LARGEST win
#
# Reporting the arms where the theory does not apply is the point, not an
# oversight: a robustness claim that holds only under ellipticity should say so.

# Measurement truncations, in the order the output columns take them:
# angle (lead2), angle_lead3, angle_lead4. See the eigen-gap table above.
ESTIM_KEEPS <- c(2L, 3L, 4L)

ESTIM_VARIANTS <- list(M1 = list(robust = TRUE,  scatter = "moment"),
                       M2 = list(robust = TRUE,  scatter = "sign"),
                       M3 = list(robust = FALSE, scatter = "moment"))

# Principal angles between two subspaces, in degrees; the largest is the
# conservative summary. On an equally spaced grid the Euclidean inner product is
# the L2 one up to a constant, and angles are invariant to that constant.
principal_angles <- function(A, B) {
  if (is.null(A) || is.null(B) || anyNA(A) || anyNA(B)) return(NA_real_)
  if (ncol(A) != ncol(B) || nrow(A) != nrow(B)) return(NA_real_)
  s <- svd(t(qr.Q(qr(A))) %*% qr.Q(qr(B)))$d
  max(acos(pmin(1, pmax(-1, s))) * 180 / pi)
}

# The fitted per-variable eigenfunctions and their grid. KFPCA and fdapace name
# these differently.
fitted_phi <- function(obj, robust) {
  if (robust) list(phi = obj$FPC_dis, grid = obj$RegGrid)
  else        list(phi = obj$phi,     grid = obj$workGrid)
}

# Stack the P elements into one multivariate eigenfunction basis, Happ and Greven
# Equation (9): psi_m^(j) = Phi_j %*% C[block_j, m].
#
# THE ELEMENTS MUST BE STACKED, not scored one at a time. A multivariate
# eigenfunction is a P-tuple and the rotation mixes information ACROSS variables
# -- that coupling is the entire purpose of the second stage. Scored per element,
# Phi_j %*% C[block_j, 1:keep] has an invertible factor and spans exactly the
# same subspace as Phi_j itself, so the measure is BLIND to the rotation: a first
# version returned angles identical to three decimals for the moment and sign
# variants in all twelve scenarios, because it was only ever measuring stage one.
#
# `blk` is the per-variable block size, passed in rather than assumed: the fitted
# side reads it from the fit, the true side uses M_true.
stack_psi <- function(getPhi, C, keep, blk) {
  off <- c(0L, cumsum(blk))
  if (nrow(C) < off[length(off)] || ncol(C) < keep) return(NULL)
  blocks <- lapply(seq_along(blk), function(j) {
    ph <- getPhi(j)
    if (is.null(ph) || ncol(ph) < blk[j]) return(NULL)
    ph[, seq_len(blk[j]), drop = FALSE] %*%
      C[(off[j] + 1L):off[j + 1L], seq_len(keep), drop = FALSE]
  })
  if (any(vapply(blocks, is.null, TRUE))) return(NULL)
  if (length(unique(vapply(blocks, nrow, 1L))) != 1L) return(NULL)
  do.call(rbind, blocks)
}

# One replicate: K classes x 3 variants, fitted on the TRUE partition.
#
#   dat  from generate_replicate()
#   cfg  the config, whose nK_fit governs the FIT and whose M_true governs the
#        generating truth. They are no longer assumed equal.
#   d    the scenario's separation, needed to rebuild the true basis
oracle_estimation <- function(dat, cfg, d) {
  K  <- as.integer(cfg$global$K_true)
  Mt <- as.integer(cfg$global$M_true)      # the TRUTH: 2
  P  <- as.integer(cfg$global$P)

  Ztrue <- cfg$Z[seq_len(P * Mt), seq_len(P * Mt), drop = FALSE]
  Ctrue <- eigen(Ztrue, symmetric = TRUE)$vectors
  blk_true <- rep(Mt, P)

  do.call(rbind, lapply(names(ESTIM_VARIANTS), function(vid) {
    v <- ESTIM_VARIANTS[[vid]]
    # Common random numbers, and the SAME key run_all_methods() uses, so the
    # oracle fits are paired with the clustering fits on this replicate.
    c_seed_fit(cfg, dat$scenario_id, dat$rep)
    optns <- build_optns(cfg, v$robust)

    do.call(rbind, lapply(seq_len(K), function(k) {
      bad <- function(msg) data.frame(
        variant = vid, class = k, angle = NA_real_, angle_lead3 = NA_real_,
        angle_lead4 = NA_real_,
        score_cc = NA_real_, nK_used = NA_integer_, error = msg,
        stringsAsFactors = FALSE)

      fit <- tryCatch(MRKCFC::MFPCA_k(k, dat$Ly, dat$Lt, dat$labels, optns,
                                     v$robust, mfpca_scatter = v$scatter),
                      error = function(e) e)
      if (inherits(fit, "error")) return(bad(conditionMessage(fit)))

      Chat    <- fit$multi_rotation
      ph_hat  <- function(j) fitted_phi(fit$ref_hatyc[[j]], v$robust)$phi
      # READ the block sizes from the fit. fdapace's selectK can differ between
      # variables, so a hard-coded stride of nK_fit would misalign the blocks --
      # silently, since the matrix arithmetic still conforms.
      blk_hat <- vapply(seq_len(P), function(j) {
        ph <- ph_hat(j); if (is.null(ph)) NA_integer_ else ncol(ph)
      }, integer(1))
      if (anyNA(blk_hat)) return(bad("a per-variable fit returned no eigenfunctions"))
      if (sum(blk_hat) != nrow(Chat))
        return(bad(sprintf("block sizes %s do not sum to rotation dim %d",
                           paste(blk_hat, collapse = "+"), nrow(Chat))))

      pg  <- fitted_phi(fit$ref_hatyc[[1]], v$robust)
      # TRUE eigenfunctions on the SAME grid the fit used, via the same builder
      # the generator used, so the corrections recorded in basis.csv (variable
      # 2's frequencies, variable 3's Gram-Schmidt) are applied here exactly as
      # they were when the data was made.
      bas <- build_basis(cfg, pg$grid, d = as.numeric(d))
      ph_true <- function(j) bas$psi[[j]][[k]]

      ang <- vapply(ESTIM_KEEPS, function(keep) principal_angles(
        stack_psi(ph_hat,  Chat,  keep, blk_hat),
        stack_psi(ph_true, Ctrue, keep, blk_true)), numeric(1))

      # SCORE CANONICAL CORRELATION -- COMPUTED, BUT NOT EVIDENCE ABOUT STAGE
      # TWO. DO NOT REPORT IT AS SUCH. Verified 2026-09-03.
      #
      # multi_score = xi %*% C with C invertible, and canonical correlation is
      # invariant to any nonsingular linear map -- so at full rank it is
      # invariant to THE ROTATION, which is the quantity under test. The M1 and
      # M2 variants share a stage one and differ only in the stage-two scatter;
      # their full-rank canonical correlations agree to three decimals (0.982
      # on A01, 0.970 on A11). A measure that cannot separate them cannot be
      # evidence about them. What it measures is how well the STAGE-ONE scores
      # predict the truth.
      #
      # THE SAME ERROR IN A DIFFERENT GUISE: an earlier angle measurement scored
      # per element, where Phi_j %*% C[idx_j, ] spans the same subspace as Phi_j
      # and the rotation was invisible. General rule: a quantity related to the
      # estimate by an invertible transform cannot measure that transform. Check
      # invariance BEFORE adopting a metric.
      #
      # Kept because it is a valid check on STAGE ONE and costs nothing.
      # Both sides take M_true columns; giving the fitted side more would
      # mechanically inflate it.
      # of either set. BOTH SIDES GET M_true COLUMNS -- giving the fitted side
      # more would mechanically inflate the correlation (more directions to match
      # against) and manufacture an nK_fit effect that is pure arithmetic.
      idxk <- which(dat$labels == k)
      rho_true <- dat$scores[idxk, , drop = FALSE] %*% Ctrue
      cc <- tryCatch({
        sh <- fit$multi_score[, seq_len(Mt), drop = FALSE]
        mean(stats::cancor(sh, rho_true[, seq_len(Mt), drop = FALSE])$cor)
      }, error = function(e) NA_real_)

      data.frame(variant = vid, class = k, angle = ang[1], angle_lead3 = ang[2],
                 angle_lead4 = ang[3],
                 score_cc = cc, nK_used = as.integer(blk_hat[1]),
                 error = NA_character_, stringsAsFactors = FALSE)
    }))
  }))
}
