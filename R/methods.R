# The eight Study A methods, plus the shared machinery they need.
#
# Every method returns the same shape: a partition, a convergence flag, and an
# error string if it failed. Failures are RECORDED, never fatal: one method
# blowing up on one replicate must not take down an array task.
#
# On MPE and RRE. Four of the eight methods cluster scores and produce no
# trajectory fit at all, so scoring each method by its own reconstruction would
# not be a like-for-like comparison. Instead every partition is scored by ONE
# common reconstruction, set by recon_robust and recon_scatter in the config.
# MPE and RRE therefore measure the PARTITION, not the fitting apparatus that
# produced it. The manuscript states this; it is a deliberate choice, since
# otherwise a reader would wonder how k-means acquired a reconstruction error.

suppressMessages({
  library(MRKCFC); library(mclust)
})

# --- shared machinery --------------------------------------------------------

build_optns <- function(cfg, robust) {
  g <- cfg$global
  interval <- c(as.numeric(g$t_min), as.numeric(g$t_max))
  if (robust) {
    list(
      interval  = interval,
      nK        = as.integer(g$nK_fit),
      bw        = as.numeric(g$bw),
      bwK       = g$bwK,
      bwmean    = g$bwmean,
      nRegGrid  = as.integer(g$nRegGrid),
      basis     = fda::create.bspline.basis(
        interval, nbasis = as.integer(g$basis_nbasis),
        norder = as.integer(g$basis_norder),
        breaks = seq(interval[1], interval[2],
                     length.out = as.integer(g$basis_nbreaks)))
    )
  } else {
    list(dataType = "Sparse", methodMuCovEst = "smooth",
         nRegGrid = as.integer(g$nRegGrid),
         methodSelectK = as.integer(g$nK_fit),
         useBinnedData = "OFF")
  }
}

# Multivariate FPC scores from ONE global MFPCA over all subjects, which is
# what the score-then-cluster comparators consume. Note this is exactly the
# handicap to state in Methods: the design gives variables 1 and 2
# class-specific eigenfunctions, so a single global basis is wrong for both
# classes. MKCFC and MRKCFC re-estimate class-specific FPCA inside their loop.
global_mfpca_scores <- function(Ly, Lt, cfg, robust) {
  n <- length(Ly[[1]])
  fit <- MRKCFC::MFPCA_k(1L, Ly, Lt, rep(1L, n), build_optns(cfg, robust),
                        robust, mfpca_scatter = "moment")
  fit$multi_score
}

# Fitted values for every observation, given any partition. This is manuscript
# Equation (14) applied class by class, and it is the only place MPE and RRE
# get their yhat.
reconstruct_given_partition <- function(Ly, Lt, labels, cfg) {
  g <- cfg$global
  robust  <- isTRUE(as.logical(g$recon_robust))
  scatter <- g$recon_scatter
  optns   <- build_optns(cfg, robust)
  P <- length(Ly); n <- length(Ly[[1]])
  K <- length(unique(labels))

  ref <- MRKCFC:::.mrkcfc_fit_multivariate_reference(Ly, Lt, K, labels, optns, robust)
  yhat <- vector("list", P)
  for (p in seq_len(P)) yhat[[p]] <- vector("list", n)
  for (i in seq_len(n)) {
    ith_yy <- lapply(Ly, function(z) z[[i]])
    ith_tt <- lapply(Lt, function(z) z[[i]])
    pr <- MRKCFC::predicted_yhat(ref[[labels[i]]], ith_yy, ith_tt,
                                robust, optns, scatter)
    for (p in seq_len(P)) yhat[[p]][[i]] <- pr$fittedX[[p]]
  }
  list(
    y        = unlist(lapply(Ly,   function(z) unlist(z, use.names = FALSE))),
    yhat     = unlist(lapply(yhat, function(z) unlist(z, use.names = FALSE))),
    variable = rep(seq_len(P), times = vapply(Ly, function(z)
      sum(lengths(z)), 1L))
  )
}

# lcmm wants one row per subject-time with all outcomes side by side. The
# design gives each subject ONE schedule shared across variables, so this is a
# clean reshape rather than a merge.
to_long <- function(Ly, Lt) {
  n <- length(Ly[[1]]); P <- length(Ly)
  do.call(rbind, lapply(seq_len(n), function(i) {
    d <- data.frame(id = i, time = Lt[[1]][[i]])
    for (p in seq_len(P)) d[[paste0("Y", p)]] <- Ly[[p]][[i]]
    d
  }))
}

# --- the eight methods -------------------------------------------------------
# Each takes (Ly, Lt, K, cfg) and returns list(labels, converged).

m_mkcfc <- function(robust, scatter) function(Ly, Lt, K, cfg) {
  g <- cfg$global
  fit <- MRKCFC::MKCFC(
    Ly = Ly, Lt = Lt, K = K, robust = robust,
    optns = build_optns(cfg, robust), mfpca_scatter = scatter,
    init.method = g$init_method, maxIter = as.integer(g$maxIter),
    parallel = FALSE, trace = 0,
    conv.prop = as.numeric(g$conv_prop),
    LOO = isTRUE(as.logical(g$LOO)))
  # `iter` separates the two ways conv can be FALSE, which need different fixes:
  # the loop exhausted maxIter, or it broke early because a group fell below ten
  # curves. Without it a non-convergence rate is not actionable.
  list(labels = as.integer(fit$res_clust$g), converged = isTRUE(fit$conv),
       iter = fit$iter)
}

m_scores_then <- function(robust, cluster_fn) function(Ly, Lt, K, cfg) {
  sc <- global_mfpca_scores(Ly, Lt, cfg, robust)
  list(labels = cluster_fn(sc, K), converged = TRUE)
}

cl_mclust <- function(sc, K) as.integer(mclust::Mclust(sc, G = K, verbose = FALSE)$classification)
cl_kmeans <- function(sc, K) as.integer(stats::kmeans(sc, centers = K, nstart = 25)$cluster)
cl_teigen <- function(sc, K) {
  fit <- teigen::teigen(sc, Gs = K, verbose = FALSE, scale = FALSE)
  # teigen returns an object with no classification when nothing converged.
  # Fail loudly here so the caller records it, rather than returning NULL and
  # having it surface later as an unrelated length error.
  if (is.null(fit$classification))
    stop("teigen: no model converged", call. = FALSE)
  as.integer(fit$classification)
}

# multlcmm does NOT converge on this design, and that is expected rather than a
# specification error. It models multiple outcomes measuring ONE common latent
# process; this design has three variables with different eigenfunctions,
# different Z blocks and only 0.2 cross-covariance, so three distinct processes.
# Raising the cap does not help: 100, 300 and 1000 iterations all end at
# conv = 2 with the log-likelihood moving only -8389.66 to -8388.32 across 700
# extra iterations. Four time specifications, linear through natural splines
# with a random spline basis, failed identically at the one-class stage.
#
# It is kept in the comparison and its non-convergence is REPORTED, since
# latent class mixed models are the natural model-based comparator and this is
# an honest boundary-of-applicability result. The `converged` column carries it; read the
# clustering metrics for this method only alongside that column.
m_multlcmm <- function(Ly, Lt, K, cfg) {
  df <- to_long(Ly, Lt)
  P <- length(Ly)
  mx <- as.integer(cfg$global$lcmm_maxiter)
  lhs <- paste(paste0("Y", seq_len(P)), collapse = " + ")
  f <- stats::as.formula(paste(lhs, "~ time"))
  # A one-class fit supplies starting values. It is kept even though it too
  # fails to converge: a cold start at ng > 1 is worse practice.
  m1 <- lcmm::multlcmm(f, random = ~ time, subject = "id", ng = 1,
                       data = df, link = "linear", maxiter = mx, verbose = FALSE)
  if (K == 1L) return(list(labels = rep(1L, length(Ly[[1]])), converged = m1$conv == 1))
  mk <- lcmm::multlcmm(f, random = ~ time, subject = "id", ng = K,
                       mixture = ~ time, data = df, link = "linear",
                       B = m1, maxiter = mx, verbose = FALSE)
  pp <- mk$pprob[order(mk$pprob$id), ]
  list(labels = as.integer(pp$class), converged = mk$conv == 1)
}

METHODS <- list(
  M1 = m_mkcfc(robust = TRUE,  scatter = "moment"),
  M2 = m_mkcfc(robust = TRUE,  scatter = "sign"),
  M3 = m_mkcfc(robust = FALSE, scatter = "moment"),
  M4 = m_scores_then(robust = FALSE, cl_mclust),
  M5 = m_scores_then(robust = FALSE, cl_kmeans),
  M6 = m_scores_then(robust = TRUE,  cl_kmeans),
  M7 = m_multlcmm,
  M8 = m_scores_then(robust = FALSE, cl_teigen)
)

# --- driver for one replicate ------------------------------------------------

run_all_methods <- function(dat, cfg, methods = names(METHODS)) {
  K <- as.integer(cfg$global$K_true)
  n <- length(dat$Ly[[1]])

  blank <- function(mid, elapsed, msg) {
    data.frame(scenario_id = dat$scenario_id, rep = dat$rep, method = mid,
               CCR = NA_real_, aRI = NA_real_, NMI = NA_real_, AMI = NA_real_,
               MPE = NA_real_, RRE = NA_real_, RRE_trim = NA_real_,
               converged = NA, iter = NA_integer_, seconds = elapsed, error = msg,
               stringsAsFactors = FALSE)
  }

  rows <- lapply(methods, function(mid) {
    # COMMON RANDOM NUMBERS. Every method starts its fit from the identical RNG
    # state rather than from whatever the previous method happened to leave
    # behind. Without this, a method's result depends on which methods ran
    # before it, so (a) toggling methods.csv's `include` silently moves every
    # method after the one toggled, (b) a recorded (seed_master, scenario, rep)
    # is not enough to reproduce a fit -- you also need the exact method list
    # and its order -- and (c) a standalone diagnostic cannot reproduce a
    # driver run. All three were observed: c4mini reproduced replicate 1 of the
    # 10-replicate round on t nu=1.5 d=3 but not replicates 2 and 3, because it
    # ran M2 alone where run_cell.R runs M1 first, and "two-staged"
    # initialisation draws on the RNG twice (k-means centres, and the random
    # fallback when k-means splits degenerately).
    #
    # The same seed for every method rather than one seed per method: the
    # methods are being COMPARED, and starting them from a common state pairs
    # them, which reduces the variance of the method-to-method differences the
    # paper reports. The cost is that a bad draw hits every method on that
    # replicate together, which is a real property of that replicate rather
    # than an artefact.
    #
    # The "|fit" tag keeps this stream distinct from the one generate_replicate()
    # used for the data itself, so the algorithm's randomness is not drawn from
    # the same stream that produced the sample it is fitting.
    set.seed(cell_seed(cfg$global$seed_master,
                       paste0(dat$scenario_id, "|fit"), dat$rep))
    t0 <- Sys.time()
    res <- tryCatch(METHODS[[mid]](dat$Ly, dat$Lt, K, cfg),
                    error = function(e) structure(list(), error = conditionMessage(e)))
    el <- function() as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!is.null(attr(res, "error"))) return(blank(mid, el(), attr(res, "error")))

    # Validate before scoring. A method can return without erroring and still
    # hand back a partition of the wrong length or with NAs; scoring that dies
    # inside the metric functions, which would otherwise take down the whole
    # cell rather than just this method.
    lab <- res$labels
    if (is.null(lab) || length(lab) != n || anyNA(lab))
      return(blank(mid, el(), sprintf(
        "malformed partition: length %s for n = %d, %d NA",
        if (is.null(lab)) "NULL" else length(lab), n, sum(is.na(lab)))))

    # Scoring is wrapped too. Everything after the fit must be non-fatal, or a
    # single bad cell kills an entire array task hours into a run.
    out <- tryCatch({
      rec <- tryCatch(reconstruct_given_partition(dat$Ly, dat$Lt, lab, cfg),
                      error = function(e) NULL)
      sc <- if (is.null(rec)) score_fit(dat$labels, lab)
            else score_fit(dat$labels, lab, rec$y, rec$yhat, rec$variable)
      data.frame(scenario_id = dat$scenario_id, rep = dat$rep, method = mid,
                 CCR = sc$CCR, aRI = sc$aRI, NMI = sc$NMI, AMI = sc$AMI,
                 MPE = sc$MPE, RRE = sc$RRE, RRE_trim = sc$RRE_trim,
                 converged = res$converged,
               iter = if (is.null(res$iter)) NA_integer_ else as.integer(res$iter),
               seconds = el(),
                 error = NA_character_, stringsAsFactors = FALSE)
    }, error = function(e) blank(mid, el(), paste("scoring failed:", conditionMessage(e))))
    out
  })
  do.call(rbind, rows)
}
