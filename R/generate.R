# Data generator for Studies A, B and C.
#
# Implements manuscript Equation (19):
#
#   Y^(p)_ij,(k) = mu^(p)_(k)(t_ij) + sum_{m=1}^{M} rho_i,(p,m) psi^(p)_m(k)(t_ij) + eps
#
# Note on the score index. Equation (19) writes the score as rho_im, indexed by
# component only, but Z is 6 x 6 = P x M_true, so rho_i is a SIX-vector and
# variable p uses its own block, elements [(p-1)*M + 1, (p-1)*M + M]. The 6 x 6
# Z is the only reading consistent with the published matrix; the manuscript's
# subscript is loose. Z's block structure confirms it: within-variable
# diagonals (9,4), (8,3.5), (8.5,3.5) and cross-variable covariance 0.2.
#
# Everything is read from study/config. Nothing here is hard-coded.

source(file.path("config", "validate_config.R"))   # load_config, eval_basis, orthonormalize_pair

# Per-cell seed. Derived from the master seed plus the cell identity, so any
# scenario x replicate reproduces in isolation regardless of execution order or
# which tasks are re-run after a failure. This is what makes an HPC array
# restartable without replaying the whole grid.
cell_seed <- function(seed_master, scenario_id, rep) {
  key <- paste(scenario_id, rep, sep = "|")
  # character codes weighted by position: cheap, deterministic, no dependency
  h <- sum(utf8ToInt(key) * seq_along(utf8ToInt(key)) * 7919)
  as.integer((as.numeric(seed_master) + h) %% 2147483647)
}

# Build the basis on a grid: mean functions and eigenfunctions per variable and
# class, applying the recorded orthonormalize flag. Returns lists indexed
# [[variable]][[class]].
build_basis <- function(cfg, tgrid, d) {
  P <- as.integer(cfg$global$P)
  K <- as.integer(cfg$global$K_true)
  b <- cfg$basis
  mu <- psi <- vector("list", P)
  for (p in seq_len(P)) {
    mu[[p]] <- psi[[p]] <- vector("list", K)
    for (k in seq_len(K)) {
      mrow <- b[b$variable == p & b$class == k & b$role == "mean", ]
      stopifnot(nrow(mrow) == 1L)
      mu[[p]][[k]] <- eval_basis(mrow$expr, tgrid, d = d)

      erows <- b[b$variable == p & b$class == k & b$role == "eigen", ]
      erows <- erows[order(erows$component), ]
      f <- lapply(erows$expr, eval_basis, t = tgrid, d = d)
      psi[[p]][[k]] <- if (any(erows$orthonormalize)) {
        orthonormalize_pair(f[[1]], f[[2]])
      } else {
        do.call(cbind, f)
      }
    }
  }
  list(mu = mu, psi = psi)
}

# Observation counts. Poisson(lambda) truncated below at obs_min: a plain draw
# gives n_i < 2 for about 0.7% of subjects, which is degenerate and which
# predicted_yhat() cannot handle. Redrawn rather than clamped, so the retained
# distribution is a genuine truncated Poisson rather than a spike at the floor.
draw_counts <- function(n, lambda, obs_min) {
  out <- stats::rpois(n, lambda)
  bad <- which(out < obs_min)
  while (length(bad)) {
    out[bad] <- stats::rpois(length(bad), lambda)
    bad <- which(out < obs_min)
  }
  out
}

# --- Study B configuration helpers -------------------------------------------
# Study B moves parameters that live in params_global.csv rather than in the
# scenario table, so each row of scenarios_b.csv carries its own values. These
# live here, not in run_cell_b.R, because test_study_b.R needs them too and a
# second copy would drift from the first -- which is exactly how
# silhouette_score() acquired a divergent twin.

# Apply one scenarios_b.csv row's overrides to a COPY of the config. Anything
# not named here keeps its global value, which is what makes Study B
# one-factor-at-a-time rather than factorial.
apply_config_b <- function(cfg, row) {
  g <- cfg$global
  g$n                  <- row$n
  g$P                  <- row$P
  g$pi                 <- row$pi
  g$obs_poisson_lambda <- row$obs_poisson_lambda
  # obs_min is per-configuration so the SPARSITY AXIS can be decided in the
  # table rather than in code. It is 4 throughout as of 2026-09-02, matching the
  # global config, but the axis is the one place a different value could be
  # argued for: at lambda = 3, obs_min = 4 redraws 65 percent of subjects and
  # the realised mean becomes 4.91, so that level is a truncated distribution
  # rather than Poisson(3). It stays at 4 because the manuscript states a
  # realised range of 4 to 12 observations, and obs_min = 2 would generate
  # subjects with two - contradicting the paper's own description of its design.
  # `realised_mu` in scenarios_b.csv records what each level ACTUALLY varies.
  #
  # Absent in the pre-2026-09-02 table, so fall back rather than fail: an old
  # scenarios_b.csv then reproduces its original behaviour instead of erroring.
  if (!is.null(row$obs_min) && !is.na(row$obs_min)) g$obs_min <- row$obs_min
  cfg$global <- g
  cfg
}

# generate_replicate() reads the scenario's identity and arm parameters from a
# one-row frame, so a Study B configuration has to be presented in the same
# column names the Study A table uses.
as_scenario_b <- function(row) data.frame(
  scenario_id = row$config_id, arm = row$arm, d = row$d,
  arm_eps = row$arm_eps, arm_k = row$arm_k, arm_nu = row$arm_nu,
  stringsAsFactors = FALSE)

# Generate one replicate of one scenario.
#
# Returns Ly and Lt in the list-of-lists shape MRKCFC expects, the true class
# labels, and the normalization attributes if normalization was applied.
generate_replicate <- function(cfg, scenario, rep, arms_env = globalenv()) {
  g   <- cfg$global
  n   <- as.integer(g$n); P <- as.integer(g$P)
  K   <- as.integer(g$K_true); M <- as.integer(g$M_true)
  pi_ <- as.numeric(strsplit(g$pi, ",")[[1]])
  stopifnot(length(pi_) == K, abs(sum(pi_) - 1) < 1e-8)

  set.seed(cell_seed(g$seed_master, scenario$scenario_id, rep))

  tgrid <- seq(as.numeric(g$t_min), as.numeric(g$t_max),
               length.out = as.integer(g$grid_n))
  bas <- build_basis(cfg, tgrid, d = as.numeric(scenario$d))

  # class membership
  labels <- sample.int(K, n, replace = TRUE, prob = pi_)

  # multivariate FPC scores, one (P * M)-vector per subject, from the arm
  # generator. Z is stored at its full 20 x 20 for Study B's P axis and
  # subsetted to the leading P * M here, so there is ONE Z in the project rather
  # than one per P. The leading 6 x 6 is bit-identical to the published matrix,
  # so Study A at P = 3 is unaffected by the extension.
  stopifnot(P * M <= nrow(cfg$Z))
  Zp <- cfg$Z[seq_len(P * M), seq_len(P * M), drop = FALSE]
  gen <- get("gen_scores", envir = arms_env)
  rho <- gen(n, scenario$arm, Zp,
             eps = if (is.na(scenario$arm_eps)) 0.05 else scenario$arm_eps,
             k   = if (is.na(scenario$arm_k))   5    else scenario$arm_k)
  stopifnot(ncol(rho) == P * M)

  # observation times, shared across variables for a subject: the manuscript
  # samples one irregular schedule per subject, not one per variable
  n_i <- draw_counts(n, as.numeric(g$obs_poisson_lambda), as.integer(g$obs_min))
  idx <- lapply(n_i, function(m) sort(sample.int(length(tgrid), m)))

  sd_eps <- as.numeric(g$error_sd)
  Ly <- Lt <- vector("list", P)
  for (p in seq_len(P)) {
    Ly[[p]] <- Lt[[p]] <- vector("list", n)
    block <- ((p - 1L) * M + 1L):(p * M)
    for (i in seq_len(n)) {
      k <- labels[i]; ii <- idx[[i]]
      fit <- bas$mu[[p]][[k]][ii] +
             as.vector(bas$psi[[p]][[k]][ii, , drop = FALSE] %*% rho[i, block])
      Ly[[p]][[i]] <- fit + stats::rnorm(length(ii), 0, sd_eps)
      Lt[[p]][[i]] <- tgrid[ii]
    }
  }

  norm_attr <- NULL
  if (isTRUE(as.logical(g$normalize))) {
    norm_attr <- vector("list", P)
    for (p in seq_len(P)) {
      pooled <- unlist(Ly[[p]], use.names = FALSE)
      z <- MRKCFC:::robust_normalize(pooled)
      norm_attr[[p]] <- list(centre = attr(z, "centre"), scale = attr(z, "scale"))
      Ly[[p]] <- lapply(Ly[[p]], function(y)
        (y - norm_attr[[p]]$centre) / norm_attr[[p]]$scale)
    }
  }

  list(Ly = Ly, Lt = Lt, labels = labels, n_i = n_i,
       scores = rho, normalization = norm_attr,
       scenario_id = scenario$scenario_id, rep = rep)
}
