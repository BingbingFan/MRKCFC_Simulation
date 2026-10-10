# Reads the configuration tables and checks they are internally consistent.
# Run from the study/ directory:
#
#   Rscript config/validate_config.R
#
# This is the only place the tables are interpreted. Simulation code should
# call load_config() rather than hard-coding any of these values, so that
# changing a parameter never means editing code.

CFG <- "config"

load_config <- function(dir = CFG) {
  g <- utils::read.csv(file.path(dir, "params_global.csv"), stringsAsFactors = FALSE)
  vals <- stats::setNames(as.list(g$value), g$param)
  list(
    global    = vals,
    meta      = g,
    scenarios = utils::read.csv(file.path(dir, "scenarios.csv"), stringsAsFactors = FALSE),
    methods   = utils::read.csv(file.path(dir, "methods.csv"), stringsAsFactors = FALSE),
    basis     = utils::read.csv(file.path(dir, "basis.csv"), stringsAsFactors = FALSE),
    Z         = local({
      z <- utils::read.csv(file.path(dir, "Z.csv"), stringsAsFactors = FALSE)
      m <- as.matrix(z[, -1]); rownames(m) <- z$row; m
    })
  )
}

# Evaluate a basis expression on a time grid. d is only used by mean functions.
eval_basis <- function(expr, t, d = 0) eval(parse(text = expr), list(t = t, d = d))

# Gram-Schmidt a pair of functions sampled on a grid, returning unit-norm,
# mutually orthogonal columns. Integration is the grid mean over [0,1].
orthonormalize_pair <- function(f1, f2) {
  ip <- function(a, b) mean(a * b)
  f1 <- f1 / sqrt(ip(f1, f1))
  f2 <- f2 - ip(f2, f1) * f1
  f2 <- f2 / sqrt(ip(f2, f2))
  cbind(f1, f2)
}

if (sys.nframe() == 0L) {
  cfg <- load_config()
  fail <- 0L
  ok <- function(lab, cond, extra = "") {
    good <- isTRUE(cond); if (!good) fail <<- fail + 1L
    cat(sprintf("%-54s %s %s\n", lab, if (good) "PASS" else "**FAIL**", extra))
  }

  cat("=== shape ===\n")
  ok("12 scenarios", nrow(cfg$scenarios) == 12L)
  ok("8 methods", nrow(cfg$methods) == 8L)
  P <- as.integer(cfg$global$P); M <- as.integer(cfg$global$M_true)
  P_MAX <- max(cfg$basis$variable)
  # Z is stored at its FULL size for Study B's P axis and subsetted to the
  # leading P * M at generation time, so it is no longer square at P * M_true.
  # The leading block must still reproduce the published 6 x 6 exactly, which is
  # what keeps Study A unaffected by the extension.
  ok(sprintf("Z is %dx%d", P_MAX * M, P_MAX * M),
     all(dim(cfg$Z) == c(P_MAX * M, P_MAX * M)))
  ok("Z symmetric", isTRUE(all.equal(cfg$Z, t(cfg$Z), check.attributes = FALSE)))
  ok("Z covers P * M_true", nrow(cfg$Z) >= P * M, sprintf("(%d x %d)", P, M))
  ok("every P in {3,5,10} gives a positive definite subset",
     all(vapply(c(3L, 5L, 10L), function(pp) {
       if (pp * M > nrow(cfg$Z)) return(NA)
       all(eigen(cfg$Z[seq_len(pp*M), seq_len(pp*M)], symmetric = TRUE)$values > 0)
     }, logical(1)), na.rm = TRUE))
  ok("published 6x6 preserved as the leading block",
     isTRUE(all.equal(unname(cfg$Z[1:6, 1:6]),
                      unname(rbind(c(9,0,.2,.2,.2,.2), c(0,4,.2,.2,.2,.2),
                                   c(.2,.2,8,0,.2,.2), c(.2,.2,0,3.5,.2,.2),
                                   c(.2,.2,.2,.2,8.5,0), c(.2,.2,.2,.2,0,3.5))))))
  # d = 1.5 and 3, chosen on measured separation rather than on the manuscript's
  # 5 and 8, which both sit at ref_ccr 1.0000 and cannot rank methods.
  ok("separations are 1.5 and 3",
     isTRUE(all.equal(sort(unique(as.numeric(cfg$scenarios$d))), c(1.5, 3))))
  # guard against the blind-substitution damage that corrupted this file once:
  # arm numbers must be 1..6 and the arm-specific parameters must survive edits
  ok("arm numbers are 1 to 6",
     isTRUE(all.equal(sort(unique(cfg$scenarios$arm)), 1:6)))
  ok("contaminated arm keeps eps = 0.05 and k = 5",
     all(cfg$scenarios$arm_eps[cfg$scenarios$arm == 4] == 0.05) &&
     all(cfg$scenarios$arm_k[cfg$scenarios$arm == 4] == 5))
  ok("t arms keep nu = 3 and nu = 1.5",
     all(cfg$scenarios$arm_nu[cfg$scenarios$arm == 5] == 3) &&
     all(cfg$scenarios$arm_nu[cfg$scenarios$arm == 6] == 1.5))
  ok("each arm appears at both separations",
     all(table(cfg$scenarios$arm) == 2L))

  cat("\n=== unresolved parameters ===\n")
  tbd <- cfg$meta$param[cfg$meta$value %in% "TBD"]
  if (length(tbd)) {
    for (p in tbd) cat("  TBD:", p, "\n")
    cat("  -> the study cannot run until these are set\n")
  } else cat("  none\n")

  cat("\n=== basis, after the recorded corrections ===\n")
  t <- seq(0, 1, length.out = 20001)
  eig <- cfg$basis[cfg$basis$role == "eigen", ]
  for (v in sort(unique(eig$variable))) for (k in sort(unique(eig$class))) {
    rows <- eig[eig$variable == v & eig$class == k, ]
    rows <- rows[order(rows$component), ]
    f <- lapply(rows$expr, eval_basis, t = t)
    B <- if (any(rows$orthonormalize)) orthonormalize_pair(f[[1]], f[[2]])
         else cbind(f[[1]], f[[2]])
    G <- crossprod(B) / length(t)
    ok(sprintf("variable %d class %d orthonormal", v, k),
       max(abs(G - diag(2))) < 1e-3,
       sprintf("norms %.4f %.4f, inner %+.4f", G[1,1], G[2,2], G[1,2]))
  }

  cat("\n=== all bases are mutually distinct ===\n")
  shape <- function(v, k = 1L) {
    rows <- eig[eig$variable == v & eig$class == k, ]
    rows <- rows[order(rows$component), ]
    f <- lapply(rows$expr, eval_basis, t = t)
    if (any(rows$orthonormalize)) orthonormalize_pair(f[[1]], f[[2]])
    else cbind(f[[1]], f[[2]])
  }
  B <- lapply(seq_len(P_MAX), shape)
  worst <- Inf; worst_lab <- ""
  for (i in seq_len(P_MAX)) for (j in seq_len(i - 1L)) {
    dd <- max(abs(abs(B[[i]]) - abs(B[[j]])))
    if (dd < worst) { worst <- dd; worst_lab <- sprintf("%d vs %d", j, i) }
  }
  ok(sprintf("all %d variable bases pairwise distinct", P_MAX), worst > 1e-3,
     sprintf("closest pair %s at %.4f", worst_lab, worst))

  # Study B's added variables must carry NO class information: same mean and
  # same basis in both classes. If this ever fails, the P axis is measuring
  # added signal rather than added noise and Delta will grow with P.
  cat("\n=== variables 4+ are uninformative (Study B's P axis) ===\n")
  if (P_MAX > 3L) {
    for (v in 4:P_MAX) {
      rows <- cfg$basis[cfg$basis$variable == v, ]
      mu <- rows[rows$role == "mean", ]
      eg <- rows[rows$role == "eigen", ]
      same_mu <- length(unique(mu$expr)) == 1L
      same_eg <- all(vapply(1:2, function(cc)
        identical(eg$expr[eg$class == 1 & eg$component == cc],
                  eg$expr[eg$class == 2 & eg$component == cc]), logical(1)))
      no_d <- !any(grepl("\\bd\\b", mu$expr))
      ok(sprintf("variable %d class-common and free of d", v),
         same_mu && same_eg && no_d,
         if (!no_d) "mean still depends on d" else "")
    }
  } else cat("  (only 3 variables configured; nothing to check)\n")

  cat("\n")
  if (fail > 0L) stop(sprintf("config validation: %d FAILED", fail), call. = FALSE)
  cat("config validation: PASS",
      if (length(tbd)) sprintf("(with %d parameters still TBD)", length(tbd)) else "", "\n")
}
