# Study D: can the ROBUST STEP ITSELF be fitted on this cohort, at what class
# size, at what bandwidth, and does transforming time to uniform design
# density change the answer?
#
#   Rscript R/realdata_bwprobe.R [outdir] [n_draw]
#
# ===========================================================================
# WHY THIS EXISTS, AND WHY IT CALLS THE REAL FUNCTION
#
# Two attempts to fix the robust branch's "No enough points in local window"
# failure were sized from a PROXY -- a count of raw pairs in each window --
# and both made things worse: GCV failed 52 of 60 robust
# fits, a fixed bwK of 0.20 failed 60 of 60, 0.15 failed 40 of 40. The proxy's
# threshold was not the smoother's real requirement, and it was sized on the
# balanced classes of n/K that the K-centres loop never produces; the loop
# drives one class to roughly 100 to 250 subjects.
#
# So this probe calls KFPCA() -- the function that throws -- exactly as
# univar_FPCA_k() calls it inside the loop, on random subsets of the cohort at
# the class sizes the loop actually reaches, across a grid of bwK including
# GCV, and records whether it succeeds and what bandwidth GCV selects.
#
# THEN THE SAME GRID WITH TIME TRANSFORMED TO UNIFORM DESIGN DENSITY. Eighty
# percent of the cohort's visits fall in the first 22 percent of the rescaled
# domain, so any fixed bandwidth is far too wide in the dense region and empty
# in the late-time corners at once. Replacing each visit time by the empirical
# distribution function of all visit times is a monotone map, back-
# transformable to months for every plot, under which a bandwidth of 0.2
# means "twenty percent of the visits" everywhere on the axis.
#
# ONE TABLE DECIDES THE RE-RUN: for each class size and bandwidth, the success
# rate of the robust step, raw against warped. If warping makes classes of 100
# to 250 estimable at a modest bandwidth, that bandwidth is read off the table
# and the arrays are re-run under it. If it does not, warping is not the fix
# and another array on it would have been wasted.
# ===========================================================================

source(file.path("R", "study_c.R"))
source(file.path("R", "realdata_load.R"))

args   <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "results_rd"
n_draw <- if (length(args) >= 2) as.integer(args[2]) else 3L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# The class sizes the loop reaches, not the balanced n/K a diagnostic would
# assume. 1405 is the global fit the combiner makes; 100 to 250 is where the
# minority class sits mid-trajectory in every fit that got that far.
N_K   <- c(1405L, 700L, 350L, 250L, 150L, 100L)
BWK   <- list("GCV", 0.10, 0.15, 0.20, 0.30, 0.40)
BWLAB <- vapply(BWK, function(b) if (is.character(b)) b else format(b), "")

cfg  <- load_config()
nK   <- as.integer(cfg$global$nK_fit)
optn <- rd_opt(nK, TRUE, cfg)          # interval, nK, bw, nRegGrid, basis; bwK replaced below

dat  <- rd_load(window = 36L)
prep <- rd_prepare(dat, log = FALSE, verbose = TRUE)
n    <- prep$n
P    <- length(prep$Ly)

# ---------------------------------------------------- the time transform ---
# Empirical distribution function of ALL visit times, mapped so the earliest
# visit is 0 and the latest is 1. Monotone, so every trajectory keeps its
# shape and the class means are plotted in months by inverting it.
allt  <- unlist(prep$Lt[[1]], use.names = FALSE)
Fhat  <- stats::ecdf(allt)
f0    <- Fhat(min(allt)); f1 <- Fhat(max(allt))
warp  <- function(t) (Fhat(t) - f0) / (f1 - f0)
Lt_w  <- lapply(prep$Lt, function(el) lapply(el, warp))
cat(sprintf("\nwarped time: deciles of visit times move from\n  raw    %s\nto\n  warped %s\n",
            paste(sprintf("%.2f", stats::quantile(allt, seq(0.1, 0.9, 0.2))), collapse = " "),
            paste(sprintf("%.2f", stats::quantile(warp(allt), seq(0.1, 0.9, 0.2))), collapse = " ")))

# ------------------------------------------------------------- one call ---
# Exactly univar_FPCA_k()'s call, so a failure here is the failure in the loop.
one_fit <- function(Lt1, Ly1, bwk) {
  t0 <- Sys.time()
  # `:::` rather than `::`: the package installed on the cluster predates the
  # NAMESPACE regeneration that exports KFPCA, and the loop itself never needs
  # the export because it calls KFPCA from inside the namespace. The first run
  # of this probe failed every cell on exactly that, before touching the data.
  r  <- tryCatch(
    MRKCFC:::KFPCA(Lt = Lt1, Ly = Ly1, interval = optn$interval, nK = optn$nK,
                 bw = optn$bw, bwK = bwk, nRegGrid = optn$nRegGrid,
                 fdParobj = optn$basis),
    error = function(e) e)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(r, "error"))
    return(list(ok = FALSE, bwK = NA_real_, secs = el,
                msg = substr(conditionMessage(r), 1, 60)))
  list(ok = TRUE, bwK = as.numeric(r$bwK), secs = el, msg = "")
}

# ---------------------------------------------------------------- the grid ---
set.seed(20260910L)
rows <- list()
for (scale in c("raw", "warped")) {
  Lt_use <- if (scale == "raw") prep$Lt else Lt_w
  cat(sprintf("\n=== %s time ===\n", scale))
  for (nk in N_K) {
    draws <- if (nk >= n) 1L else n_draw
    for (d in seq_len(draws)) {
      idx <- if (nk >= n) seq_len(n) else sort(sample.int(n, nk))
      for (b in seq_along(BWK)) {
        # every variable, since the loop fits all P and any one can throw
        res <- lapply(seq_len(P), function(p)
          one_fit(Lt_use[[p]][idx], prep$Ly[[p]][idx], BWK[[b]]))
        ok  <- all(vapply(res, `[[`, TRUE, "ok"))
        rows[[length(rows) + 1L]] <- data.frame(
          scale = scale, n_k = nk, draw = d, bwK = BWLAB[b],
          ok = ok,
          bwK_selected = if (BWLAB[b] == "GCV" && ok)
            stats::median(vapply(res, `[[`, 0, "bwK")) else NA_real_,
          secs = sum(vapply(res, `[[`, 0, "secs")),
          msg = paste(unique(vapply(res, `[[`, "", "msg")[!ok]), collapse = " | "),
          stringsAsFactors = FALSE)
      }
      cat(sprintf("  n_k=%4d draw %d: %s\n", nk, d,
                  paste(sprintf("%s:%s", BWLAB,
                                vapply(seq_along(BWK), function(b) {
                                  r <- rows[[length(rows) - length(BWK) + b]]
                                  if (r$ok) "ok" else "FAIL" }, "")),
                        collapse = "  ")))
    }
  }
}
R <- do.call(rbind, rows)
utils::write.csv(R, file.path(outdir, "rd_bwprobe.csv"), row.names = FALSE)

# ---------------------------------------------------------------- report ---
for (scale in c("raw", "warped")) {
  cat(sprintf("\n=== ROBUST STEP SUCCESS RATE, %s time (draws per cell in parentheses) ===\n",
              scale))
  cat(sprintf("%-6s %s\n", "n_k", paste(sprintf("%8s", BWLAB), collapse = "")))
  for (nk in N_K) {
    cells <- vapply(BWLAB, function(bl) {
      s <- R[R$scale == scale & R$n_k == nk & R$bwK == bl, ]
      sprintf("%3.0f%%(%d)", 100 * mean(s$ok), nrow(s))
    }, "")
    cat(sprintf("%-6d %s\n", nk, paste(sprintf("%8s", cells), collapse = "")))
  }
  g <- R[R$scale == scale & R$bwK == "GCV" & R$ok, ]
  if (nrow(g))
    cat(sprintf("GCV selected bwK, median by class size: %s\n",
                paste(vapply(N_K, function(nk) {
                  v <- g$bwK_selected[g$n_k == nk]
                  sprintf("%d:%s", nk, if (length(v)) sprintf("%.3f", stats::median(v)) else "--")
                }, ""), collapse = "  ")))
  m <- unique(R$msg[R$scale == scale & !R$ok]); m <- m[nzchar(m)]
  if (length(m)) cat(sprintf("failure messages: %s\n", paste(m, collapse = " || ")))
}
cat(sprintf("\nwrote %s\n", file.path(outdir, "rd_bwprobe.csv")))
cat("\nHOW TO READ THIS. The loop's minority class sits at 100 to 250 subjects.",
    "\nA bandwidth is usable only if it succeeds there AND at 1405, since the",
    "\ncombiner fits the whole cohort once. Compare the raw and warped tables at",
    "\nthose rows; the smallest bandwidth with 100% in both is the one to run.\n")
