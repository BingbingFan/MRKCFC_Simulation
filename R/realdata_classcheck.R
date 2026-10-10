# Study D: WHAT CLASSES DID THE LOOP FORM, AND WHY COULD THE ROBUST SMOOTHER
# NOT FIT THEM? The `error` fits of the rescaled-axis runs kept their last
# accepted partition. This reads those partitions and, for every class, reports
# the follow-up its members have and whether KFPCA() -- called exactly as the
# loop calls it -- succeeds on that class under the run's own bandwidth and
# axis, and under the warped axis at GCV. It is the evidence for the sentence
# in Section 3.2 that says WHY the time axis had to be transformed.
#
#   Rscript R/realdata_classcheck.R <outstem> <dir> [<dir> ...]
#
# <dir> may be a glob. Every rd_task_*.rds found is read. Files without a
# `status` field (an earlier driver) and fits without a partition
# (`init-only`, `error-early`) are counted and skipped. Files that carry `time`
# are handled on their own axis, so this also documents the warped runs later.
#
# Writes <outstem>_classes.csv (one row per class per axis, appended after
# every fit so a wall hit keeps what was done) and <outstem>_fits.csv, and
# prints the report.
#
# ===========================================================================
# WHAT THE SMOOTHER NEEDS, read from the code on 2026-09-10, and why one
# marker's verdict is every marker's verdict.
#
# KFPCA() has two steps that depend on the visit times ALONE:
#
# 1. GetMatList() builds, for subject i, the Nadaraya-Watson estimate of every
#    other subject j at i's visit times with the compact Epanechnikov kernel
#    and bandwidth `bw`. Where j has no visit strictly within `bw` of one of
#    i's times the estimate is 0/0, and KFPCA drops that ROW. A subject whose
#    every row is dropped is dropped ENTIRELY (KFPCA.R, `sub_id`), and so are
#    all its time pairs. An isolated late visitor never reaches the smoother.
#
# 2. fdapace::Lwls2D() smooths the raw Kendall surface over the 51 x 51 grid
#    with bandwidth `bwK`, and STOPS with "No enough points in local window"
#    at the first grid point whose window holds fewer than three pairs of
#    points at distinct locations. The input locations are the within-subject
#    visit-time pairs (both orderings) of the surviving subjects, so the late
#    corner (1, 1) needs subjects with two or more visits inside the last
#    `bwK` of the axis.
#
# 3. `bwK = "GCV"` (GCVbw2D.R) tries ten candidates from the LARGEST down,
#    r/4 with r = sqrt(2) x the class's time range (r/2 when the design's
#    minimal bandwidth already reaches r/4; see h_effective), catches each
#    failure as an infinite score, and when every candidate fails returns the
#    largest; KFPCA then calls Lwls2D at it and throws. So a GCV failure means
#    the class was infeasible even with a window covering the last 35 percent
#    of the axis, which on the rescaled axis is months 23 to 36.
#
# None of this involves the marker values, and the three markers share one
# schedule (rd_prepare()), so the exact call is made on CEA alone. The window
# count below re-implements 1 and 2 from the visit times, which lets the
# failure be LOCATED on the surface and traced along the loop's history; it is
# validated against the exact call on every class and the agreement is
# reported. The exact call is the verdict; the count is the explanation.
#
# The non-robust branch never failed because fdapace::FPCA on sparse data
# defaults to a Gaussian kernel, whose window is never empty.
# ===========================================================================

source(file.path("R", "study_c.R"))       # load_config, build_optns via methods.R
source(file.path("R", "realdata_load.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L)
  stop("usage: Rscript R/realdata_classcheck.R <outstem> <dir> [<dir> ...]",
       call. = FALSE)
stem <- args[1]
dirs <- unlist(lapply(args[-1], function(a) {
  g <- Sys.glob(a); if (length(g)) g else a }))
dirs <- unique(dirs[dir.exists(dirs)])
if (!length(dirs))
  stop("no directory found among: ", paste(args[-1], collapse = " "), call. = FALSE)
dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
# Classes above this size skip the exact call (the count still runs). Unset
# means every class; set RD_CC_MAXN in the environment if the pace is too slow.
MAX_EXACT_N <- as.integer(Sys.getenv("RD_CC_MAXN", "1000000"))

WINDOW <- 36L
cfg    <- load_config()
nK     <- as.integer(cfg$global$nK_fit)
optn   <- rd_opt(nK, TRUE, cfg)          # interval, nK, bw, nRegGrid, basis
BW_CFG <- as.numeric(optn$bw)            # the config's Nadaraya-Watson bandwidth
GRID   <- seq(optn$interval[1], optn$interval[2], length.out = optn$nRegGrid)

dat  <- rd_load(window = WINDOW)
# Both axes on both scales. The log arm could drop visits whose raw value is
# not positive, so each scale has its own schedule and subject set; a fit is
# always checked against the prep of ITS scale, matched by subject ID.
PREP <- list(
  rescaled = list(raw = rd_prepare(dat, log = FALSE, verbose = TRUE,  time = "rescaled"),
                  log = rd_prepare(dat, log = TRUE,  verbose = FALSE, time = "rescaled")),
  warped   = list(raw = rd_prepare(dat, log = FALSE, verbose = FALSE, time = "warped"),
                  log = rd_prepare(dat, log = TRUE,  verbose = FALSE, time = "warped")))
months_of <- function(p) lapply(p$Lt[[1]], function(u) u * diff(p$t_range) + p$t_range[1])
MONTHS <- lapply(PREP$rescaled, months_of)               # by scale
LAST   <- lapply(MONTHS, function(M) vapply(M, max, 0))
NVIS   <- lapply(MONTHS, lengths)
u_to_months <- function(axis, u, scale = "raw") {
  p <- PREP[[axis]][[scale]]
  if (axis == "rescaled") u * diff(p$t_range) + p$t_range[1] else rd_time_inverse(p, u)
}

cat(sprintf("\ncohort (raw scale): %d subjects; last visit median %.1f months, %.0f%% by month 8, %.0f%% by month 12\n",
            length(PREP$rescaled$raw$ids), stats::median(LAST$raw),
            100 * mean(LAST$raw <= 8), 100 * mean(LAST$raw <= 12)))
if (!identical(PREP$rescaled$log$ids, PREP$rescaled$raw$ids) ||
    !identical(PREP$rescaled$log$Lt, PREP$rescaled$raw$Lt))
  cat(sprintf("  log scale differs: %d subjects, %d visits (raw: %d subjects, %d visits)\n",
              length(PREP$rescaled$log$ids), sum(lengths(PREP$rescaled$log$Lt[[1]])),
              length(PREP$rescaled$raw$ids), sum(lengths(PREP$rescaled$raw$Lt[[1]]))))
for (mo in c(23.3, 28.8, 32.4))
  cat(sprintf("  subjects with >= 2 visits after month %.1f: %d\n", mo,
              sum(vapply(MONTHS$raw, function(m) sum(m >= mo) >= 2L, TRUE))))

# ------------------------------------------------ the count, steps 1 and 2 --
# Step 1: subject i survives iff some other member j has a visit strictly
# within bw of EVERY visit of i.
nw_survivors <- function(Lt1, bw) {
  n <- length(Lt1)
  if (n < 2L) return(integer(0))
  tt <- unlist(Lt1, use.names = FALSE)
  si <- rep.int(seq_len(n), lengths(Lt1))
  keep <- logical(n)
  for (j in seq_len(n)) {
    tj <- sort(Lt1[[j]]); nj <- length(tj)
    k  <- findInterval(tt, tj)                       # tj[k] <= tt < tj[k+1]
    lo <- ifelse(k >= 1L, tj[pmax(k, 1L)], Inf)
    hi <- ifelse(k <  nj, tj[pmin(k + 1L, nj)], Inf)
    d  <- pmin(abs(tt - lo), abs(hi - tt))
    miss <- rowsum(as.integer(!(d < bw)), si)[, 1]   # i's visits NOT covered by j
    ok <- miss == 0L; ok[j] <- FALSE
    keep <- keep | ok
  }
  which(keep)
}

# The smoother's input locations: within-subject time pairs, both orderings.
pair_locations <- function(Lt1, sub) {
  L <- lapply(Lt1[sub], function(t) {
    if (length(t) < 2L) return(NULL)
    cb <- t(utils::combn(t, 2L)); rbind(cb, cb[, 2:1, drop = FALSE]) })
  P <- do.call(rbind, L)
  if (is.null(P)) matrix(numeric(0), 0L, 2L) else P
}

# Step 2 at every grid point (x, y), x <= y: among the locations within h of it
# on both axes, at least three pairs of points at DIFFERENT locations, i.e.
# C(m,2) - sum over locations C(m_loc,2) >= 3. Three distinct first coordinates
# settle it at once, which is the common case. Returns the failing points.
infeasible_points <- function(P, h, first_only = FALSE, grid = GRID) {
  ng <- length(grid); tol <- h + 1e-6
  bad <- matrix(numeric(0), 0L, 2L)
  if (!nrow(P)) {
    all <- do.call(rbind, lapply(seq_len(ng), function(i) cbind(grid[i], grid[i:ng])))
    return(if (first_only) all[1, , drop = FALSE] else all)
  }
  o <- order(P[, 1]); s <- P[o, 1]; t <- P[o, 2]
  for (i in seq_len(ng)) {
    x <- grid[i]
    a <- findInterval(x - tol, s, left.open = TRUE) + 1L   # first s >= x - tol
    b <- findInterval(x + tol, s)                          # last  s <= x + tol
    if (b < a) {
      pts <- cbind(x, grid[i:ng])
      if (first_only) return(pts[1, , drop = FALSE])
      bad <- rbind(bad, pts); next
    }
    sw <- s[a:b]; tw <- t[a:b]                             # sw is sorted
    for (j in i:ng) {
      y <- grid[j]
      w <- which(abs(tw - y) <= tol)
      m <- length(w); ok <- FALSE
      if (m >= 3L) {
        sv <- sw[w]
        if (sum(diff(sv) > 0) >= 2L) ok <- TRUE
        else {
          tab <- table(paste(sv, tw[w]))
          ok  <- (choose(m, 2) - sum(choose(tab, 2))) >= 3
        }
      }
      if (!ok) {
        if (first_only) return(matrix(c(x, y), 1L))
        bad <- rbind(bad, c(x, y))
      }
    }
  }
  bad
}

# The smallest bandwidth at which no window fails, to 0.005. Feasibility is
# monotone in h because a wider window contains the narrower one.
h_min_feasible <- function(P, lo = 0, hi = 0.5, tol = 0.005) {
  while (nrow(infeasible_points(P, hi, TRUE))) {
    lo <- hi; hi <- hi * 1.5
    if (hi > 3) return(Inf)
  }
  while (hi - lo > tol) {
    mid <- (lo + hi) / 2
    if (nrow(infeasible_points(P, mid, TRUE))) lo <- mid else hi <- mid
  }
  hi
}

# The effective smoother bandwidth: the number given, or under GCV the largest
# candidate, which is the one a GCV failure means (point 3 above). The grid's
# top is r/4, except that when the design's own minimal bandwidth (GetMinb,
# the selector's starting point) already reaches r/4 the grid is built from
# there and tops out at r/2. That case is real: on a 33-subject class GCV
# selected 0.64 against an r/4 of 0.32 (class check of 2026-09-10, the only
# two of 328 checks where the count and the exact call disagreed).
h_effective <- function(bwK, Lt1) {
  if (is.numeric(bwK)) return(as.numeric(bwK))
  obs <- sort(unique(unlist(Lt1, use.names = FALSE)))
  r   <- diff(range(obs)) * sqrt(2)
  h0  <- tryCatch(as.numeric(MRKCFC:::GetMinb(Lt1, obsGrid = obs, dataType = "Sparse")),
                  error = function(e) NA_real_)
  if (!is.finite(h0) || h0 < r / 4) r / 4 else r / 2
}

# -------------------------------------------------------- the exact call --
# univar_FPCA_k()'s call, verbatim, on CEA. `:::` because the cluster's
# installed package predates the NAMESPACE that exports KFPCA.
exact_kfpca <- function(mem_ids, axis, scale, bw, bwK) {
  p   <- PREP[[axis]][[scale]]
  idx <- match(mem_ids, p$ids)
  Lt1 <- p$Lt[[1]][idx]
  Ly1 <- p$Ly[[1]][idx]
  t0  <- Sys.time()
  r <- tryCatch(
    MRKCFC:::KFPCA(Lt = Lt1, Ly = Ly1, interval = optn$interval, nK = optn$nK,
                  bw = bw, bwK = bwK, nRegGrid = optn$nRegGrid,
                  fdParobj = optn$basis),
    error = function(e) e)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(r, "error"))
    list(ok = FALSE, bwK = NA_real_, secs = secs,
         msg = substr(conditionMessage(r), 1, 70))
  else list(ok = TRUE, bwK = as.numeric(r$bwK), secs = secs, msg = "")
}

# ------------------------------------------------------------ the files --
# What each file was fitted under. Three generations of checkpoint exist:
#   `bw` + `bw_binds`, no `bwK`, no `time`   the 2026-09-09 GCV runs, where the
#                                             misattributed override widened
#                                             the NW bandwidth to 0.20
#   `bwK` + `bwK_binds`, no `time`            the fixed-bwK runs, NW at config
#   `bwK` + `time`                            the current driver
# Non-robust fits are checked under the config's robust settings: the question
# there is whether the robust smoother could have fitted the classes MKCFC
# formed, not what MKCFC did.
read_fit <- function(f) {
  p <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(p) || !is.character(p$status)) return(list(skip = "no status field"))
  lab <- if (!is.null(p$final_labels)) p$final_labels else p$labels
  if (is.null(lab)) return(list(skip = p$status))
  robust <- !identical(p$variant, "MKCFC")
  list(skip = "", file = basename(f), dir = basename(dirname(f)),
       scale = if (isTRUE(as.logical(p$log))) "log" else "raw",
       variant = p$variant, K = as.integer(p$K), restart = as.integer(p$restart),
       status = p$status, iter = as.integer(p$iter), chunks = as.integer(p$chunks),
       labels = as.integer(lab),
       which  = if (!is.null(p$final_labels)) "final" else "best",
       ids = p$ids, history = p$history, robust = robust,
       time = if (is.null(p$time)) "rescaled" else as.character(p$time),
       bw   = if (!robust) BW_CFG
              else if (!is.null(p$bwK)) BW_CFG
              else if (!is.null(p$bw)) as.numeric(p$bw) else BW_CFG,
       bwK  = if (!robust || is.null(p$bwK)) "GCV"
              else if (is.numeric(p$bwK)) as.numeric(p$bwK) else "GCV",
       error = if (is.null(p$error)) NA_character_ else as.character(p$error))
}

fits <- list(); skipped <- character(0)
for (d in dirs) {
  fs <- list.files(d, pattern = "^rd_task_[0-9]+[.]rds$", full.names = TRUE)
  got <- lapply(fs, read_fit)
  sk  <- vapply(got, `[[`, "", "skip")
  cat(sprintf("\n=== %s: %d files ===\n", d, length(fs)))
  if (length(fs)) {
    st <- vapply(got, function(g) if (nzchar(g$skip)) g$skip else g$status, "")
    print(table(st, dnn = NULL))
    kept <- got[!nzchar(sk)]
    if (length(kept)) {
      des <- unique(vapply(kept, function(g)
        sprintf("%s axis, %s scale, NW bw %.2f, bwK %s%s", g$time, g$scale, g$bw,
                if (is.numeric(g$bwK)) format(g$bwK) else g$bwK,
                if (g$robust) "" else " (MKCFC: config settings)"), ""))
      cat("designs among fits with a partition:\n  ", paste(des, collapse = "\n   "), "\n")
    }
    fits <- c(fits, kept)
  }
  skipped <- c(skipped, sk[nzchar(sk)])
}
if (!length(fits)) stop("no fit with a partition in ", paste(dirs, collapse = ", "), call. = FALSE)
cat(sprintf("\n%d fits with a partition; %d files skipped (%s)\n", length(fits), length(skipped),
            paste(sprintf("%s %d", names(table(skipped)), table(skipped)), collapse = ", ")))

# Robust error fits first: they are the evidence. Then the robust fits that
# finished, then MKCFC, so a wall hit loses the contrast, not the point.
pri <- vapply(fits, function(g) (if (g$robust) 0L else 10L) +
                (if (g$status %in% c("error", "degenerate")) 0L else 1L), 0L)
fits <- fits[order(pri, vapply(fits, `[[`, "", "variant"),
                   vapply(fits, `[[`, 0L, "K"), vapply(fits, `[[`, 0L, "restart"))]

# ----------------------------------------------------------- per class --
f_classes <- paste0(stem, "_classes.csv")
f_fits    <- paste0(stem, "_fits.csv")
unlink(c(f_classes, f_fits))
append_csv <- function(df, f) {
  new <- !file.exists(f)
  utils::write.table(df, f, sep = ",", row.names = FALSE, col.names = new,
                     append = !new, qmethod = "double")
}

class_check <- function(mem_ids, axis, scale, bw, bwK) {
  p   <- PREP[[axis]][[scale]]
  Lt1 <- p$Lt[[1]][match(mem_ids, p$ids)]
  h   <- h_effective(bwK, Lt1)
  sub <- nw_survivors(Lt1, bw)
  P   <- pair_locations(Lt1, sub)
  supply <- sum(vapply(Lt1[sub], function(t) sum(t >= 1 - h) >= 2L, TRUE))
  bad  <- infeasible_points(P, h)
  list(h = h, n_dropped = length(Lt1) - length(sub), supply = supply,
       corner_month = u_to_months(axis, 1 - h, scale),
       n_bad = nrow(bad),
       # the earliest time beyond which some window is empty: min over failing
       # points of the LARGER coordinate (x <= y by construction)
       first_bad_month = if (nrow(bad)) u_to_months(axis, min(bad[, 2]), scale) else NA_real_,
       h_min = h_min_feasible(P), count_ok = nrow(bad) == 0L)
}

fit_rows <- list()
for (g in fits) {
  sc  <- g$scale
  idx <- match(g$ids, PREP$rescaled[[sc]]$ids)
  if (anyNA(idx) || length(idx) != length(g$labels)) {
    cat(sprintf("%s %s: ids do not match the %s-scale cohort; skipped\n", g$dir, g$file, sc)); next
  }
  classes  <- sort(unique(g$labels))
  run_axis <- g$time
  other    <- setdiff(c("rescaled", "warped"), run_axis)
  sizes    <- as.integer(table(factor(g$labels, levels = classes)))
  cat(sprintf("\n[%s] %s K=%d r%02d %s iter %s | %s labels | %s scale | sizes %s\n",
              g$dir, g$variant, g$K, g$restart, g$status, format(g$iter), g$which, sc,
              paste(sizes, collapse = "/")))
  rows <- list()
  for (k in classes) {
    mem_ids <- g$ids[g$labels == k]
    mem  <- idx[g$labels == k]
    desc <- data.frame(dir = g$dir, file = g$file, variant = g$variant, K = g$K,
                       restart = g$restart, status = g$status, labels = g$which,
                       scale = sc, iter = g$iter, class = k, n_k = length(mem),
                       med_nvis = stats::median(NVIS[[sc]][mem]),
                       med_last_month = stats::median(LAST[[sc]][mem]),
                       p_last_le8  = mean(LAST[[sc]][mem] <= 8),
                       p_last_le12 = mean(LAST[[sc]][mem] <= 12),
                       max_last_month = max(LAST[[sc]][mem]), stringsAsFactors = FALSE)
    for (axis in c(run_axis, other)) {
      role <- if (axis == run_axis) "run" else "other"
      bw   <- if (role == "run") g$bw  else BW_CFG
      bwK  <- if (role == "run") g$bwK else "GCV"
      cc   <- class_check(mem_ids, axis, sc, bw, bwK)
      ex   <- if (length(mem) <= MAX_EXACT_N) exact_kfpca(mem_ids, axis, sc, bw, bwK)
              else list(ok = NA, bwK = NA_real_, secs = NA_real_, msg = "skipped by RD_CC_MAXN")
      rows[[length(rows) + 1L]] <- cbind(desc, data.frame(
        axis = axis, role = role, bw_nw = bw,
        bwK_arg = if (is.numeric(bwK)) format(bwK) else bwK, h_eff = cc$h,
        corner_month = cc$corner_month, nw_dropped = cc$n_dropped,
        corner_supply = cc$supply, n_bad_windows = cc$n_bad,
        first_bad_month = cc$first_bad_month, h_min = cc$h_min,
        count_ok = cc$count_ok, exact_ok = ex$ok, exact_bwK = ex$bwK,
        exact_secs = round(ex$secs, 1), exact_msg = ex$msg,
        agree = if (is.na(ex$ok)) NA else identical(cc$count_ok, ex$ok),
        stringsAsFactors = FALSE))
      cat(sprintf("  class %d n=%4d | last visit: median %4.1f mo, %3.0f%% by mo 8, max %4.1f | %-8s %-5s h=%.3f corner>=mo %4.1f: supply %3d, dropped %3d, bad windows %4d%s, h_min %s | exact %s%s (%.0fs)\n",
                  k, length(mem), desc$med_last_month, 100 * desc$p_last_le8,
                  desc$max_last_month, axis, role, cc$h, cc$corner_month,
                  cc$supply, cc$n_dropped, cc$n_bad,
                  if (cc$n_bad) sprintf(" beyond mo %.1f", cc$first_bad_month) else "",
                  if (is.finite(cc$h_min)) sprintf("%.3f", cc$h_min) else "Inf",
                  if (is.na(ex$ok)) "skipped" else if (ex$ok) "ok" else "FAIL",
                  if (isTRUE(ex$ok) && !is.numeric(bwK)) sprintf(" bwK %.3f", ex$bwK) else "",
                  if (is.na(ex$secs)) 0 else ex$secs))
    }
  }
  R <- do.call(rbind, rows)
  append_csv(R, f_classes)

  # The loop's own trajectory, for the fits that stopped on an error: the
  # corner supply of every class at each chunk boundary, on the run axis at the
  # run's effective bandwidth. Cheap (no exact call, no NW drop), and it shows
  # the class thinning toward the failure rather than starting there.
  if (g$robust && g$status %in% c("error", "degenerate") && length(g$history) > 1L) {
    p_run  <- PREP[[run_axis]][[sc]]
    Lt_run <- p_run$Lt[[1]]
    cat("  trajectory (chunk boundaries): sizes | members with >= 2 visits in the corner window | windows feasible\n")
    for (b in seq_along(g$history)) {
      hl <- as.integer(g$history[[b]])
      if (length(hl) != length(idx)) next
      cl <- sort(unique(hl))
      line <- vapply(cl, function(k) {
        mem <- match(g$ids[hl == k], p_run$ids); Lt1 <- Lt_run[mem]; h <- h_effective(g$bwK, Lt1)
        sup <- sum(vapply(Lt1, function(t) sum(t >= 1 - h) >= 2L, TRUE))
        feas <- nrow(infeasible_points(pair_locations(Lt1, seq_along(Lt1)), h, TRUE)) == 0L
        sprintf("%d:%d:%s", length(mem), sup, if (feas) "y" else "N") }, "")
      cat(sprintf("    b%-2d %s\n", b, paste(line, collapse = "  ")))
    }
  }

  run  <- R[R$role == "run", ]; oth <- R[R$role == "other", ]
  worst <- which.max(ifelse(is.finite(run$h_min), run$h_min, 99))
  fit_rows[[length(fit_rows) + 1L]] <- data.frame(
    dir = g$dir, file = g$file, variant = g$variant, K = g$K, restart = g$restart,
    status = g$status, scale = sc, iter = g$iter, sizes = paste(sizes, collapse = "/"),
    run_axis = run_axis, bwK_arg = run$bwK_arg[1], h_eff = run$h_eff[1],
    worst_class = run$class[worst], worst_n = run$n_k[worst],
    worst_med_last = run$med_last_month[worst], worst_p_le8 = run$p_last_le8[worst],
    worst_supply_run = run$corner_supply[worst], worst_h_min_run = run$h_min[worst],
    run_exact_all_ok = all(run$exact_ok, na.rm = TRUE),
    run_exact_n_fail = sum(!run$exact_ok, na.rm = TRUE),
    other_axis = other, worst_h_min_other = oth$h_min[worst],
    other_exact_all_ok = all(oth$exact_ok, na.rm = TRUE),
    other_bwK_med = if (any(!is.na(oth$exact_bwK))) stats::median(oth$exact_bwK, na.rm = TRUE) else NA_real_,
    agree = sum(R$agree, na.rm = TRUE), n_checked = sum(!is.na(R$agree)),
    error = g$error, stringsAsFactors = FALSE)
  append_csv(fit_rows[[length(fit_rows)]], f_fits)
}
Fz <- do.call(rbind, fit_rows)
Cz <- utils::read.csv(f_classes, stringsAsFactors = FALSE)

# ----------------------------------------------------------------- report --
show <- function(sub, title) {
  cat(sprintf("\n=== %s (%d fits) ===\n", title, nrow(sub)))
  if (!nrow(sub)) return(invisible())
  cat(sprintf("%-22s %-10s %-3s %-14s %-14s | %5s %7s %6s %6s %7s | %-6s | %-6s %7s %6s\n",
              "variant K restart", "status", "it", "run axis/bwK", "sizes",
              "class", "n", "medLV", "%<=8", "supply", "exact", "other", "h_minR", "h_minO"))
  for (i in seq_len(nrow(sub))) {
    z <- sub[i, ]
    cat(sprintf("%-22s %-10s %-3s %-14s %-14s | %5d %7d %6.1f %5.0f%% %7d | %-6s | %-6s %7s %6s\n",
                sprintf("%s K=%d r%02d", z$variant, z$K, z$restart), z$status,
                format(z$iter), sprintf("%s/%s", z$run_axis, z$bwK_arg), z$sizes,
                z$worst_class, z$worst_n, z$worst_med_last, 100 * z$worst_p_le8,
                z$worst_supply_run,
                if (z$run_exact_all_ok) "ok" else sprintf("FAIL%d", z$run_exact_n_fail),
                if (z$other_exact_all_ok) "ok" else "FAIL",
                if (is.finite(z$worst_h_min_run)) sprintf("%.3f", z$worst_h_min_run) else "Inf",
                if (is.finite(z$worst_h_min_other)) sprintf("%.3f", z$worst_h_min_other) else "Inf"))
  }
  cat("  class = the class with the largest minimal feasible bandwidth on the run axis; medLV = its",
      "members' median last visit in months;\n  supply = its members with >= 2 visits inside the",
      "corner window after the NW drop; exact = KFPCA on the run axis at the run's bwK,\n",
      " other = KFPCA on the other axis at GCV; h_minR/h_minO = that class's minimal feasible",
      "bandwidth on the run axis / the other axis.\n")
}
rng <- function(v, fmt = "%.1f") {
  v <- v[is.finite(v)]
  if (!length(v)) return("n/a")
  if (min(v) == max(v)) sprintf(fmt, min(v))
  else paste(sprintf(fmt, min(v)), "to", sprintf(fmt, max(v)))
}
rob <- Fz[Fz$variant != "MKCFC", ]
show(rob[rob$status %in% c("error", "degenerate"), ], "ROBUST FITS THAT STOPPED ON THE SMOOTHER")
show(rob[!rob$status %in% c("error", "degenerate"), ], "ROBUST FITS THAT FINISHED")
show(Fz[Fz$variant == "MKCFC", ], "MKCFC PARTITIONS UNDER THE ROBUST SMOOTHER (config settings)")

cat(sprintf("\nwindow count against the exact call: agrees in %d of %d class-axis checks%s\n",
            sum(Fz$agree), sum(Fz$n_checked),
            if (sum(Fz$agree) < sum(Fz$n_checked)) " -- READ THE DISAGREEMENTS in the classes csv before quoting the count" else ""))

# ------------------------------------------------------- for the paper --
e <- rob[rob$status %in% c("error", "degenerate"), ]
if (nrow(e)) {
  fails <- e[!e$run_exact_all_ok, ]
  cat(sprintf("\nFOR SECTION 3.2\n%d robust fits stopped on the smoother; the exact call fails on the saved partition of %d of them",
              nrow(e), nrow(fails)))
  if (nrow(fails) < nrow(e))
    cat(sprintf(" (the other %d saved a partition from before the failing reassignment)", nrow(e) - nrow(fails)))
  cat(".\n")
  if (nrow(fails)) {
    cat(sprintf("In those %d, the class the smoother could not fit had %s members (the smallest class in %d of %d),\n",
                nrow(fails), rng(fails$worst_n, "%d"),
                sum(vapply(seq_len(nrow(fails)), function(i)
                  fails$worst_n[i] == min(as.integer(strsplit(fails$sizes[i], "/")[[1]])), TRUE)),
                nrow(fails)))
    cat(sprintf("a median last visit of %s months, and %s members with two visits inside the smoother's corner window;\n",
                rng(fails$worst_med_last), rng(fails$worst_supply_run, "%d")))
    cat(sprintf("the smallest bandwidth that would have fitted it was %s of the %s axis against the %s used.\n",
                rng(fails$worst_h_min_run, "%.2f"), paste(unique(fails$run_axis), collapse = "/"),
                paste(unique(fails$bwK_arg), collapse = "/")))
    cat(sprintf("On the %s axis at GCV the same partitions are estimable in %d of %d (GCV selecting bwK %s), with minimal feasible bandwidth %s.\n",
                paste(unique(fails$other_axis), collapse = "/"),
                sum(fails$other_exact_all_ok), nrow(fails),
                rng(fails$other_bwK_med, "%.2f"), rng(fails$worst_h_min_other, "%.2f")))
  }
  fin <- rob[!rob$status %in% c("error", "degenerate"), ]
  if (nrow(fin))
    cat(sprintf("The %d robust fits that finished had worst-class supply %s and minimal feasible bandwidth %s.\n",
                nrow(fin), rng(fin$worst_supply_run, "%d"), rng(fin$worst_h_min_run, "%.2f")))
  mk <- Fz[Fz$variant == "MKCFC", ]
  if (nrow(mk))
    cat(sprintf("MKCFC's own partitions: the robust smoother would have failed on %d of %d of them on the %s axis at GCV, and on %d of %d on the %s axis (fdapace's Gaussian kernel has no empty window, which is why MKCFC itself never failed).\n",
                sum(!mk$run_exact_all_ok), nrow(mk), paste(unique(mk$run_axis), collapse = "/"),
                sum(!mk$other_exact_all_ok), nrow(mk), paste(unique(mk$other_axis), collapse = "/")))
}
cat(sprintf("\nwrote %s\n      %s\n", f_classes, f_fits))
