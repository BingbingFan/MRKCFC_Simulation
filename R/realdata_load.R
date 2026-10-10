# Study D: the loader for `data_yn_use_digit1_CEA_199_242.rds`.
#
# Sourced, not run. It replaces a data contract written before the file
# existed, which assumed one flat long table. The actual
# file is a list of seven, with TWO cohorts and one long frame PER MARKER.
#
# WHAT IS IN THE FILE
#
#   data_yn             the raw wide extract, 59829 x 142. NOT USED: values are
#                       factors carrying ">20000" and full-width "＞1500", and
#                       every column would need its own coercion rule.
#   CEA_199_242_12m     list of 3 long frames, CEA / CA199 / CA242 IN THAT ORDER
#                       (the list is unnamed; `marker_yn_list_12` carries the
#                       same three WITH names, and they agree).
#   dt_base_12          1233 subjects: baseline covariates, Death,
#                       Survival_Time, Recurrence, Recurrence_free_Time.
#   *_36m / dt_base_36  the same for the 36-month window, 1405 subjects.
#
# Each marker frame is `ID, time, yit, freq, yit_lg` -- the raw value and its
# log, with `freq` the subject's observation count.
#
# THREE DECISIONS THIS FILE MAKES EXPLICIT RATHER THAN BY DEFAULT
#
# 1. WHICH COHORT. 12m has 1233 subjects, 36m has 1405. They are nested, not
#    alternatives -- the 36-month window admits subjects followed longer. The
#    manuscript's application must state which it uses. `rd_load(window = 12)`.
#
# 2. RAW OR LOG. The file supplies both `yit` and `yit_lg`. THE ORIGINAL
#    ANALYSIS DID NOT LOG-TRANSFORM. The default here is
#    therefore RAW, to reproduce what was actually done. Passing `log = TRUE`
#    changes Section 3.2's numbers, not merely their scale, so report both and
#    say which is which.
#
#    SINCE 2026-09-11 THE LOG SCALE IS AN ARM OF THE APPLICATION, not a flag
#    someone might set: passing TRUE as the log argument of realdata_fit_task.R runs the
#    primary grid with every marker log-transformed BEFORE the median/IQR
#    normalisation of Equation (13), into its own directory. Every checkpoint
#    carries `log`, the guard and the combiner refuse to mix scales, and
#    `prep$norm` plus rd_value_inverse() below take a fitted value back to the
#    original units for display. The paper reports the original scale as the
#    primary analysis and the log scale as a sensitivity analysis.
#
# 3. ONE SCHEDULE PER SUBJECT. The three markers are NOT measured at identical
#    visits: on the 12-month cohort the union of (ID, time) is 8917 and the
#    intersection 8200, so 8 percent of visits are missing at least one marker.
#    The model takes one schedule per subject shared across variables, so this
#    intersects. That is a real exclusion and the count is reported.

suppressWarnings(suppressMessages(library(MRKCFC)))

RD_FILE <- "data/data_yn_use_digit1_CEA_199_242.rds"
RD_MARKERS <- c("CEA", "CA199", "CA242")   # the order of the unnamed list

rd_stop <- function(...) stop("Study D: ", ..., call. = FALSE)

# ---------------------------------------------------------------- loading ---
# SEARCH A FEW STANDARD PLACES RATHER THAN ONE. A whole array once failed in
# seconds because the file had been copied to a sibling `data/` directory while
# this defaulted to `data/` relative to the working directory. One wrong path,
# fifteen identical failures, no output. The
# candidates below cover the local tree and both plausible cluster layouts, and
# the error names every place looked instead of only the one that was tried.
rd_find <- function(path = RD_FILE) {
  cand <- unique(c(path, file.path("..", path), file.path("..", "..", path),
                   file.path("data", basename(path)),
                   file.path("..", "data", basename(path)),
                   file.path("..", "..", "data", basename(path))))
  hit <- cand[file.exists(cand)]
  if (!length(hit))
    rd_stop("data file not found. Looked in:\n  ",
            paste(normalizePath(cand, mustWork = FALSE), collapse = "\n  "),
            "\nCopy it to one of these, or pass an explicit path to rd_load().")
  normalizePath(hit[1])
}

rd_load <- function(path = RD_FILE, window = c(12, 36)) {
  window <- as.integer(match.arg(as.character(window[1]), c("12", "36")))
  path <- rd_find(path)
  cat("data: ", path, "\n", sep = "")
  D <- readRDS(path)
  mk <- D[[sprintf("CEA_199_242_%dm", window)]]
  bs <- D[[sprintf("dt_base_%d", window)]]
  if (is.null(mk) || is.null(bs)) rd_stop("window ", window, " not in the file")
  if (length(mk) != 3L) rd_stop("expected 3 marker frames, found ", length(mk))
  names(mk) <- RD_MARKERS
  for (m in RD_MARKERS) {
    need <- c("ID", "time", "yit", "yit_lg")
    miss <- setdiff(need, names(mk[[m]]))
    if (length(miss)) rd_stop(m, " is missing: ", paste(miss, collapse = ", "))
  }
  cat(sprintf("cohort: %d-month window | %d subjects in dt_base | markers %s\n",
              window, nrow(bs), paste(RD_MARKERS, collapse = ", ")))
  list(markers = mk, base = bs, window = window)
}

# ------------------------------------------------ THE FLOOR SCREEN ---------
#
# Carried as an open item since 2026-08-26 and now ANSWERABLE. `robust_normalize`
# centres on the median and scales by the IQR; if a marker is mostly at its assay
# detection floor the IQR is zero and the function silently falls back to a 10-90
# spread, putting that marker on a different scale from the others.
#
# Two bands matter and only one trips the guard: above about 75 percent at the
# floor, Q1 == Q3 and the fallback FIRES; between 50 and 75 percent nothing warns
# but the median sits on the floor, so more than half the normalized values come
# out exactly 0.
rd_floor_screen <- function(dat, log = FALSE) {
  col <- if (log) "yit_lg" else "yit"
  cat(sprintf("\n=== DETECTION-FLOOR SCREEN (%s scale) ===\n",
              if (log) "log" else "raw"))
  cat(sprintf("%-8s %7s %9s %9s %9s %9s   %s\n",
              "marker", "n", "min", "at min", "median", "IQR", "verdict"))
  out <- do.call(rbind, lapply(RD_MARKERS, function(m) {
    v <- dat$markers[[m]][[col]]; v <- v[!is.na(v)]
    at <- mean(v <= min(v)); iq <- stats::IQR(v)
    verdict <- if (iq <= 0) "FALLBACK FIRES -- different scale from the others"
      else if (at >= 0.75) "at the transition; check after any transform"
      else if (at >= 0.50) "median ON the floor: >half normalize to exactly 0"
      else "ok"
    cat(sprintf("%-8s %7d %9.4g %8.1f%% %9.4g %9.4g   %s\n",
                m, length(v), min(v), 100 * at, stats::median(v), iq, verdict))
    data.frame(marker = m, n = length(v), min = min(v), frac_at_min = at,
               median = stats::median(v), iqr = iq, verdict = verdict,
               stringsAsFactors = FALSE)
  }))
  if (all(out$verdict == "ok"))
    cat("\nNo marker is piled at a detection floor. The IQR fallback in\n",
        "robust_normalize() will not fire, and all three variables are\n",
        "normalized on the same footing. THE OPEN ITEM IS CLOSED.\n", sep = "")
  invisible(out)
}

# ----------------------------------------------- prepare for MKCFC ---------
#
# Returns Ly[[j]][[i]] and Lt[[i]] on the schedule COMMON to all three markers.
rd_prepare <- function(dat, log = FALSE, obs_min = 2L, verbose = TRUE,
                       time = c("warped", "rescaled")) {
  time <- match.arg(time)
  col <- if (log) "yit_lg" else "yit"
  keyed <- lapply(RD_MARKERS, function(m) {
    x <- dat$markers[[m]][, c("ID", "time", col)]
    names(x)[3] <- "y"
    x <- x[stats::complete.cases(x), ]
    x$key <- paste(x$ID, x$time)
    x[!duplicated(x$key), ]
  })
  names(keyed) <- RD_MARKERS

  common <- Reduce(intersect, lapply(keyed, `[[`, "key"))
  if (verbose) {
    uni <- length(Reduce(union, lapply(keyed, `[[`, "key")))
    cat(sprintf("\nvisits: %d common to all three markers, %d in the union",
                length(common), uni))
    cat(sprintf(" (%.1f%% dropped)\n", 100 * (1 - length(common) / uni)))
  }
  keyed <- lapply(keyed, function(x) x[match(common, x$key), ])

  ref <- keyed[[1]][order(keyed[[1]]$ID, keyed[[1]]$time), ]
  n_obs <- table(ref$ID)
  keep  <- names(n_obs)[n_obs >= obs_min]
  if (verbose && length(keep) < length(n_obs))
    cat(sprintf("dropped %d of %d subjects with fewer than %d visits\n",
                length(n_obs) - length(keep), length(n_obs), obs_min))
  ref <- ref[as.character(ref$ID) %in% keep, ]
  ids <- unique(ref$ID)
  ord <- order(keyed[[1]]$ID, keyed[[1]]$time)

  # Lt IS A P-LIST, NOT A FLAT n-LIST, and this is not obvious from the argument
  # name. `generate.R:151` does `Ly <- Lt <- vector("list", P)`, so the package
  # wants one schedule list PER VARIABLE even when, as here and in the
  # simulation, the schedule is shared. Passing a flat n-list gets
  # "Ly and Lt must contain the same number of outcomes" -- length(Ly) is 3 and
  # length(Lt) was 1405. That failed every task of an array in two seconds.
  # The generator, not the argument name, is the authority on this shape.
  # RESCALE TIME TO [0, 1]. The simulation's basis is defined there and the
  # cohort is observed over 0 to 36 months. Rescaling the domain is a monotone
  # reparametrisation that leaves the function space unchanged.
  #
  # IT DOES NOT PUT THE APPLICATION AT THE SIMULATION'S OPERATING POINT, and
  # this file said for three days that it did. MEASURED 2026-09-09: the
  # simulation draws Poisson(7) visits from a uniform
  # 101-point candidate grid, so its design density is flat on [0, 1]. This
  # cohort's visit times have deciles
  #
  #     50%  0.097   (month 3.5)      90%  0.408  (month 14.7)
  #     80%  0.225   (month 8.1)     100%  1.000  (month 36)
  #
  # -- EIGHTY PERCENT OF ALL VISITS IN THE FIRST EIGHT MONTHS of a 36-month
  # domain, with the upper 60 percent of the interval holding under a tenth of
  # the data. Same interval, entirely different density, so a bandwidth means
  # one thing in the simulation and another here. That is why the application
  # transforms time to uniform design density by default; see `time` below.
  t_range <- range(ref$time)
  # THE TIME AXIS. "warped" (the default since 2026-09-10) maps every visit
  # through the empirical distribution function of all common-schedule visit
  # times, so a bandwidth of 0.2 means twenty percent of the visits everywhere
  # on the axis. "rescaled" is the plain [0,1] rescale of months. Both are
  # monotone; `time_map` carries the knots so any fitted curve is displayed in
  # months through rd_time_inverse().
  #
  # WHY, measured by R/realdata_bwprobe.R on 2026-09-10, which
  # calls KFPCA() exactly as the loop does, on random subsets at the class
  # sizes the loop reaches. Success rate of the robust step, percent:
  #
  #               raw time                       warped time
  #   n_k    GCV  0.10  0.15  0.20  0.30    GCV  0.10  0.15  0.20  0.30
  #   1405   100   100   100   100   100    100   100   100   100   100
  #    350   100    33    33   100   100    100   100   100   100   100
  #    150   100    33    33    33   100    100   100   100   100   100
  #    100   100     0     0    67   100    100   100   100   100   100
  #
  # On raw time GCV survives only by adapting -- it selects 0.084 on the whole
  # cohort and 0.29 at 150 subjects -- and the loop's REAL classes are worse
  # than random subsets of the same size (GCV still failed 52 of 60 array fits):
  # a class of short-follow-up subjects has no late-time data at ANY bandwidth.
  # In warped time the late corner's window reaches back into months 8-36,
  # which every class has, and GCV selects a stable 0.10-0.19 throughout.
  if (time == "warped") {
    Fhat <- stats::ecdf(ref$time)
    f0 <- Fhat(t_range[1]); f1 <- Fhat(t_range[2])
    tw <- function(t) (Fhat(t) - f0) / (f1 - f0)
  } else {
    tw <- function(t) (t - t_range[1]) / diff(t_range)
  }
  knots <- sort(unique(ref$time))
  tmap  <- list(mode = time, months = knots, u = tw(knots), t_range = t_range)
  Lt_one <- lapply(ids, function(i) tw(ref$time[ref$ID == i]))
  Lt <- rep(list(Lt_one), length(RD_MARKERS))
  # NORMALIZE ON THE POOLED OBSERVATIONS, not per subject: per-subject
  # normalization would remove the between-subject variation the clustering
  # exists to find.
  # robust_normalize() returns the centre and scale as attributes; they are
  # kept per marker so a fitted curve can be shown in the original units.
  fitted <- lapply(RD_MARKERS, function(m) {
    x <- keyed[[m]][ord, ]
    x <- x[as.character(x$ID) %in% keep, ]
    z <- MRKCFC:::robust_normalize(x$y)
    list(Ly   = lapply(ids, function(i) as.numeric(z)[x$ID == i]),
         norm = c(centre = as.numeric(attr(z, "centre")),
                  scale  = as.numeric(attr(z, "scale"))))
  })
  names(fitted) <- RD_MARKERS
  Ly   <- lapply(fitted, `[[`, "Ly")
  norm <- lapply(fitted, `[[`, "norm")

  # Compare Lt[[1]], not Lt: both are P-lists now.
  stopifnot(identical(lengths(Lt[[1]]), lengths(Ly[[1]])))
  stopifnot(length(Lt) == length(Ly))
  if (verbose)
    cat(sprintf("prepared: %d subjects, %d markers, %d visits, %d to %d per subject (median %g)\n",
                length(ids), length(Ly), sum(lengths(Lt[[1]])),
                min(lengths(Lt[[1]])), max(lengths(Lt[[1]])),
                stats::median(lengths(Lt[[1]]))))
  # `n` IS CARRIED EXPLICITLY. Callers used `length(prep$Lt)` for the subject
  # count, which silently became 3 the moment Lt became a P-list. A count that
  # is wrong by a factor of 468 and raises no error is worse than a crash.
  if (verbose) {
    cat(sprintf("time axis: %s\n", time))
    cat(sprintf("value scale: %s, then median/IQR normalisation: %s\n",
                if (log) "log (the file's yit_lg)" else "raw (the file's yit)",
                paste(vapply(RD_MARKERS, function(m)
                  sprintf("%s centre %.4g scale %.4g", m, norm[[m]][["centre"]],
                          norm[[m]][["scale"]]), ""), collapse = "; ")))
  }
  list(Ly = Ly, Lt = Lt, ids = ids, n = length(ids), t_range = t_range, log = log,
       window = dat$window, time = time, time_map = tmap, norm = norm)
}

# Subject-level rows aligned to prep$ids, for the survival benchmark.
rd_subjects <- function(dat, ids) {
  b <- dat$base
  i <- match(ids, b$ID)
  if (anyNA(i)) rd_stop(sum(is.na(i)), " fitted subjects are absent from dt_base")
  b[i, , drop = FALSE]
}

# NO rd_optns HERE. IT WAS WRONG AND IT COST TWO SUBMISSIONS.
#
# The two branches take ENTIRELY DIFFERENT option sets, and the hand-written
# version mixed them: the robust branch wants `interval`, `nK`, `bw`, `bwK`,
# `bwmean`, `nRegGrid` and a B-SPLINE BASIS -- which was omitted altogether, so
# the multivariate stacking failed with "number of columns of matrices must
# match" -- while the non-robust branch wants `dataType`, `methodMuCovEst`,
# `nRegGrid`, `methodSelectK`, `useBinnedData` and REJECTS `nK` and `kern` by
# name.
#
# Callers use `build_optns(cfg, robust)` from `methods.R` instead: one
# definition, shared with the simulation, and it cannot drift.
#
# rd_prepare() RESCALES TIME TO [0, 1] so the simulation's B-spline basis
# applies unchanged; the original range is returned as `t_range` so results read
# back in months. That is a reparametrisation of the DOMAIN and nothing more --
# it does not equalise the DESIGN DENSITY, which is where this cohort and the
# simulation part company. See the measured deciles in rd_prepare().

# THE BANDWIDTH IS THE CONFIG'S AGAIN: `bwK = "GCV"`, as the simulation has it.
#
# Two constants sat here on 2026-09-09 and both were wrong. `RD_BW = 0.20` set
# the Nadaraya-Watson bandwidth `bw`, which is not the smoother that fails.
# `RD_BWK = 0.20` set the right smoother -- KFPCA.R:131, `fdapace::Lwls2D`
# with `bw2` from `bwK` -- to a value read off a proxy count of raw pairs per
# window, and that value failed 60 of 60 robust fits where GCV had failed 52.
# The proxy's threshold was not the smoother's requirement, and it was sized
# on balanced classes the loop never produces.
#
# What the real smoother needs was then MEASURED (R/realdata_bwprobe.R; the
# table is in rd_prepare() above): on the rescaled axis GCV succeeds only by
# adapting up to 0.29 on small classes and still fails on the loop's real,
# non-random classes; on the WARPED axis every bandwidth from 0.10 up succeeds
# at every class size and GCV selects a stable 0.10-0.19. The fix was the time
# axis, not the bandwidth. With the axis fixed, GCV is both adequate and the
# setting the simulation uses, so there is nothing to override.
#
# A numeric `bwK` remains available for the sensitivity arm, which fixes it at
# 0.20 to show that K selection does not depend on the bandwidth. It binds on
# the ROBUST branch only; MKCFC's fdapace path selects its own throughout.

# The option list, from the SIMULATION's `build_optns()` in methods.R. One
# definition, so the application and Study A cannot drift apart in anything the
# caller does not deliberately override.
#
# The caller must have sourced `study_c.R` (which sources `methods.R`) first.
# `nK_fit` is overridden HERE rather than edited into the config, so the
# reportable value in params_global.csv is never mutated to run a variant.
# `bwK` is left at the config's "GCV" unless a number is passed.
rd_opt <- function(nK = 4L, robust = TRUE, cfg = NULL, bwK = NULL) {
  if (!exists("build_optns", mode = "function"))
    rd_stop("source R/study_c.R before rd_opt(): build_optns lives in methods.R")
  if (is.null(cfg)) cfg <- load_config()
  cfg$global$nK_fit <- as.integer(nK)
  if (!is.null(bwK)) cfg$global$bwK <- as.numeric(bwK)   # numeric bypasses GCV
  build_optns(cfg, robust)
}

# Months for a point on the fitted time axis, by linear interpolation through
# the knots rd_prepare() recorded. Under "rescaled" this is the affine inverse;
# under "warped" it inverts the empirical distribution function. Every plot of
# a fitted trajectory goes through this so the axis a reader sees is months.
rd_time_inverse <- function(prep, u) {
  m <- prep$time_map
  stats::approx(m$u, m$months, xout = u, rule = 2, ties = "ordered")$y
}

# The original units for a value on the fitted scale of one marker: undo the
# median/IQR normalisation and, on the log arm, the log. The file's `yit_lg`
# is the NATURAL log of `yit` (checked on 2026-09-11 to machine precision), so
# the inverse is exp(). Every plotted trajectory goes through this as well as
# rd_time_inverse(), so the axes a reader sees are months and marker units.
rd_value_inverse <- function(prep, z, marker) {
  nm <- prep$norm[[marker]]
  if (is.null(nm)) rd_stop("no normalisation constants for marker ", marker)
  y <- MRKCFC:::back_transform(z, nm[["centre"]], nm[["scale"]])
  if (isTRUE(prep$log)) exp(y) else y
}
