# Study C3: selecting K. The only place in the project where K is estimated.
#
#   Rscript R/study_c3.R <task_id> [n_tasks] [outdir] [n_reps_override]
#
# GRID IS K = 2 TO 6, NOT 1 TO 6, and the plan's wording is amended.
# Neither criterion is defined at K = 1 -- both correctly return NA -- and in
# simulation K_true = 2 always, so the K = 1 rung is decorative: it can never be
# the right answer. Where the question is live is Study D, on real data, where
# "are there clusters at all" is genuine. STATE THE LIMITATION IN THE PAPER:
# these criteria cannot report the absence of clustering. A gap statistic would
# fix it and was declined here as work that buys nothing in simulation.
#
# TWO CRITERIA, BOTH DISTANCES. The silhouette and CH_mv (Calinski-Harabasz),
# each under Euclidean and robust Mahalanobis distance, giving four columns.
# select_K() computes them at no extra cost, and C being a sensitivity study,
# showing that the chosen K does not depend on the metric is exactly the
# evidence it exists to produce. If they disagree, that disagreement is itself
# the finding. Z's diagonal runs 9, 4, 8, 3.5, 8.5, 3.5 -- a 2.6x spread -- so
# Euclidean distance on raw scores is dominated by the high-variance
# components and the two columns are not expected to be redundant.
#
# THE SCORES ARE COMPUTED ONCE PER REPLICATE AND SHARED ACROSS K. A criterion
# compares partitions of a FIXED geometry; recomputing the score space per K
# would move the space and the partition together and the comparison would no
# longer be about K.
#
# SEVEN METHODS, all but lcmm. The MKCFC family costs 200-400 s a fit and the
# score-then-cluster comparators 1-3 s, so including them is effectively free
# and shows whether K recovery is a property of the method or of the data.
# lcmm would roughly double C3 to select a K that is close to meaningless,
# given it returns chance everywhere on this design and does not converge.

source(file.path("R", "study_c.R"))

args     <- commandArgs(trailingOnly = TRUE)
task_id  <- if (length(args) >= 1) as.integer(args[1]) else 1L
n_tasks  <- if (length(args) >= 2) as.integer(args[2]) else 1L
outdir   <- if (length(args) >= 3) args[3] else "results_c3"
rep_over <- if (length(args) >= 4) as.integer(args[4]) else NA_integer_
# Subject count override, 5th argument. For SMOKE TESTS only: it makes a cell
# cheap enough to exercise the driver without paying for a full n = 500 fit.
# Like the replicate override it is an argument rather than a config edit, so
# the reportable value in params_global.csv is never mutated to run a check.
n_over   <- if (length(args) >= 5) as.integer(args[5]) else NA_integer_

cfg     <- load_config()
if (!is.na(n_over)) {
  cfg$global$n <- n_over
  cat(sprintf("SMOKE TEST: n overridden to %d (reportable value is %s)
", n_over, load_config()$global$n))
}
configs <- load_c_configs()
n_reps  <- c_reps(cfg, rep_over)
K_GRID  <- 2:6
use_methods <- c("M1", "M2", "M3", "M4", "M5", "M6", "M8")

sel_K   <- MRKCFC:::select_K

# THE UNIT OF WORK IS ONE (configuration, replicate, K). Changed 2026-09-06.
#
# It was one (configuration, replicate) -- every method at all five K, about 15
# fits. Ten tasks then hit the 24 h wall, and `c_write()` runs only after the
# whole cell, so each timeout lost the entire task rather than part of it.
# Slicing on K makes a unit about 3 fits.
#
# THIS IS SAFE, and neither property is accidental:
#
#   `select_K()` computes CH and the silhouette PER PARTITION --
#   `vapply(partitions, ch_index, ...)` -- so one (cell, K) can be scored on its
#   own. Nothing is computed jointly across K.
#
#   `c_seed_fit()` is called before EVERY K fit and seeds on (config, rep), not
#   on position in the loop, so each K already starts from the same RNG state.
#   A K-sliced run therefore reproduces the cell-sliced one exactly. This is the
#   property C1 did NOT have, which is why its old files had to be discarded.
#
# The score space is a property of the CELL, so it is cached below rather than
# rebuilt for every K of the same cell.
grid <- expand.grid(K = K_GRID, rep = seq_len(n_reps),
                    config_id = configs$config_id, stringsAsFactors = FALSE)
grid <- grid[order(grid$config_id, grid$rep, grid$K), ]
rownames(grid) <- NULL
mine <- c_slice(grid, task_id, n_tasks)
if (!nrow(mine)) { cat("nothing to do\n"); quit(save = "no", status = 0) }
c_header("C3 selecting K", nrow(mine), nrow(grid), task_id, n_tasks)

started <- Sys.time()
out <- list()
# Cache the generated replicate and its score space by cell, so a task holding
# several K of the same cell pays for them once. Under modulo slicing a task's
# units are usually different cells, in which case this simply never hits.
cache_key <- ""; cache <- NULL
for (a in seq_len(nrow(mine))) {
  cell <- mine[a, ]
  row  <- configs[configs$config_id == cell$config_id, ]
  t0 <- Sys.time()
  cellrows <- tryCatch({
    c2  <- apply_config_c(cfg, row)
    key <- paste(cell$config_id, cell$rep, sep = "|")
    if (!identical(key, cache_key)) {
      dat <- generate_replicate(c2, as_scenario_c(row), rep = cell$rep)
      # The shared geometry. Robust scores, so the space itself survives the
      # heavy-tailed cell; the criteria's own `robust` flag is separate and
      # governs their location and scatter estimates.
      cache <- list(dat = dat,
                    scores = global_mfpca_scores(dat$Ly, dat$Lt, c2, robust = TRUE))
      cache_key <- key
    }
    dat <- cache$dat; scores <- cache$scores

    do.call(rbind, lapply(use_methods, function(mid) {
      parts <- list(); ok_K <- integer(0)
      for (KK in cell$K) {   # ONE K per unit since 2026-09-06; see the grid note
        cK <- c2; cK$global$K_true <- KK
        # Common random numbers across methods AND across K within this
        # replicate, so a difference between two K is a property of K rather
        # than of whichever fit happened to run before it.
        c_seed_fit(c2, cell$config_id, cell$rep)
        fit <- tryCatch(METHODS[[mid]](dat$Ly, dat$Lt, KK, cK),
                        error = function(e) NULL)
        if (is.null(fit) || anyNA(fit$labels) ||
            length(unique(fit$labels)) < 2L) next
        parts[[length(parts) + 1L]] <- fit$labels
        ok_K <- c(ok_K, KK)
      }
      if (!length(parts)) return(NULL)
      eu <- sel_K(scores, parts, distance = "euclidean",   robust = TRUE)
      mh <- sel_K(scores, parts, distance = "mahalanobis", robust = TRUE)
      data.frame(method = mid, K = ok_K,
                 CH_euclid = eu$CH,  sil_euclid = eu$silhouette,
                 CH_mahal  = mh$CH,  sil_mahal  = mh$silhouette,
                 # `error` IS REQUIRED, even though nothing fails on this path.
                 # The tryCatch below returns a frame WITH it, so without this
                 # the two frames are 6 and 7 columns and do.call(rbind, out)
                 # would abort the WHOLE task on the first failed cell --
                 # discarding every cell already completed, and reporting an
                 # rbind error instead of the real cause. Fixed 2026-09-04.
                 error = NA_character_,
                 stringsAsFactors = FALSE)
    }))
  }, error = function(e) data.frame(
    method = NA_character_, K = NA_integer_,
    CH_euclid = NA_real_, sil_euclid = NA_real_,
    CH_mahal = NA_real_, sil_mahal = NA_real_,
    error = paste("cell failed:", conditionMessage(e)), stringsAsFactors = FALSE))

  if (!is.null(cellrows) && nrow(cellrows)) {
    cellrows$config_id <- cell$config_id
    cellrows$role <- row$role
    cellrows$rep <- cell$rep
    cellrows$K_true <- as.integer(cfg$global$K_true)
    out[[length(out) + 1L]] <- cellrows
  }
  cat(sprintf("  [%3d/%3d] %-7s rep %3d  %6.1fs\n", a, nrow(mine),
              cell$config_id, cell$rep,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  utils::flush.console()
}

results <- do.call(rbind, out)
results$method_id <- results$method
results$method <- cfg$methods$label[match(results$method_id, cfg$methods$method_id)]

# The headline: does each criterion's argmax recover K_true? Reported here as a
# convenience; the raw curves are in the file so a reader can see how sharp the
# maximum was, which matters as much as where it fell.
cat("\n=== K selected by each criterion (argmax over K = 2..6) ===\n")
for (crit in c("CH_euclid", "sil_euclid", "CH_mahal", "sil_mahal")) {
  sp <- split(results, list(results$method, results$config_id, results$rep),
              drop = TRUE)
  kk <- vapply(sp, function(d) if (all(is.na(d[[crit]]))) NA_integer_
               else d$K[which.max(d[[crit]])], integer(1))
  cat(sprintf("  %-11s  recovered K_true in %d of %d fits (modal K = %s)\n",
              crit, sum(kk == as.integer(cfg$global$K_true), na.rm = TRUE),
              sum(!is.na(kk)),
              paste(names(sort(table(kk), decreasing = TRUE))[1])))
}

c_write(results, outdir, "studyC3", task_id, started)
