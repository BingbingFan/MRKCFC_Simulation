# Combines the per-task result files into one tidy table plus the summaries
# Study A reports.
#
#   Rscript R/combine.R [resultsdir-or-rawcsv] [outstem]
#
# The first argument is either a DIRECTORY of per-task .rds files, or an
# already-assembled RAW CSV. The second form exists so that a raw file
# assembled elsewhere can be summarised by this script rather than by a second
# copy of the summary logic. One definition, as with everything else here.
#
# Reports medians with interquartile ranges rather than means with standard
# deviations, deliberately: arms 5 and 6 have infinite excess kurtosis and arm
# 6 has infinite variance, so a mean and an sd are not defined quantities
# across the design. That is the same reason the arms are IQR-matched rather
# than variance-matched.

args   <- commandArgs(trailingOnly = TRUE)
resdir <- if (length(args) >= 1) args[1] else "results"
stem   <- if (length(args) >= 2) args[2] else file.path(resdir, "studyA")

if (grepl("[.]csv$", resdir, ignore.case = TRUE)) {
  if (!file.exists(resdir)) stop("no such file: ", resdir, call. = FALSE)
  all <- utils::read.csv(resdir, stringsAsFactors = FALSE)
  cat(sprintf("reading assembled raw file %s\n", resdir))
  if (length(args) < 2)
    stem <- sub("_raw[.]csv$|[.]csv$", "", resdir, ignore.case = TRUE)
  # A merged file already carries method labels and scenario labels from the
  # run that produced it, so the relabelling below must not be applied twice.
  # It is idempotent by construction -- an unmatched id keeps its own value --
  # but the method column would be blanked, so the ids are restored first.
  if ("method_id" %in% names(all)) all$method <- all$method_id
  if ("scenario_id_code" %in% names(all)) all$scenario_id <- all$scenario_id_code
} else {
  files <- list.files(resdir, pattern = "^study[AB]_task_.*[.]rds$", full.names = TRUE)
  if (!length(files)) stop("no result files found in ", resdir, call. = FALSE)
  cat(sprintf("combining %d task files\n", length(files)))
  all <- do.call(rbind, lapply(files, readRDS))
}
cat(sprintf("%d rows, %d scenarios, %d methods\n", nrow(all),
            length(unique(all$scenario_id)), length(unique(na.omit(all$method)))))

# Report NAMES, not codes. Nobody should have to hold a mapping from M1..M8 and
# A01..A12 in their head to read a results table, least of all months later.
.cfgdir <- if (dir.exists("config")) "config" else file.path("..", "config")
.mlab <- tryCatch(utils::read.csv(file.path(.cfgdir, "methods.csv"),
                                  stringsAsFactors = FALSE), error = function(e) NULL)
.slab <- tryCatch(utils::read.csv(file.path(.cfgdir, "scenarios.csv"),
                                  stringsAsFactors = FALSE), error = function(e) NULL)
# Study B's configurations live in their own table. Both maps are consulted so
# one combine step serves A and B; a duplicated combine_b.R would drift from
# this one, which is exactly what happened to silhouette_score().
.blab <- tryCatch(utils::read.csv(file.path(.cfgdir, "scenarios_b.csv"),
                                  stringsAsFactors = FALSE), error = function(e) NULL)
if (!is.null(.mlab)) {
  all$method_id <- all$method
  all$method <- .mlab$label[match(all$method, .mlab$method_id)]
}
.map <- character(0)
if (!is.null(.slab))
  .map <- c(.map, stats::setNames(sprintf("%s d=%g", .slab$arm_label, .slab$d),
                                  .slab$scenario_id))
if (!is.null(.blab))
  .map <- c(.map, stats::setNames(.blab$label, .blab$config_id))
if (length(.map)) {
  all$scenario_id_code <- all$scenario_id
  hit <- unname(.map[all$scenario_id])
  # An id absent from both tables keeps its own code rather than becoming NA,
  # which would silently drop it from every summary below.
  all$scenario_id <- ifelse(is.na(hit), all$scenario_id_code, hit)
}

# --- completeness, before anything is summarised -----------------------------
cfg_reps <- length(unique(all$rep))
expected <- length(unique(all$scenario_id)) * cfg_reps * length(unique(na.omit(all$method)))
cat(sprintf("\ncompleteness: %d of %d expected rows (%.1f%%)\n",
            sum(!is.na(all$method)), expected, 100 * sum(!is.na(all$method)) / expected))
errs <- all[!is.na(all$error), ]
if (nrow(errs)) {
  cat(sprintf("\n%d failed fits:\n", nrow(errs)))
  tb <- table(errs$method, substr(errs$error, 1, 60))
  print(tb)
} else cat("no failed fits\n")

nonconv <- stats::aggregate(converged ~ method, all,
                            function(x) mean(!x, na.rm = TRUE), na.action = na.omit)
names(nonconv)[2] <- "nonconvergence_rate"
cat("\nnon-convergence by method:\n"); print(nonconv, row.names = FALSE)

# --- summaries ---------------------------------------------------------------
med_iqr <- function(x) sprintf("%.3f (%.3f)", stats::median(x, na.rm = TRUE),
                               stats::IQR(x, na.rm = TRUE))

# `iter` joined the results on 2026-08-28 and is summarised alongside the rest,
# because it is what makes a non-convergence rate actionable: iter == maxIter
# means the cap bound the fit, iter < maxIter with converged = FALSE means it
# broke early on the <10-curves rule, and those need opposite remedies.
metrics <- c("CCR", "aRI", "NMI", "AMI", "MPE", "RRE", "RRE_trim",
             "seconds", "iter")
metrics <- metrics[metrics %in% names(all)]   # older result files lack `iter`
summ <- do.call(rbind, lapply(split(all, list(all$scenario_id, all$method), drop = TRUE),
  function(d) {
    r <- data.frame(scenario_id = d$scenario_id[1], method = d$method[1],
                    n_fits = sum(!is.na(d$CCR)), stringsAsFactors = FALSE)
    # Convergence belongs next to the metrics, not in a separate table: a median
    # aRI read without it is the mistake the 10-replicate round invited, where
    # MRKCFC-R2 looked worst on arm 6 purely because it converged least often.
    r$conv_rate <- if ("converged" %in% names(d))
      sprintf("%.2f", mean(d$converged, na.rm = TRUE)) else NA_character_
    for (m in metrics) r[[m]] <- med_iqr(d[[m]])
    r
  }))
rownames(summ) <- NULL
summ <- summ[order(summ$scenario_id, summ$method), ]

utils::write.csv(all,  paste0(stem, "_raw.csv"),     row.names = FALSE)
utils::write.csv(summ, paste0(stem, "_summary.csv"), row.names = FALSE)
cat(sprintf("\nwrote %s_raw.csv and %s_summary.csv\n", stem, stem))

# --- the headline table ------------------------------------------------------
cat("\n=== median aRI by scenario and method ===\n")
w <- stats::reshape(
  stats::aggregate(aRI ~ scenario_id + method, all, stats::median, na.action = na.omit),
  idvar = "scenario_id", timevar = "method", direction = "wide")
names(w) <- sub("^aRI[.]", "", names(w))
print(w, row.names = FALSE, digits = 3)

# --- separation, and what it implies ----------------------------------------
if ("delta" %in% names(all)) {
  cat("\n=== robust Mahalanobis separation, and the CCR it makes attainable ===\n")
  sep <- do.call(rbind, lapply(split(all, all$scenario_id), function(d)
    data.frame(scenario_id = d$scenario_id[1],
               delta = stats::median(d$delta, na.rm = TRUE),
               ref_ccr = stats::median(d$bayes_ccr, na.rm = TRUE),
               best_observed_CCR = max(d$CCR, na.rm = TRUE),
               stringsAsFactors = FALSE)))
  print(sep, row.names = FALSE, digits = 3)
  cat("\nref_ccr is Phi(delta/2): what the optimal LINEAR rule achieves from\n")
  cat("centroid separation alone. It is a YARDSTICK, not a bound. This design\n")
  cat("gives variables 1 and 2 class-specific eigenfunctions, so the classes\n")
  cat("differ in covariance as well as location, the true Bayes rule is\n")
  cat("quadratic, and methods can legitimately exceed ref_ccr by exploiting the\n")
  cat("shape difference; 107 to 111 percent has been observed. Read it both ways:\n")
  cat("well above means a method uses more than the centroid gap, well below\n")
  cat("means the cell is a floor. It answers R2's request for a reference value.\n")
  cat("A scenario where every method sits far above it, at CCR 1.000, cannot rank\n")
  cat("methods however many replicates are run.\n")

  cat("\n=== does the IQR-matching hold? delta should be flat ACROSS arms at fixed d ===\n")
  if (all(c("arm", "d") %in% names(cfg_scen <- tryCatch(
        utils::read.csv(file.path("config", "scenarios.csv"), stringsAsFactors = FALSE),
        error = function(e) NULL)))) {
    sep$scenario_id_code <- all$scenario_id_code[match(sep$scenario_id, all$scenario_id)]
    m <- merge(sep, cfg_scen[, c("scenario_id", "arm", "d")],
               by.x = "scenario_id_code", by.y = "scenario_id")
    w <- tapply(m$delta, list(m$d, m$arm), function(x) x[1])
    print(round(w, 2))
    cat("\nThe arms are matched on IQR so that a reader cannot attribute a\n")
    cat("difference to signal-to-noise rather than to tails. A row that is flat\n")
    cat("is evidence the matching worked; a row that is not means tail effects\n")
    cat("are confounded with separation, which is the objection the matching\n")
    cat("exists to foreclose.\n")
  }
}

cat("\nRead these against the moment ladder: arms 1 to 3 have all moments finite,\n")
cat("arm 4 is contaminated, arm 5 has infinite fourth moment, arm 6 has infinite\n")
cat("variance. The revision's claim is that the robust second stage removes a\n")
cat("failure mode down that ladder, at a small efficiency cost at the top of it.\n")
