# One array task. Computes a contiguous slice of the Study A grid and writes
# exactly one result file, so tasks never contend for a write.
#
#   Rscript R/run_cell.R <task_id> [n_tasks] [outdir]
#
# In an array job the task id comes from the scheduler, for instance
#
#   Rscript R/run_cell.R $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT results
#
# The identical code path runs locally with an explicit id, so the only surface
# that differs between a workstation and a cluster is the scheduler wrapper. That
# matters: the clusterExport defect earlier in this project existed only in
# installed mode and was invisible when sourced. An execution mode that is
# never tested hides defects.
#
# Slicing is over (scenario x replicate) cells. Each cell's seed is derived
# from the master seed and the cell identity, so a task reproduces its own
# slice regardless of execution order or which tasks are re-run after failure.

suppressMessages(library(MRKCFC))
source(file.path("R", "generate.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "methods.R"))
source(file.path("R", "separation.R"))
for (f in c("gen_scores", "vm_intermediate", "rmvn", "FLEISHMAN", "IQR_SCALE"))
  assign(f, get(f, envir = asNamespace("MRKCFC")), envir = globalenv())

args     <- commandArgs(trailingOnly = TRUE)
task_id  <- if (length(args) >= 1) as.integer(args[1]) else 1L
n_tasks  <- if (length(args) >= 2) as.integer(args[2]) else 1L
outdir   <- if (length(args) >= 3) args[3] else "results"
# Optional replicate override, for exploratory rounds. Deliberately an argument
# rather than a config edit: n_reps in the config is the REPORTABLE value and
# should not be mutated to run a quick check.
rep_over <- if (length(args) >= 4) as.integer(args[4]) else NA_integer_
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cfg <- load_config()
n_reps <- if (is.na(rep_over)) as.integer(cfg$global$n_reps) else rep_over

# Methods are selected in the config, not here. A method switched off still has
# a row in methods.csv with include = FALSE and the reason, so an exclusion is
# visible and reversible rather than silently absent.
use_methods <- names(METHODS)
if ("include" %in% names(cfg$methods)) {
  keep <- cfg$methods$method_id[toupper(as.character(cfg$methods$include)) == "TRUE"]
  use_methods <- intersect(names(METHODS), keep)
  dropped <- setdiff(names(METHODS), use_methods)
  if (length(dropped))
    cat(sprintf("excluded by config: %s\n", paste(
      cfg$methods$label[match(dropped, cfg$methods$method_id)], collapse = ", ")))
}

# the full grid, one row per cell, in a fixed order
grid <- expand.grid(rep = seq_len(n_reps),
                    scenario_id = cfg$scenarios$scenario_id,
                    stringsAsFactors = FALSE)
grid <- grid[order(grid$scenario_id, grid$rep), ]
rownames(grid) <- NULL

# contiguous slice for this task
mine <- which((seq_len(nrow(grid)) - 1L) %% n_tasks + 1L == task_id)
if (!length(mine)) {
  cat(sprintf("task %d of %d: no cells, nothing to do\n", task_id, n_tasks))
  quit(save = "no", status = 0)
}

cat(sprintf("task %d of %d: %d cells of %d total\n",
            task_id, n_tasks, length(mine), nrow(grid)))
cat(sprintf("R %s | MRKCFC %s | host %s\n",
            getRversion(), utils::packageVersion("MRKCFC"), Sys.info()[["nodename"]]))

started <- Sys.time()
out <- vector("list", length(mine))
for (a in seq_along(mine)) {
  cell <- grid[mine[a], ]
  sc <- cfg$scenarios[cfg$scenarios$scenario_id == cell$scenario_id, ]
  t0 <- Sys.time()
  res <- tryCatch({
    dat <- generate_replicate(cfg, sc, rep = cell$rep)
    # NOT named `out`: tryCatch evaluates in the calling frame, so that would
    # clobber the accumulator list this loop is filling.
    cell_rows <- run_all_methods(dat, cfg, methods = use_methods)
    # A property of the CELL, not of any method, so the same value lands on every
    # row. Carried with the results so a CCR can be read against the Bayes
    # optimum at that separation rather than in a vacuum.
    sep <- separation_of(dat, cfg)
    cell_rows$delta <- sep$delta; cell_rows$bayes_ccr <- sep$bayes
    cell_rows
  }, error = function(e) {
    # a whole cell failing must not take down the task
    data.frame(scenario_id = cell$scenario_id, rep = cell$rep, method = NA_character_,
               CCR = NA_real_, aRI = NA_real_, NMI = NA_real_, AMI = NA_real_,
               MPE = NA_real_, RRE = NA_real_, RRE_trim = NA_real_,
               # `iter` IS REQUIRED. run_all_methods() returns it on every row,
               # so omitting it here made this frame 15 columns against the
               # success path's 16 and do.call(rbind, out) would abort the WHOLE
               # task on the first failed cell -- discarding every cell already
               # completed, and reporting an rbind error instead of the real
               # cause. That is the exact failure the `error` column exists to
               # prevent. Fixed 2026-09-04; never seen to fire.
               converged = NA, iter = NA_integer_, seconds = NA_real_,
               error = paste("cell failed:", conditionMessage(e)),
               delta = NA_real_, bayes_ccr = NA_real_,
               stringsAsFactors = FALSE)
  })
  out[[a]] <- res
  cat(sprintf("  [%3d/%3d] %s rep %3d  %6.1fs%s\n", a, length(mine),
              cell$scenario_id, cell$rep,
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              if (any(!is.na(res$error))) sprintf("  (%d errors)", sum(!is.na(res$error))) else ""))
  utils::flush.console()
}

results <- do.call(rbind, out)
results$task_id <- task_id
attr(results, "config") <- cfg$global
attr(results, "session") <- list(
  R = as.character(getRversion()),
  MRKCFC = as.character(utils::packageVersion("MRKCFC")),
  host = Sys.info()[["nodename"]], started = started, finished = Sys.time())

f <- file.path(outdir, sprintf("studyA_task_%04d.rds", task_id))
saveRDS(results, f)
cat(sprintf("\nwrote %s  (%d rows, %.1f min)\n", f, nrow(results),
            as.numeric(difftime(Sys.time(), started, units = "mins"))))
