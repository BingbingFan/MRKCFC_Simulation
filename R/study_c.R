# Shared machinery for Studies C2 to C5.
#
# Sourced by study_c2.R, study_c3.R, study_c4.R and study_c5.R. It exists so
# there is ONE definition of the C configurations, the config override, and the
# result-writing convention, rather than four copies that drift -- the failure
# that gave silhouette_score() a divergent twin.

suppressMessages(library(MRKCFC))
source(file.path("R", "generate.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "methods.R"))
source(file.path("R", "separation.R"))
for (f in c("gen_scores", "vm_intermediate", "rmvn", "FLEISHMAN", "IQR_SCALE"))
  assign(f, get(f, envir = asNamespace("MRKCFC")), envir = globalenv())

# The two cells every C sub-study runs at, from config/scenarios_c.csv:
# a well-behaved reference and a deliberately hard one. C4 is the exception --
# it runs across all twelve Study A scenarios, because non-convergence is
# arm-specific and the point is to find out whether it really is confined to the
# heavy tails.
load_c_configs <- function()
  utils::read.csv(file.path("config", "scenarios_c.csv"), stringsAsFactors = FALSE)

# One row of scenarios_c.csv, shaped for generate_replicate().
as_scenario_c <- function(row) data.frame(
  scenario_id = row$config_id, arm = row$arm, d = row$d,
  arm_eps = row$arm_eps, arm_k = row$arm_k, arm_nu = row$arm_nu,
  stringsAsFactors = FALSE)

apply_config_c <- function(cfg, row) {
  cfg$global$obs_poisson_lambda <- row$obs_poisson_lambda
  cfg
}

# Study C's replicate count. 50, decided 2026-08-28: C4's headline output is a
# RATE, and 50 estimates one to about +/- 7 points where 10 would give +/- 15
# and could not separate 0.7 from 0.5.
c_reps <- function(cfg, override = NA_integer_)
  if (!is.na(override)) override else as.integer(cfg$global$n_reps_c)

# Pin the RNG state before EVERY fit.
#
# study_c2.R and study_c3.R call MKCFC() and METHODS[[mid]] directly rather than
# through run_all_methods(), so they do not inherit its re-seeding and would
# otherwise repeat the defect fixed on 2026-08-28: a method's result depending
# on which methods ran before it. `two-staged` initialisation draws on the RNG
# twice, for k-means centres and for the random fallback, so this is not
# hypothetical.
#
# Keyed on the CELL only -- configuration and replicate -- deliberately, not on
# the method or on the axis being varied. Every method, every M, every basis
# setting and every K within one replicate then starts from an identical state,
# which pairs the comparisons C2 and C3 exist to make and removes the RNG as a
# source of difference between the arms.
c_seed_fit <- function(cfg, config_id, rep) {
  set.seed(cell_seed(cfg$global$seed_master, paste0(config_id, "|fit"), rep))
  invisible(NULL)
}

# Same contract as run_cell.R: a contiguous slice, exactly one file, no
# contention between tasks.
c_slice <- function(grid, task_id, n_tasks) {
  mine <- which((seq_len(nrow(grid)) - 1L) %% n_tasks + 1L == task_id)
  grid[mine, , drop = FALSE]
}

c_write <- function(results, outdir, stem, task_id, started) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  results$task_id <- task_id
  attr(results, "session") <- list(
    R = as.character(getRversion()),
    MRKCFC = as.character(utils::packageVersion("MRKCFC")),
    host = Sys.info()[["nodename"]], started = started, finished = Sys.time())
  f <- file.path(outdir, sprintf("%s_task_%04d.rds", stem, task_id))
  saveRDS(results, f)
  cat(sprintf("\nwrote %s  (%d rows, %.1f min)\n", f, nrow(results),
              as.numeric(difftime(Sys.time(), started, units = "mins"))))
  invisible(f)
}

c_header <- function(what, mine, total, task_id, n_tasks) {
  cat(sprintf("%s | task %d of %d: %d cells of %d\n",
              what, task_id, n_tasks, mine, total))
  cat(sprintf("R %s | MRKCFC %s | host %s\n", getRversion(),
              utils::packageVersion("MRKCFC"), Sys.info()[["nodename"]]))
}
