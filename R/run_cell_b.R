# One array task of Study B. Same slicing contract as run_cell.R: compute a
# contiguous slice, write exactly one result file, never contend for a write.
#
#   Rscript R/run_cell_b.R <task_id> [n_tasks] [outdir] [n_reps_override]
#
# In an array job:
#   Rscript R/run_cell_b.R $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT results_b
#
# WHAT STUDY B VARIES, AND WHY IT IS SHAPED LIKE THIS
#
# One factor at a time from S0, which is arm 4 (contaminated normal) at d = 1.5,
# n = 500, P = 3, pi = 0.5/0.5, Poisson(7). Twelve configurations, decided
# 2026-08-28; see scenarios_b.csv, whose `axis` column records which factor each
# row moves.
#
# S0 sits at d = 1.5 deliberately. At an easier baseline a complexity effect can
# hide inside a ceiling, which is the failure that made the original d = {5, 8}
# design uninformative. d = 1.5 has measured Delta about 1.73 and realised CCR
# 0.77 to 0.87, so there is headroom to lose.
#
# THE P AXIS IS THE ONE THAT NEEDED CARE. Delta is a Mahalanobis distance in the
# P * M dimensional score space, so if every added variable carried the same mean
# shift then Delta would grow like sqrt(P) -- P = 3 to P = 10 multiplies it by
# about 1.83 -- and the axis would measure "more information is easier" rather
# than "more dimensions are harder". Variables 4 to 10 are therefore
# UNINFORMATIVE by construction: class-common mean with no d term, and
# class-common eigenfunctions, so they carry no class information through either
# location or shape. validate_config.R asserts this, because if it ever stops
# holding the axis silently changes meaning. `delta` is recorded per cell so the
# assumption is checkable in the output rather than merely asserted here.

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
outdir   <- if (length(args) >= 3) args[3] else "results_b"
rep_over <- if (length(args) >= 4) as.integer(args[4]) else NA_integer_
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cfg <- load_config()
configs <- utils::read.csv(file.path("config", "scenarios_b.csv"),
                           stringsAsFactors = FALSE)
n_reps <- if (!is.na(rep_over)) rep_over else as.integer(cfg$global$n_reps_b)

# Methods come from methods.csv's include_b column, not from this file. B runs
# MKCFC, MRKCFC (moment), MRKCFC-R2 and teigen: the three-way ablation plus
# one model-based comparator. The score-then-cluster comparators and lcmm are
# switched off there with a reason, so the exclusion is visible and reversible.
stopifnot("include_b" %in% names(cfg$methods))
keep <- cfg$methods$method_id[toupper(as.character(cfg$methods$include_b)) == "TRUE"]
use_methods <- intersect(names(METHODS), keep)
stopifnot(length(use_methods) > 0L)
cat(sprintf("methods: %s\n", paste(
  cfg$methods$label[match(use_methods, cfg$methods$method_id)], collapse = ", ")))

# apply_config_b() and as_scenario_b() are defined in R/generate.R, not here,
# because test_study_b.R needs them too and a second copy would drift.

grid <- expand.grid(rep = seq_len(n_reps), config_id = configs$config_id,
                    stringsAsFactors = FALSE)
grid <- grid[order(grid$config_id, grid$rep), ]
rownames(grid) <- NULL

mine <- which((seq_len(nrow(grid)) - 1L) %% n_tasks + 1L == task_id)
if (!length(mine)) {
  cat(sprintf("task %d of %d: no cells, nothing to do\n", task_id, n_tasks))
  quit(save = "no", status = 0)
}

cat(sprintf("task %d of %d: %d cells of %d total (%d configs x %d reps)\n",
            task_id, n_tasks, length(mine), nrow(grid), nrow(configs), n_reps))
cat(sprintf("R %s | MRKCFC %s | host %s\n", getRversion(),
            utils::packageVersion("MRKCFC"), Sys.info()[["nodename"]]))

started <- Sys.time()
out <- vector("list", length(mine))
for (a in seq_along(mine)) {
  cell <- grid[mine[a], ]
  row  <- configs[configs$config_id == cell$config_id, ]
  t0 <- Sys.time()
  res <- tryCatch({
    c2  <- apply_config_b(cfg, row)
    dat <- generate_replicate(c2, as_scenario_b(row), rep = cell$rep)
    cell_rows <- run_all_methods(dat, c2, methods = use_methods)
    # Delta is recorded per cell for two reasons here. It supplies the reference
    # CCR as in Study A, and on the P axis it is the CHECK that added variables
    # really are uninformative: if delta drifts upward with P, the axis is
    # measuring added signal and the result must not be read as a complexity
    # effect.
    sep <- separation_of(dat, c2)
    cell_rows$delta <- sep$delta; cell_rows$bayes_ccr <- sep$bayes
    cell_rows
  }, error = function(e) data.frame(
    scenario_id = cell$config_id, rep = cell$rep, method = NA_character_,
    CCR = NA_real_, aRI = NA_real_, NMI = NA_real_, AMI = NA_real_,
    MPE = NA_real_, RRE = NA_real_, RRE_trim = NA_real_,
    converged = NA, iter = NA_integer_, seconds = NA_real_,
    error = paste("cell failed:", conditionMessage(e)),
    delta = NA_real_, bayes_ccr = NA_real_, stringsAsFactors = FALSE))

  # Carried on every row so a result can be read without joining back to the
  # config table, and so `axis` can group the reporting directly.
  res$config_id <- row$config_id
  res$config_label <- row$label
  res$axis <- row$axis
  res$n <- row$n; res$P <- row$P; res$pi <- row$pi
  res$lambda <- row$obs_poisson_lambda; res$d <- row$d
  out[[a]] <- res
  cat(sprintf("  [%3d/%3d] %-4s %-14s rep %3d  %6.1fs%s\n", a, length(mine),
              row$config_id, row$label, cell$rep,
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              if (any(!is.na(res$error))) sprintf("  (%d errors)", sum(!is.na(res$error))) else ""))
  utils::flush.console()
}

results <- do.call(rbind, out)
results$task_id <- task_id
# `method` is left as the method_id on purpose. combine.R maps ids to labels for
# both studies from methods.csv, and relabelling here as well would make that
# match fail and blank the column.
attr(results, "config") <- cfg$global
attr(results, "session") <- list(
  R = as.character(getRversion()),
  MRKCFC = as.character(utils::packageVersion("MRKCFC")),
  host = Sys.info()[["nodename"]], started = started, finished = Sys.time())

f <- file.path(outdir, sprintf("studyB_task_%04d.rds", task_id))
saveRDS(results, f)
cat(sprintf("\nwrote %s  (%d rows, %.1f min)\n", f, nrow(results),
            as.numeric(difftime(Sys.time(), started, units = "mins"))))
