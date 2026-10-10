# Study C1: initialisation sensitivity.
#
#   Rscript R/study_c1.R <scenario_ids> <methods> <reps> <n_restarts> [outdir]
#   Rscript R/study_c1.R A11,A12 M2 2 5 results_c1
#
# The design: k-means on global MRFPCA scores (the default), 20 random
# restarts, a Ward hierarchical start, and an oracle true-label start; report
# the aRI distribution across starts and the pairwise agreement among the
# partitions each start produces.
#
# The algorithm is iterative and may depend on the starting labels, so the
# sensitivity to initialisation is measured directly.
#
# Two outputs, matching the two things the plan asks for:
#   1. aRI to truth for each strategy, plus the spread ACROSS the random
#      restarts, which is the quantity that says how much the start matters
#   2. pairwise aRI BETWEEN the partitions the strategies produce, which says
#      whether they are finding the same solution or different ones
#
# The restart arm selects by the algorithm's own objective, total squared
# reconstruction error, not by aRI. Selecting by aRI would use the truth to
# choose, which is unavailable in practice and would overstate the method.

suppressMessages(library(MRKCFC))
source(file.path("R", "generate.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "methods.R"))
for (f in c("gen_scores", "vm_intermediate", "rmvn", "FLEISHMAN", "IQR_SCALE"))
  assign(f, get(f, envir = asNamespace("MRKCFC")), envir = globalenv())

args      <- commandArgs(trailingOnly = TRUE)
scen_ids  <- if (length(args) >= 1) strsplit(args[1], ",")[[1]] else c("A11", "A12")
meth_ids  <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else c("M1", "M2", "M3")
n_reps    <- if (length(args) >= 3) as.integer(args[3]) else 2L
n_restart <- if (length(args) >= 4) as.integer(args[4]) else 20L
outdir    <- if (length(args) >= 5) args[5] else "results_c1"
# TASK SLICING. Added 2026-09-04, then REWRITTEN 2026-09-06 to slice on STARTS
# rather than on cells -- see the long note below the helpers. A task now takes
# a share of the (scenario, replicate, method, start) grid, so its cost is a few
# fits rather than 23. Omit both arguments and the behaviour is unchanged: one
# process, one file, the whole grid.
task_id   <- if (length(args) >= 6) as.integer(args[6]) else NA_integer_
n_tasks   <- if (length(args) >= 7) as.integer(args[7]) else NA_integer_
sliced    <- !is.na(task_id) && !is.na(n_tasks) && n_tasks > 1L
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cfg <- load_config()
K <- as.integer(cfg$global$K_true)

# method id -> the (robust, scatter) pair it fits with
mcfg <- list(M1 = list(TRUE, "moment"), M2 = list(TRUE, "sign"), M3 = list(FALSE, "moment"))

fit_with <- function(Ly, Lt, robust, scatter, init.method = "kmeans", init.labels = NULL) {
  g <- cfg$global
  MRKCFC::MKCFC(Ly = Ly, Lt = Lt, K = K, robust = robust,
               optns = build_optns(cfg, robust), mfpca_scatter = scatter,
               init.method = init.method, init.labels = init.labels,
               maxIter = as.integer(g$maxIter), parallel = FALSE, trace = 0,
               conv.prop = as.numeric(g$conv_prop),
               LOO = isTRUE(as.logical(g$LOO)))
}

# The algorithm's own objective: total squared reconstruction error under the
# partition. Lower is better, and it uses no knowledge of the truth.
objective <- function(Ly, Lt, labels) {
  r <- tryCatch(reconstruct_given_partition(Ly, Lt, labels, cfg), error = function(e) NULL)
  if (is.null(r)) return(NA_real_)
  sum((r$y - r$yhat)^2)
}

ward_start <- function(Ly, Lt, robust, scatter) {
  sc <- global_mfpca_scores(Ly, Lt, cfg, robust)
  stats::cutree(stats::hclust(stats::dist(sc), method = "ward.D2"), k = K)
}
# THE UNIT OF WORK IS ONE START, NOT ONE CELL. Changed 2026-09-06 after job
# 3286698 timed out: 13 of 20 tasks hit a 24 h wall with one cell each, and a
# cell is 23 full fits. One A01 cell -- the CHEAP arm -- did not finish in 24 h,
# so a cell is >23 h and 48 h, the partition ceiling, was an even-odds bet with
# no third attempt. The 23 starts within a cell are INDEPENDENT, so slicing on
# them makes a task ~6 fits instead of 23 and puts the wall out of reach.
#
# TWO CONSEQUENCES, both handled:
#
#   `best-of-N` MOVED TO THE COMBINER. It needs all 20 restarts and no single
#   task has them any more. `combine_c1.R` computes it from the pooled rows,
#   which is where it belonged: the driver was selecting a winner from data it
#   happened to hold.
#
#   EVERY START IS SEEDED DETERMINISTICALLY, from (scenario, replicate, method,
#   start). Previously the 20 restarts drew in sequence from one RNG state, so
#   splitting them across tasks would have changed the results -- and changed
#   them in a way that depended on the slicing, which is the worst kind of
#   irreproducibility. Now a start returns the same fit wherever it runs.
STARTS <- c("kmeans", "ward", "oracle", sprintf("random%02d", seq_len(n_restart)))

grid <- expand.grid(start = STARTS, method = meth_ids, rep = seq_len(n_reps),
                    scenario_id = scen_ids, stringsAsFactors = FALSE)
grid <- grid[order(grid$scenario_id, grid$rep, grid$method,
                   match(grid$start, STARTS)), ]
rownames(grid) <- NULL
n_all <- nrow(grid)
if (sliced) {
  grid <- grid[which(((seq_len(n_all) - 1L) %% n_tasks) + 1L == task_id), ,
               drop = FALSE]
  cat(sprintf("task %d of %d: %d of %d fits\n", task_id, n_tasks, nrow(grid), n_all))
  if (!nrow(grid)) { cat("nothing to do\n"); quit(save = "no", status = 0) }
}

# Deterministic per-start seed. Distinct across every (scenario, rep, method,
# start) and independent of how the grid is sliced.
start_seed <- function(sid, r, mid, tag)
  20260906L + 7919L * match(sid, scen_ids) + 131L * as.integer(r) +
    17L * match(mid, meth_ids) + match(tag, STARTS)

rows <- list(); parts <- list()
cur_key <- ""; dat <- NULL
for (gi in seq_len(nrow(grid))) {
  sid <- grid$scenario_id[gi]; r <- grid$rep[gi]
  mid <- grid$method[gi];      tag <- grid$start[gi]

  # regenerate only when the replicate changes; the grid is ordered so that is
  # once per cell rather than once per fit
  key <- paste(sid, r, sep = "|")
  if (!identical(key, cur_key)) {
    sc  <- cfg$scenarios[cfg$scenarios$scenario_id == sid, ]
    dat <- generate_replicate(cfg, sc, rep = r)
    cur_key <- key
  }
  rob <- mcfg[[mid]][[1]]; scat <- mcfg[[mid]][[2]]

  one <- function(tag, expr) {
    t0 <- Sys.time()
    f <- tryCatch(expr, error = function(e) conditionMessage(e))
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (is.character(f))
      return(list(row = data.frame(scenario_id = sid, rep = r, method = mid,
        start = tag, aRI = NA_real_, CCR = NA_real_, obj = NA_real_,
        converged = NA, seconds = el, error = f, stringsAsFactors = FALSE),
        labels = NULL))
    lab <- as.integer(f$res_clust$g)
    list(row = data.frame(scenario_id = sid, rep = r, method = mid, start = tag,
      aRI = ari(dat$labels, lab), CCR = ccr(dat$labels, lab),
      obj = objective(dat$Ly, dat$Lt, lab), converged = isTRUE(f$conv),
      seconds = el, error = NA_character_, stringsAsFactors = FALSE),
      labels = lab)
  }

  set.seed(start_seed(sid, r, mid, tag))
  z <- if (tag == "kmeans") {
    one(tag, fit_with(dat$Ly, dat$Lt, rob, scat, "kmeans"))
  } else if (tag == "ward") {
    one(tag, fit_with(dat$Ly, dat$Lt, rob, scat,
                      init.labels = ward_start(dat$Ly, dat$Lt, rob, scat)))
  } else if (tag == "oracle") {
    one(tag, fit_with(dat$Ly, dat$Lt, rob, scat, init.labels = dat$labels))
  } else {
    one(tag, fit_with(dat$Ly, dat$Lt, rob, scat, "random"))
  }

  rows[[length(rows) + 1L]] <- z$row
  if (!is.null(z$labels)) parts[[paste(sid, r, mid, tag, sep = "|")]] <- z$labels

  # PER-FIT PROGRESS. The previous version printed once per cell, so a 24 h task
  # emitted one line at the start and nothing until it wrote -- there was no way
  # to tell how far it had got, or to size the next run from a part-finished one.
  cat(sprintf("  [%3d/%3d] %-5s rep %2d %-3s %-10s %7.1fs%s\n", gi, nrow(grid),
              sid, r, mid, tag, z$row$seconds,
              if (!is.na(z$row$error)) "  ERROR" else ""))
  utils::flush.console()
}

all <- do.call(rbind, rows)

# A DRIVER MUST NOT SUMMARISE. Under slicing the tables below would be computed
# from one task's cells and would print something that looks complete and is
# mostly empty -- the failure `study_fve.R` and `study_estim.R` both had. When
# sliced, write only; `combine_c1.R` is the one reader.
if (sliced) {
  f <- file.path(outdir, sprintf("c1_task_%04d.rds", task_id))
  saveRDS(list(results = all, partitions = parts), f)
  cat(sprintf("\nwrote %s (%d rows)\n", f, nrow(all)))
  quit(save = "no", status = 0)
}

saveRDS(list(results = all, partitions = parts), file.path(outdir, "c1.rds"))
utils::write.csv(all, file.path(outdir, "c1_raw.csv"), row.names = FALSE)

cat("\n=== aRI to truth by starting strategy ===\n")
main <- all[!grepl("^random[0-9]", all$start), ]
w <- stats::reshape(main[, c("scenario_id", "rep", "method", "start", "aRI")],
                    idvar = c("scenario_id", "rep", "method"),
                    timevar = "start", direction = "wide")
names(w) <- sub("^aRI[.]", "", names(w))
print(w, row.names = FALSE, digits = 3)

cat("\n=== spread ACROSS random restarts (the quantity that says how much the start matters) ===\n")
rr <- all[grepl("^random[0-9]", all$start), ]
if (nrow(rr)) {
  s <- do.call(rbind, lapply(split(rr, list(rr$scenario_id, rr$rep, rr$method), drop = TRUE),
    function(d) data.frame(scenario_id = d$scenario_id[1], rep = d$rep[1], method = d$method[1],
      n_ok = sum(!is.na(d$aRI)), min = min(d$aRI, na.rm = TRUE),
      median = stats::median(d$aRI, na.rm = TRUE), max = max(d$aRI, na.rm = TRUE),
      stringsAsFactors = FALSE)))
  print(s, row.names = FALSE, digits = 3)
}

cat("\n=== pairwise agreement between the partitions each start produced ===\n")
for (key in unique(sub("[|][^|]*$", "", names(parts)))) {
  sel <- parts[startsWith(names(parts), paste0(key, "|"))]
  sel <- sel[!grepl("[|]random[0-9]", names(sel))]
  if (length(sel) < 2L) next
  tags <- sub("^.*[|]", "", names(sel))
  m <- outer(seq_along(sel), seq_along(sel), Vectorize(function(i, j) ari(sel[[i]], sel[[j]])))
  dimnames(m) <- list(tags, tags)
  cat("\n", key, "\n", sep = ""); print(round(m, 3))
}
cat(sprintf("\nwrote %s\n", file.path(outdir, "c1.rds")))
