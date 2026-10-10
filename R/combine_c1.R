# Read the C1 initialisation study. The one reader; the driver summarises
# nothing when it is sliced.
#
#   Rscript R/combine_c1.R <dir-of-c1_task_*.rds> [outstem]
#
# The outstem is a PATH, like combine_estim.R's and unlike combine_c.R's, whose
# relative stem has stranded output in code/study more than once.

args <- commandArgs(trailingOnly = TRUE)
inp  <- if (length(args) >= 1) args[1] else "results_c1"
stem <- if (length(args) >= 2) args[2] else file.path(inp, "c1")

fs <- list.files(inp, pattern = "^c1_task_[0-9]+[.]rds$", full.names = TRUE)
if (!length(fs)) stop("no c1_task_*.rds under ", inp)
objs <- lapply(fs, readRDS)
all  <- do.call(rbind, lapply(objs, `[[`, "results"))
parts <- do.call(c, lapply(objs, `[[`, "partitions"))
cat(sprintf("read %d task files, %d rows, %d partitions\n",
            length(fs), nrow(all), length(parts)))

# BEST-OF-N IS COMPUTED HERE, not in the driver. Since 2026-09-06 a task holds
# only a share of the 20 restarts, so no task can select among them -- and this
# is where the selection belonged anyway, rather than in whichever process
# happened to hold the rows.
#
# Selection is BY THE ALGORITHM'S OWN OBJECTIVE, never by aRI. Selecting by aRI
# would use the truth to choose, which a practitioner cannot do and which would
# overstate the method.
#
# Any `best-of-*` rows already present are DROPPED first. Files written before
# 2026-09-06 carry one computed inside the driver, and keeping both would double
# count the winner.
all <- all[!grepl("^best-of-", all$start), ]
rr0 <- all[grepl("^random[0-9]", all$start), ]
if (nrow(rr0)) {
  bests <- do.call(rbind, lapply(
    split(rr0, list(rr0$scenario_id, rr0$rep, rr0$method), drop = TRUE),
    function(d) {
      ok <- !is.na(d$obj)
      if (!any(ok)) return(NULL)
      b <- d[ok, ][which.min(d$obj[ok]), , drop = FALSE]
      b$start <- sprintf("best-of-%d", nrow(d))
      b
    }))
  n_r <- table(interaction(rr0$scenario_id, rr0$rep, rr0$method, drop = TRUE))
  if (length(unique(n_r)) > 1L)
    cat(sprintf("WARNING: restarts per cell differ (%s). A best-of-N selected\n",
                paste(sort(unique(as.integer(n_r))), collapse = ",")),
        "from fewer restarts is not comparable with one selected from more.\n",
        "Check for missing task files before reporting it.\n", sep = "")
  all <- rbind(all, bests)
  cat(sprintf("computed best-of-N for %d cells (%d restarts each)\n",
              nrow(bests), as.integer(n_r)[1]))
}

cells <- unique(all[, c("scenario_id", "rep")])
cat(sprintf("cells: %d | scenarios: %s | methods: %s | restarts kept: %d\n",
            nrow(cells), paste(sort(unique(all$scenario_id)), collapse = ","),
            paste(sort(unique(all$method)), collapse = ","),
            sum(grepl("^random[0-9]", all$start)) / max(1L, nrow(cells))))
if (any(!is.na(all$error)))
  cat(sprintf("errors: %d rows\n", sum(!is.na(all$error))))

# --- the named starts, side by side ------------------------------------------
cat("\n=== aRI to truth by starting strategy ===\n")
main <- all[!grepl("^random[0-9]", all$start), ]
agg <- do.call(rbind, lapply(
  split(main, list(main$scenario_id, main$method, main$start), drop = TRUE),
  function(d) data.frame(
    scenario_id = d$scenario_id[1], method = d$method[1], start = d$start[1],
    # N_OK, NOT nrow(d). A start that ABORTS still contributes a row, with aRI
    # NA, so nrow() counts attempts and the median is taken over survivors only.
    # Reporting them as one number said "A01 ward, median 0.891, n = 10" when
    # nine of the ten aborted and the median was over a single replicate. Same
    # survivorship trap as the Study A tables; do not collapse these again.
    n_attempt = nrow(d), n_ok = sum(!is.na(d$aRI)),
    aRI = stats::median(d$aRI, na.rm = TRUE),
    conv = mean(d$converged, na.rm = TRUE), stringsAsFactors = FALSE)))
rownames(agg) <- NULL
for (s in sort(unique(agg$scenario_id))) {
  cat(sprintf("\n  %s\n", s))
  cat(sprintf("    %-14s %10s %8s %8s %8s\n",
              "start", "median|ok", "conv", "n_ok", "attempts"))
  d <- agg[agg$scenario_id == s, ]
  for (i in order(-d$aRI, na.last = TRUE))
    cat(sprintf("    %-14s %10s %8.2f %8d %8d%s\n", d$start[i],
                if (is.na(d$aRI[i])) "---" else sprintf("%.3f", d$aRI[i]),
                d$conv[i], d$n_ok[i], d$n_attempt[i],
                if (d$n_ok[i] < d$n_attempt[i])
                  sprintf("   <- %d ABORTED", d$n_attempt[i] - d$n_ok[i]) else ""))
}

# --- THE QUANTITY THAT SAYS HOW MUCH THE START MATTERS -----------------------
#
# The spread ACROSS random restarts within a replicate. A wide spread means the
# loop is finding different local optima on identical data; a narrow one means
# the start does not matter. `best-of-N` is what a practitioner would actually
# report, and it is kept separate so the spread is not hidden behind it.
cat("\n\n=== spread across the random restarts, within each replicate ===\n")
rr <- all[grepl("^random[0-9]", all$start), ]
if (nrow(rr)) {
  s <- do.call(rbind, lapply(
    split(rr, list(rr$scenario_id, rr$rep, rr$method), drop = TRUE),
    function(d) data.frame(
      scenario_id = d$scenario_id[1], rep = d$rep[1], method = d$method[1],
      n_ok = sum(!is.na(d$aRI)),
      min = min(d$aRI, na.rm = TRUE), med = stats::median(d$aRI, na.rm = TRUE),
      max = max(d$aRI, na.rm = TRUE), conv = mean(d$converged, na.rm = TRUE),
      stringsAsFactors = FALSE)))
  rownames(s) <- NULL
  cat(sprintf("  %-6s %-12s %5s %8s %8s %8s %8s\n",
              "scen", "method", "n_ok", "min", "median", "max", "conv"))
  for (i in order(s$scenario_id, s$method, s$rep))
    cat(sprintf("  %-6s %-12s %5d %8.3f %8.3f %8.3f %8.2f\n", s$scenario_id[i],
                substr(s$method[i], 1, 12), s$n_ok[i], s$min[i], s$med[i],
                s$max[i], s$conv[i]))
  cat("\n  by scenario, pooled over replicates:\n")
  for (sc in sort(unique(s$scenario_id))) {
    d <- s[s$scenario_id == sc, ]
    cat(sprintf("    %-6s median spread (max - min) = %.3f, worst restart %.3f, best %.3f\n",
                sc, stats::median(d$max - d$min), min(d$min), max(d$max)))
  }
}

# --- the C1 / C4b question, stated so it cannot be read past ------------------
cat("\n\n=== IS A CHANCE-LEVEL RESULT ALWAYS A NON-CONVERGED RUN? ===\n")
cat("C1 at nK_fit = 2 found that it was, and that the only converged run was the\n")
cat("good one. C4b at nK_fit = 4 found non-convergence to be COSMETIC, with\n")
cat("partitions agreeing with the strict fit at 0.968 to 1.000. Both cannot\n")
cat("describe the current operating point. This is the table that settles it.\n\n")
ok <- !is.na(all$aRI) & !is.na(all$converged)
if (any(ok)) {
  tb <- table(converged = all$converged[ok], chance = all$aRI[ok] < 0.1)
  print(tb)
  cat("\n")
  cv <- all$converged[ok]; ch <- all$aRI[ok] < 0.1
  cat(sprintf("  converged AND chance-level: %d fits", sum(cv & ch)))
  cat(if (sum(cv & ch) > 0)
        "  -> convergence does NOT guarantee a good partition at M = 4\n"
      else "  -> none; convergence still implies a good partition\n")
  cat(sprintf("  NOT converged and GOOD (aRI > 0.5): %d fits", sum(!cv & all$aRI[ok] > 0.5)))
  cat(if (sum(!cv & all$aRI[ok] > 0.5) > 0)
        "  -> non-convergence is COSMETIC here, as C4b found\n"
      else "  -> non-convergence still means failure, as C1 found at M = 2\n")
}

summ <- agg
utils::write.csv(all,  paste0(stem, "_raw.csv"),     row.names = FALSE)
utils::write.csv(summ, paste0(stem, "_summary.csv"), row.names = FALSE)
cat(sprintf("\nwrote %s_raw.csv (%d rows) and %s_summary.csv (%d rows)\n",
            stem, nrow(all), stem, nrow(summ)))
