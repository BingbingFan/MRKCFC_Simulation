# The chunked, best-iterate K-centres loop, as one function.
#
# Sourced, not run.
#
# WHY THIS EXISTS. `realdata_fit_task.R` grew this loop inline on 2026-09-09 and
# `realdata_subsample_task.R` needs exactly the same thing. Two copies of a loop
# whose subtleties are (a) `maxIter = m + 1` to obtain m accepted reassignments,
# (b) resuming from the FINAL iterate while reporting the BEST one, and (c) the
# difference between an init-only stop and a degenerate one, is precisely how
# this project has produced plausible wrong numbers before.
#
# NOT ADOPTED BY realdata_fit_task.R, which was running when this was written
# and could not change mid-array; the loop below is a faithful extraction of
# it, used by realdata_subsample_task.R.
#
# THE LOOP DOES NOT CONVERGE ON THIS COHORT. See the note in
# realdata_fit_task.R's write_ckpt(): the objective fluctuates without trend, so
# the LAST iterate is an arbitrary point of a trajectory and the BEST one is
# what any downstream comparison must use.

# `objective` takes labels and returns a number to MINIMISE. `on_chunk`, if
# given, is called as on_chunk(status, labels, best_labels, best_obj, iter) after
# each chunk, for callers that checkpoint to disk.
fit_best_iterate <- function(Ly, Lt, K, robust, scatter, optns,
                             objective, init.labels = NULL,
                             maxIter = 100L, conv.prop = 0, LOO = FALSE,
                             ncore = 1L, chunk = 5L, trace = 1,
                             on_chunk = NULL, label = "") {
  labels <- init.labels
  total_iter <- 0L; chunks <- 0L; secs <- 0
  history <- list()
  best_obj <- Inf; best_labels <- NULL; best_iter <- NA_integer_
  trace_iter <- integer(0); trace_obj <- numeric(0)
  status <- "incomplete"

  record <- function(lab) {
    if (is.null(lab)) return(invisible(NULL))
    o <- tryCatch(objective(lab), error = function(e) NA_real_)
    if (is.finite(o)) {
      trace_iter <<- c(trace_iter, total_iter)
      trace_obj  <<- c(trace_obj, o)
      if (o < best_obj) {
        best_obj <<- o; best_labels <<- lab; best_iter <<- total_iter
      }
    }
    invisible(o)
  }

  repeat {
    remaining <- maxIter - total_iter
    if (remaining < 1L) { status <- "maxiter"; break }
    m <- min(chunk, remaining)

    t0 <- Sys.time()
    fit <- tryCatch(
      MRKCFC::MKCFC(Ly = Ly, Lt = Lt, K = K, robust = robust, optns = optns,
                   mfpca_scatter = scatter,
                   init.method = "random", init.labels = labels,
                   maxIter = m + 1L,     # m+1 yields m ACCEPTED reassignments
                   parallel = ncore > 1L,
                   nmcores = if (ncore > 1L) ncore else NULL,
                   trace = trace, conv.prop = conv.prop, LOO = LOO),
      error = function(e) e)
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    secs <- secs + el; chunks <- chunks + 1L

    if (inherits(fit, "error")) {
      # error-early, not init-only: an error inside a chunk loses whatever the
      # chunk had accepted, so nothing can be said about the start itself.
      status <- if (total_iter == 0L) "error-early" else "error"
      attr(status, "message") <- conditionMessage(fit)
      if (trace >= 1)
        cat(sprintf("%sERROR in chunk %d after %.1f min: %s\n", label, chunks,
                    el / 60, conditionMessage(fit)))
      break
    }

    gained <- as.integer(fit$iter) - 1L
    new    <- as.integer(fit$res_clust$g)
    if (trace >= 1)
      cat(sprintf("%schunk %d: +%d iter (total %d) %.1f min = %.1f min/iter | sizes %s\n",
                  label, chunks, gained, total_iter + gained, el / 60,
                  if (gained > 0) el / 60 / gained else NA_real_,
                  paste(table(new), collapse = "/")))

    if (gained == 0L) {
      # Nothing accepted. With nothing banked, `new` IS the starting partition.
      status <- if (total_iter == 0L) "init-only" else "degenerate"
      break
    }

    total_iter <- total_iter + gained
    labels     <- new
    record(labels)
    repeated <- any(vapply(history, function(h) identical(h, labels), logical(1)))
    history[[length(history) + 1L]] <- labels

    if (!is.null(on_chunk))
      on_chunk("incomplete", labels, best_labels, best_obj, total_iter)

    if (isTRUE(fit$conv) || repeated) { status <- "converged"; break }
    if (gained < m)                   { status <- "degenerate"; break }
  }

  list(labels = best_labels,
       objective = if (is.finite(best_obj)) best_obj else NA_real_,
       best_iter = best_iter, final_labels = labels,
       iter = total_iter, chunks = chunks, seconds = secs,
       status = as.character(status),
       error = if (!is.null(attr(status, "message"))) attr(status, "message")
               else NA_character_,
       obj_trace = data.frame(iter = trace_iter, objective = trace_obj))
}
