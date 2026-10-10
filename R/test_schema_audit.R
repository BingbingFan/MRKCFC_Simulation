# Self-test for schema_audit.R.
#
#   Rscript R/test_schema_audit.R
#
# A CHECKER THAT CANNOT FAIL IS WORSE THAN NONE: it converts an unexamined tree
# into a green light. This reintroduces each known defect into a COPY of R/ and
# asserts the audit flags it, and asserts the real tree passes.
#
# The unparseable case is here because the first version of schema_audit.R got
# it WRONG -- it printed a note, returned an empty column set and still exited 0,
# so a driver with a syntax error sailed through the gate. That was found by
# running this test, not by reading the code.

SRC <- if (dir.exists("R")) "R" else "."          # run from study/ or from R/
AUD <- file.path(SRC, "schema_audit.R")
stopifnot(file.exists(AUD))
RSCRIPT <- file.path(R.home("bin"), "Rscript")
TMP <- file.path(tempdir(), "audit_selftest")

run_status <- function(dir) {
  o <- suppressWarnings(system2(RSCRIPT, c(shQuote(normalizePath(AUD)), shQuote(dir)),
                                stdout = TRUE, stderr = TRUE))
  s <- attr(o, "status")
  list(status = if (is.null(s)) 0L else as.integer(s), out = o)
}

fresh <- function() {
  unlink(TMP, recursive = TRUE); dir.create(TMP, recursive = TRUE)
  file.copy(list.files(SRC, pattern = "[.]R$", full.names = TRUE), TMP)
}
# Patch through R, not the shell: a PowerShell regex replace mangled the copy's
# encoding on the first attempt and the "defect" it produced was a parse error,
# not the schema drift under test.
patch <- function(file, from, to) {
  p <- file.path(TMP, file)
  s <- paste(readLines(p, warn = FALSE), collapse = "\n")
  if (!grepl(from, s, fixed = TRUE))
    stop("anchor not found in ", file, " -- the source moved, update this test")
  writeLines(strsplit(sub(from, to, s, fixed = TRUE), "\n")[[1]], p)
}

cases <- list(
  list(name = "clean tree (must PASS)", want = 0L, patch = function() invisible(NULL)),
  # the real 2026-09-04 defect in Study A's driver; cross-file, pass two
  list(name = "run_cell.R missing iter (cross-file)", want = 1L, patch = function()
    patch("run_cell.R", "converged = NA, iter = NA_integer_, seconds = NA_real_,",
                        "converged = NA, seconds = NA_real_,")),
  # the real 2026-09-04 defect in C3; within-file, pass one
  list(name = "study_c3.R missing error (within-file)", want = 1L, patch = function()
    patch("study_c3.R", "                 error = NA_character_,\n", "")),
  list(name = "unparseable file", want = 1L, patch = function()
    cat("f <- function(  {\n", file = file.path(TMP, "run_cell.R"), append = TRUE))
)

ok <- TRUE
for (cs in cases) {
  fresh(); cs$patch()
  r <- run_status(TMP)
  pass <- identical(r$status, cs$want)
  ok <- ok && pass
  cat(sprintf("%-42s want exit %d, got %d  %s\n", cs$name, cs$want, r$status,
              if (pass) "PASS" else "*** FAIL ***"))
  if (!pass) cat(paste0("    ", r$out, collapse = "\n"), "\n")
}
unlink(TMP, recursive = TRUE)
cat(if (ok) "\nself-test passed\n" else "\nSELF-TEST FAILED\n")
quit(save = "no", status = if (ok) 0L else 1L)
