# Static check that every result frame agrees with its own failure frame.
#
#   Rscript R/schema_audit.R [dir-of-R-files]      # default: R
#
# Exits 1 if anything is found, so it can gate a submission:
#
#   Rscript R/schema_audit.R && <submit the array>
#
# WHY THIS FILE IS IN THE REPO
#
# Every driver pairs a success frame with a tryCatch failure frame, and the two
# must carry the same columns. They are written hundreds of lines apart, so when
# a column is added to one it is easy to miss the other -- which is exactly what
# happened when `angle_lead4` was added to estimation.R on 2026-09-04. Auditing
# for that one defect turned up two more that had been latent for weeks:
#
#   run_cell.R   STUDY A's driver: failure frame was missing `iter` (15 columns
#                against the success path's 16)
#   study_c3.R   success frame was missing `error` (6 against 7)
#
# THE FAILURE MODE IS NOT GRACEFUL. Drivers end with do.call(rbind, out), which
# ABORTS on mismatched columns. One failed cell would therefore take down the
# whole task, discard every cell already completed in it, and surface an rbind
# error instead of the real cause -- precisely what the `error` column exists to
# prevent. Neither defect was ever observed to fire, which is the point: this is
# the class of bug that costs a whole run the first time a cell throws.
#
# RUN IT WHENEVER A COLUMN IS ADDED, before submitting.

args <- commandArgs(trailingOnly = TRUE)
D <- if (length(args)) args[1] else "R"

# Collect the named arguments of every data.frame() call in an expression tree.
dfcalls <- function(e, acc = list()) {
  if (is.call(e)) {
    fn <- e[[1]]
    if ((is.name(fn) && as.character(fn) == "data.frame") ||
        (is.call(fn) && identical(as.character(fn[[1]]), "::") &&
         as.character(fn[[3]]) == "data.frame")) {
      nm <- names(as.list(e))[-1]
      nm <- setdiff(nm[nzchar(nm)], c("stringsAsFactors", "check.names", "row.names"))
      if (length(nm)) acc[[length(acc) + 1L]] <- nm
    }
    # Empty args (the blank in x[i, ]) are the missing symbol; BINDING one to a
    # variable and then using it errors, so index inline and skip it first.
    l <- as.list(e)
    for (i in seq_along(l)) {
      if (identical(l[[i]], quote(expr = ))) next
      if (is.call(l[[i]])) acc <- dfcalls(l[[i]], acc)
    }
  }
  acc
}

files <- sort(list.files(D, pattern = "[.]R$", full.names = TRUE))
found <- 0L

# A FILE THAT DOES NOT PARSE MUST COUNT AS A FINDING. The first version returned
# an empty set and printed a note, so an unparseable driver was silently skipped
# and the gate still exited 0 -- green-lighting a submission with broken R, which
# is worse than the schema drift this exists to catch. Caught by self-test.
unparseable <- character(0)
sets_in <- function(f) {
  ex <- tryCatch(parse(f), error = function(e) e)
  if (inherits(ex, "error")) {
    unparseable <<- union(unparseable, f)
    return(list())
  }
  unique(do.call(c, lapply(ex, function(e) dfcalls(e))))
}

# --- pass one: within a file -------------------------------------------------
# Two frames sharing most of their names are a success/failure pair, not two
# unrelated frames that happen to live in the same file.
cat("== pass one: success and failure frames within each file ==\n")
for (f in files) {
  s <- sets_in(f)
  if (length(s) < 2) next
  for (i in seq_len(length(s) - 1L)) for (j in (i + 1L):length(s)) {
    a <- s[[i]]; b <- s[[j]]
    sh <- length(intersect(a, b)); un <- length(union(a, b))
    if (sh >= 3 && sh < un && sh / un >= 0.5) {
      found <- found + 1L
      cat("\n### ", basename(f), "\n", sep = "")
      cat("  only in one: ", paste(setdiff(a, b), collapse = " "), "\n", sep = "")
      cat("  only in other: ", paste(setdiff(b, a), collapse = " "), "\n", sep = "")
    }
  }
}
if (!found) cat("  clean\n")

# --- pass two: across files --------------------------------------------------
# Pass one CANNOT see the run_cell.R defect: its success rows come from
# run_all_methods() in methods.R, so there is no second frame in the file to
# compare against. Any driver calling it without subsetting must cover all of
# its columns. study_c4.R and study_c5.R DO subset, to the same seven on both
# paths, and are correct.
cat("\n== pass two: drivers against run_all_methods() ==\n")
canon <- local({
  s <- sets_in(file.path(D, "methods.R"))
  if (!length(s)) stop("cannot read methods.R to establish the canonical columns")
  s[[which.max(vapply(s, length, 1L))]]
})
cat("  canonical (", length(canon), "): ", paste(canon, collapse = " "), "\n", sep = "")

for (f in files) {
  if (basename(f) == "methods.R") next
  src <- readLines(f, warn = FALSE)
  # Strip comments first. estimation.R only NAMES run_all_methods() in a comment
  # -- it fits via MFPCA_k and has an unrelated schema -- and matching that gave
  # a false positive, which trains the reader to ignore this report.
  code <- sub("#.*$", "", src)
  if (!any(grepl("run_all_methods\\(", code))) next
  if (any(grepl('\\[, c\\("[^)]*"\\)\\]', code))) next   # subsets to its own set
  for (s in sets_in(f)) {
    if (!("error" %in% s)) next                          # not a result frame
    miss <- setdiff(canon, s)
    if (length(miss)) {
      found <- found + 1L
      cat("\n### ", basename(f), " result frame is SHORT\n", sep = "")
      cat("  missing: ", paste(miss, collapse = " "), "\n", sep = "")
    }
  }
}

if (length(unparseable)) {
  found <- found + length(unparseable)
  cat("\n### FILES THAT DO NOT PARSE\n")
  for (f in unparseable) cat("  ", basename(f), "\n", sep = "")
}

cat("\n", if (found) sprintf("%d finding(s) -- do not submit until resolved\n", found)
         else "no findings\n", sep = "")
quit(save = "no", status = if (found) 1L else 0L)
