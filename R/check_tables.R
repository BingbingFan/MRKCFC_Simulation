# Structural check on generated LaTeX: does every body row carry as many columns
# as the tabular preamble declares?
#
#   Rscript R/check_tables.R <tex-dir>
#
# A column mismatch surfaces only when the document is compiled, long after the
# numbers have been read and trusted. It costs a second to catch here.
d <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else file.path("..","doc","tables")
bad <- 0L
for (f in list.files(d, pattern = "[.]tex$", full.names = TRUE)) {
  x  <- readLines(f, warn = FALSE)
  bt <- grep("begin{tabular}", x, fixed = TRUE, value = TRUE)[1]
  al <- strsplit(bt, "{", fixed = TRUE)[[1]]
  al <- strsplit(al[length(al)], "}", fixed = TRUE)[[1]][1]
  ncol <- nchar(gsub("[^lrc]", "", al))
  rows <- x[endsWith(x, "\\\\")]
  rows <- rows[!grepl("multicolumn", rows, fixed = TRUE)]
  amp  <- vapply(rows, function(r)
    lengths(regmatches(r, gregexpr("&", r, fixed = TRUE))) + 1L, integer(1))
  ok <- length(amp) > 0 && all(amp == ncol)
  if (!ok) bad <- bad + 1L
  cat(sprintf("%-30s cols=%-2d rows=%-3d widths=%-8s %s\n", basename(f), ncol,
      length(rows), paste(sort(unique(amp)), collapse = ","),
      if (ok) "OK" else "** MISMATCH **"))
}
cat(sprintf("\n%d file(s) with mismatched columns\n", bad))
if (bad) quit(save = "no", status = 1)
