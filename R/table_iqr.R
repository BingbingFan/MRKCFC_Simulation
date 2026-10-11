# Median (IQR) cells for the generated tables.
#
# The combined summaries carry medians; the interquartile range comes either
# from the summary itself (combine.R writes "median (IQR)" strings, which
# numify() in make_tables.R splits into <col> and <col>_iqr) or from the
# per-replicate file of the study. When it comes from the per-replicate file,
# check_median() stops the run if the median recomputed from those values
# differs from the summary's, so that the IQR describes exactly the replicates
# the printed median is taken over. Sourced by make_tables.R.

IQR_SENTENCE <- "Each cell is the median over replicates, with the interquartile range in parentheses."

iqr_of <- function(x) if (sum(!is.na(x)) < 2L) NA_real_ else stats::IQR(x, na.rm = TRUE)

check_median <- function(x, med, what, tol = 1e-6) {
  m <- stats::median(x, na.rm = TRUE)
  if (is.na(med) && is.na(m)) return(invisible(TRUE))
  if (is.na(med) || is.na(m) || abs(m - med) > tol)
    stop(sprintf("%s: the median recomputed from the per-replicate file (%s) differs from the summary (%s)",
                 what, format(m), format(med)))
  invisible(TRUE)
}

# fmt() is defined in make_tables.R; this pairs a formatted value with its IQR.
fmtq <- function(v, q, d = 3) paste0(fmt(v, d), ifelse(is.na(v) | is.na(q), "", paste0(" (", fmt(q, d), ")")))
