# Abbreviation footnote for the generated tables.
#
# abbr_footnote(lines) scans the caption, header, body and note of a table for
# the abbreviations listed below and returns one \footnotetext line in the form
# the typed tables use, "BMI indicates body mass index; CA19-9, carbohydrate
# antigen 19-9; ...", alphabetical by abbreviation; NULL when none is present.
# Expansions that themselves contain an abbreviation (FPCA) pull that entry in.
# Sourced by make_tables.R and make_rd_tables.R; run from the repository root.

ABBR <- list(
  # pattern (perl regex on the table text), printed label, expansion
  list("\\b[aA]RI\\b",   "aRI",        "adjusted Rand index"),
  list("\\bBMI\\b",      "BMI",        "body mass index"),
  list("CA19-9",         "CA19-9",     "carbohydrate antigen 19-9"),
  list("CA242",          "CA242",      "carbohydrate antigen 242"),
  list("\\bCCR\\b",      "CCR",        "correct classification rate"),
  list("\\bCEA\\b",      "CEA",        "carcinoembryonic antigen"),
  list("\\bCH\\b",       "CH",         "Calinski--Harabasz index"),
  list("\\bCI\\b",       "CI",         "confidence interval"),
  list("\\bFPCA\\b",     "FPCA",       "functional principal component analysis"),
  list("\\bFVE\\b",      "FVE",        "fraction of variance explained"),
  list("\\bHR\\b",       "HR",         "hazard ratio"),
  list("\\bIQR\\b",      "IQR",        "interquartile range"),
  list("\\bKFPCA\\b",    "KFPCA",      "Kendall's $\\tau$ functional principal component analysis"),
  list("\\blcmm\\b",     "lcmm",       "multivariate latent class mixed model"),
  list("\\bmclust\\b",   "mclust",     "Gaussian mixture on the multivariate FPCA scores"),
  list("\\bMKCFC\\b",    "MKCFC",      "multivariate K-centers functional clustering"),
  list("MRKCFC-M",       "MRKCFC-M",   "multivariate robust K-centers functional clustering with the second-moment scatter"),
  list("MRKCFC-S",       "MRKCFC-S",   "multivariate robust K-centers functional clustering with the spatial sign scatter"),
  list("R\\$k\\$-means", "R$k$-means", "$k$-means on the robust multivariate FPCA scores"),
  list("\\bRRE\\b",      "RRE",        "relative reconstruction error"),
  list("\\bSD\\b",       "SD",         "standard deviation"),
  list("\\bteigen\\b",   "teigen",     "mixture of multivariate $t$ distributions on the multivariate FPCA scores"),
  list("\\bTNM\\b",      "TNM",        "tumor-node-metastasis")
)

abbr_sentence <- function(text) {
  text <- paste(text, collapse = " ")
  hit <- vapply(ABBR, function(a) grepl(a[[1]], text, perl = TRUE), logical(1))
  # an expansion may use FPCA; list it too
  more <- vapply(ABBR, function(a)
    any(vapply(ABBR[hit], function(h) grepl(a[[1]], h[[3]], perl = TRUE), logical(1))), logical(1))
  hit <- hit | more
  if (!any(hit)) return(NULL)
  lab <- vapply(ABBR[hit], `[[`, character(1), 2)
  exp <- vapply(ABBR[hit], `[[`, character(1), 3)
  o   <- order(tolower(lab))
  lab <- lab[o]; exp <- exp[o]
  parts <- c(paste(lab[1], "indicates", exp[1]),
             if (length(lab) > 1) paste(lab[-1], exp[-1], sep = ", "))
  paste0(paste(parts, collapse = "; "), ".")
}

abbr_footnote <- function(lines) {
  s <- abbr_sentence(lines)
  if (is.null(s)) NULL else sprintf("\\footnotetext{%s}", s)
}
