# Publication-ready LaTeX tables, generated from the combiners' summary files.
#
#   Rscript R/make_tables.R <results-root> [out-dir]
#
#   Rscript R/make_tables.R $ROOT ../doc/tables
#   Rscript R/make_tables.R "C:/.../MRKCFC/HPC" "C:/.../MRKCFC/doc/tables"
#
# Requires \usepackage{booktabs}, \usepackage{amsmath} and \usepackage{rotating}
# (for the sidewaystable wrapper) in the document. Tables are emitted as
# tabular* at the full text width; see env() below.
#
# WHY THIS FILE IS IN THE REPO
#
# An earlier set of tables was produced by a script that was not kept with the
# project, so they could not be regenerated. The chain is now closed end to
# end:
#
#   fit -> <study>_task_NNNN.rds -> combine_*.R -> *_summary.csv -> here -> .tex
#
# IT READS SUMMARIES, NOT RAW FILES. The combiners already compute the medians,
# rates and completeness counts. Recomputing them here would create a second
# definition that can drift from the one the combiner prints to screen -- the
# failure that gave this project two divergent copies of silhouette_score().
# Where a table needs something the summary does not carry, add it to the
# combiner rather than reaching past it into the raw file.
#
# INPUTS ARE FOUND BY BASENAME, RECURSIVELY, so one invocation covers a results
# tree with each study in its own directory. Anything missing is SKIPPED WITH A
# MESSAGE rather than failing, so this runs while studies are still going and
# can be re-run as each lands.
#
# ---------------------------------------------------------------------------
# TABLE NUMBERING, fixed 2026-09-04. THE MANUSCRIPT AND THE SUPPLEMENT CITE
# THESE NUMBERS. Do not renumber without updating both.
#
# Ordered by where each table appears in the argument: methods provenance, then
# simulation design, then the main comparison, then the evidence for
# Proposition 1, then complexity. The supplement follows the same sequence in
# more depth.
#
#   MAIN TEXT
#   Table 1   Simulation design: arms and constants     (manuscript text, not here)
#   Table 2   Study A, median aRI                       tab_studyA_ari
#   Table 3   Estimation study, principal angles        tab_estim_angle
#   Table 4   Study B, complexity, median aRI           tab_studyB_ari
#   Table 5   Real data, Cox models                     make_rd_tables.R
#
#   SUPPLEMENT
#   Table S1  Component selection by FVE                tab_fve_selection
#   Table S2  M sweep, aRI by M                         tab_sweep_M
#   Table S3  M sweep, convergence by M                 tab_sweep_convergence
#   Table S4  Study A, correct classification rate      tab_studyA_ccr
#   Table S5  Study A, reconstruction error             tab_studyA_rre
#   Table S6  Study A, convergence rate                 tab_studyA_convergence
#   Table S7  Estimation, score canonical correlation   tab_estim_cc
#   Table S8  Study B, convergence rate                 tab_studyB_convergence
#   Table S9  C1, initialisation sensitivity            tab_c1_init
#   Table S10 C2, class-specific vs class-common basis  tab_c2_basis
#   Table S11 C3, recovery of K                         tab_c3_kselect
#   Table S12 C4b, convergence thresholds               tab_c4b_agree
#
# Every numbered simulation table has a generator here; Table 1 is manuscript
# text and the real-data tables come from make_rd_tables.R.
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else "."
OUT  <- if (length(args) >= 2) args[2] else file.path("..", "doc", "tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("results root: %s\noutput:       %s\n\n", ROOT, OUT))

grab <- function(base) {
  f <- list.files(ROOT, pattern = paste0("^", base, "$"),
                  full.names = TRUE, recursive = TRUE)
  if (!length(f)) { cat(sprintf("  SKIP  %s\n", base)); return(NULL) }
  f <- f[order(file.info(f)$mtime, decreasing = TRUE)][1]   # newest wins
  cat(sprintf("  read  %s\n", f))
  utils::read.csv(f, stringsAsFactors = FALSE)
}

slab <- c(A01 = "Gaussian, $d=1.5$",           A02 = "Gaussian, $d=3$",
          A03 = "Fleishman$(0,4)$, $d=1.5$",   A04 = "Fleishman$(0,4)$, $d=3$",
          A05 = "Fleishman$(1.5,4)$, $d=1.5$", A06 = "Fleishman$(1.5,4)$, $d=3$",
          A07 = "Contaminated, $d=1.5$",       A08 = "Contaminated, $d=3$",
          A09 = "$t_{3}$, $d=1.5$",            A10 = "$t_{3}$, $d=3$",
          A11 = "$t_{1.5}$, $d=1.5$",          A12 = "$t_{1.5}$, $d=3$")
# DISPLAY names (2026-09-13): -M and -S name the second-stage scatter,
# moment and spatial sign. The codes M1/M2 and the summaries' own labels
# ("MRKCFC-R2 (new)", config/methods.csv) are unchanged: unlabel() parses
# those, this map only writes.
mlab <- c(M1 = "MRKCFC-M", M2 = "MRKCFC-S", M3 = "MKCFC", M4 = "mclust",
          M5 = "$k$-means", M6 = "R$k$-means", M7 = "lcmm", M8 = "teigen")
EOL <- " \\\\"                     # the LaTeX row terminator, written once
source(file.path("R", "table_abbr.R"))   # abbr_footnote(): the abbreviation note
source(file.path("R", "table_iqr.R"))    # iqr_of(), check_median(), fmtq(): median (IQR) cells

wr <- function(file, lines) {
  writeLines(lines, file.path(OUT, file)); cat(sprintf("  WROTE %s\n", file))
}
fmt <- function(x, d = 3) ifelse(is.na(x), "---", formatC(x, format = "f", digits = d))
# Bold the best in a row; `lower` flips the direction.
best <- function(v, lower = FALSE, d = 3, iqr = NULL) {
  s <- fmt(v, d)
  k <- if (all(is.na(v))) NA else if (lower) which.min(v) else which.max(v)
  if (!is.na(k)) s[k] <- paste0("\\textbf{", s[k], "}")
  if (!is.null(iqr)) s <- paste0(s, ifelse(is.na(v) | is.na(iqr), "", paste0(" (", fmt(iqr, d), ")")))
  s
}
# THE FLOAT WRAPPER IS SET HERE AND NOWHERE ELSE. The manuscript \input{}s these
# files unchanged, so anything about how a table sits on the page belongs in
# this function, not in a post-processing step that a regeneration would undo.
#
# tabular* at the full width with \extracolsep{\fill}: the columns spread to
# the text width and the caption runs the same width,
# instead of the table sitting narrow in the middle of the page. Tables too
# wide for the text column are rotated, decided by label so no call site has
# to know; a rotated table takes \textheight, the width of the rotated page.
WIDE    <- c("tab:studyA-ari", "tab:studyA-ccr", "tab:studyA-rre",
             "tab:estim-angle", "tab:c4b-agree", "tab:studyB-ari",
             "tab:studyB-convergence")
SMALLER <- c("tab:estim-cc")     # overruns the text width at \small
# CAPTION AND NOTE (2026-10-01). The caption is the first sentence
# only; everything after it goes below the table as \footnotetext, which is
# the sn-jnl convention and what make_rd_tables.R writes explicitly. The
# split is at the first ". " or "? " followed by a capital letter, a
# backslash or "$"; a one-sentence caption gets no note.
split_cap <- function(cap) {
  m <- regexpr("[.?] (?=[A-Z\\\\$])", cap, perl = TRUE)
  if (m < 0) return(list(cap = cap, note = NULL))
  list(cap = substr(cap, 1, m), note = trimws(substr(cap, m + 1, nchar(cap))))
}
env <- function(cap, lab, align, head, body) {
  sc    <- split_cap(cap)
  wide  <- lab %in% WIDE
  size  <- if (lab %in% SMALLER) "\\footnotesize" else "\\small"
  width <- if (wide) "\\textheight" else "\\textwidth"
  c(paste0(if (wide) "\\begin{sidewaystable}[htbp]" else "\\begin{table}[htbp]", size),
    "\\centering", sprintf("\\caption{%s}", sc$cap), sprintf("\\label{%s}", lab),
    sprintf("\\begin{tabular*}{%s}{@{\\extracolsep{\\fill}}%s@{}}", width, align),
    "\\toprule", head, "\\midrule", body, "\\bottomrule",
    "\\end{tabular*}",
    if (!is.null(sc$note)) sprintf("\\footnotetext{%s}", sc$note),
    abbr_footnote(c(sc$cap, head, body, sc$note)),
    if (wide) "\\end{sidewaystable}" else "\\end{table}")
}
# The completeness marker. The convention: never a superscript, but
# the count written out beside the value as (N_rep = k), so a reader does not
# have to look up what a small raised number means.
nrep <- function(k) sprintf(" ($N_{\\mathrm{rep}}=%d$)", as.integer(k))
# First matching value, or NA. Conditions are passed as logical vectors.
cellf <- function(d, col, ...) {
  v <- d[[col]][Reduce(`&`, list(...))]
  if (!length(v)) NA_real_ else v[1]
}

# ACCEPT EITHER CODES OR LABELS.
#
# combine.R relabels deliberately -- `scenario_id` becomes "Gaussian d=1.5" and
# `method` becomes "MRKCFC-R2 (new)" -- because results are reported by name, not
# by code. The C, sweep and FVE combiners keep the codes. This generator keys on
# codes internally, since A01..A12 and M1..M8 are stable while labels are prose,
# so a summary arriving relabelled is mapped BACK here.
#
# Without this the Study A tables silently came out with one column and one row:
# no error, just an empty table. Any schema mismatch must fail loudly or be
# handled; it must not produce a plausible-looking empty result.
CFG <- if (dir.exists("config")) "config" else file.path("study", "config")

# combine.R's summary is a DISPLAY table, not a machine-readable one: it writes
# aRI as the string "0.000 (0.001)", the median with its spread in parentheses,
# and conv_rate through sprintf. The leading token is the median, so take it.
#
# This is deliberately narrow and it FAILS LOUDLY. If a column that should be
# numeric parses to all-NA, that means combine.R's format changed and every
# table built from it would be blank -- exactly the silent-empty-result failure
# this file exists to avoid. Better to stop here than to emit a table of dashes.
numify <- function(d, cols) {
  for (cc in intersect(cols, names(d))) {
    if (is.numeric(d[[cc]])) next
    lead <- sub("^\\s*([-+0-9.eE]+).*$", "\\1", as.character(d[[cc]]))
    v <- suppressWarnings(as.numeric(lead))
    s <- as.character(d[[cc]])
    pq <- sub("^[^(]*\\(([-+0-9.eE]+)\\).*$", "\\1", s)
    pv <- suppressWarnings(as.numeric(pq)); pv[!grepl("(", s, fixed = TRUE)] <- NA
    if (any(!is.na(pv))) d[[paste0(cc, "_iqr")]] <- pv   # "median (IQR)" strings from combine.R
    if (all(is.na(v)) && any(nzchar(as.character(d[[cc]]))))
      stop(sprintf("cannot parse a number from column '%s' (e.g. '%s'). ",
                   cc, as.character(d[[cc]])[1]),
           "combine.R's summary format has changed; fix the parse or have the ",
           "combiner emit numerics.")
    d[[cc]] <- v
  }
  d
}

unlabel <- function(d) {
  if (!dir.exists(CFG)) return(d)
  if ("method" %in% names(d) && !"method_id" %in% names(d)) d$method_id <- d$method
  if ("method_id" %in% names(d) && !any(grepl("^M[0-9]+$", d$method_id))) {
    mt <- utils::read.csv(file.path(CFG, "methods.csv"), stringsAsFactors = FALSE)
    d$method_id <- mt$method_id[match(d$method_id, mt$label)]
  }
  if ("scenario_id" %in% names(d) && !any(grepl("^A[0-9]+$", d$scenario_id))) {
    st <- utils::read.csv(file.path(CFG, "scenarios.csv"), stringsAsFactors = FALSE)
    key <- sprintf("%s d=%g", st$arm_label, st$d)          # combine.R's format
    d$scenario_id <- st$scenario_id[match(d$scenario_id, key)]
  }
  d
}

# ------------------------------------------------------------------- Study A --
cat("Study A\n")
A <- grab("studyA_summary.csv")
if (!is.null(A)) {
  # combine.R names the id column `method` (run_all_methods writes the method_id
  # into it), while the C and sweep combiners write `method_id`. Accept either
  # rather than requiring one, so a schema difference does not silently produce
  # no tables -- which is exactly what it did on the first attempt.
  A <- numify(unlabel(A),
              c("aRI", "CCR", "NMI", "AMI", "MPE", "RRE", "RRE_trim",
                "conv_rate", "seconds", "iter"))
  # No M column in combine.R's summary. The fallback is the config's operating
  # point; if Study A is ever re-run at another M this caption will be wrong, so
  # add M to the combiner rather than editing the number here.
  MM  <- if ("M" %in% names(A)) A$M[1] else 4
  scA <- intersect(names(slab), unique(A$scenario_id))
  mA  <- intersect(names(mlab), unique(A$method_id))
  gA  <- function(s, m, col) cellf(A, col, A$scenario_id == s, A$method_id == m)
  gAi <- function(s, m, col) if (paste0(col, "_iqr") %in% names(A)) gA(s, m, paste0(col, "_iqr")) else NA_real_

  # SURVIVORSHIP MARKER. Every median here is taken over the replicates on which
  # that method returned a partition, so a cell built from 89 fits is not
  # comparable with one built from 100 -- and the failures are not random, they
  # are concentrated on the hard replicates. Unmarked in round one, which let
  # teigen hold a bolded win on t(3) from about 90 fits against MRKCFC-R2's 100.
  # The count itself is the marker; no legend to look up.
  REPS <- if ("n_fits" %in% names(A)) max(A$n_fits, na.rm = TRUE) else NA_real_
  nmark <- function(s, m) {
    if (is.na(REPS)) return("")
    k <- gA(s, m, "n_fits")
    if (is.na(k) || k >= REPS) "" else nrep(k)
  }
  NOTE_N <- paste(if ("aRI_iqr" %in% names(A)) IQR_SENTENCE else "", if (is.na(REPS)) "" else sprintf(paste(
    "Where fewer than %d of the %d fits succeeded, the cell gives the number",
    "of replicates the median is taken over as $(N_{\\mathrm{rep}}=\\cdot)$;",
    "unmarked cells use all %d. Those medians are conditional on the method",
    "converging and the excluded replicates are not a random subset, so a",
    "marked cell is not directly comparable with an unmarked one."),
    REPS, REPS, REPS))

  tabA <- function(col, cap, lab, lower = FALSE) env(paste(cap, NOTE_N), lab,
    paste0("l", strrep("r", length(mA))),
    paste0("Scenario & ", paste(mlab[mA], collapse = " & "), EOL),
    vapply(scA, function(s) paste0(slab[[s]], " & ", paste(paste0(
      best(vapply(mA, function(m) gA(s, m, col), numeric(1)), lower,
           iqr = vapply(mA, function(m) gAi(s, m, col), numeric(1))),
      vapply(mA, function(m) nmark(s, m), character(1))),
      collapse = " & "), EOL), character(1)))

  # CAPTIONS ARE SHORT TITLES (2026-10-10): the design constants sit in
  # the note, and "Study A/B/C" is not in a title. Notes are plain prose: no
  # bold phrases, no dashes, no questions.
  wr("tab_studyA_ari.tex", tabA("aRI", sprintf(paste(
    "Median adjusted Rand index by scenario. $n=500$, $M=%s$,",
    "at least four visits per subject, an iteration cap of 100 and a reclassification threshold of zero;",
    "best per row in bold."), MM), "tab:studyA-ari"))
  wr("tab_studyA_ccr.tex", tabA("CCR", sprintf(paste(
    "Median correct classification rate by scenario. $n=500$, $M=%s$;",
    "best per row in bold."), MM),
    "tab:studyA-ccr"))
  if ("RRE" %in% names(A))
    wr("tab_studyA_rre.tex", tabA("RRE", sprintf(paste(
      "Median relative reconstruction error by scenario. $n=500$, $M=%s$;",
      "lower is better, best per row in bold. RRE favors methods that minimize",
      "Equation~(15) and should be read within the MKCFC family rather than",
      "across families. Every partition is scored by one common reconstruction,",
      "so the measure describes the partition and not the fitting apparatus",
      "that produced it."), MM),
      "tab:studyA-rre", lower = TRUE))

  # Convergence is REPORTED, never used as a filter, and it belongs beside the
  # accuracy it qualifies. C4b established that at M = 4 a low rate reflects a
  # few boundary subjects reassigning indefinitely rather than an unsettled
  # partition: relaxing conv.prop doubles the rate while the partitions still
  # agree with the strict fit at 0.90 to 1.00 (Table S12 body; 2026-10-09).
  if ("conv_rate" %in% names(A)) {
    mF <- intersect(c("M1", "M2", "M3"), mA)
    wr("tab_studyA_convergence.tex", env(paste(
      "Convergence rate of the K-centers loop by scenario. $n=500$, $M=4$,",
      "a reclassification threshold of zero and an iteration cap of 100. A low rate reflects a",
      "few boundary subjects reassigning indefinitely rather than an unsettled",
      "partition: relaxing the threshold raises the rate by 0.04 to 0.52",
      "wherever it was below one, and the partitions still agree with the",
      "strict fit at 0.90 to 1.00 (Supplementary Table~\\ref{tab:c4b-agree}).",
      "Rates are reported beside the accuracy they qualify and are not used to",
      "filter results."),
      "tab:studyA-convergence", paste0("l", strrep("r", length(mF))),
      paste0("Scenario & ", paste(mlab[mF], collapse = " & "), EOL),
      vapply(scA, function(s) paste0(slab[[s]], " & ", paste(fmt(
        vapply(mF, function(m) gA(s, m, "conv_rate"), numeric(1)), 2),
        collapse = " & "), EOL), character(1))))
  }
}

# ------------------------------------------------------------------ Study B --
#
# STUDY B IS A STRESS TEST OF THE DESIGN, NOT A METHOD COMPARISON. Every cell is
# arm 4 (contaminated normal) at d = 1.5, chosen so that a complexity effect is
# attributable to the factor moved rather than to a distribution already
# breaking the estimator, and low enough not to hide inside a ceiling. It
# therefore says NOTHING about distributional robustness -- that is Study A's
# job, on arms Study B does not contain. Read it as "where does the method
# break", and the answer is proportions and dimension.
#
# `unlabel()` MUST NOT BE USED HERE. Its scenario branch maps against
# `scenarios.csv`'s A-codes, and Study B's `scenario_id` holds config LABELS
# ("d = 2", "P = 10"), so every row would silently become NA and the table would
# come out empty -- the same failure that produced a one-row Study A table on
# the first attempt. Study B joins on its own config file instead.
cat("\nStudy B\n")
B <- grab("studyB_summary.csv")
if (!is.null(B)) {
  sbf <- file.path(CFG, "scenarios_b.csv")
  if (!file.exists(sbf)) stop("Study B needs ", sbf)
  sb <- utils::read.csv(sbf, stringsAsFactors = FALSE)

  # methods: combine.R writes labels; map back to codes via methods.csv
  if ("method" %in% names(B) && !("method_id" %in% names(B))) {
    mt <- utils::read.csv(file.path(CFG, "methods.csv"), stringsAsFactors = FALSE)
    B$method_id <- mt$method_id[match(B$method, mt$label)]
  }
  # configs: scenario_id is the LABEL; join to config_id and axis
  B$config_id <- sb$config_id[match(B$scenario_id, sb$label)]
  B$axis      <- sb$axis[match(B$scenario_id, sb$label)]
  if (anyNA(B$config_id) || anyNA(B$method_id))
    stop("Study B join failed: ", sum(is.na(B$config_id)), " config and ",
         sum(is.na(B$method_id)), " method labels unmatched. ",
         "A silent NA here empties the table -- fix the join, not the caption.")

  B <- numify(B, c("aRI", "CCR", "conv_rate", "n_fits"))
  mB  <- intersect(names(mlab), unique(B$method_id))
  # ORDER BY AXIS, then config id. The config file's own order puts B13 last,
  # after `dimension`, even though it is a sparsity row -- which splits the
  # sparsity block in two and makes the table read as though the axis changed
  # back. Sorting here keeps each axis contiguous.
  albl <- c(baseline = "baseline", n = "sample size", sparsity = "sparsity",
            separation = "separation", proportions = "proportions",
            dimension = "dimension")
  ord <- sb[sb$config_id %in% B$config_id, c("config_id", "axis")]
  ord <- ord[order(match(ord$axis, names(albl)), ord$config_id), ]
  cB  <- ord$config_id
  # Print the axis once per group; a label repeated down every row is noise.
  axcol <- albl[ord$axis]
  axcol[duplicated(ord$axis)] <- ""
  gB  <- function(cc, m, col) cellf(B, col, B$config_id == cc, B$method_id == m)
  gBi <- function(cc, m, col) if (paste0(col, "_iqr") %in% names(B)) gB(cc, m, paste0(col, "_iqr")) else NA_real_
  REPSB <- max(B$n_fits, na.rm = TRUE)
  nmarkB <- function(cc, m) {
    k <- gB(cc, m, "n_fits")
    if (is.na(k) || k >= REPSB) "" else nrep(k)
  }
  rowlab <- function(cc) {
    i <- match(cc, sb$config_id)
    lab <- sb$label[i]
    # The sparsity rows also lower the visit floor from 4 to 2 (obs_min in
    # scenarios_b.csv); say so in the label, as B13 already does (2026-10-09).
    if (identical(sb$axis[i], "sparsity") && as.integer(sb$obs_min[i]) == 2L && !grepl("obs_min", lab))
      lab <- paste0(lab, ", obs_min 2")
    lab <- sub("obs_min 2", "at least 2 visits", lab, fixed = TRUE)
    gsub("([%&#_])", "\\\\\\1", lab)
  }

  # The baseline cell is Study A's A07 refitted on its own 50 replicates
  # (studyB_raw.csv, reps 1-50), so the medians need not match Table 2's: on
  # 2026-10-10 MRKCFC-M read 0.698 here against 0.712 there. Said in the note.
  wr("tab_studyB_ari.tex", env(paste(
    "Median adjusted Rand index under one-factor perturbations of the design.",
    "$n=500$, $M=4$, 50 replicates. Every cell is the contaminated-normal arm",
    "at $d=1.5$, held fixed so that an effect is attributable to the factor",
    "moved rather than to a distribution already breaking the estimator, and",
    "low enough not to hide inside a ceiling. The table therefore shows where",
    "the method breaks under design complexity; distributional robustness is",
    "the subject of Table~\\ref{tab:studyA-ari}, on arms this design does not",
    "contain. The baseline row is the contaminated $d=1.5$ scenario of",
    "Table~\\ref{tab:studyA-ari} refitted on 50 replicates, so its medians can",
    "differ from that table's (MRKCFC-M 0.698 against 0.712). Where fewer than",
    sprintf("%d", REPSB), "fits succeeded,",
    "the cell gives the number of replicates the median is taken over as",
    "$(N_{\\mathrm{rep}}=\\cdot)$; those medians are conditional on convergence",
    "and the excluded replicates are not a random subset. Best per row in bold.",
    if ("aRI_iqr" %in% names(B)) IQR_SENTENCE else ""),
    "tab:studyB-ari", paste0("ll", strrep("r", length(mB))),
    paste0("Axis & Configuration & ", paste(mlab[mB], collapse = " & "), EOL),
    vapply(seq_along(cB), function(i) { cc <- cB[i]; paste0(
      axcol[i], " & ", rowlab(cc), " & ",
      paste(paste0(best(vapply(mB, function(m) gB(cc, m, "aRI"), numeric(1)), FALSE,
                        iqr = vapply(mB, function(m) gBi(cc, m, "aRI"), numeric(1))),
                   vapply(mB, function(m) nmarkB(cc, m), character(1))),
            collapse = " & "), EOL) }, character(1))))

  # K-centers methods only: methods.R sets converged = TRUE for every score
  # method that returns, so teigen's column was 1.00 by construction (2026-10-10).
  mBc <- intersect(c("M1", "M2", "M3"), mB)
  if ("conv_rate" %in% names(B))
    wr("tab_studyB_convergence.tex", env(paste(
      "Convergence rate of the K-centers loop under the one-factor perturbations of Table~\\ref{tab:studyB-ari}.",
      "$n=500$, $M=4$, 50 replicates, a reclassification threshold of zero and",
      "an iteration cap of 100. Rates are reported beside the accuracy they",
      "qualify and are not used to filter results. teigen is omitted: it has no",
      "reclassification loop, and the replicates in which it returned no model",
      "appear as $N_{\\mathrm{rep}}$ in Table~\\ref{tab:studyB-ari}."),
      "tab:studyB-convergence", paste0("ll", strrep("r", length(mBc))),
      paste0("Axis & Configuration & ", paste(mlab[mBc], collapse = " & "), EOL),
      vapply(seq_along(cB), function(i) { cc <- cB[i]; paste0(
        axcol[i], " & ", rowlab(cc), " & ",
        paste(fmt(vapply(mBc, function(m) gB(cc, m, "conv_rate"), numeric(1)), 2),
              collapse = " & "), EOL) }, character(1))))
}

# ------------------------------------------------------------------ M sweep ---
cat("\nM sweep\n")
S <- grab("sweep_summary.csv"); S <- unlabel(S)
if (!is.null(S)) {
  # lcmm (M7) uses no principal components; its sweep row was one value at M = 1
  # and dashes, bolded, so it is left out of this table (2026-10-10).
  Ms <- sort(unique(S$M)); mS <- setdiff(intersect(names(mlab), unique(S$method_id)), "M7")
  pooled <- function(m, MM, col) mean(S[[col]][S$method_id == m & S$M == MM],
                                      na.rm = TRUE)
  # Pooled over scenarios: the M effect is a property of the DESIGN, and this
  # table is the evidence that most methods gain, not only the robust one
  # (k-means gains little; robust k-means falls to chance at M >= 4).
  wr("tab_sweep_M.tex", env(paste(
    "Median adjusted Rand index by retained component count.",
    "Averaged over the twelve scenarios, $n=500$, 25 replicates; best per row",
    "in bold. Five of the seven methods improve substantially beyond $M=2$, so",
    "the component count is a property of the design rather than of any one",
    "method; $k$-means gains little and R$k$-means falls to chance at $M\\geq4$.",
    "lcmm, which uses no principal components, is not in this table.",
    "The generating truth is $M_{\\mathrm{true}}=2$; $M>2$ is deliberate",
    "over-retention, since the second stage rotates $P\\times M$ stacked scores."),
    "tab:sweep-M", paste0("l", strrep("r", length(Ms))),
    paste0("Method & ", paste(sprintf("$M=%d$", Ms), collapse = " & "), EOL),
    vapply(mS, function(m) paste0(mlab[[m]], " & ", paste(best(
      vapply(Ms, function(MM) pooled(m, MM, "aRI"), numeric(1))),
      collapse = " & "), EOL), character(1))))

  # THIS MUST SIT NEXT TO THE TABLE ABOVE. Accuracy peaks at M = 5-6 while
  # convergence goes to zero there, so the two cannot be read apart.
  if ("conv_rate" %in% names(S)) {
    mF <- intersect(c("M1", "M2", "M3"), mS)
    wr("tab_sweep_convergence.tex", env(paste(
      "Convergence rate by retained component count.",
      "Averaged over the twelve scenarios, $n=500$, 25 replicates. For the two",
      "robust methods $M\\geq5$ is not reportable: the median iteration count",
      "sits on the cap of 100 iterations, so accuracy measured there is an",
      "artifact of where the loop was stopped. $M=4$ is the largest usable value",
      "for them, and is also what the variance-explained rule selects; MKCFC",
      "converges at every $M$."),
      "tab:sweep-convergence", paste0("l", strrep("r", length(Ms))),
      paste0("Method & ", paste(sprintf("$M=%d$", Ms), collapse = " & "), EOL),
      vapply(mF, function(m) paste0(mlab[[m]], " & ", paste(fmt(
        vapply(Ms, function(MM) pooled(m, MM, "conv_rate"), numeric(1)), 2),
        collapse = " & "), EOL), character(1))))
  }
}

# ---------------------------------------------------------------------- FVE ---
cat("\nFVE\n")
F5 <- grab("fve_summary.csv"); F1 <- grab("fve_summary_n1000.csv")
if (!is.null(F5)) {
  # Worst variable AND worst branch: one nK_fit must serve all P variables and
  # every method being compared, so the most demanding of them governs.
  need <- function(d, s, thr) {
    Ms <- sort(unique(d$M))
    max(vapply(unique(d$branch), function(br)
      max(vapply(sort(unique(d$variable)), function(j) {
        v <- vapply(Ms, function(M) cellf(d, "FVE", d$scenario_id == s,
                    d$branch == br, d$variable == j, d$M == M), numeric(1))
        w <- which(v >= thr)[1]
        if (is.na(w)) max(Ms) + 1L else as.integer(Ms[w])
      }, integer(1))), integer(1)))
  }
  ff  <- function(k) if (k > 6) "$>6$" else as.character(k)
  scF <- intersect(names(slab), unique(F5$scenario_id))
  two <- !is.null(F1)
  hdr <- if (two) paste0(
      "Scenario & \\multicolumn{2}{c}{$\\geq 0.90$} & ",
      "\\multicolumn{2}{c}{$\\geq 0.95$}", EOL,
      " \\cmidrule(lr){2-3}\\cmidrule(lr){4-5} ",
      "& $n=500$ & $n=1000$ & $n=500$ & $n=1000$", EOL) else paste0(
      "Scenario & $\\geq 0.90$ & $\\geq 0.95$ & $\\geq 0.99$", EOL)
  body <- vapply(scF, function(s) {
    v <- if (two) c(need(F5, s, .90), need(F1, s, .90),
                    need(F5, s, .95), need(F1, s, .95))
         else     c(need(F5, s, .90), need(F5, s, .95), need(F5, s, .99))
    paste0(slab[[s]], " & ", paste(vapply(v, ff, character(1)), collapse = " & "), EOL)
  }, character(1))
  wr("tab_fve_selection.tex", env(paste(
    "Number of components selected by the variance-explained rule.",
    "Smallest $M$ at which the cumulative share of the leading eigenvalues",
    "reaches the threshold, from one global fit over all subjects before any",
    "clustering, which a practitioner can do. The most demanding variable and",
    "the more demanding of the two first-stage estimators govern, since a",
    "single $M$ must serve every method compared. The $\\geq 0.90$ rule selects",
    "$M=4$ on every scenario at both sample sizes, while the $M$ maximizing the",
    "adjusted Rand index is 5 for MRKCFC-S and 6 for MKCFC (3 for MRKCFC-M,",
    "Supplementary Table~\\ref{tab:sweep-M}), so the label-free rule is",
    "conservative for two of the three relative to selecting on the truth. The",
    "$\\geq 0.95$ rule selects $M=5$ at $n=1000$ and on seven of the twelve",
    "scenarios at $n=500$, where the robust loops do not converge."),
    "tab:fve-selection", paste0("l", strrep("r", if (two) 4 else 3)), hdr, body))
}

# ------------------------------------------------------------ C4 thresholds ---
cat("\nC4 thresholds\n")
C4 <- grab("c4b_summary.csv"); C4 <- unlabel(C4)
if (!is.null(C4)) {
  CP  <- sort(unique(C4$conv_prop)); scC <- sort(unique(C4$scenario_id))
  mC  <- intersect(c("M1", "M2", "M3"), unique(C4$method_id))
  gC  <- function(s, m, cp, col) cellf(C4, col, C4$scenario_id == s,
                                       C4$method_id == m, C4$conv_prop == cp)
  C4R <- grab("c4b_raw.csv")
  if (!is.null(C4R)) C4R <- unlabel(C4R)
  gCi4 <- function(s, m, cp) {
    if (is.null(C4R)) return(NA_real_)
    x <- C4R$agree_strict[C4R$scenario_id == s & C4R$method_id == m & C4R$conv_prop == cp]
    check_median(x, gC(s, m, cp, "agree_strict"), paste("thresholds", s, m, cp)); iqr_of(x)
  }
  nA  <- length(CP) - 1L            # agreement columns: the relaxed ones only
  body <- unlist(lapply(mC, function(m) vapply(scC, function(s) paste0(
    if (s == scC[1]) mlab[[m]] else "", " & ", s, " & ",
    # Agreement is measured against the strict fit, so its own column would be 1
    # by construction and is omitted.
    paste(fmtq(vapply(CP[-1], function(cp) gC(s, m, cp, "agree_strict"), numeric(1)),
               vapply(CP[-1], function(cp) gCi4(s, m, cp), numeric(1)), 3), collapse = " & "), " & ",
    paste(fmt(vapply(CP, function(cp) gC(s, m, cp, "conv_rate"),
                     numeric(1)), 2), collapse = " & "), EOL), character(1))))
  wr("tab_c4b_agree.tex", env(paste(
    "Agreement between the strict and relaxed convergence thresholds.",
    "$n=500$, $M=4$, 25 replicates. Every fit is scored against the",
    "partition fitted with a reclassification threshold of zero from the same replicate; a value near 1",
    "means relaxing the threshold left the partition unchanged. Agreement is",
    "0.90 to 1.00 throughout, while the convergence rate rises by 0.04 to 0.52",
    "at the loosest setting wherever it was below one, and the median adjusted",
    "Rand index against the truth moves by at most 0.04, except on the",
    "infinite-variance arm, where the replicate distribution is bimodal and the",
    "medians of MRKCFC-M and MKCFC shift by up to 0.19 and 0.33; the loop was",
    "reassigning a handful of boundary subjects indefinitely rather than failing",
    "to settle. A threshold of zero remains the reported setting; the relaxed",
    "columns are a sensitivity analysis, not a redefinition of convergence.",
    if (!is.null(C4R)) "For the agreement, each cell is the median over replicates, with the interquartile range in parentheses." else ""),
    "tab:c4b-agree", paste0("ll", strrep("r", nA + length(CP))),
    paste0("Method & Scenario & \\multicolumn{", nA,
           "}{c}{Agreement with the zero-threshold fit} & \\multicolumn{",
           length(CP), "}{c}{Convergence rate}", EOL,
           sprintf(" \\cmidrule(lr){3-%d}\\cmidrule(lr){%d-%d} ",
                   2L + nA, 3L + nA, 2L + nA + length(CP)),
           "& & ", paste(sprintf("$%g$", CP[-1]), collapse = " & "), " & ",
           paste(sprintf("$%g$", CP), collapse = " & "), EOL),
    body))
}

# ------------------------------------------- C1, initialisation (Table S9) ---
#
# BUILT AROUND THE ABORT COUNTS, NOT THE MEDIANS. The finding is not that the
# starting strategy shifts accuracy a little -- the medians differ by 0.015 --
# it is that the DETERMINISTIC starts fail outright while the random ones never
# do. k-means is `MKCFC()`'s default and it aborts on every replicate of the
# infinite-variance arm with "Each cluster must contain at least 10 curves".
# A table of medians alone would hide that behind four near-identical numbers.
cat("\nC1 initialization\n")
C1 <- grab("c1_summary.csv")
if (!is.null(C1)) {
  C1 <- numify(C1, c("aRI", "conv", "n_ok", "n_attempt"))
  ord <- c("best-of-20", "oracle", "kmeans", "ward")
  slb <- c(`best-of-20` = "Best of 20 random restarts", oracle = "Oracle (true labels)",
           kmeans = "$k$-means \\emph{(the default)}", ward = "Ward hierarchical")
  stC <- intersect(ord, unique(C1$start))
  scC <- intersect(names(slab), unique(C1$scenario_id))
  gC1 <- function(s, st, col) cellf(C1, col, C1$scenario_id == s, C1$start == st)
  C1R <- grab("c1_raw.csv")
  gC1i <- function(s, st) {
    if (is.null(C1R)) return(NA_real_)
    x <- C1R$aRI[C1R$scenario_id == s & C1R$start == st]
    check_median(x, gC1(s, st, "aRI"), paste("initialization", s, st)); iqr_of(x)
  }
  cellS9 <- function(s, st) {
    v <- gC1(s, st, "aRI"); ok <- gC1(s, st, "n_ok"); at <- gC1(s, st, "n_attempt")
    # No dagger: the returned/attempted count beside the median already says
    # that the median is over fewer fits than were attempted.
    paste0(if (is.na(v)) "---" else fmtq(v, gC1i(s, st), 3), " & ",
           sprintf("%.0f/%.0f", ok, at))
  }
  wr("tab_c1_init.tex", env(paste(
    "Sensitivity to the starting partition.",
    "$n=500$, $M=4$, 10 replicates per scenario, MRKCFC-S only. For each",
    "strategy the table gives the median adjusted Rand index among the fits",
    "that returned a partition, and the number that did so out of the number",
    "attempted; the second column of each pair carries the result. The",
    "deterministic starts abort rather than perform worse: every one of the 33",
    "failures raises the same condition, a class fallen below ten curves, and",
    "$k$-means, the software default, fails on every replicate of the",
    "infinite-variance scenario. Across 400 attempts the random restarts never",
    "abort. Selecting among 20 restarts by the algorithm's own objective, which",
    "uses no knowledge of the truth, also beats starting from the true labels",
    "on both scenarios.", if (!is.null(C1R)) IQR_SENTENCE else ""),
    "tab:c1-init", paste0("l", strrep("rr", length(scC))),
    paste0("Starting strategy & ",
           paste(sprintf("\\multicolumn{2}{c}{%s}", slab[scC]), collapse = " & "), EOL,
           paste(sprintf("\\cmidrule(lr){%d-%d}", 2L + 2L * (seq_along(scC) - 1L),
                         1L + 2L * seq_along(scC)), collapse = ""), " ",
           "& ", paste(rep("median & returned", length(scC)), collapse = " & "), EOL),
    vapply(stC, function(st) paste0(
      slb[[st]], " & ",
      paste(vapply(scC, function(s) cellS9(s, st), character(1)), collapse = " & "),
      EOL), character(1))))
}

# ------------------------------------------------- C3, selecting K (S11) -----
#
# THE CONTRAST IS BETWEEN CRITERIA, not between methods. The silhouette is the
# manuscript's criterion of record; Calinski-Harabasz is reported beside it
# because a selected K that depends on the criterion is itself the finding.
# CH fails completely paired with Gaussian mixture modelling on the hard cell --
# recovery 0.00 against the silhouette's 0.98 ON THE SAME FITS.
cat("\nC3 selecting K\n")
C3 <- grab("c3_summary.csv")
if (!is.null(C3)) {
  C3 <- numify(C3, c("recovery", "n_rep", "n_design"))
  cfg3 <- intersect(c("C_S0", "C_HARD"), unique(C3$config_id))
  cl3  <- c(C_S0 = "Reference", C_HARD = "Hard cell")
  crs  <- c(sil_euclid = "Silhouette", CH_euclid = "CH")
  crs  <- crs[names(crs) %in% unique(C3$criterion)]
  # This table keys on the summary's own labels rather than on method codes, so
  # the paper's names are mapped here (the same names as `mlab`), and the rows
  # follow the paper's order rather than the alphabet.
  m3disp <- c(`MRKCFC-R2 (new)` = "MRKCFC-S", `MRKCFC (current)` = "MRKCFC-M", MKCFC = "MKCFC",
              `MFPCA + mclust GMM` = "mclust", `MFPCA + k-means` = "$k$-means",
              `MRFPCA + k-means` = "R$k$-means", `t-mixture on MFPCA scores` = "teigen")
  m3   <- c(intersect(names(m3disp), unique(C3$method)),
            sort(setdiff(unique(C3$method), names(m3disp))))
  m3lab <- function(m) if (m %in% names(m3disp)) m3disp[[m]] else gsub("([%&#_])", "\\\\\\1", m)
  g3 <- function(cc, m, cr, col)
    cellf(C3, col, C3$config_id == cc, C3$method == m, C3$criterion == cr)
  cell11 <- function(cc, m, cr) {
    r <- g3(cc, m, cr, "recovery"); n <- g3(cc, m, cr, "n_rep")
    nd <- g3(cc, m, cr, "n_design")
    paste0(if (is.na(r)) "---" else fmt(r, 2),
           if (!is.na(n) && !is.na(nd) && n < nd) nrep(n) else "")
  }
  wr("tab_c3_kselect.tex", env(paste(
    "Recovery of the true number of classes by two criteria.",
    "True $K=2$; each criterion is maximized over $K=2,\\dots,6$; $n=500$,",
    "$M=4$, 50 replicates. Entries are the proportion of replicates selecting",
    "$K=2$, on the Euclidean distance. The silhouette, which the manuscript",
    "uses, is the more reliable of the two criteria throughout, and the",
    "difference is not marginal: on the hard cell Calinski--Harabasz fails",
    "completely when paired with Gaussian mixture modeling, selecting three to",
    "six classes in every replicate, while the silhouette on the same fits",
    "recovers $K=2$ in forty-nine.",
    "That is a failure of the criterion and not of the clustering. The two",
    "distances agree on 93.1 per cent of silhouette fits and 78.9 per cent of",
    "Calinski--Harabasz fits, so the metric is not a material source of",
    "variation. Neither criterion is defined at $K=1$, so the procedure answers",
    "how many classes there are given at least two, and cannot report the",
    "absence of clustering. Where fewer than 50 replicates returned a partition,",
    "the cell gives the number contributing as $(N_{\\mathrm{rep}}=\\cdot)$."),
    "tab:c3-kselect", paste0("l", strrep("rr", length(cfg3))),
    paste0("Method & ",
           paste(sprintf("\\multicolumn{2}{c}{%s}", cl3[cfg3]), collapse = " & "), EOL,
           paste(sprintf("\\cmidrule(lr){%d-%d}", 2L + 2L * (seq_along(cfg3) - 1L),
                         1L + 2L * seq_along(cfg3)), collapse = ""), " ",
           "& ", paste(rep(paste(crs, collapse = " & "), length(cfg3)),
                       collapse = " & "), EOL),
    vapply(m3, function(m) paste0(
      m3lab(m), " & ",
      paste(unlist(lapply(cfg3, function(cc)
        vapply(names(crs), function(cr) cell11(cc, m, cr), character(1)))),
        collapse = " & "), EOL), character(1))))
}

# ------------------------------------------------- C2, the basis ablation ----
#
# TABLE S10. Reports a RATE and a CONDITIONAL MEDIAN, never a pooled median.
# C_HARD is cleanly bimodal -- 17 of 50 replicates near zero, 26 above 0.75, an
# almost empty middle -- so the pooled median lands in whichever mode holds more
# than half the mass. MRKCFC-R2's 0.774 against MKCFC's 0.106 looks like a 7x
# accuracy effect and is not one: it is a 58-percent success rate against 24.
# combine_c.R computes both; this only lays them out.
cat("\nC2 basis ablation\n")
C2 <- grab("c2_summary.csv")
if (!is.null(C2)) {
  C2 <- numify(C2, c("aRI_median", "success_rate", "aRI_med_success", "n_fits"))
  mC <- intersect(names(mlab), unique(C2$method_id))
  # Reference first, hard cell second. `unique()` follows the summary's own
  # order, which is alphabetical and puts C_HARD ahead of C_S0 -- the reverse of
  # how the result reads.
  clbl <- c(C_S0 = "Reference", C_HARD = "Hard cell ($t_{1.5}$, $d=1.5$)")
  cfgC <- intersect(names(clbl), unique(C2$config_id))
  gC <- function(cc, m, b, col)
    cellf(C2, col, C2$config_id == cc, C2$method_id == m, C2$basis == b)
  C2R <- grab("c2_raw.csv")
  gCi <- function(cc, m, b) {
    if (is.null(C2R)) return(NA_real_)
    x <- C2R$aRI[C2R$config_id == cc & C2R$method_id == m & C2R$basis == b & is.na(C2R$error)]
    x <- x[!is.na(x) & x > 0.5]
    check_median(x, gC(cc, m, b, "aRI_med_success"), paste("basis", cc, m, b)); iqr_of(x)
  }
  pct <- function(x) if (is.na(x)) "---" else sprintf("%.0f\\%%", 100 * x)
  cell <- function(cc, m, b) {
    r <- gC(cc, m, b, "success_rate"); v <- gC(cc, m, b, "aRI_med_success")
    paste0(pct(r), " & ", if (is.na(v)) "---" else fmtq(v, gCi(cc, m, b), 3))
  }
  wr("tab_c2_basis.tex", env(paste(
    "Class-specific against class-common eigenfunctions.",
    "$n=500$, $M=4$, 50 replicates. For each method and basis the table gives",
    "the proportion of replicates recovering the partition (adjusted Rand",
    "index above 0.5) and the median among those replicates. The pair is",
    "reported rather than a pooled median because the hard cell is bimodal:",
    "replicates either recover the partition or fail almost completely, with",
    "very little between, so a pooled median records which side of that split",
    "holds more than half the mass and moves discontinuously with it.",
    "A class-common basis is misspecified by construction here, since the",
    "design gives two of the three variables class-specific eigenfunctions;",
    "the direction is expected and the magnitude is the result.",
    if (!is.null(C2R)) "For the medians among the recovered replicates, the interquartile range over those replicates is in parentheses." else ""),
    "tab:c2-basis", paste0("ll", strrep("r", 4)),
    paste0("Configuration & Method & \\multicolumn{2}{c}{Class-specific} & ",
           "\\multicolumn{2}{c}{Class-common}", EOL,
           "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6} & & recovered & median & ",
           "recovered & median", EOL),
    unlist(lapply(cfgC, function(cc) vapply(seq_along(mC), function(i) {
      m <- mC[i]
      paste0(if (i == 1L) clbl[[cc]] else "", " & ", mlab[[m]], " & ",
             cell(cc, m, "specific"), " & ", cell(cc, m, "common"), EOL)
    }, character(1))))))
}

# --------------------------------------------------------------- estimation ---
cat("\nEstimation study\n")
E <- grab("estim_summary.csv")
if (!is.null(E)) {
  vlab <- c(M1 = "MRKCFC-M (moment)", M2 = "MRKCFC-S (sign)", M3 = "MKCFC (moment)")
  vE  <- intersect(names(vlab), unique(E$variant))
  scE <- intersect(names(slab), unique(E$scenario_id))
  gE  <- function(s, v, col) cellf(E, col, E$scenario_id == s, E$variant == v)
  ER  <- grab("estim_raw.csv")
  gEi <- function(s, v, col) {
    if (is.null(ER)) return(NA_real_)
    x <- ER[[col]][ER$scenario_id == s & ER$variant == v]
    check_median(x, gE(s, v, col), paste("estimation", s, v, col)); iqr_of(x)
  }
  # ELLIPTICITY GATES THE CLAIM. The SSCM's eigenvectors coincide with the
  # covariance's only for elliptically symmetric distributions, and arms 2 and 3
  # (Fleishman via Vale-Maurelli) are not elliptical. Reporting them anyway is
  # the point: a robustness claim that holds only under ellipticity should say so.
  tabE <- function(col, cap, lab, lower) env(cap, lab,
    paste0("ll", strrep("r", length(vE))),
    paste0("Elliptical? & Scenario & ", paste(vlab[vE], collapse = " & "), EOL),
    vapply(scE, function(s) paste0(
      if (cellf(E, "arm", E$scenario_id == s) %in% c(2, 3)) "no" else "yes",
      " & ", slab[[s]], " & ", paste(best(
        vapply(vE, function(v) gE(s, v, col), numeric(1)), lower, 2,
        iqr = vapply(vE, function(v) gEi(s, v, col), numeric(1))),
        collapse = " & "), EOL), character(1)))

  # THE TRUNCATION IS A PROPERTY OF THE MEASUREMENT, NOT OF THE FIT. Every fit
  # here uses nK_fit = 4. `keep` is held at 3 and 4 on BOTH sides so the target
  # does not move with the setting under test. Leading-4 is shown beside the
  # headline only to answer the obvious "you fitted 4, why score 3?" -- its
  # eigen-gap is 15 per cent against 48, so it is the worse-conditioned cut.
  ANG <- c(angle_lead3 = "Leading-3 (48\\% gap)", angle_lead4 = "Leading-4 (15\\% gap)")
  ANG <- ANG[names(ANG) %in% names(E)]
  nv  <- length(vE)
  wr("tab_estim_angle.tex", env(paste(
    "Largest principal angle between the estimated and true eigenfunction",
    "subspaces.",
    "Degrees, fitted on the true partition, $n=500$, $M=4$; lower is better,",
    "best per row in bold. Fitting under the true partition removes the",
    "initializer, the K-centers loop and convergence, leaving only the question",
    "of whether the spatial sign scatter estimates the second stage better than",
    "the moment scatter. The truncation is a property of the measurement, not",
    "of the fit, and is held fixed on both sides so that the target does not",
    "move. Leading-3 is the headline because it is the only well-conditioned",
    "cut in this design: the relative eigen-gaps of $\\boldsymbol{Z}$ are 6.4 per cent below",
    "the second component, 48 per cent below the third and 15 per cent below",
    "the fourth, and by Davis--Kahan the error in an estimated eigen-subspace",
    "scales as perturbation over eigen-gap.",
    # CONDITIONAL, so an older summary without the column cannot produce a
    # caption describing a block that is not in the table. A mislabelled table
    # is worse than a missing one.
    #
    # THIS SENTENCE ORIGINALLY CLAIMED THE ORDERING WAS TRUNCATION-STABLE. It is
    # not, and the measured result replaced the claim (2026-09-04).
    if ("angle_lead4" %in% names(ANG)) paste(
      "At leading-4 every angle degrades and the three-way ranking shifts in",
      "three of the twelve scenarios, all on elliptical arms, as the poorer",
      "conditioning predicts; the sign scatter still beats the moment scatter",
      "on the same first stage in all twelve scenarios at both truncations,",
      "which is the comparison this table is for. Leading-4 is reported for",
      "completeness.") else "",
    "MRKCFC-M and MKCFC use the moment scatter; MRKCFC-S uses the spatial sign",
    "scatter, and the two blocks repeat the same three columns.",
    "The guarantee of Proposition~\\ref{prop1} applies only under elliptical",
    "symmetry, which the two Fleishman arms do not satisfy; those rows are reported",
    "rather than omitted.", if (!is.null(ER)) IQR_SENTENCE else ""), "tab:estim-angle",
    paste0("ll", strrep("r", nv * length(ANG))),
    # SHORT sub-headers. Repeating "(moment)"/"(sign)" across both blocks ran the
    # eight-column table 3.07pt over the text width; the scatter each method uses
    # is stated in the caption instead, and these names now match the Study A
    # tables. Widen this header again and the table overflows.
    paste0("Elliptical? & Scenario & ",
           paste(sprintf("\\multicolumn{%d}{c}{%s}", nv, ANG), collapse = " & "), EOL,
           paste(sprintf("\\cmidrule(lr){%d-%d}", 3L + nv * (seq_along(ANG) - 1L),
                         2L + nv * seq_along(ANG)), collapse = ""), " ",
           "& & ", paste(rep(mlab[vE], length(ANG)), collapse = " & "), EOL),
    vapply(scE, function(s) paste0(
      if (cellf(E, "arm", E$scenario_id == s) %in% c(2, 3)) "no" else "yes",
      " & ", slab[[s]], " & ",
      paste(vapply(names(ANG), function(cl) paste(best(
        vapply(vE, function(v) gE(s, v, cl), numeric(1)), TRUE, 2,
        iqr = vapply(vE, function(v) gEi(s, v, cl), numeric(1))),
        collapse = " & "), character(1)), collapse = " & "), EOL), character(1))))
  if ("score_cc" %in% names(E))
    wr("tab_estim_cc.tex", tabE("score_cc", paste(
      "Canonical correlation between the estimated and true scores.",
      "The mean canonical correlation over the components, fitted on the true partition, $n=500$, $M=4$;",
      "higher is better, best per row in bold.",
      "Canonical correlation is invariant to any nonsingular linear transform of",
      "either set, so it does not depend on resolving the near-tied eigenvalues",
      "that make the eigenfunction angles ill-conditioned.", if (!is.null(ER)) IQR_SENTENCE else ""),
      "tab:estim-cc", lower = FALSE))
}

cat("\ndone\n")
