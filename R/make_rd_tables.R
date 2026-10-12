# Study D: every table and figure of the real-data section, from the local
# copies of the pipeline outputs. Run from study/:
#
#   Rscript R/make_rd_tables.R <rd-dir> <tables-dir> <figures-dir>
#   Rscript R/make_rd_tables.R ../HPC/V3/rd ../doc/manu/tables ../doc/manu/figures
#
# READS (all under <rd-dir>): rd_combined.rds, rd_restarts.csv, rd_agreement.rds,
# rd_stability.rds, rd_survival_{OS,RFS}_adj.rds, rd_fve.csv, and the same for
# the sensitivity arms rd_log_nK3, rd_log and rd_bwK020 where present; plus the cohort file
# through realdata_load.R for the baseline tables and the figures.
#
# WRITES tab_rd_*.tex into <tables-dir> and fig_rd_*.{pdf,png} into
# <figures-dir>, and prints every number the Section 3.2 text quotes.
#
# THE COX MODELS ARE THE MAIN RESULT (2026-09-13), with the manuscript's
# three models: Model 1 unadjusted; Model 2 age, sex and the
# subject's median postoperative CEA, CA19-9 and CA242; Model 3 the full set of
# eleven covariates. Hazard ratios are for the smaller class against the
# larger. Model 3 reproduces realdata_survival.R's adjusted run exactly (same
# data, same formula); it is refitted here so that Models 1 and 2 come from
# the same code path and the table is built from one place.
#
# CLASS NAMES. The K-centres classes are numbered arbitrarily by the loop.
# They are named here by their marker level -- the class with the higher
# median postoperative CEA is "elevated", the other "low" -- which is what the
# figure shows and what the text uses. The names are DERIVED, not typed.

suppressWarnings(suppressMessages({ source(file.path("R", "realdata_load.R")); library(survival) }))
source(file.path("R", "table_abbr.R"))   # abbr_footnote(): the abbreviation note
args   <- commandArgs(trailingOnly = TRUE)
RD     <- if (length(args) >= 1) args[1] else file.path("..", "HPC", "V3", "rd")
TABDIR <- if (length(args) >= 2) args[2] else file.path("..", "doc", "manu", "tables")
FIGDIR <- if (length(args) >= 3) args[3] else file.path("..", "doc", "manu", "figures")
dir.create(TABDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

VO   <- c("MRKCFC-R2", "MRKCFC", "MKCFC")          # the paper's order: proposed first
COMP <- c("mclust", "k-means", "teigen")           # comparators on the non-robust scores
# DISPLAY NAMES (2026-09-13): the two robust constructions are named by
# their second-stage scatter, -S for the spatial sign (the proposed method,
# internally "MRKCFC-R2") and -M for the second moment (the original
# construction, internally "MRKCFC"). The internal identifiers live in every
# result file and in the drivers and are NOT renamed; this map is the only
# place the paper's names come from, together with plot_fig1_traj.R.
DISP <- c(`MRKCFC-R2` = "MRKCFC-S", MRKCFC = "MRKCFC-M", MKCFC = "MKCFC",
          mclust = "mclust", `k-means` = "k-means", `k-means/rob` = "k-means (robust scores)",
          teigen = "teigen")
COV3 <- c("Age", "Sex", "Surgical_approach", "Primary_site", "Tumor_differentiation",
          "AJCC_8th_ed_Stage", "Lymph_node_yield", "Mucinous_colloid_type",
          "Lymphovascular_invasion", "Perineural_invasion", "Adjuvant_chemotherapy")

# ---------------------------------------------------------------- helpers ---
fmt  <- function(x, d = 2) ifelse(is.na(x), "---", formatC(x, format = "f", digits = d))
fmtp <- function(p) ifelse(is.na(p), "---", ifelse(p < 0.001, "$<$0.001", formatC(p, format = "f", digits = 3)))
hrci <- function(h, lo, hi) sprintf("%s (%s, %s)", fmt(h), fmt(lo), fmt(hi))
tex  <- function(s) gsub("%", "\\\\%", gsub("&", "\\\\&", s))
env  <- function(cap, lab, align, head, body, foot = NULL, size = "\\small", wide = FALSE) {
  width <- if (wide) "\\textheight" else "\\textwidth"
  c(paste0(if (wide) "\\begin{sidewaystable}[htbp]" else "\\begin{table}[htbp]", size),
    "\\centering", sprintf("\\caption{%s}", cap), sprintf("\\label{%s}", lab),
    sprintf("\\begin{tabular*}{%s}{@{\\extracolsep{\\fill}}%s@{}}", width, align),
    "\\toprule", head, "\\midrule", body, "\\bottomrule", "\\end{tabular*}",
    if (!is.null(foot)) sprintf("\\footnotetext{%s}", foot),
    abbr_footnote(c(cap, head, body, foot)),
    if (wide) "\\end{sidewaystable}" else "\\end{table}")
}
put  <- function(lines, f) { writeLines(lines, file.path(TABDIR, f)); cat("wrote", f, "\n") }
rdf  <- function(stem, what) { f <- file.path(RD, paste0(stem, what)); if (file.exists(f)) readRDS(f) else NULL }

# ------------------------------------------------------------------ data ---
dat  <- rd_load(window = 36L)
prep <- rd_prepare(dat, log = FALSE, verbose = FALSE, time = "warped")
B    <- rd_subjects(dat, prep$ids)              # baseline rows aligned to prep$ids
n    <- prep$n
# per-subject median postoperative level of each marker, raw units, from the
# marker frames (all visits of the 36-month window, not only the common ones).
# The window INCLUDES every subject's month-0 visit, which is the preoperative
# sample (its value equals Preoperative_* in dt_base_36; found 2026-10-09). These
# summaries and Model 2 are kept as they are, and the inclusion is stated in
# the paper and in the table notes (without it: episodic vs
# stable peak CEA 2.8 vs 2.9 instead of 9.8 vs 3.7; any CEA > 5: 27 vs 19 per
# cent instead of 68 vs 32).
post <- sapply(RD_MARKERS, function(m) {
  x <- dat$markers[[m]]; v <- tapply(x$yit, x$ID, stats::median, na.rm = TRUE)
  as.numeric(v[as.character(prep$ids)]) })
colnames(post) <- paste0("post_", RD_MARKERS)
B <- cbind(B, post)
months <- lapply(prep$Lt[[1]], function(u) rd_time_inverse(prep, u))

# ------------------------------------------------------------ partitions ---
C  <- readRDS(file.path(RD, "rd_combined.rds"))
A  <- readRDS(file.path(RD, "rd_agreement.rds"))
S  <- readRDS(file.path(RD, "rd_stability.rds"))
stopifnot(identical(C$ids, prep$ids), A$n == n)
labs <- A$labels                                # includes the comparators
Kstar <- C$K_star

# CLASS NAMES, by what distinguishes the two classes of each method
# (2026-09-13). The smaller class is named; the larger is "Low"
# or, for the proposed method, "Stable".
#   MRKCFC-R2 828/577: the same median CEA, CA19-9 and CA242 at every month,
#     but the 577 have a median peak CEA of 9.8 against 3.7 ng/mL and 68
#     against 32 percent with at least one CEA above 5 ng/mL (CA19-9 31 vs 12,
#     CA242 42 vs 18 percent): EPISODIC elevation, not a difference in level.
#   MRKCFC 1248/157: the 157 are higher at every month (CEA 3.0 -> 4.2 against
#     2.0 -> 2.2; 94 percent above 5 ng/mL): persistently ELEVATED.
#   MKCFC 1394/11, k-means 7, teigen 5: a handful of EXTREME early values.
#   mclust 1166/239: persistently ELEVATED, as MRKCFC's class.
NAME_SMALL <- c(`MRKCFC-R2` = "Episodic", MRKCFC = "Elevated", MKCFC = "Extreme",
                mclust = "Elevated", `k-means` = "Extreme", `k-means/rob` = "Extreme", teigen = "Extreme")
NAME_LARGE <- c(`MRKCFC-R2` = "Stable")
class_names <- lapply(names(labs), function(v) {
  g <- as.integer(labs[[v]]); cl <- sort(unique(g)); sz <- table(g)[as.character(cl)]
  if (length(cl) != 2L) return(setNames(paste("Class", cl), cl))
  small <- cl[which.min(sz)]
  setNames(ifelse(cl == small, NAME_SMALL[[v]], if (v %in% names(NAME_LARGE)) NAME_LARGE[[v]] else "Low"), cl)
})
names(class_names) <- names(labs)
# per-subject peak and any-above-reference-limit, the facts behind the names
ULN <- c(CEA = 5, CA199 = 37, CA242 = 20)      # conventional upper reference limits
for (m in RD_MARKERS) {
  x <- dat$markers[[m]]
  B[[paste0("peak_", m)]]  <- as.numeric(tapply(x$yit, x$ID, max, na.rm = TRUE)[as.character(prep$ids)])
  B[[paste0("above_", m)]] <- as.integer(B[[paste0("peak_", m)]] > ULN[[m]])
}
# a display label per class: "Elevated (n = 577)"
class_lab <- function(v, k) sprintf("%s ($n$ = %d)", class_names[[v]][as.character(k)], sum(labs[[v]] == k))

# --------------------------------------------------------- Cox models ------
S_OS  <- Surv(as.numeric(B$Survival_Time), as.integer(B$Death))
S_RFS <- Surv(as.numeric(B$Recurrence_free_Time), as.integer(B$Recurrence))
cox3 <- function(g, Sv) {
  g <- as.integer(g); ref <- names(sort(table(g), decreasing = TRUE))[1]
  df <- cbind(B, cl = stats::relevel(factor(g), ref = ref))
  fit <- function(rhs) {
    f <- coxph(stats::as.formula(paste("Sv ~ cl", rhs)), data = df)
    sm <- summary(f); ci <- sm$conf.int[grep("^cl", rownames(sm$conf.int)), c(1, 3, 4), drop = FALSE]
    pv <- sm$coefficients[grep("^cl", rownames(sm$coefficients)), "Pr(>|z|)"]
    list(hr = ci, p = pv, zph = tryCatch(cox.zph(f)$table["GLOBAL", "p"], error = function(e) NA_real_))
  }
  list(m1 = fit(""),
       m2 = fit("+ Age + Sex + post_CEA + post_CA199 + post_CA242"),
       m3 = fit(paste("+", paste(COV3, collapse = " + "))),
       ref = ref, other = setdiff(names(table(g)), ref))
}
cox_tab <- function(Sv, methods, lab, cap, endpoint) {
  body <- character(0); zph <- character(0)
  for (v in methods) {
    r <- cox3(labs[[v]], Sv)
    body <- c(body, sprintf("%s & & & & & & \\\\", DISP[[v]]),
              sprintf("\\quad %s & 1.00 (Reference) & & 1.00 (Reference) & & 1.00 (Reference) & \\\\",
                      class_lab(v, as.integer(r$ref))))
    for (i in seq_along(r$other)) {
      k <- as.integer(r$other[i])
      body <- c(body, sprintf("\\quad %s & %s & %s & %s & %s & %s & %s \\\\", class_lab(v, k),
                              hrci(r$m1$hr[i, 1], r$m1$hr[i, 2], r$m1$hr[i, 3]), fmtp(r$m1$p[i]),
                              hrci(r$m2$hr[i, 1], r$m2$hr[i, 2], r$m2$hr[i, 3]), fmtp(r$m2$p[i]),
                              hrci(r$m3$hr[i, 1], r$m3$hr[i, 2], r$m3$hr[i, 3]), fmtp(r$m3$p[i])))
    }
    zph <- c(zph, sprintf("%s %s", DISP[[v]], fmt(r$m3$zph, 3)))
    cat(sprintf("COX %s %-10s ref %s | M1 %s p=%s | M2 %s p=%s | M3 %s p=%s | PH(M3) p=%.3f\n", endpoint, v, r$ref,
                hrci(r$m1$hr[1,1], r$m1$hr[1,2], r$m1$hr[1,3]), fmtp(r$m1$p[1]),
                hrci(r$m2$hr[1,1], r$m2$hr[1,2], r$m2$hr[1,3]), fmtp(r$m2$p[1]),
                hrci(r$m3$hr[1,1], r$m3$hr[1,2], r$m3$hr[1,3]), fmtp(r$m3$p[1]), r$m3$zph))
  }
  head <- c("& \\multicolumn{2}{c}{Model 1} & \\multicolumn{2}{c}{Model 2} & \\multicolumn{2}{c}{Model 3} \\\\",
            "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
            "Multivariate trajectory class & HR (95\\% CI) & $P$ & HR (95\\% CI) & $P$ & HR (95\\% CI) & $P$ \\\\")
  foot <- paste0("Primary analysis: original scale, transformed time axis, $K=2$. ",
                 "Model 1 is unadjusted. Model 2 adjusts for age, sex and the subject's median postoperative CEA, CA19-9 and CA242, each the median over the 36-month window with the preoperative measurement included. ",
                 "Model 3 adjusts for age, sex, surgical approach, primary site, tumor differentiation, TNM stage, lymph node yield, ",
                 "mucinous type, lymphovascular invasion, perineural invasion and adjuvant chemotherapy. Hazard ratios are for the ",
                 "class named against the reference class; $P$ is the Wald test of that contrast. Schoenfeld global test of ",
                 "proportional hazards in Model 3: ", paste(zph, collapse = "; "), ". The partitions were estimated on the same ",
                 "subjects, so the $P$ values are optimistic and support only the comparison between methods.")
  put(env(cap, lab, "lcccccc", head, body, foot, size = "\\footnotesize", wide = TRUE),
      sprintf("tab_rd_cox_%s.tex", tolower(endpoint)))
}
# Captions are short titles (2026-10-10); what they used to carry is
# the first sentence of each note.
cox_tab(S_OS,  c(VO, COMP), "table4", "Trajectory classes and overall survival", "OS")
# The recurrence-free table is no longer part of the manuscript (author, 2026-10-11);
# the call stays here so it can be switched back on:
# cox_tab(S_RFS, c(VO, COMP), "tab:rd-cox-rfs", "Trajectory classes and recurrence-free survival", "RFS")

# Kaplan-Meier facts for the text: deaths, median survival, 5-year survival, log-rank
km_facts <- function(v) {
  g <- as.integer(labs[[v]]); cl <- sort(unique(g))
  sf <- survfit(S_OS ~ factor(g)); s5 <- summary(sf, times = 60, extend = TRUE)
  lr <- survdiff(S_OS ~ factor(g))$pvalue
  for (i in seq_along(cl)) {
    k <- cl[i]; sel <- g == k
    cat(sprintf("KM %-10s %-8s n=%4d deaths %3d (%.1f%%) median OS %.1f mo | 5-yr survival %.1f%% (%.1f, %.1f) | log-rank p=%.2g\n",
                v, class_names[[v]][as.character(k)], sum(sel), sum(B$Death[sel]), 100 * mean(B$Death[sel]),
                stats::median(B$Survival_Time[sel]), 100 * s5$surv[i], 100 * s5$lower[i], 100 * s5$upper[i], lr))
  }
}
for (v in c(VO, "mclust")) km_facts(v)

# ---------------------------------------------------- baseline tables -------
num_row <- function(label, x, g, kind = c("meansd", "mediqr")) {
  kind <- match.arg(kind); cl <- sort(unique(g))
  cells <- vapply(cl, function(k) { v <- x[g == k]
    if (kind == "meansd") sprintf("%.1f (%.1f)", mean(v, na.rm = TRUE), stats::sd(v, na.rm = TRUE))
    else sprintf("%.1f [%.1f, %.1f]", stats::median(v, na.rm = TRUE), stats::quantile(v, .25, na.rm = TRUE), stats::quantile(v, .75, na.rm = TRUE)) }, "")
  p <- if (kind == "meansd") stats::t.test(x ~ g)$p.value else stats::wilcox.test(x ~ g)$p.value
  sprintf("%s & %s & %s \\\\", label, paste(cells, collapse = " & "), fmtp(p))
}
LEVEL_LAB <- c(Well_moderate = "Well or moderate", Poor_undifferentiated = "Poor or undifferentiated",
               Laparoscopic_resection = "Laparoscopic resection", Open_resection = "Open resection")
lab_of <- function(lv) if (lv %in% names(LEVEL_LAB)) LEVEL_LAB[[lv]] else tex(gsub("_", " ", lv))
cat_rows <- function(label, x, g, level = NULL) {
  cl <- sort(unique(g)); x <- factor(x)
  p <- tryCatch(stats::chisq.test(table(x, g))$p.value, error = function(e) NA_real_)
  if (!is.null(level)) {           # a yes/no variable: one row, the "yes" level
    cells <- vapply(cl, function(k) sprintf("%d (%.1f)", sum(x[g == k] == level), 100 * mean(x[g == k] == level)), "")
    return(sprintf("%s & %s & %s \\\\", label, paste(cells, collapse = " & "), fmtp(p)))
  }
  out <- sprintf("%s & %s & %s \\\\", label, paste(rep("", length(cl)), collapse = " & "), fmtp(p))
  for (lv in levels(x)) {
    cells <- vapply(cl, function(k) sprintf("%d (%.1f)", sum(x[g == k] == lv), 100 * mean(x[g == k] == lv)), "")
    out <- c(out, sprintf("\\quad %s & %s & \\\\", lab_of(lv), paste(cells, collapse = " & ")))
  }
  out
}
baseline_tab <- function(v, lab, cap) {
  g <- as.integer(labs[[v]]); cl <- sort(unique(g))
  head <- sprintf("Variables & %s & $P$ \\\\", paste(vapply(cl, function(k) class_lab(v, k), ""), collapse = " & "))
  body <- c(sprintf("$N$ & %s & \\\\", paste(table(g)[as.character(cl)], collapse = " & ")),
    num_row("Age, years", B$Age, g, "meansd"),
    cat_rows("Male, n (\\%)", B$Sex, g, level = "Male"),
    num_row("BMI, kg/m$^2$", B$BMI, g, "meansd"),
    num_row("Preoperative CEA, ng/mL", B$Preoperative_CEA, g, "mediqr"),
    num_row("Preoperative CA19-9, U/mL", B$Preoperative_CA199, g, "mediqr"),
    num_row("Preoperative CA242, U/mL", B$Preoperative_CA242, g, "mediqr"),
    num_row("Postoperative CEA, ng/mL", B$post_CEA, g, "mediqr"),
    num_row("Postoperative CA19-9, U/mL", B$post_CA199, g, "mediqr"),
    num_row("Postoperative CA242, U/mL", B$post_CA242, g, "mediqr"),
    num_row("Peak postoperative CEA, ng/mL", B$peak_CEA, g, "mediqr"),
    num_row("Peak postoperative CA19-9, U/mL", B$peak_CA199, g, "mediqr"),
    num_row("Peak postoperative CA242, U/mL", B$peak_CA242, g, "mediqr"),
    cat_rows("Any CEA $>$ 5 ng/mL, n (\\%)", B$above_CEA, g, level = "1"),
    cat_rows("Any CA19-9 $>$ 37 U/mL, n (\\%)", B$above_CA199, g, level = "1"),
    cat_rows("Any CA242 $>$ 20 U/mL, n (\\%)", B$above_CA242, g, level = "1"),
    cat_rows("Death, n (\\%)", B$Death, g, level = "1"),
    num_row("Survival time, months", B$Survival_Time, g, "mediqr"),
    cat_rows("Primary site, n (\\%)", B$Primary_site, g),
    cat_rows("Surgical approach, n (\\%)", B$Surgical_approach, g),
    cat_rows("Tumor differentiation, n (\\%)", B$Tumor_differentiation, g),
    cat_rows("Perineural invasion, n (\\%)", B$Perineural_invasion, g, level = "Yes"),
    cat_rows("TNM stage, n (\\%)", B$AJCC_8th_ed_Stage, g),
    cat_rows("Lymph node yield $>$12, n (\\%)", B$Lymph_node_yield, g, level = "More_than_12"),
    cat_rows("Mucinous (colloid) type, n (\\%)", B$Mucinous_colloid_type, g, level = "Yes"),
    cat_rows("Lymphovascular invasion, n (\\%)", B$Lymphovascular_invasion, g, level = "Yes"),
    cat_rows("Adjuvant chemotherapy, n (\\%)", B$Adjuvant_chemotherapy, g, level = "Yes"))
  foot <- paste("The multivariate trajectory classes of the primary analysis, $K=2$.",
                "Data are mean (SD) or median [IQR] for numeric variables and n (\\%) for categorical variables;",
                "postoperative levels are each subject's median, and peaks each subject's maximum, over all measurements in the 36-month window, the preoperative measurement at month 0 included;",
                "the reference limits are 5 ng/mL for CEA, 37 U/mL for CA19-9 and 20 U/mL for CA242. $P$ from the $t$ test (mean (SD) rows),",
                "the Wilcoxon rank-sum test (median [IQR] rows) or the $\\chi^2$ test.")
  put(env(cap, lab, paste0("l", strrep("c", length(cl)), "c"), head, body, foot, size = "\\footnotesize"),
      sprintf("tab_rd_baseline_%s.tex", gsub("[^A-Za-z0-9]", "", v)))
}
baseline_tab("MRKCFC-R2", "table_s2", "Cohort characteristics by MRKCFC-S class")
baseline_tab("MRKCFC",    "table_s3", "Cohort characteristics by MRKCFC-M class")

# ------------------------------------------------ K selection + outcomes ---
R  <- read.csv(file.path(RD, "rd_restarts.csv"), stringsAsFactors = FALSE)
sel <- C$selection
body <- character(0)
for (v in VO) for (K in sort(unique(sel$K))) {
  s <- sel[sel$variant == v & sel$K == K, ]; r <- R[R$variant == v & R$K == K, ]
  st <- table(factor(r$status, levels = c("converged", "maxiter", "degenerate", "error")))
  sz <- r$sizes[r$restart == s$restart]
  body <- c(body, sprintf("%s & %d & %d / %d / %d / %d & %s & %s & %s & %s \\\\", if (K == min(sel$K)) DISP[[v]] else "", K,
                          st["converged"], st["maxiter"], st["degenerate"], st["error"], gsub("/", " / ", sz),
                          fmt(s$silhouette, 3), fmt(s$silhouette_nonrobust, 3), fmt(s$CH, 3)))
}
put(env("Restart outcomes and the $K$-selection criteria",
        "tab:rd-kselect", "llcccccc",
        c("Method & $K$ & converged / cap / guard / smoother & class sizes of the best restart & silhouette & silhouette (non-robust scores) & CH \\\\"),
        body,
        paste("Primary analysis, ten random restarts per cell.",
              "Outcomes: converged, a repeated labeling; cap, the 100-iteration cap without one; guard, a class fell below ten curves;",
              "smoother, the covariance smoother could not be fitted to a class. The best restart minimizes the objective of Equation~(\\ref{EQ15}).",
              "The silhouette and CH are computed on the scores of one global robust multivariate FPCA (the paper's criterion) and, for comparison,",
              "on the scores of the non-robust global fit; $K$ is chosen by the silhouette on the former."),
        size = "\\footnotesize", wide = TRUE), "tab_rd_kselect.tex")

# ------------------------------------------------------------- stability ---
body <- character(0)
for (v in VO) {
  p <- S$per_variant[[v]]; if (is.null(p)) next
  jac <- vapply(seq_len(ncol(p$jaccard)), function(k) fmt(p$jaccard_mean[k]), "")
  sub <- if (!is.null(S$subsample) && v %in% S$subsample$variant) {
    ss <- S$subsample[S$subsample$variant == v, ]; paste(fmt(ss$jaccard_mean), collapse = " / ") } else "---"
  # Agreement with the reported partition over the OTHER restarts: the reported
  # partition is one of the ten and would otherwise contribute its own 1.00.
  tr <- p$ari_to_reported; self <- which(tr >= 1 - 1e-9)[1]; if (!is.na(self)) tr <- tr[-self]
  body <- c(body, sprintf("%s & %s [%s, %s] & %s [%s, %s] & %d\\%% & %s & %s \\\\", DISP[[v]],
                          fmt(stats::median(p$pairwise_ari)), fmt(min(p$pairwise_ari)), fmt(max(p$pairwise_ari)),
                          fmt(stats::median(tr)), fmt(min(tr)), fmt(max(tr)),
                          round(100 * p$modal_share), paste(jac, collapse = " / "), sub))
  cat(sprintf("STAB %-10s pairwise ARI median %.3f [%.3f, %.3f]; vs reported %.3f; modal %d%%; Jaccard %s; subsample %s\n", v,
              stats::median(p$pairwise_ari), min(p$pairwise_ari), max(p$pairwise_ari), stats::median(p$ari_to_reported),
              round(100 * p$modal_share), paste(jac, collapse = "/"), sub))
}
put(env("Stability of the $K=2$ partitions",
        "tab:rd-stability", "lccccc",
        "Method & Pairwise aRI, median [min, max] & aRI to reported, median [min, max] & Modal share & Jaccard, restarts (class 1 / 2) & Jaccard, subsamples (class 1 / 2) \\\\",
        body,
        paste("Primary analysis: ten random initializations and 100 subsamples of 80 per cent of the patients.",
              "The agreement with the reported partition is over the other nine restarts.",
              "Cluster-wise Jaccard follows Hennig (2007): the mean, over restarts or over subsample draws, of the",
              "Jaccard similarity between each reported class and its best-matching class in the refit. The subsample refits start from the",
              "reported partition restricted to the draw, so they remove initialization noise and are optimistic relative to a full bootstrap;",
              "they were run for the proposed method only. Classes are numbered as in the reported partition."),
        size = "\\footnotesize", wide = TRUE), "tab_rd_stability.tex")

# ------------------------------------------------------------- agreement ---
M <- A$agreement; nm <- rownames(M); dn <- DISP[nm]
head <- sprintf("& %s & own restarts \\\\", paste(dn, collapse = " & "))
body <- vapply(seq_along(nm), function(i) {
  own <- if (nm[i] %in% names(S$per_variant)) fmt(stats::median(S$per_variant[[nm[i]]]$pairwise_ari)) else "---"
  cells <- vapply(seq_along(nm), function(j) if (j < i) fmt(M[i, j]) else if (j == i) "1" else "", "")
  sprintf("%s & %s & %s \\\\", dn[i], paste(cells, collapse = " & "), own) }, "")
put(env("Agreement between methods",
        "tab:rd-agreement", paste0("l", strrep("c", length(nm)), "c"), head, body,
        paste("Adjusted Rand index between the $K=2$ partitions of the primary analysis; the last column is each method's median agreement",
              "with its own random restarts. mclust, $k$-means and teigen are fitted to the scores of the non-robust global multivariate FPCA, as in the simulation study;",
              "$k$-means (robust scores) to the scores of the robust global fit. lcmm is not included: its multivariate form needs a per-marker",
              "specification that was not undertaken, and it did not converge on the simulation design (Section~\\ref{subsec31}).",
              "The matrix is descriptive: methods sharing a first stage agree partly because they share it."),
        size = "\\footnotesize", wide = TRUE), "tab_rd_agreement.tex")

# ------------------------------------------------------------------- FVE ---
Fv <- read.csv(file.path(RD, "rd_fve.csv"), stringsAsFactors = FALSE)
body <- character(0)
for (sc in unique(Fv$scale)) for (br in unique(Fv$branch)) for (m in RD_MARKERS) {
  v <- Fv$FVE[Fv$scale == sc & Fv$branch == br & Fv$marker == m]
  body <- c(body, sprintf("%s & %s & %s & %s \\\\", sc, br, sub("CA199", "CA19-9", m), paste(fmt(v, 3), collapse = " & ")))
}
put(env("Fraction of variance explained on the cohort",
        "tab:rd-fve", "lllcccccc", "Scale & Branch & Marker & $M=1$ & $M=2$ & $M=3$ & $M=4$ & $M=5$ & $M=6$ \\\\", body,
        paste("Cumulative fraction of variance explained by the first $M$ components, per marker, scale and branch.",
              "Eigenvalues of the estimated covariance surface with negative values set to zero, as in Supplementary Table~\\ref{tab:fve-selection}.",
              "The rule retains the smallest $M$ with cumulative FVE $\\ge 0.90$ on the worst marker and branch: $M=4$ on the original scale, $M=3$ on the log scale."),
        size = "\\footnotesize"), "tab_rd_fve.tex")

# ---------------------------------------------------- sensitivity arms -----
sens_row <- function(stem, name) {
  Cs <- rdf(stem, "_combined.rds"); Ss <- rdf(stem, "_stability.rds"); Vs <- rdf(stem, "_survival_OS_adj.rds")
  if (is.null(Cs)) return(NULL)
  out <- character(0)
  for (v in VO) {
    ks <- Cs$selection[Cs$selection$variant == v, ]; kstar <- ks$K[which.max(ks$silhouette)]
    sz <- Cs$info$sizes[Cs$info$variant == v & Cs$info$K == Cs$K_star & Cs$info$restart == Cs$best$restart[Cs$best$variant == v & Cs$best$K == Cs$K_star]]
    ari <- if (!is.null(Ss$per_variant[[v]])) { pa <- Ss$per_variant[[v]]$pairwise_ari
      sprintf("%s [%s, %s]", fmt(stats::median(pa)), fmt(min(pa)), fmt(max(pa))) } else "---"
    # Wald P of the class term from the unrounded estimate and CI, as Table 5
    # prints it; p_cluster in the rds is the likelihood-ratio test (2026-10-10).
    hr <- if (!is.null(Vs)) { t <- Vs$table[Vs$table$method == v, ]
      z  <- log(t$hr[1]) / ((log(t$hr_hi[1]) - log(t$hr_lo[1])) / (2 * stats::qnorm(0.975)))
      pw <- 2 * stats::pnorm(-abs(z))
      sprintf("%s, %s", hrci(t$hr[1], t$hr_lo[1], t$hr_hi[1]),
              if (pw < 0.001) "$P<0.001$" else sprintf("$P$ = %s", formatC(pw, format = "f", digits = 3))) } else "---"
    out <- c(out, sprintf("%s & %s & %d & %s & %s & %s \\\\", if (v == VO[1]) name else "", DISP[[v]], kstar, gsub("/", " / ", sz), ari, hr))
    cat(sprintf("SENS %-22s %-10s K*=%d sizes %s ARI %s HR(M3) %s\n", name, v, kstar, sz, ari, hr))
  }
  out
}
# Row order as the text presents them (2026-09-14): the log scale at
# the M the rule selects there first, then one component more, then the bandwidth.
SENS_ARMS <- list(c("rd_log_nK3", "Log scale, $M=3$"), c("rd_log", "Log scale, $M=4$"),
                  c("rd_bwK020", "Fixed bandwidth 0.20"))
body <- unlist(lapply(SENS_ARMS, function(s) sens_row(s[1], s[2])))
put(env("Sensitivity analyses",
        "tab:rd-sens", "llcccc",
        "Analysis & Method & selected $K$ & class sizes at $K=2$ & pairwise aRI among restarts, median [min, max] & adjusted HR (95\\% CI), Model 3, overall survival \\\\",
        body,
        paste("The log-transformed markers at two component counts, and a fixed covariance bandwidth.",
              "Each analysis repeats the primary procedure with one change: every marker log-transformed before normalization (at the $M=3$ the",
              "variance-explained rule selects on that scale, and at $M=4$, one component more), or the covariance bandwidth fixed at 0.20 of the transformed axis instead of",
              "chosen by generalized cross-validation. The hazard ratio is for the smaller class against the larger, adjusted as Model 3 of Table~\\ref{table4};",
              "$P$ is the Wald test of the class term, as in that table."),
        size = "\\footnotesize", wide = TRUE), "tab_rd_sens.tex")

# ------------------------------------------------- cross-arm overlap -------
# What each sensitivity arm's classes have to do with the primary analysis'
# classes, patient by patient: the adjusted Rand index between the arm's
# reported K = 2 partition and the primary one for the same method, and how
# many of the arm's smaller class belong to the primary's smaller class. The
# Section 3.2 sensitivity paragraph quotes these ("MRKCFC-S's class of 407
# contains 332 of the 577 episodic patients"; "MKCFC's class of 485, which
# has two patients in common with its eleven"). Same ARI function as the
# agreement tables (metrics.R: mclust::adjustedRandIndex). An arm whose
# combined object lists the subjects in a different order is skipped, not
# realigned.
xarm <- function(stem, name) {
  Cs <- rdf(stem, "_combined.rds"); if (is.null(Cs)) return(invisible(NULL))
  if (!identical(Cs$ids, C$ids)) { cat(sprintf("XARM %-20s subject order differs from the primary; skipped\n", name)); return(invisible(NULL)) }
  smaller <- function(g) { g <- as.integer(g); tb <- table(g); g == as.integer(names(tb)[which.min(tb)]) }
  for (v in VO) {
    a <- as.integer(labs[[v]]); b <- as.integer(Cs$labels[[v]]); sa <- smaller(a); sb <- smaller(b)
    cat(sprintf("XARM %-20s %-10s ARI with primary %5.2f | smaller class %4d, of whom %4d in the primary's smaller class of %4d (%3.0f%%)\n",
                name, v, mclust::adjustedRandIndex(a, b), sum(sb), sum(sa & sb), sum(sa), 100 * sum(sa & sb) / sum(sb)))
  }
}
for (s in SENS_ARMS) xarm(s[1], s[2])

# --------------------------------------------------------------- figures ---
col2 <- c(Low = "#1f77b4", Stable = "#1f77b4", Elevated = "#d62728", Episodic = "#d62728", Extreme = "#d62728")
traj_fig <- function(file, methods) {
  grid_mo <- seq(0, 36, by = 1)
  dev_fun <- if (grepl("pdf$", file)) function() grDevices::pdf(file, width = 9, height = 2.9 * length(methods))
             else function() grDevices::png(file, width = 9, height = 2.9 * length(methods), units = "in", res = 200)
  dev_fun(); on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(length(methods), 3), mar = c(4, 4.2, 2.5, 1), oma = c(0, 0, 0, 0))
  for (v in methods) {
    g <- as.integer(labs[[v]]); cl <- sort(unique(g))
    for (m in RD_MARKERS) {
      x <- dat$markers[[m]]; x <- x[x$ID %in% prep$ids, ]
      x$cl <- g[match(x$ID, prep$ids)]
      # pointwise median and IQR of the raw values in 3-month bins, per class
      x$bin <- pmin(35, floor(x$time / 3) * 3) + 1.5
      # pointwise quartiles and the 90th percentile per class: the proposed
      # method's classes differ in their excursions, which only the upper tail shows
      st <- lapply(cl, function(k) { xx <- x[x$cl == k, ]
        b <- sort(unique(xx$bin)); q <- t(vapply(b, function(bb) stats::quantile(xx$yit[xx$bin == bb], c(.25, .5, .75, .9), na.rm = TRUE), numeric(4)))
        list(b = b, q = q) })
      ylim <- range(unlist(lapply(st, function(s) s$q)), na.rm = TRUE); ylim[1] <- max(ylim[1], 0.1)
      graphics::plot(NA, xlim = c(0, 36), ylim = ylim, log = "y", xlab = "Months after surgery",
                     ylab = sprintf("%s (%s)", sub("CA199", "CA19-9", m), if (m == "CEA") "ng/mL" else "U/mL"),
                     main = sprintf("%s", DISP[[v]]), las = 1)
      for (i in seq_along(cl)) { s <- st[[i]]; nmk <- class_names[[v]][as.character(cl[i])]
        cc <- if (nmk %in% names(col2)) col2[[nmk]] else grDevices::rainbow(length(cl))[i]
        graphics::polygon(c(s$b, rev(s$b)), c(s$q[, 1], rev(s$q[, 3])), col = grDevices::adjustcolor(cc, 0.18), border = NA)
        graphics::lines(s$b, s$q[, 2], col = cc, lwd = 2.2)
        graphics::lines(s$b, s$q[, 4], col = cc, lwd = 1.2, lty = 2) }
      if (m == RD_MARKERS[1]) graphics::legend("topleft", bty = "n", lwd = 2.2, cex = 0.85,
        col = vapply(cl, function(k) { nmk <- class_names[[v]][as.character(k)]; if (nmk %in% names(col2)) col2[[nmk]] else "grey40" }, ""),
        legend = vapply(cl, function(k) sprintf("%s (n = %d)", class_names[[v]][as.character(k)], sum(g == k)), ""))
    }
  }
}
km_fig <- function(file, methods) {
  dev_fun <- if (grepl("pdf$", file)) function() grDevices::pdf(file, width = 3.1 * length(methods), height = 3.7)
             else function() grDevices::png(file, width = 3.1 * length(methods), height = 3.7, units = "in", res = 200)
  dev_fun(); on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(1, length(methods)), mar = c(4.2, 4.2, 2.6, 0.8),
                cex.lab = 1.1, cex.axis = 1, cex.main = 1.15)
  for (v in methods) {
    g <- as.integer(labs[[v]]); cl <- sort(unique(g))
    sf <- survfit(S_OS ~ factor(g)); lr <- survdiff(S_OS ~ factor(g))$pvalue
    cols <- vapply(cl, function(k) { nmk <- class_names[[v]][as.character(k)]; if (nmk %in% names(col2)) col2[[nmk]] else "grey40" }, "")
    graphics::plot(sf, col = cols, lwd = 2.2, xlab = "Months after surgery", ylab = "Overall survival probability",
                   main = DISP[[v]], las = 1, xlim = c(0, 96), ylim = c(0.5, 1), mark.time = FALSE)
    graphics::legend("bottomleft", bty = "n", cex = 0.95,
      legend = c(vapply(cl, function(k) sprintf("%s (n = %d)", class_names[[v]][as.character(k)], sum(g == k)), ""),
                 sprintf("log-rank P %s", if (lr < 0.001) "< 0.001" else sprintf("= %.3f", lr))),
      col = c(cols, NA), lwd = c(rep(2.2, length(cl)), NA))
  }
}
# The diagnostic panel Section 3.2 cites: the marginal distribution of every
# observed value, original and log scale, with sample skewness and excess
# kurtosis. Marginal only, and said so in the text.
diag_fig <- function(file) {
  dev_fun <- if (grepl("pdf$", file)) function() grDevices::pdf(file, width = 9, height = 5.6)
             else function() grDevices::png(file, width = 9, height = 5.6, units = "in", res = 200)
  dev_fun(); on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(2, 3), mar = c(4, 4.2, 2.5, 1))
  sk <- function(v) mean((v - mean(v))^3) / stats::sd(v)^3
  ku <- function(v) mean((v - mean(v))^4) / stats::sd(v)^4 - 3
  for (row in c("original", "log")) for (m in RD_MARKERS) {
    x <- dat$markers[[m]]; v <- x$yit[x$ID %in% prep$ids & !is.na(x$yit)]
    if (row == "log") v <- log(v)
    graphics::hist(v, breaks = 60, freq = FALSE, col = "grey85", border = "grey60", las = 1,
                   main = sprintf("%s, %s scale", sub("CA199", "CA19-9", m), row),
                   xlab = if (row == "log") sprintf("log %s", sub("CA199", "CA19-9", m)) else sprintf("%s (%s)", sub("CA199", "CA19-9", m), if (m == "CEA") "ng/mL" else "U/mL"))
    graphics::curve(stats::dnorm(x, mean(v), stats::sd(v)), add = TRUE, col = "#d62728", lwd = 1.5)
    graphics::legend("topright", bty = "n", cex = 0.85,
      legend = c(sprintf("n = %d", length(v)), sprintf("skewness %.1f", sk(v)), sprintf("excess kurtosis %.0f", ku(v))))
    cat(sprintf("DIAG %-6s %-8s n=%d skewness %.2f excess kurtosis %.1f\n", m, row, length(v), sk(v), ku(v)))
  }
}
# FIGURE 1 AS THE PAPER SHOWS IT (2026-09-13): one penalised-spline
# smooth per class, no points, log axis, one row per method with its legend at
# the right. Defined ONCE, in R/plot_fig1_traj.R, which also runs standalone;
# the median-and-quantile version above is kept as fig_rd_traj_quantiles for
# the record and for the text's statement that the proposed method's classes
# share the same median at every month.
source(file.path("R", "plot_fig1_traj.R"))
stopifnot(identical(derive_class_names(labs[VO]), class_names[VO]))   # one naming rule
for (ext in c("pdf", "png")) {
  traj_fig(file.path(FIGDIR, paste0("fig_rd_traj_quantiles.", ext)), VO)
  km_fig(file.path(FIGDIR, paste0("fig_rd_km.", ext)), c("MRKCFC-R2", "MRKCFC", "MKCFC"))
  diag_fig(file.path(FIGDIR, paste0("fig_rd_diag.", ext)))
}
fig1_gam(fig1_data(dat$markers, prep$ids, labs, class_names, VO), VO,
         file = file.path(FIGDIR, "fig_rd_traj"))
cat("figures written to", FIGDIR, "\n")
cat("\nclass names by median postoperative CEA:\n"); for (v in names(class_names)) cat(sprintf("  %-12s %s\n", v, paste(sprintf("%s=%s", names(class_names[[v]]), class_names[[v]]), collapse = " ")))
