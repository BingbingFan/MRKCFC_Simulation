# Study D3: the EXTERNAL benchmark. Does a partition discriminate survival?
# This is the most useful question to ask of the partitions, because the
# evidence for preferring one partition over another is otherwise internal fit,
# and internal fit cannot distinguish a clinically meaningful grouping from a
# well-fitted one.
#
#   Rscript R/realdata_survival.R <results-dir> <combined-stem> \
#                                 [endpoint] [t1] [t2] [covariates]
#
# endpoint    OS (default) or RFS
# t1, t2      times for the time-dependent AUC, IN THE DATA'S OWN UNITS
# covariates  '+'-separated column names from dt_base; default NONE
#
# e.g.  Rscript R/realdata_survival.R .../results_rd .../rd OS 12 36
#
# THE MAIN RESULT IS THE ADJUSTED COX MODEL (2026-09-13), with the
# manuscript's covariate set, which maps onto dt_base_36 as
#
#   Age+Sex+Surgical_approach+Primary_site+Tumor_differentiation+AJCC_8th_ed_Stage+Lymph_node_yield+Mucinous_colloid_type+Lymphovascular_invasion+Perineural_invasion+Adjuvant_chemotherapy
#
# (age, sex, surgical approach, primary site, tumour differentiation, TNM
# stage, lymph node yield, mucinous type, lymphovascular invasion, perineural
# invasion, adjuvant chemotherapy; all complete, all factors but Age). An
# adjusted run writes `<stem>_survival_<endpoint>_adj.*`, reports the
# unadjusted hazard ratio beside the adjusted one as the original Table 4
# did, tests the cluster term by likelihood ratio GIVEN the covariates, and
# prints the Schoenfeld global test.
#
# ===========================================================================
# REWRITTEN 2026-09-09. The analysis is unchanged; the plumbing was all wrong.
#
# This file was first written on 2026-09-04 against a data contract invented
# before the file existed. Every accessor it used is gone:
# `rd_read`, `rd_validate`, `rd_subject_table`, and the `RD` list of column
# names -- including `RD$surv_time`, `RD$surv_event` and `RD$covariates`, none
# of which corresponded to anything in the delivered data.
#
# It now uses `rd_load()` and `rd_subjects()` from realdata_load.R, and takes
# the partitions from the combined array output.
#
# TWO THINGS ARE DELIBERATELY NOT GUESSED:
#
#   COVARIATES default to NONE. `dt_base_36` carries baseline columns whose
#   names and meanings this code has never seen, and inventing an adjustment set
#   is exactly the blind-writing that cost three submissions. The available
#   columns are PRINTED; choose them and pass them.
#
#   t1 and t2 default to quantiles of the observed event times, not to 3 and 5.
#   The previous version's 3 and 5 assumed YEARS; if `Survival_Time` is in
#   months -- as the 12/36-month window naming suggests -- those would have
#   evaluated discrimination at three and five MONTHS and reported it as
#   three- and five-year AUC. The follow-up range is printed so a unit mismatch
#   is visible rather than silent.
#
# THREE MEASURES, because they fail differently:
#   adjusted hazard ratio   effect size. Interpretable, but a single number and
#                           dependent on the proportional-hazards assumption.
#   Harrell's C             concordance over the whole follow-up. Assumption
#                           free, but averages over time and can hide a
#                           partition that separates early and converges later.
#   time-dependent AUC      discrimination AT t1 and t2. The one a clinician
#                           reads.
#
# THIS IS NOT A HYPOTHESIS TEST OF THE CLUSTERING. The partition was estimated
# from the same subjects, so the p-values are optimistic and must not be
# reported as evidence that clusters exist. What the comparison supports is
# RELATIVE: given that each method proposes a partition, whose separates
# survival best. Say that, and attach no significance claim to it.
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L)
  stop("usage: realdata_survival.R <results-dir> <combined-stem> [endpoint] [t1] [t2] [covariates]",
       call. = FALSE)
inp      <- args[1]
stem     <- args[2]
endpoint <- if (length(args) >= 3) toupper(args[3]) else "OS"
t1_arg   <- if (length(args) >= 4) as.numeric(args[4]) else NA_real_
t2_arg   <- if (length(args) >= 5) as.numeric(args[5]) else NA_real_
covars   <- if (length(args) >= 6 && nzchar(args[6]))
              trimws(strsplit(args[6], "[+,]")[[1]]) else character(0)
# Covariates are '+'-separated so that the list survives a job wrapper that
# splits its arguments on commas.

source(file.path("R", "study_c.R"))
source(file.path("R", "realdata_load.R"))
if (!requireNamespace("survival", quietly = TRUE))
  stop("the survival package is required", call. = FALSE)
# ATTACHED, not merely installed: timeROC::timeROC() calls Surv() unqualified,
# so with survival only installed every AUC came back NA, and the tryCatch
# around the call hid the reason for two runs (2026-09-12).
suppressMessages(library(survival))

cf <- paste0(stem, "_combined.rds")
if (!file.exists(cf))
  stop("run realdata_combine.R first; no ", cf, call. = FALSE)
C <- readRDS(cf)
if (is.na(C$K_star)) stop("the combine did not determine K.", call. = FALSE)

EP <- list(OS  = c(time = "Survival_Time",         event = "Death"),
           RFS = c(time = "Recurrence_free_Time",  event = "Recurrence"))
if (!endpoint %in% names(EP))
  stop("endpoint must be OS or RFS, not ", endpoint, call. = FALSE)
tcol <- EP[[endpoint]][["time"]]; ecol <- EP[[endpoint]][["event"]]

dat <- rd_load(window = C$window)
sub <- rd_subjects(dat, C$ids)
cat(sprintf("\navailable columns in dt_base_%d:\n  %s\n", C$window,
            paste(names(sub), collapse = ", ")))
miss <- setdiff(c(tcol, ecol, covars), names(sub))
if (length(miss))
  stop("not in dt_base: ", paste(miss, collapse = ", "),
       "\n  Pick from the list above.", call. = FALSE)

tt <- as.numeric(sub[[tcol]]); ee <- as.integer(sub[[ecol]])
keep <- is.finite(tt) & tt >= 0 & !is.na(ee)
if (any(!keep))
  cat(sprintf("dropping %d subject(s) with missing or negative %s\n",
              sum(!keep), tcol))
sub <- sub[keep, , drop = FALSE]; tt <- tt[keep]; ee <- ee[keep]
S <- survival::Surv(tt, ee)

cat(sprintf("\nStudy D3: endpoint %s (%s / %s), K = %d, bwK = %s\n",
            endpoint, tcol, ecol, C$K_star, format(C$bwK)))
cat(sprintf("  %d subjects, %d events (%.1f%%), follow-up %.3g to %.3g\n",
            length(tt), sum(ee), 100 * mean(ee), min(tt), max(tt)))
# t1/t2 FROM THE DATA, not from an assumed unit. Printed loudly: if these look
# like the wrong scale, the units are not what the caller assumed.
qt <- stats::quantile(tt[ee == 1], c(0.25, 0.50), na.rm = TRUE)
t1 <- if (is.na(t1_arg)) unname(qt[1]) else t1_arg
t2 <- if (is.na(t2_arg)) unname(qt[2]) else t2_arg
cat(sprintf("  time-dependent AUC at t = %.3g and %.3g%s\n", t1, t2,
            if (is.na(t1_arg)) "  [DEFAULTED to the 25th/50th percentile of event times -- pass clinically meaningful times instead]" else ""))
cat(sprintf("  adjustment: %s\n", if (length(covars))
            paste(covars, collapse = ", ") else "NONE (unadjusted)"))

# ------------------------------------------------------- the partitions -----
# Every partition available: the array's reported ones, plus each comparator if
# D2 has run, so all are scored on the same external criterion.
parts <- as.list(C$labels[setdiff(names(C$labels), "ID")])
agf <- file.path(dirname(stem), paste0(basename(stem), "_agreement.rds"))
if (file.exists(agf)) {
  parts <- readRDS(agf)$labels
  cat(sprintf("  scoring %d partitions from %s\n", length(parts), basename(agf)))
} else {
  cat(sprintf("  scoring %d partition(s) from the combine; run realdata_agreement.R\n",
              length(parts)),
      "  first to include the comparators.\n", sep = "")
}
parts <- lapply(parts, function(l) as.integer(l)[keep])

cov_rhs <- if (length(covars))
  paste("+", paste(sprintf("`%s`", covars), collapse = " + ")) else ""

score_one <- function(lab) {
  if (length(unique(lab)) < 2L) return(NULL)
  # Largest cluster as reference, so hazard ratios are comparable across methods
  # whose label numbering is arbitrary.
  ref <- names(sort(table(lab), decreasing = TRUE))[1]
  df  <- cbind(sub, cl = stats::relevel(factor(lab), ref = ref))
  fit <- tryCatch(survival::coxph(stats::as.formula(paste("S ~ cl", cov_rhs)),
                                  data = df), error = function(e) e)
  if (inherits(fit, "error")) return(list(error = conditionMessage(fit)))
  sm <- summary(fit)
  hr <- sm$conf.int[grep("^cl", rownames(sm$conf.int)), c(1, 3, 4), drop = FALSE]
  # Proportional hazards: Schoenfeld global test on the fitted model. Reported,
  # not acted on; under non-proportional hazards a Cox HR is a time average.
  ph_p <- tryCatch(survival::cox.zph(fit)$table["GLOBAL", "p"],
                   error = function(e) NA_real_)
  # The test of the PARTITION. With covariates, the model's overall likelihood
  # ratio test would test the covariates as well, so the cluster term is tested
  # by likelihood ratio against the covariate-only model.
  p_cl <- if (length(covars)) tryCatch({
    fit0 <- survival::coxph(stats::as.formula(paste("S ~", sub("^\\+ ", "", cov_rhs))),
                            data = df)
    stats::anova(fit0, fit)[2, "Pr(>|Chi|)"]
  }, error = function(e) NA_real_) else unname(sm$logtest["pvalue"])

  # Harrell's C on the CLUSTER VARIABLE ALONE, so it measures the partition
  # rather than being dominated by shared covariates. The same fit gives the
  # UNADJUSTED hazard ratio, reported beside the adjusted one as the original
  # Table 4 did (there, MKCFC's association held only unadjusted).
  cfit <- survival::coxph(S ~ cl, data = df)
  cc   <- survival::concordance(cfit)
  smu  <- summary(cfit)
  hr_unadj <- smu$conf.int[grep("^cl", rownames(smu$conf.int)), c(1, 3, 4), drop = FALSE]
  n_cl <- as.integer(table(df$cl)[sub("^cl", "", rownames(hr))])

  auc <- c(NA_real_, NA_real_)
  if (requireNamespace("timeROC", quietly = TRUE)) {
    tr <- tryCatch(timeROC::timeROC(T = tt, delta = ee,
                                    marker = stats::predict(cfit, type = "lp"),
                                    cause = 1, times = c(t1, t2), iid = FALSE),
                   error = function(e) {
                     cat(sprintf("  timeROC failed: %s\n", conditionMessage(e)))
                     NULL })
    if (!is.null(tr)) auc <- as.numeric(tr$AUC)
  }
  list(hr = hr, hr_unadj = hr_unadj, n_cl = n_cl, n_used = fit$n, ph_p = ph_p,
       c_index = cc$concordance, c_se = sqrt(cc$var),
       auc_t1 = auc[1], auc_t2 = auc[2],
       lrt_p = p_cl, k = length(unique(lab)))
}

res <- lapply(parts, score_one)

cat("\n=== survival discrimination by partition ===\n")
cat(sprintf("%-12s %3s %9s %9s %9s %9s\n", "method", "K", "C-index",
            sprintf("AUC(%.3g)", t1), sprintf("AUC(%.3g)", t2),
            if (length(covars)) "p(cl|cov)" else "LRT p"))
for (nm in names(res)) {
  r <- res[[nm]]
  if (is.null(r) || !is.null(r$error)) {
    cat(sprintf("%-12s  --  %s\n", nm,
                if (is.null(r)) "single cluster" else r$error))
    next
  }
  cat(sprintf("%-12s %3d %9.3f %9s %9s %9.2g\n", nm, r$k, r$c_index,
              if (is.na(r$auc_t1)) "---" else sprintf("%.3f", r$auc_t1),
              if (is.na(r$auc_t2)) "---" else sprintf("%.3f", r$auc_t2),
              r$lrt_p))
}
if (!requireNamespace("timeROC", quietly = TRUE))
  cat("\nNOTE: timeROC is not installed, so time-dependent AUC is not reported.\n",
      "Install it or drop the AUC claim -- do not substitute a\n",
      "different measure and call it AUC.\n", sep = "")

cat(sprintf("\n=== hazard ratios, largest cluster as reference (%s) ===\n",
            if (length(covars)) paste("adjusted for",
                                      paste(covars, collapse = ", "))
            else "UNADJUSTED"))
for (nm in names(res)) {
  r <- res[[nm]]
  if (is.null(r) || !is.null(r$error)) next
  cat(sprintf("  %s%s\n", nm,
              if (length(covars)) sprintf("   [n = %d; Schoenfeld global p = %s]", r$n_used,
                                          if (is.na(r$ph_p)) "NA" else sprintf("%.3f", r$ph_p))
              else ""))
  for (i in seq_len(nrow(r$hr)))
    cat(sprintf("    %-8s n=%4d  HR %6.2f  (%.2f, %.2f)%s\n", rownames(r$hr)[i], r$n_cl[i],
                r$hr[i, 1], r$hr[i, 2], r$hr[i, 3],
                if (length(covars)) sprintf("   unadjusted %6.2f  (%.2f, %.2f)",
                                            r$hr_unadj[i, 1], r$hr_unadj[i, 2], r$hr_unadj[i, 3])
                else ""))
}

cat("\nThe partitions were estimated on these same subjects, so the p-values are\n")
cat("optimistic and are NOT evidence that clusters exist. The comparison that is\n")
cat("supported is relative: whose proposed partition separates survival best.\n")

tab <- do.call(rbind, lapply(names(res), function(nm) {
  r <- res[[nm]]
  if (is.null(r) || !is.null(r$error)) return(NULL)
  data.frame(method = nm, K = r$k, endpoint = endpoint,
             cluster = rownames(r$hr), n_cluster = r$n_cl,
             hr = r$hr[, 1], hr_lo = r$hr[, 2], hr_hi = r$hr[, 3],
             hr_unadj = r$hr_unadj[, 1], hr_unadj_lo = r$hr_unadj[, 2],
             hr_unadj_hi = r$hr_unadj[, 3], ph_p = r$ph_p, n_used = r$n_used,
             c_index = r$c_index, c_se = r$c_se,
             auc_t1 = r$auc_t1, auc_t2 = r$auc_t2, t1 = t1, t2 = t2,
             p_cluster = r$lrt_p, adjusted = length(covars) > 0,
             stringsAsFactors = FALSE)
}))
# Adjusted runs get their own files, so the unadjusted benchmark is not
# overwritten by the Cox model the paper reports as its main result.
sfx <- if (length(covars)) "_adj" else ""
out <- file.path(dirname(stem), paste0(basename(stem), "_survival_", endpoint, sfx, ".rds"))
saveRDS(list(results = res, table = tab, endpoint = endpoint, t1 = t1, t2 = t2,
             n = length(tt), n_events = sum(ee), covariates = covars,
             bwK = C$bwK, K = C$K_star, when = Sys.time()), out)
if (!is.null(tab))
  utils::write.csv(tab, sprintf("%s_survival_%s%s.csv", stem, endpoint, sfx),
                   row.names = FALSE)
cat(sprintf("\nwrote %s\n", out))
