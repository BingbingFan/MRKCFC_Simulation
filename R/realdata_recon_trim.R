# Study D: the reconstruction error of the three reported K = 2 partitions,
# observation by observation, under the COMMON reconstruction of Equation (14).
#
#   Rscript R/realdata_recon_trim.R <rd-dir> <out-dir> [smoke]
#
#   <rd-dir>   the directory holding rd_combined.rds of the primary arm (the
#              combine writes it beside its other outputs; pass that directory)
#   <out-dir>  receives rd_recon_trim.csv and rd_recon_trim_<variant>.rds
#   smoke      optional: a 150-subject subset, plumbing only, no validation
#
# WHY THIS EXISTS (2026-10-10). Section 3.2 says the total squared
# reconstruction error is lowest for MKCFC's 1394/11 partition because 87 to
# 91 per cent of it comes from the one per cent of observations with the
# largest errors, and that with those excluded the mean squared error is 2.10,
# 1.61 and 1.90 (MKCFC, MRKCFC-S, MRKCFC-M) and the ordering reverses. Those
# figures were first computed on 2026-09-14 and not written to a file; the task
# checkpoints hold labels and the total, not the observations.
# This script recomputes them from the same reconstruction the checkpoints
# used, so that every number in that sentence has a file behind it.
#
# THE RECONSTRUCTION IS THE CHECKPOINTS' OWN. realdata_fit_task.R's
# write_ckpt() scores every partition with
#   r <- reconstruct_given_partition(prep$Ly, prep$Lt, labels, cfg)
#   sum((r$y - r$yhat)^2)
# which methods.R defines under the config's recon_robust = FALSE and
# recon_scatter = "moment": one non-robust first stage and one moment second
# stage for every variant alike (Equation (14) applied class by class). The
# total this script gets must therefore equal the objective realdata_combine.R
# selected on, and that equality is printed as the validation.
#
# THE TRIM. "Top q" means the ceiling(q * N) observations with the largest
# squared error, N the number of observations over all markers (the three
# markers share one schedule, so N = 3 x visits). The share is their sum over
# the total; the trimmed MSE is the mean over the others. q = 0.01 is what the
# text quotes; 0.005, 0.02 and 0.05 are printed beside it so the reader of the
# log can see how sharp the 1 per cent is. The file also records how many
# subjects the top 1 per cent come from and how many of those sit in the
# smaller class, which is the mechanism the sentence describes.
#
# Three partitions, three independent reconstructions: mclapply over them
# with up to three cores (SLURM_CPUS_PER_TASK, or 1). Nothing here can use more.
#
# For example, with the combined file in results/rd:
#
#   Rscript R/realdata_recon_trim.R results/rd results/rd

source(file.path("R", "study_c.R"))       # load_config, build_optns and, via
                                          # methods.R, reconstruct_given_partition
source(file.path("R", "realdata_load.R"))

args  <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: Rscript R/realdata_recon_trim.R <rd-dir> <out-dir> [smoke]")
RD    <- args[1]; OUT <- args[2]
smoke <- length(args) >= 3L && identical(args[3], "smoke")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

C <- readRDS(file.path(RD, "rd_combined.rds"))
if (!identical(as.integer(C$K_star), 2L)) stop("rd_combined.rds selects K = ", C$K_star, ", not 2")
if (isTRUE(C$log)) stop("rd_combined.rds is a log arm; the sentence is about the primary analysis")
if (!identical(as.character(C$time), "warped")) stop("rd_combined.rds is not on the warped axis")
VARIANTS <- c("MKCFC", "MRKCFC", "MRKCFC-R2")
if (!all(VARIANTS %in% names(C$labels))) stop("rd_combined.rds$labels lacks a variant column")

dat  <- rd_load(window = as.integer(C$window))
prep <- rd_prepare(dat, log = FALSE, verbose = TRUE, time = "warped")
if (!identical(as.integer(prep$ids), as.integer(C$ids)) ||
    !identical(as.integer(C$labels$ID), as.integer(prep$ids)))
  stop("the prepared cohort does not match rd_combined.rds (ids differ)")
cfg <- load_config(); cfg$global$nK_fit <- as.integer(C$nK)
g <- cfg$global
cat(sprintf("cohort n=%d | nK=%d | common reconstruction: recon_robust=%s recon_scatter=%s\n",
            prep$n, as.integer(C$nK), as.character(g$recon_robust), as.character(g$recon_scatter)))

# one subject id and one marker index per observation, in the order
# reconstruct_given_partition() unlists them (marker-major, then subject)
subj_of <- function(Ly) unlist(lapply(Ly, function(z) rep(seq_along(z), lengths(z))), use.names = FALSE)

ncore <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
if (is.na(ncore) || ncore < 1L) ncore <- 1L
QS <- c(0.005, 0.01, 0.02, 0.05)

one <- function(v) {
  labels <- as.integer(C$labels[[v]])
  Ly <- prep$Ly; Lt <- prep$Lt; keep <- seq_len(prep$n)
  if (smoke) {
    # both classes represented: the whole smaller class plus the first 100 others
    small <- as.integer(names(which.min(table(labels))))
    keep  <- sort(c(which(labels == small), head(which(labels != small), 100L)))
    Ly <- lapply(Ly, function(z) z[keep]); Lt <- lapply(Lt, function(z) z[keep])
    labels <- labels[keep]
  }
  t0 <- Sys.time()
  r  <- reconstruct_given_partition(Ly, Lt, labels, cfg)
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  e2 <- (r$y - r$yhat)^2
  N  <- length(e2); sse <- sum(e2)
  subj  <- subj_of(Ly)
  cls   <- labels[subj]
  small <- as.integer(names(which.min(table(labels))))
  best  <- C$best$objective[C$best$variant == v & C$best$K == C$K_star]
  best  <- if (length(best) == 1L) best else NA_real_
  row <- data.frame(variant = v, n_subjects = length(labels), n_obs = N,
                    sizes = paste(table(labels), collapse = "/"),
                    sse = sse, objective_combined = if (smoke) NA_real_ else best,
                    mse_all = sse / N, minutes = el, stringsAsFactors = FALSE)
  for (q in QS) {
    m   <- ceiling(q * N)
    top <- rank(-e2, ties.method = "first") <= m
    tag <- sub("^0\\.", "", formatC(q, format = "f", digits = 3))
    row[[paste0("share_top", tag)]]    <- sum(e2[top]) / sse
    row[[paste0("mse_trim", tag)]]     <- mean(e2[!top])
    row[[paste0("n_subj_top", tag)]]   <- length(unique(subj[top]))
    row[[paste0("n_small_top", tag)]]  <- length(unique(subj[top & cls == small]))
  }
  for (p in seq_along(Ly))
    row[[paste0("share_marker", p)]] <- sum(e2[r$variable == p]) / sse
  saveRDS(list(variant = v, smoke = smoke, ids = prep$ids[keep], labels = labels,
               subject = subj, variable = r$variable, class = cls,
               y = r$y, yhat = r$yhat, when = Sys.time()),
          file.path(OUT, sprintf("rd_recon_trim_%s%s.rds", gsub("[^A-Za-z0-9]", "", v), if (smoke) "_smoke" else "")))
  cat(sprintf("%-10s sse %.1f (combined objective %s) | top 1%%: %.1f%% of sse from %d obs of %d subjects, %d in the smaller class | mse all %.3f, trimmed %.3f | %.1f min\n",
              v, sse, if (is.na(row$objective_combined)) "NA" else sprintf("%.1f", row$objective_combined),
              100 * row$share_top010, ceiling(0.01 * N), row$n_subj_top010, row$n_small_top010,
              row$mse_all, row$mse_trim010, el))
  row
}

rows <- if (ncore > 1L) parallel::mclapply(VARIANTS, one, mc.cores = min(3L, ncore)) else lapply(VARIANTS, one)
bad <- vapply(rows, function(x) inherits(x, "try-error") || !is.data.frame(x), logical(1))
if (any(bad)) stop("reconstruction failed for: ", paste(VARIANTS[bad], collapse = ", "))
res <- do.call(rbind, rows)
res$sse_matches_objective <- ifelse(is.na(res$objective_combined), NA,
                                    abs(res$sse - res$objective_combined) < 1e-6 * pmax(1, res$objective_combined))
f <- file.path(OUT, if (smoke) "rd_recon_trim_smoke.csv" else "rd_recon_trim.csv")
write.csv(res, f, row.names = FALSE)
cat("\nwrote", f, "\n")
print(res[, c("variant", "sizes", "n_obs", "sse", "objective_combined", "sse_matches_objective",
              "share_top010", "mse_all", "mse_trim010", "share_top005", "share_top020", "share_top050")],
      row.names = FALSE, digits = 4)
if (!smoke && !all(res$sse_matches_objective %in% TRUE))
  cat("\nWARNING: a total does not match the objective realdata_combine.R selected on; see the sse and objective_combined columns.\n")
