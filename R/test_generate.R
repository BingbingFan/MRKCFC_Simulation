# Checks the generator against the design it is supposed to implement.
# Run from study/:   Rscript R/test_generate.R

suppressMessages(library(MRKCFC))
source(file.path("R", "generate.R"))
# gen_scores and friends are internal to MRKCFC; expose them for the generator
for (f in c("gen_scores", "vm_intermediate", "rmvn", "FLEISHMAN", "IQR_SCALE"))
  assign(f, get(f, envir = asNamespace("MRKCFC")), envir = globalenv())

cfg <- load_config()
fail <- 0L
ok <- function(lab, cond, extra = "") {
  good <- isTRUE(cond); if (!good) fail <<- fail + 1L
  cat(sprintf("%-56s %s %s\n", lab, if (good) "PASS" else "**FAIL**", extra))
}

sc <- cfg$scenarios[cfg$scenarios$scenario_id == "A01", ]
d1 <- generate_replicate(cfg, sc, rep = 1L)

cat("=== shape ===\n")
n <- as.integer(cfg$global$n); P <- as.integer(cfg$global$P)
ok("P variables", length(d1$Ly) == P && length(d1$Lt) == P)
ok("n subjects in every variable", all(vapply(d1$Ly, length, 1L) == n))
ok("labels are 1..K", all(sort(unique(d1$labels)) == c(1L, 2L)))
ok("scores are n x (P * M_true)",
   all(dim(d1$scores) == c(n, P * as.integer(cfg$global$M_true))))
ok("times match values elementwise",
   all(vapply(seq_len(P), function(p)
     all(lengths(d1$Ly[[p]]) == lengths(d1$Lt[[p]])), TRUE)))
ok("one shared schedule per subject across variables",
   identical(d1$Lt[[1]][[7]], d1$Lt[[2]][[7]]) &&
   identical(d1$Lt[[1]][[7]], d1$Lt[[3]][[7]]))

cat("\n=== observation counts ===\n")
obs_min <- as.integer(cfg$global$obs_min)
ok(sprintf("no subject below obs_min = %d", obs_min), min(d1$n_i) >= obs_min)
cat(sprintf("   realised range %d to %d, mean %.2f (manuscript reports 4 to 12)\n",
            min(d1$n_i), max(d1$n_i), mean(d1$n_i)))
ok("times lie on the 101-point grid",
   all(unlist(d1$Lt[[1]]) %in% seq(0, 1, length.out = as.integer(cfg$global$grid_n))))
ok("times are sorted within subject",
   all(vapply(d1$Lt[[1]], function(x) !is.unsorted(x), TRUE)))

cat("\n=== reproducibility ===\n")
d1b <- generate_replicate(cfg, sc, rep = 1L)
ok("same cell reproduces exactly", identical(d1, d1b))
d2 <- generate_replicate(cfg, sc, rep = 2L)
ok("a different replicate differs", !identical(d1$labels, d2$labels))
sc2 <- cfg$scenarios[cfg$scenarios$scenario_id == "A02", ]
ok("a different scenario differs",
   !identical(generate_replicate(cfg, sc2, rep = 1L)$Ly[[1]][[1]], d1$Ly[[1]][[1]]))

cat("\n=== normalization (Equation 13) ===\n")
if (isTRUE(as.logical(cfg$global$normalize))) {
  for (p in seq_len(P)) {
    pooled <- unlist(d1$Ly[[p]], use.names = FALSE)
    ok(sprintf("variable %d centred at 0 with unit IQR", p),
       abs(stats::median(pooled)) < 1e-8 && abs(stats::IQR(pooled) - 1) < 1e-8,
       sprintf("median %+.2e IQR %.6f", stats::median(pooled), stats::IQR(pooled)))
  }
  ok("normalization constants were retained for back-transform",
     length(d1$normalization) == P && all(vapply(d1$normalization,
       function(z) is.finite(z$centre) && z$scale > 0, TRUE)))
}

cat("\n=== separation acts on the mean, as d should ===\n")
sepd <- function(id) {
  s <- cfg$scenarios[cfg$scenarios$scenario_id == id, ]
  g <- generate_replicate(cfg, s, rep = 1L)
  # raw scale: undo normalization so the comparison is in original units
  nb <- g$normalization
  raw <- lapply(seq_len(P), function(p)
    vapply(seq_len(n), function(i)
      mean(g$Ly[[p]][[i]]) * (if (is.null(nb)) 1 else nb[[p]]$scale) +
        (if (is.null(nb)) 0 else nb[[p]]$centre), 0))
  vapply(raw, function(v) abs(mean(v[g$labels == 1]) - mean(v[g$labels == 2])), 0)
}
g5 <- sepd("A01"); g8 <- sepd("A02")
cat(sprintf("   d=5 per-variable class mean gap: %s\n", paste(round(g5, 3), collapse = ", ")))
cat(sprintf("   d=8 per-variable class mean gap: %s\n", paste(round(g8, 3), collapse = ", ")))
ok("d = 8 separates more than d = 5 on every variable", all(g8 > g5))

cat("\n=== arms differ in tail behaviour, IQR-matched ===\n")
iqr1 <- NULL
for (id in c("A01", "A03", "A05", "A07", "A09", "A11")) {
  s <- cfg$scenarios[cfg$scenarios$scenario_id == id, ]
  g <- generate_replicate(cfg, s, rep = 1L)
  q <- stats::IQR(g$scores[, 1]); if (is.null(iqr1)) iqr1 <- q
  cat(sprintf("   %s arm %d: score IQR ratio %.3f, excess kurtosis %8.2f\n",
              id, s$arm, q / iqr1,
              mean((g$scores[,1] - mean(g$scores[,1]))^4) /
                stats::sd(g$scores[,1])^4 - 3))
}

cat("\n")
if (fail > 0L) stop(sprintf("generator tests: %d FAILED", fail), call. = FALSE)
cat("generator tests: PASS\n")
