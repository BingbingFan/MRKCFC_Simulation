# Figure 1 of the manuscript: class-wise trajectories of CEA, CA19-9 and CA242
# for the K = 2 partitions of MRKCFC-R2, MRKCFC and MKCFC, primary analysis.
#
# One penalised-spline smooth per class of ALL observed values of the class
# (ggplot2::geom_smooth with mgcv's gam, shrinkage cubic regression spline basis (bs = "cs") of
# dimension 8), no points, on a log axis; because the axis is logarithmic the
# smooth is fitted to log10 of the value, so every curve is a GEOMETRIC-MEAN
# trajectory with its 95 per cent confidence band. One row per method, each
# with its own legend at the right, stacked with patchwork so the panel
# columns align. Time is months after surgery.
#
# Chosen on 2026-09-13 over the median-and-quantile version
# (make_rd_tables.R, traj_fig(), kept as fig_rd_traj_quantiles).
#
# STANDALONE, from study/:
#
#   Rscript R/plot_fig1_traj.R <rd-dir> <out-stem>
#   Rscript R/plot_fig1_traj.R ../HPC/V3/rd ../doc/manu/figures/fig_rd_traj
#
# writes <out-stem>.pdf and <out-stem>.png. <rd-dir> must hold
# rd_agreement.rds (the reported partitions; the agreement file rather than
# the combine because it carries the same labels under the same names) and
# the cohort file is found through realdata_load.R.
#
# AS A LIBRARY: make_rd_tables.R sources this file and calls fig1_gam() with
# the data it has already prepared, so the manuscript figure and this script
# cannot drift apart.

suppressWarnings(suppressMessages({ library(ggplot2); library(patchwork) }))

# The class names, derived from what distinguishes the two classes of each
# method (2026-09-13): the smaller class is named, the larger is
# "Low", or "Stable" for the proposed method whose classes share the same
# median and differ in excursions.
FIG1_NAME_SMALL <- c(`MRKCFC-R2` = "Episodic", MRKCFC = "Elevated", MKCFC = "Extreme",
                     mclust = "Elevated", `k-means` = "Extreme", `k-means/rob` = "Extreme",
                     teigen = "Extreme")
FIG1_NAME_LARGE <- c(`MRKCFC-R2` = "Stable")
# Display names (2026-09-13): -S, spatial-sign second stage (internally
# "MRKCFC-R2"); -M, second-moment second stage (internally "MRKCFC"). The same
# map as DISP in make_rd_tables.R; internal identifiers are never renamed.
FIG1_DISP <- c(`MRKCFC-R2` = "MRKCFC-S", MRKCFC = "MRKCFC-M", MKCFC = "MKCFC")

derive_class_names <- function(labels) {
  out <- lapply(names(labels), function(v) {
    g <- as.integer(labels[[v]]); cl <- sort(unique(g)); sz <- table(g)[as.character(cl)]
    if (length(cl) != 2L) return(setNames(paste("Class", cl), cl))
    small <- cl[which.min(sz)]
    setNames(ifelse(cl == small, FIG1_NAME_SMALL[[v]],
                    if (v %in% names(FIG1_NAME_LARGE)) FIG1_NAME_LARGE[[v]] else "Low"), cl)
  })
  names(out) <- names(labels)
  out
}

# Long data frame of every observed value with its method, marker and class.
#   markers  the three long frames from rd_load()$markers (ID, time in months, yit)
#   ids      the fitted subjects, in the order of the label vectors
#   labels   named list of integer label vectors, one per method
fig1_data <- function(markers, ids, labels, class_names, methods, disp = FIG1_DISP) {
  rows <- list()
  for (v in methods) {
    g <- as.integer(labels[[v]]); cl <- sort(unique(g)); small_cl <- cl[which.min(table(g))]
    for (m in names(markers)) {
      x <- markers[[m]]; x <- x[x$ID %in% ids, ]
      k <- g[match(x$ID, ids)]
      rows[[length(rows) + 1L]] <- data.frame(
        method = v, marker = sub("CA199", "CA19-9", m), months = x$time, value = x$yit,
        class = sprintf("%s (n = %d)", class_names[[v]][as.character(k)], as.integer(table(g)[as.character(k)])),
        small = k == small_cl, stringsAsFactors = FALSE)
    }
  }
  D <- do.call(rbind, rows)
  D$marker <- factor(D$marker, levels = c("CEA", "CA19-9", "CA242"),
                     labels = c("CEA (ng/mL)", "CA19-9 (U/mL)", "CA242 (U/mL)"))
  D
}

# The figure. Returns the patchwork object; writes it if `file` is given.
fig1_gam <- function(D, methods, disp = FIG1_DISP, file = NULL,
                     width = 7.6, height = 6.4, k = 8) {
  plots <- lapply(seq_along(methods), function(i) {
    v <- methods[i]; Dv <- D[D$method == v, ]
    lev <- unique(Dv$class[order(!Dv$small)])           # smaller class first -> red
    Dv$class <- factor(Dv$class, levels = lev)
    cols <- setNames(c("#d62728", "#1f77b4")[seq_along(lev)], lev)
    ggplot(Dv, aes(months, value, colour = class, fill = class)) +
      geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = k), se = TRUE,
                  linewidth = 0.9, alpha = 0.22) +
      facet_wrap(~ marker, nrow = 1, scales = "free_y") +
      scale_y_log10() + scale_x_continuous(breaks = seq(0, 36, 6)) +
      scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
      labs(x = if (i == length(methods)) "Months after surgery" else NULL, y = NULL,
           colour = disp[[v]], fill = disp[[v]]) +
      theme_bw(base_size = 9) +
      theme(legend.position = "right", legend.title = element_text(face = "bold", size = 8.5),
            legend.text = element_text(size = 8), strip.background = element_blank(),
            panel.grid.minor = element_blank(),
            strip.text.x = if (i == 1) element_text(size = 8.5, face = "bold") else element_blank())
  })
  p <- Reduce(`/`, plots) + plot_layout(guides = "keep")
  if (!is.null(file)) {
    for (ext in c("pdf", "png"))
      ggsave(paste0(file, ".", ext), p, width = width, height = height, units = "in", dpi = 300)
  }
  invisible(p)
}

# ------------------------------------------------------------ standalone ---
if (sys.nframe() == 0L) {
  suppressWarnings(suppressMessages(source(file.path("R", "realdata_load.R"))))
  args <- commandArgs(trailingOnly = TRUE)
  rd   <- if (length(args) >= 1) args[1] else file.path("..", "HPC", "V3", "rd")
  stem <- if (length(args) >= 2) args[2] else file.path("..", "doc", "manu", "figures", "fig_rd_traj")
  dat  <- rd_load(window = 36L)
  prep <- rd_prepare(dat, log = FALSE, verbose = FALSE, time = "warped")
  A    <- readRDS(file.path(rd, "rd_agreement.rds"))
  stopifnot(A$n == prep$n)
  methods <- c("MRKCFC-R2", "MRKCFC", "MKCFC")
  cn <- derive_class_names(A$labels[methods])
  D  <- fig1_data(dat$markers, prep$ids, A$labels, cn, methods)
  fig1_gam(D, methods, file = stem)
  cat(sprintf("wrote %s.pdf and %s.png (%d observed values)\n", stem, stem, nrow(D)))
}
