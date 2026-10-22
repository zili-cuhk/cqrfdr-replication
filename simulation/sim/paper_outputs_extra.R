## Usage: Rscript sim/paper_outputs_extra.R <extra summary.csv> <main summary.csv> [outdir]
## Figure and table for the additional studies (CQR_SIM_SET=extra):
##   fig_signal.pdf    FDR and power against the signal strength |beta_j| (A1)
##   table_sparsity.tex FDR and power for s = 5, 10, 20, 30 active covariates (B1)
## The baseline of the main study (A_p1000_c30_gauss: |beta_j| = 0.8, s = 10)
## is added as the last point / the s = 10 row.  Style as in paper_outputs.R.

args <- commandArgs(TRUE)
ex <- read.csv(args[1], stringsAsFactors = FALSE)
mn <- read.csv(args[2], stringsAsFactors = FALSE)
outdir <- if (length(args) >= 3) args[3] else "results/extra/paper"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
base <- mn[mn$scenario == "A_p1000_c30_gauss", ]
base$scenario <- "BASE"
s <- rbind(ex[, intersect(names(ex), names(base))], base[, intersect(names(ex), names(base))])
get <- function(sc, m, what) { v <- s[s$scenario == sc & s$method == m, what]; if (length(v)) v[1] else NA }

spec <- list(
  list(m = "scqr_scad",    lab = "SCQR-SCAD",              col = "#0072B2", pch = 6,  lty = "44",   lwd = 1.2),
  list(m = "scqr_lasso",   lab = "SCQR-Lasso",             col = "#CC79A7", pch = 5,  lty = "1343", lwd = 1.2),
  list(m = "lhdcqr_z18",   lab = "L-HDCQR",                col = "#009E73", pch = 0,  lty = "13",   lwd = 1.3),
  list(m = "mds_dai_q0.2",  lab = "MDS (q = 0.2)",          col = "#E69F00", pch = 2,  lty = "4212", lwd = 1.3),
  list(m = "ds_q0.2",      lab = "DS (q = 0.2)",           col = "#D55E00", pch = 4,  lty = "2262", lwd = 1.3),
  list(m = "ssms_q0.2",    lab = "SS-MS, CV screen (q = 0.2)", col = "grey50", pch = 15, lty = "22", lwd = 1.3),
  list(m = "ssms_path_q0.1", lab = "SS-MS (q = 0.1)",      col = "black",   pch = 1,  lty = "32",   lwd = 1.4),
  list(m = "ssms_path_q0.2", lab = "SS-MS (q = 0.2)",      col = "black",   pch = 16, lty = "solid", lwd = 1.6))

## ---- A1 figure: x = |beta_j| ----
bs <- sort(unique(ex$beta_level[grepl("^S_beta", ex$scenario)]))
rows <- c(sprintf("S_beta%03d", round(1000 * bs)), "BASE"); x <- c(bs, 0.8)
panel <- function(what, ylab) {
  plot.new(); plot.window(xlim = c(0, 0.85), ylim = c(-0.02, 1.02), xaxs = "i", yaxs = "i")
  abline(h = seq(0, 1, 0.2), col = "grey90", lwd = 0.6); abline(v = x, col = "grey94", lwd = 0.6)
  if (what == "FDR") { abline(h = 0.2, col = "grey35", lty = "55", lwd = 0.9); abline(h = 0.1, col = "grey35", lty = "13", lwd = 0.9) }
  for (sp in spec) {
    y <- vapply(rows, function(r) get(r, sp$m, what), 0)
    lines(x, y, col = sp$col, lty = sp$lty, lwd = sp$lwd)
    points(x, y, col = sp$col, pch = sp$pch, cex = 0.85, lwd = 1.1, bg = "white")
  }
  box(col = "grey20", lwd = 0.8)
  axis(1, at = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.65, 0.8), tcl = -0.25, cex.axis = 0.8, mgp = c(0, 0.3, 0))
  axis(2, at = seq(0, 1, 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"), las = 1, tcl = -0.25,
       cex.axis = 0.8, mgp = c(0, 0.4, 0))
  mtext(expression("Signal strength " * "|" * beta[j] * "|"), side = 1, line = 1.7, cex = 0.75)
  mtext(ylab, side = 2, line = 2.0, cex = 0.8)
}
STYLE <- Sys.getenv("CQR_FIG_STYLE", "sun")
if (STYLE == "sun") {
  ## Sun, McCaw and Lin (2025) style: lines only, colour and line type per
  ## method, ggplot hue palette (4th black, 5th blue), grey dashed nominal
  ## levels, bold panel letters, axis lines left and bottom only.
  gg_hue <- function(n) { h <- seq(15, 375, length = n + 1); hcl(h = h, l = 65, c = 100)[1:n] }
  sun_cols <- gg_hue(8); sun_cols[4] <- "black"; sun_cols[5] <- "blue"
  ## same methods, colours and line types as the main figures (paper_outputs.R)
  sspec <- list(list(m = "ssms_path_q0.2", lab = "SS-MS (q = 0.2)", col = sun_cols[1], lty = "solid"),
                list(m = "ssms_path_q0.1", lab = "SS-MS (q = 0.1)", col = sun_cols[1], lty = "22"),
                list(m = "msms_q0.2",      lab = "MS-MS (q = 0.2)", col = sun_cols[3], lty = "solid"),
                list(m = "msms_q0.1",      lab = "MS-MS (q = 0.1)", col = sun_cols[3], lty = "22"),
                list(m = "ds_q0.2",        lab = "DS (q = 0.2)",    col = sun_cols[4], lty = "dotted"),
                list(m = "mds_dai_q0.2",   lab = "MDS (q = 0.2)",   col = sun_cols[5], lty = "dotdash"),
                list(m = "lhdcqr_z18",     lab = "L-HDCQR",         col = sun_cols[6], lty = "longdash"),
                list(m = "scqr_lasso",     lab = "SCQR-Lasso",      col = sun_cols[7], lty = "twodash"),
                list(m = "scqr_scad",      lab = "SCQR-SCAD",       col = sun_cols[8], lty = "1343"),
                list(m = "ssms_q0.2",      lab = "SS-MS, CV screen (q = 0.2)", col = "grey55", lty = "13"))
  span <- function(what, ylab, letter) {
    plot.new(); plot.window(xlim = c(0.08, 0.82), ylim = c(0, 1))
    if (what == "FDR") abline(h = c(0.1, 0.2), col = "grey", lty = 2, lwd = 1.3)
    for (sp in rev(sspec)) lines(x, vapply(rows, function(r) get(r, sp$m, what), 0), col = sp$col, lty = sp$lty, lwd = 2.2)
    axis(1, at = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.65, 0.8), lwd = 0, lwd.ticks = 1, tcl = -0.3, cex.axis = 0.95, mgp = c(0, 0.45, 0))
    axis(2, at = seq(0, 1, 0.25), labels = c("0.00", "0.25", "0.50", "0.75", "1.00"), las = 1, lwd = 0,
         lwd.ticks = 1, tcl = -0.3, cex.axis = 0.95, mgp = c(0, 0.5, 0))
    u <- par("usr"); segments(u[1], u[3], u[2], u[3], lwd = 1.2, xpd = NA); segments(u[1], u[3], u[1], u[4], lwd = 1.2, xpd = NA)
    title(xlab = expression("Signal strength " * "|" * beta[j] * "|"), line = 2.1, cex.lab = 1.08)
    title(ylab = ylab, line = 2.6, cex.lab = 1.08)
    mtext(letter, side = 3, adj = -0.18, line = 0.4, font = 2, cex = 1.15)
  }
  pdf(file.path(outdir, "fig_signal.pdf"), width = 10, height = 5.1, family = "Helvetica", pointsize = 14)
  layout(rbind(1:2, c(3, 3)), heights = c(1, 0.34))
  par(mar = c(3.4, 4.0, 1.6, 0.8))
  span("FDR", "FDR", "A"); span("power", "Power", "B")
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", title = "Method", legend = sapply(sspec, `[[`, "lab"), col = sapply(sspec, `[[`, "col"),
         lty = sapply(sspec, `[[`, "lty"), lwd = 2.2, ncol = 4, bty = "n", cex = 0.85, seg.len = 3.2, x.intersp = 0.5)
  invisible(dev.off())
} else {
pdf(file.path(outdir, "fig_signal.pdf"), width = 7.4, height = 3.6, family = "Helvetica")
layout(rbind(1:2, c(3, 3)), heights = c(1, 0.3))
par(mar = c(3, 3.2, 0.8, 0.6))
panel("FDR", "FDR"); panel("power", "Power")
par(mar = c(0, 0, 0, 0)); plot.new()
o <- rev(seq_along(spec))
legend("center", legend = sapply(spec[o], `[[`, "lab"), col = sapply(spec[o], `[[`, "col"),
       pch = sapply(spec[o], `[[`, "pch"), lty = sapply(spec[o], `[[`, "lty"), lwd = sapply(spec[o], `[[`, "lwd"),
       ncol = 4, bty = "n", cex = 0.78, seg.len = 2.4, x.intersp = 0.5, pt.bg = "white")
invisible(dev.off())
}

## ---- B1 table: s = 5, 10, 20, 30 ----
srow <- c(S_s05 = "5", BASE = "10", S_s20 = "20", S_s30 = "30")
cols <- c(ssms_path_q0.1 = "\\shortstack{SS-MS\\\\($q=0.1$)}", msms_q0.1 = "\\rev{\\shortstack{MS-MS\\\\($q=0.1$)}}",
          ssms_path_q0.2 = "\\shortstack{SS-MS\\\\($q=0.2$)}", msms_q0.2 = "\\rev{\\shortstack{MS-MS\\\\($q=0.2$)}}",
          ds_q0.2 = "\\rev{\\shortstack{DS\\\\($q=0.2$)}}", mds_dai_q0.2 = "\\rev{\\shortstack{MDS\\\\($q=0.2$)}}",
          lhdcqr_z18 = "\\rev{L-HDCQR}", scqr_lasso = "SCQR-Lasso", scqr_scad = "SCQR-SCAD")
qn <- function(m) if (grepl("q0.1$", m)) 0.1 else 0.2
cell <- function(sc, m) {
  f <- get(sc, m, "FDR"); se <- get(sc, m, "FDR_mcse"); pw <- get(sc, m, "power")
  if (is.na(f)) return("-- & --")
  fs <- sprintf("%.3f", f)
  if (!m %in% c("lhdcqr", "lhdcqr_z18", "scqr_lasso", "scqr_scad") && f > qn(m) + 2 * se ||
      m %in% c("lhdcqr", "lhdcqr_z18", "scqr_lasso", "scqr_scad") && f > 0.2 + 2 * se) fs <- paste0("\\textbf{", fs, "}")
  paste(fs, sprintf("%.3f", pw), sep = " & ")
}
tex <- c("\\begin{table}[!htbp]", "\\centering\\footnotesize", "\\setlength{\\tabcolsep}{2.2pt}",
  "\\caption{Empirical FDR and power by the number $s$ of active covariates ($n=800$, $p=1000$, $30\\%$ censoring, $|\\beta_j|=0.8$).}",
  "\\label{tab:sim-sparsity}", sprintf("\\resizebox{\\linewidth}{!}{\\begin{tabular}{l%s}", strrep("cc", length(cols))), "\\toprule",
  paste0(" & ", paste(sprintf("\\multicolumn{2}{c}{%s}", cols), collapse = " & "), " \\\\"),
  paste0(paste(sprintf("\\cmidrule(lr){%d-%d}", 2 * seq_along(cols), 2 * seq_along(cols) + 1), collapse = ""), ""),
  paste0("$s$ & ", paste(rep("FDR & Power", length(cols)), collapse = " & "), " \\\\"), "\\midrule",
  vapply(names(srow), function(sc) paste0(srow[[sc]], " & ", paste(vapply(names(cols), function(m) cell(sc, m), ""), collapse = " & "), " \\\\"), ""),
  "\\bottomrule", "\\end{tabular}}",
  "\\par\\smallskip\\parbox{\\linewidth}{\\scriptsize Entries are averages over 100 replications; $s=10$ is the baseline of Table~\\ref{tab:sim-factorial}. Bold: FDR above the nominal level of the column (0.2 for L-HDCQR, SCQR-Lasso and SCQR-SCAD) by more than two Monte Carlo standard errors. \\rev{MS-MS, DS,} MDS, L-HDCQR, SCQR-Lasso and SCQR-SCAD as in Table~\\ref{tab:sim-factorial}.}",
  "\\end{table}")
writeLines(tex, file.path(outdir, "table_sparsity.tex"))
cat("wrote", file.path(outdir, c("fig_signal.pdf", "table_sparsity.tex")), sep = "\n")

## ---- screening ablation (supplement): CV screen vs path-ordered screen ----
both <- rbind(ex[, intersect(names(ex), names(mn))], mn[, intersect(names(ex), names(mn))])
gb <- function(sc, m, what) { v <- both[both$scenario == sc & both$method == m, what]; if (length(v)) v[1] else NA }
abl <- c(setNames(sprintf("$|\\beta_j|=%s$", formatC(bs, format = "fg")), sprintf("S_beta%03d", round(1000 * bs))),
         A_p1000_c30_gauss = "$|\\beta_j|=0.8$ (baseline)",
         D_hetero_tau025 = "Low SNR, $\\tau=0.25$", D_hetero_tau050 = "Low SNR, $\\tau=0.50$",
         D_hetero_tau060 = "Low SNR, $\\tau=0.60$", Dm_hetero_tau075 = "Heteroscedastic, $\\tau=0.75$")
f3 <- function(v) ifelse(is.na(v), "--", sprintf("%.3f", v))
f1 <- function(v) ifelse(is.na(v), "--", sprintf("%.1f", v))
rowtex <- vapply(names(abl), function(sc) paste(abl[[sc]],
  f3(gb(sc, "ssms_q0.2", "screen_recall")), f1(gb(sc, "ssms_q0.2", "screen_size")),
  f3(gb(sc, "ssms_q0.2", "FDR")), f3(gb(sc, "ssms_q0.2", "power")),
  f3(gb(sc, "ssms_path_q0.2", "screen_recall")), f1(gb(sc, "ssms_path_q0.2", "screen_size")),
  f3(gb(sc, "ssms_path_q0.2", "FDR")), f3(gb(sc, "ssms_path_q0.2", "power")), sep = " & "), "")
tex <- c("\\begin{table}[!htbp]", "\\centering\\footnotesize", "\\setlength{\\tabcolsep}{3pt}",
  "\\caption{First-stage screen of SS-MS: cross-validated penalty only (CV) versus the path-ordered screen of Algorithm~\\ref{alg:ssms} ($q=0.2$, $n=800$, $p=1000$).}",
  "\\label{tab:sim-screen}", "\\begin{tabular}{lcccccccc}", "\\toprule",
  " & \\multicolumn{4}{c}{CV screen} & \\multicolumn{4}{c}{Path-ordered screen} \\\\",
  "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}",
  "Setting & Recall & Size & FDR & Power & Recall & Size & FDR & Power \\\\", "\\midrule",
  paste0(rowtex, " \\\\"), "\\bottomrule", "\\end{tabular}",
  "\\par\\smallskip\\parbox{\\linewidth}{\\scriptsize Recall: share of active covariates in the candidate set; Size: candidate-set size (at most $\\lfloor d_2/\\log d_2\\rfloor$). Averages over 100 replications.}",
  "\\end{table}")
writeLines(tex, file.path(outdir, "table_screen.tex"))
cat("wrote", file.path(outdir, "table_screen.tex"), "\n")
