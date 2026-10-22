## Usage: Rscript sim/paper_outputs.R [summary.csv] [outdir]
## Builds the LaTeX tables and PDF figures of Section 4 directly from
## results/main/summary.csv, so every number in the paper is traceable.

args <- commandArgs(TRUE)
infile <- if (length(args) >= 1) args[1] else "results/main/summary.csv"
outdir <- if (length(args) >= 2) args[2] else "results/main/paper"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
s <- read.csv(infile, stringsAsFactors = FALSE)

## methods shown: SS-MS and MS-MS at q = 0.1 and 0.2, DS and MDS at q = 0.2.
## SS-MS = Algorithm 1 with the path-ordered screen (method "ssms_path").
M <- c(ssms01 = "ssms_path_q0.1", msms01 = "msms_q0.1", ssms02 = "ssms_path_q0.2", msms02 = "msms_q0.2",
       ds = "ds_q0.2", mds = "mds_dai_q0.2", lhd = "lhdcqr_z18", las = "scqr_lasso", scad = "scqr_scad")
lab <- c(ssms01 = "\\shortstack{SS-MS\\\\($q=0.1$)}", msms01 = "\\rev{\\shortstack{MS-MS\\\\($q=0.1$)}}",
         ssms02 = "\\shortstack{SS-MS\\\\($q=0.2$)}", msms02 = "\\rev{\\shortstack{MS-MS\\\\($q=0.2$)}}",
         ds = "\\rev{\\shortstack{DS\\\\($q=0.2$)}}", mds = "\\rev{\\shortstack{MDS\\\\($q=0.2$)}}",
         lhd = "\\rev{L-HDCQR}", las = "SCQR-Lasso", scad = "SCQR-SCAD")
## \rev{} (blue) marks the columns added or recomputed in the revision
get <- function(sc, m, what) { v <- s[s$scenario == sc & s$method == m, what]; if (length(v)) v else NA }

## an FDR entry is flagged when it exceeds the nominal level by more than 2 MCSE
## (nominal 0.1 for SS-MS q = 0.1, 0.2 for every other column)
cell <- function(sc, key) {
  m <- M[[key]]; f <- get(sc, m, "FDR"); se <- get(sc, m, "FDR_mcse"); pw <- get(sc, m, "power")
  q <- if (key %in% c("ssms01", "msms01")) 0.1 else 0.2
  fs <- sprintf("%.3f", f); if (!is.na(f) && f > q + 2 * se) fs <- paste0("\\textbf{", fs, "}")
  paste(fs, sprintf("%.3f", pw), sep = " & ")
}
table_tex <- function(rows, labels, caption, label, group = NULL) {
  keys <- names(M)
  head <- c("\\begin{table}[!htbp]", "\\centering\\footnotesize", "\\setlength{\\tabcolsep}{2.2pt}",
            sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
            sprintf("\\resizebox{\\linewidth}{!}{\\begin{tabular}{l%s}", strrep("cc", length(keys))), "\\toprule",
            paste0(" & ", paste(sprintf("\\multicolumn{2}{c}{%s}", lab[keys]), collapse = " & "), " \\\\"),
            paste(sprintf("\\cmidrule(lr){%d-%d}", 2 * seq_along(keys), 2 * seq_along(keys) + 1), collapse = ""),
            paste0("Setting & ", paste(rep("FDR & Power", length(keys)), collapse = " & "), " \\\\"), "\\midrule")
  body <- character()
  for (i in seq_along(rows)) {
    if (!is.null(group) && !is.na(group[i])) body <- c(body, sprintf("\\multicolumn{%d}{l}{\\emph{%s}} \\\\", 1 + 2 * length(keys), group[i]))
    body <- c(body, paste0(labels[i], " & ", paste(vapply(keys, function(k) cell(rows[i], k), ""), collapse = " & "), " \\\\"))
  }
  mcse <- s[s$scenario %in% rows & s$method %in% M, ]
  foot <- c("\\bottomrule", "\\end{tabular}}",
            sprintf("\\par\\smallskip\\parbox{\\linewidth}{\\scriptsize Entries are averages over 100 replications. The largest Monte Carlo standard error is %.3f for FDR and %.3f for power. The column label gives the nominal level of SS-MS, \\rev{MS-MS, DS} and MDS; L-HDCQR, SCQR-Lasso and SCQR-SCAD do not target an FDR level. Bold: FDR above the nominal level of the column by more than two Monte Carlo standard errors (0.2 for L-HDCQR, SCQR-Lasso and SCQR-SCAD).}",
                    max(mcse$FDR_mcse, na.rm = TRUE), max(mcse$power_mcse, na.rm = TRUE)),
            "\\end{table}")
  c(head, body, foot)
}

## Table 1: factorial block
fa <- expand.grid(x = c("gauss", "t5"), c = c(10, 30, 50), p = c(100, 1000), stringsAsFactors = FALSE)
fa <- fa[order(fa$p, fa$x, fa$c), ]
rows1 <- sprintf("A_p%d_c%02d_%s", fa$p, fa$c, fa$x)
labs1 <- sprintf("%d\\%%", fa$c)
grp1 <- rep(NA, nrow(fa)); first <- !duplicated(paste(fa$p, fa$x))
grp1[first] <- sprintf("$p=%d$, %s covariates", fa$p[first], ifelse(fa$x[first] == "t5", "multivariate $t_5$", "Gaussian"))
t1 <- table_tex(rows1, labs1, "Empirical FDR and power by dimension, censoring rate and covariate distribution ($n=800$, $\\tau=0.5$, $\\rho=0.5$, normal errors).", "tab:sim-factorial", grp1)

## Table 2: one-factor block around the baseline (p = 1000, 30% censoring)
rows2 <- c("B_tau025", "B_tau040", "A_p1000_c30_gauss", "B_tau060", "B_tau070",
           "B_rho01", "B_rho03", "A_p1000_c30_gauss", "B_rho07",
           "A_p1000_c30_gauss", "B_err_t5", "B_err_exp")
labs2 <- c("$\\tau=0.25$", "$\\tau=0.40$", "$\\tau=0.50$ (baseline)", "$\\tau=0.60$", "$\\tau=0.70$",
           "$\\rho=0.1$", "$\\rho=0.3$", "$\\rho=0.5$ (baseline)", "$\\rho=0.7$",
           "Normal (baseline)", "$t_5$", "Exponential")
grp2 <- c("Target quantile", NA, NA, NA, NA, "Covariate correlation", NA, NA, NA, "Error distribution", NA, NA)
t2 <- table_tex(rows2, labs2, "Empirical FDR and power when one factor is changed from the baseline ($n=800$, $p=1000$, 30\\% censoring, Gaussian covariates).", "tab:sim-onefactor", grp2)

## Table 3: heteroscedastic model
rows3 <- c("Dm_hetero_tau025", "Dm_hetero_tau050", "Dm_hetero_tau075")
labs3 <- c("$\\tau=0.25$", "$\\tau=0.50$", "$\\tau=0.75$")
t3 <- table_tex(rows3, labs3, "Empirical FDR and power under the heteroscedastic model, in which the eleventh covariate is inactive at the median and active at $\\tau=0.25$ and $0.75$ ($n=800$, $p=1000$, 30\\% censoring).", "tab:sim-hetero")

## Supplementary table: low signal-to-noise heteroscedastic design
rows4 <- c("D_hetero_tau025", "D_hetero_tau050", "D_hetero_tau060")
labs4 <- c("$\\tau=0.25$", "$\\tau=0.50$", "$\\tau=0.60$")
t4 <- table_tex(rows4, labs4, "Heteroscedastic model with a larger noise scale ($\\sigma=0.5$, $\\gamma=3$; average noise scale 2.5 times that of the location model).", "tab:sim-hetero-lowsnr")

writeLines(c(t1, "", t2, "", t3), file.path(outdir, "tables_main.tex"))
writeLines(t4, file.path(outdir, "table_supp_lowsnr.tex"))

## ---------------------------------------------------------------- figures ----
## Style of the simulation figures in Candes et al. (2018, JRSSB), Barber and
## Candes (2015, AoS) and Dai et al. (2023, JASA): a grid of framed panels with
## grey strip labels (columns = settings, rows = FDR / power), the true x values
## on the horizontal axis, one line + symbol per method, the nominal FDR levels
## as dashed / dotted reference lines and one legend row below the grid.
## Competitors use open symbols, so coincident values (e.g. power 1) remain
## visible as superimposed outlines; SS-MS is black and drawn last.
## Okabe-Ito colours; colour, symbol and line type all differ between methods.
fig_spec <- list(
  list(m = "scqr_scad",  lab = "SCQR-SCAD",       col = "#0072B2", pch = 6,  lty = "44",   lwd = 1.2, cex = 1.0),
  list(m = "scqr_lasso", lab = "SCQR-Lasso",      col = "#CC79A7", pch = 5,  lty = "1343", lwd = 1.2, cex = 1.0),
  list(m = "lhdcqr_z18", lab = "L-HDCQR",         col = "#009E73", pch = 0,  lty = "13",   lwd = 1.3, cex = 1.0),
  list(m = "mds_dai_q0.2",   lab = "MDS (q = 0.2)",   col = "#E69F00", pch = 2,  lty = "4212", lwd = 1.3, cex = 1.0),
  list(m = "ds_q0.2",       lab = "DS (q = 0.2)",    col = "#D55E00", pch = 4,  lty = "2262", lwd = 1.3, cex = 1.0),
  list(m = "msms_q0.2",     lab = "MS-MS (q = 0.2)", col = "grey45",  pch = 17, lty = "62",   lwd = 1.5, cex = 0.8),
  list(m = "ssms_path_q0.1",  lab = "SS-MS (q = 0.1)", col = "black",   pch = 1,  lty = "32",   lwd = 1.4, cex = 0.9),
  list(m = "ssms_path_q0.2",  lab = "SS-MS (q = 0.2)", col = "black",   pch = 16, lty = "solid", lwd = 1.6, cex = 0.8))
leg_order <- c(8, 7, 6, 5, 4, 3, 2, 1) # legend: SS-MS first
FONT <- "Helvetica"

## one framed panel; strip = column label drawn above the frame (top row only)
fpanel <- function(x, rows, what, xlim, strip = NULL, xaxis = FALSE, yaxis = FALSE, ref = FALSE,
                   xlabels = x) {
  plot.new(); plot.window(xlim = xlim, ylim = c(-0.04, 1.04), xaxs = "i", yaxs = "i")
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col = "white", border = NA)
  abline(h = seq(0, 1, 0.2), col = "grey90", lwd = 0.6)
  abline(v = x, col = "grey94", lwd = 0.6)
  if (ref) {
    abline(h = 0.2, col = "grey35", lty = "55", lwd = 0.9)
    abline(h = 0.1, col = "grey35", lty = "13", lwd = 0.9)
  }
  for (sp in fig_spec) {
    y <- vapply(rows, function(r) get(r, sp$m, what), 0)
    lines(x, y, col = sp$col, lty = sp$lty, lwd = sp$lwd)
    points(x, y, col = sp$col, pch = sp$pch, cex = sp$cex, lwd = 1.1, bg = "white")
  }
  box(col = "grey20", lwd = 0.8)
  if (!is.null(strip)) {
    h <- 0.13 * diff(usr[3:4])
    rect(usr[1], usr[4], usr[2], usr[4] + h, col = "grey88", border = "grey20", lwd = 0.8, xpd = NA)
    text(mean(usr[1:2]), usr[4] + h / 2, strip, cex = 0.85, xpd = NA)
  }
  if (xaxis) axis(1, at = x, labels = xlabels, tcl = -0.25, cex.axis = 0.78, mgp = c(0, 0.3, 0), col = "grey20", lwd = 0.8)
  if (yaxis) axis(2, at = seq(0, 1, 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"),
                  las = 1, tcl = -0.25, cex.axis = 0.78, mgp = c(0, 0.4, 0), col = "grey20", lwd = 0.8)
}
## grid: rows FDR / power, one column per entry of `cols`
fgrid <- function(file, cols, width, height, xtitle_line = 1.6) {
  K <- length(cols)
  pdf(file, width = width, height = height, family = FONT)
  layout(rbind(seq_len(K), K + seq_len(K), rep(2 * K + 1, K)), heights = c(1, 1, 0.22))
  par(oma = c(0, 2.6, 1.6, 0.4), mar = c(0.5, 0.5, 0.5, 0.5))
  for (what in c("FDR", "power")) for (k in seq_len(K)) {
    cl <- cols[[k]]
    xl <- range(cl$x) + c(-1, 1) * 0.07 * diff(range(cl$x))
    if (what == "power") par(mar = c(2.8, 0.5, 0.9, 0.5)) else par(mar = c(0.9, 0.5, 0.5, 0.5))
    fpanel(cl$x, cl$rows, what, xl, strip = if (what == "FDR") cl$strip, xaxis = what == "power",
           yaxis = k == 1, ref = what == "FDR", xlabels = cl$xlabels)
    if (what == "power") mtext(cl$xtitle, side = 1, line = xtitle_line, cex = 0.72)
    if (k == 1) mtext(if (what == "FDR") "FDR" else "Power", side = 2, line = 2.0, cex = 0.8)
  }
  par(mar = c(0, 0, 0, 0)); plot.new()
  o <- leg_order
  legend("center", legend = sapply(fig_spec[o], `[[`, "lab"), col = sapply(fig_spec[o], `[[`, "col"),
         pch = sapply(fig_spec[o], `[[`, "pch"), lty = sapply(fig_spec[o], `[[`, "lty"),
         lwd = sapply(fig_spec[o], `[[`, "lwd"), pt.cex = 0.95, horiz = TRUE, bty = "n", cex = 0.86,
         seg.len = 2.6, x.intersp = 0.5, text.width = NA)
  invisible(dev.off())
}

## ---- alternative style: Sun, McCaw and Lin (2025, JASA; csmGmm reproduce code) ----
## ggplot2 + cowplot::theme_cowplot look, rebuilt in base R so no packages are
## needed: white background, no grid or frame, left and bottom axis lines only,
## one line per method (colour AND line type 1:6, no symbols, no error bars),
## grey dashed line at the nominal FDR, panels lettered A, B, ... in bold, the
## setting written into the y-axis label, one horizontal legend below.
## Colours: ggplot's default hue palette for six methods with the 4th set to
## black and the 5th to blue, exactly as in Sun's plot_fig1.R / plot_fig2.R.
## Default style; CQR_FIG_STYLE=journal gives the framed-panel style instead.
gg_hue <- function(n) { h <- seq(15, 375, length = n + 1); hcl(h = h, l = 65, c = 100)[1:n] }
sun_cols <- gg_hue(8); sun_cols[4] <- "black"; sun_cols[5] <- "blue"
## SS-MS and MS-MS: one colour per method, solid at q = 0.2 and dashed at q = 0.1
sun_spec <- list(   # legend / colour order: proposed methods first, as in Sun's figures
  list(m = "ssms_path_q0.2", lab = "SS-MS (q = 0.2)", col = sun_cols[1], lty = "solid"),
  list(m = "ssms_path_q0.1", lab = "SS-MS (q = 0.1)", col = sun_cols[1], lty = "22"),
  list(m = "msms_q0.2",      lab = "MS-MS (q = 0.2)", col = sun_cols[3], lty = "solid"),
  list(m = "msms_q0.1",      lab = "MS-MS (q = 0.1)", col = sun_cols[3], lty = "22"),
  list(m = "ds_q0.2",        lab = "DS (q = 0.2)",    col = sun_cols[4], lty = "dotted"),
  list(m = "mds_dai_q0.2",   lab = "MDS (q = 0.2)",   col = sun_cols[5], lty = "dotdash"),
  list(m = "lhdcqr_z18",     lab = "L-HDCQR",         col = sun_cols[6], lty = "longdash"),
  list(m = "scqr_lasso",     lab = "SCQR-Lasso",      col = sun_cols[7], lty = "twodash"),
  list(m = "scqr_scad",      lab = "SCQR-SCAD",       col = sun_cols[8], lty = "1343"))
sun_panel <- function(x, rows, what, xlab, ylab, letter, xlabels = x, ref = what == "FDR") {
  xl <- range(x) + c(-0.03, 0.03) * diff(range(x))
  plot.new(); plot.window(xlim = xl, ylim = c(0, 1))
  if (ref) abline(h = 0.1, col = "grey", lty = 2, lwd = 1.3)
  if (ref) abline(h = 0.2, col = "grey", lty = 2, lwd = 1.3)
  for (k in rev(seq_along(sun_spec))) {           # SS-MS drawn last, on top
    sp <- sun_spec[[k]]
    y <- vapply(rows, function(r) get(r, sp$m, what), 0)
    lines(x, y, col = sp$col, lty = sp$lty, lwd = 2.2)
  }
  axis(1, at = x, labels = xlabels, lwd = 0, lwd.ticks = 1, tcl = -0.3, cex.axis = 0.95, mgp = c(0, 0.45, 0))
  axis(2, at = seq(0, 1, 0.25), labels = c("0.00", "0.25", "0.50", "0.75", "1.00"), las = 1,
       lwd = 0, lwd.ticks = 1, tcl = -0.3, cex.axis = 0.95, mgp = c(0, 0.5, 0))
  u <- par("usr"); segments(u[1], u[3], u[2], u[3], lwd = 1.2, xpd = NA); segments(u[1], u[3], u[1], u[4], lwd = 1.2, xpd = NA)
  title(xlab = xlab, line = 2.0, cex.lab = 1.08); title(ylab = ylab, line = 2.6, cex.lab = 1.08)
  mtext(letter, side = 3, adj = -0.22, line = 0.4, font = 2, cex = 1.15)
}
sun_grid <- function(file, cols, width, height, ps = 12) {
  K <- length(cols)
  pdf(file, width = width, height = height, family = "Helvetica", pointsize = ps)
  layout(rbind(seq_len(K), K + seq_len(K), rep(2 * K + 1, K)), heights = c(1, 1, 0.3))
  par(mar = c(3.4, 4.0, 1.6, 0.8), oma = c(0, 0, 0, 0))
  L <- 0L
  for (what in c("FDR", "power")) for (cl in cols) {
    L <- L + 1L
    yl <- if (is.null(cl$ytag)) (if (what == "FDR") "FDR" else "Power")
          else bquote(.(if (what == "FDR") "FDR" else "Power") ~ "(" * .(cl$ytag) * ")")
    sun_panel(cl$x, cl$rows, what, cl$xtitle, yl, LETTERS[L], cl$xlabels)
  }
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", title = "Method", legend = sapply(sun_spec, `[[`, "lab"),
         col = sapply(sun_spec, `[[`, "col"), lty = sapply(sun_spec, `[[`, "lty"), lwd = 2.2,
         ncol = 5, bty = "n", cex = 0.95, seg.len = 3.2, x.intersp = 0.5)
  invisible(dev.off())
}
STYLE <- Sys.getenv("CQR_FIG_STYLE", "sun")

## Figure 1: factorial; columns = (p, covariates), x = censoring rate
cfg <- list(c(100, "gauss"), c(100, "t5"), c(1000, "gauss"), c(1000, "t5"))
cols1 <- lapply(cfg, function(cf) list(
  x = c(10, 30, 50), xlabels = c("10%", "30%", "50%"),
  rows = sprintf("A_p%s_c%02d_%s", cf[1], c(10, 30, 50), cf[2]),
  ytag = if (cf[2] == "gauss") bquote(italic(p) * " = " * .(cf[1]) * ", Gaussian")
         else bquote(italic(p) * " = " * .(cf[1]) * ", " * italic(t)[5]),
  strip = if (cf[2] == "gauss") bquote(italic(p) * " = " * .(cf[1]) * ", Gaussian")
          else bquote(italic(p) * " = " * .(cf[1]) * ", " * italic(t)[5]),
  xtitle = "Censoring rate"))
if (STYLE == "sun") {
  lapply(seq_along(cols1), function(k) cols1[[k]]$xtitle <<- "Censoring rate")
  sun_grid(file.path(outdir, "fig_factorial.pdf"), cols1, width = 12, height = 7.4, ps = 18)
} else fgrid(file.path(outdir, "fig_factorial.pdf"), cols1, width = 7.4, height = 4.3)

## Figure 2: target quantile and covariate correlation (baseline at 0.5)
cols2 <- list(
  list(x = c(0.25, 0.4, 0.5, 0.6, 0.7), xlabels = c("0.25", "0.4", "0.5", "0.6", "0.7"),
       rows = c("B_tau025", "B_tau040", "A_p1000_c30_gauss", "B_tau060", "B_tau070"),
       strip = expression("Target quantile " * tau), xtitle = expression("Target quantile" ~ tau), ytag = "location"),
  list(x = c(0.1, 0.3, 0.5, 0.7), xlabels = c("0.1", "0.3", "0.5", "0.7"),
       rows = c("B_rho01", "B_rho03", "A_p1000_c30_gauss", "B_rho07"),
       strip = expression("Covariate correlation " * rho), xtitle = expression("Correlation" ~ rho), ytag = "location"),
  list(x = c(0.25, 0.5, 0.75), xlabels = c("0.25", "0.5", "0.75"), rows = rows3,
       strip = expression("Heteroscedastic model, " * tau), xtitle = expression("Target quantile" ~ tau), ytag = "heteroscedastic"))
## one figure: location model (tau, rho) and heteroscedastic model (tau)
if (STYLE == "sun") sun_grid(file.path(outdir, "fig_onefactor.pdf"), cols2, width = 10, height = 7, ps = 15) else
fgrid(file.path(outdir, "fig_onefactor.pdf"), cols2, width = 7.4, height = 4.3)
cat("wrote", list.files(outdir, full.names = TRUE), sep = "\n")
