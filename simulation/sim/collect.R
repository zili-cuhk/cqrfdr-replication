## Usage: [CQR_SIM_SET=main|df_pilot] Rscript sim/collect.R
## Reads results/<set>/*/rep*.rds and writes into results/<set>/
##   replications.csv  one row per scenario x replication x method
##   results/summary.csv       FDR, power, MCSE, failure counts and bounds
##   results/summary.tex       LaTeX table body (FDR / power, valid counts)
##   results/figures.pdf       FDR and power by scenario block

SIM_ROOT <- normalizePath(Sys.getenv("CQR_SIM_ROOT", "."))
source(file.path(SIM_ROOT, "sim", "config.R"))
files <- list.files(results_dir(SIM_ROOT), pattern = "^rep[0-9]+\\.rds$",
                    recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("no result files found")
rep <- do.call(rbind, lapply(files, function(f) readRDS(f)$metrics))
write.csv(rep, file.path(results_dir(SIM_ROOT), "replications.csv"), row.names = FALSE)

sc <- make_scenarios()
summ <- do.call(rbind, lapply(split(rep, list(rep$scenario, rep$method), drop = TRUE), function(g) {
  ok <- !g$failed; n <- nrow(g); m <- sum(ok)
  mc <- function(x) if (length(x) > 1L) sd(x) / sqrt(length(x)) else NA
  data.frame(scenario = g$scenario[1], method = g$method[1], base_method = g$base_method[1],
             q = g$q[1], reps = n, valid = m,
             failed = n - m,
             FDR = if (m) mean(g$FDP[ok]) else NA, FDR_mcse = mc(g$FDP[ok]),
             power = if (m) mean(g$TPP[ok]) else NA, power_mcse = mc(g$TPP[ok]),
             ## failure-inclusive bounds: failed replications counted as 0 or 1
             FDR_lo = sum(g$FDP[ok]) / n, FDR_hi = (sum(g$FDP[ok]) + n - m) / n,
             power_lo = sum(g$TPP[ok]) / n, power_hi = (sum(g$TPP[ok]) + n - m) / n,
             mean_R = if (m) mean(g$R[ok]) else NA,
             no_threshold = sum(g$status == "no_threshold"),
             screen_recall = mean(g$screen_recall[ok], na.rm = TRUE),
             screen_size = mean(g$screen_size[ok], na.rm = TRUE),
             median_seconds = median(g$seconds),
             top_failure = if (m < n) names(sort(table(g$status[!ok]), decreasing = TRUE))[1] else "",
             stringsAsFactors = FALSE)
}))
summ <- merge(sc, summ, by = "scenario")
summ <- summ[order(summ$scenario, summ$method), ]
write.csv(summ, file.path(results_dir(SIM_ROOT), "summary.csv"), row.names = FALSE)

## LaTeX body: one row per scenario, FDR & power per method (valid count as superscript)
## method labels in the order of SIM_CONFIG$methods, then by q
lab <- unique(summ[, c("method", "base_method", "q")])
lab <- lab[order(match(lab$base_method, SIM_CONFIG$methods), lab$q), ]
meths <- lab$method
fmt <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x))
tex <- vapply(unique(summ$scenario), function(s) {
  cells <- vapply(meths, function(mm) {
    r <- summ[summ$scenario == s & summ$method == mm, ]
    if (!nrow(r)) return("-- & --")
    sup <- if (r$valid < r$reps) sprintf("$^{(%d)}$", r$valid) else ""
    paste0(fmt(r$FDR), sup, " & ", fmt(r$power))
  }, character(1))
  paste0(gsub("_", "\\\\_", s), " & ", paste(cells, collapse = " & "), " \\\\")
}, character(1))
writeLines(tex, file.path(results_dir(SIM_ROOT), "summary.tex"))

## Figures (base graphics; no extra packages)
pdf(file.path(results_dir(SIM_ROOT), "figures.pdf"), width = 9, height = 6)
cols <- c(ssms = "#0072B2", ssms_size = "#56B4E9", mds = "#009E73", lhdcqr = "#D55E00",
          scqr_lasso = "#CC79A7", scqr_scad = "#E69F00")
for (blk in unique(substr(summ$scenario, 1, 1))) {
  d <- summ[substr(summ$scenario, 1, 1) == blk, ]
  scen <- unique(d$scenario)
  op <- par(mfrow = c(1, 2), mar = c(9, 4, 3, 1))
  for (what in c("FDR", "power")) {
    plot(NA, xlim = c(0.5, length(scen) + 0.5), ylim = c(0, 1), xaxt = "n",
         xlab = "", ylab = what, main = paste("Block", blk, "-", what))
    axis(1, at = seq_along(scen), labels = scen, las = 2, cex.axis = 0.7)
    if (what == "FDR") abline(h = SIM_CONFIG$q, lty = 2, col = "grey40")
    off <- seq(-0.25, 0.25, length.out = length(meths))
    for (j in seq_along(meths)) {
      r <- d[d$method == meths[j], ]
      x <- match(r$scenario, scen) + off[j]
      y <- r[[what]]; se <- r[[paste0(what, "_mcse")]]
      cl <- cols[lab$base_method[j]]
      points(x, y, pch = 14 + j, col = cl)
      segments(x, y - se, x, y + se, col = cl)
    }
    legend("topright", legend = meths, col = cols[lab$base_method], pch = 14 + seq_along(meths),
           cex = 0.6, bg = "white")
  }
  par(op)
}
invisible(dev.off())
cat("wrote", file.path(results_dir(SIM_ROOT), c("summary.csv", "summary.tex", "figures.pdf")), sep = "\n")
