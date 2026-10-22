## Usage: Rscript sim/paper_table_stability.R [results folder] [outdir]
## Supplementary table on split-to-split variability: SS-MS with the
## path-ordered screen (results/stability_path, CQR_STAB_RULE=path) and the
## Jaccard index of SS-MS with the CV screen (results/stability, the default
## rule), both produced by sim/stability.R and sim/stability_collect.R.
args <- commandArgs(TRUE)
res <- if (length(args) >= 1) args[1] else "results"
outdir <- if (length(args) >= 2) args[2] else file.path(res, "paper")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
pa <- read.csv(file.path(res, "stability_path", "stability_summary.csv"))
cv <- read.csv(file.path(res, "stability", "stability_summary.csv"))
lab <- c(T_p1000_strong = "$p=1000$, $|\\beta_j|=0.8$", T_p100_strong = "$p=100$, $|\\beta_j|=0.8$",
         T_p1000_weak = "$p=1000$, $|\\beta_j|=0.3$")
rows <- character()
for (sc in names(lab)) for (q in c(0.1, 0.2)) {
  a <- pa[pa$scenario == sc & pa$q == q, ]; b <- cv[cv$scenario == sc & cv$q == q, ]
  rows <- c(rows, sprintf("%s & %.1f & %.2f & %.2f & %.2f (%.2f) & %.2f & %.2f & %.2f \\\\",
                          if (q == 0.1) lab[[sc]] else "", q, a$jaccard_mean, a$fdp_sd_within, a$power_mean,
                          a$power_sd_within, a$signal_sel_min, a$null_sel_max, b$jaccard_mean))
}
tex <- c("\\begin{tabular}{llcccccc}", "\\toprule",
         " & & \\multicolumn{5}{c}{SS-MS} & SS-MS, CV screen \\\\",
         "\\cmidrule(lr){3-7}\\cmidrule(lr){8-8}",
         "Setting & $q$ & Jaccard & SD(FDP) & Power (SD) & Min signal & Max null & Jaccard \\\\", "\\midrule",
         rows, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(outdir, "table_stability.tex"))
cat(tex, sep = "\n")
