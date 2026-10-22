## Usage: Rscript sim/stability_collect.R
## Summarises results/stability/*/data###.rds into results/stability/
##   stability_by_data.csv  one row per scenario x data set x q
##   stability_summary.csv  one row per scenario x q
##   table_stability.tex    LaTeX table for the paper / supplement
## Jaccard index J(A, B) = |A & B| / |A | B|, with J(empty, empty) = 1.

SIM_ROOT <- normalizePath(Sys.getenv("CQR_SIM_ROOT", "."))
STAB_RULE <- Sys.getenv("CQR_STAB_RULE", "cv")
dir <- file.path(SIM_ROOT, "results", paste0("stability", if (STAB_RULE == "cv") "" else paste0("_", STAB_RULE)))
files <- list.files(dir, "^data[0-9]+\\.rds$", recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("no stability results found")
## optional two-run stability of MS-MS / MDS (results/stability_multi; not run for the paper, columns NA)
mdir <- file.path(SIM_ROOT, "results", "stability_multi")
multi <- function(scen, d, m, j, act) {
  f <- file.path(mdir, scen, sprintf("data%03d_%s.rds", d, m))
  if (!file.exists(f)) return(c(NA, NA, NA))
  r <- readRDS(f)$runs
  if (any(vapply(r, `[[`, logical(1), "failed"))) return(c(NA, NA, NA))
  A <- r[[1]]$by_q[[j]]$selected; B <- r[[2]]$by_q[[j]]$selected
  c(jac(A, B), mean(c(fdp(A, act), fdp(B, act))), mean(c(tpp(A, act), tpp(B, act))))
}
jac <- function(a, b) { u <- length(union(a, b)); if (u == 0) 1 else length(intersect(a, b)) / u }
fdp <- function(S, act) if (length(S)) mean(!S %in% act) else 0
tpp <- function(S, act) mean(act %in% S)

rows <- list()
for (f in files) {
  x <- readRDS(f); act <- x$active
  ok <- !vapply(x$splits, `[[`, logical(1), "failed")
  for (j in seq_along(x$splits[[which(ok)[1]]]$by_q)) {
    q <- x$splits[[which(ok)[1]]]$by_q[[j]]$q
    S <- lapply(x$splits[ok], function(s) s$by_q[[j]]$selected)
    pr <- utils::combn(length(S), 2)
    J <- apply(pr, 2, function(ix) jac(S[[ix[1]]], S[[ix[2]]]))
    F <- vapply(S, fdp, 0, act = act); P <- vapply(S, tpp, 0, act = act)
    nulls_ever <- setdiff(unique(unlist(S)), act)
    m <- x$mds
    if (is.null(m)) m <- list(list(failed = TRUE), list(failed = TRUE))
    mds_ok <- !vapply(m, `[[`, logical(1), "failed")
    Jm <- if (all(mds_ok)) jac(m[[1]]$by_q[[j]]$selected, m[[2]]$by_q[[j]]$selected) else NA
    Fm <- if (mds_ok[1]) fdp(m[[1]]$by_q[[j]]$selected, act) else NA
    rows[[length(rows) + 1L]] <- data.frame(
      scenario = x$scenario, data = x$data, q = q, splits = sum(ok),
      jaccard_mean = mean(J), jaccard_q10 = unname(quantile(J, 0.1)),
      identical_share = mean(J == 1),
      fdp_mean = mean(F), fdp_sd = sd(F), fdp_max = max(F),
      power_mean = mean(P), power_sd = sd(P),
      signal_sel_min = min(vapply(act, function(a) mean(vapply(S, function(s) a %in% s, TRUE)), 0)),
      nulls_ever = length(nulls_ever),
      null_sel_max = if (length(nulls_ever)) max(vapply(nulls_ever, function(a) mean(vapply(S, function(s) a %in% s, TRUE)), 0)) else 0,
      mds_jaccard = Jm, mds_fdp = Fm,
      msms_jaccard = multi(x$scenario, x$data, "msms", j, act)[1],
      msms_fdp = multi(x$scenario, x$data, "msms", j, act)[2],
      msms_power = multi(x$scenario, x$data, "msms", j, act)[3],
      mdsdai_jaccard = multi(x$scenario, x$data, "mds_dai", j, act)[1],
      stringsAsFactors = FALSE)
  }
}
by_data <- do.call(rbind, rows)
write.csv(by_data, file.path(dir, "stability_by_data.csv"), row.names = FALSE)
agg <- function(g) data.frame(
  scenario = g$scenario[1], q = g$q[1], datasets = nrow(g),
  jaccard_mean = mean(g$jaccard_mean), jaccard_q10 = mean(g$jaccard_q10),
  identical_share = mean(g$identical_share),
  fdp_mean = mean(g$fdp_mean), fdp_sd_within = mean(g$fdp_sd), fdp_max = max(g$fdp_max),
  power_mean = mean(g$power_mean), power_sd_within = mean(g$power_sd),
  signal_sel_min = mean(g$signal_sel_min), null_sel_max = mean(g$null_sel_max),
  mds_jaccard = mean(g$mds_jaccard, na.rm = TRUE), mds_fdp = mean(g$mds_fdp, na.rm = TRUE),
  msms_datasets = sum(!is.na(g$msms_jaccard)), msms_jaccard = mean(g$msms_jaccard, na.rm = TRUE),
  msms_fdp = mean(g$msms_fdp, na.rm = TRUE), msms_power = mean(g$msms_power, na.rm = TRUE),
  mdsdai_jaccard = mean(g$mdsdai_jaccard, na.rm = TRUE))
summ <- do.call(rbind, lapply(split(by_data, list(by_data$scenario, by_data$q), drop = TRUE), agg))
summ <- summ[order(summ$scenario, summ$q), ]
write.csv(summ, file.path(dir, "stability_summary.csv"), row.names = FALSE)

lab <- c(T_p1000_strong = "$p=1000$, $|\\beta_j|=0.8$", T_p1000_weak = "$p=1000$, $|\\beta_j|=0.3$",
         T_p100_strong = "$p=100$, $|\\beta_j|=0.8$")
f3 <- function(v) sprintf("%.3f", v)
body <- vapply(seq_len(nrow(summ)), function(i) with(summ[i, ], paste(
  if (scenario %in% names(lab)) lab[[scenario]] else gsub("_", "\\\\_", scenario), sprintf("%.1f", q),
  f3(jaccard_mean), f3(identical_share), f3(fdp_mean), f3(fdp_sd_within), f3(power_mean),
  f3(null_sel_max), f3(mds_jaccard), sep = " & ")), character(1))
tex <- c("\\begin{table}[!htbp]", "\\centering\\footnotesize",
  "\\caption{Split-to-split stability of SS-MS: each data set is analysed with 50 random splits (averages over data sets).}",
  "\\label{tab:sim-stability}", "\\begin{tabular}{llrrrrrrr}", "\\toprule",
  "Setting & $q$ & Jaccard & Identical & FDP & SD(FDP) & Power & Max null freq. & MDS Jaccard \\\\", "\\midrule",
  paste0(body, " \\\\"), "\\bottomrule", "\\end{tabular}",
  "\\par\\smallskip\\parbox{\\linewidth}{\\scriptsize Jaccard: mean pairwise Jaccard index of the selected sets across splits of the same data set. Identical: share of split pairs with identical selections. SD(FDP): standard deviation of the false discovery proportion across splits. Max null freq.: largest share of splits selecting any one null covariate. MDS Jaccard: Jaccard index between two MDS runs ($B=50$ each) with different split seeds.}",
  "\\end{table}")
writeLines(tex, file.path(dir, "table_stability.tex"))
print(summ, digits = 3)
