## ---------------------------------------------------------------------------
## Section 5 of the paper: every table and every number in the text, from the
## collected result files (results_paper/, the output of R/04_collect.R and
## plasmode/plasmode_collect.R on the HPC runs).
## ---------------------------------------------------------------------------
## Usage (from the realdata/ folder):  Rscript R/11_paper_tables.R [results folder]
## Prints Tables 1-3, the numbers quoted for METABRIC (targeted) and the two
## supplementary tables, and writes them as CSV files to <results>/paper_tables/.
##
## Method labels in the result files -> names in the paper:
##   "MS-MS", "SS-MS"             -> MS-MS, SS-MS
##   "MDS-Dai"                    -> MDS   ("MDS" = an earlier implementation, not used)
##   "L-HDCQR-Z18"                -> L-HDCQR (Zheng et al. 2018)  ("L-HDCQR" = single
##                                   cross-validated penalty, not used)
##   "DS", "SCQR-Lasso", "SCQR-SCAD", "Cox-lasso"
##   plasmode: "SS-MS (path, U+)" -> SS-MS

args <- commandArgs(TRUE)
RES <- if (length(args)) args[1] else "results_paper"
OUT <- file.path(RES, "paper_tables"); dir.create(OUT, showWarnings = FALSE)
rd <- function(...) utils::read.csv(file.path(RES, ...), stringsAsFactors = FALSE)
show <- function(title, df) { cat("\n==== ", title, " ====\n", sep = ""); print(df, row.names = FALSE) }
save <- function(df, name) utils::write.csv(df, file.path(OUT, paste0(name, ".csv")), row.names = FALSE)
split_genes <- function(s) if (is.na(s) || !nzchar(s)) character(0) else strsplit(s, ";")[[1]]
jac <- function(a, b) { u <- length(union(a, b)); if (!u) 1 else length(intersect(a, b)) / u }
r3 <- function(x) round(x, 3)

## ---- Table 1: plasmode on the SUPPORT design ------------------------------------
pm <- rd("plasmode", "plasmode_summary_SUPPORT2.csv")
pm$q <- suppressWarnings(as.numeric(pm$q))
cell <- function(meth, qq) {
  z <- pm[pm$method == meth & (if (is.na(qq)) is.na(pm$q) else (!is.na(pm$q) & pm$q == qq)), ]
  z[order(z$s, z$b), c("s", "b", "FDR", "power", "FDR_se", "power_se")]
}
t1 <- cell("SS-MS (path, U+)", 0.1)[, c("s", "b")]
for (k in list(c("SS-MS (path, U+)", 0.1, "SSMS_q01"), c("MS-MS", 0.1, "MSMS_q01"),
               c("SS-MS (path, U+)", 0.2, "SSMS_q02"), c("MS-MS", 0.2, "MSMS_q02"),
               c("L-HDCQR-Z18", NA, "LHDCQR"))) {
  z <- cell(k[1], as.numeric(k[2]))
  t1[[paste0(k[3], "_FDR")]] <- z$FDR; t1[[paste0(k[3], "_power")]] <- z$power
}
used <- pm[pm$method %in% c("SS-MS (path, U+)", "MS-MS", "L-HDCQR-Z18"), ]
show("Table 1: plasmode simulation on the SUPPORT design (tau = 0.5)", t1)
cat(sprintf("largest MCSE: FDR %.3f, power %.3f\n", max(used$FDR_se), max(used$power_se)))
sc <- rd("plasmode", "plasmode_screen_SUPPORT2.csv"); sc <- sc[order(sc$s, sc$b), ]
cat("screen recall (path):", paste(sc$path_recall, collapse = ", "), "\n")
save(t1, "table1_plasmode_support")

## ---- Table 2: SUPPORT selections and diagnostics -----------------------------------
ms <- rd("tables_support", "main_selection.csv")
nsel <- function(tau, m, q = NA) {
  z <- ms[ms$tau == tau & ms$method == m & (if (is.na(q)) is.na(ms$q) else (!is.na(ms$q) & ms$q == q)), ]
  if (nrow(z)) z$n_selected[1] else NA
}
genes <- function(tau, m, q = NA) {
  z <- ms[ms$tau == tau & ms$method == m & (if (is.na(q)) is.na(ms$q) else (!is.na(ms$q) & ms$q == q)), ]
  if (nrow(z)) split_genes(z$genes[1]) else character(0)
}
mx <- rd("tables_support", "multi_split.csv")
st <- rd("tables_support", "stability.csv")
dc <- rd("tables_support", "decoy.csv")
ho <- rd("tables_support", "holdout.csv")
cols <- expand.grid(q = c(0.1, 0.2), tau = c(0.25, 0.5))
row_q <- function(f) vapply(seq_len(nrow(cols)), function(i) f(cols$tau[i], cols$q[i]), 0)
pick <- function(df, tau, m, q = NA, col) {
  z <- df[df$tau == tau & df$method == m & (if (is.na(q)) is.na(df$q) else (!is.na(df$q) & df$q == q)), col]
  if (length(z)) z[1] else NA
}
t2 <- rbind(
  `n: SS-MS`        = row_q(function(t, q) nsel(t, "SS-MS", q)),
  `n: MS-MS run 1`  = row_q(function(t, q) mx$msms[mx$tau == t & mx$q == q]),
  `n: MS-MS run 2`  = row_q(function(t, q) mx$msms_run2[mx$tau == t & mx$q == q]),
  `n: DS`           = row_q(function(t, q) nsel(t, "DS", q)),
  `n: MDS`          = row_q(function(t, q) nsel(t, "MDS-Dai", q)),
  `n: L-HDCQR`      = row_q(function(t, q) nsel(t, "L-HDCQR-Z18")),
  `n: SCQR-Lasso`   = row_q(function(t, q) nsel(t, "SCQR-Lasso")),
  `n: SCQR-SCAD`    = row_q(function(t, q) nsel(t, "SCQR-SCAD")),
  `n: Cox lasso`    = row_q(function(t, q) nsel(t, "Cox-lasso")),
  `Jaccard: SS-MS 100 splits` = row_q(function(t, q) st$jaccard[st$tau == t & st$q == q]),
  `Jaccard: MS-MS two runs`   = row_q(function(t, q) mx$msms_jaccard_runs[mx$tau == t & mx$q == q]),
  `MS-MS share of empty splits` = row_q(function(t, q) mx$msms_share_empty_splits[mx$tau == t & mx$q == q]),
  `decoys: SS-MS`   = row_q(function(t, q) pick(dc, t, "SS-MS", q, "decoy_share")),
  `decoys: MS-MS`   = row_q(function(t, q) pick(dc, t, "MS-MS", q, "decoy_share")),
  `decoys: L-HDCQR` = row_q(function(t, q) pick(dc, t, "L-HDCQR-Z18", NA, "decoy_share")),
  `decoys: Cox`     = row_q(function(t, q) pick(dc, t, "Cox-lasso", NA, "decoy_share")),
  `C: SS-MS`        = row_q(function(t, q) pick(ho, t, "SS-MS", q, "C")),
  `C: MS-MS`        = row_q(function(t, q) pick(ho, t, "MS-MS", q, "C")),
  `C: L-HDCQR`      = row_q(function(t, q) pick(ho, t, "L-HDCQR-Z18 (refit)", NA, "C")),
  `C: Cox`          = row_q(function(t, q) pick(ho, t, "Cox-lasso", NA, "C")),
  `n in hold-out: SS-MS`   = row_q(function(t, q) pick(ho, t, "SS-MS", q, "mean_selected")),
  `n in hold-out: MS-MS`   = row_q(function(t, q) pick(ho, t, "MS-MS", q, "mean_selected")),
  `n in hold-out: L-HDCQR` = row_q(function(t, q) pick(ho, t, "L-HDCQR-Z18 (refit)", NA, "mean_selected")),
  `n in hold-out: Cox`     = row_q(function(t, q) pick(ho, t, "Cox-lasso", NA, "mean_selected")))
colnames(t2) <- sprintf("tau=%.2f,q=%.1f", cols$tau, cols$q)
t2 <- data.frame(row = rownames(t2), round(t2, 3), check.names = FALSE)
show("Table 2: SUPPORT data (L-HDCQR, SCQR and Cox do not depend on q)", t2)
save(t2, "table2_support")
cat("\nSD of the hold-out C-index at tau = 0.25, q = 0.1: SS-MS",
    pick(ho, 0.25, "SS-MS", 0.1, "C_sd"), " MS-MS", pick(ho, 0.25, "MS-MS", 0.1, "C_sd"), "\n")
cat("any decoy selected (share of replications): Cox",
    paste(dc$any_decoy[dc$method == "Cox-lasso"], collapse = "/"), "\n")
cat("\nAgreement of MS-MS (q = 0.1) with other methods:\n")
for (t in c(0.25, 0.5)) {
  a <- genes(t, "MS-MS", 0.1)
  cat(sprintf("  tau %.2f: Jaccard with DS %.2f, with MDS %.2f; %d of L-HDCQR's %d also in MS-MS\n", t,
              jac(a, genes(t, "DS", 0.1)), jac(a, genes(t, "MDS-Dai", 0.1)),
              sum(genes(t, "L-HDCQR-Z18") %in% a), length(genes(t, "L-HDCQR-Z18"))))
}
cat("\nRun time per tau (minutes; SS-MS, DS, L-HDCQR, SCQR, MS-MS on 4 cores, MDS on 10):\n")
tm <- ms[ms$method %in% c("SS-MS", "MS-MS", "DS", "MDS-Dai", "L-HDCQR-Z18", "SCQR-Lasso", "SCQR-SCAD") &
         (is.na(ms$q) | ms$q == 0.1), c("tau", "method", "seconds")]
tm$minutes <- round(tm$seconds / 60, 1); print(tm[, c("tau", "method", "minutes")], row.names = FALSE)

## ---- Table 3: main effects and squares selected by MS-MS at q = 0.1 ------------------
a <- genes(0.25, "MS-MS", 0.1); b <- genes(0.5, "MS-MS", 0.1)
is_main <- function(v) !grepl(":", v)
t3 <- list(both = intersect(a, b)[is_main(intersect(a, b))],
           tau025_only = setdiff(a, b)[is_main(setdiff(a, b))],
           tau050_only = setdiff(b, a)[is_main(setdiff(b, a))])
cat("\n==== Table 3: MS-MS (q = 0.1), main effects and squared terms ====\n")
for (k in names(t3)) cat(sprintf("%-12s %s\n", k, paste(t3[[k]], collapse = ", ")))
cat(sprintf("pairwise products: %d at tau = 0.25, %d at tau = 0.5, %d at both; Jaccard between the levels %.2f\n",
            sum(!is_main(a)), sum(!is_main(b)), sum(!is_main(intersect(a, b))), jac(a, b)))
save(data.frame(set = rep(names(t3), lengths(t3)), variable = unlist(t3)), "table3_support_msms_variables")

## ---- METABRIC, targeted (Section 5.4) -------------------------------------------------
mt <- rd("tables_targeted", "main_selection.csv")
tt <- stats::aggregate(n_selected ~ method + q, data = transform(mt, q = ifelse(is.na(q), "-", q)),
                       FUN = function(v) paste(v, collapse = "/"))
show("METABRIC targeted: number selected at tau = 0.25 / 0.5", tt)
cat("MS-MS selections (q = 0.2):", unique(mt$genes[mt$method == "MS-MS" & mt$q == 0.2]), "\n")
dt <- rd("tables_targeted", "decoy.csv"); ht <- rd("tables_targeted", "holdout.csv")
show("METABRIC targeted: decoys", dt[, c("tau", "method", "q", "decoy_share", "any_decoy")])
show("METABRIC targeted: hold-out", ht[ht$method %in% c("MS-MS", "L-HDCQR-Z18 (refit)"), c("tau", "method", "q", "C", "mean_selected")])

## ---- Supplementary: plasmode on three designs ------------------------------------------
sp <- do.call(rbind, lapply(c("BRCA_METABRIC_TARGETED", "BRCA_METABRIC", "GLIOMA_TCGA"), function(nm) {
  x <- rd("plasmode", paste0("plasmode_summary_", nm, ".csv"))
  r <- rd("plasmode", paste0("plasmode_screen_", nm, ".csv"))
  g <- unique(x[, c("s", "b")]); g <- g[order(g$s, g$b), ]
  do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
    f <- function(m, q, w) { z <- x[x$s == g$s[i] & x$b == g$b[i] & x$method == m & x$q == q, w]; if (length(z)) z else NA }
    data.frame(design = nm, s = g$s[i], b = g$b[i],
               recall = r$path_recall[r$s == g$s[i] & r$b == g$b[i]],
               SSMS_q01 = sprintf("%.3f/%.3f", f("SS-MS (path, U+)", "0.1", "FDR"), f("SS-MS (path, U+)", "0.1", "power")),
               MSMS_q01 = sprintf("%.3f/%.3f", f("MS-MS", "0.1", "FDR"), f("MS-MS", "0.1", "power")),
               SSMS_q02 = sprintf("%.3f/%.3f", f("SS-MS (path, U+)", "0.2", "FDR"), f("SS-MS (path, U+)", "0.2", "power")),
               MSMS_q02 = sprintf("%.3f/%.3f", f("MS-MS", "0.2", "FDR"), f("MS-MS", "0.2", "power")))
  }))
}))
show("Supplementary table: plasmode on three designs (FDR/power)", sp)
save(sp, "supp_plasmode_three_designs")

## ---- Supplementary: genome-wide cohorts ------------------------------------------------
mg <- rbind(rd("tables_genome", "multi_split.csv"), rd("tables_glioma", "multi_split.csv"))
el <- rd("tables_genome", "eligibility.csv")
gd <- rd("tables_glioma", "glioma_design.csv")
el <- rbind(el[, c("dataset", "n", "events", "screen_cap")], gd[, c("dataset", "n", "events", "screen_cap")])
mg <- mg[mg$q == 0.2, ]
sg <- do.call(rbind, lapply(unique(mg$dataset), function(ds) {
  z <- mg[mg$dataset == ds, ]; e <- el[el$dataset == ds, ]
  f <- function(t, w) z[z$tau == t, w]
  data.frame(cohort = ds, n = if (nrow(e)) e$n else NA, events = if (nrow(e)) e$events else NA,
             cap = if (nrow(e)) e$screen_cap else NA,
             SSMS_025 = f(0.25, "ssms"), MSMS_025 = f(0.25, "msms"), empty_025 = r3(f(0.25, "msms_share_empty_splits")),
             J_025 = r3(f(0.25, "msms_jaccard_runs")),
             SSMS_050 = f(0.5, "ssms"), MSMS_050 = f(0.5, "msms"), empty_050 = r3(f(0.5, "msms_share_empty_splits")),
             J_050 = r3(f(0.5, "msms_jaccard_runs")))
}))
sg <- sg[order(-sg$events), ]
show("Supplementary table: genome-wide cohorts (q = 0.2); J = 1 means both runs empty", sg)
save(sg, "supp_genome_wide")
cat(sprintf("\nq = 0.1: largest number selected by SS-MS / MS-MS over all cohorts: %d / %d\n",
            max(rbind(rd("tables_genome", "multi_split.csv"), rd("tables_glioma", "multi_split.csv"))[, "ssms"][
              rbind(rd("tables_genome", "multi_split.csv"), rd("tables_glioma", "multi_split.csv"))$q == 0.1]),
            max(rbind(rd("tables_genome", "multi_split.csv"), rd("tables_glioma", "multi_split.csv"))[, "msms"][
              rbind(rd("tables_genome", "multi_split.csv"), rd("tables_glioma", "multi_split.csv"))$q == 0.1])))
cat("\nwritten:", paste(list.files(OUT), collapse = ", "), "\n")
