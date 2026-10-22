## ---------------------------------------------------------------------------
## Collect all task results into tables
## ---------------------------------------------------------------------------
## Usage: Rscript R/04_collect.R
## tables/
##   eligibility.csv          (from 02_prepare.R; includes the pre-specified role)
##   main_selection.csv       selected genes per data set, tau, method, q
##   tau_overlap.csv          SS-MS selections at tau = 0.25 vs 0.5
##   stability.csv            split-to-split stability of SS-MS
##   decoy.csv                permuted-gene (null) selections
##   holdout.csv              70/30 prediction: C-index and held-out Peng-Huang loss
##   external.csv             external validation (refit on discovery, C-index on validation)
##   summary_by_dataset.csv   one row per data set
##   supp_realdata.tex        LaTeX table for the supplement
##   progress.csv             finished / expected tasks

source(file.path(Sys.getenv("CQR_RD_ROOT", "."), "R", "config.R"))
suppressPackageStartupMessages({ library(cqrfdr); library(survival) })
TB <- TABLES_DIR; dir.create(TB, showWarnings = FALSE, recursive = TRUE)
el <- utils::read.csv(file.path(RD_ROOT, "tables", "eligibility.csv"), stringsAsFactors = FALSE)
tasks <- utils::read.csv(TASKS_FILE, stringsAsFactors = FALSE)
rpath <- function(ds, tau, type, rep) file.path(RD_ROOT, "results", ds, sprintf("tau%.2f", tau),
                                               sprintf("%s_%03d.rds", type, rep))
tasks$done <- file.exists(mapply(rpath, tasks$dataset, tasks$tau, tasks$type, tasks$rep))
prog <- stats::aggregate(done ~ dataset + tau + type, tasks, function(v) sprintf("%d/%d", sum(v), length(v)))
write.csv(prog, file.path(TB, "progress.csv"), row.names = FALSE)
cat(sprintf("finished %d of %d tasks\n", sum(tasks$done), nrow(tasks)))
read_all <- function(ds, tau, type) {
  t <- tasks[tasks$dataset == ds & tasks$tau == tau & tasks$type == type & tasks$done, ]
  lapply(seq_len(nrow(t)), function(i) readRDS(rpath(ds, tau, type, t$rep[i])))
}
qs <- format(CFG$q)
dt <- unique(tasks[, c("dataset", "tau")])
sel_ok <- function(r) !is.null(r) && is.null(r$error) && !is.null(r$selected) && !identical(r$status, "glmnet_not_installed")
jac <- function(a, b) { u <- length(union(a, b)); if (u == 0) 1 else length(intersect(a, b)) / u }
g <- function(v) paste(v, collapse = ";")

## 1. main selections ---------------------------------------------------------------
ms <- list(); main <- list()
for (i in seq_len(nrow(dt))) {
  ds <- dt$dataset[i]; tau <- dt$tau[i]
  r <- read_all(ds, tau, "main"); if (!length(r)) next
  r <- r[[1L]]
  rm_ <- read_all(ds, tau, "mds"); if (length(rm_)) r$mds <- rm_[[1L]]$mds
  rr_ <- Filter(function(z) z$task$rep == 1L, read_all(ds, tau, "msms")); if (length(rr_)) r$msms <- rr_[[1L]]$msms
  rd_ <- read_all(ds, tau, "mds_dai"); if (length(rd_)) r$mds_dai <- rd_[[1L]]$mds_dai
  rc_ <- read_all(ds, tau, "cmp"); if (length(rc_)) for (m in grep("^cmp_", names(rc_[[1L]]), value = TRUE)) r[[m]] <- rc_[[1L]][[m]]
  main[[paste(ds, tau)]] <- r
  add <- function(method, q, genes, sec) ms[[length(ms) + 1L]] <<- data.frame(
    dataset = ds, tau = tau, method = method, q = q, n_selected = length(genes),
    genes = g(genes), seconds = round(sec), stringsAsFactors = FALSE)
  lab_m <- c(ssms = "SS-MS", msms = "MS-MS", mds = "MDS", mds_dai = "MDS-Dai")   # MDS = earlier implementation
  for (m in names(lab_m)) if (sel_ok(r[[m]])) for (q in qs)
    add(lab_m[[m]], q, r[[m]]$selected[[q]], r[[m]]$seconds)
  if (sel_ok(r$hdcqr)) add("L-HDCQR", NA, r$hdcqr$selected, r$hdcqr$seconds)
  ## comparators of the simulation study (task "cmp"); L-HDCQR-Z18 = Zheng et al. (2018)
  if (sel_ok(r$cmp_ds)) for (q in qs) add("DS", q, r$cmp_ds$selected[[q]], r$cmp_ds$seconds)
  cm <- c(lhdcqr = "L-HDCQR-Z18", scqr_lasso = "SCQR-Lasso", scqr_scad = "SCQR-SCAD")
  for (k in names(cm)) { z <- r[[paste0("cmp_", k)]]; if (sel_ok(z)) add(cm[[k]], NA, z$selected, z$seconds) }
  if (sel_ok(r$cox)) add("Cox-lasso", NA, r$cox$selected, r$cox$seconds)
}
ms <- if (length(ms)) do.call(rbind, ms) else data.frame()
write.csv(ms, file.path(TB, "main_selection.csv"), row.names = FALSE)

## 2. tau overlap (SS-MS) ------------------------------------------------------------------
ov <- list()
for (ds in unique(dt$dataset)) for (q in qs) {
  a <- main[[paste(ds, 0.25)]]$ssms$selected[[q]]; b <- main[[paste(ds, 0.5)]]$ssms$selected[[q]]
  if (is.null(main[[paste(ds, 0.25)]]) || is.null(main[[paste(ds, 0.5)]])) next
  ov[[length(ov) + 1L]] <- data.frame(dataset = ds, q = q, n_tau025 = length(a), n_tau050 = length(b),
    both = length(intersect(a, b)), only_tau025 = g(setdiff(a, b)), only_tau050 = g(setdiff(b, a)),
    jaccard = round(jac(a, b), 3), stringsAsFactors = FALSE)
}
ov <- if (length(ov)) do.call(rbind, ov) else data.frame()
write.csv(ov, file.path(TB, "tau_overlap.csv"), row.names = FALSE)

## 3. stability --------------------------------------------------------------------------------
st <- list()
for (i in seq_len(nrow(dt))) {
  ds <- dt$dataset[i]; tau <- dt$tau[i]
  rr <- Filter(function(r) sel_ok(r$ssms), read_all(ds, tau, "stab")); if (length(rr) < 2L) next
  for (q in qs) {
    S <- lapply(rr, function(r) r$ssms$selected[[q]])
    pr <- utils::combn(length(S), 2L)
    J <- apply(pr, 2L, function(ix) jac(S[[ix[1]]], S[[ix[2]]]))
    fr <- table(unlist(S)) / length(S)
    msel <- main[[paste(ds, tau)]]$ssms$selected[[q]]
    st[[length(st) + 1L]] <- data.frame(dataset = ds, tau = tau, q = q, splits = length(S),
      mean_selected = round(mean(lengths(S)), 1), share_empty = round(mean(lengths(S) == 0), 3),
      jaccard = round(mean(J), 3), genes_freq50 = sum(fr >= 0.5),
      main_genes_median_freq = if (length(msel)) round(stats::median(as.numeric(fr[msel]), na.rm = TRUE), 3) else NA,
      top_genes = g(names(sort(-fr))[seq_len(min(15L, length(fr)))]), stringsAsFactors = FALSE)
  }
}
st <- if (length(st)) do.call(rbind, st) else data.frame()
write.csv(st, file.path(TB, "stability.csv"), row.names = FALSE)

## 4. decoys --------------------------------------------------------------------------------------
dc <- list()
for (i in seq_len(nrow(dt))) {
  ds <- dt$dataset[i]; tau <- dt$tau[i]
  rr <- read_all(ds, tau, "decoy"); if (!length(rr)) next
  one <- function(sets, method, q) {
    sets <- Filter(Negate(is.null), sets); if (!length(sets)) return(NULL)
    V <- vapply(sets, function(s) sum(startsWith(s, "DECOY_")), 0); R <- lengths(sets)
    data.frame(dataset = ds, tau = tau, method = method, q = q, reps = length(sets),
               mean_selected = round(mean(R), 1), mean_decoys = round(mean(V), 2),
               decoy_share = round(mean(V / pmax(R, 1)), 3), any_decoy = round(mean(V > 0), 3))
  }
  for (q in qs) dc[[length(dc) + 1L]] <- one(lapply(rr, function(r) if (sel_ok(r$ssms)) r$ssms$selected[[q]]), "SS-MS", q)
  dc[[length(dc) + 1L]] <- one(lapply(rr, function(r) if (sel_ok(r$hdcqr)) r$hdcqr$selected), "L-HDCQR", NA)
  rz <- read_all(ds, tau, "decoy_z18")
  if (length(rz)) dc[[length(dc) + 1L]] <- one(lapply(rz, function(r) if (sel_ok(r$lhdcqr_z18)) r$lhdcqr_z18$selected), "L-HDCQR-Z18", NA)
  rm2 <- read_all(ds, tau, "decoy_msms")
  if (length(rm2)) for (q in qs) dc[[length(dc) + 1L]] <- one(lapply(rm2, function(r) if (sel_ok(r$msms)) r$msms$selected[[q]]), "MS-MS", q)
  dc[[length(dc) + 1L]] <- one(lapply(rr, function(r) if (sel_ok(r$cox)) r$cox$selected), "Cox-lasso", NA)
}
dc <- if (length(dc)) do.call(rbind, dc) else data.frame()
write.csv(dc, file.path(TB, "decoy.csv"), row.names = FALSE)

## 5. hold-out prediction -------------------------------------------------------------------
ho <- list()
for (i in seq_len(nrow(dt))) {
  ds <- dt$dataset[i]; tau <- dt$tau[i]
  rr <- read_all(ds, tau, "holdout"); if (!length(rr)) next
  one <- function(v, method, q) {
    v <- Filter(function(z) is.numeric(z) && !is.null(names(z)), v); if (!length(v)) return(NULL)
    M <- do.call(rbind, lapply(v, function(z) z[c("n_sel", "C", "loss")]))
    data.frame(dataset = ds, tau = tau, method = method, q = q, reps = nrow(M),
               mean_selected = round(mean(M[, "n_sel"]), 1),
               C = round(mean(M[, "C"], na.rm = TRUE), 3), C_sd = round(stats::sd(M[, "C"], na.rm = TRUE), 3),
               loss = round(mean(M[, "loss"], na.rm = TRUE), 4))
  }
  for (q in qs) ho[[length(ho) + 1L]] <- one(lapply(rr, function(r) r$ssms[[q]]), "SS-MS", q)
  ho[[length(ho) + 1L]] <- one(lapply(rr, `[[`, "hdcqr"), "L-HDCQR (penalized)", NA)
  ho[[length(ho) + 1L]] <- one(lapply(rr, function(r) if (is.numeric(r$hdcqr))
    c(n_sel = r$hdcqr[["n_sel"]], C = r$hdcqr[["C_refit"]], loss = r$hdcqr[["loss_refit"]])), "L-HDCQR (refit)", NA)
  ho[[length(ho) + 1L]] <- one(lapply(rr, `[[`, "cox"), "Cox-lasso", NA)
  rz <- read_all(ds, tau, "holdout_z18")
  if (length(rz)) ho[[length(ho) + 1L]] <- one(lapply(rz, `[[`, "lhdcqr_z18"), "L-HDCQR-Z18 (refit)", NA)
  rm2 <- read_all(ds, tau, "holdout_msms")
  if (length(rm2)) for (q in qs) ho[[length(ho) + 1L]] <- one(lapply(rm2, function(r) r$msms[[q]]), "MS-MS", q)
  ho[[length(ho) + 1L]] <- one(lapply(rr, `[[`, "base"), "intercept only", NA)
}
ho <- if (length(ho)) do.call(rbind, ho) else data.frame()
write.csv(ho, file.path(TB, "holdout.csv"), row.names = FALSE)

## 6. external validation ------------------------------------------------------------------------
## Genes selected on the full discovery data (main task) are refitted without
## penalty on the discovery data (genes z-scored within each cohort); the
## predicted tau-quantile is evaluated by the C-index in each validation cohort.
zs <- function(x) scale(x)
load_raw <- function(nm) {
  f <- file.path(RD_ROOT, "data", "raw", paste0(nm, ".rds")); if (!file.exists(f)) return(NULL)
  r <- readRDS(f); x <- r$x
  if (anyNA(x)) for (j in which(colSums(is.na(x)) > 0)) x[is.na(x[, j]), j] <- stats::median(x[, j], na.rm = TRUE)
  r$x <- x; r
}
ex <- list()
for (k in seq_len(nrow(dt))) {
  ds <- dt$dataset[k]; tau <- dt$tau[k]; r <- main[[paste(ds, tau)]]
  vals <- VALIDATION[[ds]]; if (is.null(r) || !length(vals)) next
  D <- readRDS(file.path(RD_ROOT, "data", "clean", paste0(ds, ".rds")))
  sets <- list()
  if (sel_ok(r$ssms)) for (q in qs) sets[[paste0("SS-MS q=", q)]] <- r$ssms$selected[[q]]
  if (sel_ok(r$mds)) for (q in qs) sets[[paste0("MDS q=", q)]] <- r$mds$selected[[q]]
  if (sel_ok(r$msms)) for (q in qs) sets[[paste0("MS-MS q=", q)]] <- r$msms$selected[[q]]
  if (sel_ok(r$mds_dai)) for (q in qs) sets[[paste0("MDS-Dai q=", q)]] <- r$mds_dai$selected[[q]]
  if (sel_ok(r$hdcqr)) sets[["L-HDCQR"]] <- r$hdcqr$selected
  if (sel_ok(r$cmp_lhdcqr)) sets[["L-HDCQR-Z18"]] <- r$cmp_lhdcqr$selected
  if (sel_ok(r$cox)) sets[["Cox-lasso"]] <- r$cox$selected
  for (v in vals) {
    V <- load_raw(v); if (is.null(V)) next
    for (nm in names(sets)) {
      gs <- intersect(sets[[nm]], colnames(V$x))
      res <- data.frame(discovery = ds, validation = v, tau = tau, method = nm,
                        n_selected = length(sets[[nm]]), n_in_validation = length(gs),
                        C = NA_real_, status = "no genes", stringsAsFactors = FALSE)
      if (length(gs)) {
        taus <- seq(CFG$nu, tau, length.out = CFG$n_grid)
        f <- hdcqr_path(D$y, D$delta, cbind(1, zs(D$x[, gs, drop = FALSE])), taus, 0)
        res$status <- f$status
        if (f$status == "ok") {
          pred <- drop(cbind(1, zs(V$x[, gs, drop = FALSE])) %*% f$B[, length(taus)])
          res$C <- round(unname(survival::concordance(survival::Surv(V$time, V$status) ~ pred)$concordance), 3)
        }
      }
      ex[[length(ex) + 1L]] <- res
    }
  }
}
ex <- if (length(ex)) do.call(rbind, ex) else data.frame()
write.csv(ex, file.path(TB, "external.csv"), row.names = FALSE)

## 7. one row per data set ------------------------------------------------------------------------
pick <- function(df, ...) { cond <- list(...); for (n in names(cond)) df <- df[!is.na(df[[n]]) & df[[n]] == cond[[n]], , drop = FALSE]; df }
sm <- list()
for (i in which(el$analysed)) {
  ds <- el$dataset[i]
  for (tau in CFG$taus[c(el$tau025[i], el$tau050[i])]) {
    nsel <- function(m, q = NA) { z <- if (is.na(q)) pick(ms, dataset = ds, tau = tau, method = m) else pick(ms, dataset = ds, tau = tau, method = m, q = q); if (nrow(z)) z$n_selected[1] else NA }
    sm[[length(sm) + 1L]] <- data.frame(
      dataset = ds, role = el$role[i], n = el$n[i], events = el$events[i], tau = tau,
      ssms_q01 = nsel("SS-MS", "0.1"), ssms_q02 = nsel("SS-MS", "0.2"),
      mds_q01 = nsel("MDS", "0.1"), lhdcqr = nsel("L-HDCQR"), cox = nsel("Cox-lasso"),
      jaccard_q01 = { z <- pick(st, dataset = ds, tau = tau, q = "0.1"); if (nrow(z)) z$jaccard else NA },
      decoy_share_q01 = { z <- pick(dc, dataset = ds, tau = tau, method = "SS-MS", q = "0.1"); if (nrow(z)) z$decoy_share else NA },
      decoy_share_lhdcqr = { z <- pick(dc, dataset = ds, tau = tau, method = "L-HDCQR"); if (nrow(z)) z$decoy_share else NA },
      C_ssms_q01 = { z <- pick(ho, dataset = ds, tau = tau, method = "SS-MS", q = "0.1"); if (nrow(z)) z$C else NA },
      C_lhdcqr = { z <- pick(ho, dataset = ds, tau = tau, method = "L-HDCQR (penalized)"); if (nrow(z)) z$C else NA },
      stringsAsFactors = FALSE)
  }
}
sm <- if (length(sm)) do.call(rbind, sm) else data.frame()
write.csv(sm, file.path(TB, "summary_by_dataset.csv"), row.names = FALSE)

## 8. LaTeX table for the supplement ------------------------------------------------------------------
if (nrow(sm)) {
  esc <- function(s) gsub("_", "\\\\_", s)
  f <- function(v, d = 2) ifelse(is.na(v), "--", formatC(v, format = "f", digits = d))
  fi <- function(v) ifelse(is.na(v), "--", as.character(v))
  body <- with(sm, paste(esc(dataset), n, events, f(tau), fi(ssms_q01), fi(ssms_q02), fi(mds_q01),
                         fi(lhdcqr), fi(cox), f(jaccard_q01), f(decoy_share_q01), f(decoy_share_lhdcqr),
                         f(C_ssms_q01), f(C_lhdcqr), sep = " & "))
  tex <- c("\\begin{table}[!htbp]", "\\centering\\scriptsize",
    "\\caption{All eligible data sets, analysed with the same pre-specified pipeline. Numbers of selected genes on the full data; split-to-split Jaccard index of SS-MS; share of permuted (decoy) genes among the selections; hold-out C-index (70/30 splits).}",
    "\\label{tab:realdata-all}", "\\resizebox{\\linewidth}{!}{%", "\\begin{tabular}{lrrr rrrrr rrr rr}", "\\toprule",
    "& & & & \\multicolumn{5}{c}{Selected genes} & & \\multicolumn{2}{c}{Decoy share} & \\multicolumn{2}{c}{C-index}\\\\",
    "\\cmidrule(lr){5-9}\\cmidrule(lr){11-12}\\cmidrule(lr){13-14}",
    "Data set & $n$ & Events & $\\tau$ & SS-MS$_{0.1}$ & SS-MS$_{0.2}$ & MDS$_{0.1}$ & L-HDCQR & Cox & Jaccard & SS-MS$_{0.1}$ & L-HDCQR & SS-MS$_{0.1}$ & L-HDCQR \\\\",
    "\\midrule", paste0(body, " \\\\"), "\\bottomrule", "\\end{tabular}}", "\\end{table}")
  writeLines(tex, file.path(TB, "supp_realdata.tex"))
}
cat("\nmain-text data set (pre-specified rule):", el$dataset[which(el$role == "MAIN TEXT")], "\n")
print(sm, row.names = FALSE)

## 9. LaTeX table: all eligible data sets and the pre-specified analysis rule ----------------------
ee <- el[el$eligible, ]
if (nrow(ee)) {
  esc <- function(s) gsub("_", "\\\\_", s)
  body <- with(ee, paste(esc(dataset), esc(cancer), n, events, sprintf("%.2f", censoring),
                         screen_cap, ifelse(analysed, "yes", "no"),
                         ifelse(role == "MAIN TEXT", "main text", ifelse(analysed, ifelse(low_cap, "supplement (low cap)", "supplement"), "--")), sep = " & "))
  tex <- c("\\begin{table}[!htbp]", "\\centering\\scriptsize",
    sprintf("\\caption{Eligible data sets (n $\\geq$ %d, at least %d events, $\\tau=0.25$ identifiable). Every eligible data set with both $\\tau=0.25$ and $\\tau=0.5$ identifiable is analysed with the same pipeline. Cap: screening cap $\\lfloor d_2/\\log d_2\\rfloor$ with $d_2$ half the number of events; with a cap below %d, a selection at $q=0.1$ (at least 10 variables) needs nearly half of the screened set to be clear signals, so little power is expected. The main-text data set has the most events among the analysed data sets with an external validation cohort; the rule was fixed before any selection was computed.}",
            CFG$min_n, CFG$min_events, CFG$low_cap),
    "\\label{tab:realdata-eligibility}", "\\begin{tabular}{llrrrrll}", "\\toprule",
    "Data set & Cancer & $n$ & Events & Censoring & Cap & Analysed & Role \\\\", "\\midrule",
    paste0(body, " \\\\"), "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(tex, file.path(TB, "supp_eligibility.tex"))
}

## 10. single versus multiple splitting, all data sets of this task table ----------------------------
## Selections of SS-MS, MS-MS and MDS (Dai design; earlier implementation in
## brackets), run-to-run agreement of MS-MS (rep 1 vs rep 2), share of MS-MS
## splits whose own SS-MS selection was empty, and overlap with SS-MS.
mm <- list()
for (i in seq_len(nrow(dt))) {
  ds <- dt$dataset[i]; tau <- dt$tau[i]; r <- main[[paste(ds, tau)]]
  rr <- read_all(ds, tau, "msms")
  m1 <- Filter(function(z) z$task$rep == 1L, rr); m2 <- Filter(function(z) z$task$rep == 2L, rr)
  m1 <- if (length(m1)) m1[[1L]]$msms; m2 <- if (length(m2)) m2[[1L]]$msms
  for (q in qs) {
    gsel <- function(z) if (sel_ok(z)) z$selected[[q]] else NULL
    nz <- function(z) if (is.null(z)) NA else length(z)
    a <- gsel(r$ssms); b <- gsel(m1); b2 <- gsel(m2)
    mm[[length(mm) + 1L]] <- data.frame(
      dataset = ds, tau = tau, q = q, ssms = nz(a), msms = nz(b), msms_run2 = nz(b2),
      msms_jaccard_runs = if (!is.null(b) && !is.null(b2)) round(jac(b, b2), 3) else NA,
      msms_share_empty_splits = if (sel_ok(m1)) round(unname(m1$share_empty[q]), 3) else NA,
      msms_median_split_size = if (sel_ok(m1)) stats::median(m1$split_sizes[, q]) else NA,
      ssms_and_msms = if (!is.null(a) && !is.null(b)) length(intersect(a, b)) else NA,
      mds_dai = nz(gsel(r$mds_dai)), mds_earlier = nz(gsel(r$mds)),
      msms_genes = if (is.null(b)) "" else g(b), stringsAsFactors = FALSE)
  }
}
mm <- if (length(mm)) do.call(rbind, mm) else data.frame()
write.csv(mm, file.path(TB, "multi_split.csv"), row.names = FALSE)
if (nrow(mm)) { cat("\nSingle vs multiple splitting (numbers of selected variables):\n")
  print(mm[, setdiff(names(mm), "msms_genes")], row.names = FALSE) }
