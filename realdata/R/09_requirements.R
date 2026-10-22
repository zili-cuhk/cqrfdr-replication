## ---------------------------------------------------------------------------
## Requirement check for every analysed data set (run before and after jobs)
## ---------------------------------------------------------------------------
## Usage: Rscript R/09_requirements.R
## Reads every data/clean/*.rds and checks, without looking at any selection
## result:
##   A. data integrity : x numeric, finite, no constant or duplicated columns,
##                       unique gene names; time > 0 and finite; delta in {0,1};
##                       y = log(time); dimensions agree
##   B. eligibility    : n >= min_n, events >= min_events, tau = 0.25 identifiable
##                       (the pre-specified rule of config.R)
##   C. both quantiles : tau = 0.25 and 0.5 identifiable (KM reaches tau + margin,
##                       >= min_at_risk at risk at the KM tau-quantile)
##   D. screening cap  : floor(d2 / log d2) >= low_cap (25), d2 = events / 2, so
##                       that U+ at q = 0.1 (>= 10 selections) is attainable within
##                       the candidate set; also reports the cap relative to p
## and, per task table, which analyses are finished (main, mds, msms 1/2, mds_dai).
## Writes tables/requirements.csv and prints PASS / FAIL per check.

source(file.path(Sys.getenv("CQR_RD_ROOT", "."), "R", "config.R"))
suppressPackageStartupMessages(library(survival))
dir.create(file.path(RD_ROOT, "tables"), showWarnings = FALSE)

tau_ok <- function(time, status, tau) {
  km <- survfit(Surv(time, status) ~ 1); Fh <- 1 - km$surv
  if (max(Fh) < tau + CFG$tau_margin) return(FALSE)
  sum(time >= km$time[which(Fh >= tau)[1L]]) >= CFG$min_at_risk
}
capf <- function(ev) { d2 <- floor(ev / 2); if (d2 > 2) floor(d2 / log(d2)) else 0 }

## finished analyses, from all task tables
tfiles <- list.files(RD_ROOT, "^tasks.*\\.csv$", full.names = TRUE)
tk <- if (length(tfiles)) do.call(rbind, lapply(tfiles, function(f) {
  z <- utils::read.csv(f, stringsAsFactors = FALSE)[, c("dataset", "tau", "type", "rep")]
  z$file <- basename(f); z })) else NULL
done <- function(ds, type, rep = 1L) {
  if (is.null(tk)) return(NA_character_)
  z <- unique(tk[tk$dataset == ds & tk$type == type & tk$rep == rep, c("tau")])
  if (!length(z)) return("no task")
  ok <- file.exists(file.path(RD_ROOT, "results", ds, sprintf("tau%.2f", z), sprintf("%s_%03d.rds", type, rep)))
  sprintf("%d/%d", sum(ok), length(ok))
}

files <- list.files(file.path(RD_ROOT, "data", "clean"), "\\.rds$", full.names = TRUE)
if (!length(files)) stop("no data/clean/*.rds found")
rows <- list()
for (f in files) {
  D <- readRDS(f); x <- D$x; nm <- if (!is.null(D$name)) D$name else sub("\\.rds$", "", basename(f))
  n <- NROW(x); p <- NCOL(x); ev <- sum(D$delta == 1, na.rm = TRUE)
  ## A. integrity
  prob <- character()
  if (!is.matrix(x) || !is.numeric(x)) prob <- c(prob, "x not a numeric matrix")
  if (any(!is.finite(x))) prob <- c(prob, sprintf("%d non-finite x", sum(!is.finite(x))))
  sds <- apply(x, 2L, stats::sd)
  if (any(!is.finite(sds) | sds < 1e-8)) prob <- c(prob, sprintf("%d constant columns", sum(!is.finite(sds) | sds < 1e-8)))
  if (is.null(colnames(x)) || anyDuplicated(colnames(x))) prob <- c(prob, "missing/duplicated column names")
  if (p > 1 && p <= 5000) {
    cc <- abs(stats::cor(x)); diag(cc) <- 0
    if (any(cc > 1 - 1e-8, na.rm = TRUE)) prob <- c(prob, sprintf("%d identical column pairs", sum(cc[upper.tri(cc)] > 1 - 1e-8, na.rm = TRUE)))
  }
  if (length(D$time) != n || length(D$delta) != n || length(D$y) != n) prob <- c(prob, "length mismatch")
  if (any(!is.finite(D$time) | D$time <= 0)) prob <- c(prob, "time not positive/finite")
  if (!all(D$delta %in% 0:1)) prob <- c(prob, "delta not 0/1")
  if (max(abs(D$y - log(D$time))) > 1e-8) prob <- c(prob, "y != log(time)")
  ## B-D
  okt <- vapply(CFG$taus, function(t) tau_ok(D$time, D$delta, t), TRUE)
  cap <- capf(ev)
  rows[[length(rows) + 1L]] <- data.frame(
    dataset = nm, n = n, p = p, events = ev, censoring = round(1 - ev / n, 3),
    d2 = floor(ev / 2), cap = cap, cap_over_p = round(cap / p, 3),
    A_integrity = if (length(prob)) "FAIL" else "PASS", integrity_note = paste(prob, collapse = "; "),
    B_eligible = if (n >= CFG$min_n && ev >= CFG$min_events && okt[1L]) "PASS" else "FAIL",
    C_tau025 = okt[1L], C_tau050 = okt[2L], C_both = if (all(okt)) "PASS" else "FAIL",
    D_cap25 = if (cap >= CFG$low_cap) "PASS" else "FAIL",
    all_pass = length(prob) == 0 && n >= CFG$min_n && ev >= CFG$min_events && all(okt) && cap >= CFG$low_cap,
    main = done(nm, "main"), mds = done(nm, "mds"), msms_1 = done(nm, "msms", 1L),
    msms_2 = done(nm, "msms", 2L), mds_dai = done(nm, "mds_dai"),
    stringsAsFactors = FALSE)
}
R <- do.call(rbind, rows)
R <- R[order(!R$all_pass, -R$events), ]
utils::write.csv(R, file.path(RD_ROOT, "tables", "requirements.csv"), row.names = FALSE)
options(width = 220)
print(R[, c("dataset", "n", "p", "events", "cap", "A_integrity", "B_eligible", "C_both", "D_cap25",
            "all_pass", "main", "mds", "msms_1", "msms_2", "mds_dai")], row.names = FALSE)
bad <- R[R$A_integrity == "FAIL", ]
if (nrow(bad)) { cat("\nIntegrity problems:\n"); print(bad[, c("dataset", "integrity_note")], row.names = FALSE) }
cat(sprintf("\n%d of %d data sets meet all requirements (A-D).\n", sum(R$all_pass), nrow(R)))
