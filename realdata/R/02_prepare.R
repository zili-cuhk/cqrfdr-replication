## ---------------------------------------------------------------------------
## Eligibility, outcome-blind gene filter, task table
## ---------------------------------------------------------------------------
## Usage: Rscript R/02_prepare.R
## Writes tables/eligibility.csv, data/clean/<name>.rds and tasks.csv.
## The gene filter uses only the expression matrix (MAD), never the outcome:
## a filter that looks at survival on the full data would let both halves of
## every split see the outcome and void the FDR guarantee.

source(file.path(Sys.getenv("CQR_RD_ROOT", "."), "R", "config.R"))
suppressPackageStartupMessages(library(survival))
dirs()
files <- list.files(file.path(RD_ROOT, "data", "raw"), "\\.rds$", full.names = TRUE)
if (!length(files)) stop("no raw data sets; run R/01_download.R first")

## KM-based identifiability of tau: the tau-quantile must be reached before the
## last event time with a margin, and enough subjects must remain at risk there.
tau_ok <- function(time, status, tau) {
  km <- survfit(Surv(time, status) ~ 1)
  Fhat <- 1 - km$surv
  if (max(Fhat) < tau + CFG$tau_margin) return(FALSE)
  t_tau <- km$time[which(Fhat >= tau)[1L]]
  sum(time >= t_tau) >= CFG$min_at_risk
}

rows <- list()
for (f in files) {
  r <- readRDS(f)
  x <- r$x; n <- nrow(x)
  ## genes: < 20% missing, impute by gene median, non-constant
  miss <- colMeans(is.na(x)); x <- x[, miss < 0.2, drop = FALSE]
  if (anyNA(x)) for (j in which(colSums(is.na(x)) > 0)) x[is.na(x[, j]), j] <- stats::median(x[, j], na.rm = TRUE)
  mad <- apply(x, 2L, stats::mad)
  x <- x[, mad > 1e-8, drop = FALSE]; mad <- mad[mad > 1e-8]
  keep <- order(-mad)[seq_len(min(CFG$p_filter, ncol(x)))]
  x <- x[, sort(keep), drop = FALSE]
  km <- survfit(Surv(r$time, r$status) ~ 1)
  med <- unname(summary(km)$table["median"])
  ok_tau <- vapply(CFG$taus, function(t) tau_ok(r$time, r$status, t), TRUE)
  ev <- sum(r$status)
  eligible <- n >= CFG$min_n && ev >= CFG$min_events && ok_tau[1L]
  vcand <- VALIDATION[[r$name]]
  vavail <- vcand[file.exists(file.path(RD_ROOT, "data", "raw", paste0(vcand, ".rds")))]
  rows[[length(rows) + 1L]] <- data.frame(
    dataset = r$name, cancer = r$cancer, source = r$source, platform = r$platform,
    endpoint = r$endpoint, n = n, events = ev, censoring = round(1 - ev / n, 3),
    median_fu_years = round(unname(summary(survfit(Surv(r$time, 1 - r$status) ~ 1))$table["median"]), 2),
    km_median_years = round(med, 2), genes_total = ncol(r$x), p = ncol(x),
    tau025 = ok_tau[1L], tau050 = ok_tau[2L], eligible = eligible,
    validation = paste(vavail, collapse = ";"), stringsAsFactors = FALSE)
  if (eligible)
    save_atomic(list(name = r$name, cancer = r$cancer, x = x, y = log(r$time), time = r$time,
                     delta = r$status, taus = CFG$taus[ok_tau]),
                file.path(RD_ROOT, "data", "clean", paste0(r$name, ".rds")))
}
el <- do.call(rbind, rows)

## pre-specified rules (config.R)
el$both_taus <- el$tau025 & el$tau050
el$has_validation <- nzchar(el$validation)
d2 <- floor(el$events / 2)
el$screen_cap <- ifelse(d2 > 2, floor(d2 / log(pmax(d2, 3))), 0)
el$analysed <- el$eligible & el$both_taus & el$screen_cap >= CFG$min_cap
el$low_cap <- el$screen_cap < CFG$low_cap
o <- order(!el$analysed, !el$eligible, !el$has_validation, -el$events)
el <- el[o, ]
el$rank <- ifelse(el$analysed, cumsum(el$analysed), NA)
el$role <- ifelse(!el$eligible, "not eligible",
             ifelse(!el$analysed, sprintf("listed only (cap < %d)", CFG$min_cap),
               ifelse(el$rank == 1L, "MAIN TEXT",
                 ifelse(el$low_cap, sprintf("supplement (cap < %d)", CFG$low_cap), "supplement"))))
## OV_TCGA (curatedOvarianData, Affymetrix) and TCGA_OV (Xena, RNA-seq) are the
## same patients: report both, never count them as independent cohorts
el$note <- ifelse(el$dataset %in% c("OV_TCGA", "TCGA_OV"),
                  "same TCGA-OV patients on two platforms", "")
write.csv(el, file.path(RD_ROOT, "tables", "eligibility.csv"), row.names = FALSE)
print(el[el$eligible, c("dataset", "n", "events", "censoring", "screen_cap", "tau025", "tau050",
                        "validation", "role")], row.names = FALSE)

## task table: one row per (data set, tau, type, replicate)
tk <- list()
for (i in which(el$analysed)) {
  d <- el$dataset[i]; taus <- CFG$taus[c(el$tau025[i], el$tau050[i])]
  for (t in taus) tk[[length(tk) + 1L]] <- rbind(
    data.frame(dataset = d, tau = t, type = "main", rep = 1L),
    data.frame(dataset = d, tau = t, type = "mds", rep = 1L),
    data.frame(dataset = d, tau = t, type = "stab", rep = seq_len(CFG$n_stab)),
    data.frame(dataset = d, tau = t, type = "decoy", rep = seq_len(CFG$n_decoy)),
    data.frame(dataset = d, tau = t, type = "holdout", rep = seq_len(CFG$n_holdout)))
}
tk <- do.call(rbind, tk)
## main and MDS tasks first (main-text data set first within them)
tk <- tk[order(!tk$type %in% c("main", "mds"), match(tk$dataset, el$dataset), tk$type, tk$tau), ]
tk$id <- seq_len(nrow(tk))
write.csv(tk, file.path(RD_ROOT, "tasks.csv"), row.names = FALSE)
cat(sprintf("\n%d eligible, %d analysed data sets; %d tasks written to tasks.csv\n",
            sum(el$eligible), sum(el$analysed), nrow(tk)))
