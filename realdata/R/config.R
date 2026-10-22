## ---------------------------------------------------------------------------
## Real-data analysis: configuration and pre-specified rules
## ---------------------------------------------------------------------------
## Everything that decides WHICH data set goes to the main text is fixed here,
## before any selection result is seen.  The rule uses only properties of the
## data (events, identifiable quantiles, external validation), never the
## behaviour of SS-MS or of a competitor.

RD_ROOT <- normalizePath(Sys.getenv("CQR_RD_ROOT", "."), mustWork = FALSE)
RD_TEST <- nzchar(Sys.getenv("CQR_RD_TEST"))     # small settings for a dry run

CFG <- list(
  endpoint     = Sys.getenv("CQR_RD_ENDPOINT", "OS"),  # TCGA endpoint (OS or PFI)
  p_filter     = 2000L,         # genes kept by an outcome-blind MAD filter
  min_events   = 150L,          # eligibility: at least this many events
  min_n        = 200L,
  min_cap      = 0L,            # analysed only if floor(d2 / log d2) >= this, d2 = events / 2
                                # (0: every eligible cohort with both taus identifiable is analysed)
  low_cap      = 25L,           # cohorts with a smaller cap are flagged as low-power: a selection
                                # at q = 0.1 needs >= 10 variables out of fewer than 25 candidates
  taus         = c(0.25, 0.5),  # target quantiles
  tau_margin   = 0.05,          # tau identifiable if KM F(last event time) >= tau + margin
  min_at_risk  = 20L,           # ... and at least this many at risk at the KM tau-quantile
  nu           = 0.1, n_grid = 40L, q = c(0.1, 0.2), nfolds = 5L,
  n_stab       = 100L,          # random splits for split-to-split stability
  n_decoy      = 20L,           # decoy (permuted-gene) replicates
  n_decoy_vars = 200L,          # decoys appended per replicate
  n_holdout    = 20L,           # random 70/30 train/test splits
  holdout_frac = 0.7,
  mds_B        = 50L,
  cores        = as.integer(Sys.getenv("CQR_RD_CORES", "4"))
)
if (RD_TEST) CFG[c("p_filter", "min_events", "min_n", "n_grid", "nfolds", "n_stab",
                   "n_decoy", "n_decoy_vars", "n_holdout", "mds_B", "cores")] <-
  list(150L, 60L, 100L, 10L, 3L, 3L, 2L, 20L, 2L, 3L, 1L)
if (RD_TEST) CFG$low_cap <- 20L

## TCGA cohorts pulled from UCSC Xena (TCGA-CDR survival).  OV and BRCA are
## included too: TCGA-OV via Xena RNA-seq is a second OV data set, TCGA-BRCA
## serves as the external cohort for METABRIC.
TCGA_COHORTS <- c("BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC",
                  "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD",
                  "SARC", "SKCM", "STAD", "UCEC")

## External validation pairs: discovery -> independent cohorts of the same
## cancer from a different study.  Fixed before the analysis.
VALIDATION <- list(
  OV_TCGA       = c("OV_GSE9891", "OV_GSE32062.GPL6480", "OV_GSE26712"),
  BRCA_METABRIC = c("TCGA_BRCA"),
  OV_GSE9891    = c("OV_TCGA", "OV_GSE32062.GPL6480"),
  TCGA_OV       = c("OV_GSE9891", "OV_GSE32062.GPL6480")
)

## Pre-specified rules (fixed before any selection result was seen):
##  eligible : n >= min_n, events >= min_events, tau = 0.25 identifiable;
##  analysed : eligible and both taus identifiable (the first version also
##             required floor(d2 / log d2) >= 25, d2 = events / 2; it was dropped
##             so that every eligible cohort is analysed; such cohorts are flagged);
##  main text: among the analysed, external validation available, then most events.
## Eligible but not analysed data sets are listed (n, events, cap) in the supplement.

## task table and output tables (the targeted analysis uses its own pair)
TASKS_FILE <- file.path(RD_ROOT, Sys.getenv("CQR_RD_TASKS", "tasks.csv"))
TABLES_DIR <- file.path(RD_ROOT, Sys.getenv("CQR_RD_TABLES", "tables"))
VALIDATION$BRCA_METABRIC_TARGETED <- "TCGA_BRCA"

dirs <- function() {
  d <- file.path(RD_ROOT, c("data/raw", "data/clean", "results", "tables", "logs"))
  for (x in d) dir.create(x, recursive = TRUE, showWarnings = FALSE)
  invisible(d)
}

## deterministic seed from strings
str_seed <- function(...) {
  h <- 0
  for (v in utf8ToInt(paste(..., sep = "|"))) h <- (h * 131 + v) %% 2147483629
  as.integer(h %% 2147483000 + 1)
}

save_atomic <- function(obj, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dirname(path), fileext = ".tmp")
  saveRDS(obj, tmp)
  if (!file.rename(tmp, path)) stop("could not move ", tmp, " to ", path)
  invisible(path)
}
