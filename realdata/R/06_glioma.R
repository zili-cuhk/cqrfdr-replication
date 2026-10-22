## ---------------------------------------------------------------------------
## Strong-signal benchmark: TCGA diffuse glioma (LGG + GBM, WHO grade II-IV)
## ---------------------------------------------------------------------------
## Usage: Rscript R/06_glioma.R
## Needs data/raw/TCGA_LGG.rds and data/raw/TCGA_GBM.rds (R/01_download.R tcga).
## The Xena clinical matrices (age, sex) are downloaded here if missing.
## Writes data/clean/GLIOMA_TCGA.rds, tasks_glioma.csv and
## tables_glioma/glioma_design.csv.
##
## Why this data set (fixed before any analysis of it): in diffuse glioma,
## survival is known to be driven by strong, well-replicated factors (grade,
## IDH mutation, 1p/19q co-deletion, age), with large expression differences
## between the molecular subtypes (Brat et al. 2015 NEJM; Ceccarelli et al.
## 2016 Cell).  It is the standard "strong signal" setting in TCGA and is
## used here to show the method when signals are strong, next to METABRIC,
## where they are weak.  LGG and GBM are combined as in Ceccarelli et al.;
## neither is eligible alone (LGG has too few deaths, GBM too few patients).
##
## Design, same as the other analyses:
##   * genes: top 2000 by MAD on the combined expression matrix (outcome-blind);
##   * clinical: age at diagnosis, male sex, GBM (grade IV) indicator;
##   * no outcome information is used to build the candidate set.

source(file.path(Sys.getenv("CQR_RD_ROOT", "."), "R", "config.R"))
suppressPackageStartupMessages(library(survival))
dirs()
NAME <- "GLIOMA_TCGA"
TBG <- file.path(RD_ROOT, "tables_glioma"); dir.create(TBG, showWarnings = FALSE)
RAW <- file.path(RD_ROOT, "data", "raw")
DL  <- file.path(RD_ROOT, "data", "download")
xena <- "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"

rd <- lapply(c(LGG = "TCGA_LGG", GBM = "TCGA_GBM"), function(nm) {
  f <- file.path(RAW, paste0(nm, ".rds"))
  if (!file.exists(f)) stop(f, " not found; run: Rscript R/01_download.R tcga")
  readRDS(f)
})
genes <- Reduce(intersect, lapply(rd, function(r) colnames(r$x)))
x <- do.call(rbind, lapply(rd, function(r) r$x[, genes, drop = FALSE]))
time <- unlist(lapply(rd, `[[`, "time"), use.names = FALSE)
delta <- unlist(lapply(rd, `[[`, "status"), use.names = FALSE)
gbm <- rep(c(0L, 1L), vapply(rd, function(r) nrow(r$x), 1L))
pat <- substr(rownames(x), 1L, 12L)
if (anyDuplicated(pat)) {                      # a patient in both cohorts: keep the first
  k <- !duplicated(pat); x <- x[k, ]; time <- time[k]; delta <- delta[k]; gbm <- gbm[k]; pat <- pat[k]
}

## clinical: age and sex from the Xena clinical matrices ----------------------
clin <- do.call(rbind, lapply(c("LGG", "GBM"), function(co) {
  f <- file.path(DL, paste0(co, "_clinicalMatrix.txt"))
  if (!file.exists(f))
    utils::download.file(paste0(xena, "TCGA.", co, ".sampleMap%2F", co, "_clinicalMatrix"), f, mode = "wb")
  d <- utils::read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  data.frame(patient = substr(d[[1L]], 1L, 12L),
             age = suppressWarnings(as.numeric(d$age_at_initial_pathologic_diagnosis)),
             male = ifelse(toupper(d$gender) %in% c("MALE", "FEMALE"), as.integer(toupper(d$gender) == "MALE"), NA))
}))
clin <- clin[!is.na(clin$age) & !duplicated(clin$patient), ]
m <- match(pat, clin$patient)
C <- cbind(age = clin$age[m], male = clin$male[m], GBM = gbm)
keep <- stats::complete.cases(C)
x <- x[keep, , drop = FALSE]; C <- C[keep, , drop = FALSE]
time <- time[keep]; delta <- delta[keep]

## outcome-blind gene filter, as in R/02_prepare.R ------------------------------
miss <- colMeans(is.na(x)); x <- x[, miss < 0.2, drop = FALSE]
if (anyNA(x)) for (j in which(colSums(is.na(x)) > 0)) x[is.na(x[, j]), j] <- stats::median(x[, j], na.rm = TRUE)
mad <- apply(x, 2L, stats::mad); x <- x[, mad > 1e-8, drop = FALSE]; mad <- mad[mad > 1e-8]
x <- x[, sort(order(-mad)[seq_len(min(CFG$p_filter, ncol(x)))]), drop = FALSE]
X <- cbind(x, C)

tau_ok <- function(time, status, tau) {
  km <- survfit(Surv(time, status) ~ 1); Fh <- 1 - km$surv
  if (max(Fh) < tau + CFG$tau_margin) return(FALSE)
  sum(time >= km$time[which(Fh >= tau)[1L]]) >= CFG$min_at_risk
}
taus <- CFG$taus[vapply(CFG$taus, function(t) tau_ok(time, delta, t), TRUE)]
if (!length(taus)) stop("neither tau is identifiable")
save_atomic(list(name = NAME, cancer = "GLIOMA", x = X, y = log(time), time = time, delta = delta,
                 taus = taus, genes = colnames(x), clinical = colnames(C)),
            file.path(RD_ROOT, "data", "clean", paste0(NAME, ".rds")))

d2 <- floor(sum(delta) / 2)
km <- survfit(Surv(time, delta) ~ 1)
design <- data.frame(dataset = NAME, n = nrow(X), n_LGG = sum(C[, "GBM"] == 0), n_GBM = sum(C[, "GBM"] == 1),
                     events = sum(delta), censoring = round(1 - mean(delta), 3),
                     km_median_years = round(unname(summary(km)$table["median"]), 2),
                     genes = ncol(x), clinical = ncol(C), p = ncol(X),
                     screen_cap = floor(d2 / log(d2)), taus = paste(taus, collapse = ";"),
                     dropped_no_clinical = sum(!keep))
write.csv(design, file.path(TBG, "glioma_design.csv"), row.names = FALSE)
print(t(design))
if (floor(d2 / log(d2)) < CFG$low_cap)
  cat("\nNOTE: screening cap below", CFG$low_cap, "- little power expected at q = 0.1\n")

tk <- do.call(rbind, lapply(taus, function(t) rbind(
  data.frame(dataset = NAME, tau = t, type = "main", rep = 1L),
  data.frame(dataset = NAME, tau = t, type = "mds", rep = 1L),
  data.frame(dataset = NAME, tau = t, type = "stab", rep = seq_len(CFG$n_stab)),
  data.frame(dataset = NAME, tau = t, type = "decoy", rep = seq_len(CFG$n_decoy)),
  data.frame(dataset = NAME, tau = t, type = "holdout", rep = seq_len(CFG$n_holdout)))))
tk <- tk[order(!tk$type %in% c("main", "mds"), tk$type, tk$tau), ]
tk$id <- seq_len(nrow(tk))
write.csv(tk, file.path(RD_ROOT, "tasks_glioma.csv"), row.names = FALSE)
cat(sprintf("\n%d tasks written to tasks_glioma.csv\n", nrow(tk)))
