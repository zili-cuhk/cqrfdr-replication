## ---------------------------------------------------------------------------
## Targeted METABRIC analysis: published prognostic genes + clinical factors
## ---------------------------------------------------------------------------
## Usage: Rscript R/05_targeted.R
## Needs data/raw/BRCA_METABRIC.rds (R/01_download.R) and the METABRIC clinical
## files; data_clinical_sample.txt is downloaded here if it is missing.
## Writes data/clean/BRCA_METABRIC_TARGETED.rds, tasks_targeted.csv and
## tables_targeted/targeted_design.csv.
##
## The candidate set is fixed in advance from published signatures, not from
## the METABRIC outcome:
##   * the 50 PAM50 genes (Parker et al., 2009, J Clin Oncol);
##   * the 16 cancer-related genes of the Oncotype DX recurrence score
##     (Paik et al., 2004, N Engl J Med), 11 of which are new;
##   * ten routine clinical factors: age, tumour size, positive lymph nodes,
##     grade, ER, PR and HER2 status, chemotherapy, hormone therapy, radiotherapy.
## The question: which of these have an independent conditional effect on the
## 0.25 and 0.5 quantiles of overall survival time, given all the others?

source(file.path(Sys.getenv("CQR_RD_ROOT", "."), "R", "config.R"))
suppressPackageStartupMessages(library(survival))
dirs()
NAME <- "BRCA_METABRIC_TARGETED"
TBT <- file.path(RD_ROOT, "tables_targeted"); dir.create(TBT, showWarnings = FALSE)

PAM50 <- c("ACTR3B", "ANLN", "BAG1", "BCL2", "BIRC5", "BLVRA", "CCNB1", "CCNE1", "CDC20",
           "CDC6", "CDH3", "CENPF", "CEP55", "CXXC5", "EGFR", "ERBB2", "ESR1", "EXO1",
           "FGFR4", "FOXA1", "FOXC1", "GPR160", "GRB7", "KIF2C", "KRT14", "KRT17", "KRT5",
           "MAPT", "MDM2", "MELK", "MIA", "MKI67", "MLPH", "MMP11", "MYBL2", "MYC", "NAT1",
           "NDC80", "NUF2", "ORC6", "PGR", "PHGDH", "PTTG1", "RRM2", "SFRP1", "SLC39A6",
           "TMEM45B", "TYMS", "UBE2C", "UBE2T")
ONCOTYPE <- c("MKI67", "AURKA", "BIRC5", "CCNB1", "MYBL2", "GRB7", "ERBB2", "ESR1", "PGR",
              "BCL2", "SCUBE2", "MMP11", "CTSV", "GSTM1", "CD68", "BAG1")
ALIAS <- list(ORC6 = "ORC6L", CTSV = "CTSL2", AURKA = "STK15", NDC80 = "KNTC2",
              ACTR3B = "ACTR3BP", CENPF = "CENPF")
GENES <- unique(c(PAM50, ONCOTYPE))

r <- readRDS(file.path(RD_ROOT, "data", "raw", "BRCA_METABRIC.rds"))
X <- r$x
pick <- vapply(GENES, function(g) {
  cand <- c(g, ALIAS[[g]])
  hit <- cand[cand %in% colnames(X)]
  if (length(hit)) hit[1] else NA_character_
}, "")
missing_genes <- GENES[is.na(pick)]
G <- X[, pick[!is.na(pick)], drop = FALSE]
colnames(G) <- names(pick)[!is.na(pick)]
for (j in which(colSums(is.na(G)) > 0)) G[is.na(G[, j]), j] <- stats::median(G[, j], na.rm = TRUE)

## clinical data ------------------------------------------------------------
dl <- file.path(RD_ROOT, "data", "download", "brca_metabric")
fs <- list.files(dl, recursive = TRUE, full.names = TRUE)
f_pat <- fs[basename(fs) == "data_clinical_patient.txt"][1]
f_smp <- fs[basename(fs) == "data_clinical_sample.txt"][1]
if (is.na(f_smp)) {
  f_smp <- file.path(dl, "data_clinical_sample.txt")
  utils::download.file(paste0("https://media.githubusercontent.com/media/cBioPortal/datahub/master/",
                              "public/brca_metabric/data_clinical_sample.txt"), f_smp, mode = "wb")
}
pat <- utils::read.delim(f_pat, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
smp <- utils::read.delim(f_smp, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
cl <- merge(pat, smp, by = "PATIENT_ID", all.x = TRUE, suffixes = c("", ".smp"))
col <- function(nm) { for (n in nm) if (n %in% names(cl)) return(cl[[n]]); rep(NA, nrow(cl)) }
yes <- function(v) ifelse(is.na(v) | v == "", NA, as.integer(toupper(v) %in% c("YES", "POSITIVE", "POSITVE", "+")))
num <- function(v) suppressWarnings(as.numeric(v))
C <- data.frame(
  age        = num(col("AGE_AT_DIAGNOSIS")),
  log_size   = log(num(col("TUMOR_SIZE"))),
  log1p_nodes = log1p(num(col("LYMPH_NODES_EXAMINED_POSITIVE"))),
  grade      = num(col("GRADE")),
  ER_pos     = yes(col(c("ER_STATUS", "ER_IHC"))),
  PR_pos     = yes(col("PR_STATUS")),
  HER2_pos   = yes(col(c("HER2_STATUS", "HER2_SNP6"))),
  chemo      = yes(col("CHEMOTHERAPY")),
  hormone    = yes(col("HORMONE_THERAPY")),
  radio      = yes(col("RADIO_THERAPY")),
  row.names  = cl$PATIENT_ID)
## HER2_SNP6 is GAIN/NEUTRAL/LOSS/UNDEF; use it only if HER2_STATUS is absent
if (!"HER2_STATUS" %in% names(cl) && "HER2_SNP6" %in% names(cl))
  C$HER2_pos <- ifelse(cl$HER2_SNP6 %in% c("GAIN", "NEUTRAL", "LOSS"), as.integer(cl$HER2_SNP6 == "GAIN"), NA)
bad <- vapply(C, function(v) all(is.na(v)) || stats::sd(v, na.rm = TRUE) == 0, TRUE)
if (any(bad)) { cat("dropping clinical variables not found or constant:", names(C)[bad], "\n"); C <- C[, !bad, drop = FALSE] }

## align and keep complete cases --------------------------------------------
ids <- intersect(rownames(G), rownames(C))
keep <- ids[stats::complete.cases(C[ids, , drop = FALSE])]
m <- match(keep, rownames(G))
x <- cbind(G[m, , drop = FALSE], as.matrix(C[keep, , drop = FALSE]))
time <- r$time[m]; delta <- r$status[m]
tau_ok <- function(time, status, tau) {
  km <- survfit(Surv(time, status) ~ 1); Fh <- 1 - km$surv
  if (max(Fh) < tau + CFG$tau_margin) return(FALSE)
  sum(time >= km$time[which(Fh >= tau)[1L]]) >= CFG$min_at_risk
}
taus <- CFG$taus[vapply(CFG$taus, function(t) tau_ok(time, delta, t), TRUE)]
save_atomic(list(name = NAME, cancer = "BRCA", x = x, y = log(time), time = time, delta = delta,
                 taus = taus, genes = colnames(G), clinical = colnames(C)),
            file.path(RD_ROOT, "data", "clean", paste0(NAME, ".rds")))
d2 <- floor(sum(delta) / 2)
design <- data.frame(dataset = NAME, n = nrow(x), events = sum(delta),
                     censoring = round(1 - mean(delta), 3), genes = ncol(G), clinical = ncol(C),
                     p = ncol(x), screen_cap = floor(d2 / log(d2)),
                     taus = paste(taus, collapse = ";"),
                     genes_missing = paste(missing_genes, collapse = ";"),
                     dropped_incomplete = length(ids) - length(keep))
write.csv(design, file.path(TBT, "targeted_design.csv"), row.names = FALSE)
print(t(design))

## task table -----------------------------------------------------------------
tk <- do.call(rbind, lapply(taus, function(t) rbind(
  data.frame(dataset = NAME, tau = t, type = "main", rep = 1L),
  data.frame(dataset = NAME, tau = t, type = "mds", rep = 1L),
  data.frame(dataset = NAME, tau = t, type = "stab", rep = seq_len(CFG$n_stab)),
  data.frame(dataset = NAME, tau = t, type = "decoy", rep = seq_len(CFG$n_decoy)),
  data.frame(dataset = NAME, tau = t, type = "holdout", rep = seq_len(CFG$n_holdout)))))
tk <- tk[order(!tk$type %in% c("main", "mds"), tk$type, tk$tau), ]
tk$id <- seq_len(nrow(tk))
write.csv(tk, file.path(RD_ROOT, "tasks_targeted.csv"), row.names = FALSE)
cat(sprintf("\n%d tasks written to tasks_targeted.csv\n", nrow(tk)))
