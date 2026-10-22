## ---------------------------------------------------------------------------
## Download all candidate data sets and store them in one common format
## ---------------------------------------------------------------------------
## Usage (on the LOGIN node, which has internet access):
##   Rscript R/01_download.R            # all sources
##   Rscript R/01_download.R ov tcga    # only some sources: ov, metabric, dlbcl, tcga
## Each data set is saved as data/raw/<name>.rds, a list with
##   name, source, cancer, platform,
##   x      (n x G expression matrix, columns = gene symbols, one row per patient),
##   time   (years), status (1 = event), endpoint, note.
## A source that fails is reported and skipped; re-running only fetches what is missing.

source(file.path(Sys.getenv("CQR_RD_ROOT", "."), "R", "config.R"))
dirs()
RAW <- file.path(RD_ROOT, "data", "raw")
DL  <- file.path(RD_ROOT, "data", "download"); dir.create(DL, showWarnings = FALSE, recursive = TRUE)
options(timeout = max(3600, getOption("timeout")))
args <- tolower(commandArgs(TRUE))
want <- function(s) !length(args) || s %in% args
log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

## collapse probes to genes: keep, for each symbol, the probe with largest IQR
collapse_genes <- function(E, sym) {           # E: probes x samples
  sym <- trimws(sub("///.*$", "", as.character(sym)))
  ok <- !is.na(sym) & nzchar(sym) & sym != "---"
  E <- E[ok, , drop = FALSE]; sym <- sym[ok]
  iqr <- apply(E, 1L, stats::IQR, na.rm = TRUE)
  o <- order(sym, -iqr); keep <- o[!duplicated(sym[o])]
  out <- E[keep, , drop = FALSE]; rownames(out) <- sym[keep]
  out
}
to_log2 <- function(E) if (max(E, na.rm = TRUE) > 50) log2(pmax(E, 1)) else E

save_raw <- function(name, x, time, status, source, cancer, platform, endpoint = "OS", note = "") {
  ok <- is.finite(time) & time > 0 & status %in% 0:1
  if (sum(ok) < 50) stop(sprintf("only %d usable patients (check sample matching)", sum(ok)))
  x <- x[ok, , drop = FALSE]
  obj <- list(name = name, source = source, cancer = cancer, platform = platform,
              x = x, time = time[ok], status = as.integer(status[ok]),
              endpoint = endpoint, note = note)
  save_atomic(obj, file.path(RAW, paste0(name, ".rds")))
  log_msg("saved %-24s n = %4d, events = %4d, genes = %5d", name, nrow(x), sum(obj$status), ncol(x))
}
have <- function(name) file.exists(file.path(RAW, paste0(name, ".rds")))
attempt <- function(label, expr) {
  r <- tryCatch({ force(expr); TRUE }, error = function(e) {
    log_msg("FAILED %s: %s", label, conditionMessage(e)); FALSE })
  invisible(r)
}

## 1. Ovarian cancer: curatedOvarianData (Bioconductor) ---------------------------
## Every expression set with overall survival.  TCGA RNA-seq / miRNA sets are
## skipped (same patients as TCGA_eset).
if (want("ov")) attempt("curatedOvarianData", {
  suppressPackageStartupMessages(library(curatedOvarianData))
  sets <- data(package = "curatedOvarianData")$results[, "Item"]
  sets <- sets[grepl("_eset$", sets) & !grepl("mirna|RNASeq", sets, ignore.case = TRUE)]
  for (s in sets) {
    nm <- paste0("OV_", sub("_eset$", "", s))
    if (have(nm)) next
    attempt(nm, {
      e <- new.env(); data(list = s, package = "curatedOvarianData", envir = e)
      es <- get(s, envir = e); pd <- Biobase::pData(es)
      if (!all(c("vital_status", "days_to_death") %in% names(pd))) stop("no survival columns")
      keep <- !is.na(pd$vital_status) & !is.na(pd$days_to_death)
      if ("sample_type" %in% names(pd)) keep <- keep & pd$sample_type %in% c("tumor", NA)
      if (sum(keep) < 50) stop("fewer than 50 tumours with survival")
      E <- Biobase::exprs(es)[, keep, drop = FALSE]
      ## duplicated patients (if a patient id column exists)
      pid <- if ("unique_patient_ID" %in% names(pd)) pd$unique_patient_ID[keep] else colnames(E)
      d <- !duplicated(pid) | is.na(pid); E <- E[, d, drop = FALSE]
      st <- as.integer(pd$vital_status[keep][d] == "deceased")
      tm <- as.numeric(pd$days_to_death[keep][d]) / 365.25
      save_raw(nm, t(to_log2(E)), tm, st, "curatedOvarianData", "OV",
               Biobase::annotation(es), note = s)
    })
  }
})

## 2. Breast cancer: METABRIC from the cBioPortal datahub ---------------------------
if (want("metabric") && !have("BRCA_METABRIC")) attempt("METABRIC", {
  ## sources, in order: an already extracted folder or an archive placed in
  ## data/download by hand (e.g. the zip from https://www.cbioportal.org/datasets),
  ## the S3 tarball, the individual files of the GitHub datahub (Git LFS)
  ex <- file.path(DL, "brca_metabric")
  tgz <- file.path(DL, "brca_metabric.tar.gz"); zp <- file.path(DL, "brca_metabric.zip")
  if (!dir.exists(ex) && file.exists(zp)) utils::unzip(zp, exdir = DL)
  if (!dir.exists(ex) && !file.exists(tgz))
    tryCatch(utils::download.file("https://cbioportal-datahub.s3.amazonaws.com/brca_metabric.tar.gz",
                                  tgz, mode = "wb"),
             error = function(e) { unlink(tgz); log_msg("S3 tarball not available: %s", conditionMessage(e)) },
             warning = function(w) { unlink(tgz); log_msg("S3 tarball not available: %s", conditionMessage(w)) })
  if (!dir.exists(ex) && file.exists(tgz)) utils::untar(tgz, exdir = DL)
  if (!dir.exists(ex)) {
    dir.create(ex)
    gh <- "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_metabric/"
    for (f in c("data_clinical_patient.txt", "data_mrna_illumina_microarray.txt", "data_expression_median.txt")) {
      dest <- file.path(ex, f)
      ok <- tryCatch({ utils::download.file(paste0(gh, f), dest, mode = "wb"); TRUE },
                     error = function(e) FALSE, warning = function(w) FALSE)
      if (!ok || file.size(dest) < 10000 || grepl("^version https://git-lfs", readLines(dest, 1L, warn = FALSE))) {
        unlink(dest); log_msg("GitHub datahub: %s not available", f)
      } else log_msg("GitHub datahub: got %s", f)
    }
  }
  fs <- list.files(ex, full.names = TRUE, recursive = TRUE)
  fcl <- fs[basename(fs) == "data_clinical_patient.txt"]
  fex <- fs[grepl("^data_(mrna|expression).*\\.txt$", basename(fs)) &
            !grepl("zscore|Zscore|seq", basename(fs))][1L]
  if (!length(fcl) || is.na(fex)) stop("clinical or expression file not found in ", ex)
  cl <- utils::read.delim(fcl, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
  ed <- data.table::fread(fex, data.table = FALSE)
  sym <- ed$Hugo_Symbol
  E <- as.matrix(ed[, setdiff(names(ed), c("Hugo_Symbol", "Entrez_Gene_Id")), drop = FALSE])
  E <- collapse_genes(E, sym)
  st <- as.integer(substr(cl$OS_STATUS, 1L, 1L) == "1")
  tm <- as.numeric(cl$OS_MONTHS) / 12
  m <- match(colnames(E), cl$PATIENT_ID); ok <- !is.na(m) & !is.na(st[m]) & !is.na(tm[m])
  save_raw("BRCA_METABRIC", t(E[, ok]), tm[m[ok]], st[m[ok]],
           "cBioPortal brca_metabric", "BRCA", basename(fex))
})

## 3. DLBCL: GSE10846 (Lenz et al., 2008; CHOP and R-CHOP) ---------------------------
if (want("dlbcl") && !have("DLBCL_GSE10846")) attempt("GSE10846", {
  suppressPackageStartupMessages(library(GEOquery))
  g <- GEOquery::getGEO("GSE10846", GSEMatrix = TRUE, getGPL = TRUE, destdir = DL)[[1L]]
  pd <- Biobase::pData(g)
  ch <- apply(pd[, grepl("^characteristics", names(pd)), drop = FALSE], 1L, paste, collapse = " | ")
  grab <- function(pat) sub(paste0(".*", pat, "\\s*([^|]+).*"), "\\1", ifelse(grepl(pat, ch, ignore.case = TRUE), ch, NA))
  st_raw <- toupper(trimws(grab("Follow up status:")))
  tm <- suppressWarnings(as.numeric(trimws(grab("Follow up years:"))))
  st <- ifelse(st_raw == "DEAD", 1L, ifelse(st_raw == "ALIVE", 0L, NA))
  fd <- Biobase::fData(g)
  scol <- grep("^Gene.Symbol$|^Gene Symbol$|^gene_assignment$", names(fd), value = TRUE)[1L]
  E <- collapse_genes(to_log2(Biobase::exprs(g)), fd[[scol]])
  ok <- !is.na(st) & !is.na(tm)
  save_raw("DLBCL_GSE10846", t(E[, ok]), tm[ok], st[ok], "GEO GSE10846", "DLBCL",
           Biobase::annotation(g), note = "CHOP and R-CHOP combined")
})

## 4. TCGA cohorts from UCSC Xena (HiSeqV2 = log2(RSEM + 1); TCGA-CDR survival) -----------
xena <- "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
if (want("tcga")) for (co in TCGA_COHORTS) {
  nm <- paste0("TCGA_", co)
  if (have(nm)) next
  attempt(nm, {
    fe <- file.path(DL, paste0("TCGA.", co, ".HiSeqV2.gz"))
    fs <- file.path(DL, paste0(co, "_survival.txt"))
    if (!file.exists(fe)) utils::download.file(paste0(xena, "TCGA.", co, ".sampleMap%2FHiSeqV2.gz"), fe, mode = "wb")
    if (!file.exists(fs)) utils::download.file(paste0(xena, "survival%2F", co, "_survival.txt"), fs, mode = "wb")
    ed <- data.table::fread(fe, data.table = FALSE)
    E <- as.matrix(ed[, -1L]); rownames(E) <- ed[[1L]]
    sv <- utils::read.delim(fs, check.names = FALSE, stringsAsFactors = FALSE)
    ep <- CFG$endpoint
    ## one tumour sample per patient: primary (01) > primary blood (03) > metastatic (06)
    code <- substr(colnames(E), 14L, 15L); pat <- substr(colnames(E), 1L, 12L)
    pr <- match(code, c("01", "03", "06")); keep <- !is.na(pr)
    o <- order(pat[keep], pr[keep]); idx <- which(keep)[o]
    idx <- idx[!duplicated(pat[idx])]
    ## match on the patient barcode (survival is patient-level); fall back to the sample barcode
    m <- if ("_PATIENT" %in% names(sv)) match(pat[idx], sv[["_PATIENT"]]) else match(colnames(E)[idx], sv$sample)
    if (all(is.na(m))) m <- match(substr(colnames(E)[idx], 1L, 15L), substr(sv$sample, 1L, 15L))
    st <- suppressWarnings(as.integer(sv[[ep]][m])); tm <- suppressWarnings(as.numeric(sv[[paste0(ep, ".time")]][m])) / 365.25
    ok <- !is.na(m) & !is.na(st) & !is.na(tm)
    save_raw(nm, t(E[, idx[ok]]), tm[ok], st[ok], "UCSC Xena TCGA hub", co, "IlluminaHiSeq RNASeqV2",
             endpoint = ep)
  })
}
log_msg("done; raw data sets: %d", length(list.files(RAW, "\\.rds$")))
