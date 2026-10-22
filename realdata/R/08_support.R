## ---------------------------------------------------------------------------
## SUPPORT (support2): many events, moderately correlated clinical covariates
## ---------------------------------------------------------------------------
## Usage: Rscript R/08_support.R          (login node: downloads the data)
## Writes data/clean/SUPPORT2.rds, tasks_support.csv and
## tables_support/support_design.csv.
##
## Why this data set (fixed from the plasmode/ablation results, before any
## analysis of it): SS-MS needs a screen that covers the true signals and at
## least about 1/q of them.  The ablation showed both fail in gene-expression
## data with a few hundred events and strongly correlated genes.  SUPPORT has
## 9105 seriously ill hospitalised adults, most of whom died during follow-up
## (up to 5.6 years), and about 30 baseline covariates with moderate
## correlation and many established independent prognostic effects (Knaus et
## al. 1995, Ann Intern Med).  The question: which main effects, curvatures
## and pairwise interactions act on the 0.25 and 0.5 quantiles of survival?
##
## Candidate set, built without the outcome:
##   * baseline covariates (no post-baseline or outcome-derived variables:
##     sps, aps, surv2m/6m, prg2m/6m, dnr, avtisst, charges, slos excluded);
##   * missing laboratory values imputed with Harrell's "normal" values, plus a
##     missing-value indicator for each laboratory variable with >5% missing;
##   * skewed laboratory values on the log scale;
##   * squares of the centred physiology variables (U-shaped risks);
##   * all pairwise products of the centred base covariates.

source(file.path(Sys.getenv("CQR_RD_ROOT", "."), "R", "config.R"))
suppressPackageStartupMessages(library(survival))
dirs()
NAME <- "SUPPORT2"
TBS <- file.path(RD_ROOT, "tables_support"); dir.create(TBS, showWarnings = FALSE)
DL <- file.path(RD_ROOT, "data", "download", "support2"); dir.create(DL, recursive = TRUE, showWarnings = FALSE)

## download ----------------------------------------------------------------------
get_support <- function() {
  csv <- list.files(DL, "support2.*\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(csv)) return(utils::read.csv(csv[1], stringsAsFactors = FALSE))
  tries <- list(
    list(url = "https://archive.ics.uci.edu/static/public/880/support2.zip", zip = TRUE),
    list(url = "https://hbiostat.org/data/repo/support2csv.zip", zip = TRUE))
  for (t in tries) {
    f <- file.path(DL, basename(t$url))
    ok <- tryCatch({ utils::download.file(t$url, f, mode = "wb", quiet = TRUE); TRUE }, error = function(e) FALSE)
    if (ok && file.size(f) > 1000) {
      utils::unzip(f, exdir = DL)
      csv <- list.files(DL, "\\.csv$", recursive = TRUE, full.names = TRUE)
      if (length(csv)) { cat("downloaded", t$url, "\n"); return(utils::read.csv(csv[1], stringsAsFactors = FALSE)) }
    }
  }
  if (requireNamespace("Hmisc", quietly = TRUE)) {
    e <- new.env(); Hmisc::getHdata(support2, where = e)
    return(as.data.frame(lapply(e$support2, function(v) if (is.factor(v)) as.character(v) else as.vector(v))))
  }
  stop("could not download support2; put support2.csv into ", DL)
}
S <- get_support()
names(S) <- gsub("_", ".", tolower(names(S)))
cat("raw data:", nrow(S), "rows,", ncol(S), "columns\n")

num <- function(v) suppressWarnings(as.numeric(v))
low <- function(v) tolower(trimws(as.character(v)))
time <- num(S$d.time); delta <- as.integer(num(S$death))
ok <- is.finite(time) & time > 0 & delta %in% 0:1
S <- S[ok, ]; time <- time[ok]; delta <- delta[ok]

## base covariates -----------------------------------------------------------------
dz <- low(S$dzgroup)
dz_lev <- c(copd = "copd", chf = "chf", cirrhosis = "cirrhosis", coma = "coma",
            colon_ca = "colon", lung_ca = "lung", mosf_malig = "malig")      # reference: ARF/MOSF w/ sepsis
D <- data.frame(age = num(S$age), male = as.integer(low(S$sex) == "male"))
for (k in names(dz_lev)) D[[paste0("dz_", k)]] <- as.integer(grepl(dz_lev[[k]], dz))
D$num_co   <- num(S$num.co)
D$edu      <- num(S$edu)
D$scoma    <- num(S$scoma)
D$race_black <- as.integer(low(S$race) == "black")
D$race_other <- as.integer(low(S$race) %in% c("asian", "hispanic", "other"))
D$diabetes <- num(S$diabetes)
D$dementia <- num(S$dementia)
D$ca_yes   <- as.integer(low(S$ca) == "yes")
D$ca_meta  <- as.integer(low(S$ca) == "metastatic")
D$hday     <- log1p(num(S$hday))
adl <- if ("adlsc" %in% names(S)) num(S$adlsc) else num(S$adls)
D$adl      <- adl

## physiology on day 3; Harrell's normal values for missing laboratory results
normal <- c(meanbp = NA, hrt = NA, resp = NA, temp = NA, sod = NA,
            pafi = 333.3, alb = 3.5, bili = 1.01, crea = 1.01, bun = 6.51,
            wblc = 9, urine = 2502, ph = 7.4, glucose = NA)
for (v in names(normal)) {
  x <- num(S[[v]])
  if (mean(is.na(x)) > 0.05) D[[paste0("miss_", v)]] <- as.integer(is.na(x))
  fill <- if (is.na(normal[[v]])) stats::median(x, na.rm = TRUE) else normal[[v]]
  x[is.na(x)] <- fill
  D[[v]] <- if (v %in% c("bili", "crea", "bun", "urine", "glucose", "wblc")) log1p(pmax(x, 0)) else x
}
## remaining missing values (few, non-laboratory): median, outcome-blind
for (v in names(D)) if (anyNA(D[[v]])) D[[v]][is.na(D[[v]])] <- stats::median(D[[v]], na.rm = TRUE)
D <- D[, vapply(D, function(v) stats::sd(v) > 0, TRUE), drop = FALSE]

## expanded candidate set -------------------------------------------------------------
miss_cols <- grep("^miss_", names(D), value = TRUE)
base <- setdiff(names(D), miss_cols)
Bc <- scale(as.matrix(D[, base]), center = TRUE, scale = TRUE)
phys <- intersect(c("age", "meanbp", "hrt", "resp", "temp", "sod", "ph", "pafi", "alb",
                    "bili", "crea", "bun", "wblc", "urine", "glucose", "scoma", "adl"), base)
SQ <- Bc[, phys, drop = FALSE]^2; colnames(SQ) <- paste0(phys, "^2")
pr <- utils::combn(base, 2L)
dzv <- grep("^dz_", base, value = TRUE); rcv <- grep("^race_", base, value = TRUE); cav <- grep("^ca_", base, value = TRUE)
same_factor <- function(a, b) (a %in% dzv && b %in% dzv) || (a %in% rcv && b %in% rcv) || (a %in% cav && b %in% cav)
is_bin <- vapply(base, function(v) all(D[[v]] %in% 0:1), TRUE)
rare <- function(a, b) is_bin[[a]] && is_bin[[b]] && sum(D[[a]] == 1 & D[[b]] == 1) < 50   # too few joint cases
keep <- !vapply(seq_len(ncol(pr)), function(j) same_factor(pr[1, j], pr[2, j]) || rare(pr[1, j], pr[2, j]), TRUE)
pr <- pr[, keep, drop = FALSE]
INT <- vapply(seq_len(ncol(pr)), function(j) Bc[, pr[1, j]] * Bc[, pr[2, j]], numeric(nrow(Bc)))
colnames(INT) <- paste(pr[1, ], pr[2, ], sep = ":")
X <- cbind(as.matrix(D[, base]), as.matrix(D[, miss_cols, drop = FALSE]), SQ, INT)
X <- X[, apply(X, 2L, stats::sd) > 1e-8, drop = FALSE]
time <- time / 365.25

tau_ok <- function(time, status, tau) {
  km <- survfit(Surv(time, status) ~ 1); Fh <- 1 - km$surv
  if (max(Fh) < tau + CFG$tau_margin) return(FALSE)
  sum(time >= km$time[which(Fh >= tau)[1L]]) >= CFG$min_at_risk
}
taus <- CFG$taus[vapply(CFG$taus, function(t) tau_ok(time, delta, t), TRUE)]
if (!length(taus)) stop("neither tau is identifiable")
save_atomic(list(name = NAME, cancer = "SUPPORT", x = X, y = log(time), time = time, delta = delta,
                 taus = taus, genes = colnames(X), clinical = base),
            file.path(RD_ROOT, "data", "clean", paste0(NAME, ".rds")))

d2 <- floor(sum(delta) / 2)
km <- survfit(Surv(time, delta) ~ 1)
q25 <- unname(quantile(km, 0.25)$quantile); q50 <- unname(quantile(km, 0.5)$quantile)
R <- stats::cor(X[, c(base, miss_cols)])
design <- data.frame(dataset = NAME, n = nrow(X), events = sum(delta), censoring = round(1 - mean(delta), 3),
                     km_q25_days = round(q25 * 365.25), km_median_days = round(q50 * 365.25),
                     base = length(base), missing_ind = length(miss_cols), squares = ncol(SQ),
                     interactions = sum(grepl(":", colnames(X))), p = ncol(X),
                     screen_cap = floor(d2 / log(d2)), taus = paste(taus, collapse = ";"),
                     max_abs_cor_base = round(max(abs(R[upper.tri(R)])), 2),
                     median_abs_cor_base = round(stats::median(abs(R[upper.tri(R)])), 3))
write.csv(design, file.path(TBS, "support_design.csv"), row.names = FALSE)
print(t(design))

tk <- do.call(rbind, lapply(taus, function(t) rbind(
  data.frame(dataset = NAME, tau = t, type = "main", rep = 1L),
  data.frame(dataset = NAME, tau = t, type = "mds", rep = 1L),
  data.frame(dataset = NAME, tau = t, type = "stab", rep = seq_len(CFG$n_stab)),
  data.frame(dataset = NAME, tau = t, type = "decoy", rep = seq_len(CFG$n_decoy)),
  data.frame(dataset = NAME, tau = t, type = "holdout", rep = seq_len(CFG$n_holdout)))))
tk <- tk[order(!tk$type %in% c("main", "mds"), tk$type, tk$tau), ]
tk$id <- seq_len(nrow(tk))
write.csv(tk, file.path(RD_ROOT, "tasks_support.csv"), row.names = FALSE)
cat(sprintf("\n%d tasks written to tasks_support.csv\n", nrow(tk)))
