## ---------------------------------------------------------------------------
## Method-agnostic signal strength of every analysed data set
## ---------------------------------------------------------------------------
## Usage: Rscript R/07_signal.R            (about 1-3 minutes per data set)
## Writes tables/signal_strength.csv.
##
## Two measures that do not involve SS-MS, MDS or L-HDCQR:
##   * n_bh05 : number of candidate variables with univariate Cox BH-adjusted
##              p < 0.05 (how many variables carry a marginal signal);
##   * cv_C   : 5-fold cross-validated C-index of a ridge Cox model on all
##              candidates (how much the candidates predict jointly).
## The table documents, for every data set, how strong the signal is, so the
## choice of the strong-signal example can be checked against a rule that
## does not depend on the proposed method.

source(file.path(Sys.getenv("CQR_RD_ROOT", "."), "R", "config.R"))
suppressPackageStartupMessages({ library(survival); library(glmnet) })
fs <- list.files(file.path(RD_ROOT, "data", "clean"), "\\.rds$", full.names = TRUE)
a <- commandArgs(TRUE); if (length(a)) fs <- fs[sub("\\.rds$", "", basename(fs)) %in% a]

one <- function(f) {
  D <- readRDS(f); X <- scale(D$x); S <- Surv(D$time, D$delta)
  p_uni <- apply(X, 2L, function(v) {
    fit <- tryCatch(coxph(S ~ v), error = function(e) NULL)
    if (is.null(fit)) NA_real_ else summary(fit)$coefficients[1L, 5L] })
  set.seed(str_seed("signal", D$name))
  fold <- sample(rep_len(seq_len(5L), nrow(X)))
  lp <- numeric(nrow(X))
  for (k in 1:5) {
    tr <- fold != k
    cv <- cv.glmnet(X[tr, ], S[tr], family = "cox", alpha = 0, nfolds = 5)
    lp[!tr] <- drop(predict(cv, X[!tr, , drop = FALSE], s = "lambda.min"))
  }
  cC <- concordance(S ~ lp, reverse = TRUE)$concordance
  data.frame(dataset = D$name, n = nrow(X), events = sum(D$delta), p = ncol(X),
             n_bh05 = sum(stats::p.adjust(p_uni, "BH") < 0.05, na.rm = TRUE),
             share_bh05 = round(mean(stats::p.adjust(p_uni, "BH") < 0.05, na.rm = TRUE), 3),
             cv_C = round(cC, 3))
}
out <- do.call(rbind, lapply(fs, function(f) {
  cat(basename(f), "... "); r <- tryCatch(one(f), error = function(e) { cat("error:", conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) cat(sprintf("cv_C = %.3f, n_bh05 = %d\n", r$cv_C, r$n_bh05)); r }))
out <- out[order(-out$cv_C), ]
dir.create(file.path(RD_ROOT, "tables"), showWarnings = FALSE)
write.csv(out, file.path(RD_ROOT, "tables", "signal_strength.csv"), row.names = FALSE)
print(out, row.names = FALSE)
