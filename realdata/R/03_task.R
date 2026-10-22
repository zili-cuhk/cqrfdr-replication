## ---------------------------------------------------------------------------
## One analysis task (one row of tasks.csv)
## ---------------------------------------------------------------------------
## Usage: Rscript R/03_task.R <task id> [<number of consecutive tasks>]
## Types
##   main    : SS-MS (one split), L-HDCQR and Cox-lasso on the full data
##   mds     : MDS (B splits) on the full data, earlier implementation (penalty
##             tuned once, capped screen); kept so that old results stay reproducible
##   msms    : MS-MS (SS-MS on B splits, aggregated by inclusion rates) on the full
##             data; rep 1 is the analysis, rep 2 an independent run with other
##             split seeds (run-to-run agreement)
##   mds_dai : MDS as in Dai et al. (2023): every split a DS fit with its own
##             10-fold CV, no screening cap, their tie rule
##   cmp     : the comparators of the simulation study, each following its own
##             publication (cqr_comparators/): L-HDCQR of Zheng et al. (2018),
##             DS of Dai et al. (2023), SCQR-Lasso and SCQR-SCAD of He et al. (2022)
##   decoy_z18, holdout_z18 : L-HDCQR (Zheng et al.) in the decoy and hold-out
##             experiments, with the same decoys / splits as "decoy" / "holdout"
##   decoy_msms, holdout_msms : MS-MS in the same decoy and hold-out experiments
##   stab    : SS-MS with split seed = rep (split-to-split stability)
##   decoy   : n_decoy_vars permuted genes appended (null by construction);
##             SS-MS and L-HDCQR; counts how many decoys are selected
##   holdout : random 70/30 split; selection on the training part, prediction
##             (C-index, held-out Peng-Huang loss) on the test part
## Results: results/<dataset>/tau<tau>/<type>_<rep>.rds; finished tasks are skipped.

source(file.path(Sys.getenv("CQR_RD_ROOT", "."), "R", "config.R"))
suppressPackageStartupMessages({ library(cqrfdr); library(survival) })
if (utils::packageVersion("cqrfdr") < "0.2.0") stop("cqrfdr >= 0.2.0 is needed (msms, Dai-design mds); run R/00_install.R")
has_glmnet <- requireNamespace("glmnet", quietly = TRUE)

args <- commandArgs(TRUE)
if (!length(args)) stop("usage: Rscript R/03_task.R <task id> [count]")
first <- as.integer(args[1L]); cnt <- if (length(args) > 1L) as.integer(args[2L]) else 1L
tasks <- utils::read.csv(TASKS_FILE, stringsAsFactors = FALSE)

out_path <- function(tk) file.path(RD_ROOT, "results", tk$dataset, sprintf("tau%.2f", tk$tau),
                                   sprintf("%s_%03d.rds", tk$type, tk$rep))
timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(expr, error = function(e) list(error = conditionMessage(e)))
  r$seconds <- proc.time()[["elapsed"]] - t0
  r
}

## methods ----------------------------------------------------------------------
fit_ssms <- function(x, y, d, tau, seed) {
  set.seed(seed)
  f <- ssms(x, y, d, tau = tau, q = CFG$q, nu = CFG$nu, n_grid = CFG$n_grid,
            nfolds = CFG$nfolds, cores = CFG$cores)
  list(status = f$status, selected = lapply(f$selected, names),
       screen = names(f$screen), M = f$M[f$screen], b1 = f$b1[f$screen], b2 = f$b2[f$screen],
       threshold = f$threshold, mult = f$tuning$mult, cap = f$cap)
}
fit_mds <- function(x, y, d, tau, seed) {             # earlier implementation (cqrfdr 0.1.0 default)
  set.seed(seed)
  f <- mds(x, y, d, tau = tau, q = CFG$q, B = CFG$mds_B, nu = CFG$nu, n_grid = CFG$n_grid,
           nfolds = CFG$nfolds, tune = "cv_once", screen_cap = "d2", ties = "random",
           cores = CFG$cores)
  list(status = f$status, selected = lapply(f$selected, names),
       inclusion = f$inclusion[rowSums(f$inclusion) > 0, , drop = FALSE],
       n_failed_splits = f$n_failed_splits, mult = f$mult)
}
fit_multi <- function(fun, x, y, d, tau, seed, ...) {
  set.seed(seed)
  f <- fun(x, y, d, tau = tau, q = CFG$q, B = CFG$mds_B, nu = CFG$nu, n_grid = CFG$n_grid,
           cores = CFG$cores, ...)
  list(status = f$status, selected = lapply(f$selected, names),
       inclusion = f$inclusion[rowSums(f$inclusion) > 0, , drop = FALSE],
       n_failed_splits = f$n_failed_splits, share_empty = f$share_empty,
       split_sizes = f$split_sizes)
}
fit_msms    <- function(x, y, d, tau, seed) fit_multi(msms, x, y, d, tau, seed, nfolds = CFG$nfolds)
fit_mds_dai <- function(x, y, d, tau, seed) fit_multi(mds, x, y, d, tau, seed)   # Dai design: package default
## comparators of the simulation study, loaded into their own environment (the
## bundle defines ssms()/mds() that must not mask the cqrfdr package)
CMP <- local({
  e <- new.env()
  e$CQR_COMP_ROOT <- normalizePath(Sys.getenv("CQR_COMP_ROOT", file.path(RD_ROOT, "cqr_comparators")))
  for (f in c("dgp.R", "solver.R", "methods.R", "utils.R"))
    sys.source(file.path(e$CQR_COMP_ROOT, "R", f), envir = e)
  e
})
fit_cmp <- function(which, x, y, d, tau, seed) {
  set.seed(seed)
  taus <- taus_for(tau)
  r <- switch(which,
    lhdcqr = CMP$lhdcqr_z18(y, d, x, taus, K = CFG$nfolds, cores = CFG$cores),
    ds     = CMP$ds(y, d, x, taus, q = CFG$q, K = 10L, loss = "ph", cores = CFG$cores),
    scqr_lasso = CMP$scqr(y, d, x, taus, penalty = "lasso", K = CFG$nfolds),
    scqr_scad  = CMP$scqr(y, d, x, taus, penalty = "scad", K = CFG$nfolds))
  nm <- function(idx) colnames(x)[idx]
  if (!is.null(r$by_q))
    list(status = r$status, selected = stats::setNames(lapply(r$by_q, function(l) nm(l$selected)),
                                                       vapply(r$by_q, function(l) format(l$q), "")))
  else list(status = r$status, selected = nm(r$selected), failed = isTRUE(r$failed))
}
fit_hdcqr <- function(x, y, d, tau, seed) {
  set.seed(seed)
  f <- hdcqr(x, y, d, tau = tau, nu = CFG$nu, n_grid = CFG$n_grid, nfolds = CFG$nfolds,
             cores = CFG$cores)
  list(status = "ok", selected = names(f$selected), mult = f$mult, fit = f)
}
fit_cox <- function(x, time, d, seed) {
  if (!has_glmnet) return(list(status = "glmnet_not_installed", selected = character(0)))
  set.seed(seed)
  extra <- if ("cox.ties" %in% names(formals(glmnet::glmnet))) list(cox.ties = "breslow") else list()
  cv <- do.call(glmnet::cv.glmnet, c(list(x = x, y = survival::Surv(time, d), family = "cox",
                                          nfolds = CFG$nfolds), extra))
  b <- as.numeric(stats::coef(cv, s = "lambda.min"))
  list(status = "ok", selected = colnames(x)[b != 0], beta = stats::setNames(b, colnames(x)),
       lambda = cv$lambda.min)
}

## prediction on a test set ----------------------------------------------------------
## Unpenalized Peng-Huang refit on the training data with the selected genes;
## returns the predicted tau-quantile of log time on the test data and the
## held-out Peng-Huang loss (smaller is better).
taus_for <- function(tau) seq(CFG$nu, tau, length.out = CFG$n_grid)
refit_predict <- function(xtr, ytr, dtr, xte, yte, dte, genes, tau) {
  taus <- taus_for(tau); m <- length(taus)
  ctr <- colMeans(xtr[, genes, drop = FALSE]); scl <- apply(xtr[, genes, drop = FALSE], 2L, stats::sd)
  Ztr <- cbind(1, scale(xtr[, genes, drop = FALSE], ctr, scl))
  Zte <- cbind(1, scale(xte[, genes, drop = FALSE], ctr, scl))
  f <- hdcqr_path(ytr, dtr, Ztr, taus, 0)
  if (f$status != "ok") return(list(status = f$status, pred = rep(NA_real_, nrow(xte)), loss = NA_real_))
  list(status = "ok", pred = drop(Zte %*% f$B[, m]),
       loss = cqrfdr:::cv_loss(yte, dte, Zte, f$B, taus, "ph"))
}
cidx <- function(time, d, pred, reverse = FALSE) {
  if (all(!is.finite(pred)) || stats::sd(pred, na.rm = TRUE) == 0) return(0.5)
  unname(survival::concordance(survival::Surv(time, d) ~ pred, reverse = reverse)$concordance)
}

## run ---------------------------------------------------------------------------------
cache <- new.env()
get_data <- function(nm) {
  if (is.null(cache[[nm]])) cache[[nm]] <- readRDS(file.path(RD_ROOT, "data", "clean", paste0(nm, ".rds")))
  cache[[nm]]
}
for (id in first:(first + cnt - 1L)) {
  if (id > nrow(tasks)) break
  tk <- tasks[tasks$id == id, ]
  path <- out_path(tk)
  if (file.exists(path)) { cat("task", id, "done already\n"); next }
  D <- get_data(tk$dataset); x <- D$x; y <- D$y; d <- D$delta; tm <- D$time; tau <- tk$tau
  sd0 <- function(method) str_seed(tk$dataset, tau, tk$type, tk$rep, method)
  cat(sprintf("task %d: %s tau=%.2f %s rep %d (n=%d, p=%d)\n", id, tk$dataset, tau, tk$type,
              tk$rep, nrow(x), ncol(x)))
  res <- list(task = tk)
  if (tk$type == "main") {
    res$ssms   <- timed(fit_ssms(x, y, d, tau, sd0("ssms")))
    res$hdcqr  <- timed({ h <- fit_hdcqr(x, y, d, tau, sd0("hdcqr")); h$fit <- NULL; h })
    res$cox    <- timed(fit_cox(x, tm, d, sd0("cox")))
  } else if (tk$type == "mds") {
    res$mds    <- timed(fit_mds(x, y, d, tau, sd0("mds")))
  } else if (tk$type == "msms") {
    res$msms <- timed(fit_msms(x, y, d, tau, sd0("msms")))
  } else if (tk$type == "mds_dai") {
    res$mds_dai <- timed(fit_mds_dai(x, y, d, tau, sd0("mds_dai")))
  } else if (tk$type == "cmp") {
    for (m in c("lhdcqr", "ds", "scqr_lasso", "scqr_scad"))
      res[[paste0("cmp_", m)]] <- timed(fit_cmp(m, x, y, d, tau, sd0(m)))
  } else if (tk$type == "decoy_z18") {
    set.seed(str_seed(tk$dataset, tau, "decoy", tk$rep, "decoy"))   # the decoys of task "decoy"
    src <- sample.int(ncol(x), min(CFG$n_decoy_vars, ncol(x)))
    xd <- x[sample.int(nrow(x)), src, drop = FALSE]
    colnames(xd) <- paste0("DECOY_", colnames(x)[src])
    res$lhdcqr_z18 <- timed(fit_cmp("lhdcqr", cbind(x, xd), y, d, tau, sd0("lhdcqr")))
    res$n_decoys <- ncol(xd)
  } else if (tk$type == "holdout_z18") {
    set.seed(str_seed(tk$dataset, tau, "holdout", tk$rep, "split"))  # the split of task "holdout"
    tr <- integer(0)
    for (lev in 0:1) { ii <- which(d == lev); tr <- c(tr, ii[sample.int(length(ii), round(CFG$holdout_frac * length(ii)))]) }
    tr <- sort(tr); te <- setdiff(seq_len(nrow(x)), tr)
    h <- timed(fit_cmp("lhdcqr", x[tr, ], y[tr], d[tr], tau, sd0("lhdcqr")))
    res$lhdcqr_z18 <- if (is.null(h$error)) {
      r <- refit_predict(x[tr, , drop = FALSE], y[tr], d[tr], x[te, , drop = FALSE], y[te], d[te], h$selected, tau)
      c(n_sel = length(h$selected), C = cidx(tm[te], d[te], r$pred), loss = r$loss)
    } else h
  } else if (tk$type == "decoy_msms") {
    set.seed(str_seed(tk$dataset, tau, "decoy", tk$rep, "decoy"))   # the decoys of task "decoy"
    src <- sample.int(ncol(x), min(CFG$n_decoy_vars, ncol(x)))
    xd <- x[sample.int(nrow(x)), src, drop = FALSE]
    colnames(xd) <- paste0("DECOY_", colnames(x)[src])
    res$msms <- timed(fit_msms(cbind(x, xd), y, d, tau, sd0("msms")))
    res$n_decoys <- ncol(xd)
  } else if (tk$type == "holdout_msms") {
    set.seed(str_seed(tk$dataset, tau, "holdout", tk$rep, "split"))  # the split of task "holdout"
    tr <- integer(0)
    for (lev in 0:1) { ii <- which(d == lev); tr <- c(tr, ii[sample.int(length(ii), round(CFG$holdout_frac * length(ii)))]) }
    tr <- sort(tr); te <- setdiff(seq_len(nrow(x)), tr)
    ev <- function(genes) {
      r <- refit_predict(x[tr, , drop = FALSE], y[tr], d[tr], x[te, , drop = FALSE], y[te], d[te], genes, tau)
      c(n_sel = length(genes), C = cidx(tm[te], d[te], r$pred), loss = r$loss)
    }
    m <- timed(fit_msms(x[tr, ], y[tr], d[tr], tau, sd0("msms")))
    res$msms <- if (is.null(m$error)) lapply(m$selected, ev) else m
  } else if (tk$type == "stab") {
    res$ssms <- timed(fit_ssms(x, y, d, tau, sd0("ssms")))
  } else if (tk$type == "decoy") {
    set.seed(sd0("decoy"))
    src <- sample.int(ncol(x), min(CFG$n_decoy_vars, ncol(x)))
    xd <- x[sample.int(nrow(x)), src, drop = FALSE]          # one row permutation: decoys keep
    colnames(xd) <- paste0("DECOY_", colnames(x)[src])       # their mutual correlation but are
    xa <- cbind(x, xd)                                        # independent of the outcome
    res$ssms  <- timed(fit_ssms(xa, y, d, tau, sd0("ssms")))
    res$hdcqr <- timed({ h <- fit_hdcqr(xa, y, d, tau, sd0("hdcqr")); h$fit <- NULL; h })
    res$cox   <- timed(fit_cox(xa, tm, d, sd0("cox")))
    res$n_decoys <- ncol(xd)
  } else if (tk$type == "holdout") {
    set.seed(sd0("split"))
    tr <- integer(0)
    for (lev in 0:1) { ii <- which(d == lev); tr <- c(tr, ii[sample.int(length(ii), round(CFG$holdout_frac * length(ii)))]) }
    tr <- sort(tr); te <- setdiff(seq_len(nrow(x)), tr)
    ev <- function(genes) {
      r <- refit_predict(x[tr, , drop = FALSE], y[tr], d[tr], x[te, , drop = FALSE], y[te], d[te], genes, tau)
      c(n_sel = length(genes), C = cidx(tm[te], d[te], r$pred), loss = r$loss)
    }
    base <- ev(character(0))
    s <- timed(fit_ssms(x[tr, ], y[tr], d[tr], tau, sd0("ssms")))
    res$ssms <- if (is.null(s$error)) lapply(s$selected, ev) else s
    h <- timed(fit_hdcqr(x[tr, ], y[tr], d[tr], tau, sd0("hdcqr")))
    if (is.null(h$error)) {
      taus <- taus_for(tau)
      Zte <- cbind(1, scale(x[te, , drop = FALSE], h$fit$center, h$fit$scale))
      er <- ev(h$selected)
      res$hdcqr <- c(n_sel = length(h$selected),
                     C = cidx(tm[te], d[te], predict(h$fit, x[te, , drop = FALSE])),
                     loss = cqrfdr:::cv_loss(y[te], d[te], Zte, h$fit$coef_std, taus, "ph"),
                     C_refit = unname(er["C"]), loss_refit = unname(er["loss"]))
    } else res$hdcqr <- h
    cx <- timed(fit_cox(x[tr, ], tm[tr], d[tr], sd0("cox")))
    res$cox <- if (identical(cx$status, "ok")) {
      er <- ev(cx$selected)
      c(n_sel = length(cx$selected), C = cidx(tm[te], d[te], drop(x[te, ] %*% cx$beta), reverse = TRUE),
        loss = NA_real_, C_refit = unname(er["C"]), loss_refit = unname(er["loss"]))
    } else cx
    res$base <- base
  }
  save_atomic(res, path)
  cat("  saved", path, "\n")
}
