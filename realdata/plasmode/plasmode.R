## ---------------------------------------------------------------------------
## Plasmode (semi-synthetic) simulation on the real METABRIC design
## ---------------------------------------------------------------------------
## The covariate matrix is the real METABRIC expression matrix (MAD-filtered,
## 2000 genes, n = 1979), so the correlation structure is the real one.  For
## each replicate, s genes are drawn at random as true signals with effects
## +-b (genes standardized), and
##     log T = X_S beta + eps,  eps ~ N(0, 1),
##     log C = c0 + N(0, 1),    c0 calibrated to the METABRIC censoring rate.
## In this location model the active set is the same at every quantile, so
## FDR and power at tau are known exactly.
##
## Usage: Rscript plasmode/plasmode.R <task id>
## Results: results/plasmode/task_<id>.rds
## Methods (CQR_PM_METHODS, comma-separated; default all):
##   ssms     SS-MS, path and CV screens, cutoffs U+ and U
##   msms     MS-MS (B = mds_B splits)
##   mds_dai  MDS as in Dai et al. (2023) (cqrfdr >= 0.2.0 default)
##   mds      MDS, earlier implementation (penalty tuned once, capped screen)
##   hdcqr    L-HDCQR with one cross-validated penalty level (cqrfdr::hdcqr)
##   lhdcqr_z18  L-HDCQR as in Zheng, Peng and He (2018) (cqr_comparators/)
## If the result file exists, only the requested methods that are missing from
## it are run and appended (same data: the seeds depend on s, b and rep only).

source(file.path(Sys.getenv("CQR_RD_ROOT", "."), "R", "config.R"))
suppressPackageStartupMessages(library(cqrfdr))
PM_TEST <- nzchar(Sys.getenv("CQR_PM_TEST"))

PM <- list(
  data   = Sys.getenv("CQR_PM_DATA", "BRCA_METABRIC"),
  cens   = NA_real_,                   # censoring rate of the real data (0.42 for METABRIC)
  tau    = 0.5,
  grid   = expand.grid(s = as.integer(strsplit(Sys.getenv("CQR_PM_S", "10,30"), ",")[[1]]),
                       b = as.numeric(strsplit(Sys.getenv("CQR_PM_B", "0.1,0.2"), ",")[[1]])),
  reps   = 50L,
  q      = c(0.1, 0.2),
  n_grid = 40L, nfolds = 5L, mds_B = 50L,
  cores  = CFG$cores,
  methods = strsplit(Sys.getenv("CQR_PM_METHODS", "ssms,msms,mds_dai,mds,hdcqr"), ",")[[1]]
)
if (utils::packageVersion("cqrfdr") < "0.2.0") stop("cqrfdr >= 0.2.0 is needed; run R/00_install.R")
if (PM_TEST) PM[c("reps", "n_grid", "nfolds", "mds_B", "cores")] <- list(2L, 10L, 3L, 3L, 1L)
n_tasks <- function() nrow(PM$grid) * PM$reps
PM_DIR <- file.path(RD_ROOT, "results",
                    if (PM$data == "BRCA_METABRIC") "plasmode" else paste0("plasmode_", PM$data))

run_task <- function(id) {
  out <- file.path(PM_DIR, sprintf("task_%04d.rds", id))
  old <- if (file.exists(out)) readRDS(out) else NULL
  label_of <- c(ssms = "SS-MS (path, U+)", msms = "MS-MS", mds_dai = "MDS-Dai", mds = "MDS", hdcqr = "L-HDCQR",
                lhdcqr_z18 = "L-HDCQR-Z18")
  todo <- PM$methods
  if (!is.null(old)) todo <- todo[!label_of[todo] %in% unique(old$rows$method)]
  if (!length(todo)) { cat("task", id, "done already\n"); return(invisible()) }
  g <- (id - 1L) %/% PM$reps + 1L; rep <- (id - 1L) %% PM$reps + 1L
  s <- PM$grid$s[g]; b <- PM$grid$b[g]
  D <- readRDS(file.path(RD_ROOT, "data", "clean", paste0(PM$data, ".rds")))
  X <- scale(D$x); n <- nrow(X); p <- ncol(X)
  cens <- if (PM$data == "BRCA_METABRIC") 0.42 else 1 - mean(D$delta)
  set.seed(str_seed("plasmode", s, b, rep))
  active <- sort(sample.int(p, s))
  beta <- b * sample(c(-1, 1), s, replace = TRUE)
  logT <- drop(X[, active, drop = FALSE] %*% beta) + stats::rnorm(n)
  gC <- stats::rnorm(n)
  c0 <- stats::uniroot(function(c) mean(c + gC < logT) - cens, range(logT) + c(-10, 10))$root
  logC <- c0 + gC
  y <- pmin(logT, logC); d <- as.integer(logT <= logC)
  cat(sprintf("task %d: s=%d b=%.2f rep %d  n=%d p=%d events=%d\n", id, s, b, rep, n, p, sum(d)))

  res <- if (is.null(old)) list(id = id, s = s, b = b, rep = rep, active = active, events = sum(d), rows = list())
         else { if (!identical(old$active, active)) stop("task ", id, ": regenerated data differ from the stored result")
                o <- old; o$rows <- list(old$rows); o }
  cat("  running:", paste(todo, collapse = ", "), "\n")
  add <- function(method, q, sel, sec) res$rows[[length(res$rows) + 1L]] <<- data.frame(
    s = s, b = b, rep = rep, method = method, q = q, n_sel = length(sel),
    FDP = fdp_power(sel, active)[["FDP"]], power = fdp_power(sel, active)[["power"]],
    seconds = sec)
  tm <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- tryCatch(expr, error = function(e) e)
                         list(v = v, sec = proc.time()[["elapsed"]] - t0) }
  sd_split <- str_seed("plasmode-split", s, b, rep)
  if ("ssms" %in% todo) for (sc in c("path", "cv")) {
    r <- tm({ set.seed(sd_split)              # same split and folds for both screens
              ssms(X, y, d, tau = PM$tau, q = PM$q, n_grid = PM$n_grid, nfolds = PM$nfolds,
                   screen = sc, cores = PM$cores) })
    if (inherits(r$v, "error")) { cat("  ssms", sc, "error:", conditionMessage(r$v), "\n"); next }
    f <- r$v; Ms <- f$M[f$screen]
    for (q in PM$q) {
      add(sprintf("SS-MS (%s, U+)", sc), q, selected(f, q, names = FALSE), r$sec)
      thr0 <- if (length(Ms)) mirror_threshold(Ms, q, 0) else Inf
      add(sprintf("SS-MS (%s, U)", sc), q, which(f$M >= thr0), r$sec)
    }
    res[[paste0("screen_", sc)]] <- list(size = length(f$screen), recall = mean(active %in% f$screen))
  }
  if ("msms" %in% todo) {
    r <- tm({ set.seed(str_seed("plasmode-msms", s, b, rep))
              msms(X, y, d, tau = PM$tau, q = PM$q, B = PM$mds_B, n_grid = PM$n_grid,
                   nfolds = PM$nfolds, cores = PM$cores) })
    if (!inherits(r$v, "error")) {
      for (q in PM$q) add("MS-MS", q, selected(r$v, q, names = FALSE), r$sec)
      res$msms_share_empty <- r$v$share_empty
    } else cat("  msms error:", conditionMessage(r$v), "\n")
  }
  if ("mds_dai" %in% todo) {
    r <- tm({ set.seed(str_seed("plasmode-mds-dai", s, b, rep))
              mds(X, y, d, tau = PM$tau, q = PM$q, B = PM$mds_B, n_grid = PM$n_grid, cores = PM$cores) })
    if (!inherits(r$v, "error")) for (q in PM$q) add("MDS-Dai", q, selected(r$v, q, names = FALSE), r$sec)
    else cat("  mds_dai error:", conditionMessage(r$v), "\n")
  }
  if ("mds" %in% todo) {
    r <- tm({ set.seed(str_seed("plasmode-mds", s, b, rep))
              mds(X, y, d, tau = PM$tau, q = PM$q, B = PM$mds_B, n_grid = PM$n_grid,
                  nfolds = PM$nfolds, tune = "cv_once", screen_cap = "d2", ties = "random",
                  cores = PM$cores) })
    if (!inherits(r$v, "error")) for (q in PM$q) add("MDS", q, selected(r$v, q, names = FALSE), r$sec)
  }
  if ("lhdcqr_z18" %in% todo) {
    CMP <- new.env()
    CMP$CQR_COMP_ROOT <- normalizePath(Sys.getenv("CQR_COMP_ROOT", file.path(RD_ROOT, "cqr_comparators")))
    for (f in c("dgp.R", "solver.R", "methods.R", "utils.R"))
      sys.source(file.path(CMP$CQR_COMP_ROOT, "R", f), envir = CMP)
    r <- tm({ set.seed(str_seed("plasmode-lhdcqr-z18", s, b, rep))
              CMP$lhdcqr_z18(y, d, X, seq(0.1, PM$tau, length.out = PM$n_grid), K = PM$nfolds, cores = PM$cores) })
    if (!inherits(r$v, "error") && !isTRUE(r$v$failed)) add("L-HDCQR-Z18", NA, r$v$selected, r$sec)
    else cat("  lhdcqr_z18 failed:", if (inherits(r$v, "error")) conditionMessage(r$v) else r$v$status, "\n")
  }
  if ("hdcqr" %in% todo) {
  r <- tm({ set.seed(str_seed("plasmode-hdcqr", s, b, rep))
            hdcqr(X, y, d, tau = PM$tau, n_grid = PM$n_grid, nfolds = PM$nfolds, cores = PM$cores) })
  if (!inherits(r$v, "error")) add("L-HDCQR", NA, selected(r$v, names = FALSE), r$sec)
  }
  res$rows <- do.call(rbind, res$rows)
  save_atomic(res, out)
  print(res$rows[, c("method", "q", "n_sel", "FDP", "power")], row.names = FALSE)
}

if (sys.nframe() == 0L) {
  a <- commandArgs(TRUE)
  if (!length(a)) { cat("number of tasks:", n_tasks(), "\n"); quit("no") }
  run_task(as.integer(a[1L]))
}
