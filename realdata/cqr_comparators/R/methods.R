## ---------------------------------------------------------------------------
## Selection methods
## ---------------------------------------------------------------------------
## Every method returns a list with
##   status    "ok" | "empty_screen" | "no_threshold" | <failure code>
##   selected  integer indices in 1..p (empty on failure)
##   failed    TRUE if the method could not produce a selection
## plus method-specific diagnostics.

## Mirror statistic and data-driven threshold --------------------------------
mirror_stat <- function(b1, b2) sign(b1 * b2) * (abs(b1) + abs(b2))

## smallest t > 0 with (offset + #{M <= -t}) / max(#{M >= t}, 1) <= q
## offset = 1: U_+ (finite-sample bound, Theorem 1); offset = 0: U
mirror_threshold <- function(M, q, offset = 1) {
  ts <- sort(unique(abs(M[M != 0])))
  if (!length(ts)) return(Inf)
  neg <- vapply(ts, function(t) sum(M <= -t), numeric(1))
  pos <- vapply(ts, function(t) sum(M >= t), numeric(1))
  ok <- which((offset + neg) / pmax(pos, 1) <= q)
  if (length(ok)) ts[ok[1L]] else Inf
}

## std = "sd": mean / standard deviation (the draft); "robust": median / MAD,
## which stays stable when covariates have no second moment (t_1, t_2).
standardize_by <- function(X, rows, std = c("sd", "robust")) {
  std <- match.arg(std)
  Xr <- X[rows, , drop = FALSE]
  if (std == "sd") {
    ctr <- colMeans(Xr)
    scl <- sqrt(colSums(sweep(Xr, 2L, ctr)^2) / (length(rows) - 1))
  } else {
    ctr <- apply(Xr, 2L, stats::median)
    scl <- apply(Xr, 2L, stats::mad)
  }
  if (any(!is.finite(scl) | scl <= 1e-10)) stop("degenerate covariate on the fitting subset")
  cbind(`(Intercept)` = 1, sweep(sweep(X, 2L, ctr), 2L, scl, "/"))
}

fail <- function(status, ...) c(list(status = status, selected = integer(0), failed = TRUE), list(...))

## Single split with mirror statistics (Algorithm 1) -------------------------
## screen_cap: "theory" keeps at most floor(d2 / log d2) candidates with the
##   largest |b1|, d2 = number of events in split 2 (Condition 5(b) with the
##   effective sample size); an integer sets the cap; NULL keeps the full
##   support.  Without a cap the refit fails whenever the candidate set is
##   large relative to d2 (this is what produced the draft's 56 SS-MS
##   "numerical failures" at 50% censoring).
ssms <- function(y, delta, X, taus, q = 0.2, offset = 1,
                 lambda_rule = "cv", mult = NULL, fixed_mult = 0.5, K = 5L,
                 loss = "ph", screen_cap = "theory", coef_tol = 1e-6,
                 cores = 1L, split = NULL, std = "sd", foldid = NULL,
                 screen_rule = c("cv", "size", "path")) {
  screen_rule <- match.arg(screen_rule)
  n <- length(y); p <- ncol(X); m <- length(taus)
  i1 <- if (is.null(split)) sort(sample.int(n, floor(n / 2))) else sort(split)
  i2 <- setdiff(seq_len(n), i1)
  Z <- standardize_by(X, i1, std)                # split-1 centring/scaling
  targs <- list(y = y[i1], delta = delta[i1], Z = Z[i1, , drop = FALSE], taus = taus,
                rule = lambda_rule, fixed_mult = fixed_mult, K = K, loss = loss,
                cores = cores)
  ## shared partition restricted to split 1 (if a partition was supplied)
  if (!is.null(foldid)) targs$foldid <- foldid[i1]
  if (!is.null(mult)) targs$mult <- mult
  tun <- tryCatch(do.call(tune_lambda, targs), error = function(e) e)
  if (inherits(tun, "error")) return(fail("tuning_failed", message = conditionMessage(tun)))
  f1 <- hdcqr_path(y[i1], delta[i1], Z[i1, , drop = FALSE], taus, tun$lambda)
  if (f1$status != "ok") return(fail(paste0("stage1_", f1$status), tau_fail = taus[f1$k]))
  b1 <- f1$B[-1L, m]
  screen <- which(abs(b1) > coef_tol)
  d2 <- sum(delta[i2])
  cap <- if (identical(screen_cap, "theory")) floor(d2 / log(d2)) else
    if (identical(screen_cap, "none")) NULL else screen_cap
  ## screen_rule = "size": if the CV fit keeps fewer than `cap` candidates,
  ## move down the multiplier grid (smaller penalties) until at least `cap`
  ## candidates are nonzero, then keep the `cap` largest |b1|.  This uses
  ## the full candidate budget allowed by Condition 5(b) to protect sure
  ## screening (Condition 1); nothing depends on split 2.
  mult_screen <- tun$mult
  if (screen_rule == "size" && !is.null(cap) && length(screen) < cap && !is.null(tun$grid)) {
    for (mm in rev(tun$grid[tun$grid < tun$mult])) {
      f <- hdcqr_path(y[i1], delta[i1], Z[i1, , drop = FALSE], taus, lambda_grid(mm * tun$pivotal, taus))
      if (f$status != "ok") break
      f1 <- f; b1 <- f1$B[-1L, m]; mult_screen <- mm
      screen <- which(abs(b1) > coef_tol)
      if (length(screen) >= cap) break
    }
  }
  ## screen_rule = "path": as "size", but candidates are ranked by the order in
  ## which they enter the lasso path as the multiplier decreases from the CV
  ## choice (ties: larger |coefficient| first), and each candidate keeps its
  ## coefficient from the fit at which it entered.  The CV support is always
  ## kept (when it fits under the cap), so the screen can only gain recall
  ## relative to "cv"; the late, nearly unpenalized fits at upper quantiles,
  ## whose largest coefficients are mostly noise, no longer decide the screen.
  ## b1 stays D1-measurable, so the mirror-statistic argument is unchanged.
  if (screen_rule == "path" && !is.null(cap) && length(screen) < cap && !is.null(tun$grid)) {
    ord <- screen[order(-abs(b1[screen]))]
    b1e <- numeric(p); b1e[screen] <- b1[screen]
    for (mm in rev(tun$grid[tun$grid < tun$mult])) {
      f <- hdcqr_path(y[i1], delta[i1], Z[i1, , drop = FALSE], taus, lambda_grid(mm * tun$pivotal, taus))
      if (f$status != "ok") break
      bb <- f$B[-1L, m]; mult_screen <- mm
      new <- setdiff(which(abs(bb) > coef_tol), ord)
      if (length(new)) {
        new <- new[order(-abs(bb[new]))]
        b1e[new] <- bb[new]; ord <- c(ord, new)
      }
      if (length(ord) >= cap) break
    }
    screen <- sort(ord[seq_len(min(cap, length(ord)))]); b1 <- b1e
  }
  full_screen <- length(screen)
  if (!is.null(cap) && length(screen) > cap)
    screen <- sort(screen[order(-abs(b1[screen]))[seq_len(cap)]])
  out <- list(status = "ok", selected = integer(0), failed = FALSE,
              screen = screen, screen_full_size = full_screen,
              mult = tun$mult, mult_screen = mult_screen, threshold = Inf, M = NULL)
  if (!length(screen)) { out$status <- "empty_screen"; return(out) }
  f2 <- hdcqr_path(y[i2], delta[i2], Z[i2, c(1L, screen + 1L), drop = FALSE], taus, 0)
  if (f2$status != "ok") return(fail(paste0("stage2_", f2$status), tau_fail = taus[f2$k],
                                     screen = screen))
  b2 <- f2$B[-1L, m]
  M <- numeric(p); M[screen] <- mirror_stat(b1[screen], b2)
  out$M <- M[screen]
  ## one selection per nominal level: only the threshold depends on q
  out$by_q <- lapply(q, function(qq) {
    thr <- mirror_threshold(M, qq, offset)
    list(q = qq, selected = which(M >= thr), threshold = thr,
         status = if (is.finite(thr)) "ok" else "no_threshold")
  })
  out$selected <- out$by_q[[1L]]$selected
  out$threshold <- out$by_q[[1L]]$threshold
  out$status <- out$by_q[[1L]]$status
  out
}

## Multiple data splitting (Dai et al., 2023) ---------------------------------
## Each split runs Algorithm 1 (threshold offset 0, as in Dai et al.) for every
## q; the inclusion rates are aggregated separately for each q.  Splits that
## fail are dropped and counted (n_failed_splits).
## tune = "cv_once": the penalty multiplier is chosen once per data set by the
##   same K-fold CV as SS-MS (on one random half, with the shared folds) and
##   reused, with each split's own pivotal level, in all B splits; "pivotal"
##   uses fixed_mult (the earlier, cheaper choice).
## ties = "random": the aggregation removes exactly the l smallest inclusion
##   rates, ties broken at random.  With ties = "dai" the rule
##   {j: I_j > I_(l)} removes a whole block of tied rates; when every split
##   selects all s signals their rates are tied, and if the null mass is below
##   q the block of signals is removed together (selection becomes empty).
mds <- function(y, delta, X, taus, q = 0.2, B = 50L, offset = 0,
                tune = c("cv_once", "pivotal", "cv_each"), ties = c("random", "dai"),
                lambda_rule = "pivotal", fixed_mult = 0.5, K = 5L, loss = "ph",
                foldid = NULL, cores = 1L, ...) {
  tune <- match.arg(tune); ties <- match.arg(ties)
  p <- ncol(X); n <- length(y)
  mult_used <- fixed_mult
  if (tune == "cv_once") {
    i1 <- sort(sample.int(n, floor(n / 2)))
    Z <- standardize_by(X, i1, "sd")
    targs <- list(y = y[i1], delta = delta[i1], Z = Z[i1, , drop = FALSE], taus = taus,
                  rule = "cv", K = K, loss = loss, cores = cores)
    if (!is.null(foldid)) targs$foldid <- foldid[i1]
    tun <- tryCatch(do.call(tune_lambda, targs), error = function(e) e)
    if (inherits(tun, "error")) return(fail("tuning_failed", message = conditionMessage(tun)))
    mult_used <- tun$mult
  }
  incl <- matrix(0, p, length(q)); nf <- 0L
  for (b in seq_len(B)) {
    ## tune = "cv_each": every split is a complete single-split fit with its own
    ## K-fold CV on its first half (Dai et al.: MDS repeats DS)
    r <- if (tune == "cv_each")
           ssms(y, delta, X, taus, q = q, offset = offset, lambda_rule = "cv", K = K,
                loss = loss, cores = cores, foldid = NULL, ...)
         else ssms(y, delta, X, taus, q = q, offset = offset, lambda_rule = "pivotal",
                   fixed_mult = mult_used, cores = cores, ...)
    if (isTRUE(r$failed)) { nf <- nf + 1L; next }
    if (is.null(r$by_q)) next                       # empty screen: selects nothing
    for (j in seq_along(q)) {
      sel <- r$by_q[[j]]$selected
      if (length(sel)) incl[sel, j] <- incl[sel, j] + 1 / length(sel)
    }
  }
  if (nf == B) return(fail("all_splits_failed", n_failed_splits = nf))
  incl <- incl / (B - nf)
  by_q <- lapply(seq_along(q), function(j) {
    I <- incl[, j]
    o <- if (ties == "dai") order(I) else order(I, runif(p))
    cs <- cumsum(I[o]); l <- max(c(0L, which(cs <= q[j])))
    sel <- if (ties == "dai") {
      cut <- if (l == 0L) -Inf else I[o][l]
      which(I > cut & I > 0)
    } else {
      keep <- if (l < p) o[(l + 1L):p] else integer(0)
      sort(keep[I[keep] > 0])
    }
    list(q = q[j], selected = sel, status = "ok")
  })
  list(status = "ok", selected = by_q[[1L]]$selected, failed = FALSE,
       n_failed_splits = nf, by_q = by_q, mult = mult_used)
}

## L-HDCQR of Zheng et al. (2018) on the full sample ---------------------------
## flat = TRUE, loss = "deviance_z18": the procedure of Zheng et al. (2018):
## one penalty level for all grid points, chosen by K-fold CV of the
## integrated absolute deviance residual (Fei et al., 2021, Remark 3).  The
## solver uses the augmented LAD problem at level 0.5 at every grid point.
lhdcqr <- function(y, delta, X, taus, lambda_rule = "cv", mult = NULL, K = 5L,
                   loss = "ph", coef_tol = 1e-6, cores = 1L, std = "sd",
                   foldid = NULL, flat = FALSE, ...) {
  m <- length(taus)
  Z <- standardize_by(X, seq_along(y), std)
  targs <- list(y = y, delta = delta, Z = Z, taus = taus, rule = lambda_rule,
                K = K, loss = loss, cores = cores, foldid = foldid, flat = flat)
  if (!is.null(mult)) targs$mult <- mult
  tun <- tryCatch(do.call(tune_lambda, targs), error = function(e) e)
  if (inherits(tun, "error")) return(fail("tuning_failed", message = conditionMessage(tun)))
  f <- hdcqr_path(y, delta, Z, taus, tun$lambda)
  if (f$status != "ok") return(fail(f$status, tau_fail = taus[f$k]))
  list(status = "ok", selected = which(abs(f$B[-1L, m]) > coef_tol), failed = FALSE,
       mult = tun$mult)
}

## L-HDCQR as specified by Zheng, Peng and He (2018, Section 5.1) -------------
## Sequential L1-penalized estimation on the tau grid with a penalty that grows
## along the grid:
##   * lambda_0: conservative pivotal choice at tau_0 (Belloni and
##     Chernozhukov, 2011);
##   * k = 1, ..., m: lambda_k = lambda_{k-1} + increment, the increment chosen
##     by 5-fold CV from a grid c * sqrt(n log(p v n)) * eps_n, eps_n the grid
##     spacing; the held-out error at tau_k is, by default (crit = "ph"), the
##     held-out Peng-Huang objective at tau_k (Zheng et al. do not state the
##     CV criterion); crit = "deviance" uses the sum of absolute deviance
##     residuals with martingale residual M_i = N_i - sum_r 1{at risk} dH - tau_0
##     (Fei et al., 2021, Remark 3), which over-penalizes in our experience;
##   * the compensator starts at tau_0, as in their equation (7).
## The constants (their a, b) are not reported; the increment grid below spans
## them on the log scale and includes 0.  Each LAD-lasso step is solved at
## level 0.5 of the augmented problem (their Section 3.1).
z18_state <- function(y, delta, Z, rows, a0) {
  ev <- rows[delta[rows] == 1]
  list(rows = rows, Ze = Z[ev, , drop = FALSE], ye = y[ev],
       zminus = -colSums(Z[ev, , drop = FALSE]), a = rep(a0, length(rows)), b = NULL)
}
z18_solve <- function(st, Z, lam, R0 = 1e4, R_growth = 10, R_tries = 5L) {
  P <- ncol(Z)
  zplus <- 2 * colSums(Z[st$rows, , drop = FALSE] * st$a)
  pen <- c(0, rep(lam, P - 1L))
  init <- if (is.null(st$b)) integer(0) else which(st$b != 0)
  R <- R0
  for (t in seq_len(R_tries)) {
    fit <- tryCatch(
      if (P > 100L) lad_lasso_ws(rbind(st$Ze, st$zminus, zplus), c(st$ye, R, R), pen, init = init)
      else lad_lasso(rbind(st$Ze, st$zminus, zplus), c(st$ye, R, R), pen),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    b <- fit$coef
    if (all(R - c(sum(st$zminus * b), sum(zplus * b)) > 0)) return(b)
    R <- R * R_growth
  }
  NULL
}
## advance a state to the next grid point after committing b at the current one
z18_advance <- function(st, y, Z, b, dH) {
  st$a <- st$a + at_risk(y[st$rows], drop(Z[st$rows, , drop = FALSE] %*% b)) * dH
  st$b <- b
  st
}
lhdcqr_z18 <- function(y, delta, X, taus, K = 5L, foldid = NULL, std = "sd",
                       inc_grid = c(0, exp(seq(log(0.02), log(20), length.out = 12))),
                       nsim = 500L, coef_tol = 1e-6, cores = 1L,
                       crit = c("ph", "deviance"), ...) {
  crit <- match.arg(crit)
  n <- length(y); m <- length(taus)
  Z <- standardize_by(X, seq_len(n), std)
  p <- ncol(Z) - 1L
  fold <- if (is.null(foldid)) make_folds(delta, K) else match(foldid, sort(unique(foldid)))
  K <- max(fold)
  dH <- diff(H(taus)); a0 <- taus[1L]
  unit <- sqrt(n * log(max(p, n))) * (if (m > 1L) taus[2L] - taus[1L] else 0)
  lam <- numeric(m)
  lam[1L] <- pivotal_lambda(Z[, -1L, drop = FALSE], taus[1L], nsim)
  full <- z18_state(y, delta, Z, seq_len(n), a0)
  tr <- lapply(seq_len(K), function(k) z18_state(y, delta, Z, which(fold != k), a0))
  ## held-out compensators (one vector per fold, over the held-out rows)
  te_rows <- lapply(seq_len(K), function(k) which(fold == k))
  te_comp <- lapply(te_rows, function(r) rep(a0, length(r)))
  dev_err <- function(k, b, kk) {               # held-out error of fold kk at grid point k
    r <- te_rows[[kk]]; fit <- drop(Z[r, , drop = FALSE] %*% b)
    if (crit == "ph")                              # held-out Peng-Huang objective at tau_k
      return(sum(delta[r] * abs(y[r] - fit) + (delta[r] - 2 * te_comp[[kk]]) * fit))
    M <- delta[r] * had_event_by(y[r], fit) - te_comp[[kk]]
    D <- -2 * (M + delta[r] * log(pmax(delta[r] - M, 1e-300)))
    sum(sqrt(pmax(D, 0)))
  }
  B <- matrix(NA_real_, ncol(Z), m)
  choice <- rep(NA_real_, m)
  for (k in seq_len(m)) {
    if (k == 1L) {
      cand <- lam[1L]; sols <- lapply(tr, function(st) list(z18_solve(st, Z, cand)))
      best <- 1L
    } else {
      cand <- lam[k - 1L] + inc_grid * unit
      jobs <- expand.grid(f = seq_len(K), c = seq_along(cand))
      run <- function(j) z18_solve(tr[[jobs$f[j]]], Z, cand[jobs$c[j]])
      res <- if (cores > 1L) parallel::mclapply(seq_len(nrow(jobs)), run, mc.cores = cores,
                                                mc.set.seed = FALSE)
             else lapply(seq_len(nrow(jobs)), run)
      err <- vapply(seq_along(cand), function(cc) {
        e <- vapply(seq_len(K), function(f) {
          b <- res[[which(jobs$f == f & jobs$c == cc)]]
          if (is.null(b)) Inf else dev_err(k, b, f) }, 0)
        sum(e) }, 0)
      if (all(!is.finite(err))) return(fail("no_root", tau_fail = taus[k]))
      best <- which(err == min(err))[1L]           # smallest increment on ties
      sols <- lapply(seq_len(K), function(f) list(res[[which(jobs$f == f & jobs$c == best)]]))
      choice[k] <- inc_grid[best]
    }
    lam[k] <- cand[best]
    bf <- z18_solve(full, Z, lam[k])
    if (is.null(bf)) return(fail("no_root", tau_fail = taus[k]))
    B[, k] <- bf
    if (k < m) {
      full <- z18_advance(full, y, Z, bf, dH[k])
      for (f in seq_len(K)) {
        bt <- sols[[f]][[1L]]
        if (is.null(bt)) bt <- bf                    # failed fold: continue from the full fit
        te_comp[[f]] <- te_comp[[f]] + at_risk(y[te_rows[[f]]], drop(Z[te_rows[[f]], , drop = FALSE] %*% bt)) * dH[k]
        tr[[f]] <- z18_advance(tr[[f]], y, Z, bt, dH[k])
      }
    }
  }
  list(status = "ok", selected = which(abs(B[-1L, m]) > coef_tol), failed = FALSE,
       lambda = lam, increment = choice, mult = mean(choice, na.rm = TRUE))
}

## SCQR-Lasso / SCQR-SCAD of He et al. (author C++ code) ----------------------
.scqr_env <- new.env()
## The authors' file is not redistributed in the public repository (their
## repository has no license); if it is missing it is downloaded from the
## authors' repository at a fixed commit, and the C++11 plugin line is changed
## to C++17 (RcppArmadillo >= 14 requires C++14 or later).  Nothing else changes.
HDSCQR_URL <- paste0("https://raw.githubusercontent.com/XiaoouPan/scqr/",
                     "e69daceb4509fa6e79f81b65c801ad10bc95f207/src/hdscqr.cpp")
fetch_hdscqr <- function(path, url = HDSCQR_URL) {
  src <- readLines(url, warn = FALSE)
  i <- grep("Rcpp::plugins(cpp11)", src, fixed = TRUE)
  if (length(i) != 1L) stop("unexpected content of ", url)
  src <- append(src[-i], c("// (changed from cpp11: RcppArmadillo >= 14 requires C++14 or later)",
                           "// [[Rcpp::plugins(cpp17)]]"), after = i - 1L)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(src, path)
  message("downloaded hdscqr.cpp (He, Pan, Tan and Zhou, 2022) to ", path)
  invisible(path)
}
load_scqr <- function(path = file.path(CQR_COMP_ROOT, "vendor", "hdscqr.cpp")) {
  if (is.null(.scqr_env$f)) {
    if (!file.exists(path)) fetch_hdscqr(path)
    e <- new.env()
    Rcpp::sourceCpp(path, env = e)
    .scqr_env$f <- e
  }
  .scqr_env$f
}
scqr <- function(y, delta, X, taus, penalty = c("lasso", "scad"), K = 5L,
                 lambda = exp(seq(log(0.01), log(0.20), length.out = 50L)),
                 coef_tol = 1e-6, foldid = NULL, ...) {
  penalty <- match.arg(penalty)
  e <- load_scqr()
  n <- length(y); p <- ncol(X); m <- length(taus)
  fold <- if (is.null(foldid)) make_folds(delta, K) else as.integer(foldid)
  K <- max(fold)
  h <- 0.5 * (log(p) / n)^(1 / 4)
  r <- tryCatch(
    if (penalty == "lasso") e$cvSqrLassoGrow(X, delta, y, lambda, fold, taus, K, h)$beta
    else e$cvSqrScadGrow(X, delta, y, lambda, fold, taus, K, h),
    error = function(err) err)
  if (inherits(r, "error")) return(fail("scqr_error", message = conditionMessage(r)))
  b <- as.matrix(r)[-1L, m]
  if (any(!is.finite(b))) return(fail("scqr_nonfinite"))
  list(status = "ok", selected = which(abs(b) > coef_tol), failed = FALSE)
}

## DS: single data splitting of Dai et al. (2023) --------------------------------
## Their Algorithm 1 transferred to censored quantile regression: half-half
## split; lasso on the first half with the penalty chosen by 10-fold CV (their
## simulations); the candidate set is the CV support, without a cap; the
## unpenalized estimator on the second half restricted to that set (their OLS
## step); mirror statistics; cutoff #{M_j < -t} / max(#{M_j > t}, 1) <= q
## without the +1.  The CV folds are drawn within the first half (Dai et al.
## tune the lasso on that half only), so the shared full-sample partition is
## not used.
ds <- function(y, delta, X, taus, q = 0.2, K = 10L, loss = "ph", cores = 1L, ...) {
  ssms(y, delta, X, taus, q = q, offset = 0, lambda_rule = "cv", K = K, loss = loss,
       screen_cap = NULL, screen_rule = "cv", cores = cores, foldid = NULL)
}

## Cox-lasso (glmnet), used as a reference in the real-data analyses -----------
## time on the original scale; penalty chosen by K-fold CV (lambda.min).
cox_lasso <- function(time, delta, X, K = 5L, foldid = NULL) {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("install.packages('glmnet')")
  extra <- if ("cox.ties" %in% names(formals(glmnet::glmnet))) list(cox.ties = "breslow") else list()
  args <- c(list(x = X, y = survival::Surv(time, delta), family = "cox", nfolds = K), extra)
  if (!is.null(foldid)) args$foldid <- foldid
  cv <- do.call(glmnet::cv.glmnet, args)
  b <- as.numeric(stats::coef(cv, s = "lambda.min"))
  list(status = "ok", selected = which(b != 0), failed = FALSE, beta = b, lambda = cv$lambda.min)
}
