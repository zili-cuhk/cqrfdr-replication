## ---------------------------------------------------------------------------
## Simulation configuration: edit this file only.
## ---------------------------------------------------------------------------
## Two scenario sets, chosen with the environment variable CQR_SIM_SET:
##   main      (default) the simulation reported in the paper
##   df_pilot  short pilot to fix the degrees of freedom of the t covariates
## Results go to results/<set>/, calibration to sim/calibration_<set>.csv.

SIM_SET <- Sys.getenv("CQR_SIM_SET", "main")

## Degrees of freedom of the heavy-tailed covariates / errors in the main set.
## Chosen after the df_pilot (t3, t5, t10; 20 replications): t5.
XDF <- 5
EDF <- 5
## Signal-strength grid of the "extra" set (A1); the baseline 0.8 is in "main"
BETA_GRID <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.65)

SIM_CONFIG <- list(
  master_seed = 20261014L,
  reps        = 100L,          # replications per scenario
  q           = c(0.1, 0.2),   # nominal FDR levels (SS-MS and MDS report one row per level)
  tau_lower   = 0.1,           # first grid point nu
  n_grid      = 40L,           # grid points from nu to the target tau
  methods     = c("ssms", "ssms_size", "ssms_path", "ds", "msms", "mds", "mds_cv", "mds_dai", "lhdcqr", "lhdcqr_z18", "scqr_lasso", "scqr_scad"),
  cv_folds    = 5L,            # same K and the same (event-stratified) fold partition for every method
  cores       = as.integer(Sys.getenv("LSB_DJOB_NUMPROC", "1")),  # LSF sets this from -n

  ## proposed single split (Algorithm 1)
  ssms   = list(lambda_rule = "cv", loss = "ph", screen_cap = "theory", offset = 1),
  ## multiple splitting (Dai et al. 2023); fixed pivotal penalty keeps 50 splits affordable
  mds    = list(B = 50L, tune = "pivotal", ties = "dai", fixed_mult = 0.5, offset = 0,
                screen_cap = "theory"),
  ## corrected MDS (replaces "mds" in the paper): penalty multiplier chosen by
  ## the same 5-fold CV as SS-MS (once per data set), ties broken at random
  mds_cv = list(B = 50L, tune = "cv_once", ties = "random", offset = 0,
                screen_cap = "theory"),
  ## single data splitting of Dai et al. (2023), their tuning: 10-fold CV on
  ## the first half, CV support without a cap, cutoff without the +1
  ds     = list(K = 10L, loss = "ph"),
  ## MDS with the design of Dai et al. (2023): B = 50 complete DS fits, each
  ## with its own 10-fold CV on its first half, CV support without a cap,
  ## cutoff U, their aggregation rule including ties
  mds_dai = list(B = 50L, tune = "cv_each", ties = "dai", offset = 0, K = 10L,
                 screen_cap = "none", screen_rule = "cv", own_folds = TRUE),
  ## MS-MS: SS-MS (path-ordered screen, cap, cutoff U+) on B = 50 splits,
  ## aggregated by inclusion rates as MDS aggregates DS; multiplier by the
  ## same 5-fold CV as SS-MS, chosen once per data set; ties broken at random
  msms   = list(B = 50L, tune = "cv_once", ties = "random", offset = 1,
                screen_cap = "theory", screen_rule = "path"),
  ## comparators
  lhdcqr = list(lambda_rule = "cv", loss = "ph"),
  ## L-HDCQR as specified by Zheng, Peng and He (2018, Section 5.1): penalty
  ## increased along the grid, each increment chosen by 5-fold CV of the
  ## held-out Peng-Huang objective (crit = "ph"); lambda_0 pivotal at tau_0
  lhdcqr_z18 = list(),
  scqr_lasso = list(),
  scqr_scad  = list(),
  ## SS-MS with the size-calibrated screen (d_pilot only)
  ssms_size = list(lambda_rule = "cv", loss = "ph", screen_cap = "theory", offset = 1),
  ## SS-MS with the path-ordered screen (candidates ranked by entry along the
  ## lasso path below the CV multiplier, up to the cap floor(d2 / log d2))
  ssms_path = list(lambda_rule = "cv", loss = "ph", screen_cap = "theory", offset = 1)
)

## Scenario tables ------------------------------------------------------------
## Common: n = 800, s = 10 active covariates (first 10), Toeplitz rho = 0.5,
## target tau = 0.5, normal errors (sd 0.5), log-normal censoring calibrated
## to the target rate, unless a scenario says otherwise.
##
## main (27 scenarios; n = 800 throughout)
##   A  full factorial  p {100, 1000} x censoring {10%, 30%, 50%}
##                      x covariates {Gaussian, multivariate t_XDF}       (12)
##      baseline = A_p1000_c30_gauss
##   B  one factor at a time around the baseline:
##      tau {0.25, 0.4, 0.6, 0.7}; rho {0.1, 0.3, 0.7};
##      errors t_EDF and centred exponential                                (9)
##   D  heteroscedastic model (covariate 11 null at the median, active at
##      other quantiles), tau {0.25, 0.5, 0.6}                              (3)
## extra (11 scenarios, 100 replications; added after review of the main set)
##   S_beta*  signal strength: coefficient function scaled by b / 0.8 (A1)
##   S_s*     number of active covariates s = 5, 20, 30 (B1)
## df_pilot (6 scenarios, 20 replications)
##      p {100, 1000} x covariates t_3, t_5, t_10, 30% censoring
make_scenarios <- function(set = SIM_SET) {
  base <- list(n = 800L, p = 1000L, rho = 0.5, tau = 0.5, xdist = "gaussian",
               edist = "normal", cens = 0.30, model = "location",
               xdf = XDF, edf = EDF, err_scale = "sd", sigma = 0.5, gamma = 3)
  if (set == "extra") base <- c(base, list(s = 10L, beta_level = 0.8, beta_slope = 0.8))
  rows <- list()
  add <- function(.label, ...) {
    r <- base; a <- list(...); r[names(a)] <- a
    rows[[length(rows) + 1L]] <<- as.data.frame(c(list(scenario = .label), r), stringsAsFactors = FALSE)
  }
  if (set == "main") {
    for (p in c(100L, 1000L)) for (cr in c(0.1, 0.3, 0.5)) for (xd in c("gaussian", "t"))
      add(sprintf("A_p%d_c%02d_%s", p, round(100 * cr), if (xd == "t") paste0("t", XDF) else "gauss"),
          p = p, cens = cr, xdist = xd)
    for (t in c(0.25, 0.4, 0.6, 0.7)) add(sprintf("B_tau%03d", round(100 * t)), tau = t)
    for (r in c(0.1, 0.3, 0.7)) add(sprintf("B_rho%02d", round(10 * r)), rho = r)
    add(sprintf("B_err_t%d", EDF), edist = "t")
    add("B_err_exp", edist = "exp")
    ## D  (supplementary, low signal-to-noise) sigma 0.5, gamma 3: mean noise
    ##    scale 1.25 = 2.5 x the location model
    for (t in c(0.25, 0.5, 0.6)) add(sprintf("D_hetero_tau%03d", round(100 * t)), model = "hetero", tau = t)
    ## Dm (main heteroscedastic block, signal-to-noise matched) sigma 0.25,
    ##    gamma 4: mean noise scale 0.75; coefficient of Z11 is q(tau):
    ##    -0.674, 0, +0.674 at tau = 0.25, 0.5, 0.75.  Appended at the end so
    ##    the task ids of the first 24 scenarios do not change.
    for (t in c(0.25, 0.5, 0.75)) add(sprintf("Dm_hetero_tau%03d", round(100 * t)), model = "hetero",
                                      tau = t, sigma = 0.25, gamma = 4)
  } else if (set == "df_pilot") {
    for (p in c(100L, 1000L)) for (df in c(3, 5, 10))
      add(sprintf("P_p%d_t%d", p, df), p = p, xdist = "t", xdf = df)
  } else if (set == "d_pilot") {
    ## heteroscedastic model: original (sigma 0.5, gamma 3; mean noise scale
    ## 1.25 = 2.5 x the location model) vs signal-to-noise matched
    ## (sigma 0.25, gamma 4; mean scale 0.75, same Z11 coefficient q(tau))
    add("D_hetero_tau025", model = "hetero", tau = 0.25)
    for (t in c(0.25, 0.5, 0.6))
      add(sprintf("Dm_hetero_tau%03d", round(100 * t)), model = "hetero", tau = t, sigma = 0.25, gamma = 4)
  } else if (set == "extra") {
    ## A1  signal strength: beta_j(u) = b {sign_j + (u - 0.5) slope_j}, i.e. the
    ##     baseline coefficient function scaled by b / 0.8; |beta_j(0.5)| = b,
    ##     baseline otherwise (p = 1000, 30% censoring, Gaussian, tau = 0.5)
    for (b in BETA_GRID) add(sprintf("S_beta%03d", round(1000 * b)), beta_level = b, beta_slope = b)
    ## B1  number of active covariates (strong signals, |beta_j| = 0.8)
    for (ss in c(5L, 20L, 30L)) add(sprintf("S_s%02d", ss), s = ss)
  } else stop("unknown CQR_SIM_SET: ", set)
  do.call(rbind, rows)
}
if (SIM_SET == "df_pilot") SIM_CONFIG$reps <- 20L
if (SIM_SET == "extra") SIM_CONFIG$methods <- c("ssms", "ssms_size", "ssms_path", "ds", "msms", "mds_cv", "mds_dai", "lhdcqr", "lhdcqr_z18", "scqr_lasso", "scqr_scad")
if (SIM_SET == "d_pilot") {
  SIM_CONFIG$reps <- 20L
  SIM_CONFIG$methods <- c("ssms", "ssms_size", "scqr_scad")
}

## Smoke mode: CQR_SIM_SMOKE=1 shrinks everything so the whole pipeline
## (calibrate -> run -> collect) finishes in a few minutes on a laptop.
if (identical(Sys.getenv("CQR_SIM_SMOKE"), "1")) {
  SIM_CONFIG$reps <- 2L; SIM_CONFIG$n_grid <- 15L
  SIM_CONFIG$cv_folds <- 3L; SIM_CONFIG$ds$K <- 3L; SIM_CONFIG$mds$B <- 3L; SIM_CONFIG$mds_cv$B <- 3L; SIM_CONFIG$msms$B <- 3L; SIM_CONFIG$mds_dai$B <- 3L; SIM_CONFIG$mds_dai$K <- 3L
  .full_scenarios <- make_scenarios
  make_scenarios <- function(set = SIM_SET) {
    s <- .full_scenarios(set)
    s <- s[s$scenario %in% c("A_p1000_c30_gauss", "A_p100_c50_t5", "B_tau025", "D_hetero_tau025",
                             "P_p100_t3", "P_p1000_t10", "Dm_hetero_tau025"), ]
    s$n <- 300L; s$p <- pmin(s$p, 80L)
    s
  }
}

results_dir <- function(root) file.path(root, "results", SIM_SET)
calibration_file <- function(root) file.path(root, "sim", sprintf("calibration_%s.csv", SIM_SET))
