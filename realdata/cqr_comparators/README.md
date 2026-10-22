# Variable selection for censored quantile regression: stand-alone method code

The code of the simulation study (`simulation/R`), packaged so that every method can be called on any data set with the settings of the paper. It contains the proposed methods (SS-MS, MS-MS) and all comparators (DS, MDS, L-HDCQR, SCQR-Lasso, SCQR-SCAD, Cox-lasso). Each comparator follows its original publication, including its tuning.

## Use

```r
source("cqr_comparators/load.R")          # loads the solver and all methods
taus <- seq(0.1, 0.5, length.out = 40)    # quantile grid; the last point is the target tau
fit  <- cqr_msms(y, delta, X, taus, q = c(0.1, 0.2))
cqr_selected(fit, q = 0.1)                # selected column indices
```

- `y` = log(min(T, C)), `delta` = event indicator (1 = event), `X` = n × p numeric matrix (standardized internally).
- Every wrapper returns a list with `status`, `selected` (column indices) and `failed`; with several q there is also `by_q`.
- `Rscript demo.R` (from this folder) runs every method on a small simulated data set (a few minutes) and prints the number selected, FDP and power; `demo_output_example.txt` is a reference output.

Needs `survival`, `quantreg`, `parallel`; SCQR needs `Rcpp`, `RcppArmadillo` and a C++ compiler (the C++ is compiled at the first call, about a minute); Cox-lasso needs `glmnet`.

## Methods

| Function | Method | Procedure | Reference |
|---|---|---|---|
| `cqr_ssms` | SS-MS (proposed) | half split; lasso with the pivotal penalty × a multiplier chosen by 5-fold CV on the first half; path-ordered screen capped at ⌊d₂/log d₂⌋; unpenalized refit on the second half; mirror statistics, cutoff U₊ (offset 1) | this paper |
| `cqr_msms` | MS-MS (proposed) | SS-MS on B = 50 splits; multiplier tuned once; inclusion-rate aggregation of Dai et al.; ties broken at random | this paper |
| `cqr_ds` | DS | half split; lasso tuned by 10-fold CV within the first half; candidate set = CV support (no cap); Peng–Huang refit on the second half (in place of their OLS refit); cutoff U (no +1) | Dai, Lin, Xing and Liu (2023, JASA) |
| `cqr_mds` | MDS | 50 complete DS fits, each with its own 10-fold CV; inclusion-rate aggregation with their rule {j : I_j > I_(ℓ)} | same |
| `cqr_lhdcqr` | L-HDCQR | sequential L1-penalized estimation on the grid; λ₀ pivotal at τ₀; at every grid point an increment chosen by 5-fold CV from {0} ∪ 12 log-spaced values × √(n log(p∨n))·εₙ; each step solved as an augmented median (LAD) problem | Zheng, Peng and He (2018, AoS), Section 5.1 |
| `cqr_scqr_lasso`, `cqr_scqr_scad` | SCQR | the authors' C++ code; 50 log-spaced penalties in [0.01, 0.2]; bandwidth 0.5 (log p / n)^{1/4}; the authors' CV | He, Pan, Tan and Zhou (2022, AoS) |
| `cqr_cox_lasso` | Cox-lasso | glmnet, 5-fold CV, lambda.min; takes time on the original scale | reference (not a quantile method) |

## Departures from the original publications

- **DS / MDS.** Dai et al. refit by OLS; here the refit is the unpenalized Peng–Huang estimator, its censored-data counterpart. Everything else (10-fold CV, no cap, cutoff U, aggregation rule) is as published.
- **L-HDCQR.** Zheng et al. do not state the CV criterion. The default is the held-out Peng–Huang objective (`crit = "ph"`); `crit = "deviance"` uses the deviance residual of Fei et al. (2021, Remark 3), which over-penalizes in our experience. The authors' routine `quantproc()` (`reference/zheng2018_author_quantproc.R`) calls `rq(..., tau = JJ[s])` at every grid point, but the augmented problem must be solved at τ = 0.5; `reference/check_author_quantproc.R` shows this on uncensored data, where the correct path must equal ordinary quantile regression.
- **SCQR.** `vendor/hdscqr.cpp` is the file of the authors' repository (github.com/XiaoouPan/scqr) with one change: `Rcpp::plugins(cpp11)` → `cpp17` (current RcppArmadillo needs C++14 or later). The repository has no license file, so the public GitHub copy of this code does not contain the file: `load_scqr()` downloads it at the first SCQR call from the authors' repository at a fixed commit and makes the same change (`fetch_hdscqr()` in `R/methods.R`). Please cite He et al. (2022) when using it.

## Files

```
load.R              loads the code and defines the cqr_* wrappers (settings of the paper)
demo.R              small example with every method
R/solver.R          Peng–Huang path solver (augmented LAD-lasso, interior point), pivotal penalty, CV
R/methods.R         ssms(), mds(), ds(), lhdcqr_z18(), scqr(), cox_lasso() and other low-level functions
R/dgp.R             simulation models (location and heteroscedastic)
R/utils.R           seeds, folds, selection_metrics() and other helpers
vendor/hdscqr.cpp   SCQR author C++ code (downloaded at first use if missing)
reference/          the L-HDCQR author routine and the check described above (the routine is not in the public GitHub copy)
```

The low-level functions take more arguments than the wrappers, e.g. `mds(tune = "pivotal")` uses a fixed multiplier and `ssms(screen_rule = "cv")` screens with the CV support only; these were used for the ablations. The R package `cqrfdr` (0.3.0) gives the same results through `ssms()`, `msms()`, `ds()`, `mds()` and `hdcqr(lambda = "zheng")`, with the interface `f(x, y, delta, tau = ...)`.
