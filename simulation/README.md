# Simulation study (Section 4 and Supplement)

## Layout

```
R/          dgp.R      data-generating models, censoring calibration
            solver.R   Peng–Huang path (penalized / unpenalized), LP solvers
            methods.R  SS-MS, MS-MS, DS, MDS, L-HDCQR, SCQR wrappers
            utils.R    seeds, folds, mirror statistics, cutoffs, load_sources()
sim/        config.R   scenarios and method settings (the only file to edit)
            calibrate.R, run.R, submit_lsf.sh, collect.R          main pipeline
            stability.R, submit_stability.sh, stability_collect.R split-to-split stability
            paper_outputs.R, paper_outputs_extra.R, paper_table_stability.R  paper tables / figures
vendor/     hdscqr.cpp  He, Pan, Tan and Zhou (2022) author code for SCQR-Lasso / SCQR-SCAD
            (downloaded at the first SCQR call if missing; see the top-level README)
tests/      test_core.R            unit and end-to-end checks (about a minute)
            check_author_lhdcqr.R  check of the L-HDCQR author routine (see below)
            ref/                   author functions used by that check
results/    summaries of the runs reported in the paper (see "Stored results")
```

All scripts are run from this folder (`simulation/`), or with `CQR_SIM_ROOT` set to it.

## Quick check

```bash
Rscript tests/test_core.R
```
Every line should start with `PASS`. Needs `quantreg`, `survival`, `glmnet`; SCQR needs `Rcpp`, `RcppArmadillo` and a C++ compiler.

## Scenario sets

`CQR_SIM_SET` chooses the set (defined in `sim/config.R`, `make_scenarios()`). Common settings: n = 800, s = 10, Toeplitz ρ = 0.5, τ = 0.5, normal errors, log-normal censoring calibrated to the target rate, 100 replications, q ∈ {0.1, 0.2}, 40 grid points on [0.1, τ], 5-fold CV.

| Set | Scenarios | Paper |
|---|---|---|
| `main` | A: p ∈ {100, 1000} × censoring {10, 30, 50}% × Gaussian / t₅ covariates (12); B: τ ∈ {0.25, 0.4, 0.6, 0.7}, ρ ∈ {0.1, 0.3, 0.7}, t₅ and exponential errors (9); Dm: heteroscedastic, τ ∈ {0.25, 0.5, 0.75} (3); D: heteroscedastic with larger noise (3, Supplement) | Tables 1–3, Figures 1–2; Supplement (larger noise scale) |
| `extra` | signal strength b ∈ {0.1, …, 0.65} (8); s ∈ {5, 20, 30} (3) | Table 4, Figure 3; Supplement (screen) |
| `df_pilot` | t₃, t₅, t₁₀ covariates, p ∈ {100, 1000} (20 replications) | pilot that fixed df = 5; not reported |

## Rerunning (LSF)

```bash
bash sim/submit_lsf.sh                          # main: calibration (once), then the job array
CQR_SIM_SET=extra bash sim/submit_lsf.sh        # extra set
CQR_SIM_ONLY=msms bash sim/submit_lsf.sh        # add or re-run one method in existing results
Rscript sim/collect.R                           # results/main/summary.csv (and replications.csv)
CQR_SIM_SET=extra Rscript sim/collect.R         # results/extra/summary.csv

bash sim/submit_stability.sh                    # stability, SS-MS with the CV screen
CQR_STAB_RULE=path bash sim/submit_stability.sh # stability, SS-MS with the path-ordered screen
Rscript sim/stability_collect.R
CQR_STAB_RULE=path Rscript sim/stability_collect.R
```

- One task = one replication of one scenario with all methods (`run.R <task_id>`); 4 cores, 24 h limit, 8 GB.
- `calibrate.R` writes `sim/calibration_<set>.csv` (censoring location per scenario); `submit_lsf.sh` runs it first if the file is missing.
- Seeds are derived from `master_seed`, scenario, replication and method, so each method has its own random stream: re-running one method (`CQR_SIM_ONLY`) or resubmitting after a crash gives the same numbers as an uninterrupted run. Finished replications are skipped.
- Without LSF: loop `Rscript sim/run.R $t` over task ids 1 … (scenarios × 100), after `Rscript sim/calibrate.R`.
- Cost: about 38,000 core-hours for `main` + `extra`, about 33,000 of them for MDS (`mds_dai`, which runs 10-fold CV in each of its 50 splits).

## Methods (`base_method` column of summary.csv; `method` adds `_q0.1` / `_q0.2` for the FDR methods)

| Code | Paper | Settings (`sim/config.R`) |
|---|---|---|
| `ssms_path` | SS-MS | half split; CV multiplier of the pivotal penalty; path-ordered screen capped at ⌊d₂/log d₂⌋; unpenalized refit; cutoff U₊ |
| `msms` | MS-MS | SS-MS on B = 50 splits; multiplier chosen once; inclusion-rate aggregation; random ties |
| `ds` | DS (Dai et al., 2023) | 10-fold CV on the first half; CV support without a cap; cutoff U |
| `mds_dai` | MDS (Dai et al., 2023) | B = 50 DS fits, each with its own 10-fold CV; their aggregation and tie rule |
| `lhdcqr_z18` | L-HDCQR (Zheng, Peng and He, 2018) | pivotal λ₀ at τ₀; increment at each grid point chosen by 5-fold CV |
| `scqr_lasso`, `scqr_scad` | SCQR-Lasso, SCQR-SCAD (He et al., 2022) | author C++ code, `vendor/hdscqr.cpp` |
| `ssms` | SS-MS with the CV screen | Figure 3 and the screen table of the Supplement |

`ssms_size`, `mds`, `mds_cv` and `lhdcqr` are development variants kept in the code; they are not reported. `tests/check_author_lhdcqr.R` shows why L-HDCQR is our implementation rather than the author routine `quantproc`: that routine solves the augmented problem at τ = τ_s instead of 0.5, so it does not reproduce ordinary quantile regression without censoring.

## Stored results

| File | Contents |
|---|---|
| `results/main/summary.csv` | one row per scenario × method (× q): FDR and power with Monte Carlo SEs and intervals, mean selection size, screen recall and size, failures, median run time |
| `results/extra/summary.csv` | same for the `extra` set |
| `results/stability/stability_{summary,by_data}.csv` | split-to-split Jaccard, FDP, power (CV screen) |
| `results/stability_path/stability_{summary,by_data}.csv` | same, path-ordered screen |

The per-replication `.rds` files (several GB) are not included; `collect.R` rebuilds the summaries from them.

## Tables and figures

```bash
Rscript sim/paper_outputs.R       results/main/summary.csv results/paper
Rscript sim/paper_outputs_extra.R results/extra/summary.csv results/main/summary.csv results/paper
Rscript sim/paper_table_stability.R results results/paper
```
Written to `results/paper/`: `tables_main.tex`, `table_supp_lowsnr.tex`, `fig_factorial.pdf`, `fig_onefactor.pdf` (main set); `table_sparsity.tex`, `table_screen.tex`, `fig_signal.pdf` (extra set); `table_stability.tex`. New columns relative to the first submission (MS-MS, DS, MDS, L-HDCQR) are wrapped in `\rev{}`; define `\newcommand{\rev}[1]{#1}` to print them normally.
