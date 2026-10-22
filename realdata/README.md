# Real data analysis (Section 5 and Supplement)

## Layout

```
R/config.R             settings and the pre-specified eligibility rules (fixed before any selection was seen)
R/00_install.R         installs cqrfdr and the CRAN / Bioconductor packages
R/01_download.R        downloads the gene-expression cohorts (login node, needs internet)
R/02_prepare.R         eligibility, outcome-blind gene filter, task table tasks.csv
R/05_targeted.R        METABRIC targeted design (PAM50 + Oncotype DX genes + 10 clinical factors)
R/06_glioma.R          TCGA diffuse glioma (LGG + GBM)
R/08_support.R         SUPPORT design (main effects, squares, pairwise products)
R/10_add_tasks.R       appends task types (MS-MS, MDS, comparators, decoy / hold-out) to a task table
R/03_task.R            runs one task (one row of a task table)
R/04_collect.R         collects the task results into tables
R/07_signal.R          method-free signal strength of each cohort (univariate Cox BH, ridge-Cox C-index)
R/09_requirements.R    requirement checks A–D for every data set and task completion
R/11_paper_tables.R    all tables and numbers of Section 5 from results_paper/
plasmode/              plasmode simulation on the real designs (plasmode.R, submit_plasmode.sh, plasmode_collect.R)
cqr_comparators/       stand-alone code of L-HDCQR, DS and SCQR (the simulation code; see its README)
submit_realdata.sh     LSF job array over a task table
results_paper/         collected results of the runs reported in the paper
```

All scripts are run from this folder (`realdata/`), or with `CQR_RD_ROOT` set to it.

## Tables from the stored results

```bash
Rscript R/11_paper_tables.R results_paper
```
Prints Tables 1–3 of Section 5, the METABRIC (targeted) numbers in the text, the supplementary plasmode table (three designs) and the genome-wide cohorts table, and writes them as CSV files to `results_paper/paper_tables/`. It needs base R only and takes a few seconds.

| Folder in `results_paper/` | Produced by | Contents |
|---|---|---|
| `tables_support/` | `04_collect.R` on `tasks_support.csv` | SUPPORT: selections, run-to-run agreement (`multi_split.csv`), decoys, hold-out C-index, design |
| `tables_targeted/` | `04_collect.R` on `tasks_targeted.csv` | METABRIC targeted |
| `tables_glioma/` | `04_collect.R` on `tasks_glioma.csv` | TCGA glioma |
| `tables_genome/` | `02_prepare.R`, `04_collect.R`, `07_signal.R`, `09_requirements.R` on `tasks.csv` | 13 genome-wide cohorts: eligibility, requirements, signal strength, selections |
| `plasmode/` | `plasmode/plasmode_collect.R` | FDR and power (`plasmode_summary_<design>.csv`) and screen size and recall (`plasmode_screen_<design>.csv`) |

## Rerunning the analyses (LSF)

### 1. Software and data (login node)

```bash
module load R/4.5.1
Rscript R/00_install.R ../package/cqrfdr_0.3.0.tar.gz   # or: Rscript R/00_install.R github
Rscript R/01_download.R          # all gene-expression sources; or e.g. "Rscript R/01_download.R metabric tcga"
```
`01_download.R` writes `data/raw/<name>.rds` (expression matrix, time in years, status) for METABRIC (cBioPortal), the TCGA cohorts (UCSC Xena, TCGA-CDR survival), the ovarian cohorts (`curatedOvarianData`) and DLBCL (GEO GSE10846). Re-running fetches only what is missing.

**SUPPORT** is downloaded by `08_support.R` from the UCI repository or hbiostat.org. If the cluster cannot reach them, download `support2csv.zip` from https://hbiostat.org/data (or the UCI repository, data set 880) and unzip it so that `data/download/support2/support2.csv` exists.

### 2. Designs and task tables

```bash
Rscript R/02_prepare.R           # genome-wide cohorts -> tasks.csv
Rscript R/05_targeted.R          # METABRIC targeted   -> tasks_targeted.csv
Rscript R/06_glioma.R            # TCGA glioma         -> tasks_glioma.csv
Rscript R/08_support.R           # SUPPORT             -> tasks_support.csv
Rscript R/10_add_tasks.R                                   # msms (2 runs), mds_dai, all tables
Rscript R/10_add_tasks.R cmp,decoy_z18,holdout_z18 SUPPORT2,BRCA_METABRIC_TARGETED
Rscript R/10_add_tasks.R decoy_msms,holdout_msms   SUPPORT2,BRCA_METABRIC_TARGETED
Rscript R/09_requirements.R      # checks A-D; PASS / FAIL per data set
```
`10_add_tasks.R` prints the id range of each new type, for `ONLY=` below. New rows are appended, so the ids of existing tasks do not change.

### 3. Jobs

```bash
TASKS=tasks_support.csv bash submit_realdata.sh                    # all SUPPORT tasks
TASKS=tasks_support.csv CORES=10 ONLY="<mds_dai ids>" bash submit_realdata.sh
TASKS=tasks_targeted.csv bash submit_realdata.sh
TASKS=tasks_glioma.csv   bash submit_realdata.sh
bash submit_realdata.sh                                            # genome-wide (tasks.csv)
```
Each task writes `results/<dataset>/tau<tau>/<type>_<rep>.rds` and is skipped if that file exists. Defaults: 4 cores, 32 GB (`MEM=`), 24 h. Without LSF: `Rscript R/03_task.R <id>` with `CQR_RD_TASKS=<table>`.

| Type | What it runs |
|---|---|
| `main` | SS-MS, L-HDCQR with one CV penalty (not reported) and Cox-lasso on the full data |
| `msms` | MS-MS, B = 50; rep 1 is the analysis, rep 2 an independent run with other split seeds (run-to-run Jaccard) |
| `mds_dai` | MDS as in Dai et al. (2023): every split a DS fit with its own 10-fold CV (`CORES=10`) |
| `cmp` | L-HDCQR (Zheng et al., 2018), DS, SCQR-Lasso, SCQR-SCAD from `cqr_comparators/` |
| `decoy`, `decoy_z18`, `decoy_msms` | 20 replicates; min(200, p) row-permuted copies of candidate columns appended (null by construction); `decoy` runs SS-MS and Cox-lasso (and the unreported L-HDCQR variant), `_z18` L-HDCQR, `_msms` MS-MS, all on the same decoys |
| `holdout`, `holdout_z18`, `holdout_msms` | 20 event-stratified 70/30 splits; selection on 70%, unpenalized refit, C-index on 30%; `holdout` runs SS-MS and Cox-lasso (and the unreported L-HDCQR variant), `_z18` L-HDCQR, `_msms` MS-MS, all on the same splits |
| `stab` | SS-MS with 100 split seeds (split-to-split Jaccard) |
| `mds` | earlier MDS implementation, kept so old results stay reproducible; not reported |

### 4. Collect

```bash
CQR_RD_TASKS=tasks_support.csv  CQR_RD_TABLES=tables_support  Rscript R/04_collect.R
CQR_RD_TASKS=tasks_targeted.csv CQR_RD_TABLES=tables_targeted Rscript R/04_collect.R
CQR_RD_TASKS=tasks_glioma.csv   CQR_RD_TABLES=tables_glioma   Rscript R/04_collect.R
Rscript R/04_collect.R; Rscript R/07_signal.R; Rscript R/09_requirements.R   # genome-wide -> tables/
```
Copy the resulting `tables*/` folders into `results_paper/` (`tables/` as `tables_genome/`) and run `R/11_paper_tables.R`.

## Plasmode simulation

The real covariate matrix is kept, s covariates are drawn at random as signals with effects ±b, log T = X_S β + N(0, 1), and censoring is calibrated to the observed rate; 50 replicates per (s, b).

```bash
PM_DATA=SUPPORT2 PM_S=10,30 PM_B=0.05,0.1 PM_METHODS=ssms,msms,hdcqr,lhdcqr_z18 bash plasmode/submit_plasmode.sh
PM_DATA=BRCA_METABRIC                         bash plasmode/submit_plasmode.sh   # default s = 10, 30; b = 0.1, 0.2
PM_DATA=BRCA_METABRIC_TARGETED PM_S=5,15 PM_B=0.15,0.3 bash plasmode/submit_plasmode.sh
PM_DATA=GLIOMA_TCGA            PM_S=5,15 PM_B=0.15,0.3 bash plasmode/submit_plasmode.sh
CQR_PM_DATA=SUPPORT2 Rscript plasmode/plasmode_collect.R      # -> tables/plasmode_{summary,screen}_SUPPORT2.csv
```
4 settings × 50 replicates = 200 tasks per design. Seeds depend on (s, b, replicate) only, so `PM_METHODS=<method>` adds a method to finished tasks on the same simulated data. In the summaries, `SS-MS (path, U+)` is SS-MS as defined in the paper; the `cv` screen and the cutoff `U` are ablations, and `L-HDCQR` (one CV penalty) and `MDS` are not reported.

## Method labels in the result files

| File label | Paper |
|---|---|
| `SS-MS`, `MS-MS`, `DS`, `SCQR-Lasso`, `SCQR-SCAD`, `Cox-lasso` | same |
| `MDS-Dai` | MDS |
| `L-HDCQR-Z18` | L-HDCQR |
| `MDS`, `L-HDCQR` | earlier variants, not reported |

## Pre-specified rules

`R/config.R` fixes, before any selection result was seen, which data sets are analysed and the requirements checked by `09_requirements.R`: A, data integrity; B, n ≥ 200, at least 150 events, τ = 0.25 identifiable; C, τ = 0.25 and 0.5 both identifiable (Kaplan–Meier reaches τ with a margin and enough subjects remain at risk); D, screening cap ⌊d₂/log d₂⌋ ≥ 25 (d₂ = events in the second half), so that U₊ at q = 0.1 can be attained. The gene filter in `02_prepare.R` uses the expression matrix only (MAD), never the outcome.
