# Replication materials

**FDR-Controlled Variable Selection for High-Dimensional Censored Quantile Regression**
Jinzhao Yu and Zili Liu

This archive reproduces every table and figure in Section 4 (simulation studies), Section 5 (real data analysis) and the corresponding parts of the Supplementary Material.

```
package/       cqrfdr_0.3.0.tar.gz: R package implementing SS-MS, MS-MS and the comparators DS, MDS, L-HDCQR
               (also at https://github.com/zili-cuhk/cqrfdr)
simulation/    code for Section 4, and the summaries of the runs reported in the paper
realdata/      code for Section 5, and the collected results reported in the paper
```

## Two levels of reproduction

1. **Tables and figures from the stored results (minutes, any computer).** The summaries of the HPC runs are included, so all tables and figures are rebuilt without rerunning the analyses:
   ```bash
   cd simulation
   Rscript sim/paper_outputs.R       results/main/summary.csv results/paper
   Rscript sim/paper_outputs_extra.R results/extra/summary.csv results/main/summary.csv results/paper
   Rscript sim/paper_table_stability.R results results/paper
   cd ../realdata
   Rscript R/11_paper_tables.R results_paper
   ```
2. **Rerunning the analyses (HPC).** The simulation and the real-data analyses were run on an LSF cluster (MD Anderson "seadragon", R 4.5.1), as job arrays with 4 cores per task. The submission scripts use `bsub`; for another scheduler only the `bsub` header lines need to change. Every task has a deterministic seed, writes one result file and is skipped when its file exists, so interrupted runs can be resubmitted. Approximate cost, from the recorded run times: simulation about 38,000 core-hours, of which about 33,000 are MDS, which cross-validates the penalty in each of its 50 splits (all other methods together need about 5,000); real data about 4,000 core-hours, mostly MS-MS and MDS on the SUPPORT data.

## Software

R 4.5.1. CRAN: `survival`, `quantreg`, `glmnet`, `Rcpp`, `RcppArmadillo`, `data.table`; Bioconductor: `Biobase`, `curatedOvarianData`, `GEOquery` (real-data downloads only). Install the package with
```r
# from GitHub
install.packages("remotes")
remotes::install_github("zili-cuhk/cqrfdr")

# or from the tarball in this archive (identical to the GitHub version 0.3.0)
install.packages("package/cqrfdr_0.3.0.tar.gz", repos = NULL, type = "source")
```
`packageVersion("cqrfdr")` should return 0.3.0. The package source, documentation and vignette are at https://github.com/zili-cuhk/cqrfdr. cqrfdr is written in base R and needs no compiler. SCQR-Lasso and SCQR-SCAD use the authors' C++ code (He, Pan, Tan and Zhou, 2022; github.com/XiaoouPan/scqr, file `hdscqr.cpp`, unchanged except for the C++17 plugin line required by current RcppArmadillo); it is compiled with `Rcpp::sourceCpp` at run time and needs a C++ compiler.

**Third-party code.** The authors' SCQR repository has no license file, so `hdscqr.cpp` is not redistributed in the public GitHub repository of these materials: at the first SCQR call it is downloaded from the authors' repository at a fixed commit (`fetch_hdscqr()` in `R/methods.R`) and the plugin line is changed as above. The same applies to the L-HDCQR author routine `quantproc()` used only by the two check scripts (`simulation/tests/check_author_lhdcqr.R`, `realdata/cqr_comparators/reference/check_author_quantproc.R`). Both files are included in the copy of these materials submitted as Supplementary Material.

Our code is released under GPL (≥ 3); see `LICENSE`.

## Where each result comes from

| Paper | Produced by | Input |
|---|---|---|
| Section 4, Tables 1–3, Figures 1–2 | `simulation/sim/paper_outputs.R` | `simulation/results/main/summary.csv` |
| Section 4, Table 4, Figure 3 | `simulation/sim/paper_outputs_extra.R` | `simulation/results/extra/summary.csv` (+ main for the baseline) |
| Supplement: larger noise scale; first-stage screen | `paper_outputs.R`; `paper_outputs_extra.R` | as above |
| Supplement: split-to-split variability | `simulation/sim/paper_table_stability.R` | `simulation/results/stability_path/`, `simulation/results/stability/` |
| Section 5, Tables 1–3 and all numbers in the text | `realdata/R/11_paper_tables.R` | `realdata/results_paper/` |
| Supplement: plasmode on three designs; genome-wide cohorts | `realdata/R/11_paper_tables.R` | `realdata/results_paper/` |

The pipelines that produce these inputs are described in `simulation/README.md` and `realdata/README.md`.

## Method names

| Paper | Simulation code (`base_method` in summary.csv) | Real-data results |
|---|---|---|
| SS-MS | `ssms_path` | `SS-MS` |
| MS-MS | `msms` | `MS-MS` |
| DS | `ds` | `DS` |
| MDS | `mds_dai` | `MDS-Dai` |
| L-HDCQR | `lhdcqr_z18` | `L-HDCQR-Z18` |
| SCQR-Lasso, SCQR-SCAD | `scqr_lasso`, `scqr_scad` | `SCQR-Lasso`, `SCQR-SCAD` |
| SS-MS with the CV screen (Figure 3, screen table) | `ssms` | — |

Other labels in the result files (`ssms_size`, `mds`, `mds_cv`, `lhdcqr`; real data `MDS`, `L-HDCQR`) are earlier variants that were run during development and are not reported in the paper.

## Data

The simulation needs no data. The real-data analyses use public data that we do not redistribute; `realdata/R/01_download.R` downloads them, except SUPPORT, which is downloaded by hand (see `realdata/README.md`):

- SUPPORT (Knaus et al., 1995): `support2.csv`, Vanderbilt Department of Biostatistics data repository (https://hbiostat.org/data) or UCI Machine Learning Repository.
- METABRIC (Curtis et al., 2012): cBioPortal datahub, study `brca_metabric`.
- TCGA cohorts: UCSC Xena TCGA hub (HiSeqV2 expression) with the TCGA-CDR survival data.
- Ovarian cancer cohorts: Bioconductor package `curatedOvarianData`. DLBCL: GEO GSE10846.

## Code and package

The simulation was run with the code in `simulation/R` and the real-data analyses with the package (SS-MS, MS-MS, MDS) and with `realdata/cqr_comparators/` (L-HDCQR, DS, SCQR), which is the simulation code. The package implements the same algorithms; for example `hdcqr(lambda = "zheng")` returns the same fit as `lhdcqr_z18()` in the simulation code on the same seed and folds.
