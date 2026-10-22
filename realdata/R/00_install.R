## Install the packages needed by the real-data pipeline (run once on the login node).
## Usage: Rscript R/00_install.R ../package/cqrfdr_0.3.0.tar.gz   (cqrfdr from the tarball in this archive)
##        Rscript R/00_install.R github                          (cqrfdr from github.com/zili-cuhk/cqrfdr)
## Packages go to the personal library R_LIBS_USER (created if missing), so no
## interactive "use a personal library?" question can stop the installation.
lib <- Sys.getenv("R_LIBS_USER")
if (!nzchar(lib)) lib <- file.path(Sys.getenv("HOME"), "R", "library")
lib <- path.expand(strsplit(lib, .Platform$path.sep)[[1L]][1L])
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
cat("installing into", lib, "\n")
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 3600)
need <- c("data.table", "glmnet", "survival", "quantreg", "Rcpp", "RcppArmadillo", "BiocManager")
miss <- need[!vapply(need, requireNamespace, TRUE, quietly = TRUE)]
if (length(miss)) install.packages(miss, lib = lib, Ncpus = 4)
bioc <- c("Biobase", "curatedOvarianData", "GEOquery")
miss <- bioc[!vapply(bioc, requireNamespace, TRUE, quietly = TRUE)]
if (length(miss)) BiocManager::install(miss, lib = lib, update = FALSE, ask = FALSE, Ncpus = 4)
tgz <- commandArgs(TRUE)[1L]
if (identical(tgz, "github")) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", lib = lib)
  remotes::install_github("zili-cuhk/cqrfdr", lib = lib, upgrade = "never")
} else if (!is.na(tgz)) install.packages(tgz, lib = lib, repos = NULL, type = "source")
cat("\n")
for (p in c(need, bioc, "cqrfdr"))
  cat(sprintf("%-20s %s\n", p, if (requireNamespace(p, quietly = TRUE)) "OK" else "MISSING"))
