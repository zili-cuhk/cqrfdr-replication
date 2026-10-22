## ---------------------------------------------------------------------------
## Append the MS-MS and Dai-design MDS tasks to the existing task tables
## ---------------------------------------------------------------------------
## Usage: Rscript R/10_add_tasks.R [types] [datasets]   (default: msms,mds_dai; all data sets)
##   e.g. Rscript R/10_add_tasks.R cmp,decoy_z18,holdout_z18 SUPPORT2,BRCA_METABRIC_TARGETED
##        Rscript R/10_add_tasks.R decoy_msms,holdout_msms SUPPORT2,BRCA_METABRIC_TARGETED
## For every task table present (tasks.csv, tasks_targeted.csv, tasks_glioma.csv,
## tasks_support.csv) and every data set x tau with a "main" task, adds
##   msms    rep 1 and rep 2   (rep 2: independent run, run-to-run agreement)
##   mds_dai rep 1
## New rows get ids after the existing ones, so finished results and the ids of
## old tasks are unchanged.  Rows already present are not added twice.  Prints
## the id range of each new type for submit_realdata.sh (ONLY=...).
## Run it again after any of 02_prepare.R / 05 / 06 / 08 rewrites a task table.

source(file.path(Sys.getenv("CQR_RD_ROOT", "."), "R", "config.R"))
a <- commandArgs(TRUE)
types <- if (length(a)) strsplit(a[1], ",")[[1]] else c("msms", "mds_dai")
only_ds <- if (length(a) > 1L) strsplit(a[2], ",")[[1]] else NULL
reps <- c(msms = 2L, mds_dai = 1L, cmp = 1L, decoy_z18 = CFG$n_decoy, holdout_z18 = CFG$n_holdout,
          decoy_msms = CFG$n_decoy, holdout_msms = CFG$n_holdout)
if (!all(types %in% names(reps))) stop("unknown task type: ", paste(setdiff(types, names(reps)), collapse = ", "))
files <- c("tasks.csv", "tasks_targeted.csv", "tasks_glioma.csv", "tasks_support.csv")
for (f in files) {
  path <- file.path(RD_ROOT, f)
  if (!file.exists(path)) { cat(sprintf("%-20s not found (skipped)\n", f)); next }
  tk <- utils::read.csv(path, stringsAsFactors = FALSE)
  dt <- unique(tk[tk$type == "main", c("dataset", "tau")])
  if (!is.null(only_ds)) dt <- dt[dt$dataset %in% only_ds, , drop = FALSE]
  if (!nrow(dt)) { cat(sprintf("%-20s no matching data set\n", f)); next }
  new <- do.call(rbind, lapply(types, function(ty) do.call(rbind, lapply(seq_len(reps[[ty]]), function(r)
    data.frame(dataset = dt$dataset, tau = dt$tau, type = ty, rep = r)))))
  key <- function(z) paste(z$dataset, z$tau, z$type, z$rep)
  new <- new[!key(new) %in% key(tk), , drop = FALSE]
  if (!nrow(new)) { cat(sprintf("%-20s nothing to add (%d tasks)\n", f, nrow(tk))); next }
  new <- new[order(match(new$type, types), new$rep, match(new$dataset, unique(tk$dataset)), new$tau), ]
  new$id <- max(tk$id) + seq_len(nrow(new))
  for (nm in setdiff(names(tk), names(new))) new[[nm]] <- NA
  tk <- rbind(tk, new[, names(tk)])
  utils::write.csv(tk, path, row.names = FALSE)
  cat(sprintf("%-20s added %d tasks (now %d)\n", f, nrow(new), nrow(tk)))
  for (ty in types) {
    if (reps[[ty]] <= 2L) for (r in seq_len(reps[[ty]])) {
      ids <- new$id[new$type == ty & new$rep == r]
      if (length(ids)) cat(sprintf("    %-11s rep %d : ids %d-%d (%d tasks)\n", ty, r, min(ids), max(ids), length(ids)))
    } else {
      ids <- new$id[new$type == ty]
      if (length(ids)) cat(sprintf("    %-11s       : ids %d-%d (%d tasks)\n", ty, min(ids), max(ids), length(ids)))
    }
  }
}
