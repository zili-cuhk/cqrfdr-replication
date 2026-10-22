## Summarise the plasmode simulation: FDR (mean FDP) and power with Monte Carlo SE.
## Usage: CQR_PM_DATA=<design> Rscript plasmode/plasmode_collect.R   (default BRCA_METABRIC)
## Writes tables/plasmode_summary_<design>.csv and tables/plasmode_screen_<design>.csv.
source(file.path(Sys.getenv("CQR_RD_ROOT", "."), "R", "config.R"))
PM_DATA <- Sys.getenv("CQR_PM_DATA", "BRCA_METABRIC")
tag <- if (PM_DATA == "BRCA_METABRIC") "" else paste0("_", PM_DATA)
fs <- list.files(file.path(RD_ROOT, "results", paste0("plasmode", tag)), "^task_.*\\.rds$", full.names = TRUE)
if (!length(fs)) stop("no plasmode results yet")
R <- do.call(rbind, lapply(fs, function(f) readRDS(f)$rows))
R$q <- ifelse(is.na(R$q), "-", format(R$q))      # L-HDCQR has no q; keep it in the grouping
sc <- do.call(rbind, lapply(fs, function(f) { r <- readRDS(f)
  data.frame(s = r$s, b = r$b, path_size = r$screen_path$size, path_recall = r$screen_path$recall,
             cv_size = r$screen_cv$size, cv_recall = r$screen_cv$recall) }))
key <- interaction(R$s, R$b, R$method, R$q, drop = TRUE, lex.order = TRUE)
S <- do.call(rbind, lapply(split(R, key), function(g) data.frame(
  s = g$s[1], b = g$b[1], method = g$method[1], q = g$q[1], reps = nrow(g),
  FDR = round(mean(g$FDP), 3), FDR_se = round(stats::sd(g$FDP) / sqrt(nrow(g)), 3),
  power = round(mean(g$power), 3), power_se = round(stats::sd(g$power) / sqrt(nrow(g)), 3),
  mean_selected = round(mean(g$n_sel), 1))))
S <- S[order(S$s, S$b, S$q, S$method), ]
dir.create(file.path(RD_ROOT, "tables"), showWarnings = FALSE)
write.csv(S, file.path(RD_ROOT, "tables", paste0("plasmode_summary_", PM_DATA, ".csv")), row.names = FALSE)
cat(sprintf("%d tasks finished\n\n", length(fs)))
print(S, row.names = FALSE)
cat("\nScreen size and recall of the true signals:\n")
SC <- stats::aggregate(cbind(path_size, path_recall, cv_size, cv_recall) ~ s + b, sc, function(v) round(mean(v), 2))
print(SC)
write.csv(SC, file.path(RD_ROOT, "tables", paste0("plasmode_screen_", PM_DATA, ".csv")), row.names = FALSE)
