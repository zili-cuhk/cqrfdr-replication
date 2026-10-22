## ---------------------------------------------------------------------------
## Seeds, metrics, I/O helpers
## ---------------------------------------------------------------------------

## Deterministic 31-bit seed from (master, scenario name, replication, role).
## Every method gets its own stream, so adding/removing a method or changing
## the execution order never changes the data or another method's split.
str_hash <- function(s) {
  h <- 0
  for (v in utf8ToInt(s)) h <- (h * 131 + v) %% 2147483629
  h
}
make_seed <- function(master, scenario, rep, role) {
  as.integer((str_hash(paste(master, scenario, rep, role, sep = "|")) %% 2147483000) + 1)
}
with_seed <- function(seed, expr) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  expr
}

selection_metrics <- function(selected, active, p) {
  V <- length(setdiff(selected, active)); R <- length(selected)
  T <- length(intersect(selected, active))
  c(R = R, V = V, T = T, FDP = V / max(R, 1), TPP = T / length(active))
}

## atomic save: write to a temporary file in the same directory, then rename
save_atomic <- function(obj, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dirname(path), fileext = ".tmp")
  saveRDS(obj, tmp)
  if (!file.rename(tmp, path)) stop("could not move ", tmp, " to ", path)
  invisible(path)
}

grid_taus <- function(tau, cfg) seq(cfg$tau_lower, tau, length.out = cfg$n_grid)

load_sources <- function(root) {
  for (f in c("dgp.R", "solver.R", "methods.R", "utils.R"))
    sys.source(file.path(root, "R", f), envir = globalenv())
}
