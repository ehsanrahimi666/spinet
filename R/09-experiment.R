# ---------------------------------------------------------------------------
# Module 8b: virtual experiments
# ---------------------------------------------------------------------------

#' Run a factorial virtual experiment
#'
#' Builds a design from a list of factor levels, runs one complete
#' current-to-future analysis per design row, and returns a tidy table of
#' responses. Each run is independently seeded, so the whole experiment is
#' exactly reproducible and can be split across machines.
#'
#' @param design A named list of factor levels, or a `data.frame` in which each
#'   row is one run. Recognised factors are
#'   \describe{
#'     \item{`richness`}{number of species per level.}
#'     \item{`connectance`}{connectance of the constraint layer.}
#'     \item{`architecture`}{`"random"`, `"nested"`, `"modular"`,
#'       `"heterogeneous"`.}
#'     \item{`scenario`}{climate scenario, see [si_simulate_climate()].}
#'     \item{`severity`}{magnitude of climate change.}
#'     \item{`shift_sd`}{among-species standard deviation of the phenological
#'       shift. Zero means every species shifts alike, which permits no
#'       rewiring.}
#'     \item{`dynamic_phenology`}{logical; recompute the constraint layer from
#'       the future phenology.}
#'     \item{`threshold`}{suitability threshold for the binary pathway.}
#'     \item{`pathway`}{`"binary"` or `"continuous"`.}
#'     \item{`grid`}{side length of the simulation landscape.}
#'   }
#'   Unspecified factors take their defaults.
#' @param replicates Number of independent replicates per design row.
#' @param metrics Character vector of network metrics to record.
#' @param robustness Logical. Also compute \eqn{R} and \eqn{R_w} on the
#'   pooled network. Adds appreciable cost.
#' @param cells Integer. Number of grid cells sampled per run for the per-cell
#'   metrics; `NULL` uses all of them. Sampling is what keeps large designs
#'   tractable.
#' @param cores Number of parallel workers.
#' @param seed Base random seed.
#' @param progress Logical. Print a progress bar.
#'
#' @return A `data.frame` of class `si_experiment`: one row per run, carrying
#'   the design columns and the responses.
#'
#' @examples
#' # a small design: does constraint density govern the network response?
#' des <- list(richness = 20, connectance = c(1, 0.5, 0.15),
#'             scenario = "decline", severity = 0.4, grid = 14)
#' ex <- si_experiment(des, replicates = 2, cells = 40, seed = 1)
#' ex[, c("connectance", "d_connectance", "d_NODF", "beta_WN", "beta_OS")]
#'
#' # rewiring appears only when the phenological shift varies among species
#' des2 <- list(richness = 20, connectance = 0.3, shift_sd = c(0, 30),
#'              dynamic_phenology = TRUE, severity = 0.3, grid = 14)
#' ex2 <- si_experiment(des2, replicates = 2, cells = 40, seed = 2)
#' ex2[, c("shift_sd", "beta_WN", "beta_OS", "rewiring_gains")]
#' @seealso [si_simulate_climate()], [si_beta()], [si_sensitivity()]
#' @export
si_experiment <- function(design, replicates = 1,
                          metrics = c("links", "connectance", "NODF",
                                      "H2prime", "evenness"),
                          robustness = FALSE, cells = 100, cores = 1,
                          seed = 1, progress = FALSE) {

  defaults <- list(richness = 25, connectance = 0.25, architecture = "random",
                   scenario = "decline", severity = 0.3, shift_sd = 0,
                   dynamic_phenology = FALSE, threshold = 0.3,
                   pathway = "binary", grid = 20, distance = 4)
  if (is.list(design) && !is.data.frame(design)) {
    for (nm in names(defaults)) if (is.null(design[[nm]])) design[[nm]] <- defaults[[nm]]
    grid <- expand.grid(design, stringsAsFactors = FALSE,
                        KEEP.OUT.ATTRS = FALSE)
  } else {
    grid <- as.data.frame(design, stringsAsFactors = FALSE)
    for (nm in names(defaults)) if (is.null(grid[[nm]])) grid[[nm]] <- defaults[[nm]]
  }
  grid <- grid[rep(seq_len(nrow(grid)), each = replicates), , drop = FALSE]
  grid$replicate <- rep(seq_len(replicates), times = nrow(grid) / replicates)
  grid$run <- seq_len(nrow(grid))
  rownames(grid) <- NULL

  runner <- function(i) .si_run_one(grid[i, ], metrics, robustness, cells,
                                    seed = seed + i)
  res <- if (cores > 1 && .Platform$OS.type == "unix") {
    parallel::mclapply(seq_len(nrow(grid)), runner, mc.cores = cores)
  } else if (progress) {
    pb <- utils::txtProgressBar(0, nrow(grid), style = 3)
    o <- vector("list", nrow(grid))
    for (i in seq_len(nrow(grid))) { o[[i]] <- runner(i); utils::setTxtProgressBar(pb, i) }
    close(pb); o
  } else lapply(seq_len(nrow(grid)), runner)

  out <- cbind(grid, do.call(rbind, lapply(res, function(z) as.data.frame(t(z)))))
  structure(out, class = c("si_experiment", "data.frame"))
}

.si_run_one <- function(p, metrics, robustness, cells, seed) {
  set.seed(seed)
  S <- as.integer(p$richness); g <- as.integer(p$grid)

  A1 <- si_simulate_species(S, g, g, prefix = "A", seed = seed)
  B1 <- si_simulate_species(S, g, g, prefix = "B", seed = seed + 1e6)
  A2 <- si_simulate_climate(A1, p$scenario, severity = p$severity,
                            distance = p$distance, seed = seed)
  B2 <- si_simulate_climate(B1, p$scenario, severity = p$severity,
                            distance = p$distance, seed = seed + 1)

  ph <- si_phenology_simulate(S, S, shift_mean = -10,
                              shift_sd = p$shift_sd, seed = seed + 2)
  base <- si_simulate_metaweb(S, S, connectance = p$connectance,
                              architecture = p$architecture, seed = seed + 3,
                              names_A = A1$names, names_B = B1$names)
  if (isTRUE(p$dynamic_phenology)) {
    # constraint layer = structural metaweb AND phenological overlap
    thr <- stats::quantile(si_phenology_overlap(ph, "current"),
                           1 - p$connectance)
    M1 <- base$matrix * (si_phenology_overlap(ph, "current") > thr)
    M2 <- base$matrix * (si_phenology_overlap(ph, "future")  > thr)
  } else {
    M1 <- M2 <- base$matrix
  }
  mw1 <- si_metaweb(M1, names_A = A1$names, names_B = B1$names)
  mw2 <- si_metaweb(M2, names_A = A1$names, names_B = B1$names)

  if (identical(p$pathway, "binary")) {
    A1 <- si_threshold(A1, "fixed", p$threshold); B1 <- si_threshold(B1, "fixed", p$threshold)
    A2 <- si_threshold(A2, "fixed", p$threshold); B2 <- si_threshold(B2, "fixed", p$threshold)
    meth <- "binary"; fl <- 0
  } else { meth <- "product"; fl <- 0.05 }

  f1 <- si_overlap(A1, B1, mw1, method = meth, floor = fl,
                   scenario = "current", check = FALSE)
  f2 <- si_overlap(A2, B2, mw2, method = meth, floor = fl,
                   scenario = "future", check = FALSE)

  k <- if (is.null(cells) || cells >= f1$n_cells) seq_len(f1$n_cells) else
    sample.int(f1$n_cells, cells)
  m1 <- si_metrics(f1, what = metrics, cells = k)
  m2 <- si_metrics(f2, what = metrics, cells = k)
  bb <- suppressWarnings(si_beta(f1, f2, cells = k))
  rw <- si_rewiring(f1, f2, cells = k)
  ex <- si_network_extinction(f1, f2)

  cur <- vapply(metrics, function(v) mean(m1[[v]], na.rm = TRUE), numeric(1))
  fut <- vapply(metrics, function(v) mean(m2[[v]], na.rm = TRUE), numeric(1))
  out <- c(stats::setNames(cur, paste0("cur_", metrics)),
           stats::setNames(fut, paste0("fut_", metrics)),
           stats::setNames(fut - cur, paste0("d_", metrics)),
           beta_S  = mean(bb$beta_S,  na.rm = TRUE),
           beta_WN = mean(bb$beta_WN, na.rm = TRUE),
           beta_ST = mean(bb$beta_ST, na.rm = TRUE),
           beta_OS = mean(bb$beta_OS, na.rm = TRUE),
           rewiring_gains  = sum(rw$realised, na.rm = TRUE),
           rewiring_losses = sum(rw$lost, na.rm = TRUE),
           n_extinct = sum(ex$status == "extinct"),
           n_persist = sum(ex$status == "persists"),
           realised_C = mean(m1$connectance, na.rm = TRUE))
  if (robustness) {
    pooled1 <- .pool_web(f1, k)
    r0 <- si_robustness(pooled1, n_sim = 25)
    r1 <- si_robustness(pooled1, rewiring = mw1$matrix, n_sim = 25)
    out <- c(out, R = r0$robustness, Rw = r1$robustness,
             D = r1$robustness - r0$robustness)
  }
  out
}

.pool_web <- function(field, cells) {
  W <- matrix(0, length(field$names_A), length(field$names_B))
  for (k in cells) W <- W + .local_matrix(field, k, drop = FALSE)
  dimnames(W) <- list(field$names_A, field$names_B)
  rk <- rowSums(W) > 0; ck <- colSums(W) > 0
  W[rk, ck, drop = FALSE]
}

#' @export
print.si_experiment <- function(x, ...) {
  cat("<si_experiment>  runs: ", nrow(x), "\n", sep = "")
  fac <- intersect(c("richness", "connectance", "architecture", "scenario",
                     "severity", "shift_sd", "pathway"), names(x))
  for (f in fac) {
    u <- unique(x[[f]])
    if (length(u) > 1)
      cat("  ", f, ": ", paste(utils::head(u, 6), collapse = ", "),
          if (length(u) > 6) " ..." else "", "\n", sep = "")
  }
  invisible(x)
}

#' Variance decomposition of an experiment
#'
#' Fits a linear model of a response on the design factors and reports the
#' proportion of variance each factor explains. This is how a factorial
#' experiment answers "which driver matters most" rather than only "does this
#' driver matter".
#'
#' @param experiment An [si_experiment()] result.
#' @param response Character name of the response column.
#' @param factors Character vector of design columns; defaults to those that
#'   actually vary.
#' @param interactions Logical. Include all two-way interactions.
#' @return A `data.frame` of variance components, sorted by importance.
#' @examples
#' des <- list(richness = c(15, 30), connectance = c(0.15, 0.5),
#'             severity = c(0.2, 0.5), grid = 12)
#' ex <- si_experiment(des, replicates = 2, cells = 30, seed = 3)
#' si_sensitivity(ex, "d_connectance")
#' @export
si_sensitivity <- function(experiment, response, factors = NULL,
                           interactions = TRUE) {
  d <- as.data.frame(experiment)
  if (!response %in% names(d))
    stop("`response` '", response, "' not found. Available: ",
         paste(utils::head(setdiff(names(d), c("run", "replicate")), 12),
               collapse = ", "), call. = FALSE)
  cand <- c("richness", "connectance", "architecture", "scenario", "severity",
            "shift_sd", "dynamic_phenology", "threshold", "pathway", "grid")
  if (is.null(factors))
    factors <- cand[cand %in% names(d) &
                      vapply(cand, function(f) f %in% names(d) &&
                               length(unique(d[[f]])) > 1, logical(1))]
  if (!length(factors))
    stop("No design factor varies; nothing to decompose.", call. = FALSE)
  d <- d[stats::complete.cases(d[, c(response, factors)]), , drop = FALSE]
  rhs <- if (interactions && length(factors) > 1)
    paste0("(", paste(factors, collapse = " + "), ")^2") else
      paste(factors, collapse = " + ")
  fit <- stats::lm(stats::as.formula(paste(response, "~", rhs)), data = d)
  a <- stats::anova(fit)
  ss <- a[["Sum Sq"]]
  out <- data.frame(term = rownames(a), variance_explained = ss / sum(ss),
                    row.names = NULL)
  out <- out[order(-out$variance_explained), ]
  rownames(out) <- NULL
  out
}

#' @export
`[.si_experiment` <- function(x, ...) {
  y <- NextMethod("[")
  if (is.data.frame(y)) class(y) <- "data.frame"
  y
}
