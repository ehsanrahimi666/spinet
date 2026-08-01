# ---------------------------------------------------------------------------
# Module 4: network metrics, computed natively, with degeneracy detection
# ---------------------------------------------------------------------------

#' Detect structurally degenerate interaction fields
#'
#' Reports which network metrics are mathematically constant for a given field,
#' and therefore carry no information about space or about environmental
#' change. This diagnostic is the reason `spinet` exists, and it runs by
#' default inside [si_overlap()].
#'
#' @param field An [si_field()].
#' @param warn Logical. Emit a warning when a degeneracy is found.
#'
#' @return A `data.frame` with one row per issue, invisibly, with columns
#'   `issue`, `affects` and `explanation`. A zero-row result means no
#'   degeneracy was detected.
#'
#' @details
#' Three distinct failure modes are detected.
#'
#' **Rank-one degeneracy.** With no constraint layer and a multiplicative
#' combination rule, \eqn{W_k = p_k \otimes a_k} has rank one. Then
#' connectance is identically 1, NODF identically 0, and
#' \eqn{H_2' \equiv 0} exactly, because the maximum-entropy matrix with the
#' same marginals *is* \eqn{p_k \otimes a_k}. All 'topological' variation
#' across such a field is variation in the marginal suitabilities, not in
#' network structure.
#'
#' **Complete-bipartite degeneracy.** With no constraint layer and
#' `method = "binary"`, every local network is the complete bipartite graph
#' \eqn{K_{n,m}}. Every binary metric is a deterministic function of
#' \eqn{(n, m)} alone.
#'
#' **Floorless-continuous degeneracy.** Correlative suitability surfaces are
#' strictly positive almost everywhere, so with `floor = 0` the support of
#' \eqn{W_k} equals the metaweb in every cell. Presence-based metrics then
#' return the metaweb itself, invariant in space *and* time. Fix with
#' [si_prune()].
#'
#' A fourth, non-local degeneracy is reported by [si_beta()]: with a
#' time-invariant constraint layer, interaction rewiring \eqn{\beta_{OS}} is
#' identically zero.
#'
#' @examples
#' # unconstrained continuous field: rank one
#' P <- matrix(runif(100 * 8), 100, 8); A <- matrix(runif(100 * 6), 100, 6)
#' f_bad <- si_overlap(P, A, metaweb = NULL, method = "product", check = FALSE)
#' si_degenerate(f_bad, warn = FALSE)
#'
#' # empirically: connectance is 1 and NODF is 0 in every cell
#' m <- si_local(f_bad, 1, drop = FALSE)
#' c(connectance = mean(m > 0), NODF = si_nodf(m > 0))
#'
#' # constrained field: no degeneracy
#' nrow(si_degenerate(si_example_field(), warn = FALSE))
#' @export
si_degenerate <- function(field, warn = TRUE) {
  stopifnot(inherits(field, "si_field"))
  unc  <- isTRUE(attr(field$F, "unconstrained"))
  rank1 <- field$method %in% c("product", "geometric", "harmonic")
  cont  <- !field$A$binary || !field$B$binary
  out <- data.frame(issue = character(), affects = character(),
                    explanation = character(), stringsAsFactors = FALSE)
  add <- function(i, a, e) out[nrow(out) + 1L, ] <<- list(i, a, e)

  if (unc && rank1)
    add("rank-one",
        "connectance, NODF, WNODF, H2prime, modularity",
        paste("No constraint layer with a multiplicative rule gives",
              "W = p (x) a, rank 1. Connectance == 1, NODF == 0 and",
              "H2' == 0 exactly. Supply a metaweb."))
  if (unc && field$method == "binary")
    add("complete-bipartite",
        "all binary metrics",
        paste("Every local network is the complete bipartite graph K(n,m);",
              "every binary metric is a function of richness alone.",
              "Supply a metaweb."))
  if (!unc && cont && field$floor <= 0 && field$method != "binary")
    add("floorless-continuous",
        "connectance, links, NODF, degree, modularity",
        paste("Continuous suitabilities are positive almost everywhere, so",
              "the support of W equals the metaweb in every cell. Binary",
              "metrics will be constant in space and time.",
              "Use si_prune() or set floor > 0."))
  if (warn && nrow(out))
    warning("si_field is structurally degenerate:\n",
            paste0("  * ", out$issue, " -> ", out$affects, collapse = "\n"),
            "\n  See ?si_degenerate.", call. = FALSE)
  invisible(out)
}


# --- single-matrix metric primitives ---------------------------------------

#' Nestedness (NODF) of a binary matrix
#'
#' Nestedness based on Overlap and Decreasing Fill (Almeida-Neto and others,
#' 2008), implemented with matrix products so that it is fast enough to run on
#' hundreds of thousands of local networks.
#'
#' @param m A binary (or coercible) matrix.
#' @return Numeric NODF on 0--100, or `NA` for degenerate matrices.
#' @examples
#' set.seed(1)
#' perfect <- outer(1:6, 1:6, function(i, j) as.numeric(j <= 7 - i))
#' si_nodf(perfect)          # a perfectly nested matrix
#' si_nodf(matrix(1, 6, 6))  # complete graph: no nestedness by definition
#' @references
#' Almeida-Neto, M., Guimaraes, P., Guimaraes, P.R., Loyola, R.D. & Ulrich, W.
#' (2008) A consistent metric for nestedness analysis in ecological systems.
#' *Oikos*, 117, 1227--1239.
#' @export
si_nodf <- function(m) {
  b <- (as.matrix(m) > 0) * 1
  if (nrow(b) < 2 || ncol(b) < 2 || sum(b) == 0) return(NA_real_)
  half <- function(x) {
    d <- rowSums(x)
    if (length(d) < 2) return(c(0, 0))
    O <- tcrossprod(x)
    dm <- outer(d, d, pmin)
    neq <- outer(d, d, "!=") & dm > 0
    ct <- ifelse(neq, 100 * O / pmax(dm, 1e-12), 0)
    iu <- upper.tri(ct)
    c(sum(ct[iu]), sum(iu))
  }
  r <- half(b); cc <- half(t(b))
  tot <- r[2] + cc[2]
  if (tot == 0) return(NA_real_)
  (r[1] + cc[1]) / tot
}

#' Weighted nestedness (WNODF)
#'
#' The quantitative extension of NODF (Almeida-Neto & Ulrich, 2011): paired
#' overlap counts only those cells in which the less-generalist species has a
#' strictly smaller weight than the more-generalist one.
#'
#' @param m A weighted matrix.
#' @return Numeric WNODF on 0--100.
#' @examples
#' w <- matrix(c(5, 3, 1, 0,  4, 2, 0, 0,  3, 0, 0, 0,  1, 0, 0, 0), 4, 4,
#'             byrow = TRUE)
#' si_wnodf(w)
#' @references
#' Almeida-Neto, M. & Ulrich, W. (2011) A straightforward computational
#' approach for measuring nestedness using quantitative matrices.
#' *Environmental Modelling & Software*, 26, 173--178.
#' @export
si_wnodf <- function(m) {
  w <- as.matrix(m)
  if (nrow(w) < 2 || ncol(w) < 2 || sum(w > 0) == 0) return(NA_real_)
  half <- function(x) {
    d <- rowSums(x > 0); n <- nrow(x); tot <- 0; cnt <- 0
    if (n < 2) return(c(0, 0))
    for (i in 1:(n - 1)) for (j in (i + 1):n) {
      cnt <- cnt + 1
      if (d[i] == d[j] || d[i] == 0 || d[j] == 0) next
      hi <- if (d[i] > d[j]) i else j
      lo <- if (d[i] > d[j]) j else i
      sel <- x[lo, ] > 0
      if (!any(sel)) next
      tot <- tot + 100 * sum(x[lo, sel] < x[hi, sel]) / sum(sel)
    }
    c(tot, cnt)
  }
  r <- half(w); cc <- half(t(w))
  tot <- r[2] + cc[2]
  if (tot == 0) return(NA_real_)
  (r[1] + cc[1]) / tot
}

#' Network-level specialisation H2'
#'
#' The standardised two-dimensional Shannon entropy of Bluthgen and others
#' (2006). `H2max` is the entropy of the maximum-entropy matrix with the same
#' marginal totals, which is exactly the outer product of the marginals;
#' `H2min` is obtained by the standard greedy concentration algorithm.
#'
#' @param m A weighted matrix.
#' @param components Logical. Return `H2`, `H2max`, `H2min` and `H2prime`
#'   rather than `H2prime` alone. The raw difference `H2max - H2` is the useful
#'   quantity when comparing fields, because it is zero exactly when the matrix
#'   is a pure outer product.
#' @return A number, or a named numeric vector when `components = TRUE`.
#' @examples
#' # a pure outer product has H2' exactly 0 -- this is the degeneracy
#' p <- runif(20); a <- runif(15)
#' si_h2prime(outer(p, a), components = TRUE)
#'
#' # a constrained matrix does not
#' set.seed(2)
#' mw <- si_simulate_metaweb(20, 15, connectance = 0.2, seed = 2)
#' si_h2prime(outer(p, a) * mw$matrix, components = TRUE)
#' @references
#' Bluthgen, N., Menzel, F. & Bluthgen, N. (2006) Measuring specialization in
#' species interaction networks. *BMC Ecology*, 6, 9.
#' @export
si_h2prime <- function(m, components = FALSE) {
  w <- as.matrix(m); w[is.na(w)] <- 0
  tot <- sum(w)
  if (tot <= 0 || nrow(w) < 2 || ncol(w) < 2) {
    v <- c(H2 = NA, H2max = NA, H2min = NA, H2prime = NA)
    return(if (components) v else unname(v["H2prime"]))
  }
  ent <- function(x) { x <- x[x > 0]; -sum(x * log(x)) }
  H2  <- ent(w / tot)
  rs <- rowSums(w); cs <- colSums(w)
  H2max <- ent(outer(rs, cs) / tot^2)
  H2min <- ent(.greedy_concentrate(rs, cs) / tot)
  den <- H2max - H2min
  H2p <- if (abs(den) < .Machine$double.eps^0.5) 0 else (H2max - H2) / den
  v <- c(H2 = H2, H2max = H2max, H2min = H2min,
         H2prime = max(0, min(1, H2p)))
  if (components) v else unname(v["H2prime"])
}

.greedy_concentrate <- function(rs, cs) {
  # maximally concentrated matrix with the given marginals (minimum entropy)
  M <- matrix(0, length(rs), length(cs))
  r <- rs; cc <- cs
  repeat {
    i <- which.max(r); j <- which.max(cc)
    v <- min(r[i], cc[j])
    if (!is.finite(v) || v <= 0) break
    M[i, j] <- M[i, j] + v
    r[i] <- r[i] - v; cc[j] <- cc[j] - v
    if (sum(r) <= .Machine$double.eps) break
  }
  M
}

#' Shannon interaction evenness
#' @param m A weighted matrix.
#' @return Evenness on 0--1.
#' @examples
#' si_evenness(matrix(c(4, 1, 1, 4), 2, 2))
#' si_evenness(matrix(1, 4, 4))   # perfectly even
#' @export
si_evenness <- function(m) {
  w <- as.matrix(m); p <- w[w > 0]
  if (length(p) < 2) return(NA_real_)
  p <- p / sum(p)
  -sum(p * log(p)) / log(length(p))
}

#' Modularity of a bipartite network
#'
#' Barber's bipartite modularity for the partition returned by a fast greedy
#' community search. Requires the optional `igraph` package; returns `NA` if it
#' is not installed.
#'
#' @param m A weighted matrix.
#' @return Numeric modularity, or `NA`.
#' @examples
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   set.seed(4)
#'   modular <- si_simulate_metaweb(20, 20, connectance = 0.15,
#'                                  architecture = "modular", n_modules = 4,
#'                                  seed = 4)
#'   random  <- si_simulate_metaweb(20, 20, connectance = 0.15,
#'                                  architecture = "random", seed = 4)
#'   c(modular = si_modularity(modular$matrix),
#'     random  = si_modularity(random$matrix))
#' }
#' @export
si_modularity <- function(m) {
  if (!requireNamespace("igraph", quietly = TRUE)) return(NA_real_)
  w <- as.matrix(m)
  if (nrow(w) < 2 || ncol(w) < 2 || sum(w > 0) < 2) return(NA_real_)
  nr <- nrow(w); nc <- ncol(w)
  full <- matrix(0, nr + nc, nr + nc)
  full[1:nr, (nr + 1):(nr + nc)] <- w
  full[(nr + 1):(nr + nc), 1:nr] <- t(w)
  g <- igraph::graph_from_adjacency_matrix(full, mode = "undirected",
                                           weighted = TRUE, diag = FALSE)
  if (igraph::gorder(g) < 2 || igraph::gsize(g) < 1) return(NA_real_)
  cl <- try(igraph::cluster_fast_greedy(g), silent = TRUE)
  if (inherits(cl, "try-error")) return(NA_real_)
  as.numeric(igraph::modularity(g, igraph::membership(cl),
                                weights = igraph::E(g)$weight))
}

#' All network metrics for one matrix
#'
#' @param m A weighted or binary matrix.
#' @param what Character vector of metric names, or `"all"`.
#' @return A named numeric vector.
#' @examples
#' si_web_metrics(si_local(si_example_field(), 3))
#' @export
si_web_metrics <- function(m, what = "all") {
  w <- as.matrix(m); w[is.na(w)] <- 0
  b <- w > 0
  nr <- nrow(w); nc <- ncol(w); L <- sum(b)
  if (nr == 0 || nc == 0 || L == 0) {
    v <- c(n_A = nr, n_B = nc, links = 0, connectance = NA, weighted_links = 0,
           linkage_density = NA, generality = NA, vulnerability = NA,
           NODF = NA, WNODF = NA, H2prime = NA, H2_gap = NA,
           evenness = NA, modularity = NA, web_asymmetry = NA)
    if (!identical(what, "all")) v <- v[intersect(what, names(v))]
    return(v)
  }
  h2 <- si_h2prime(w, components = TRUE)
  gen <- mean(rowSums(b)[rowSums(b) > 0])
  vul <- mean(colSums(b)[colSums(b) > 0])
  v <- c(n_A = nr, n_B = nc, links = L, connectance = L / (nr * nc),
         weighted_links = sum(w),
         linkage_density = L / (nr + nc),
         generality = gen, vulnerability = vul,
         NODF = si_nodf(b), WNODF = si_wnodf(w),
         H2prime = unname(h2["H2prime"]),
         H2_gap = unname(h2["H2max"] - h2["H2"]),
         evenness = si_evenness(w),
         modularity = si_modularity(w),
         web_asymmetry = (nc - nr) / (nc + nr))
  if (!identical(what, "all")) v <- v[intersect(what, names(v))]
  v
}


#' Per-cell network metrics for a field
#'
#' @param field An [si_field()].
#' @param what Character vector of metric names (see [si_web_metrics()]) or
#'   `"all"`. Restricting `what` is worthwhile: `NODF`, `WNODF` and
#'   `modularity` dominate the cost.
#' @param cells Optional integer vector of cells to evaluate; defaults to all.
#' @param min_species Integer. Cells with fewer than this many species in
#'   either level return `NA` rather than a degenerate value.
#' @param cores Integer. Number of parallel workers (forked; 1 on Windows).
#' @param progress Logical. Print a progress bar.
#' @return An object of class `si_metrics`: a `data.frame` with one row per
#'   cell, carrying the field's template as an attribute so that
#'   [si_metric_map()] can rasterise it.
#' @examples
#' f <- si_example_field(n_cells = 60)
#' mt <- si_metrics(f, what = c("links", "connectance", "NODF", "H2prime"))
#' head(mt)
#' summary(mt)
#' @seealso [si_metric_map()], [si_degenerate()]
#' @export
si_metrics <- function(field, what = "all", cells = NULL, min_species = 2,
                       cores = 1, progress = FALSE) {
  stopifnot(inherits(field, "si_field"))
  si_degenerate(field, warn = TRUE)
  k <- if (is.null(cells)) seq_len(field$n_cells) else as.integer(cells)
  one <- function(i) {
    w <- .local_matrix(field, i, drop = TRUE)
    if (nrow(w) < min_species || ncol(w) < min_species)
      return(si_web_metrics(matrix(0, 0, 0), what))
    si_web_metrics(w, what)
  }
  res <- if (cores > 1 && .Platform$OS.type == "unix") {
    parallel::mclapply(k, one, mc.cores = cores)
  } else if (progress) {
    pb <- utils::txtProgressBar(0, length(k), style = 3)
    o <- vector("list", length(k))
    for (ii in seq_along(k)) { o[[ii]] <- one(k[ii]); utils::setTxtProgressBar(pb, ii) }
    close(pb); o
  } else lapply(k, one)

  df <- as.data.frame(do.call(rbind, res))
  df <- cbind(cell = k, df)
  structure(df, class = c("si_metrics", "data.frame"),
            template = field$template, field_cells = field$cells[k],
            scenario = field$scenario)
}

#' @export
print.si_metrics <- function(x, ...) {
  cat("<si_metrics>  scenario: ", attr(x, "scenario"), "   cells: ", nrow(x),
      "\n", sep = "")
  print(utils::head(as.data.frame(x), 6))
  if (nrow(x) > 6) cat("  ...\n")
  invisible(x)
}

#' @export
summary.si_metrics <- function(object, ...) {
  d <- as.data.frame(object)[, -1, drop = FALSE]
  data.frame(metric = names(d),
             mean = sapply(d, mean, na.rm = TRUE),
             sd   = sapply(d, stats::sd, na.rm = TRUE),
             min  = sapply(d, min, na.rm = TRUE),
             max  = sapply(d, max, na.rm = TRUE),
             n_na = sapply(d, function(z) sum(is.na(z))),
             row.names = NULL)
}

#' Rasterise per-cell metrics
#'
#' @param metrics An [si_metrics()] result.
#' @param what Character vector of metric names, or `"all"`.
#' @return A `SpatRaster`, one layer per metric.
#' @examples
#' library(terra)
#' r <- rast(nrows = 12, ncols = 12, nlyrs = 6); values(r) <- runif(ncell(r) * 6)
#' names(r) <- paste0("p", 1:6)
#' r2 <- rast(nrows = 12, ncols = 12, nlyrs = 5); values(r2) <- runif(ncell(r2) * 5)
#' names(r2) <- paste0("a", 1:5)
#' f <- si_overlap(si_stack(r), si_stack(r2),
#'                 si_simulate_metaweb(6, 5, 0.3, seed = 1,
#'                                     names_A = paste0("p", 1:6),
#'                                     names_B = paste0("a", 1:5)),
#'                 floor = 0.05, check = FALSE)
#' m <- si_metrics(f, what = c("links", "connectance"))
#' plot(si_metric_map(m))
#' @export
si_metric_map <- function(metrics, what = "all") {
  tmpl <- si_grid_rast(attr(metrics, "template"))
  if (is.null(tmpl))
    stop("These metrics carry no raster template (the field was built from a ",
         "matrix).", call. = FALSE)
  d <- as.data.frame(metrics)
  nm <- setdiff(names(d), "cell")
  if (!identical(what, "all")) nm <- intersect(what, nm)
  cells <- attr(metrics, "field_cells")
  lst <- lapply(nm, function(v) {
    r <- terra::rast(tmpl); terra::values(r) <- NA_real_
    r[cells] <- d[[v]]
    names(r) <- v
    r
  })
  terra::rast(lst)
}
