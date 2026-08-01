# ---------------------------------------------------------------------------
# Module 4b: null models, species-level metrics, interoperability
# ---------------------------------------------------------------------------

#' Null models and standardised effect sizes
#'
#' Compares an observed network with an ensemble of randomisations that
#' preserve its marginal totals, and returns standardised effect sizes. This
#' matters especially in spatial projections: a raw metric map largely tracks
#' local species richness, whereas a z-score map shows where the network is
#' structured differently from what its own marginals imply.
#'
#' @param m An interaction matrix, or an [si_field()] with `cell` supplied.
#' @param cell Integer cell index when `m` is a field.
#' @param what Character vector of metrics (see [si_web_metrics()]).
#' @param model Null model. `"patefield"` redistributes interaction weights at
#'   random while fixing both sets of marginal totals (`stats::r2dtable`);
#'   `"curveball"` is the fixed-fixed binary swap algorithm of Strona and
#'   others (2014), which preserves both degree sequences exactly;
#'   `"shuffle"` preserves only the number of links.
#' @param n Number of randomisations.
#' @param seed Optional random seed.
#'
#' @return A `data.frame` with the observed value, the null mean and standard
#'   deviation, the z-score and a two-sided p-value for each metric.
#'
#' @details
#' A structurally degenerate field cannot produce a meaningful effect size: if
#' the observed matrix is the outer product of its own marginals, every
#' marginal-preserving randomisation reproduces it, the null standard deviation
#' is zero and the z-score is undefined. `si_null()` returns `NA` in that case
#' rather than a spurious number, which is itself a useful diagnostic.
#'
#' @examples
#' set.seed(1)
#' mw <- si_simulate_metaweb(20, 15, 0.25, architecture = "nested", seed = 1)
#' w  <- mw$matrix * matrix(rpois(300, 4) + 1, 20, 15)
#' si_null(w, what = c("NODF", "H2prime"), n = 99, seed = 1)
#'
#' # a matrix that is exactly the outer product of its marginals has no
#' # specialisation to detect: the observed value sits on the null mean
#' op <- round(outer(runif(12) * 20, runif(10) * 20)) + 1
#' si_null(op, what = "H2prime", n = 99, seed = 2)
#' @references
#' Strona, G., Nappo, D., Boccacci, F., Fattorini, S. & San-Miguel-Ayanz, J.
#' (2014) A fast and unbiased procedure to randomize ecological binary
#' matrices with fixed row and column totals. *Nature Communications*, 5, 4114.
#'
#' Patefield, W.M. (1981) Algorithm AS 159: an efficient method of generating
#' random R x C tables with given row and column totals. *Applied Statistics*,
#' 30, 91--97.
#' @seealso [si_web_metrics()], [si_degenerate()]
#' @export
si_null <- function(m, cell = NULL, what = c("NODF", "H2prime", "evenness"),
                    model = c("patefield", "curveball", "shuffle"),
                    n = 199, seed = NULL) {
  model <- match.arg(model)
  if (inherits(m, "si_field")) {
    if (is.null(cell)) stop("Supply `cell` when `m` is an si_field.", call. = FALSE)
    m <- .local_matrix(m, cell, drop = TRUE)
  }
  W <- as.matrix(m); W[is.na(W)] <- 0
  if (nrow(W) < 2 || ncol(W) < 2 || sum(W) == 0)
    return(data.frame(metric = what, observed = NA_real_, null_mean = NA_real_,
                      null_sd = NA_real_, z = NA_real_, p = NA_real_))
  if (!is.null(seed)) set.seed(seed)
  obs <- si_web_metrics(W, what)
  gen <- switch(model,
    patefield = {
      Wi <- round(W / min(W[W > 0]))            # integerise, preserving shape
      r <- rowSums(Wi); cc <- colSums(Wi)
      if (sum(r) != sum(cc) || sum(r) == 0) NULL else
        function() stats::r2dtable(1, r, cc)[[1]]
    },
    curveball = function() .curveball(W > 0, n_iter = 5 * sum(W > 0)) * 1,
    shuffle   = function() {
      z <- numeric(length(W)); z[sample.int(length(W), sum(W > 0))] <- W[W > 0]
      matrix(z, nrow(W), ncol(W))
    })
  if (is.null(gen))
    return(data.frame(metric = what, observed = obs, null_mean = NA_real_,
                      null_sd = NA_real_, z = NA_real_, p = NA_real_,
                      row.names = NULL))
  sim <- vapply(seq_len(n), function(i) si_web_metrics(gen(), what),
                numeric(length(obs)))
  if (is.null(dim(sim))) sim <- matrix(sim, nrow = 1)
  mu <- rowMeans(sim, na.rm = TRUE)
  sd_ <- apply(sim, 1, stats::sd, na.rm = TRUE)
  z <- ifelse(sd_ > .Machine$double.eps^0.5, (obs - mu) / sd_, NA_real_)
  p <- vapply(seq_along(obs), function(i) {
    if (is.na(obs[i])) return(NA_real_)
    2 * min(mean(sim[i, ] >= obs[i], na.rm = TRUE),
            mean(sim[i, ] <= obs[i], na.rm = TRUE))
  }, numeric(1))
  data.frame(metric = names(obs), observed = unname(obs), null_mean = mu,
             null_sd = sd_, z = z, p = pmin(1, p), row.names = NULL)
}

.curveball <- function(b, n_iter = 1000L) {
  b <- b > 0
  nr <- nrow(b); nc <- ncol(b)
  hp <- lapply(seq_len(nr), function(i) which(b[i, ]))
  for (k in seq_len(n_iter)) {
    ab <- sample.int(nr, 2L)
    A <- hp[[ab[1]]]; B <- hp[[ab[2]]]
    ints <- intersect(A, B); li <- length(ints)
    if (length(A) == li || length(B) == li) next
    pool <- sample(setdiff(union(A, B), ints))
    nA <- length(A) - li
    hp[[ab[1]]] <- c(ints, pool[seq_len(nA)])
    hp[[ab[2]]] <- c(ints, pool[-seq_len(nA)])
  }
  out <- matrix(FALSE, nr, nc)
  for (i in seq_len(nr)) out[i, hp[[i]]] <- TRUE
  out
}


#' Species-level network metrics
#'
#' Per-species position in the local network: degree, normalised degree,
#' species strength sensu Bascompte and others (2006) -- the sum of the
#' dependences that partners place on this species --
#' the paired differences index of specialisation, and the raw
#' Kullback-Leibler specialisation `d`.
#'
#' @param m An interaction matrix, or an [si_field()] with `cell`.
#' @param cell Integer cell index when `m` is a field.
#' @param level `"A"` (rows) or `"B"` (columns).
#' @return A `data.frame`, one row per species.
#' @details
#' `d_raw` is the unnormalised Kullback-Leibler divergence between a species'
#' interaction distribution and the overall availability of partners
#' (Blüthgen and others, 2006). It is reported unnormalised because the
#' standardisation to \eqn{d'} requires a constrained optimisation that
#' `bipartite::dfun()` implements; use that package if the normalised index is
#' required.
#' @examples
#' f <- si_example_field()
#' head(si_species_metrics(f, cell = 3, level = "A"))
#' @references
#' Bascompte, J., Jordano, P. & Olesen, J.M. (2006) Asymmetric coevolutionary
#' networks facilitate biodiversity maintenance. *Science*, 312, 431--433.
#'
#' Poisot, T., Canard, E., Mouquet, N. & Hochberg, M.E. (2012) A comparative
#' study of ecological specialization estimators. *Methods in Ecology and
#' Evolution*, 3, 537--544.
#' @export
si_species_metrics <- function(m, cell = NULL, level = c("A", "B")) {
  level <- match.arg(level)
  if (inherits(m, "si_field")) {
    if (is.null(cell)) stop("Supply `cell` when `m` is an si_field.", call. = FALSE)
    m <- .local_matrix(m, cell, drop = TRUE)
  }
  W <- as.matrix(m); W[is.na(W)] <- 0
  if (level == "B") W <- t(W)
  if (!nrow(W) || !ncol(W) || sum(W) == 0)
    return(data.frame(species = rownames(W), degree = integer(0)))
  A <- rowSums(W); tot <- sum(W); q <- colSums(W) / tot
  dep <- W / pmax(A, 1e-12)                       # dependence of i on j
  pdi <- apply(W, 1, function(x) {
    x <- sort(x[x > 0], decreasing = TRUE)
    if (length(x) < 2) return(1)
    sum(x[1] - x[-1]) / ((length(x) - 1) * x[1])
  })
  dkl <- vapply(seq_len(nrow(W)), function(i) {
    p <- dep[i, ]; k <- p > 0 & q > 0
    if (!any(k)) return(NA_real_)
    sum(p[k] * log(p[k] / q[k]))
  }, numeric(1))
  data.frame(species = rownames(W),
             degree = rowSums(W > 0),
             normalised_degree = rowSums(W > 0) / ncol(W),
             species_strength = rowSums(sweep(W, 2, pmax(colSums(W), 1e-12), "/")),
             pdi = pdi, d_raw = dkl, row.names = NULL)
}


#' Convert a local network to igraph or to a bipartite-style matrix
#'
#' Interoperability with the established network toolchain. `spinet` computes
#' spatial fields; `igraph`, `tidygraph` and `bipartite` remain the right tools
#' for single-network analysis and plotting.
#'
#' @param field An [si_field()], or a matrix.
#' @param cell Cell index when `field` is a field.
#' @param drop Remove unlinked species.
#' @return `si_to_igraph()` returns an `igraph` object with a logical `type`
#'   vertex attribute (bipartite convention) and `weight` edge attribute;
#'   `si_to_bipartite()` returns a plain matrix with level A in rows, the
#'   orientation expected by the `bipartite` package.
#' @examples
#' f <- si_example_field()
#' w <- si_to_bipartite(f, cell = 2)
#' dim(w)
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   g <- si_to_igraph(f, cell = 2)
#'   c(vertices = igraph::gorder(g), edges = igraph::gsize(g))
#' }
#' @export
si_to_bipartite <- function(field, cell = 1, drop = TRUE) {
  if (!inherits(field, "si_field")) return(as.matrix(field))
  .local_matrix(field, cell, drop = drop)
}

#' @rdname si_to_bipartite
#' @export
si_to_igraph <- function(field, cell = 1, drop = TRUE) {
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Package 'igraph' is required for si_to_igraph().", call. = FALSE)
  W <- si_to_bipartite(field, cell, drop)
  nr <- nrow(W); nc <- ncol(W)
  if (!nr || !nc) stop("Empty network in this cell.", call. = FALSE)
  el <- which(W > 0, arr.ind = TRUE)
  g <- igraph::make_empty_graph(n = nr + nc, directed = FALSE)
  igraph::V(g)$name <- c(rownames(W), colnames(W))
  igraph::V(g)$type <- c(rep(FALSE, nr), rep(TRUE, nc))
  g <- igraph::add_edges(g, as.vector(t(cbind(el[, 1], el[, 2] + nr))))
  igraph::E(g)$weight <- W[el]
  g
}


#' Links that are newly possible in the future
#'
#' Identifies pairs permitted by the future constraint layer but not by the
#' present one, and those permitted at both times but only co-occurring in the
#' future. The first is genuine novelty in the interaction rules; the second is
#' novelty only in geography.
#'
#' @param current,future [si_field()] objects.
#' @param cells Optional cell subset.
#' @return A `data.frame` with per-cell counts of `novel_rule`,
#'   `novel_spatial` and `lost` links.
#' @examples
#' set.seed(3)
#' A1 <- si_simulate_species(10, 14, 14, prefix = "A", seed = 3)
#' B1 <- si_simulate_species(8, 14, 14, prefix = "B", seed = 4)
#' A2 <- si_simulate_climate(A1, "shift", distance = 3)
#' B2 <- si_simulate_climate(B1, "shift", distance = 3)
#' ph <- si_phenology_simulate(10, 8, shift_mean = -10, shift_sd = 40, seed = 3)
#' m1 <- si_metaweb(si_rule_phenology(ph, "current", 0.3)$matrix,
#'                  names_A = A1$names, names_B = B1$names)
#' m2 <- si_metaweb(si_rule_phenology(ph, "future", 0.3)$matrix,
#'                  names_A = A1$names, names_B = B1$names)
#' th <- function(s) si_threshold(s, "fixed", 0.3)
#' f1 <- si_overlap(th(A1), th(B1), m1, method = "binary", check = FALSE)
#' f2 <- si_overlap(th(A2), th(B2), m2, method = "binary", check = FALSE)
#' colSums(si_novelty(f1, f2)[, -1])
#' @references
#' Vizentin-Bugoni, J. and others (2019) Structure, spatial dynamics, and
#' stability of novel seed dispersal mutualistic networks in Hawai'i.
#' *Science*, 364, 78--82.
#' @export
si_novelty <- function(current, future, cells = NULL) {
  stopifnot(inherits(current, "si_field"), inherits(future, "si_field"))
  k <- if (is.null(cells)) seq_len(current$n_cells) else as.integer(cells)
  F1 <- current$F > 0; F2 <- future$F > 0
  rule_new <- F2 & !F1
  res <- vapply(k, function(i) {
    w1 <- .local_matrix(current, i, drop = FALSE) > 0
    w2 <- .local_matrix(future,  i, drop = FALSE) > 0
    c(novel_rule    = sum(w2 & !w1 & rule_new),
      novel_spatial = sum(w2 & !w1 & !rule_new),
      lost          = sum(w1 & !w2))
  }, numeric(3))
  df <- cbind(cell = k, as.data.frame(t(res)))
  structure(df, class = c("si_metrics", "data.frame"),
            template = current$template, field_cells = current$cells[k],
            scenario = "novelty")
}


#' Aggregate species into functional or taxonomic groups
#'
#' Collapses a stack, and optionally a constraint layer, from species to
#' groups. Grouping changes measured specialisation systematically, so any
#' analysis conducted at group level should state the rule used.
#'
#' @param x An [si_stack()] or [si_metaweb()].
#' @param groups Character or factor vector, one entry per species.
#' @param fun How to combine within a group: `"max"` (the default; a group is
#'   present where any member is), `"mean"` or `"sum"`.
#' @return An object of the same class with one column per group.
#' @examples
#' st <- si_stack(matrix(runif(200 * 6), 200, 6,
#'                       dimnames = list(NULL, paste0("sp", 1:6))))
#' g  <- c("Apidae", "Apidae", "Syrphidae", "Syrphidae", "Syrphidae", "Muscidae")
#' si_aggregate(st, g)
#' @references
#' Rahimi, E. & Jung, C. (2025) Impact of taxonomic and functional grouping on
#' specialization in plant-pollinator networks. *Entomological Research*, 55,
#' e70036.
#' @export
si_aggregate <- function(x, groups, fun = c("max", "mean", "sum")) {
  fun <- match.arg(fun)
  f <- switch(fun, max = function(z) apply(z, 1, max),
              mean = rowMeans, sum = rowSums)
  if (inherits(x, "si_stack")) {
    if (length(groups) != ncol(x$values))
      stop("`groups` must have one entry per species (", ncol(x$values), ").",
           call. = FALSE)
    lv <- unique(as.character(groups))
    v <- vapply(lv, function(g)
      f(x$values[, groups == g, drop = FALSE]), numeric(nrow(x$values)))
    x$values <- v; x$names <- lv; colnames(x$values) <- lv
    x$binary <- .is_binary(v)
    return(x)
  }
  if (inherits(x, "si_metaweb")) {
    M <- x$matrix
    lv <- unique(as.character(groups))
    if (length(groups) == nrow(M)) {
      M <- t(vapply(lv, function(g) {
        s <- M[groups == g, , drop = FALSE]
        if (fun == "max") apply(s, 2, max) else if (fun == "mean") colMeans(s) else colSums(s)
      }, numeric(ncol(M))))
      rownames(M) <- lv
    } else if (length(groups) == ncol(M)) {
      M <- vapply(lv, function(g) {
        s <- M[, groups == g, drop = FALSE]
        if (fun == "max") apply(s, 1, max) else if (fun == "mean") rowMeans(s) else rowSums(s)
      }, numeric(nrow(M)))
      colnames(M) <- lv
    } else stop("`groups` must match one dimension of the metaweb.", call. = FALSE)
    return(new_si_metaweb(M, all(M == 1)))
  }
  stop("`x` must be an si_stack or an si_metaweb.", call. = FALSE)
}
