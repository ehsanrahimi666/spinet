# ---------------------------------------------------------------------------
# Module 5: change, turnover and rewiring
# ---------------------------------------------------------------------------

#' Partition interaction turnover into species turnover and rewiring
#'
#' Implements the additive partition
#' \deqn{\beta_{WN} = \beta_{ST} + \beta_{OS}}
#' of Poisot and others (2012), applied cell by cell between two time slices.
#' \eqn{\beta_{OS}} is interaction rewiring: dissimilarity of the links among
#' species present in both slices. \eqn{\beta_{ST}} is the part of interaction
#' turnover attributable to species turnover.
#'
#' @param current,future Two [si_field()] objects covering the same cells,
#'   species and constraint layer, differing in their suitability surfaces (and
#'   optionally in their constraint layer).
#' @param cells Optional integer vector of cells to evaluate.
#' @param index Dissimilarity index. `"whittaker"` (the default, used by Poisot
#'   and by CaraDonna and others) or `"sorensen"`.
#' @param cores Integer number of parallel workers.
#'
#' @return A `data.frame` of class `si_beta` with per-cell `beta_S`,
#'   `beta_WN`, `beta_OS`, `beta_ST`, the proportion `rewiring_share`, and link
#'   counts. The raster template is retained so the result can be mapped with
#'   [si_metric_map()].
#'
#' @section The zero-rewiring theorem:
#' If the constraint layer is the same in both slices, then two species present
#' in a cell at both times are linked at both times whenever the layer permits
#' it, so the link set among shared species is identical and
#' \eqn{\beta_{OS} = 0} exactly, in every cell, for every climate scenario.
#' Rewiring is therefore not detectable from range shifts alone: it requires a
#' constraint layer that itself changes, for example one recomputed from a
#' shifted phenology (see [si_rule_phenology()]). `si_beta()` warns when the two
#' fields share a constraint layer.
#'
#' @examples
#' set.seed(6)
#' st_now <- si_simulate_species(12, nrow = 20, ncol = 20, seed = 6)
#' st_fut <- si_simulate_climate(st_now, "decline", severity = 0.3)
#' po_now <- si_simulate_species(10, nrow = 20, ncol = 20, prefix = "poll",
#'                               seed = 7)
#' po_fut <- si_simulate_climate(po_now, "decline", severity = 0.3)
#' mw <- si_simulate_metaweb(12, 10, 0.25, seed = 6,
#'                           names_A = st_now$names, names_B = po_now$names)
#'
#' f1 <- si_overlap(st_now, po_now, mw, floor = 0.05, scenario = "current",
#'                  check = FALSE)
#' f2 <- si_overlap(st_fut, po_fut, mw, floor = 0.05, scenario = "future",
#'                  check = FALSE)
#' b <- si_beta(f1, f2)
#' colMeans(b[, c("beta_S", "beta_WN", "beta_ST", "beta_OS")], na.rm = TRUE)
#' # beta_OS is exactly zero: a static constraint layer cannot rewire
#' @references
#' Poisot, T., Canard, E., Mouillot, D., Mouquet, N. & Gravel, D. (2012) The
#' dissimilarity of species interaction networks. *Ecology Letters*, 15,
#' 1353--1361.
#'
#' CaraDonna, P.J. and others (2017) Interaction rewiring and the rapid
#' turnover of plant-pollinator networks. *Ecology Letters*, 20, 385--394.
#' @seealso [si_rewiring()], [si_network_extinction()]
#' @export
si_beta <- function(current, future, cells = NULL,
                    index = c("whittaker", "sorensen"), cores = 1) {
  index <- match.arg(index)
  stopifnot(inherits(current, "si_field"), inherits(future, "si_field"))
  if (current$n_cells != future$n_cells)
    stop("The two fields cover different numbers of cells.", call. = FALSE)
  if (!identical(current$names_A, future$names_A) ||
      !identical(current$names_B, future$names_B))
    stop("The two fields must contain the same species in the same order.",
         call. = FALSE)
  static <- identical(unname(current$F), unname(future$F))
  if (static)
    warning("Both fields share one constraint layer, so beta_OS (rewiring) ",
            "is identically zero by construction. See ?si_beta.",
            call. = FALSE)

  k <- if (is.null(cells)) seq_len(current$n_cells) else as.integer(cells)
  bfun <- if (index == "whittaker") .beta_whittaker else .beta_sorensen

  one <- function(i) {
    Lc <- .link_set(current, i); Lf <- .link_set(future, i)
    nc <- length(Lc$links); nf <- length(Lf$links)
    if (nc == 0 && nf == 0)
      return(c(beta_S = NA, beta_WN = NA, beta_ST = NA, beta_OS = NA,
               links_current = 0, links_future = 0, shared_links = 0))
    sc <- Lc$species; sf <- Lf$species
    bS <- bfun(sc, sf)
    bWN <- bfun(Lc$links, Lf$links)
    shp <- intersect(Lc$A, Lf$A); shb <- intersect(Lc$B, Lf$B)
    ck <- .cross_keys(shp, shb, Lc$nB)
    Oc <- Lc$links[Lc$links %in% ck]; Of <- Lf$links[Lf$links %in% ck]
    # No links among shared species at either time: nothing can rewire, so
    # all interaction turnover is species turnover (beta_OS = 0).
    bOS <- if (is.na(bWN)) NA_real_ else if (!length(Oc) && !length(Of)) 0 else bfun(Oc, Of)
    bST <- if (is.na(bWN)) NA_real_ else bWN - bOS
    c(beta_S = bS, beta_WN = bWN, beta_ST = bST, beta_OS = bOS,
      links_current = nc, links_future = nf,
      shared_links = length(intersect(Lc$links, Lf$links)))
  }
  res <- if (cores > 1 && .Platform$OS.type == "unix")
    parallel::mclapply(k, one, mc.cores = cores) else lapply(k, one)
  df <- as.data.frame(do.call(rbind, res))
  df$rewiring_share <- with(df, ifelse(is.na(beta_WN) | beta_WN == 0, NA,
                                       beta_OS / beta_WN))
  df <- cbind(cell = k, df)
  structure(df, class = c("si_beta", "data.frame"),
            template = current$template,
            field_cells = current$cells[k],
            static_constraint = static)
}

#' @export
print.si_beta <- function(x, ...) {
  cat("<si_beta>   cells: ", nrow(x), "\n", sep = "")
  m <- colMeans(as.data.frame(x)[, c("beta_S", "beta_WN", "beta_ST",
                                     "beta_OS")], na.rm = TRUE)
  cat(sprintf("  beta_S  %.4f   beta_WN %.4f   beta_ST %.4f   beta_OS %.4f\n",
              m[1], m[2], m[3], m[4]))
  if (isTRUE(attr(x, "static_constraint")))
    cat("  ! constraint layer is static: beta_OS == 0 by construction\n")
  invisible(x)
}

# Link and species sets are encoded as integers rather than strings: a link
# (i, j) becomes (i - 1) * nB + j, and a species becomes +i for level A and
# -j for level B. On a dense field this is roughly an order of magnitude
# faster and far lighter than pasted character keys.
.link_set <- function(field, k) {
  w <- .local_matrix(field, k, drop = FALSE)
  ij <- which(w > 0, arr.ind = TRUE)
  nB <- ncol(w)
  if (!nrow(ij))
    return(list(links = integer(0), species = integer(0),
                A = integer(0), B = integer(0), nB = nB))
  A <- sort(unique(ij[, 1])); B <- sort(unique(ij[, 2]))
  list(links = (ij[, 1] - 1L) * nB + ij[, 2],
       species = c(A, -B), A = A, B = B, nB = nB)
}

.cross_keys <- function(A, B, nB) {
  if (!length(A) || !length(B)) return(integer(0))
  as.vector(outer((A - 1L) * nB, B, "+"))
}

.beta_whittaker <- function(X, Y) {
  a <- length(intersect(X, Y)); b <- length(setdiff(Y, X)); cc <- length(setdiff(X, Y))
  s <- a + b + cc
  if (s == 0) return(NA_real_)
  s / ((2 * a + b + cc) / 2) - 1
}
.beta_sorensen <- function(X, Y) {
  a <- length(intersect(X, Y)); b <- length(setdiff(Y, X)); cc <- length(setdiff(X, Y))
  if (2 * a + b + cc == 0) return(NA_real_)
  (b + cc) / (2 * a + b + cc)
}


#' Rewiring diagnostics
#'
#' Summarises how much partner switching occurs between two slices and how much
#' *could* occur given the constraint layer. The ratio of the two is the
#' rewiring capacity that is actually used.
#'
#' @param current,future [si_field()] objects.
#' @param cells Optional cell subset.
#' @return A `data.frame` of class `si_rewiring` with per-cell `realised`
#'   (links gained among shared species), `lost` (links lost among shared
#'   species), `potential` (permitted but unrealised links among shared
#'   species in the current slice) and `capacity_used`.
#' @examples
#' set.seed(8)
#' ph  <- si_phenology_simulate(12, 10, shift_mean = -10, shift_sd = 30, seed = 8)
#' A1  <- si_simulate_species(12, 18, 18, prefix = "plant", seed = 8)
#' B1  <- si_simulate_species(10, 18, 18, prefix = "poll",  seed = 9)
#' A2  <- si_simulate_climate(A1, "decline", severity = 0.2)
#' B2  <- si_simulate_climate(B1, "decline", severity = 0.2)
#' mw1 <- si_metaweb(si_rule_phenology(ph, "current", 0.3)$matrix,
#'                   names_A = A1$names, names_B = B1$names)
#' mw2 <- si_metaweb(si_rule_phenology(ph, "future", 0.3)$matrix,
#'                   names_A = A1$names, names_B = B1$names)
#' f1 <- si_overlap(A1, B1, mw1, floor = 0.05, check = FALSE)
#' f2 <- si_overlap(A2, B2, mw2, floor = 0.05, check = FALSE)
#' colMeans(si_rewiring(f1, f2)[, -1], na.rm = TRUE)
#' @export
si_rewiring <- function(current, future, cells = NULL) {
  stopifnot(inherits(current, "si_field"), inherits(future, "si_field"))
  k <- if (is.null(cells)) seq_len(current$n_cells) else as.integer(cells)
  res <- lapply(k, function(i) {
    Lc <- .link_set(current, i); Lf <- .link_set(future, i)
    shp <- intersect(Lc$A, Lf$A); shb <- intersect(Lc$B, Lf$B)
    keys <- .cross_keys(shp, shb, Lc$nB)
    Oc <- Lc$links[Lc$links %in% keys]; Of <- Lf$links[Lf$links %in% keys]
    pot <- length(keys) - length(Oc)
    gained <- length(setdiff(Of, Oc)); lost <- length(setdiff(Oc, Of))
    c(realised = gained, lost = lost, potential = pot,
      capacity_used = if (pot > 0) gained / pot else NA_real_,
      shared_A = length(shp), shared_B = length(shb))
  })
  df <- cbind(cell = k, as.data.frame(do.call(rbind, res)))
  structure(df, class = c("si_rewiring", "data.frame"),
            template = current$template, field_cells = current$cells[k])
}


#' Locate networks that disappear
#'
#' Identifies cells that hold at least one interaction in the current slice and
#' none in the future one. In the published Chilean analyses this was the
#' dominant effect of climate change, and it remains meaningful even in fields
#' whose topological metrics are degenerate.
#'
#' @param current,future [si_field()] objects.
#' @param min_links Integer. A network is considered present when it has at
#'   least this many links.
#' @return A `data.frame` of class `si_extinction` with per-cell link counts and
#'   a `status` factor: `"persists"`, `"extinct"`, `"colonised"` or `"absent"`.
#' @examples
#' f <- si_example_field(n_cells = 80)
#' g <- f; g$A$values <- g$A$values * 0.5
#' tab <- si_network_extinction(f, g)
#' table(tab$status)
#' @export
si_network_extinction <- function(current, future, min_links = 1) {
  stopifnot(inherits(current, "si_field"), inherits(future, "si_field"))
  k <- seq_len(current$n_cells)
  lc <- vapply(k, function(i) sum(.local_matrix(current, i, drop = FALSE) > 0),
               numeric(1))
  lf <- vapply(k, function(i) sum(.local_matrix(future, i, drop = FALSE) > 0),
               numeric(1))
  status <- ifelse(lc >= min_links & lf >= min_links, "persists",
            ifelse(lc >= min_links & lf < min_links, "extinct",
            ifelse(lc < min_links & lf >= min_links, "colonised", "absent")))
  df <- data.frame(cell = k, links_current = lc, links_future = lf,
                   change = lf - lc,
                   status = factor(status, levels = c("persists", "extinct",
                                                      "colonised", "absent")))
  structure(df, class = c("si_extinction", "data.frame"),
            template = current$template, field_cells = current$cells)
}


#' Compare two fields metric by metric
#'
#' @param current,future [si_field()] objects.
#' @param what Metrics to compute, passed to [si_metrics()].
#' @param cores Parallel workers.
#' @return A `data.frame` of class `si_compare` with `_current`, `_future` and
#'   `_change` columns. The sign of `_change` is multiplied by the
#'   interaction's `overlap_sign`, so that for antagonistic or competitive
#'   networks a positive change consistently means deterioration.
#' @examples
#' f1 <- si_example_field(n_cells = 50, seed = 1)
#' f2 <- f1; f2$A$values <- f2$A$values * 0.6; f2$scenario <- "future"
#' head(si_compare(f1, f2, what = c("links", "connectance", "NODF")), 3)
#' @export
si_compare <- function(current, future, what = "all", cores = 1) {
  a <- si_metrics(current, what = what, cores = cores)
  b <- si_metrics(future,  what = what, cores = cores)
  vars <- setdiff(names(a), "cell")
  out <- data.frame(cell = a$cell)
  sgn <- current$interaction$overlap_sign
  for (v in vars) {
    out[[paste0(v, "_current")]] <- a[[v]]
    out[[paste0(v, "_future")]]  <- b[[v]]
    out[[paste0(v, "_change")]]  <- sgn * (b[[v]] - a[[v]])
  }
  structure(out, class = c("si_compare", "data.frame"),
            template = current$template, field_cells = current$cells,
            overlap_sign = sgn)
}

#' @export
`[.si_beta` <- function(x, ...) {
  y <- NextMethod("["); if (is.data.frame(y)) class(y) <- "data.frame"; y
}
#' @export
`[.si_metrics` <- function(x, ...) {
  y <- NextMethod("["); if (is.data.frame(y)) class(y) <- "data.frame"; y
}
