# ---------------------------------------------------------------------------
# Module 2: the spatial interaction engine
# ---------------------------------------------------------------------------

#' Build a spatial interaction field
#'
#' The central function of the package. Combines two suitability stacks with a
#' constraint layer into a *lazy* representation of one interaction network per
#' grid cell. The local matrices are never all held in memory; they are
#' reconstructed on demand by [si_local()].
#'
#' @param A,B [si_stack()] objects, or anything [si_stack()] accepts, giving
#'   the two levels of the network.
#' @param metaweb An [si_metaweb()] or [si_forbidden()] object, a matrix, or
#'   `NULL`. `NULL` means "every pair may interact", which produces a
#'   structurally degenerate field; the function warns in that case.
#' @param method Character. How the two suitabilities are combined into a local
#'   interaction weight:
#'   \describe{
#'     \item{`"product"`}{\eqn{w = p_i a_j}. Probability of joint occurrence
#'       under independence. Follows the continuous framework of Rahimi & Jung.}
#'     \item{`"binary"`}{\eqn{w = 1} where both exceed zero. The presence-absence
#'       framework.}
#'     \item{`"min"`}{\eqn{w = \min(p_i, a_j)}. A Liebig-style limiting rule:
#'       the interaction is capped by the scarcer partner.}
#'     \item{`"geometric"`}{\eqn{w = \sqrt{p_i a_j}}, on the same scale as the
#'       inputs, unlike `"product"`.}
#'     \item{`"harmonic"`}{\eqn{w = 2 p_i a_j / (p_i + a_j)}. Penalises
#'       asymmetry more strongly than `"min"`.}
#'     \item{`"difference"`}{\eqn{w = 1 - |p_i - a_j|}. A matching rule for
#'       niche-overlap framing.}
#'   }
#'   `"product"`, `"geometric"` and `"harmonic"` are all rank-one transforms
#'   and share the degeneracy properties described in [si_degenerate()];
#'   `"min"`, `"difference"` and `"binary"` are not.
#' @param interaction An [si_interaction()] object. Defaults to a
#'   plant-pollinator mutualism.
#' @param scenario Character label for this time slice, e.g. `"current"`.
#' @param floor Numeric. Interaction weights below this value are set to zero.
#'   For continuous SDM output this is essential: MaxEnt-type surfaces are
#'   never exactly zero, so without a floor every permitted link is "present"
#'   in every cell and all presence-based metrics collapse to the metaweb. See
#'   [si_prune()].
#' @param check Logical. Run [si_degenerate()] and warn about metrics that will
#'   be constant.
#'
#' @return An object of class `si_field`.
#'
#' @details
#' The local interaction matrix for cell \eqn{k} is
#' \deqn{W_k = f(p_k, a_k) \circ F}
#' where \eqn{f} is the chosen `method`, \eqn{\circ} is the Hadamard product
#' and \eqn{F} is the constraint layer. This is the spatially explicit,
#' climate-projectable form of the spatial-overlap term \eqn{S} in the
#' probabilistic framework of Vazquez and others (2009).
#'
#' Memory: a field stores \eqn{O(n_{cells}(n_A + n_B) + n_A n_B)} numbers
#' rather than \eqn{O(n_{cells} n_A n_B)}. For 2800 cells and 187 x 171 species
#' that is roughly 8 MB rather than 740 MB.
#'
#' @examples
#' set.seed(42)
#' # two levels and a sparse metaweb
#' P  <- si_stack(matrix(runif(300 * 12), 300, 12,
#'                dimnames = list(NULL, paste0("p", 1:12))), level = "plant")
#' A  <- si_stack(matrix(runif(300 * 10), 300, 10,
#'                dimnames = list(NULL, paste0("a", 1:10))), level = "pollinator")
#' mw <- si_simulate_metaweb(12, 10, connectance = 0.25, seed = 1)
#'
#' fld <- si_overlap(P, A, metaweb = mw, method = "product", floor = 0.05)
#' fld
#'
#' # the local network in cell 7
#' si_local(fld, 7)
#'
#' # what happens without a constraint layer
#' bad <- si_overlap(P, A, metaweb = NULL, method = "product", check = FALSE)
#' si_degenerate(bad)
#' @references
#' Vazquez, D.P., Chacoff, N.P. & Cagnolo, L. (2009) Evaluating multiple
#' determinants of the structure of plant-animal mutualistic networks.
#' *Ecology*, 90, 2039--2046.
#' @seealso [si_local()], [si_metrics()], [si_degenerate()], [si_prune()]
#' @export
si_overlap <- function(A, B, metaweb = NULL,
                       method = c("product", "binary", "min", "geometric",
                                  "harmonic", "difference"),
                       interaction = si_interaction("pollinatedBy",
                                                    level_A = "plant",
                                                    level_B = "pollinator"),
                       scenario = "current", floor = 0, check = TRUE) {

  method <- match.arg(method)
  if (!inherits(A, "si_stack")) A <- si_stack(A, level = interaction$level_A)
  if (!inherits(B, "si_stack")) B <- si_stack(B, level = interaction$level_B)
  if (!inherits(interaction, "si_interaction"))
    stop("`interaction` must be an si_interaction object; see ?si_interaction",
         call. = FALSE)
  if (!identical(as.integer(A$cells), as.integer(B$cells))) {
    same_grid <- !is.null(A$template) && !is.null(B$template) &&
      identical(A$template$dim, B$template$dim) &&
      isTRUE(all.equal(A$template$ext, B$template$ext))
    if (!same_grid)
      stop("The two stacks cover different cells (", nrow(A$values), " vs ",
           nrow(B$values), ") and do not share a raster grid, so their cells ",
           "cannot be aligned.\n  Build every stack with the same `mask` so ",
           "that cell i means the same place in each.", call. = FALSE)
    # Same grid, different valid cells: keep the cells present in both, in
    # raster order, so that row k refers to the same place in each stack.
    common <- sort(intersect(A$cells, B$cells))
    if (!length(common))
      stop("The two stacks share a grid but have no cells in common.",
           call. = FALSE)
    message("Aligning stacks on ", length(common), " shared cells (",
            nrow(A$values) - length(common), " level-A and ",
            nrow(B$values) - length(common), " level-B cells dropped).")
    A$values <- A$values[match(common, A$cells), , drop = FALSE]; A$cells <- common
    B$values <- B$values[match(common, B$cells), , drop = FALSE]; B$cells <- common
  }
  if (!is.numeric(floor) || length(floor) != 1L || floor < 0)
    stop("`floor` must be a single non-negative number.", call. = FALSE)

  F <- .resolve_metaweb(metaweb, A$names, B$names)

  fld <- structure(list(
    A = A, B = B, F = F,
    method = method, floor = floor,
    interaction = interaction, scenario = scenario,
    names_A = A$names, names_B = B$names,
    n_cells = nrow(A$values),
    template = si_grid(if (!is.null(A$template)) A$template else B$template),
    cells = A$cells
  ), class = "si_field")

  if (check) {
    d <- si_degenerate(fld, warn = TRUE)
    fld$degeneracy <- d
  }
  fld
}

.resolve_metaweb <- function(metaweb, nA, nB) {
  if (is.null(metaweb)) {
    m <- matrix(1, length(nA), length(nB), dimnames = list(nA, nB))
    attr(m, "unconstrained") <- TRUE
    return(m)
  }
  if (inherits(metaweb, "si_forbidden")) metaweb <- metaweb$metaweb
  m <- .constraint_matrix(metaweb)
  storage.mode(m) <- "double"
  if (is.null(rownames(m))) rownames(m) <- nA
  if (is.null(colnames(m))) colnames(m) <- nB
  rownames(m) <- .norm_names(rownames(m)); colnames(m) <- .norm_names(colnames(m))
  ri <- match(nA, rownames(m)); ci <- match(nB, colnames(m))
  if (anyNA(ri) || anyNA(ci)) {
    # If the metaweb has the right shape but shares no names at all, it was
    # almost certainly built without names; positional matching is unambiguous.
    dims_ok <- nrow(m) == length(nA) && ncol(m) == length(nB)
    no_overlap <- !any(nA %in% rownames(m)) && !any(nB %in% colnames(m))
    if (dims_ok && no_overlap) {
      message("Metaweb species names do not match the stacks but the ",
              "dimensions do; matching positionally.")
      dimnames(m) <- list(nA, nB)
      attr(m, "unconstrained") <- all(m == 1)
      return(m)
    }
    ma <- si_match_names(nA, rownames(m)); mb <- si_match_names(nB, colnames(m))
    stop("Metaweb does not cover every species in the stacks.\n",
         "  missing level-A: ", length(ma$only_a), " e.g. ",
         paste(utils::head(ma$only_a, 3), collapse = ", "), "\n",
         "  missing level-B: ", length(mb$only_a), " e.g. ",
         paste(utils::head(mb$only_a, 3), collapse = ", "), "\n",
         "  Rebuild it with si_metaweb(x, names_A = ..., names_B = ...).",
         call. = FALSE)
  }
  m <- m[ri, ci, drop = FALSE]
  attr(m, "unconstrained") <- all(m == 1)
  m
}

#' @export
print.si_field <- function(x, ...) {
  cat("<si_field>  scenario: ", x$scenario, "\n", sep = "")
  cat("  interaction : ", x$interaction$type, " (", x$interaction$sign, ")\n", sep = "")
  cat("  levels      : ", length(x$names_A), " ", x$interaction$level_A,
      "  x  ", length(x$names_B), " ", x$interaction$level_B, "\n", sep = "")
  cat("  cells       : ", x$n_cells, "\n", sep = "")
  cat("  method      : ", x$method,
      if (x$floor > 0) paste0("   floor = ", x$floor) else "", "\n", sep = "")
  cat("  constraint  : ", if (isTRUE(attr(x$F, "unconstrained")))
    "NONE (unconstrained - see ?si_degenerate)" else
      sprintf("metaweb, connectance %.4f", mean(x$F > 0)), "\n", sep = "")
  invisible(x)
}

#' Extract the local network of one or more cells
#'
#' Materialises the interaction matrix \eqn{W_k} for specific cells. This is
#' the only point at which a full local matrix exists in memory.
#'
#' @param field An [si_overlap()] field.
#' @param cell Integer vector of cell indices (positions within the field, not
#'   raster cell numbers), or a two-column matrix / `SpatVector` of coordinates.
#' @param drop Logical. Remove species with no links from the returned matrix.
#'   Network metrics are conventionally computed on the reduced matrix.
#' @param simplify Logical. When one cell is requested, return the matrix
#'   itself rather than a list of length one.
#' @return A matrix, or a named list of matrices.
#' @examples
#' set.seed(3)
#' f <- si_example_field()
#' w <- si_local(f, 5)
#' dim(w)
#' round(w[1:4, 1:4], 3)
#'
#' # several cells at once
#' length(si_local(f, c(1, 2, 3)))
#' @export
si_local <- function(field, cell, drop = TRUE, simplify = TRUE) {
  stopifnot(inherits(field, "si_field"))
  if (is.matrix(cell) || inherits(cell, "SpatVector")) cell <- .cells_from_xy(field, cell)
  cell <- as.integer(cell)
  bad <- cell < 1L | cell > field$n_cells
  if (any(bad))
    stop("Cell index out of range: ", paste(utils::head(cell[bad], 3),
         collapse = ", "), " (field has ", field$n_cells, " cells).",
         call. = FALSE)
  out <- lapply(cell, function(k) .local_matrix(field, k, drop))
  names(out) <- paste0("cell_", cell)
  if (simplify && length(out) == 1L) out[[1]] else out
}

.local_matrix <- function(field, k, drop = TRUE) {
  p <- field$A$values[k, ]
  a <- field$B$values[k, ]
  w <- switch(field$method,
    product    = outer(p, a, "*"),
    binary     = outer(as.numeric(p > 0), as.numeric(a > 0), "*"),
    min        = outer(p, a, pmin),
    geometric  = sqrt(outer(p, a, "*")),
    harmonic   = {
      s <- outer(p, a, "+"); pr <- outer(p, a, "*")
      ifelse(s > 0, 2 * pr / s, 0)
    },
    difference = pmax(0, 1 - abs(outer(p, a, "-"))) *
                 outer(as.numeric(p > 0), as.numeric(a > 0), "*")
  )
  w <- w * field$F
  if (field$floor > 0) w[w < field$floor] <- 0
  dimnames(w) <- list(field$names_A, field$names_B)
  if (drop) {
    rk <- rowSums(w > 0) > 0; ck <- colSums(w > 0) > 0
    w <- w[rk, ck, drop = FALSE]
  }
  w
}

.cells_from_xy <- function(field, xy) {
  if (is.null(field$template))
    stop("This field has no raster template; supply integer cell indices.",
         call. = FALSE)
  cn <- terra::cellFromXY(field$template, if (inherits(xy, "SpatVector"))
    terra::crds(xy) else xy)
  k <- match(cn, field$cells)
  if (anyNA(k))
    stop("Some coordinates fall outside the analysed cells.", call. = FALSE)
  k
}


#' Threshold a continuous stack or field
#'
#' Converts continuous suitability into presence-absence. Because the choice of
#' threshold is subjective and materially changes network results, the function
#' always reports the resulting prevalence so that the sensitivity of downstream
#' conclusions can be checked.
#'
#' @param x An [si_stack()] or [si_field()].
#' @param rule One of `"fixed"`, `"quantile"`, `"prevalence"` or `"mean"`.
#' @param value Numeric. The threshold for `"fixed"`, the quantile for
#'   `"quantile"`, or the target prevalence for `"prevalence"`. Ignored for
#'   `"mean"` (which uses each species' own mean suitability).
#' @param per_species Logical. Apply the rule to each species separately
#'   (the default) or to the whole stack at once.
#' @return An object of the same class, with binary values.
#' @examples
#' st <- si_stack(matrix(runif(500 * 5), 500, 5), level = "plant")
#' b1 <- si_threshold(st, "fixed", 0.5)
#' b2 <- si_threshold(st, "quantile", 0.9)
#' c(fixed = mean(b1$values), quantile90 = mean(b2$values))
#' @export
si_threshold <- function(x, rule = c("fixed", "quantile", "prevalence", "mean"),
                         value = 0.5, per_species = TRUE) {
  rule <- match.arg(rule)
  if (inherits(x, "si_field")) {
    x$A <- si_threshold(x$A, rule, value, per_species)
    x$B <- si_threshold(x$B, rule, value, per_species)
    x$method <- "binary"
    return(x)
  }
  stopifnot(inherits(x, "si_stack"))
  v <- x$values
  thr <- switch(rule,
    fixed      = rep(value, ncol(v)),
    mean       = colMeans(v),
    quantile   = if (per_species) apply(v, 2, stats::quantile, probs = value)
                 else rep(stats::quantile(v, value), ncol(v)),
    prevalence = if (per_species) apply(v, 2, stats::quantile, probs = 1 - value)
                 else rep(stats::quantile(v, 1 - value), ncol(v)))
  # Presence requires positive suitability. Without this, a quantile- or
  # prevalence-based threshold for a species absent from most of the grid
  # resolves to 0 and every cell, including those with zero suitability,
  # would be scored as present.
  out <- (v >= matrix(thr, nrow(v), ncol(v), byrow = TRUE) & v > 0) * 1
  zero_thr <- x$names[thr <= 0]
  if (length(zero_thr) && rule %in% c("quantile", "prevalence"))
    warning(length(zero_thr), " species have a ", rule, " threshold of zero ",
            "because most of their cells have zero suitability (e.g. ",
            paste(utils::head(zero_thr, 3), collapse = ", "),
            "); presence was restricted to cells with positive suitability.",
            call. = FALSE)
  dim(out) <- dim(v)
  dimnames(out) <- list(NULL, x$names)
  x$values <- out
  x$binary <- TRUE
  x
}


#' Remove negligible interaction weights
#'
#' Sets to zero every interaction weight below a floor. For continuous SDM
#' output this is not cosmetic. Correlative suitability surfaces are strictly
#' positive almost everywhere, so without a floor `W > 0` is true wherever the
#' metaweb permits a link, and every presence-based metric returns the metaweb
#' itself, identically in every cell and every time slice.
#'
#' @param field An [si_field()].
#' @param floor Numeric floor, or `"quantile"` with `q` to use a data-driven
#'   floor.
#' @param q Numeric quantile used when `floor = "quantile"`.
#' @return The field with an updated floor.
#' @examples
#' f  <- si_example_field(floor = 0)
#' f2 <- si_prune(f, 0.05)
#' c(before = mean(si_local(f, 1, drop = FALSE) > 0),
#'   after  = mean(si_local(f2, 1, drop = FALSE) > 0))
#' @export
si_prune <- function(field, floor = 0.05, q = 0.5) {
  stopifnot(inherits(field, "si_field"))
  if (identical(floor, "quantile")) {
    smp <- sample.int(field$n_cells, min(field$n_cells, 200))
    w <- unlist(lapply(smp, function(k) .local_matrix(field, k, drop = FALSE)))
    floor <- unname(stats::quantile(w[w > 0], q))
  }
  field$floor <- as.numeric(floor)
  field
}


#' A small ready-made field for examples and tests
#'
#' @param n_cells,n_A,n_B Dimensions.
#' @param connectance Metaweb connectance.
#' @param method,floor Passed to [si_overlap()].
#' @param seed Random seed.
#' @return An [si_field()].
#' @examples
#' si_example_field()
#' @export
si_example_field <- function(n_cells = 200, n_A = 15, n_B = 12,
                             connectance = 0.25, method = "product",
                             floor = 0.02, seed = 1) {
  set.seed(seed)
  A <- si_stack(matrix(stats::runif(n_cells * n_A), n_cells, n_A,
                       dimnames = list(NULL, paste0("plant", seq_len(n_A)))),
                level = "plant")
  B <- si_stack(matrix(stats::runif(n_cells * n_B), n_cells, n_B,
                       dimnames = list(NULL, paste0("poll", seq_len(n_B)))),
                level = "pollinator")
  mw <- si_simulate_metaweb(n_A, n_B, connectance = connectance, seed = seed,
                            names_A = A$names, names_B = B$names)
  si_overlap(A, B, mw, method = method, floor = floor, check = FALSE)
}
