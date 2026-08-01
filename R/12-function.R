# ---------------------------------------------------------------------------
# Module 7: interaction-based ecosystem function, redundancy and uncertainty
# ---------------------------------------------------------------------------

#' Interaction-based ecosystem function
#'
#' Converts a spatial interaction field into a per-cell measure of delivered
#' function. The default is the summed interaction weight, which for a
#' pollination network is a relative index of pollination supply; supplying
#' `effectiveness` and `dependence` turns it into a weighted service estimate.
#'
#' @param field An [si_field()].
#' @param effectiveness Optional matrix of per-link effectiveness, the same
#'   shape as the constraint layer. Not all visitors are equal pollinators, and
#'   this is where that is encoded.
#' @param dependence Optional numeric vector, one value per level-A species,
#'   giving its dependence on the interaction (for crops, the pollinator
#'   dependence of Klein and others 2007).
#' @param per_species Logical. Return the level-A species-by-cell matrix rather
#'   than the cell totals.
#' @param cells Optional cell subset.
#'
#' @return A `data.frame` of class `si_metrics` with per-cell `function_index`,
#'   `n_partners` and `redundancy`, or a matrix if `per_species = TRUE`.
#'
#' @details
#' The index is
#' \deqn{\Phi_k = \sum_{i} d_i \sum_{j} W_k[i,j] \, E[i,j],}
#' the sum over level-A species of their dependence-weighted realised
#' interaction strength. Values are relative, not absolute: they compare cells
#' and scenarios within one analysis and should not be read as visitation rates.
#'
#' @examples
#' f <- si_example_field(n_cells = 60)
#' fx <- si_function(f)
#' head(fx)
#'
#' # weight by crop dependence on animal pollination
#' dep <- c(rep(0.95, 5), rep(0.25, 5), rep(0.05, 5))
#' head(si_function(f, dependence = dep))
#' @references
#' Klein, A.-M., Vaissiere, B.E., Cane, J.H., Steffan-Dewenter, I.,
#' Cunningham, S.A., Kremen, C. & Tscharntke, T. (2007) Importance of
#' pollinators in changing landscapes for world crops. *Proceedings of the
#' Royal Society B*, 274, 303--313.
#' @seealso [si_service_change()], [si_metric_map()]
#' @export
si_function <- function(field, effectiveness = NULL, dependence = NULL,
                        per_species = FALSE, cells = NULL) {
  stopifnot(inherits(field, "si_field"))
  nA <- length(field$names_A)
  if (is.null(dependence)) dependence <- rep(1, nA)
  if (length(dependence) != nA)
    stop("`dependence` must have one value per ", field$interaction$level_A,
         " (", nA, ").", call. = FALSE)
  if (!is.null(effectiveness)) {
    effectiveness <- as.matrix(effectiveness)
    if (!identical(dim(effectiveness), dim(field$F)))
      stop("`effectiveness` must have the same dimensions as the constraint ",
           "layer (", nrow(field$F), " x ", ncol(field$F), ").", call. = FALSE)
  }
  k <- if (is.null(cells)) seq_len(field$n_cells) else as.integer(cells)
  sp <- matrix(0, length(k), nA, dimnames = list(NULL, field$names_A))
  npart <- red <- numeric(length(k))
  for (ii in seq_along(k)) {
    W <- .local_matrix(field, k[ii], drop = FALSE)
    if (!is.null(effectiveness)) W <- W * effectiveness
    sp[ii, ] <- rowSums(W) * dependence
    b <- W > 0
    npart[ii] <- mean(rowSums(b)[rowSums(b) > 0])
    p <- W / pmax(rowSums(W), 1e-12)
    red[ii] <- mean(apply(p, 1, function(z) {
      z <- z[z > 0]; if (length(z) < 1) NA_real_ else 1 / sum(z^2)
    }), na.rm = TRUE)
  }
  if (per_species) return(sp)
  df <- data.frame(cell = k, function_index = rowSums(sp),
                   n_partners = ifelse(is.finite(npart), npart, NA_real_),
                   redundancy = red)
  structure(df, class = c("si_metrics", "data.frame"),
            template = field$template, field_cells = field$cells[k],
            scenario = field$scenario)
}


#' Change in ecosystem function between two scenarios
#'
#' @param current,future [si_field()] objects.
#' @param ... Passed to [si_function()].
#' @return A `data.frame` of class `si_metrics` with present and future
#'   function, absolute and relative change. The sign convention follows the
#'   interaction's `overlap_sign`, so that for antagonistic networks an
#'   increase in interaction strength is reported as a loss.
#' @examples
#' f1 <- si_example_field(n_cells = 40, seed = 1)
#' f2 <- f1; f2$A$values <- f2$A$values * 0.6; f2$scenario <- "future"
#' head(si_service_change(f1, f2))
#' @export
si_service_change <- function(current, future, ...) {
  a <- si_function(current, ...); b <- si_function(future, ...)
  sgn <- current$interaction$overlap_sign
  df <- data.frame(cell = a$cell,
                   f_current = a$function_index, f_future = b$function_index,
                   change = sgn * (b$function_index - a$function_index),
                   rel_change = sgn * (b$function_index - a$function_index) /
                     pmax(a$function_index, 1e-12),
                   redundancy_change = b$redundancy - a$redundancy)
  structure(df, class = c("si_metrics", "data.frame"),
            template = attr(a, "template"), field_cells = attr(a, "field_cells"),
            scenario = "change")
}


#' Functional redundancy and response diversity
#'
#' Redundancy is the effective number of partners per level-A species
#' (the inverse Simpson index of its dependence distribution). Response
#' diversity is the dispersion of projected suitability change among a
#' species' partners: a plant whose pollinators all respond to climate in the
#' same direction has low response diversity and is vulnerable even if it has
#' many partners.
#'
#' @param field An [si_field()].
#' @param change_B Optional numeric vector of projected suitability change for
#'   each level-B species; required for response diversity.
#' @param cells Optional cell subset.
#' @return A `data.frame` of class `si_metrics`.
#' @examples
#' f <- si_example_field(n_cells = 40)
#' set.seed(1)
#' head(si_redundancy(f, change_B = rnorm(12, -0.3, 0.25)))
#' @references
#' Elmqvist, T., Folke, C., Nystrom, M., Peterson, G., Bengtsson, J.,
#' Walker, B. & Norberg, J. (2003) Response diversity, ecosystem change, and
#' resilience. *Frontiers in Ecology and the Environment*, 1, 488--494.
#' @export
si_redundancy <- function(field, change_B = NULL, cells = NULL) {
  stopifnot(inherits(field, "si_field"))
  nB <- length(field$names_B)
  if (!is.null(change_B) && length(change_B) != nB)
    stop("`change_B` must have one value per ", field$interaction$level_B,
         " (", nB, ").", call. = FALSE)
  k <- if (is.null(cells)) seq_len(field$n_cells) else as.integer(cells)
  out <- t(vapply(k, function(i) {
    W <- .local_matrix(field, i, drop = FALSE)
    rs <- rowSums(W)
    keep <- rs > 0
    if (!any(keep)) return(c(redundancy = NA_real_, response_diversity = NA_real_))
    p <- W[keep, , drop = FALSE] / rs[keep]
    redu <- mean(apply(p, 1, function(z) 1 / sum(z[z > 0]^2)))
    rd <- if (is.null(change_B)) NA_real_ else
      mean(apply(p, 1, function(z) {
        m <- sum(z * change_B)
        sqrt(sum(z * (change_B - m)^2))
      }))
    c(redundancy = redu, response_diversity = rd)
  }, numeric(2)))
  df <- cbind(cell = k, as.data.frame(out))
  structure(df, class = c("si_metrics", "data.frame"),
            template = field$template, field_cells = field$cells[k],
            scenario = field$scenario)
}


#' Uncertainty across an ensemble of fields
#'
#' Summarises a metric across an ensemble of interaction fields built from
#' different distribution models, thresholds, climate models or constraint
#' rules, producing per-cell mean, standard deviation, coefficient of
#' variation and agreement on the direction of change. Uncertainty maps should
#' accompany every projection, because the choice of threshold and of
#' constraint rule moves network results at least as much as the climate
#' scenario does.
#'
#' @param fields A list of [si_field()] objects covering the same cells.
#' @param what A single metric name, or `"function"` for [si_function()].
#' @param baseline Optional list of matching present-day fields; when supplied,
#'   the direction of change is assessed and `agreement` reports the proportion
#'   of ensemble members agreeing with the majority direction.
#' @param cells Optional cell subset.
#' @return A `data.frame` of class `si_metrics`.
#' @examples
#' set.seed(5)
#' base <- si_example_field(n_cells = 40, seed = 5)
#' ens <- lapply(c(0.5, 0.7, 0.9), function(s) {
#'   g <- base; g$A$values <- g$A$values * s; g
#' })
#' head(si_uncertainty(ens, what = "links", baseline = rep(list(base), 3)))
#' @export
si_uncertainty <- function(fields, what = "links", baseline = NULL,
                           cells = NULL) {
  if (!is.list(fields) || !length(fields))
    stop("`fields` must be a non-empty list of si_field objects.", call. = FALSE)
  if (!all(vapply(fields, inherits, logical(1), "si_field")))
    stop("Every element of `fields` must be an si_field.", call. = FALSE)
  k <- if (is.null(cells)) seq_len(fields[[1]]$n_cells) else as.integer(cells)
  grab <- function(f) {
    if (identical(what, "function")) si_function(f, cells = k)$function_index
    else suppressWarnings(si_metrics(f, what = what, cells = k))[[what]]
  }
  M <- vapply(fields, grab, numeric(length(k)))
  mu <- rowMeans(M, na.rm = TRUE)
  sd_ <- apply(M, 1, stats::sd, na.rm = TRUE)
  df <- data.frame(cell = k, mean = mu, sd = sd_,
                   cv = sd_ / pmax(abs(mu), 1e-12),
                   min = apply(M, 1, min, na.rm = TRUE),
                   max = apply(M, 1, max, na.rm = TRUE))
  if (!is.null(baseline)) {
    if (length(baseline) != length(fields))
      stop("`baseline` must have the same length as `fields`.", call. = FALSE)
    B <- vapply(baseline, grab, numeric(length(k)))
    D <- M - B
    df$mean_change <- rowMeans(D, na.rm = TRUE)
    df$agreement <- apply(D, 1, function(z) {
      z <- z[!is.na(z)]
      if (!length(z)) return(NA_real_)
      max(mean(z > 0), mean(z < 0), mean(z == 0))
    })
  }
  structure(df, class = c("si_metrics", "data.frame"),
            template = fields[[1]]$template,
            field_cells = fields[[1]]$cells[k], scenario = "ensemble")
}
