# ---------------------------------------------------------------------------
# Module 6: robustness and coextinction
# ---------------------------------------------------------------------------

#' Network robustness to sequential species loss
#'
#' Computes robustness \eqn{R}, the area under the attack-tolerance curve, and
#' optionally \eqn{R_w}, its extension in which a species that loses a partner
#' may rewire to a surviving one with a user-specified probability. The
#' cascade rule is taken from the interaction type, so the same function
#' handles mutualistic coextinction, consumer-only extinction in antagonistic
#' webs, and competitive release.
#'
#' @param web A weighted or binary interaction matrix, or an [si_field()] with
#'   `cell` supplied.
#' @param cell Integer cell index when `web` is a field.
#' @param remove Which level to remove species from: `"A"` (rows) or `"B"`
#'   (columns). Robustness is then measured on the other level.
#' @param sequence Extinction order. `"random"`, `"degree"` (most generalised
#'   first), `"abundance"` (fewest interactions first), or a numeric vector
#'   giving a custom order. `"climate"` requires `suitability_change`.
#' @param suitability_change Numeric vector, one value per species in the
#'   removed level, giving the projected change in climatic suitability.
#'   Species are then removed from the largest decline to the smallest, which
#'   is the climate-informed extinction sequence rather than an arbitrary one.
#' @param rewiring Rewiring probability matrix of the same dimensions as `web`,
#'   giving the probability that a species rewires to each potential partner.
#'   `NULL` disables rewiring and returns the classic \eqn{R}. A convenient
#'   source is an [si_forbidden()] layer built from traits and phenology.
#' @param attempts Integer. Number of rewiring attempts allowed per lost link.
#' @param threshold Numeric on 0--1. Proportion of a species' original
#'   interaction strength that must be lost before it goes extinct. The
#'   classic model uses `1` (all partners lost); empirical evidence suggests
#'   lower values are more realistic.
#' @param n_sim Number of simulation replicates.
#' @param interaction An [si_interaction()] controlling the cascade rule.
#' @param seed Optional random seed.
#'
#' @return An object of class `si_robustness` with the mean robustness, the
#'   per-replicate values, and the mean attack-tolerance curve.
#'
#' @details
#' \eqn{R} is the area under the curve relating the surviving fraction of one
#' level to the removed fraction of the other. \eqn{R_w} follows the same
#' construction but allows partner switching, following Vizentin-Bugoni and
#' others (2020). The effect size \eqn{D = R_w - R} quantifies how much
#' rewiring buffers the network. Setting `rewiring` to a trait- or
#' phenology-constrained matrix rather than a matrix of ones is what
#' distinguishes attainable from unconstrained rewiring.
#'
#' @examples
#' set.seed(11)
#' mw <- si_simulate_metaweb(25, 20, connectance = 0.2, architecture = "nested",
#'                           seed = 11)
#' w  <- mw$matrix * matrix(runif(500, 1, 10), 25, 20)
#'
#' # classic robustness, no rewiring
#' r0 <- si_robustness(w, remove = "A", sequence = "random", n_sim = 50)
#' r0
#'
#' # unconstrained rewiring: the optimistic bound
#' r1 <- si_robustness(w, remove = "A", rewiring = matrix(1, 25, 20),
#'                     n_sim = 50, seed = 1)
#'
#' # trait-constrained rewiring: what is actually attainable
#' rew <- si_rule_morphology(runif(25, 5, 40), runif(20, 3, 35),
#'                           rule = "gower")$matrix
#' r2 <- si_robustness(w, remove = "A", rewiring = rew, n_sim = 50, seed = 1)
#'
#' c(R = r0$robustness, Rw_free = r1$robustness, Rw_constrained = r2$robustness,
#'   D_free = r1$robustness - r0$robustness,
#'   D_constrained = r2$robustness - r0$robustness)
#'
#' # an antagonistic web: the resource level is released, not coextinct
#' si_robustness(w, remove = "B", interaction = si_interaction("+-"),
#'               n_sim = 20)$robustness
#' @references
#' Burgos, E. and others (2007) Why nestedness in mutualistic networks?
#' *Journal of Theoretical Biology*, 249, 307--313.
#'
#' Memmott, J., Waser, N.M. & Price, M.V. (2004) Tolerance of pollination
#' networks to species extinctions. *Proceedings of the Royal Society B*, 271,
#' 2605--2611.
#'
#' Vizentin-Bugoni, J., Debastiani, V.J., Bastazini, V.A.G., Maruyama, P.K. &
#' Sperry, J.H. (2020) Including rewiring in the estimation of the robustness
#' of mutualistic networks. *Methods in Ecology and Evolution*, 11, 106--116.
#' @seealso [si_robustness_map()], [si_interaction()]
#' @export
si_robustness <- function(web, cell = NULL, remove = c("A", "B"),
                          sequence = c("random", "degree", "abundance",
                                       "climate"),
                          suitability_change = NULL,
                          rewiring = NULL, attempts = 1, threshold = 1,
                          n_sim = 100,
                          interaction = si_interaction("pollinatedBy"),
                          seed = NULL) {
  remove <- match.arg(remove)
  if (inherits(web, "si_field")) {
    if (is.null(cell)) stop("Supply `cell` when `web` is an si_field.",
                            call. = FALSE)
    interaction <- web$interaction
    web <- .local_matrix(web, cell, drop = TRUE)
  }
  W <- as.matrix(web); W[is.na(W)] <- 0
  if (!is.character(sequence)) {
    order_fixed <- as.integer(sequence); sequence <- "custom"
  } else {
    sequence <- match.arg(sequence); order_fixed <- NULL
  }
  if (!is.null(seed)) set.seed(seed)

  if (remove == "B") { W <- t(W); if (!is.null(rewiring)) rewiring <- t(rewiring) }
  nR <- nrow(W); nC <- ncol(W)
  if (nR < 2 || nC < 2)
    return(structure(list(robustness = NA_real_, sims = NA_real_,
                          curve = NULL, remove = remove,
                          rewiring = !is.null(rewiring)),
                     class = "si_robustness"))

  rule <- si_cascade_rule(interaction)
  victim <- if (remove == "A") rule[["B"]] else rule[["A"]]
  if (victim != "extinct")
    return(structure(list(robustness = 1, sims = rep(1, n_sim), curve = NULL,
                          remove = remove, rewiring = !is.null(rewiring),
                          note = paste0("Cascade rule '", interaction$cascade,
                                        "': the surviving level is ", victim,
                                        ", so robustness is 1 by definition.")),
                     class = "si_robustness"))

  base_strength <- colSums(W)
  seq_order <- switch(sequence,
    degree     = order(rowSums(W > 0), decreasing = TRUE),
    abundance  = order(rowSums(W), decreasing = FALSE),
    climate    = {
      if (is.null(suitability_change))
        stop("`sequence = \"climate\"` needs `suitability_change`.",
             call. = FALSE)
      if (length(suitability_change) != nR)
        stop("`suitability_change` must have one value per removed species (",
             nR, ").", call. = FALSE)
      order(suitability_change, decreasing = FALSE)
    },
    custom     = order_fixed,
    random     = NULL)

  reps <- if (sequence %in% c("random")) n_sim else
    if (is.null(rewiring)) 1L else n_sim
  curves <- matrix(NA_real_, reps, nR + 1L)
  Rvals <- numeric(reps)

  for (s in seq_len(reps)) {
    ord <- if (is.null(seq_order)) sample.int(nR) else seq_order
    M <- W
    alive <- rep(TRUE, nC)
    curve <- numeric(nR + 1L); curve[1] <- 1
    for (t in seq_len(nR)) {
      i <- ord[t]
      lost <- M[i, ]
      M[i, ] <- 0
      if (!is.null(rewiring) && any(lost > 0)) {
        for (j in which(lost > 0 & alive)) {
          amount <- lost[j]
          for (att in seq_len(attempts)) {
            cand <- which(rowSums(M > 0) > 0)
            cand <- setdiff(cand, ord[seq_len(t)])
            if (!length(cand)) break
            p <- rewiring[cand, j]
            if (sum(p, na.rm = TRUE) <= 0) break
            pick <- cand[sample.int(length(cand), 1L, prob = p)]
            if (stats::runif(1) < rewiring[pick, j]) {
              M[pick, j] <- M[pick, j] + amount
              break
            }
          }
        }
      }
      remaining <- colSums(M)
      alive <- alive & (remaining > (1 - threshold) * base_strength |
                          remaining > 0 & threshold < 1) &
        remaining > 0
      curve[t + 1L] <- sum(alive) / nC
    }
    curves[s, ] <- curve
    Rvals[s] <- mean((curve[-1] + curve[-length(curve)]) / 2)
  }
  structure(list(robustness = mean(Rvals), sims = Rvals,
                 curve = colMeans(curves), remove = remove,
                 sequence = sequence, rewiring = !is.null(rewiring),
                 n_removed = nR, n_target = nC),
            class = "si_robustness")
}

#' @export
print.si_robustness <- function(x, ...) {
  cat("<si_robustness>\n")
  cat("  removing level : ", x$remove, "   sequence: ",
      if (is.null(x$sequence)) "-" else x$sequence, "\n", sep = "")
  cat("  rewiring       : ", if (isTRUE(x$rewiring)) "yes (Rw)" else "no (R)",
      "\n", sep = "")
  cat("  robustness     : ", format(round(x$robustness, 4), nsmall = 4),
      if (length(x$sims) > 1)
        sprintf("   (sd %.4f, n = %d)", stats::sd(x$sims), length(x$sims))
      else "", "\n", sep = "")
  if (!is.null(x$note)) cat("  note           : ", x$note, "\n", sep = "")
  invisible(x)
}

#' @export
plot.si_robustness <- function(x, ...) {
  if (is.null(x$curve)) { message("No attack-tolerance curve to plot."); return(invisible(x)) }
  n <- length(x$curve) - 1
  graphics::plot(seq(0, 1, length.out = n + 1), x$curve, type = "l", lwd = 2,
                 xlab = paste0("fraction of level ", x$remove, " removed"),
                 ylab = "fraction of partners surviving",
                 ylim = c(0, 1), ...)
  graphics::abline(0, -1, lty = 3, col = "grey60")
  graphics::legend("bottomleft", bty = "n",
                   legend = sprintf("%s = %.3f",
                                    if (isTRUE(x$rewiring)) "Rw" else "R",
                                    x$robustness))
  invisible(x)
}


#' Robustness across every cell of a field
#'
#' @param field An [si_field()].
#' @param cells Optional cell subset.
#' @param n_sim Replicates per cell (keep small: cost is cells x replicates).
#' @param ... Passed to [si_robustness()].
#' @return A `data.frame` of class `si_metrics`, mappable with
#'   [si_metric_map()].
#' @examples
#' f <- si_example_field(n_cells = 30)
#' rb <- si_robustness_map(f, n_sim = 10)
#' summary(rb)
#' @export
si_robustness_map <- function(field, cells = NULL, n_sim = 20, ...) {
  stopifnot(inherits(field, "si_field"))
  k <- if (is.null(cells)) seq_len(field$n_cells) else as.integer(cells)
  v <- vapply(k, function(i) {
    r <- si_robustness(field, cell = i, n_sim = n_sim,
                       interaction = field$interaction, ...)
    r$robustness
  }, numeric(1))
  df <- data.frame(cell = k, robustness = v)
  structure(df, class = c("si_metrics", "data.frame"),
            template = field$template, field_cells = field$cells[k],
            scenario = field$scenario)
}


#' Climate-informed coextinction
#'
#' Removes species in order of projected loss of climatic suitability rather
#' than at random or by degree, and reports the resulting cascade. This is the
#' extinction sequence implied by species distribution models, and it is what
#' makes a coextinction simulation a projection rather than a thought
#' experiment.
#'
#' @param web An interaction matrix.
#' @param change_A,change_B Numeric vectors of projected suitability change for
#'   each level. Species with the largest decline are removed first.
#' @param threshold Proportion of interaction strength that must be lost to
#'   trigger secondary extinction.
#' @param rewiring Optional rewiring probability matrix.
#' @param n_sim Replicates.
#' @param interaction An [si_interaction()].
#' @param seed Optional seed.
#' @return A `data.frame` giving, at each step, the number of primary
#'   extinctions and secondary coextinctions in each level.
#' @examples
#' set.seed(12)
#' mw <- si_simulate_metaweb(20, 16, 0.2, seed = 12)
#' w  <- mw$matrix * matrix(runif(320, 1, 5), 20, 16)
#' dA <- rnorm(20, -0.3, 0.2); dB <- rnorm(16, -0.2, 0.2)
#' cx <- si_coextinction(w, dA, dB, threshold = 0.5, n_sim = 20)
#' head(cx)
#' tail(cx, 1)
#' @references
#' Schleuning, M. and others (2016) Ecological networks are more sensitive to
#' plant than to animal extinction under climate change. *Nature
#' Communications*, 7, 13965.
#' @export
si_coextinction <- function(web, change_A, change_B, threshold = 1,
                            rewiring = NULL, n_sim = 50,
                            interaction = si_interaction("pollinatedBy"),
                            seed = NULL) {
  W0 <- as.matrix(web); W0[is.na(W0)] <- 0
  nA <- nrow(W0); nB <- ncol(W0)
  stopifnot(length(change_A) == nA, length(change_B) == nB)
  if (!is.null(seed)) set.seed(seed)
  rule <- si_cascade_rule(interaction)

  ordA <- order(change_A); ordB <- order(change_B)
  steps <- max(nA, nB)
  acc <- matrix(0, steps + 1L, 4,
                dimnames = list(NULL, c("primary_A", "primary_B",
                                        "coextinct_A", "coextinct_B")))
  for (s in seq_len(n_sim)) {
    W <- W0
    baseA <- rowSums(W0); baseB <- colSums(W0)
    aliveA <- rep(TRUE, nA); aliveB <- rep(TRUE, nB)
    for (t in seq_len(steps)) {
      pa <- ceiling(t * nA / steps); pb <- ceiling(t * nB / steps)
      killA <- ordA[seq_len(pa)]; killB <- ordB[seq_len(pb)]
      aliveA[killA] <- FALSE; aliveB[killB] <- FALSE
      W <- W0; W[!aliveA, ] <- 0; W[, !aliveB] <- 0
      if (!is.null(rewiring)) W <- .reallocate(W0, W, aliveA, aliveB, rewiring)
      if (rule[["A"]] == "extinct")
        aliveA <- aliveA & (rowSums(W) > (1 - threshold) * baseA) & rowSums(W) > 0
      if (rule[["B"]] == "extinct")
        aliveB <- aliveB & (colSums(W) > (1 - threshold) * baseB) & colSums(W) > 0
      acc[t + 1L, ] <- acc[t + 1L, ] +
        c(pa, pb, sum(!aliveA) - pa, sum(!aliveB) - pb)
    }
  }
  acc <- acc / n_sim
  data.frame(step = 0:steps, acc,
             surviving_A = nA - acc[, "primary_A"] - acc[, "coextinct_A"],
             surviving_B = nB - acc[, "primary_B"] - acc[, "coextinct_B"])
}

.reallocate <- function(W0, W, aliveA, aliveB, rewiring) {
  lostB <- colSums(W0) - colSums(W)
  for (j in which(aliveB & lostB > 0)) {
    cand <- which(aliveA & rewiring[, j] > 0)
    if (!length(cand)) next
    p <- rewiring[cand, j] / sum(rewiring[cand, j])
    W[cand, j] <- W[cand, j] + lostB[j] * p
  }
  W
}
