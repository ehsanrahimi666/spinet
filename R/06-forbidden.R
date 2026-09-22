# ---------------------------------------------------------------------------
# Module 3: forbidden links
# ---------------------------------------------------------------------------

#' Assemble a constraint layer from forbidden-link rules
#'
#' Combines any number of rule matrices into the constraint layer used by
#' [si_overlap()]. Each rule is a matrix on 0--1 of the same dimensions, in
#' which zero means the link is impossible and one means it is unconstrained;
#' intermediate values express graded plausibility.
#'
#' @param ... Rule matrices, or objects returned by the `si_rule_*` functions.
#'   Named arguments are retained for attribution by [si_forbidden_partition()].
#' @param observed Optional record of documented interactions (an
#'   [si_metaweb()] or a matrix), aligned to the same species as the rules.
#' @param observed_action How `observed` is used. `"override"` (the default)
#'   applies it *after* the rules: a documented pair is permitted whatever the
#'   rules say, and undocumented pairs are left to the rules. This is the right
#'   choice when the record is trusted evidence and the rules are inference,
#'   because a trait rule that contradicts an observation is wrong about that
#'   pair. `"filter"` instead treats the record as one more rule, so a pair
#'   that was never recorded is forbidden however plausible the traits make it;
#'   use this only when the record is close to complete, since absence of
#'   evidence is otherwise read as evidence of absence.
#' @param combine How to combine rules: `"product"` (the Hadamard product used
#'   by Vazquez and others, 2009, and the default), `"min"` (the most
#'   restrictive rule wins) or `"mean"`.
#' @param binarise Logical. Threshold the combined layer at `threshold`.
#' @param threshold Numeric cut used when `binarise = TRUE`.
#'
#' @return An object of class `si_forbidden`, carrying the combined `metaweb`,
#'   the individual `rules`, per-rule statistics, and — when `observed` is
#'   supplied — the number of documented links and, under `"override"`, the
#'   number of them that the rules would otherwise have forbidden
#'   (`rescued_links`).
#'
#' @details
#' Rules multiply, so a link survives only if every rule permits it. This is
#' the standard treatment of multiple interaction determinants, in which
#' abundance, phenological overlap, spatial overlap and trait matching each
#' contribute a probability matrix that is combined element-wise.
#'
#' Documented interactions are handled separately from rules, because they are
#' evidence rather than inference. Under the default `observed_action =
#' "override"` they are applied after the rules have been combined and, if
#' requested, binarised, so a recorded interaction is never removed by a trait
#' or phenology rule that disagrees with it. `si_forbidden()` reports how many
#' links this rescued, which is a useful diagnostic: a large number means the
#' rules and the record disagree badly and the rules should be revisited.
#'
#' @examples
#' set.seed(1)
#' nA <- 20; nB <- 15
#' ph  <- si_phenology_simulate(nA, nB, seed = 1)
#' tr  <- data.frame(A = runif(nA, 5, 40))        # corolla depth
#' tb  <- data.frame(B = runif(nB, 3, 35))        # proboscis length
#'
#' fl <- si_forbidden(
#'   phenology  = si_rule_phenology(ph, threshold = 0.3),
#'   morphology = si_rule_morphology(tr$A, tb$B, rule = "barrier"),
#'   combine    = "product", binarise = TRUE)
#' fl
#' si_forbidden_partition(fl)
#'
#' # documented interactions, some of which the trait rule would forbid
#' obs <- matrix(0, nA, nB)
#' obs[cbind(sample(nA, 6), sample(nB, 6))] <- 1
#'
#' # default: the record wins where it disagrees with the rules
#' ov <- si_forbidden(morphology = si_rule_morphology(tr$A, tb$B), observed = obs)
#' c(rules_only = sum(si_rule_morphology(tr$A, tb$B)$matrix > 0),
#'   with_record = sum(ov$matrix > 0), rescued = ov$rescued_links)
#'
#' # filter: only documented pairs survive
#' ft <- si_forbidden(morphology = si_rule_morphology(tr$A, tb$B), observed = obs,
#'                    observed_action = "filter")
#' sum(ft$matrix > 0)
#' @references
#' Vazquez, D.P., Blüthgen, N., Cagnolo, L. & Chacoff, N.P. (2009) Uniting
#' pattern and process in plant-animal mutualistic networks: a review.
#' *Annals of Botany*, 103, 1445--1457.
#' @seealso [si_rule_phenology()], [si_rule_morphology()],
#'   [si_rule_elevation()], [si_rule_taxonomy()], [si_rule_custom()]
#' @export
si_forbidden <- function(..., observed = NULL,
                         observed_action = c("override", "filter"),
                         combine = c("product", "min", "mean"),
                         binarise = FALSE, threshold = 0.5) {
  combine <- match.arg(combine)
  observed_action <- match.arg(observed_action)
  rules <- list(...)
  obs <- NULL
  if (!is.null(observed)) {
    om <- .constraint_matrix(observed)
    obs <- (om > 0) * 1
    # "filter" makes the record a rule like any other, so an undocumented pair
    # is forbidden. "override" applies it after the rules, so a documented pair
    # is permitted whatever the rules say, and undocumented pairs are left to
    # the rules.
    if (observed_action == "filter") rules$observed <- obs
  }
  if (!length(rules))
    stop("Supply at least one rule matrix (see ?si_rule_phenology).",
         call. = FALSE)
  rules <- lapply(rules, function(r)
    if (inherits(r, "si_rule")) r$matrix else as.matrix(r))
  d <- vapply(rules, dim, integer(2))
  if (any(d[1, ] != d[1, 1]) || any(d[2, ] != d[2, 1]))
    stop("All rules must have the same dimensions. Got: ",
         paste(apply(d, 2, paste, collapse = "x"), collapse = ", "),
         call. = FALSE)
  if (is.null(names(rules)) || any(names(rules) == ""))
    names(rules) <- paste0("rule", seq_along(rules))

  # Rules built by different functions need not list species in the same
  # order. When every rule carries species names, align them to the first
  # rule by name rather than by position; refuse to combine rules whose
  # species sets differ.
  named <- vapply(rules, function(r) !is.null(rownames(r)) && !is.null(colnames(r)),
                  logical(1))
  if (all(named) && length(rules) > 1) {
    rn <- .norm_names(rownames(rules[[1]])); cn <- .norm_names(colnames(rules[[1]]))
    for (k in seq_along(rules)[-1]) {
      ri <- match(rn, .norm_names(rownames(rules[[k]])))
      ci <- match(cn, .norm_names(colnames(rules[[k]])))
      if (anyNA(ri) || anyNA(ci))
        stop("Rule '", names(rules)[k], "' does not contain the same species as ",
             "rule '", names(rules)[1], "'.", call. = FALSE)
      rules[[k]] <- rules[[k]][ri, ci, drop = FALSE]
    }
    if (!is.null(obs) && !is.null(rownames(obs)) && !is.null(colnames(obs))) {
      ri <- match(rn, .norm_names(rownames(obs))); ci <- match(cn, .norm_names(colnames(obs)))
      if (!anyNA(ri) && !anyNA(ci)) obs <- obs[ri, ci, drop = FALSE]
      if (observed_action == "filter") rules$observed <- obs
    }
  } else if (length(rules) > 1 && any(named) && !all(named)) {
    warning("Some rules carry species names and others do not; rules were ",
            "combined by position.", call. = FALSE)
  }

  arr <- simplify2array(rules)
  M <- switch(combine,
    product = apply(arr, 1:2, prod),
    min     = apply(arr, 1:2, min),
    mean    = apply(arr, 1:2, mean))
  if (binarise) M <- (M >= threshold) * 1
  dimnames(M) <- dimnames(rules[[1]])

  rescued <- 0L
  if (!is.null(obs) && observed_action == "override") {
    if (!identical(dim(obs), dim(M)))
      stop("`observed` is ", nrow(obs), " x ", ncol(obs), " but the rules are ",
           nrow(M), " x ", ncol(M), ".", call. = FALSE)
    rescued <- sum(obs > 0 & M == 0)
    M <- pmax(M, obs)
    dimnames(M) <- dimnames(rules[[1]])
  }

  stats <- data.frame(
    rule = names(rules),
    permitted = vapply(rules, function(r) mean(r > 0), numeric(1)),
    forbidden = vapply(rules, function(r) mean(r == 0), numeric(1)),
    row.names = NULL)

  structure(list(metaweb = new_si_metaweb(M, all(M == 1)),
                 matrix = M, rules = rules, stats = stats,
                 combine = combine, binarised = binarise,
                 observed_action = if (is.null(obs)) NA_character_ else observed_action,
                 observed_links = if (is.null(obs)) 0L else sum(obs > 0),
                 rescued_links = rescued),
            class = "si_forbidden")
}

#' @export
print.si_forbidden <- function(x, ...) {
  cat("<si_forbidden>  combine = ", x$combine,
      if (x$binarised) "  (binarised)" else "", "\n", sep = "")
  cat("  dimensions : ", nrow(x$matrix), " x ", ncol(x$matrix), "\n", sep = "")
  cat("  permitted  : ", sum(x$matrix > 0), " links  (connectance ",
      format(round(mean(x$matrix > 0), 4), nsmall = 4), ")\n", sep = "")
  cat("  rules      :\n")
  s <- x$stats
  for (i in seq_len(nrow(s)))
    cat(sprintf("     %-14s forbids %5.1f%% of all pairs\n",
                s$rule[i], 100 * s$forbidden[i]))
  if (!is.na(x$observed_action)) {
    cat("  observed   : ", x$observed_links, " documented links, applied as \"",
        x$observed_action, "\"\n", sep = "")
    if (identical(x$observed_action, "override"))
      cat("               ", x$rescued_links,
          " of them would have been forbidden by the rules\n", sep = "")
  }
  invisible(x)
}

new_si_rule <- function(m, type, detail = list()) {
  structure(list(matrix = m, type = type, detail = detail), class = "si_rule")
}

#' @export
print.si_rule <- function(x, ...) {
  cat("<si_rule: ", x$type, ">  ", nrow(x$matrix), " x ", ncol(x$matrix),
      "   forbids ", sprintf("%.1f%%", 100 * mean(x$matrix == 0)),
      " of pairs\n", sep = "")
  invisible(x)
}

#' Attribute forbidden links to individual rules
#'
#' Reports how many links each rule removes on its own, how many it removes
#' uniquely (that no other rule would have removed), and the marginal effect of
#' dropping it. Useful for reporting which constraint actually drives a result.
#'
#' @param forbidden An [si_forbidden()] object.
#' @return A `data.frame`.
#' @examples
#' set.seed(2)
#' fl <- si_forbidden(
#'   r1 = matrix(rbinom(300, 1, 0.7), 20, 15),
#'   r2 = matrix(rbinom(300, 1, 0.6), 20, 15))
#' si_forbidden_partition(fl)
#' @export
si_forbidden_partition <- function(forbidden) {
  stopifnot(inherits(forbidden, "si_forbidden"))
  R <- forbidden$rules; nm <- names(R)
  zero <- lapply(R, function(r) r == 0)
  anyz <- Reduce(`|`, zero)
  out <- data.frame(
    rule = nm,
    forbids = vapply(zero, sum, numeric(1)),
    forbids_pct = 100 * vapply(zero, mean, numeric(1)),
    unique_to_rule = vapply(seq_along(zero), function(i) {
      others <- if (length(zero) > 1) Reduce(`|`, zero[-i]) else
        matrix(FALSE, nrow(zero[[i]]), ncol(zero[[i]]))
      sum(zero[[i]] & !others)
    }, numeric(1)),
    row.names = NULL)
  out$pct_of_all_forbidden <- 100 * out$forbids / max(sum(anyz), 1)
  out
}


# --- individual rules -------------------------------------------------------

#' Forbidden links from phenological non-overlap
#'
#' Two species cannot interact if they are never active at the same time. This
#' is the rule whose time-dependence generates rewiring: if the constraint
#' layer is recomputed from a shifted future phenology, shared species can
#' change partners, and interaction rewiring becomes possible.
#'
#' @param phenology An [si_phenology_simulate()] object, a list with `mu_A`,
#'   `sd_A`, `mu_B`, `sd_B`, or a pair of observed activity matrices supplied as
#'   `list(A = , B = )`. Activity matrices must have **time periods in rows and
#'   species in columns** (the layout returned by
#'   [si_phenology_from_records()]), and both must have the same number of
#'   rows. The returned rule is always level-A by level-B, matching every other
#'   `si_rule_*` function.
#' @param slice `"current"` or `"future"`.
#' @param threshold Numeric on 0--1. Overlap at or below this is forbidden.
#' @param graded Logical. Return the overlap itself as a graded rule rather
#'   than a hard 0/1 mask.
#' @return An `si_rule`.
#' @examples
#' ph <- si_phenology_simulate(15, 12, shift_mean = -10, shift_sd = 20, seed = 3)
#' now <- si_rule_phenology(ph, "current", threshold = 0.3)
#' fut <- si_rule_phenology(ph, "future",  threshold = 0.3)
#' c(now = mean(now$matrix > 0), future = mean(fut$matrix > 0))
#' # links that open and close
#' c(opened = sum(fut$matrix > 0 & now$matrix == 0),
#'   closed = sum(fut$matrix == 0 & now$matrix > 0))
#'
#' # observed activity matrices: months in rows, species in columns
#' A <- matrix(0L, 12, 5, dimnames = list(month.abb, paste0("plant", 1:5)))
#' B <- matrix(0L, 12, 3, dimnames = list(month.abb, paste0("bird", 1:3)))
#' A[1:6, ] <- 1L; B[4:9, ] <- 1L          # every pair shares Apr-Jun
#' obs <- si_rule_phenology(list(A = A, B = B), threshold = 0.1)
#' dim(obs$matrix)                          # 5 x 3, level A by level B
#' all(obs$matrix == 1)                     # TRUE: every pair overlaps
#' @export
si_rule_phenology <- function(phenology, slice = c("current", "future"),
                              threshold = 0.1, graded = FALSE) {
  slice <- match.arg(slice)
  O <- if (is.list(phenology) && !is.null(phenology$A) && is.matrix(phenology$A)) {
    a <- phenology$A > 0; b <- phenology$B > 0
    if (nrow(a) != nrow(b))
      stop("`phenology$A` and `phenology$B` must have the same number of rows ",
           "(time periods). Got ", nrow(a), " and ", nrow(b), ".\n",
           "  Both must be period-by-species matrices, periods in rows.",
           call. = FALSE)
    # Szymkiewicz-Simpson overlap: shared periods / periods of the less active
    # partner. crossprod(a, b) is already level-A by level-B; do not transpose.
    ov <- crossprod(a, b) / pmax(1, outer(colSums(a), colSums(b), pmin))
    dimnames(ov) <- list(colnames(phenology$A), colnames(phenology$B))
    ov
  } else {
    si_phenology_overlap(phenology, slice)
  }
  m <- if (graded) O else (O > threshold) * 1
  new_si_rule(m, "phenology", list(slice = slice, threshold = threshold))
}

#' Forbidden links from morphological mismatch
#'
#' @param trait_A,trait_B Numeric trait vectors, e.g. effective corolla depth
#'   and proboscis or bill length.
#' @param rule One of:
#'   \describe{
#'     \item{`"barrier"`}{a hard constraint: the consumer cannot exploit a
#'       resource deeper than its own reach (`trait_B >= trait_A - tolerance`).}
#'     \item{`"gower"`}{graded matching by the modified Gower similarity used
#'       for rewiring probabilities: closer traits are more likely partners.}
#'     \item{`"ratio"`}{a body-size ratio window `[lower, upper]`, appropriate
#'       for predator-prey and host-parasitoid systems.}
#'   }
#' @param tolerance Numeric slack added to the barrier rule.
#' @param lower,upper Ratio limits for `rule = "ratio"`.
#' @return An `si_rule`.
#' @examples
#' set.seed(4)
#' corolla   <- runif(12, 5, 40)
#' proboscis <- runif(10, 3, 35)
#' b <- si_rule_morphology(corolla, proboscis, rule = "barrier")
#' g <- si_rule_morphology(corolla, proboscis, rule = "gower")
#' c(barrier_permits = mean(b$matrix > 0), gower_mean = mean(g$matrix))
#'
#' # predator-prey mass ratio window
#' si_rule_morphology(runif(8, 1, 100), runif(6, 10, 1000),
#'                    rule = "ratio", lower = 2, upper = 100)
#' @references
#' Vizentin-Bugoni, J., Maruyama, P.K. & Sazima, M. (2014) Processes
#' entangling interactions in communities: forbidden links are more important
#' than abundance in a hummingbird-plant network. *Proceedings of the Royal
#' Society B*, 281, 20132397.
#' @export
si_rule_morphology <- function(trait_A, trait_B,
                               rule = c("barrier", "gower", "ratio"),
                               tolerance = 0, lower = 1, upper = Inf) {
  rule <- match.arg(rule)
  m <- switch(rule,
    barrier = outer(trait_A, trait_B, function(a, b) as.numeric(b + tolerance >= a)),
    gower = {
      rng <- diff(range(c(trait_A, trait_B)))
      if (rng == 0) rng <- 1
      1 - abs(outer(trait_A, trait_B, "-")) / rng
    },
    ratio = outer(trait_A, trait_B, function(a, b) {
      r <- b / a
      as.numeric(r >= lower & r <= upper)
    }))
  new_si_rule(m, paste0("morphology:", rule),
              list(rule = rule, tolerance = tolerance))
}

#' Forbidden links from non-overlapping elevational ranges
#'
#' @param lo_A,hi_A,lo_B,hi_B Numeric vectors of elevational limits.
#' @param min_overlap Numeric. Minimum overlap in metres required.
#' @param na_permit Logical. Treat species with missing limits as unconstrained
#'   (`TRUE`, the default) rather than forbidden.
#' @return An `si_rule`.
#' @examples
#' set.seed(5)
#' loA <- runif(10, 0, 2000); hiA <- loA + runif(10, 200, 1500)
#' loB <- runif(8, 0, 2000);  hiB <- loB + runif(8, 200, 1500)
#' si_rule_elevation(loA, hiA, loB, hiB)
#' @export
si_rule_elevation <- function(lo_A, hi_A, lo_B, hi_B, min_overlap = 0,
                              na_permit = TRUE) {
  # limits supplied in the wrong order are swapped rather than silently
  # producing an empty range
  swapA <- !is.na(lo_A) & !is.na(hi_A) & lo_A > hi_A
  swapB <- !is.na(lo_B) & !is.na(hi_B) & lo_B > hi_B
  if (any(swapA) || any(swapB)) {
    warning(sum(swapA) + sum(swapB), " elevational range(s) had the lower limit ",
            "above the upper limit and were swapped.", call. = FALSE)
    tmp <- lo_A[swapA]; lo_A[swapA] <- hi_A[swapA]; hi_A[swapA] <- tmp
    tmp <- lo_B[swapB]; lo_B[swapB] <- hi_B[swapB]; hi_B[swapB] <- tmp
  }
  ov <- outer(hi_A, lo_B, "-")
  ov2 <- outer(lo_A, hi_B, "-")
  m <- (ov >= min_overlap & -ov2 >= min_overlap) * 1
  if (na_permit) {
    miss <- outer(is.na(lo_A) | is.na(hi_A), is.na(lo_B) | is.na(hi_B), "|")
    m[miss] <- 1
  }
  m[is.na(m)] <- if (na_permit) 1 else 0
  new_si_rule(m, "elevation", list(min_overlap = min_overlap))
}

#' Forbidden links from taxonomic or functional specificity
#'
#' Encodes host specificity: a consumer restricted to one plant family cannot
#' use resources outside it. Also used for functional-group aggregation.
#'
#' @param group_A Character or factor vector giving each level-A species'
#'   group (family, functional group).
#' @param allowed A named list mapping each level-B species to the groups it
#'   can use, or a data frame with columns `species` and `group`.
#' @param names_B Character vector of level-B species names, in the order they
#'   appear in the stack.
#' @param default Numeric 0 or 1. What to do for level-B species with no entry
#'   in `allowed`. Default `1` (unconstrained).
#' @return An `si_rule`.
#' @examples
#' gA <- rep(c("Asteraceae", "Fabaceae", "Myrtaceae"), each = 4)
#' allowed <- list(bee1 = "Asteraceae", bee2 = c("Fabaceae", "Myrtaceae"))
#' si_rule_taxonomy(gA, allowed, names_B = c("bee1", "bee2", "fly1"))
#' @export
si_rule_taxonomy <- function(group_A, allowed, names_B, default = 1) {
  if (is.data.frame(allowed))
    allowed <- split(as.character(allowed$group), as.character(allowed$species))
  m <- matrix(default, length(group_A), length(names_B),
              dimnames = list(NULL, names_B))
  for (j in seq_along(names_B)) {
    a <- allowed[[names_B[j]]]
    if (!is.null(a)) m[, j] <- as.numeric(as.character(group_A) %in% a)
  }
  new_si_rule(m, "taxonomy", list(n_specified = sum(names_B %in% names(allowed))))
}

#' A user-defined forbidden-link rule
#'
#' The extension point of the forbidden-link framework. Any function of the two
#' species pools that returns a matrix on 0--1 can be used as a rule.
#'
#' @param f A function taking `(A, B)` and returning a matrix, or a matrix.
#' @param A,B Arguments passed to `f`.
#' @param name Character label used in [si_forbidden_partition()].
#' @return An `si_rule`.
#' @examples
#' # nocturnal pollinators cannot visit diurnal flowers
#' diurnal_flower <- c(TRUE, TRUE, FALSE, TRUE)
#' nocturnal_bee  <- c(FALSE, TRUE, FALSE)
#' si_rule_custom(function(a, b) outer(a, b, function(x, y) as.numeric(!(x & y))),
#'                diurnal_flower, nocturnal_bee, name = "activity_period")
#' @export
si_rule_custom <- function(f, A = NULL, B = NULL, name = "custom") {
  m <- if (is.function(f)) f(A, B) else as.matrix(f)
  if (!is.matrix(m)) stop("`f` must return a matrix.", call. = FALSE)
  if (any(m < 0 | m > 1, na.rm = TRUE))
    stop("Rule values must lie on [0, 1]. Got range [",
         paste(round(range(m, na.rm = TRUE), 3), collapse = ", "), "].",
         call. = FALSE)
  new_si_rule(m, name)
}

#' Neutral (abundance-based) interaction propensity
#'
#' Not a forbidden-link rule but its neutral counterpart: the probability that
#' two species meet, given their relative abundances. Included so that
#' niche-based and neutral determinants can be compared within one framework.
#'
#' @param abund_A,abund_B Numeric abundance vectors.
#' @param normalise Logical. Rescale so the maximum is one.
#' @return An `si_rule`.
#' @examples
#' si_rule_abundance(c(10, 1, 5), c(3, 20))
#' @export
si_rule_abundance <- function(abund_A, abund_B, normalise = TRUE) {
  m <- outer(abund_A / sum(abund_A), abund_B / sum(abund_B))
  if (normalise) m <- m / max(m)
  new_si_rule(m, "abundance")
}
