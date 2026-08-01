# ---------------------------------------------------------------------------
# Module 8: simulation engine
# ---------------------------------------------------------------------------

#' Simulate an interaction constraint layer with controlled architecture
#'
#' Generates metawebs of a target connectance and a target architecture. Because
#' architecture and connectance can be varied independently, this is the tool
#' that lets a virtual experiment separate the effect of *how many* links are
#' permitted from the effect of *which* links are permitted.
#'
#' @param n_A,n_B Numbers of species in the two levels.
#' @param connectance Target proportion of realised links.
#' @param architecture One of:
#'   \describe{
#'     \item{`"random"`}{Erdos-Renyi: links placed uniformly.}
#'     \item{`"nested"`}{generalists interact with subsets of each other's
#'       partners, produced by thresholding an outer product of ranked
#'       propensities.}
#'     \item{`"modular"`}{`n_modules` blocks, with `p_between` of the links
#'       leaking across block boundaries.}
#'     \item{`"heterogeneous"`}{power-law degree distribution with exponent
#'       `gamma`.}
#'   }
#' @param n_modules Integer, for `"modular"`.
#' @param p_between Numeric, proportion of links placed between modules.
#' @param gamma Numeric exponent, for `"heterogeneous"`.
#' @param names_A,names_B Optional species names.
#' @param min_degree Integer. Guarantee every species at least this many links,
#'   so that no species is unconnected in the metaweb.
#' @param seed Optional random seed.
#'
#' @return An [si_metaweb()].
#'
#' @examples
#' set.seed(1)
#' archs <- c("random", "nested", "modular", "heterogeneous")
#' t(sapply(archs, function(a) {
#'   mw <- si_simulate_metaweb(30, 25, connectance = 0.15,
#'                             architecture = a, seed = 1)
#'   c(connectance = si_connectance(mw),
#'     NODF = si_nodf(mw$matrix),
#'     modularity = si_modularity(mw$matrix))
#' }))
#' @seealso [si_metaweb()], [si_forbidden()]
#' @export
si_simulate_metaweb <- function(n_A, n_B, connectance = 0.2,
                                architecture = c("random", "nested",
                                                 "modular", "heterogeneous"),
                                n_modules = 4, p_between = 0.1, gamma = 2.5,
                                names_A = NULL, names_B = NULL,
                                min_degree = 1, seed = NULL) {
  architecture <- match.arg(architecture)
  if (!is.null(seed)) set.seed(seed)
  if (connectance <= 0 || connectance > 1)
    stop("`connectance` must be in (0, 1].", call. = FALSE)
  target <- max(1L, round(connectance * n_A * n_B))

  P <- switch(architecture,
    random = matrix(stats::runif(n_A * n_B), n_A, n_B),
    nested = {
      ra <- sort(stats::runif(n_A), decreasing = TRUE)
      rb <- sort(stats::runif(n_B), decreasing = TRUE)
      outer(ra, rb) + stats::runif(n_A * n_B, 0, 0.02)
    },
    modular = {
      ga <- rep(seq_len(n_modules), length.out = n_A)
      gb <- rep(seq_len(n_modules), length.out = n_B)
      within <- outer(ga, gb, "==")
      m <- matrix(stats::runif(n_A * n_B), n_A, n_B)
      m + ifelse(within, 1, p_between)
    },
    heterogeneous = {
      da <- (seq_len(n_A))^(-1 / (gamma - 1))
      db <- (seq_len(n_B))^(-1 / (gamma - 1))
      outer(sample(da), sample(db)) * matrix(stats::runif(n_A * n_B), n_A, n_B)
    })

  thr <- sort(as.vector(P), decreasing = TRUE)[target]
  M <- (P >= thr) * 1
  # guarantee a minimum degree so that no species is isolated in the metaweb
  if (min_degree > 0) {
    for (i in which(rowSums(M) < min_degree))
      M[i, sample.int(n_B, min_degree)] <- 1
    for (j in which(colSums(M) < min_degree))
      M[sample.int(n_A, min_degree), j] <- 1
  }
  if (is.null(names_A)) names_A <- paste0("A", seq_len(n_A))
  if (is.null(names_B)) names_B <- paste0("B", seq_len(n_B))
  dimnames(M) <- list(names_A, names_B)
  mw <- new_si_metaweb(M, FALSE)
  mw$architecture <- architecture
  mw
}


#' Simulate virtual species suitability surfaces
#'
#' Generates spatially autocorrelated habitat-suitability surfaces on a regular
#' grid. Each species is given a Gaussian response to two latent environmental
#' gradients, which are themselves spatially smooth. This follows the logic of
#' virtual-species simulation (Leroy and others, 2016) but is implemented
#' natively so that hundreds of thousands of surfaces can be generated without
#' an external dependency.
#'
#' @param n_species Number of species.
#' @param nrow,ncol Grid dimensions.
#' @param niche_breadth Numeric, or a length-2 vector giving the range from
#'   which each species' breadth is drawn. Small values give specialists with
#'   small ranges.
#' @param autocorrelation Numeric > 0. Spatial smoothing of the environmental
#'   gradients; larger values give coarser, more clumped environments.
#' @param prefix Character prefix for species names.
#' @param as_raster Logical. Return a `SpatRaster` rather than an
#'   [si_stack()].
#' @param seed Optional random seed.
#'
#' @return An [si_stack()], or a `SpatRaster` if `as_raster = TRUE`, with the
#'   niche parameters kept in the `niche` attribute so that a matching future
#'   projection can be produced by [si_simulate_climate()].
#'
#' @examples
#' set.seed(1)
#' st <- si_simulate_species(6, nrow = 30, ncol = 30, seed = 1)
#' st
#' head(summary(st), 3)
#'
#' # specialists occupy less of the grid than generalists
#' sp <- si_simulate_species(40, nrow = 25, ncol = 25,
#'                           niche_breadth = c(0.05, 0.5), seed = 2)
#' cor(attr(sp, "niche")$breadth, colMeans(sp$values), method = "spearman")
#' @references
#' Leroy, B., Meynard, C.N., Bellard, C. & Courchamp, F. (2016)
#' virtualspecies, an R package to generate virtual species distributions.
#' *Ecography*, 39, 599--607.
#' @export
si_simulate_species <- function(n_species, nrow = 40, ncol = 40,
                                niche_breadth = c(0.08, 0.35),
                                autocorrelation = 3, prefix = "sp",
                                as_raster = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  env <- lapply(1:2, function(i) .smooth_field(nrow, ncol, autocorrelation))
  if (length(niche_breadth) == 1L) niche_breadth <- rep(niche_breadth, 2)
  opt1 <- stats::runif(n_species); opt2 <- stats::runif(n_species)
  br <- stats::runif(n_species, niche_breadth[1], niche_breadth[2])
  v <- vapply(seq_len(n_species), function(i) {
    exp(-((env[[1]] - opt1[i])^2 + (env[[2]] - opt2[i])^2) / (2 * br[i]^2))
  }, numeric(nrow * ncol))
  colnames(v) <- paste0(prefix, seq_len(n_species))
  niche <- data.frame(species = colnames(v), opt1 = opt1, opt2 = opt2,
                      breadth = br)
  if (as_raster) {
    r <- terra::rast(nrows = nrow, ncols = ncol, nlyrs = n_species,
                     xmin = 0, xmax = ncol, ymin = 0, ymax = nrow)
    terra::values(r) <- v
    names(r) <- colnames(v)
    attr(r, "niche") <- niche
    attr(r, "env") <- env
    return(r)
  }
  st <- si_stack(v, level = prefix)
  attr(st, "niche") <- niche
  attr(st, "env") <- env
  attr(st, "grid_dim") <- c(nrow, ncol)
  st
}

.smooth_field <- function(nr, nc, sigma) {
  z <- matrix(stats::rnorm(nr * nc), nr, nc)
  k <- max(1L, round(sigma))
  # separable box smoothing repeated twice approximates a Gaussian
  sm <- function(x, k) {
    n <- nrow(x)
    cs <- rbind(0, apply(x, 2, cumsum))
    lo <- pmax(1, seq_len(n) - k); hi <- pmin(n, seq_len(n) + k)
    (cs[hi + 1, , drop = FALSE] - cs[lo, , drop = FALSE]) / (hi - lo + 1)
  }
  for (i in 1:2) { z <- sm(z, k); z <- t(sm(t(z), k)) }
  z <- as.vector(z)
  (z - min(z)) / diff(range(z))
}


#' Project a suitability stack under a climate-change scenario
#'
#' Applies one of several transformations to a current suitability surface to
#' produce a future one. Both non-spatial scenarios (which rescale suitability
#' in place) and spatially explicit ones (which displace the surface) are
#' provided, because the two make different predictions about interaction loss.
#'
#' @param stack An [si_stack()], ideally one produced by
#'   [si_simulate_species()] so that its `niche` and `grid_dim` attributes are
#'   available.
#' @param scenario One of:
#'   \describe{
#'     \item{`"decline"`}{all species lose suitability by `severity`.}
#'     \item{`"increase"`}{all species gain suitability.}
#'     \item{`"mixed"`}{half decline, half increase.}
#'     \item{`"random"`}{each species' change drawn from
#'       `N(-severity, severity)`.}
#'     \item{`"trait_dependent"`}{narrow-niche species decline most, following
#'       the empirical pattern that biotic specialists have narrow climatic
#'       niches.}
#'     \item{`"shift"`}{spatially explicit poleward displacement of the
#'       surface by `distance` cells.}
#'     \item{`"upslope"`}{spatially explicit displacement along an elevation
#'       gradient, implemented as displacement plus contraction.}
#'     \item{`"contraction"`}{range contraction toward each species'
#'       optimum, without displacement.}
#'   }
#' @param severity Numeric on 0--1. Proportional change in suitability for the
#'   non-spatial scenarios.
#' @param distance Numeric. Displacement in cells for `"shift"` and
#'   `"upslope"`.
#' @param dims Integer vector `c(nrow, ncol)`, required for the spatial
#'   scenarios if the stack carries no `grid_dim` attribute.
#' @param seed Optional random seed.
#'
#' @return An [si_stack()] with the same species and cells.
#'
#' @examples
#' st <- si_simulate_species(10, nrow = 30, ncol = 30, seed = 5)
#' scen <- c("decline", "trait_dependent", "shift", "contraction")
#' round(sapply(scen, function(s)
#'   mean(si_simulate_climate(st, s, severity = 0.3, distance = 5,
#'                            seed = 1)$values)), 4)
#' round(mean(st$values), 4)   # current, for comparison
#' @export
si_simulate_climate <- function(stack,
                                scenario = c("decline", "increase", "mixed",
                                             "random", "trait_dependent",
                                             "shift", "upslope",
                                             "contraction"),
                                severity = 0.3, distance = 5, dims = NULL,
                                seed = NULL) {
  scenario <- match.arg(scenario)
  if (!inherits(stack, "si_stack")) stop("`stack` must be an si_stack.",
                                         call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  v <- stack$values; n <- ncol(v)

  f <- switch(scenario,
    decline  = rep(-severity, n),
    increase = rep(severity, n),
    mixed    = rep(c(-severity, severity), length.out = n),
    random   = stats::rnorm(n, -severity, severity),
    trait_dependent = {
      br <- attr(stack, "niche")$breadth
      if (is.null(br)) br <- apply(v, 2, stats::sd)
      r <- rank(br) / length(br)                # 1 = broadest niche
      -severity * 2 * (1 - r)                   # narrow niches decline most
    },
    NULL)

  if (!is.null(f)) {
    out <- sweep(v, 2, 1 + f, "*")
    out[out < 0] <- 0; out[out > 1] <- 1
    stack$values <- out
    colnames(stack$values) <- stack$names
    return(stack)
  }

  d <- if (!is.null(dims)) dims else attr(stack, "grid_dim")
  if (is.null(d))
    stop("Spatial scenario '", scenario, "' needs grid dimensions.\n",
         "  Supply `dims = c(nrow, ncol)` or use a stack from ",
         "si_simulate_species().", call. = FALSE)
  nr <- d[1]; nc <- d[2]
  if (nr * nc != nrow(v))
    stop("`dims` (", nr, " x ", nc, ") does not match the ", nrow(v),
         " cells in the stack.", call. = FALSE)

  shift_rows <- function(mat, k) {
    a <- array(mat, c(nr, nc, ncol(mat)))
    k <- round(k)
    if (k == 0) return(mat)
    out <- array(0, dim(a))
    if (k > 0) out[(k + 1):nr, , ] <- a[1:(nr - k), , , drop = FALSE]
    else out[1:(nr + k), , ] <- a[(1 - k):nr, , , drop = FALSE]
    matrix(out, nr * nc, ncol(mat))
  }

  out <- switch(scenario,
    shift = shift_rows(v, distance),
    upslope = {
      s <- shift_rows(v, distance)
      sweep(s, 2, 1 - severity, "*")
    },
    contraction = {
      m <- apply(v, 2, max)
      s <- sweep(v, 2, pmax(m, 1e-9), "/")
      sweep(s^(1 / max(1e-6, 1 - severity)), 2, m, "*")
    })
  out[out < 0] <- 0; out[out > 1] <- 1
  stack$values <- out
  colnames(stack$values) <- stack$names
  stack
}


#' Simulate species phenologies and their shift under climate change
#'
#' Produces a Gaussian activity phenology for each species and, optionally, a
#' shifted future phenology. The among-species standard deviation of the shift
#' (`shift_sd`) is the key parameter: with `shift_sd = 0` every species
#' advances by the same amount, phenological overlap is unchanged, and no
#' rewiring is possible. Rewiring is generated by *variance* in the shift, not
#' by its mean.
#'
#' @param n_A,n_B Numbers of species in the two levels.
#' @param season Numeric length of the annual cycle (365 for days, 12 for
#'   months).
#' @param duration Length-2 numeric giving the range of activity durations
#'   (standard deviations) drawn per species.
#' @param shift_mean,shift_sd Mean and among-species standard deviation of the
#'   phenological shift, in the same units as `season`. Negative `shift_mean`
#'   is an advance.
#' @param seed Optional random seed.
#'
#' @return A list of class `si_phenology` with `current` and `future`
#'   components, each a list of `mu` and `sigma` for both levels.
#'
#' @examples
#' ph <- si_phenology_simulate(20, 15, shift_mean = -10, shift_sd = 15, seed = 1)
#' str(ph$current$mu_A[1:5])
#'
#' # the overlap matrices differ only when shift_sd > 0
#' o_now <- si_phenology_overlap(ph, "current")
#' o_fut <- si_phenology_overlap(ph, "future")
#' mean(abs(o_now - o_fut))
#' @seealso [si_rule_phenology()], [si_phenology_from_records()]
#' @export
si_phenology_simulate <- function(n_A, n_B, season = 365,
                                  duration = c(20, 60),
                                  shift_mean = -10, shift_sd = 0,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  mk <- function(n) list(mu = stats::runif(n, 0, season),
                         sigma = stats::runif(n, duration[1], duration[2]))
  a <- mk(n_A); b <- mk(n_B)
  cur <- list(mu_A = a$mu, sd_A = a$sigma, mu_B = b$mu, sd_B = b$sigma)
  fut <- cur
  fut$mu_A <- (cur$mu_A + stats::rnorm(n_A, shift_mean, shift_sd)) %% season
  fut$mu_B <- (cur$mu_B + stats::rnorm(n_B, shift_mean, shift_sd)) %% season
  structure(list(current = cur, future = fut, season = season,
                 shift_mean = shift_mean, shift_sd = shift_sd),
            class = "si_phenology")
}

#' Phenological overlap matrix
#'
#' @param phenology An [si_phenology_simulate()] object, or a list with `mu_A`,
#'   `sd_A`, `mu_B`, `sd_B`.
#' @param slice `"current"` or `"future"`.
#' @return A matrix of overlap on 0--1.
#' @examples
#' ph <- si_phenology_simulate(8, 6, seed = 2)
#' round(si_phenology_overlap(ph)[1:4, 1:4], 3)
#' @export
si_phenology_overlap <- function(phenology, slice = c("current", "future")) {
  slice <- match.arg(slice)
  p <- if (inherits(phenology, "si_phenology")) phenology[[slice]] else phenology
  season <- if (inherits(phenology, "si_phenology")) phenology$season else 365
  d <- abs(outer(p$mu_A, p$mu_B, "-"))
  d <- pmin(d, season - d)                      # circular distance
  s <- sqrt(outer(p$sd_A^2, p$sd_B^2, "+"))
  exp(-0.5 * (d / s)^2)
}


#' Derive species phenology from dated interaction records
#'
#' Extracts a month-by-species activity matrix from an interaction table that
#' carries an ISO 8601 `eventDate` field, as used by Darwin Core and Global
#' Biotic Interactions. Records whose date span covers a full year or more are
#' uninformative and are discarded; the remainder contribute the months in
#' their window to both partners.
#'
#' The result is an upper bound on activity, because the recorded window is the
#' study's sampling period rather than the species' phenology. It is
#' nonetheless internally checkable: an observed interaction implies that both
#' partners were active simultaneously, so a high proportion of observed links
#' with zero month overlap would indicate that the inference has failed.
#' [si_phenology_check()] reports that proportion.
#'
#' @param records A `data.frame` of interaction records.
#' @param species_A,species_B Column names holding the two species names.
#' @param date Column name holding the date or date range.
#' @param names_A,names_B Optional species pools to align the output to.
#' @param max_span Integer. Discard windows spanning at least this many months.
#' @return A list of class `si_phenology_records` with month-by-species
#'   matrices `A` and `B` and a `coverage` summary.
#' @examples
#' recs <- data.frame(
#'   plant  = c("Sp_a", "Sp_a", "Sp_b", "Sp_c"),
#'   animal = c("An_x", "An_y", "An_y", "An_z"),
#'   eventDate = c("2014-09/2015-03", "2012-11", "2015-08/2015-11",
#'                 "1991-01/1997-12"))
#' ph <- si_phenology_from_records(recs, "plant", "animal", "eventDate")
#' ph$coverage
#' rowSums(t(ph$A))          # months of record per plant
#' @seealso [si_phenology_simulate()], [si_rule_phenology()]
#' @export
si_phenology_from_records <- function(records, species_A, species_B,
                                      date = "eventDate",
                                      names_A = NULL, names_B = NULL,
                                      max_span = 12L) {
  for (v in c(species_A, species_B, date))
    if (!v %in% names(records))
      stop("Column '", v, "' not found in `records`.", call. = FALSE)
  mo <- lapply(as.character(records[[date]]), .months_from_span, max_span = max_span)
  keep <- !vapply(mo, is.null, logical(1))
  a <- .norm_names(records[[species_A]])[keep]
  b <- .norm_names(records[[species_B]])[keep]
  mo <- mo[keep]
  if (is.null(names_A)) names_A <- sort(unique(a[!is.na(a) & a != "NA"]))
  if (is.null(names_B)) names_B <- sort(unique(b[!is.na(b) & b != "NA"]))
  A <- matrix(0L, 12L, length(names_A), dimnames = list(month.abb, names_A))
  B <- matrix(0L, 12L, length(names_B), dimnames = list(month.abb, names_B))
  for (i in seq_along(mo)) {
    ia <- match(a[i], names_A); ib <- match(b[i], names_B)
    if (!is.na(ia)) A[mo[[i]], ia] <- 1L
    if (!is.na(ib)) B[mo[[i]], ib] <- 1L
  }
  cov <- data.frame(
    level = c("A", "B"),
    n_species = c(length(names_A), length(names_B)),
    n_dated = c(sum(colSums(A) > 0), sum(colSums(B) > 0)),
    median_months = c(stats::median(colSums(A)[colSums(A) > 0]),
                      stats::median(colSums(B)[colSums(B) > 0])))
  structure(list(A = A, B = B, coverage = cov, n_records = sum(keep)),
            class = "si_phenology_records")
}

.months_from_span <- function(s, max_span = 12L) {
  s <- trimws(s)
  if (is.na(s) || !nzchar(s)) return(NULL)
  rng <- regmatches(s, regexec(
    "^(\\d{4})-(\\d{2})(?:-\\d{2})?/(\\d{4})-(\\d{2})(?:-\\d{2})?$", s))[[1]]
  if (length(rng) == 5L) {
    y1 <- as.integer(rng[2]); m1 <- as.integer(rng[3])
    y2 <- as.integer(rng[4]); m2 <- as.integer(rng[5])
    span <- (y2 - y1) * 12L + (m2 - m1) + 1L
    if (span >= max_span || span < 1L) return(NULL)
    return(((m1 - 1L + seq_len(span) - 1L) %% 12L) + 1L)
  }
  one <- regmatches(s, regexec("^(\\d{4})-(\\d{2})$", s))[[1]]
  if (length(one) == 3L) return(as.integer(one[3]))
  NULL
}

#' @export
print.si_phenology_records <- function(x, ...) {
  cat("<si_phenology_records>  informative records: ", x$n_records, "\n", sep = "")
  print(x$coverage)
  invisible(x)
}

#' Internal consistency check for a derived phenology
#'
#' An observed interaction requires that both partners were active at the same
#' time. This function reports the proportion of observed links whose partners
#' have no month in common, which should be close to zero if the derived
#' phenology is sound.
#'
#' @param phenology An [si_phenology_from_records()] result.
#' @param metaweb The observed interaction matrix (an [si_metaweb()] or matrix)
#'   aligned to the same species.
#' @return A list with the number of testable links, the number of violations
#'   and the violation rate.
#' @examples
#' recs <- data.frame(
#'   plant  = c("Sp_a", "Sp_a", "Sp_b"),
#'   animal = c("An_x", "An_y", "An_y"),
#'   eventDate = c("2014-09/2015-03", "2012-11", "2015-08/2015-11"))
#' ph <- si_phenology_from_records(recs, "plant", "animal", "eventDate")
#' mw <- si_metaweb(recs[, 1:2], names_A = colnames(ph$A),
#'                  names_B = colnames(ph$B))
#' si_phenology_check(ph, mw)
#' @export
si_phenology_check <- function(phenology, metaweb) {
  M <- if (inherits(metaweb, "si_metaweb")) metaweb$matrix else as.matrix(metaweb)
  A <- phenology$A; B <- phenology$B
  ij <- which(M > 0, arr.ind = TRUE)
  dated <- colSums(A) > 0; datedB <- colSums(B) > 0
  ok <- dated[ij[, 1]] & datedB[ij[, 2]]
  ij <- ij[ok, , drop = FALSE]
  if (!nrow(ij))
    return(list(testable = 0L, violations = 0L, violation_rate = NA_real_))
  viol <- sum(vapply(seq_len(nrow(ij)), function(k)
    sum(A[, ij[k, 1]] & B[, ij[k, 2]]) == 0, logical(1)))
  list(testable = nrow(ij), violations = viol,
       violation_rate = viol / nrow(ij))
}
