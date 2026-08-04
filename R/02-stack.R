# ---------------------------------------------------------------------------
# Module 1b: suitability stacks and constraint layers
# ---------------------------------------------------------------------------

#' Read a set of per-species suitability surfaces
#'
#' Builds the validated, memory-efficient representation of one level of a
#' network (for example, all modelled plants). `spinet` never stores the full
#' raster stack: it extracts a cells-by-species matrix once and keeps a
#' template raster for writing results back to disk.
#'
#' @param x One of: a directory path containing one raster per species; a
#'   character vector of raster file paths; a `SpatRaster` with one layer per
#'   species; or a numeric matrix with cells in rows and species in columns.
#' @param band Integer or character. Which band of each file to read when there
#'   is one file per species. Many SDM workflows write a continuous suitability
#'   band and a thresholded binary band into the same GeoTIFF; use `band = 1`
#'   (or `"max"`) for the continuous surface and `band = 2` (or
#'   `"sensitivity"`) for the binary one. Ignored when `x` is a matrix, or when
#'   `x` is a single multi-layer raster whose layers are the species.
#' @param mask Optional. A `SpatVector`, `SpatRaster` or logical vector used to
#'   restrict the analysis to a study area. Cells outside the mask are dropped
#'   entirely rather than set to `NA`, which is what makes continental
#'   analyses tractable.
#' @param names Optional character vector of species names. Defaults to the
#'   file names without extension, or the layer names.
#' @param level Character label for this level, used in messages.
#' @param na_value Numeric. Value substituted for `NA` cells. The default `0`
#'   treats "no prediction" as "not suitable", which is correct for masked SDM
#'   output.
#' @param pattern Regular expression used to select files when `x` is a
#'   directory. Default matches GeoTIFFs.
#' @param grid Optional [si_grid()] (or `SpatRaster`) describing the analysis
#'   grid. Only needed for matrix input that should still be mappable, for
#'   example the bundled [chile_pp] data.
#' @param cells Optional integer vector of raster cell indices matching the
#'   rows of a matrix input. Supply together with `grid`.
#' @param quiet Logical. Suppress progress messages.
#'
#' @return An object of class `si_stack`: a list with the cells-by-species
#'   matrix `values`, species `names`, the retained `cells` indices, and a
#'   zero-layer `template` (a `SpatRaster`) or `NULL` for matrix input.
#'
#' @details
#' Two file conventions are supported and distinguished automatically. If `x`
#' is a directory, or a vector of several paths, each file is taken to be one
#' species and `band` selects which variant to read from it. If `x` is a single
#' raster with more than one layer, the layers are taken to be the species and
#' `band` is ignored.
#'
#' All rasters must share one grid. The function checks extent, resolution and
#' CRS across every file and fails with an informative message naming the first
#' offender rather than silently recycling values.
#'
#' Species names are normalised (whitespace to underscore, trailing
#' underscores removed) so that a matrix read from a spreadsheet can be matched
#' against file names. [si_match_names()] reports the reconciliation.
#'
#' @examples
#' # From a matrix (no raster backend needed)
#' set.seed(1)
#' m <- matrix(runif(200 * 6), nrow = 200,
#'             dimnames = list(NULL, paste0("plant_", 1:6)))
#' st <- si_stack(m, level = "plant")
#' st
#'
#' # Mappable matrix input: the bundled Chilean subset
#' data(chile_pp)
#' P <- si_stack(chile_pp$plants_current_bin, level = "plant",
#'               grid = chile_pp$grid, cells = chile_pp$cells)
#' P
#'
#' # From a SpatRaster
#' library(terra)
#' r <- rast(nrows = 20, ncols = 20, nlyrs = 4)
#' values(r) <- runif(terra::ncell(r) * 4)
#' names(r) <- paste0("sp", 1:4)
#' si_stack(r, level = "pollinator")
#' @seealso [si_overlap()], [si_match_names()]
#' @export
si_stack <- function(x, band = 1, mask = NULL, names = NULL, level = "A",
                     na_value = 0, pattern = "\\.(tif|tiff|grd|img|asc)$",
                     grid = NULL, cells = NULL, quiet = FALSE) {

  if (is.matrix(x) || is.data.frame(x)) {
    v <- as.matrix(x)
    storage.mode(v) <- "double"
    nm <- if (!is.null(names)) names else colnames(v)
    if (is.null(nm)) nm <- paste0(level, "_", seq_len(ncol(v)))
    v[is.na(v)] <- na_value
    if (!is.null(cells) && length(cells) != nrow(v))
      stop("`cells` has length ", length(cells), " but the matrix has ",
           nrow(v), " rows.", call. = FALSE)
    return(new_si_stack(v, .norm_names(nm),
                        if (is.null(cells)) seq_len(nrow(v)) else as.integer(cells),
                        si_grid(grid), level))
  }

  if (inherits(x, "SpatRaster")) {
    r <- x
  } else if (is.character(x)) {
    files <- if (length(x) == 1L && dir.exists(x)) {
      f <- list.files(x, pattern = pattern, full.names = TRUE,
                      ignore.case = TRUE)
      if (!length(f))
        stop("No rasters matching '", pattern, "' found in: ", x, call. = FALSE)
      sort(f)
    } else {
      miss <- x[!file.exists(x)]
      if (length(miss)) {
        if (length(x) == 1L)
          stop("Path not found: '", x, "'\n",
               "  `x` should be a directory containing one raster per species, ",
               "a vector of raster file paths, a SpatRaster, or a matrix.\n",
               "  Current working directory is: ", getwd(), "\n",
               "  For a runnable example with no external files, see ",
               "?si_stack or data(chile_pp).", call. = FALSE)
        stop("File(s) not found: ", paste(utils::head(miss, 3), collapse = ", "),
             if (length(miss) > 3) " ..." else "", call. = FALSE)
      }
      x
    }
    # Two conventions are supported and disambiguated here:
    #   many files, one per species  -> `band` selects the variant in each
    #   one file, many layers        -> the layers are the species
    single_multilayer <- length(files) == 1L && terra::nlyr(terra::rast(files)) > 1L
    if (single_multilayer) {
      r <- terra::rast(files)
      if (!identical(band, 1) && !quiet)
        message("`x` is a single multi-layer raster, so its layers are taken ",
                "as species and `band` is ignored.")
      if (!quiet) message("Reading ", terra::nlyr(r), " ", level,
                          " layers from one raster ...")
      if (is.null(names)) names <- base::names(r)
    } else {
      if (!quiet) message("Reading ", length(files), " ", level, " rasters ...")
      r <- .read_band_stack(files, band)
      if (is.null(names)) names <- tools::file_path_sans_ext(basename(files))
    }
  } else {
    stop("`x` must be a directory, a vector of file paths, a SpatRaster, ",
         "or a matrix.", call. = FALSE)
  }

  nm <- if (!is.null(names)) names else base::names(r)
  if (length(nm) != terra::nlyr(r))
    stop("`names` has length ", length(nm), " but the stack has ",
         terra::nlyr(r), " layers.", call. = FALSE)

  keep <- .resolve_mask(r, mask)
  v <- terra::values(r, mat = TRUE)[keep, , drop = FALSE]
  v[is.na(v)] <- na_value
  storage.mode(v) <- "double"
  new_si_stack(v, .norm_names(nm), which(keep), si_grid(r[[1]]), level)
}

new_si_stack <- function(values, names, cells, template, level) {
  colnames(values) <- names
  structure(list(values = values, names = names, cells = cells,
                 template = template, level = level,
                 binary = .is_binary(values)),
            class = "si_stack")
}

#' @export
print.si_stack <- function(x, ...) {
  cat("<si_stack>  level: ", x$level, "\n", sep = "")
  cat("  species  : ", length(x$names), "  (",
      paste(utils::head(x$names, 3), collapse = ", "),
      if (length(x$names) > 3) ", ..." else "", ")\n", sep = "")
  cat("  cells    : ", nrow(x$values), "\n", sep = "")
  cat("  values   : ", if (x$binary) "binary (0/1)" else
    sprintf("continuous [%.3f, %.3f]", min(x$values), max(x$values)), "\n",
    sep = "")
  cat("  raster   : ", if (is.null(x$template)) "none (matrix input)" else
    paste0(x$template$dim[1], " x ", x$template$dim[2], ", ",
           x$template$crs_code), "\n", sep = "")
  invisible(x)
}

#' @export
summary.si_stack <- function(object, ...) {
  occ <- colSums(object$values > 0)
  data.frame(species = object$names,
             occupied_cells = occ,
             prevalence = occ / nrow(object$values),
             mean_suitability = colMeans(object$values),
             row.names = NULL)
}


#' Reconcile species names between two sources
#'
#' Interaction matrices assembled in a spreadsheet and raster file names
#' assembled by a modelling pipeline almost never agree exactly. This function
#' reports the reconciliation explicitly instead of letting a silent `reindex`
#' or `match` drop species.
#'
#' @param a,b Character vectors of species names.
#' @param normalise Logical. Apply the standard normalisation (trim, spaces to
#'   underscores, drop trailing underscores, case-fold) before matching.
#' @return A list with `matched`, `only_a`, `only_b` and a `map` data frame
#'   giving the correspondence for matched names.
#' @examples
#' si_match_names(c("Villa", "Apis mellifera"), c("Villa_", "Apis_mellifera"))
#' @export
si_match_names <- function(a, b, normalise = TRUE) {
  ka <- if (normalise) .norm_names(a) else a
  kb <- if (normalise) .norm_names(b) else b
  la <- tolower(ka); lb <- tolower(kb)
  m <- match(la, lb)
  list(matched = a[!is.na(m)],
       only_a  = a[is.na(m)],
       only_b  = b[!(lb %in% la)],
       map     = data.frame(a = a[!is.na(m)], b = b[m[!is.na(m)]],
                            row.names = NULL))
}


#' Build an interaction constraint layer (metaweb)
#'
#' The metaweb records which pairs of species *can* interact, independently of
#' where they occur. It is the object that prevents [si_overlap()] from
#' producing structurally degenerate networks (see the package overview).
#'
#' @param x A binary or weighted matrix with level-A species in rows and
#'   level-B species in columns; or a two-column `data.frame` of interaction
#'   records (edge list); or `NULL` for an unconstrained metaweb.
#' @param names_A,names_B Character vectors giving the full species pools. When
#'   `x` is an edge list these are required; when `x` is a matrix they are used
#'   to reorder and pad it, so that a metaweb assembled from a partial species
#'   list can be aligned to the modelled species.
#' @param weighted Logical. Keep the supplied weights rather than binarising.
#' @param symmetric Logical. For unipartite (competition) networks, force
#'   symmetry and zero the diagonal.
#'
#' @return An object of class `si_metaweb`.
#'
#' @examples
#' # from an edge list
#' el <- data.frame(plant = c("p1", "p1", "p2", "p3"),
#'                  poll  = c("a1", "a2", "a2", "a3"))
#' mw <- si_metaweb(el, names_A = paste0("p", 1:4), names_B = paste0("a", 1:3))
#' mw
#' si_connectance(mw)
#'
#' # an unconstrained metaweb: this is what makes networks degenerate
#' si_connectance(si_metaweb(NULL, names_A = paste0("p", 1:4),
#'                                 names_B = paste0("a", 1:3)))
#' @seealso [si_forbidden()], [si_simulate_metaweb()]
#' @export
si_metaweb <- function(x = NULL, names_A = NULL, names_B = NULL,
                       weighted = FALSE, symmetric = FALSE) {

  if (is.null(x)) {
    if (is.null(names_A) || is.null(names_B))
      stop("`names_A` and `names_B` are required when `x` is NULL.",
           call. = FALSE)
    m <- matrix(1, length(names_A), length(names_B),
                dimnames = list(.norm_names(names_A), .norm_names(names_B)))
    return(new_si_metaweb(m, TRUE))
  }

  if (is.data.frame(x) && ncol(x) >= 2L && !is.numeric(x[[1]])) {
    if (is.null(names_A)) names_A <- unique(x[[1]])
    if (is.null(names_B)) names_B <- unique(x[[2]])
    nA <- .norm_names(names_A); nB <- .norm_names(names_B)
    m <- matrix(0, length(nA), length(nB), dimnames = list(nA, nB))
    i <- match(.norm_names(x[[1]]), nA)
    j <- match(.norm_names(x[[2]]), nB)
    ok <- !is.na(i) & !is.na(j)
    if (!any(ok)) stop("No edge-list rows matched `names_A`/`names_B`.",
                       call. = FALSE)
    w <- if (ncol(x) >= 3L && is.numeric(x[[3]])) x[[3]][ok] else 1
    m[cbind(i[ok], j[ok])] <- w
  } else {
    m <- as.matrix(x)
    storage.mode(m) <- "double"
    if (is.null(rownames(m))) rownames(m) <- names_A
    if (is.null(colnames(m))) colnames(m) <- names_B
    rownames(m) <- .norm_names(rownames(m))
    colnames(m) <- .norm_names(colnames(m))
    if (!is.null(names_A) || !is.null(names_B)) {
      nA <- if (is.null(names_A)) rownames(m) else .norm_names(names_A)
      nB <- if (is.null(names_B)) colnames(m) else .norm_names(names_B)
      out <- matrix(0, length(nA), length(nB), dimnames = list(nA, nB))
      ri <- match(nA, rownames(m)); ci <- match(nB, colnames(m))
      okr <- !is.na(ri); okc <- !is.na(ci)
      out[okr, okc] <- m[ri[okr], ci[okc], drop = FALSE]
      m <- out
    }
  }
  m[is.na(m)] <- 0
  if (!weighted) m[] <- as.numeric(m > 0)
  if (symmetric) {
    m <- pmax(m, t(m)); diag(m) <- 0
  }
  new_si_metaweb(m, all(m == 1))
}

new_si_metaweb <- function(m, unconstrained) {
  structure(list(matrix = m,
                 names_A = rownames(m), names_B = colnames(m),
                 unconstrained = isTRUE(unconstrained)),
            class = "si_metaweb")
}

#' @export
print.si_metaweb <- function(x, ...) {
  m <- x$matrix
  cat("<si_metaweb>\n")
  cat("  dimensions : ", nrow(m), " x ", ncol(m), "\n", sep = "")
  cat("  links      : ", sum(m > 0), "\n", sep = "")
  cat("  connectance: ", format(round(mean(m > 0), 4), nsmall = 4), "\n", sep = "")
  if (x$unconstrained)
    cat("  ! UNCONSTRAINED: every pair is permitted. Networks built on this\n",
        "    metaweb are rank-one and structurally degenerate. See ?si_degenerate.\n",
        sep = "")
  invisible(x)
}

#' Connectance of a metaweb or matrix
#' @param x An [si_metaweb()], a matrix, or an [si_field()].
#' @return Numeric proportion of realised links.
#' @examples
#' si_connectance(matrix(c(1, 0, 1, 1), 2, 2))
#' @export
si_connectance <- function(x) {
  m <- if (inherits(x, "si_metaweb")) x$matrix else as.matrix(x)
  mean(m > 0)
}


# ---------------------------------------------------------------------------
# internal helpers
# ---------------------------------------------------------------------------

.norm_names <- function(x) {
  x <- as.character(x)
  # Real ecological data frequently carries invalid bytes and non-breaking
  # spaces from spreadsheets. Drop invalid bytes first so that downstream
  # regular expressions cannot fail on them.
  x <- iconv(x, from = "UTF-8", to = "UTF-8", sub = "")
  x[is.na(x)] <- ""
  x <- gsub("\xc2\xa0", " ", x, useBytes = TRUE)   # NBSP
  x <- gsub("\xe2\x80\x87|\xe2\x80\xaf|\xe2\x81\x9f|\xe3\x80\x80",
            " ", x, useBytes = TRUE)                # figure/narrow/ideographic
  x <- trimws(x)
  x <- gsub("[[:space:]]+", "_", x)
  sub("_+$", "", x)
}

.is_binary <- function(v) {
  u <- unique(as.vector(v[seq_len(min(length(v), 1e5))]))
  all(u %in% c(0, 1))
}

.read_band_stack <- function(files, band) {
  first <- terra::rast(files[1])
  b <- .resolve_band(first, band)
  out <- vector("list", length(files))
  ref <- NULL
  for (i in seq_along(files)) {
    r <- terra::rast(files[i])
    if (terra::nlyr(r) < b)
      stop("File '", basename(files[i]), "' has ", terra::nlyr(r),
           " band(s); band ", b, " requested.", call. = FALSE)
    ri <- r[[b]]
    if (is.null(ref)) {
      ref <- ri
    } else if (!terra::compareGeom(ref, ri, stopOnError = FALSE)) {
      stop("Raster '", basename(files[i]), "' does not share the grid of '",
           basename(files[1]), "'.\n  All suitability surfaces must have the ",
           "same extent, resolution and CRS.", call. = FALSE)
    }
    out[[i]] <- ri
  }
  terra::rast(out)
}

.resolve_band <- function(r, band) {
  if (is.numeric(band)) return(as.integer(band))
  nm <- base::names(r)
  b <- match(band, nm)
  if (is.na(b))
    stop("Band '", band, "' not found. Available: ",
         paste(nm, collapse = ", "), call. = FALSE)
  b
}

.resolve_mask <- function(r, mask) {
  n <- terra::ncell(r)
  if (is.null(mask)) {
    v <- terra::values(r[[1]], mat = FALSE)
    return(!is.na(v))
  }
  # An explicit mask defines the analysis grid exactly. It is deliberately NOT
  # intersected with the first layer's non-NA cells: different species have
  # different NA patterns, and intersecting would give each stack a different
  # set of cells, so fields could not be compared. NA cells inside the mask are
  # replaced by `na_value`.
  if (is.logical(mask)) {
    if (length(mask) != n)
      stop("Logical `mask` must have length ncell(r) = ", n, ".", call. = FALSE)
    return(mask)
  }
  if (inherits(mask, "SpatVector")) {
    mr <- terra::rasterize(terra::project(mask, terra::crs(r)), r[[1]],
                           field = 1, background = NA)
    return(!is.na(terra::values(mr, mat = FALSE)))
  }
  if (inherits(mask, "SpatRaster")) {
    v <- terra::values(mask[[1]], mat = FALSE)
    return(!is.na(v) & v > 0)
  }
  stop("`mask` must be NULL, a logical vector, a SpatVector or a SpatRaster.",
       call. = FALSE)
}


#' Build a common analysis mask across several raster sets
#'
#' Species distribution projections for different time slices are often written
#' with slightly different `NA` masks, because a future climate layer may be
#' undefined where the baseline is defined. If the current and future analyses
#' are then run on different cell sets, cells that merely lack a future
#' prediction are silently counted as cells that have lost their network.
#'
#' `si_common_mask()` returns the intersection of valid cells across every
#' supplied source and reports how many cells each one contributes, so that the
#' discrepancy is visible rather than absorbed into the results.
#'
#' @param ... Directories, file paths or `SpatRaster` objects. One
#'   representative raster per set is enough.
#' @param region Optional `SpatVector` or `SpatRaster` study-area mask.
#' @param band Band to inspect.
#' @param quiet Logical. Suppress the report.
#' @return A logical vector over the cells of the reference grid, with a
#'   `report` attribute.
#' @examples
#' library(terra)
#' r1 <- rast(nrows = 10, ncols = 10); values(r1) <- 1
#' r2 <- rast(nrows = 10, ncols = 10); values(r2) <- 1; r2[1:15] <- NA
#' m <- si_common_mask(r1, r2)
#' attr(m, "report")
#' @seealso [si_stack()]
#' @export
si_common_mask <- function(..., region = NULL, band = 1, quiet = FALSE) {
  src <- list(...)
  if (!length(src)) stop("Supply at least one raster source.", call. = FALSE)
  reps <- lapply(src, function(s) {
    if (inherits(s, "SpatRaster")) return(s[[min(band, terra::nlyr(s))]])
    f <- if (length(s) == 1L && dir.exists(s))
      sort(list.files(s, pattern = "\\.(tif|tiff|grd|img|asc)$",
                      full.names = TRUE, ignore.case = TRUE))[1] else s[1]
    if (is.na(f) || !file.exists(f))
      stop("No raster found for source: ", s, call. = FALSE)
    terra::rast(f)[[band]]
  })
  ref <- reps[[1]]
  valid <- lapply(reps, function(r) !is.na(terra::values(r, mat = FALSE)))
  reg <- if (is.null(region)) rep(TRUE, terra::ncell(ref)) else
    .resolve_mask(ref, region)
  keep <- Reduce(`&`, valid) & reg
  rep_tab <- data.frame(
    source = vapply(src, function(s)
      if (inherits(s, "SpatRaster")) "SpatRaster" else basename(s[1]),
      character(1)),
    valid_in_region = vapply(valid, function(v) sum(v & reg), numeric(1)),
    dropped_vs_common = vapply(valid, function(v) sum(v & reg) - sum(keep),
                               numeric(1)))
  if (!quiet) {
    message(sprintf("common mask: %d cells", sum(keep)))
    if (any(rep_tab$dropped_vs_common > 0)) {
      message("  sources differ in their valid cells:")
      for (i in which(rep_tab$dropped_vs_common > 0))
        message(sprintf("    %-28s %d valid, %d not shared",
                        rep_tab$source[i], rep_tab$valid_in_region[i],
                        rep_tab$dropped_vs_common[i]))
      message("  Cells outside the common mask would otherwise be scored as ",
              "lost networks.")
    }
  }
  structure(keep, report = rep_tab)
}


#' A serialisable description of a raster grid
#'
#' `terra` objects hold external pointers and therefore do not survive
#' `saveRDS()`. `spinet` stores the analysis grid as a small description that
#' can be saved, shared and rebuilt, so that fields, metrics and experiment
#' results are fully serialisable.
#'
#' @param x A `SpatRaster`, or an `si_grid` (returned unchanged).
#' @return An object of class `si_grid`.
#' @examples
#' library(terra)
#' g <- si_grid(rast(nrows = 10, ncols = 12))
#' g
#' si_grid_rast(g)
#' # unlike a SpatRaster, this survives a round trip
#' f <- tempfile(); saveRDS(g, f); identical(readRDS(f)$dim, g$dim)
#' @export
si_grid <- function(x) {
  if (inherits(x, "si_grid")) return(x)
  if (is.null(x)) return(NULL)
  stopifnot(inherits(x, "SpatRaster"))
  structure(list(dim = c(terra::nrow(x), terra::ncol(x)),
                 ext = as.vector(terra::ext(x)),
                 wkt = terra::crs(x),
                 crs_code = {
                   d <- terra::crs(x, describe = TRUE)$code
                   if (is.na(d)) "unknown CRS" else d
                 }),
            class = "si_grid")
}

#' Rebuild a SpatRaster from a grid description
#' @param g An [si_grid()].
#' @return An empty single-layer `SpatRaster`.
#' @examples
#' library(terra)
#' si_grid_rast(si_grid(rast(nrows = 5, ncols = 5)))
#' @export
si_grid_rast <- function(g) {
  if (is.null(g)) return(NULL)
  if (inherits(g, "SpatRaster")) return(g)
  r <- terra::rast(nrows = g$dim[1], ncols = g$dim[2],
                   xmin = g$ext[1], xmax = g$ext[2],
                   ymin = g$ext[3], ymax = g$ext[4], crs = g$wkt)
  terra::values(r) <- NA_real_
  r
}

#' @export
print.si_grid <- function(x, ...) {
  cat("<si_grid> ", x$dim[1], " x ", x$dim[2], "   ", x$crs_code, "\n", sep = "")
  invisible(x)
}


#' Paths to the bundled example rasters
#'
#' The package ships a small multi-band raster subset of the Chilean system:
#' 20 plants and 20 pollinators, present and 2070, each as a continuous
#' suitability surface and a thresholded binary surface, plus the study-area
#' polygon and the documented metaweb. These files exist so that the raster
#' reading path can be demonstrated and tested without external data;
#' [chile_pp] holds the same information as matrices.
#'
#' @param what Which file to locate: `"plants_current"`, `"plants_future"`,
#'   `"pollinators_current"`, `"pollinators_future"` (multi-layer GeoTIFFs of
#'   continuous suitability), the corresponding `"_binary"` versions,
#'   `"region"` (the study-area polygon) or `"metaweb"` (a CSV). `"all"`
#'   returns every path.
#' @return A file path, or a named character vector when `what = "all"`.
#' @examples
#' si_example_rasters("all")
#'
#' library(terra)
#' r <- rast(si_example_rasters("plants_current"))
#' c(layers = terra::nlyr(r), cells = terra::ncell(r))
#'
#' # the raster workflow, end to end, with no external files
#' region <- terra::vect(si_example_rasters("region"))
#' P <- si_stack(si_example_rasters("plants_current_binary"),
#'               mask = region, level = "plant", quiet = TRUE)
#' A <- si_stack(si_example_rasters("pollinators_current_binary"),
#'               mask = region, level = "pollinator", quiet = TRUE)
#' mw <- si_metaweb(as.matrix(read.csv(si_example_rasters("metaweb"),
#'                                     row.names = 1)))
#' f <- si_overlap(P, A, mw, method = "binary", check = FALSE)
#' f
#' @seealso [chile_pp], [si_stack()], [si_common_mask()]
#' @export
si_example_rasters <- function(what = "all") {
  files <- c(
    plants_current             = "plants_current_suitability.tif",
    plants_future              = "plants_future_suitability.tif",
    pollinators_current        = "pollinators_current_suitability.tif",
    pollinators_future         = "pollinators_future_suitability.tif",
    plants_current_binary      = "plants_current_binary.tif",
    plants_future_binary       = "plants_future_binary.tif",
    pollinators_current_binary = "pollinators_current_binary.tif",
    pollinators_future_binary  = "pollinators_future_binary.tif",
    region                     = "chile.gpkg",
    metaweb                    = "chile_metaweb.csv")
  paths <- vapply(files, function(f)
    system.file("extdata", f, package = "spinet"), character(1))
  if (identical(what, "all")) return(paths)
  if (!what %in% names(paths))
    stop("Unknown file '", what, "'. Available: ",
         paste(names(paths), collapse = ", "), call. = FALSE)
  unname(paths[what])
}
