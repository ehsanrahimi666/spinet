#' spinet: Spatial Interaction Networks Under Environmental Change
#'
#' @description
#' `spinet` projects ecological interaction networks through space and time by
#' coupling per-species habitat suitability surfaces (from any species
#' distribution model) with an interaction constraint layer. It is deliberately
#' *downstream* of distribution modelling: it consumes SDM output and never
#' fits distribution models itself.
#'
#' @section Why a constraint layer is mandatory:
#' If a local interaction network is built from co-occurrence alone, the local
#' interaction matrix is the outer product of the two suitability vectors,
#' \eqn{W_k = p_k \otimes a_k}. Such a matrix has rank one, which has three
#' exact consequences:
#' \itemize{
#'   \item in the binary pathway every local network is the complete bipartite
#'         graph, so connectance is identically 1 and nestedness identically 0;
#'   \item in the continuous pathway \eqn{H_2'} is identically 0, because the
#'         null expectation of \eqn{H_2'} *is* the outer product of the
#'         marginal totals;
#'   \item interaction rewiring \eqn{\beta_{OS}} is identically 0 whenever the
#'         constraint layer is time-invariant, because species shared between
#'         two realisations necessarily interact identically.
#' }
#' `spinet` therefore treats the constraint layer as the primary scientific
#' object rather than as a refinement, and [si_degenerate()] warns whenever a
#' requested metric is mathematically constant for the field supplied.
#'
#' @section Main workflow:
#' \enumerate{
#'   \item [si_stack()] -- read per-species suitability rasters.
#'   \item [si_metaweb()] / [si_forbidden()] -- build the constraint layer.
#'   \item [si_overlap()] -- construct the spatial interaction field.
#'   \item [si_metrics()] -- per-cell network metrics.
#'   \item [si_beta()], [si_rewiring()] -- change between two time slices.
#'   \item [si_robustness()] -- coextinction robustness with rewiring.
#'   \item [si_experiment()] -- factorial virtual experiments.
#' }
#'
#' @references
#' Poisot, T., Canard, E., Mouillot, D., Mouquet, N. & Gravel, D. (2012)
#' The dissimilarity of species interaction networks. *Ecology Letters*,
#' 15, 1353--1361. \doi{10.1111/ele.12002}
#'
#' Vazquez, D.P., Chacoff, N.P. & Cagnolo, L. (2009) Evaluating multiple
#' determinants of the structure of plant-animal mutualistic networks.
#' *Ecology*, 90, 2039--2046. \doi{10.1890/08-1837.1}
#'
#' Vizentin-Bugoni, J., Debastiani, V.J., Bastazini, V.A.G., Maruyama, P.K. &
#' Sperry, J.H. (2020) Including rewiring in the estimation of the robustness
#' of mutualistic networks. *Methods in Ecology and Evolution*, 11, 106--116.
#' \doi{10.1111/2041-210X.13306}
#'
#' @keywords internal
"_PACKAGE"

# ---------------------------------------------------------------------------
# Interaction-type vocabulary
# ---------------------------------------------------------------------------

#' Relations Ontology interaction vocabulary
#'
#' The lookup table used by [si_interaction()] to translate a Relations
#' Ontology (RO) term, or a common English label, into the sign structure and
#' cascade rule of an ecological interaction. The vocabulary follows the terms
#' used by Global Biotic Interactions (GloBI), so that any GloBI-standard
#' interaction table can be read by `spinet` without a bespoke mapping.
#'
#' @format A `data.frame` with columns
#' \describe{
#'   \item{label}{canonical English label}
#'   \item{ro_id}{Relations Ontology IRI, where one exists}
#'   \item{sign_A}{effect of level B on level A: `"+"`, `"-"` or `"0"`}
#'   \item{sign_B}{effect of level A on level B}
#'   \item{sign}{the two-character sign structure, e.g. `"++"`}
#'   \item{mode}{`"bipartite"` or `"unipartite"`}
#'   \item{cascade}{how partner loss propagates; see [si_interaction()]}
#' }
#' @examples
#' head(si_vocabulary, 12)
#' # all the antagonistic terms
#' si_vocabulary[si_vocabulary$sign == "+-", c("label", "ro_id")]
#' @seealso [si_interaction()]
#' @export
si_vocabulary <- data.frame(
  label = c(
    "pollinatedBy", "hasVector", "mutualistOf", "hasHost",
    "hasFlowersVisitedBy", "visitedBy", "commensalistOf", "epiphyteOf",
    "eatenBy", "preyedUponBy", "parasitizedBy", "hostOf", "pathogenOf",
    "amensalistOf", "allelopathOf",
    "competesWith",
    "neutralWith", "interactsWith"),
  ro_id = c(
    "RO_0002456", "RO_0002460", "RO_0002442", "RO_0002454",
    "RO_0002623", "RO_0002623", "RO_0002441", "RO_0008501",
    "RO_0002471", "RO_0002458", "RO_0002456", "RO_0002453", "RO_0002556",
    NA_character_, NA_character_,
    "RO_0002448",
    NA_character_, "RO_0002437"),
  sign_A = c("+", "+", "+", "+",
             "0", "0", "0", "0",
             "-", "-", "-", "-", "-",
             "-", "-",
             "-",
             "0", NA),
  sign_B = c("+", "+", "+", "+",
             "+", "+", "+", "+",
             "+", "+", "+", "+", "+",
             "0", "0",
             "-",
             "0", NA),
  mode = c(rep("bipartite", 15), "unipartite", "bipartite", "bipartite"),
  cascade = c(rep("mutual", 4),
              rep("one_way", 4),
              rep("consumer", 5),
              rep("release", 2),
              "competitive_release",
              "none", "mutual"),
  stringsAsFactors = FALSE
)
si_vocabulary$sign <- paste0(si_vocabulary$sign_A, si_vocabulary$sign_B)
si_vocabulary$sign[is.na(si_vocabulary$sign_A)] <- NA_character_
si_vocabulary <- si_vocabulary[, c("label", "ro_id", "sign_A", "sign_B",
                                   "sign", "mode", "cascade")]


#' Define an ecological interaction type
#'
#' Creates the object that tells every other `spinet` function what kind of
#' interaction it is handling. The interaction type controls three things that
#' cannot be inferred from the data: the sign structure (who benefits), the
#' direction in which the loss of a partner propagates, and whether increasing
#' spatial overlap is good news or bad news for the species involved.
#'
#' @param type Character. Either a canonical label or a Relations Ontology
#'   identifier from [si_vocabulary] (e.g. `"pollinatedBy"`, `"RO_0002456"`),
#'   or one of the six shorthand sign structures `"++"`, `"+-"`, `"+0"`,
#'   `"-0"`, `"--"`, `"00"`.
#' @param level_A,level_B Character. Names for the two levels, used in
#'   printing and in output column names. Defaults `"A"` and `"B"`; for a
#'   pollination network you would typically use `"plant"` and `"pollinator"`.
#' @param cascade Character, optionally overriding the vocabulary default.
#'   One of:
#'   \describe{
#'     \item{`"mutual"`}{a species in either level is lost when it loses all
#'       partners (the classic mutualistic coextinction rule).}
#'     \item{`"consumer"`}{only the consumer (level B) coextincts; the resource
#'       (level A) is *released* by the loss of its consumers.}
#'     \item{`"one_way"`}{only the beneficiary coextincts; the other level is
#'       indifferent (commensalism).}
#'     \item{`"release"`}{the harmed level is released; nothing coextincts
#'       (amensalism).}
#'     \item{`"competitive_release"`}{loss of a partner increases the
#'       persistence of the survivor (competition).}
#'     \item{`"none"`}{no propagation (neutralism); used as a null benchmark.}
#'   }
#' @param overlap_sign Numeric, `+1` or `-1`, optionally overriding the
#'   default. `+1` means that greater spatial overlap implies a stronger or
#'   more beneficial interaction (mutualism, commensalism); `-1` means that
#'   greater overlap implies greater pressure (competition, amensalism, and
#'   the resource side of an antagonism). This flips the interpretation of
#'   every change map produced by [si_compare()].
#'
#' @return An object of class `si_interaction`.
#'
#' @details
#' The default `overlap_sign` is `+1` for `"++"`, `"+0"`, and `"+-"` (where
#' the analysis is normally framed from the consumer's perspective) and `-1`
#' for `"--"` and `"-0"`. Set it explicitly when the framing matters: for a
#' crop-pest analysis, for example, you want `overlap_sign = -1` so that
#' increasing future overlap is reported as a deterioration.
#'
#' @examples
#' # a plant-pollinator mutualism
#' pol <- si_interaction("pollinatedBy", level_A = "plant",
#'                       level_B = "pollinator")
#' pol
#'
#' # the same thing addressed by its ontology identifier
#' identical(si_interaction("RO_0002456")$sign, si_interaction("++")$sign)
#'
#' # a crop-pest antagonism, framed so that more overlap is worse
#' pest <- si_interaction("eatenBy", level_A = "crop", level_B = "pest",
#'                        overlap_sign = -1)
#' pest
#'
#' # competition is a one-mode problem
#' si_interaction("competesWith")$mode
#' @seealso [si_vocabulary], [si_cascade_rule()]
#' @export
si_interaction <- function(type = "pollinatedBy",
                           level_A = "A", level_B = "B",
                           cascade = NULL, overlap_sign = NULL) {
  if (!is.character(type) || length(type) != 1L)
    stop("`type` must be a single character string.", call. = FALSE)
  v <- si_vocabulary
  key <- trimws(type)
  idx <- which(tolower(v$label) == tolower(key))
  if (!length(idx)) idx <- which(!is.na(v$ro_id) & toupper(v$ro_id) == toupper(sub("^.*/", "", key)))
  if (!length(idx) && nchar(key) == 2L && all(strsplit(key, "")[[1]] %in% c("+", "-", "0"))) {
    rev_key <- paste0(rev(strsplit(key, "")[[1]]), collapse = "")
    hit <- which(!is.na(v$sign) & (v$sign == key | v$sign == rev_key))
    idx <- hit[1L]
  }
  if (!length(idx) || all(is.na(idx))) {
    stop("Unknown interaction type: '", type, "'.\n",
         "  Use one of the labels or RO ids in `si_vocabulary`, or a sign ",
         "structure: \"++\", \"+-\", \"+0\", \"-0\", \"--\", \"00\".",
         call. = FALSE)
  }
  idx <- idx[1L]
  row <- v[idx, ]
  if (is.na(row$sign))
    stop("`", type, "` is a generic term with no defined sign structure. ",
         "Supply a specific interaction type.", call. = FALSE)

  cas <- if (is.null(cascade)) row$cascade else match.arg(
    cascade, c("mutual", "consumer", "one_way", "release",
               "competitive_release", "none"))
  osign <- if (!is.null(overlap_sign)) {
    if (!overlap_sign %in% c(-1, 1))
      stop("`overlap_sign` must be -1 or 1.", call. = FALSE)
    as.integer(overlap_sign)
  } else {
    if (row$sign %in% c("--", "-0")) -1L else 1L
  }

  structure(list(
    type         = row$label,
    ro_id        = row$ro_id,
    sign         = row$sign,
    sign_A       = row$sign_A,
    sign_B       = row$sign_B,
    mode         = row$mode,
    cascade      = cas,
    overlap_sign = osign,
    level_A      = level_A,
    level_B      = level_B
  ), class = "si_interaction")
}

#' @export
print.si_interaction <- function(x, ...) {
  cat("<si_interaction>\n")
  cat("  type         : ", x$type,
      if (!is.na(x$ro_id)) paste0("  (", x$ro_id, ")") else "", "\n", sep = "")
  cat("  sign         : ", x$sign,
      "   [", x$level_A, " ", x$sign_A, " / ", x$level_B, " ", x$sign_B, "]\n", sep = "")
  cat("  mode         : ", x$mode, "\n", sep = "")
  cat("  cascade      : ", x$cascade, "\n", sep = "")
  cat("  overlap sign : ", if (x$overlap_sign > 0)
    "+1  (more overlap = stronger/beneficial)" else
      "-1  (more overlap = more pressure)", "\n", sep = "")
  invisible(x)
}

#' Cascade rule for an interaction type
#'
#' Returns, for each level, whether losing all partners causes loss of the
#' species (`"extinct"`), no change (`"neutral"`), or a benefit
#' (`"released"`). This is the single place in `spinet` where the sign
#' structure is translated into extinction dynamics, and it is what makes the
#' package applicable to antagonistic and competitive networks rather than
#' only to mutualisms.
#'
#' @param interaction An [si_interaction()] object.
#' @return A named character vector of length two, `A` and `B`.
#' @examples
#' si_cascade_rule(si_interaction("++"))   # both coextinct
#' si_cascade_rule(si_interaction("+-"))   # consumer dies, resource released
#' si_cascade_rule(si_interaction("+0"))   # only the beneficiary dies
#' si_cascade_rule(si_interaction("--"))   # both released
#' @export
si_cascade_rule <- function(interaction) {
  stopifnot(inherits(interaction, "si_interaction"))
  switch(interaction$cascade,
    mutual              = c(A = "extinct",  B = "extinct"),
    consumer            = c(A = "released", B = "extinct"),
    one_way             = c(A = "neutral",  B = "extinct"),
    release             = c(A = "released", B = "neutral"),
    competitive_release = c(A = "released", B = "released"),
    none                = c(A = "neutral",  B = "neutral"))
}


#' Chilean plant-pollinator subset
#'
#' A 20 x 20 subset of the Chilean plant-pollinator system used to demonstrate
#' the package: the twenty best-connected plants and the twenty best-connected
#' pollinators from the catalogue of Muschett and Fonturbel (2022), with
#' habitat suitability projected for the present and for 2070 under SSP585.
#'
#' @format A list with:
#' \describe{
#'   \item{plants_current, plants_future}{2800 x 20 matrices of continuous
#'     habitat suitability (MaxEnt cloglog).}
#'   \item{pollinators_current, pollinators_future}{2800 x 20 matrices.}
#'   \item{plants_current_bin, plants_future_bin,
#'     pollinators_current_bin, pollinators_future_bin}{the same surfaces
#'     thresholded at maximum sensitivity.}
#'   \item{metaweb}{20 x 20 binary matrix of documented interactions.}
#'   \item{phenology_plants, phenology_pollinators}{12 x 20 month-by-species
#'     activity matrices derived from the catalogue's dated records.}
#'   \item{cells}{raster cell indices of the 2800 analysed cells.}
#'   \item{grid}{an [si_grid()] describing the 10 arc-minute WGS84 grid.}
#' }
#' @source Species distribution models fitted with MaxEnt to GBIF occurrences;
#'   interactions from Muschett, G. & Fonturbel, F.E. (2022) A comprehensive
#'   catalogue of plant-pollinator interactions for Chile. *Scientific Data*,
#'   9, 78. \doi{10.1038/s41597-022-01195-8}
#' @examples
#' data(chile_pp)
#' str(chile_pp$metaweb)
#' P <- si_stack(chile_pp$plants_current_bin, level = "plant")
#' A <- si_stack(chile_pp$pollinators_current_bin, level = "pollinator")
#' P$template <- chile_pp$grid; P$cells <- chile_pp$cells
#' A$template <- chile_pp$grid; A$cells <- chile_pp$cells
#' mw <- si_metaweb(chile_pp$metaweb)
#' f <- si_overlap(P, A, mw, method = "binary", check = FALSE)
#' f
#' summary(si_metrics(f, what = c("links", "connectance", "NODF")))
"chile_pp"
