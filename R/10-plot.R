# ---------------------------------------------------------------------------
# Module 9: visualisation
# ---------------------------------------------------------------------------

#' Plot a local interaction network
#'
#' Draws a bipartite interaction matrix as a two-row graph, with node width
#' proportional to marginal interaction strength and link width proportional to
#' interaction weight.
#'
#' @param x An [si_field()], or an interaction matrix.
#' @param cell Cell index when `x` is a field.
#' @param sort Logical. Order species by marginal strength, which makes
#'   nestedness visible.
#' @param col_A,col_B,col_link Colours.
#' @param labels Logical. Draw species labels.
#' @param ... Passed to `plot`.
#' @return Invisibly, the plotted matrix.
#' @examples
#' f <- si_example_field()
#' si_plot_network(f, cell = 4)
#' @export
si_plot_network <- function(x, cell = 1, sort = TRUE,
                            col_A = "#2E7D32", col_B = "#EF6C00",
                            col_link = "#37474F", labels = TRUE, ...) {
  W <- if (inherits(x, "si_field")) .local_matrix(x, cell, drop = TRUE) else
    as.matrix(x)
  if (!nrow(W) || !ncol(W)) { graphics::plot.new(); graphics::title("empty network"); return(invisible(W)) }
  if (sort) W <- W[order(-rowSums(W)), order(-colSums(W)), drop = FALSE]
  nA <- nrow(W); nB <- ncol(W)
  xa <- seq(0, 1, length.out = nA + 2)[2:(nA + 1)]
  xb <- seq(0, 1, length.out = nB + 2)[2:(nB + 1)]
  op <- graphics::par(mar = c(3, 1, 3, 1)); on.exit(graphics::par(op))
  graphics::plot(NA, xlim = c(-0.05, 1.05), ylim = c(-0.25, 1.25),
                 axes = FALSE, xlab = "", ylab = "", ...)
  mx <- max(W)
  ij <- which(W > 0, arr.ind = TRUE)
  for (r in seq_len(nrow(ij))) {
    i <- ij[r, 1]; j <- ij[r, 2]
    graphics::segments(xa[i], 1, xb[j], 0,
                       lwd = 0.4 + 3 * W[i, j] / mx,
                       col = grDevices::adjustcolor(col_link, 0.45))
  }
  ra <- rowSums(W) / max(rowSums(W)); rb <- colSums(W) / max(colSums(W))
  graphics::points(xa, rep(1, nA), pch = 15, cex = 0.8 + 2 * ra, col = col_A)
  graphics::points(xb, rep(0, nB), pch = 15, cex = 0.8 + 2 * rb, col = col_B)
  if (labels && nA <= 30)
    graphics::text(xa, 1.10, rownames(W), srt = 90, adj = 0, cex = 0.55)
  if (labels && nB <= 30)
    graphics::text(xb, -0.10, colnames(W), srt = 90, adj = 1, cex = 0.55)
  invisible(W)
}

#' @export
plot.si_field <- function(x, cell = 1, ...) si_plot_network(x, cell = cell, ...)

#' Plot a matrix as a shaded grid
#'
#' @param m A matrix.
#' @param sort Logical. Sort rows and columns by marginal totals.
#' @param col Colour ramp function.
#' @param main Title.
#' @param ... Passed to `image`.
#' @return Invisibly, the plotted matrix.
#' @examples
#' mw <- si_simulate_metaweb(30, 25, 0.2, architecture = "nested", seed = 1)
#' si_plot_matrix(mw$matrix, main = "nested metaweb")
#' @export
si_plot_matrix <- function(m, sort = TRUE, col = NULL, main = "", ...) {
  W <- as.matrix(m)
  if (sort) W <- W[order(-rowSums(W)), order(-colSums(W)), drop = FALSE]
  if (is.null(col)) col <- grDevices::colorRampPalette(
    c("#FFFFFF", "#BBDEFB", "#1E88E5", "#0D47A1"))(64)
  op <- graphics::par(mar = c(3, 3, 3, 1)); on.exit(graphics::par(op))
  graphics::image(t(W[nrow(W):1, , drop = FALSE]), col = col, axes = FALSE,
                  main = main, ...)
  graphics::box(col = "grey70")
  invisible(W)
}

#' Plot the response surface of an experiment
#'
#' @param experiment An [si_experiment()] result.
#' @param response Column name of the response.
#' @param x Design factor for the horizontal axis.
#' @param group Optional design factor mapped to colour.
#' @param fun Summary function applied within each combination.
#' @param ... Passed to `plot`.
#' @return Invisibly, the summarised table.
#' @examples
#' des <- list(richness = 20, connectance = c(1, 0.6, 0.3, 0.1),
#'             severity = c(0.2, 0.5), grid = 12)
#' ex <- suppressWarnings(si_experiment(des, replicates = 2, cells = 30, seed = 4))
#' si_plot_experiment(ex, "beta_WN", x = "connectance", group = "severity")
#' @export
si_plot_experiment <- function(experiment, response, x = "connectance",
                               group = NULL, fun = mean, ...) {
  d <- as.data.frame(experiment)
  if (!response %in% names(d)) stop("Unknown response '", response, "'.",
                                    call. = FALSE)
  if (is.null(group)) {
    agg <- stats::aggregate(d[[response]], by = list(x = d[[x]]), FUN = fun,
                            na.rm = TRUE)
    names(agg)[2] <- "y"
    graphics::plot(agg$x, agg$y, type = "b", pch = 19, xlab = x,
                   ylab = response, ...)
    return(invisible(agg))
  }
  agg <- stats::aggregate(d[[response]],
                          by = list(x = d[[x]], g = d[[group]]), FUN = fun,
                          na.rm = TRUE)
  names(agg)[3] <- "y"
  gs <- unique(agg$g)
  cols <- grDevices::hcl.colors(max(2, length(gs)), "Dark 3")
  graphics::plot(range(agg$x), range(agg$y, na.rm = TRUE), type = "n",
                 xlab = x, ylab = response, ...)
  for (i in seq_along(gs)) {
    s <- agg[agg$g == gs[i], ]
    s <- s[order(s$x), ]
    graphics::lines(s$x, s$y, type = "b", pch = 19, col = cols[i], lwd = 2)
  }
  graphics::legend("topright", legend = paste0(group, " = ", gs),
                   col = cols[seq_along(gs)], lwd = 2, bty = "n", cex = 0.8)
  invisible(agg)
}

#' A perceptually ordered palette for change maps
#'
#' Diverging palette centred on zero, for mapping current-to-future change.
#'
#' @param n Number of colours.
#' @return A character vector of hex colours.
#' @examples
#' si_palette_change(5)
#' @export
si_palette_change <- function(n = 11) {
  grDevices::colorRampPalette(
    c("#762A83", "#9970AB", "#C2A5CF", "#E7D4E8", "#F7F7F7",
      "#D9F0D3", "#A6DBA0", "#5AAE61", "#1B7837"))(n)
}

#' A sequential palette for network metrics
#' @param n Number of colours.
#' @return A character vector of hex colours.
#' @examples
#' si_palette_metric(5)
#' @export
si_palette_metric <- function(n = 11) {
  grDevices::colorRampPalette(
    c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59",
      "#EF6548", "#D7301F", "#B30000", "#7F0000"))(n)
}
