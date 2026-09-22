## Supplementary Code S3. Scripts that produce Figures 1-6. Set FIG_DATA to
## Supplementary_Data_S2, FIG_OUT to an output folder and CHILE_OUTLINE to a
## national-boundary polygon (optional, used for map outlines).
## ---------------------------------------------------------------------------
## Figures 1-6 for the spinet manuscript, drawn at final print size
## (6.5 in wide, 9-pt base font) so that text is legible on the page.
## Figure 5 = sign structures and simulation engine; Figure 6 = Chilean example.
## ---------------------------------------------------------------------------
invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))
suppressMessages({library(spinet); library(terra)})
DATA <- Sys.getenv("FIG_DATA", "data")
FIG  <- Sys.getenv("FIG_OUT", "figures"); dir.create(FIG, showWarnings = FALSE)
OUTLINE <- Sys.getenv("CHILE_OUTLINE", "")

GREY <- "#3C3C3C"; LG <- "#D9D9D9"
PAL  <- c("#1F4E79", "#C0392B", "#1E8449", "#B7950B", "#7D3C98", "#566573")
dev_open <- function(f, w, h, ps = 9)
  png(file.path(FIG, f), width = w, height = h, units = "in", res = 600,
      pointsize = ps, type = "cairo")
base_par <- function(...) par(family = "sans", las = 1, bty = "n", tcl = -0.25,
                              mgp = c(2.0, 0.5, 0), cex.axis = 0.95, cex.lab = 1,
                              col.axis = GREY, fg = GREY, ...)
title_panel <- function(lab, txt, cex = 1)
  mtext(bquote(bold(.(lab)) ~~ .(txt)), side = 3, adj = 0, line = 0.6,
        cex = cex, col = "black")
sem  <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

## ===========================================================================
## Figure 1 - workflow
## ===========================================================================
rrect <- function(x0, y0, x1, y1, col, border = NA, r = 1.2, lty = 1, lwd = 1) {
  t <- seq(0, pi / 2, length.out = 8)
  xs <- c(x1 - r + r * cos(t), x0 + r - r * sin(t), x0 + r - r * cos(t), x1 - r + r * sin(t))
  ys <- c(y1 - r + r * sin(t), y1 - r + r * cos(t), y0 + r - r * sin(t), y0 + r - r * cos(t))
  polygon(xs, ys, col = col, border = border, lty = lty, lwd = lwd)
}
wbox <- function(x0, y0, x1, y1, col, num, head, sub = NULL, tc = "white", hc = 0.92, sc = 0.80) {
  rrect(x0, y0, x1, y1, col)
  cy <- if (is.null(sub)) (y0 + y1) / 2 else y1 - 2.4
  text((x0 + x1) / 2, cy, paste0(num, "  ", head), font = 2, cex = hc, col = tc)
  if (!is.null(sub)) for (i in seq_along(sub))
    text((x0 + x1) / 2, cy - 2.6 * i, sub[i], cex = sc, col = tc)
}
arr <- function(x0, y0, x1, y1, col = GREY, lwd = 1.1)
  arrows(x0, y0, x1, y1, length = 0.05, col = col, lwd = lwd)

dev_open("Figure1.png", 6.5, 5.4, ps = 8.5)
par(mar = c(0.2, 0.2, 0.2, 0.2), family = "sans")
plot(NA, xlim = c(0, 100), ylim = c(0, 100), axes = FALSE, xlab = "", ylab = "")
EXT <- "#DCEBF7"; PREP <- "#1F4E79"; CON <- "#C0392B"; SIMC <- "#7D3C98"
FLD <- "#1E8449"; DIA <- "#7B241C"; OUT <- "#566573"
MOD <- c("#2E75B6", "#CA6F1E", "#148F77", "#9A7D0A")
# 1 inputs
rrect(2, 84, 74, 99, EXT, border = "#7FA7C9", lty = 2, lwd = 1.2)
text(38, 96.7, "1  Inputs (external data, not part of spinet)", font = 2, cex = 0.92, col = PREP)
ib <- list(c(3.5, 19.5, "Occurrences", "GBIF, surveys"),
           c(22.5, 37.5, "SDM outputs", "e.g. flexsdm, biomod2"),
           c(40.5, 55.5, "Interaction records", "catalogues, GloBI"),
           c(58, 72.5, "Traits, phenology", "field or literature data"))
for (b in ib) {
  x0 <- as.numeric(b[1]); x1 <- as.numeric(b[2])
  rrect(x0, 85.3, x1, 94, "white", border = "#9DBBD6")
  text((x0 + x1) / 2, 91.5, b[3], font = 2, cex = 0.8, col = PREP)
  text((x0 + x1) / 2, 88.0, b[4], cex = 0.7, col = GREY)
}
arr(19.7, 89.6, 22.3, 89.6)
# legend
rrect(77, 70.5, 99, 99, "white", border = LG)
text(88, 96.8, "Legend", font = 2, cex = 0.82, col = GREY)
leg <- list(c(EXT, "External inputs"), c(PREP, "Data preparation"), c(CON, "Constraint layer"),
            c(SIMC, "Virtual simulation"), c(FLD, "Interaction field"), c(DIA, "Diagnostic"),
            c(MOD[1], "Analysis modules"), c(OUT, "Outputs"))
for (i in seq_along(leg)) {
  y <- 94.2 - (i - 1) * 2.95
  rect(78.6, y - 1.0, 81.6, y + 1.0, col = leg[[i]][1], border = if (i == 1) "#7FA7C9" else NA)
  text(82.6, y, leg[[i]][2], adj = 0, cex = 0.76, col = GREY)
}
# 2-5 constructors
cx <- list(c(2, 19), c(21, 38), c(40, 57), c(59, 76))
lab <- list(c("2", "si_stack()", "suitability rasters,", "one layer per species"),
            c("3", "si_metaweb()", "documented", "interactions"),
            c("4", "si_rule_*()", "phenology, traits,", "elevation, taxonomy"),
            c("5", "si_interaction()", "sign structure", "(RO / GloBI terms)"))
for (i in 1:4) wbox(cx[[i]][1], 68, cx[[i]][2], 80.5, PREP, lab[[i]][1], lab[[i]][2],
                   lab[[i]][3:4], hc = 0.86, sc = 0.76)
arr(30, 85.0, 11.5, 80.8); arr(48, 85.0, 30.5, 80.8); arr(65.2, 85.0, 49.5, 80.8)
# 6 constraint, 7 simulation
wbox(21, 51, 57, 63.5, CON, "6", "si_forbidden()",
    c("rules combined by product; documented", "records override or filter the rules"), hc = 0.9, sc = 0.76)
arr(29.5, 67.7, 33, 63.8); arr(48.5, 67.7, 45, 63.8); arr(67.5, 67.7, 55, 63.8, col = "#8C8C8C")
wbox(60, 51, 99, 63.5, SIMC, "7", "virtual systems",
    c("si_simulate_species(), si_simulate_climate(),", "si_simulate_metaweb(), si_phenology_simulate()"),
    hc = 0.9, sc = 0.72)
# 8 field, 9 diagnostic
rrect(2, 33, 60, 46.5, FLD)
text(31, 43.6, "8  si_overlap()  \u2192  si_field", font = 2, cex = 0.95, col = "white")
text(31, 39.6, "lazy spatial interaction field", cex = 0.8, col = "white")
text(31, 36.0, expression(italic(W)[italic(k)] == italic(f)(italic(p)[italic(k)] * "," ~ italic(a)[italic(k)]) ~ "\u2218" ~ bold(F)),
     cex = 0.9, col = "white")
arr(10.5, 67.7, 10.5, 46.8); arr(39, 50.7, 31, 46.8); arr(79.5, 50.7, 55, 46.8, col = "#8C8C8C")
wbox(66, 33, 99, 46.5, DIA, "9", "si_degenerate()",
    c("flags metrics that are", "constant by construction"), hc = 0.9, sc = 0.78)
arr(60.3, 39.8, 65.7, 39.8, col = DIA, lwd = 1.4)
# 10-13 modules
mx <- list(c(1, 25), c(26, 50), c(51, 75), c(76, 99.5))
ml <- list(c("10", "si_metrics()", "per-cell network metrics", "si_null(): effect sizes"),
           c("11", "si_beta()", "", "si_rewiring(), si_novelty()"),
           c("12", "si_robustness()", "R and Rw, sign-aware", "si_coextinction()"),
           c("13", "si_function()", "service, redundancy", "si_uncertainty()"))
for (i in 1:4) {
  wbox(mx[[i]][1], 14.5, mx[[i]][2], 27.5, MOD[i], ml[[i]][1], ml[[i]][2], ml[[i]][3:4],
      hc = 0.86, sc = 0.74)
  arr(31, 32.7, mean(mx[[i]]), 27.8, col = "#8C8C8C")
}
text(38, 22.4, expression(beta[WN] == beta[ST] + beta[OS]), cex = 0.8, col = "white")
# 14-15 outputs
wbox(1, 1, 49, 11.5, OUT, "14", "si_metric_map()  \u2192  SpatRaster", "maps of every response", hc = 0.88, sc = 0.78)
wbox(51, 1, 99.5, 11.5, OUT, "15", "si_experiment()  \u2192  data frame",
    "factorial designs; si_sensitivity()", hc = 0.88, sc = 0.78)
for (i in 1:2) arr(mean(mx[[i]]), 14.2, 25, 11.8, col = "#8C8C8C")
for (i in 3:4) arr(mean(mx[[i]]), 14.2, 75, 11.8, col = "#8C8C8C")
dev.off()

## ===========================================================================
## Figure 2 - computational design
## ===========================================================================
mem <- readRDS(file.path(DATA, "bench_mem.rds")); tim <- readRDS(file.path(DATA, "bench_time.rds"))
dev_open("Figure2.png", 6.5, 2.75)
base_par(mfrow = c(1, 3), mar = c(3.4, 3.6, 2.2, 0.6))
plot(mem$S, mem$dense_MB, type = "b", pch = 16, log = "xy", col = PAL[2], lwd = 1.8,
     xlab = "Species per level", ylab = "Memory (MB)",
     ylim = range(c(mem$lazy_MB, mem$dense_MB)) * c(0.8, 1.4))
lines(mem$S, mem$lazy_MB, type = "b", pch = 17, col = PAL[1], lwd = 1.8)
legend("topleft", bty = "n", cex = 0.85, lwd = 1.8, pch = c(16, 17), col = PAL[2:1],
       legend = c("One matrix per cell", "si_field"))
text(max(mem$S), mem$dense_MB[nrow(mem)] * 0.45, sprintf("%.0f\u00d7", max(mem$ratio)),
     col = PAL[2], font = 2, cex = 1, pos = 2)
title_panel("a", "Memory, 2,800 cells")
plot(tim$S, tim$full, type = "b", pch = 16, col = PAL[1], lwd = 1.8,
     xlab = "Species per level", ylab = "Seconds (1,000 cells)", ylim = c(0, max(tim$full) * 1.25))
lines(tim$S, tim$basic, type = "b", pch = 17, col = PAL[3], lwd = 1.8)
lines(tim$S, tim$beta,  type = "b", pch = 15, col = PAL[4], lwd = 1.8)
legend("topleft", bty = "n", cex = 0.82, lwd = 1.8, pch = c(16, 17, 15), col = PAL[c(1, 3, 4)],
       legend = c("All metrics", "Links, connectance", "\u03b2 partition"))
title_panel("b", "Runtime, one core")
par(mar = c(0.6, 0.6, 2.2, 0.6))
plot(NA, xlim = c(0, 10), ylim = c(0, 10), axes = FALSE, xlab = "", ylab = "")
rect(0.3, 6.2, 3.5, 9.6, col = PAL[1], border = NA); text(1.9, 7.9, "P\ncells \u00d7 nA", col = "white", cex = 0.9, font = 2)
rect(3.9, 6.2, 7.1, 9.6, col = PAL[1], border = NA); text(5.5, 7.9, "A\ncells \u00d7 nB", col = "white", cex = 0.9, font = 2)
rect(7.5, 6.2, 9.7, 9.6, col = PAL[2], border = NA); text(8.6, 7.9, "F\nnA \u00d7 nB", col = "white", cex = 0.9, font = 2)
text(5, 5.35, "stored", font = 3, cex = 0.9, col = GREY)
arrows(5, 4.9, 5, 3.9, length = 0.06, col = GREY, lwd = 1.3)
rect(1.2, 0.5, 8.8, 3.6, col = NA, border = PAL[3], lwd = 1.6, lty = 2)
text(5, 2.65, expression(italic(W)[italic(k)] == italic(f)(italic(p)[italic(k)] * "," ~ italic(a)[italic(k)]) ~ "\u2218" ~ bold(F)), cex = 1)
text(5, 1.3, "rebuilt on demand by si_local()", cex = 0.85, col = GREY)
title_panel("c", "Contents of an si_field")
dev.off()

## ===========================================================================
## Figure 3 - degeneracy diagnostic
## ===========================================================================
set.seed(7)
P <- matrix(runif(400 * 25), 400, 25, dimnames = list(NULL, paste0("A", 1:25)))
A <- matrix(runif(400 * 20), 400, 20, dimnames = list(NULL, paste0("B", 1:20)))
mwg <- si_simulate_metaweb(25, 20, 0.2, seed = 7, names_A = colnames(P), names_B = colnames(A))
f_bad  <- si_overlap(P, A, NULL, method = "product", check = FALSE)
f_good <- si_overlap(P, A, mwg, method = "product", floor = 0.05, check = FALSE)
m_bad  <- suppressWarnings(si_metrics(f_bad,  what = c("connectance", "NODF", "H2prime")))
m_good <- suppressWarnings(si_metrics(f_good, what = c("connectance", "NODF", "H2prime")))
mat_plot <- function(M, lab, txt) {
  M <- M[order(-rowSums(M)), order(-colSums(M))]
  image(t(M[nrow(M):1, ]), col = colorRampPalette(c("white", "#BBDEFB", "#1E88E5", "#0D47A1"))(64),
        axes = FALSE)
  box(col = LG); mtext("Species B", 1, line = 0.5, cex = 0.85, col = GREY)
  mtext("Species A", 2, line = 0.5, cex = 0.85, col = GREY, las = 0)
  title_panel(lab, txt)
}
dev_open("Figure3.png", 6.5, 4.6)
base_par(mfrow = c(2, 3), mar = c(2.2, 2.2, 2.2, 0.8))
mat_plot(si_local(f_bad, 1, drop = FALSE), "a", "Co-occurrence only")
mat_plot(mwg$matrix, "b", "Constraint layer F")
mat_plot(si_local(f_good, 1, drop = FALSE), "c", "Co-occurrence \u2218 F")
par(mar = c(3.4, 2.6, 2.2, 0.8))
labs <- c(connectance = "Connectance", NODF = "Nestedness (NODF)", H2prime = "H2\u2032")
for (v in names(labs)) {
  xb <- m_bad[[v]][is.finite(m_bad[[v]])]; xg <- m_good[[v]][is.finite(m_good[[v]])]
  rng <- range(c(xb, xg)); if (diff(rng) == 0) rng <- rng + c(-0.1, 0.1)
  plot(NA, xlim = rng, ylim = c(0, 1.12), yaxt = "n", xlab = labs[v], ylab = "")
  mtext("Scaled density", 2, line = 0.6, cex = 0.85, col = GREY, las = 0)
  d <- density(xg); polygon(d$x, d$y / max(d$y), col = adjustcolor(PAL[1], .25), border = PAL[1], lwd = 1.6)
  if (sd(xb) < 1e-12) segments(xb[1], 0, xb[1], 1, col = PAL[2], lwd = 3)
  else { d <- density(xb); polygon(d$x, d$y / max(d$y), col = adjustcolor(PAL[2], .25), border = PAL[2], lwd = 1.6) }
  if (v == "connectance")
    legend("top", bty = "n", cex = 0.82, lwd = c(3, 1.6), col = PAL[2:1], horiz = FALSE,
           legend = c("No constraint (SD = 0)", "With constraint"))
  title_panel(c(connectance = "d", NODF = "e", H2prime = "f")[[v]], "Across 400 cells")
}
dev.off()

## ===========================================================================
## Figure 4 - forbidden-link module
## ===========================================================================
set.seed(11); nA <- 30; nB <- 25
ph  <- si_phenology_simulate(nA, nB, seed = 11)
r_ph <- si_rule_phenology(ph, "current", threshold = 0.35)
r_mo <- si_rule_morphology(runif(nA, 5, 40), runif(nB, 3, 35), rule = "barrier")
elA <- runif(nA, 0, 2500); elB <- runif(nB, 0, 2500)
r_el <- si_rule_elevation(elA, elA + 900, elB, elB + 900)
fl   <- si_forbidden(phenology = r_ph, morphology = r_mo, elevation = r_el, binarise = TRUE)
part <- si_forbidden_partition(fl)
rule_plot <- function(M, lab, txt) {
  image(t(M[nrow(M):1, ]), col = c("#8C2D19", "#F4F4F4"), axes = FALSE, zlim = c(0, 1))
  box(col = LG); mtext("Species B", 1, line = 0.5, cex = 0.85, col = GREY)
  mtext("Species A", 2, line = 0.5, cex = 0.85, col = GREY, las = 0)
  title_panel(lab, txt)
}
dev_open("Figure4.png", 6.5, 4.7)
base_par(mfrow = c(2, 3), mar = c(2.2, 2.2, 2.4, 0.8))
rule_plot(r_ph$matrix, "a", sprintf("Phenology (%.0f%% forbidden)", 100 * mean(r_ph$matrix == 0)))
rule_plot(r_mo$matrix, "b", sprintf("Morphology (%.0f%%)", 100 * mean(r_mo$matrix == 0)))
rule_plot(r_el$matrix, "c", sprintf("Elevation (%.0f%%)", 100 * mean(r_el$matrix == 0)))
rule_plot(fl$matrix, "d", sprintf("Combined F (C = %.2f)", mean(fl$matrix > 0)))
par(mar = c(3.4, 3.4, 2.4, 0.8))
b <- barplot(part$forbids_pct, names.arg = c("Phenology", "Morphology", "Elevation"),
             col = PAL[c(1, 2, 3)], border = NA, ylab = "% of pairs forbidden",
             ylim = c(0, max(part$forbids_pct) * 1.25), cex.names = 0.9)
text(b, part$forbids_pct, labels = paste(part$unique_to_rule, "unique"), pos = 3, cex = 0.82, col = GREY)
title_panel("e", "Attribution by rule")
cs <- sapply(seq(0.1, 0.9, 0.1), function(t)
  mean(si_forbidden(phenology = si_rule_phenology(ph, "current", threshold = t),
                    morphology = r_mo, binarise = TRUE)$matrix > 0))
plot(seq(0.1, 0.9, 0.1), cs, type = "b", pch = 16, col = PAL[2], lwd = 1.8,
     xlab = "Phenological overlap threshold", ylab = "Connectance of F")
title_panel("f", "Sensitivity to threshold")
dev.off()

## ===========================================================================
## Figure 5 - sign structures and simulation engine
## ===========================================================================
sg <- readRDS(file.path(DATA, "signs.rds")); e1 <- readRDS(file.path(DATA, "exp1.rds"))
dev_open("Figure5.png", 6.5, 4.9)
base_par(mfrow = c(2, 3), mar = c(3.4, 3.4, 2.6, 0.8))
m <- rbind(sg$R_removeA, sg$R_removeB)
barplot(m, beside = TRUE, names.arg = sg$sign, col = c(PAL[1], PAL[3]), border = NA,
        ylab = "Robustness", ylim = c(0, 1.32), yaxt = "n", cex.names = 0.95)
axis(2, at = seq(0, 1, 0.2))
legend("top", horiz = TRUE, bty = "n", cex = 0.82, fill = c(PAL[1], PAL[3]), border = NA,
       legend = c("Remove A", "Remove B"), inset = c(0, -0.02))
title_panel("a", "Cascades by sign")
par(mar = c(0.6, 0.4, 2.6, 0.4))
plot(NA, xlim = c(0, 10), ylim = c(0, 7.4), axes = FALSE, xlab = "", ylab = "")
text(c(0.2, 1.7, 6.1, 8.3), 7.0, c("Sign", "Term", "Fate A", "Fate B"), adj = 0, font = 2, cex = 0.85)
fcol <- c(extinct = "#A93226", released = "#1E8449", neutral = GREY)
for (i in 1:6) {
  y <- 6.0 - (i - 1) * 1.05
  text(0.2, y, sg$sign[i], adj = 0, font = 2, cex = 0.95, col = PAL[1])
  text(1.7, y, sg$type[i], adj = 0, cex = 0.8, col = GREY)
  text(6.1, y, sg$fate_A[i], adj = 0, cex = 0.82, col = fcol[sg$fate_A[i]])
  text(8.3, y, sg$fate_B[i], adj = 0, cex = 0.82, col = fcol[sg$fate_B[i]])
}
title_panel("b", "si_cascade_rule()")
par(mar = c(2.4, 1.2, 2.6, 3.2))
st <- si_simulate_species(4, 40, 40, seed = 3)
im <- matrix(st$values[, 1], 40, 40)
image(im, col = hcl.colors(50, "YlOrRd", rev = TRUE), axes = FALSE, zlim = c(0, 1)); box(col = LG)
usr <- par("usr"); xs <- usr[2] + 0.04; par(xpd = NA)
yy <- seq(usr[3], usr[4], length.out = 51)
rect(xs, yy[-51], xs + 0.06, yy[-1], col = hcl.colors(50, "YlOrRd", rev = TRUE), border = NA)
text(xs + 0.08, c(usr[3], usr[4]), c("0", "1"), adj = 0, cex = 0.8)
text(xs + 0.03, usr[4] + 0.07, "Suitability", cex = 0.78); par(xpd = FALSE)
title_panel("c", "si_simulate_species()")
par(mar = c(5.2, 3.4, 2.6, 0.8))
scen <- c("decline", "trait_dependent", "contraction", "shift", "upslope", "increase")
v <- sapply(scen, function(s) mean(si_simulate_climate(st, s, severity = 0.35, distance = 6, seed = 1)$values))
cols <- c(PAL[1], PAL[1], PAL[1], PAL[2], PAL[2], PAL[3])
barplot(v, names.arg = c("Decline", "Trait-dep.", "Contract.", "Shift", "Upslope", "Increase"),
        col = cols, border = NA, las = 2, ylab = "Mean suitability", cex.names = 0.88,
        ylim = c(0, max(v, mean(st$values)) * 1.3))
abline(h = mean(st$values), lty = 2, col = GREY)
text(6.9, mean(st$values), "Present", pos = 3, cex = 0.78, col = GREY, offset = 0.2)
legend("topleft", bty = "n", cex = 0.78, fill = PAL[c(1, 2, 3)], border = NA,
       legend = c("Non-spatial decline", "Spatial displacement", "Increase"))
title_panel("d", "si_simulate_climate()")
par(mar = c(3.4, 3.4, 2.6, 0.8))
mm <- tapply(e1$d_connectance, e1$connectance, mean)
plot(as.numeric(names(mm)), mm, type = "b", pch = 16, col = PAL[1], lwd = 1.8,
     xlab = "Connectance of F", ylab = "\u0394 connectance")
abline(h = 0, lty = 3, col = LG)
title_panel("e", "si_experiment()")
sv <- si_sensitivity(e1, "beta_WN")
sv <- sv[sv$term != "Residuals", ][1:3, ]
par(mar = c(3.4, 7.6, 2.6, 0.8))
barplot(rev(sv$variance_explained), horiz = TRUE,
        names.arg = rev(c(severity = "Severity", connectance = "Connectance",
                          "connectance:severity" = "Interaction")[sv$term]),
        col = PAL[2], border = NA, las = 1, xlab = "Variance explained", cex.names = 0.9,
        xlim = c(0, 1.18))
text(rev(sv$variance_explained), seq(0.7, by = 1.2, length.out = 3),
     sprintf("%.3f", rev(sv$variance_explained)), pos = 4, cex = 0.8, col = GREY)
title_panel("f", "si_sensitivity()")
dev.off()

## ===========================================================================
## Figure 6 - worked example on the bundled Chilean subset
## ===========================================================================
data(chile_pp)
mk <- function(v) si_stack(v, grid = chile_pp$grid, cells = chile_pp$cells)
f1 <- si_overlap(mk(chile_pp$plants_current_bin), mk(chile_pp$pollinators_current_bin),
                 si_metaweb(chile_pp$metaweb), method = "binary", scenario = "current", check = FALSE)
f2 <- si_overlap(mk(chile_pp$plants_future_bin), mk(chile_pp$pollinators_future_bin),
                 si_metaweb(chile_pp$metaweb), method = "binary", scenario = "future", check = FALSE)
mc <- suppressWarnings(si_metrics(f1, what = c("links", "connectance")))
ex <- si_network_extinction(f1, f2)
bx <- ext(-76.2, -66.2, -56.2, -17.3)
grid_r <- si_grid_rast(chile_pp$grid)
lay <- function(vals) { r <- grid_r; r[chile_pp$cells] <- vals; crop(r, bx) }
r_links <- lay(mc$links)
r_conn  <- lay(ifelse(mc$links > 0, mc$connectance, NA))
r_chg   <- lay(ex$links_future - ex$links_current)
r_fate  <- lay(as.integer(ex$status))
r_land  <- lay(rep(1, length(chile_pp$cells)))
outl <- if (nzchar(OUTLINE) && file.exists(OUTLINE)) crop(vect(OUTLINE), bx) else NULL

classify_plot <- function(r, brks, cols, zero_col = "#E3E3E3") {
  plot(r_land, col = zero_col, legend = FALSE, axes = FALSE, mar = NA, box = FALSE)
  v <- values(r, mat = FALSE); k <- cut(v, brks, include.lowest = TRUE, labels = FALSE)
  rr <- r; values(rr) <- k
  plot(rr, col = cols, breaks = seq(0.5, length(cols) + 0.5, 1), legend = FALSE,
       axes = FALSE, add = TRUE)
  if (!is.null(outl)) lines(outl, col = "#6E6E6E", lwd = 0.35)
}
cbar <- function(cols, labs, title) {
  par(mar = c(1.6, 0.35, 1.1, 0.35))
  n <- length(cols)
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  if (length(labs) == n) {                       # categorical
    for (i in seq_len(n)) {
      yy <- 0.95 - (i - 1) * 0.24
      rect(0.02, yy - 0.16, 0.2, yy, col = cols[i], border = NA)
      text(0.25, yy - 0.08, labs[i], adj = 0, cex = 0.78)
    }
  } else {                                       # classed ramp
    x <- seq(0.08, 0.92, length.out = n + 1)
    rect(x[-(n + 1)], 0.45, x[-1], 0.85, col = cols, border = "white", lwd = 0.4)
    text(x, 0.26, labs, cex = 0.7, xpd = NA)
  }
  mtext(title, side = 3, line = 0.05, cex = 0.8, font = 2)
}

dev_open("Figure6.png", 6.5, 6.1, ps = 8.5)
layout(matrix(c(1, 3, 4, 5, 6,
                1, 3, 4, 5, 6,
                2, 3, 4, 5, 6,
                2, 7, 8, 9, 10), nrow = 4, byrow = TRUE),
       widths = c(1.75, 1, 1, 1, 1), heights = c(1, 1, 1, 0.52))
par(family = "sans", col.axis = GREY, fg = GREY)
# a metaweb
par(mar = c(2.2, 2.4, 2.2, 0.6))
M <- chile_pp$metaweb; M <- M[order(-rowSums(M)), order(-colSums(M))]
image(t(M[nrow(M):1, ]), col = c("white", PAL[1]), axes = FALSE); box(col = LG)
mtext("Pollinators (20)", 1, line = 0.6, cex = 0.8); mtext("Plants (20)", 2, line = 0.6, cex = 0.8, las = 0)
title_panel("a", sprintf("Metaweb, %d links", sum(chile_pp$metaweb > 0)), cex = 0.95)
# b richest local network
par(mar = c(1.2, 0.6, 2.2, 0.6))
si_plot_network(f1, cell = which.max(mc$links), labels = FALSE)
title_panel("b", sprintf("Richest cell, %d links", max(mc$links)), cex = 0.95)
# maps
lk_br <- c(0.5, 2.5, 5.5, 10.5, 20.5, 40.5, Inf)
lk_col <- hcl.colors(6, "viridis")
par(mar = c(0.2, 0.2, 2.2, 0.2)); classify_plot(r_links, lk_br, lk_col); title_panel("c", "Links", cex = 0.95)
cn_br <- c(0, 0.4, 0.5, 0.6, 0.7, 0.85, 1.0001)
cn_col <- hcl.colors(6, "Plasma")
par(mar = c(0.2, 0.2, 2.2, 0.2)); classify_plot(r_conn, cn_br, cn_col); title_panel("d", "Connectance", cex = 0.95)
ch_br <- c(-Inf, -20.5, -10.5, -3.5, -0.5, 0.5, 3.5, 10.5, Inf)
ch_col <- c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#E3E3E3", "#D1E5F0", "#92C5DE", "#2166AC")
par(mar = c(0.2, 0.2, 2.2, 0.2)); classify_plot(r_chg, ch_br, ch_col); title_panel("e", "Change", cex = 0.95)
fate_col <- c("#2C7FB8", "#D7301F", "#31A354", "#E3E3E3")
par(mar = c(0.2, 0.2, 2.2, 0.2)); classify_plot(r_fate, c(0.5, 1.5, 2.5, 3.5, 4.5), fate_col)
title_panel("f", "Fate by 2070", cex = 0.95)
# colour bars
cbar(lk_col, c("1", "3", "6", "11", "21", "41", ""), "links per cell")
cbar(cn_col, c("0.25", "0.4", "0.5", "0.6", "0.7", "0.85", "1"), "connectance")
cbar(ch_col[c(1:4, 6:8)], c("", "-20", "-10", "-3", "+3", "+10", "", ""), "change in links")
cbar(fate_col, c("persists", "lost", "gained", "none"), "network fate")
dev.off()
cat("figures written:", paste(list.files(FIG, pattern = "^Figure"), collapse = ", "), "\n")
