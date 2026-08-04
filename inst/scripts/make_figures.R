## ---------------------------------------------------------------------------
## spinet paper: figures
## ---------------------------------------------------------------------------
suppressMessages({library(spinet); library(terra)})
## Output and input directories. Override with environment variables or edit.
FIG <- Sys.getenv("SPINET_FIG", unset = "figures")
OUT <- Sys.getenv("SPINET_OUT", unset = "spinet_output")
SIM <- Sys.getenv("SPINET_SIM", unset = "simulation_output")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
for (d in c(OUT, SIM)) if (!dir.exists(d))
  stop("Directory not found: '", d, "'. Run chile_analysis.R and the ",
       "simulation script first, or set SPINET_OUT / SPINET_SIM.")

## house style ---------------------------------------------------------------
PAL   <- c("#1B3A5C", "#C1553B", "#3F8F6E", "#B8912F", "#7A5C8E", "#5A6B75")
GREY  <- "#4A4A4A"; LGREY <- "#D8D8D8"
pdev <- function(f, w, h) grDevices::png(file.path(FIG, f), width = w, height = h,
                                         units = "in", res = 300, type = "cairo")
par_base <- function(...) par(family = "sans", mgp = c(2.1, 0.6, 0), tcl = -0.25,
                             cex.axis = 0.82, cex.lab = 0.95, las = 1,
                             bty = "n", col.axis = GREY, fg = GREY, ...)
panel <- function(l) mtext(l, side = 3, adj = -0.13, line = 0.9, font = 2,
                           cex = 0.95, col = "black")
sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
band <- function(x, m, s, col) polygon(c(x, rev(x)), c(m - s, rev(m + s)),
                                       col = adjustcolor(col, 0.18), border = NA)

## ===========================================================================
## FIGURE 1 -- the degeneracy, on real data
## ===========================================================================
AB <- readRDS(file.path(OUT, "AB.rds")); CD <- readRDS(file.path(OUT, "CD.rds"))$res
lab <- c("A\nbinary\nno metaweb", "B\ncontinuous\nno metaweb",
         "C\nbinary\n\u00d7 metaweb", "D\ncontinuous\n\u00d7 metaweb")
get <- function(p, s, v) { z <- if (p %in% c("A","B")) AB[[p]] else CD[[p]]; z[[s]][[v]] }

pdev("Fig1_degeneracy.png", 9.2, 6.6)
par_base(mfrow = c(2, 3), mar = c(4.4, 4.2, 2.6, 1))
for (v in c("connectance", "NODF", "H2_gap")) {
  ms <- sapply(c("A","B","C","D"), function(p) mean(get(p,"m1",v), na.rm = TRUE))
  ss <- sapply(c("A","B","C","D"), function(p) sd(get(p,"m1",v), na.rm = TRUE))
  ylab <- switch(v, connectance = "connectance", NODF = "nestedness (NODF)",
                 H2_gap = expression(H[2]^{"max"} - H[2]))
  b <- barplot(ms, names.arg = lab, ylab = ylab, col = c(LGREY, LGREY, PAL[1], PAL[2]),
               border = NA, ylim = c(0, max(ms + ss, 0.05) * 1.25), cex.names = 0.62)
  arrows(b, ms - ss, b, ms + ss, code = 3, angle = 90, length = 0.03, col = GREY, lwd = 1)
  for (i in 1:2) text(b[i], max(ms + ss) * 0.12,
                      sprintf("%.3g\nsd = 0", ms[i]), cex = 0.62, col = PAL[2], font = 2)
  panel(c(connectance = "a", NODF = "b", H2_gap = "c")[v])
}
## per-cell distributions
for (v in c("connectance", "NODF", "H2_gap")) {
  xs <- lapply(c("A","C","D"), function(p) get(p, "m1", v))
  xs <- lapply(xs, function(z) z[is.finite(z)])
  rng <- range(unlist(xs)); if (diff(rng) == 0) rng <- rng + c(-0.05, 0.05)
  plot(NA, xlim = rng, ylim = c(0, 1), xlab = switch(v, connectance="connectance",
       NODF="nestedness (NODF)", H2_gap=expression(H[2]^{"max"}-H[2])),
       ylab = "scaled density", yaxt = "n")
  cols <- c(GREY, PAL[1], PAL[2])
  for (i in seq_along(xs)) {
    if (sd(xs[[i]]) < 1e-12) {
      segments(xs[[i]][1], 0, xs[[i]][1], 1, col = cols[i], lwd = 3)
    } else {
      d <- density(xs[[i]], na.rm = TRUE); lines(d$x, d$y/max(d$y), col = cols[i], lwd = 2)
      polygon(d$x, d$y/max(d$y), col = adjustcolor(cols[i], 0.16), border = NA)
    }
  }
  if (v == "connectance") legend("topright", legend = c("A / B  (no metaweb)",
      "C  (binary \u00d7 metaweb)", "D  (continuous \u00d7 metaweb)"),
      col = cols, lwd = 2, bty = "n", cex = 0.66)
  panel(c(connectance = "d", NODF = "e", H2_gap = "f")[v])
}
dev.off()

## ===========================================================================
## FIGURE 2 -- Chile maps
## ===========================================================================
mp <- rast(file.path(OUT, "chile_metrics_current.tif"))
st <- rast(file.path(OUT, "chile_status.tif"))
ex <- CD$C$ex
tm <- readRDS(file.path(OUT, "CD.rds"))
chg <- rast(mp[[1]]); values(chg) <- NA
cells <- readRDS(file.path(OUT, "stacks.rds"))$S$Pc_c$cells
chg[cells] <- ex$links_future - ex$links_current

pdev("Fig2_chile_maps.png", 10.2, 6.2)
par_base(mfrow = c(1, 4), mar = c(2, 1.4, 3, 4.2))
plot(mp[["links"]], col = si_palette_metric(60), axes = FALSE, mar = c(2,1.4,3,4.2),
     main = "", plg = list(title = "links", cex = 0.7)); panel("a  local links, present")
plot(mp[["connectance"]], col = si_palette_metric(60), axes = FALSE, mar = c(2,1.4,3,4.2),
     main = "", plg = list(title = "C", cex = 0.7)); panel("b  connectance")
plot(chg, col = si_palette_change(60), axes = FALSE, mar = c(2,1.4,3,4.2),
     main = "", plg = list(title = "\u0394 links", cex = 0.7)); panel("c  change by 2070")
stf <- as.factor(st)
plot(st, col = c("#BDBDBD", "#8C2D19", "#2E6E4E", "#F2F2F2"), axes = FALSE,
     mar = c(2,1.4,3,4.2), main = "", type = "classes",
     levels = c("persists","extinct","colonised","no network"),
     plg = list(cex = 0.65)); panel("d  network fate")
dev.off()

## ===========================================================================
## FIGURE 3 -- constraint density governs the topological response
## ===========================================================================
e1 <- readRDS(file.path(SIM, "exp1.rds"))
pdev("Fig3_constraint_density.png", 9.6, 3.5)
par_base(mfrow = c(1, 3), mar = c(4.2, 4.4, 2.6, 1))
sev <- sort(unique(e1$severity))
for (v in c("d_connectance", "d_H2prime", "beta_WN")) {
  ylab <- switch(v, d_connectance = expression(Delta*" connectance"),
                 d_H2prime = expression(Delta*H[2]*"'"), beta_WN = expression(beta[WN]))
  yl <- range(tapply(e1[[v]], list(e1$connectance, e1$severity), mean), na.rm = TRUE)
  plot(NA, xlim = rev(range(e1$connectance)), ylim = yl,
       xlab = "constraint layer connectance", ylab = ylab)
  abline(h = 0, col = LGREY, lwd = 1)
  for (i in seq_along(sev)) {
    s <- e1[e1$severity == sev[i], ]
    m <- tapply(s[[v]], s$connectance, mean, na.rm = TRUE)
    e <- tapply(s[[v]], s$connectance, sem)
    x <- as.numeric(names(m)); o <- order(x)
    band(x[o], m[o], e[o], PAL[i]); lines(x[o], m[o], col = PAL[i], lwd = 2.1)
    points(x[o], m[o], col = PAL[i], pch = 16, cex = 0.8)
  }
  if (v == "d_connectance") legend("topright", legend = paste("severity", sev),
      col = PAL[seq_along(sev)], lwd = 2, bty = "n", cex = 0.7)
  if (v != "beta_WN") text(1, 0, "0 exactly", pos = 3, cex = 0.68, col = PAL[2], font = 2)
  panel(c(d_connectance="a", d_H2prime="b", beta_WN="c")[v])
}
dev.off()

## ===========================================================================
## FIGURE 4 -- rewiring: mechanism and its irrelevance to persistence
## ===========================================================================
e2 <- readRDS(file.path(SIM, "exp2.rds"))
pdev("Fig4_rewiring.png", 9.6, 6.6)
par_base(mfrow = c(2, 2), mar = c(4.3, 4.6, 2.6, 1))
sd_ <- sort(unique(e2$shift_sd)); sev <- sort(unique(e2$severity))

## a: beta_OS vs shift_sd
m <- tapply(e2$beta_OS, e2$shift_sd, mean); e <- tapply(e2$beta_OS, e2$shift_sd, sem)
plot(sd_, m, type = "n", xlab = "among-species SD of phenological shift (days)",
     ylab = expression("interaction rewiring  "*beta[OS]), ylim = c(0, max(m+e)*1.1))
band(sd_, m, e, PAL[1]); lines(sd_, m, col = PAL[1], lwd = 2.4); points(sd_, m, pch = 16, col = PAL[1])
points(0, m[1], pch = 21, bg = "white", col = PAL[2], cex = 1.6, lwd = 2)
text(0, m[1], "  0 exactly", pos = 4, cex = 0.72, col = PAL[2], font = 2)
panel("a")

## b: partition of turnover
ms <- t(sapply(sd_, function(s) c(ST = mean(e2$beta_ST[e2$shift_sd == s], na.rm = TRUE),
                                  OS = mean(e2$beta_OS[e2$shift_sd == s], na.rm = TRUE))))
b <- barplot(t(ms), beside = FALSE, names.arg = sd_, col = c(PAL[6], PAL[2]), border = NA,
             xlab = "among-species SD of phenological shift (days)",
             ylab = expression(beta[WN]*"  =  "*beta[ST]*"  +  "*beta[OS]))
legend("topleft", legend = c(expression(beta[ST]*"  species turnover"),
                             expression(beta[OS]*"  rewiring")),
       fill = c(PAL[6], PAL[2]), border = NA, bty = "n", cex = 0.75)
panel("b")

## c: rewiring gains
m <- tapply(e2$rewiring_gains, e2$shift_sd, mean); e <- tapply(e2$rewiring_gains, e2$shift_sd, sem)
plot(sd_, m, type = "n", xlab = "among-species SD of phenological shift (days)",
     ylab = "new links among shared species", ylim = c(0, max(m+e)*1.12))
band(sd_, m, e, PAL[3]); lines(sd_, m, col = PAL[3], lwd = 2.4); points(sd_, m, pch = 16, col = PAL[3])
panel("c")

## d: extinctions do NOT respond to rewiring
plot(NA, xlim = range(sd_), ylim = range(tapply(e2$n_extinct, list(e2$shift_sd, e2$severity), mean)),
     xlab = "among-species SD of phenological shift (days)", ylab = "networks lost (cells)")
for (i in seq_along(sev)) {
  s <- e2[e2$severity == sev[i], ]
  m <- tapply(s$n_extinct, s$shift_sd, mean); e <- tapply(s$n_extinct, s$shift_sd, sem)
  band(sd_, m, e, PAL[i]); lines(sd_, m, col = PAL[i], lwd = 2.1); points(sd_, m, pch = 16, col = PAL[i], cex = 0.8)
  ft <- lm(n_extinct ~ shift_sd, data = s)
  text(max(sd_), m[length(m)], sprintf("  p = %.2f", summary(ft)$coefficients[2,4]),
       pos = 4, cex = 0.62, col = PAL[i], xpd = NA)
}
legend("topleft", legend = paste("severity", sev), col = PAL[seq_along(sev)],
       lwd = 2, bty = "n", cex = 0.72)
panel("d")
dev.off()

## ===========================================================================
## FIGURE 5 -- architecture, richness, scenarios
## ===========================================================================
e3 <- readRDS(file.path(SIM, "exp3.rds")); e4 <- readRDS(file.path(SIM, "exp4.rds"))
e5 <- readRDS(file.path(SIM, "exp5.rds"))
pdev("Fig5_architecture_richness_scenarios.png", 9.8, 6.6)
par_base(mfrow = c(2, 2), mar = c(4.4, 4.6, 2.6, 1))

## a architecture
ar <- c("random","nested","modular","heterogeneous")
m <- sapply(ar, function(a) mean(e3$n_extinct[e3$architecture == a]))
e <- sapply(ar, function(a) sem(e3$n_extinct[e3$architecture == a]))
b <- barplot(m, names.arg = ar, col = PAL[1:4], border = NA, ylab = "networks lost (cells)",
             ylim = c(0, max(m+e)*1.2), cex.names = 0.78)
arrows(b, m-e, b, m+e, code = 3, angle = 90, length = 0.03, col = GREY)
kt <- kruskal.test(n_extinct ~ factor(architecture), data = e3)
mtext(sprintf("Kruskal-Wallis  chi2 = %.1f,  p = %.3f", kt$statistic, kt$p.value),
      side = 3, line = -0.4, cex = 0.66, col = GREY)
panel("a  metaweb architecture (connectance matched)")

## b richness scaling
rr <- sort(unique(e4$richness))
abs_ <- tapply(e4$cur_links - e4$fut_links, e4$richness, mean)
rel_ <- tapply((e4$cur_links - e4$fut_links)/pmax(e4$cur_links,1e-9), e4$richness, mean)
par(mar = c(4.4, 4.6, 2.6, 4.4))
plot(rr, abs_, type = "b", pch = 16, col = PAL[1], lwd = 2, log = "xy",
     xlab = "species per level", ylab = "absolute link loss")
par(new = TRUE)
plot(rr, rel_, type = "b", pch = 17, col = PAL[2], lwd = 2, log = "x",
     axes = FALSE, xlab = "", ylab = "", ylim = c(0.5, 0.85))
axis(4, col = PAL[2], col.axis = PAL[2]); mtext("relative link loss", 4, line = 2.5, col = PAL[2], las = 0, cex = 0.9)
legend("bottomright", legend = c("absolute (left)", "relative (right)"),
       col = PAL[1:2], pch = c(16,17), lwd = 2, bty = "n", cex = 0.72)
panel("b  richness scaling")

## c networks lost vs richness
par(mar = c(4.4, 4.6, 2.6, 1))
m <- tapply(e4$n_extinct, e4$richness, mean); e <- tapply(e4$n_extinct, e4$richness, sem)
plot(rr, m, type = "n", xlab = "species per level", ylab = "networks lost (cells)", log = "x",
     ylim = c(0, max(m+e)*1.1))
band(rr, m, e, PAL[3]); lines(rr, m, col = PAL[3], lwd = 2.4); points(rr, m, pch = 16, col = PAL[3])
panel("c  small networks are the vulnerable ones")

## d scenarios
sc <- c("increase","trait_dependent","contraction","mixed","random","decline","shift","upslope")
m <- sapply(sc, function(s) mean(e5$n_extinct[e5$scenario == s]))
e <- sapply(sc, function(s) sem(e5$n_extinct[e5$scenario == s]))
cols <- ifelse(sc %in% c("shift","upslope"), PAL[2], PAL[6])
b <- barplot(m, names.arg = sc, col = cols, border = NA, ylab = "networks lost (cells)",
             las = 2, cex.names = 0.68, ylim = c(0, max(m+e)*1.2))
arrows(b, m-e, b, m+e, code = 3, angle = 90, length = 0.03, col = GREY)
legend("topleft", legend = c("non-spatial", "spatially explicit"), fill = c(PAL[6], PAL[2]),
       border = NA, bty = "n", cex = 0.72)
panel("d  climate scenario")
dev.off()

## ===========================================================================
## FIGURE 6 -- robustness: constrained vs unconstrained rewiring
## ===========================================================================
set.seed(99)
archs <- c("random","nested","modular","heterogeneous")
rb <- do.call(rbind, lapply(archs, function(a) {
  do.call(rbind, lapply(1:12, function(i) {
    mw <- si_simulate_metaweb(30, 25, 0.2, architecture = a, seed = i)
    w  <- mw$matrix * matrix(runif(750, 1, 10), 30, 25)
    tr <- si_rule_morphology(runif(30, 5, 40), runif(25, 3, 35), rule = "gower")$matrix
    r0 <- si_robustness(w, remove = "A", n_sim = 20, seed = i)
    r1 <- si_robustness(w, remove = "A", rewiring = matrix(1, 30, 25), n_sim = 20, seed = i)
    r2 <- si_robustness(w, remove = "A", rewiring = tr, n_sim = 20, seed = i)
    data.frame(architecture = a, R = r0$robustness, Rw_free = r1$robustness,
               Rw_constrained = r2$robustness)
  }))
}))
saveRDS(rb, file.path(SIM, "robustness.rds"))

pdev("Fig6_robustness.png", 9.4, 3.6)
par_base(mfrow = c(1, 3), mar = c(4.4, 4.6, 2.6, 1))
## a curves
mw <- si_simulate_metaweb(30, 25, 0.2, architecture = "nested", seed = 1)
w  <- mw$matrix * matrix(runif(750, 1, 10), 30, 25)
tr <- si_rule_morphology(runif(30, 5, 40), runif(25, 3, 35), rule = "gower")$matrix
cs <- list(si_robustness(w, remove="A", n_sim=60, seed=1),
           si_robustness(w, remove="A", rewiring=tr, n_sim=60, seed=1),
           si_robustness(w, remove="A", rewiring=matrix(1,30,25), n_sim=60, seed=1))
plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "fraction of plants removed",
     ylab = "fraction of pollinators surviving")
abline(0, -1, lty = 3, col = LGREY)
for (i in 1:3) lines(seq(0,1,length.out = length(cs[[i]]$curve)), cs[[i]]$curve,
                     col = PAL[c(6,1,3)][i], lwd = 2.4)
legend("bottomleft", bty = "n", cex = 0.72, lwd = 2.4, col = PAL[c(6,1,3)],
       legend = sprintf(c("R = %.3f  (no rewiring)","Rw = %.3f  (trait-constrained)",
                          "Rw = %.3f  (unconstrained)"), sapply(cs, `[[`, "robustness")))
panel("a")
## b effect sizes
d1 <- rb$Rw_constrained - rb$R; d2 <- rb$Rw_free - rb$R
bp <- boxplot(list(`trait-constrained` = d1, unconstrained = d2), col = c(PAL[1], PAL[3]),
              border = GREY, ylab = expression("rewiring effect size  D = "*R[w]*" - R"),
              cex.axis = 0.78)
tt <- t.test(d2, d1, paired = TRUE)
mtext(sprintf("paired t = %.1f,  p = %.1e,  inflation = %.0f%%", tt$statistic, tt$p.value,
              100*(mean(d2)/mean(d1) - 1)), side = 3, line = -0.4, cex = 0.62, col = GREY)
panel("b")
## c by architecture
m1 <- tapply(rb$Rw_constrained - rb$R, rb$architecture, mean)[archs]
m2 <- tapply(rb$Rw_free - rb$R, rb$architecture, mean)[archs]
b <- barplot(rbind(m1, m2), beside = TRUE, names.arg = archs, col = c(PAL[1], PAL[3]),
             border = NA, ylab = "rewiring effect size D", cex.names = 0.7, las = 2)
legend("topleft", legend = c("trait-constrained","unconstrained"), fill = c(PAL[1], PAL[3]),
       border = NA, bty = "n", cex = 0.7)
panel("c")
dev.off()

cat("figures written to", FIG, "\n"); print(list.files(FIG))
