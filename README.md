# spinet

**Spatial interaction networks under environmental change**

`spinet` projects potential ecological interaction networks through space and
time. It combines species distribution model (SDM) output for two groups of
species, such as plants and pollinators, with a constraint layer that records
which pairs can interact, and computes network structure in every grid cell for
present and future climates. It works downstream of SDM software such as
flexsdm or biomod2 and does not fit distribution models itself.

Networks produced by `spinet` describe **potential** interactions: two species
are linked in a cell when both are predicted to occur there and the pair is
permitted by the constraint layer.

## Installation

`spinet` needs R (>= 4.1) and the `terra` package. It contains no compiled
code, so Rtools is not required on Windows.

```r
# From GitHub, latest version
install.packages("remotes")
remotes::install_github("ehsanrahimi666/spinet")

# A specific release, e.g. the version described in the paper
remotes::install_github("ehsanrahimi666/spinet@v0.1.2")

# From a downloaded source archive
install.packages("spinet_0.1.2.tar.gz", repos = NULL, type = "source")

# Optional packages used by some functions
install.packages(c("igraph", "bipartite", "sf", "vegan"))
```

Check that it works:

```r
library(spinet)
packageVersion("spinet")
```

On Linux, `terra` needs the GDAL, GEOS and PROJ system libraries; see
<https://github.com/rspatial/terra> if its installation fails.

## Quick start (runs as-is, no files needed)

The package includes a 20 × 20 subset of a Chilean plant–pollinator system,
with present-day and 2070 suitability surfaces and the documented interactions.

```r
library(spinet)
data(chile_pp)

# 1. Suitability surfaces (binary presence), one stack per level
plants <- si_stack(chile_pp$plants_current_bin, level = "plant",
                   grid = chile_pp$grid, cells = chile_pp$cells)
polls  <- si_stack(chile_pp$pollinators_current_bin, level = "pollinator",
                   grid = chile_pp$grid, cells = chile_pp$cells)

# 2. Constraint layer from documented interactions
mw <- si_metaweb(chile_pp$metaweb)
mw                                   # 114 links, connectance 0.285

# 3. Spatial interaction field and per-cell network metrics
f_now <- si_overlap(plants, polls, mw, method = "binary", scenario = "present")
m_now <- si_metrics(f_now, what = c("links", "connectance", "NODF"))
summary(m_now)

# 4. The same for 2070, then change between periods
plants_f <- si_stack(chile_pp$plants_future_bin, level = "plant",
                     grid = chile_pp$grid, cells = chile_pp$cells)
polls_f  <- si_stack(chile_pp$pollinators_future_bin, level = "pollinator",
                     grid = chile_pp$grid, cells = chile_pp$cells)
f_fut <- si_overlap(plants_f, polls_f, mw, method = "binary", scenario = "2070")

b <- si_beta(f_now, f_fut)           # interaction turnover per cell
round(colMeans(b[, c("beta_WN", "beta_ST", "beta_OS")], na.rm = TRUE), 3)
table(si_network_extinction(f_now, f_fut)$status)
```

Output of the last two lines:

```
beta_WN beta_ST beta_OS
  0.584   0.584   0.000

 persists   extinct colonised    absent
      485       110        23      2182
```

`si_beta()` also warns that rewiring (`beta_OS`) is zero by construction here,
because both periods use the same documented metaweb: shared species keep the
same potential partners. Rewiring can only be estimated when the constraint
layer changes between periods, for example through projected phenology.

## Example 1: from GeoTIFF files to maps

The package also ships a 20 × 20 subset as GeoTIFF files, to show the raster
path. This is the workflow you would use with your own SDM output.

```r
library(terra)
si_example_rasters("all")            # file names of the example data

region <- vect(si_example_rasters("region"))   # study area, used as mask
P_now <- si_stack(si_example_rasters("plants_current_binary"),      mask = region, level = "plant")
A_now <- si_stack(si_example_rasters("pollinators_current_binary"), mask = region, level = "pollinator")
P_fut <- si_stack(si_example_rasters("plants_future_binary"),       mask = region, level = "plant")
A_fut <- si_stack(si_example_rasters("pollinators_future_binary"),  mask = region, level = "pollinator")

web <- as.matrix(read.csv(si_example_rasters("metaweb"), row.names = 1))
mw2 <- si_metaweb(web, names_A = P_now$names, names_B = A_now$names)

f1 <- si_overlap(P_now, A_now, mw2, method = "binary", scenario = "present")
f2 <- si_overlap(P_fut, A_fut, mw2, method = "binary", scenario = "2070")
m1 <- si_metrics(f1, what = c("links", "connectance", "NODF"))
m2 <- si_metrics(f2, what = c("links", "connectance", "NODF"))

# Maps of network metrics
maps <- si_metric_map(m1, c("links", "connectance", "NODF"))
plot(maps)

# Change in the number of potential links
change <- si_metric_map(m2, "links") - si_metric_map(m1, "links")
plot(change, main = "Change in potential links, present to 2070")

# The local network of the richest cell
si_plot_network(f1, cell = which.max(m1$links))

# Save the maps as a multi-layer GeoTIFF
writeRaster(maps, "spinet_metrics_present.tif", overwrite = TRUE)
```

## Example 2: continuous suitability

With continuous suitability, each potential link is weighted by the product of
the two species' suitabilities. Because correlative suitability is rarely
exactly zero, set a floor with `si_prune()` so that presence-based metrics
remain informative.

```r
pc <- si_stack(chile_pp$plants_current, level = "plant",
               grid = chile_pp$grid, cells = chile_pp$cells)
ac <- si_stack(chile_pp$pollinators_current, level = "pollinator",
               grid = chile_pp$grid, cells = chile_pp$cells)

f_cont <- si_prune(si_overlap(pc, ac, mw, method = "product", check = FALSE),
                   floor = 0.05)
summary(si_metrics(f_cont, what = c("links", "connectance", "H2prime")))
```

## Example 3: forbidden links from traits and phenology

Biological rules estimate which links are impossible without relying on
interaction records. Here, phenological overlap and a morphological barrier
(proboscis shorter than corolla) are combined, then documented interactions
are added.

```r
set.seed(1)
nA <- 20; nB <- 15
ph        <- si_phenology_simulate(nA, nB, seed = 1)   # activity periods
corolla   <- runif(nA, 5, 40)                          # corolla depth, mm
proboscis <- runif(nB, 3, 35)                          # proboscis length, mm

rules <- list(
  phenology  = si_rule_phenology(ph, threshold = 0.3),
  morphology = si_rule_morphology(corolla, proboscis, rule = "barrier"))

fl <- do.call(si_forbidden, c(rules, binarise = TRUE))
fl                             # permitted pairs and share forbidden by each rule
si_forbidden_partition(fl)     # pairs forbidden by each rule alone

# Documented interactions: 10 pairs the rules permit, 4 they forbid
M   <- fl$metaweb$matrix
obs <- matrix(0, nA, nB)
obs[sample(which(M > 0), 10)] <- 1
obs[sample(which(M == 0), 4)] <- 1

# "override": documented pairs are always permitted (default)
# "filter":   only documented pairs that the rules also permit are kept
over <- do.call(si_forbidden, c(rules, list(observed = obs, binarise = TRUE,
                                            observed_action = "override")))
filt <- do.call(si_forbidden, c(rules, list(observed = obs, binarise = TRUE,
                                            observed_action = "filter")))
c(rules_only = si_connectance(fl), override = si_connectance(over),
  filter = si_connectance(filt))

# Use the rule-based layer with simulated suitability surfaces
P <- si_simulate_species(nA, nrow = 30, ncol = 30, prefix = "plant", seed = 1)
A <- si_simulate_species(nB, nrow = 30, ncol = 30, prefix = "bee",   seed = 2)
f_rules <- si_prune(si_overlap(P, A, over, method = "product", check = FALSE),
                    floor = 0.05)
summary(si_metrics(f_rules, what = c("links", "connectance", "H2prime")))
```

## Example 4: the degeneracy check

Without a constraint layer, every plant in a cell is linked to every pollinator
there. Each local network is then complete, and its metrics describe species
richness only. `si_overlap()` detects this and warns before any metric is
computed:

```r
f_bad <- si_overlap(plants, polls, NULL, method = "binary")
#> Warning: si_field is structurally degenerate:
#>   * complete-bipartite -> all binary metrics
```

## Example 5: robustness and interaction sign

```r
web <- si_simulate_metaweb(25, 20, connectance = 0.2,
                           architecture = "nested", seed = 11)
si_robustness(web$matrix, remove = "A", sequence = "random",
              n_sim = 100, seed = 1)
#> robustness : 0.7571   (sd 0.0645, n = 100)

# How partner loss propagates depends on the interaction type
si_cascade_rule(si_interaction("pollinatedBy"))   # mutualism: both levels coextinct
si_cascade_rule(si_interaction("eatenBy"))         # antagonism: resources released
```

## Example 6: virtual species and climate scenarios

```r
st <- si_simulate_species(10, nrow = 30, ncol = 30, seed = 5)
round(sapply(c("decline", "shift", "upslope"), function(s)
  mean(si_simulate_climate(st, s, severity = 0.3, distance = 5,
                           seed = 1)$values)), 3)
#> decline   shift upslope
#>   0.170   0.204   0.143
```

`si_experiment()` runs factorial designs over such scenarios and constraint
layers, and `si_sensitivity()` partitions the variance among design factors.

## Using your own species distribution models

Organise one folder per network level and period, with one raster per species
named after the species. All rasters must share one grid (extent, resolution
and coordinate reference system).

```
my_project/
  sdm/plants_current/        Adesmia_conferta.tif, Alstroemeria_ligtu.tif, ...
  sdm/plants_2070/
  sdm/pollinators_current/
  sdm/pollinators_2070/
  data/interactions.csv      plant names in column 1, pollinator names as header
  gis/study_area.shp
```

```r
# Template: replace the paths with your own
library(spinet); library(terra)
area <- vect("gis/study_area.shp")

# band 1 = continuous suitability, band 2 = binary presence (flexsdm/MaxEnt
# convention). Use the same mask for every stack.
P_now <- si_stack("sdm/plants_current",      band = 2, mask = area, level = "plant")
A_now <- si_stack("sdm/pollinators_current", band = 2, mask = area, level = "pollinator")
P_fut <- si_stack("sdm/plants_2070",         band = 2, mask = area, level = "plant")
A_fut <- si_stack("sdm/pollinators_2070",    band = 2, mask = area, level = "pollinator")

web <- as.matrix(read.csv("data/interactions.csv", row.names = 1, check.names = FALSE))
mw  <- si_metaweb(web, names_A = P_now$names, names_B = A_now$names)

f_now <- si_overlap(P_now, A_now, mw, method = "binary", scenario = "present")
f_fut <- si_overlap(P_fut, A_fut, mw, method = "binary", scenario = "2070")

what <- c("links", "connectance", "NODF")
writeRaster(si_metric_map(si_metrics(f_now, what = what), what),
            "metrics_present.tif", overwrite = TRUE)
writeRaster(si_metric_map(si_metrics(f_fut, what = what), what),
            "metrics_2070.tif", overwrite = TRUE)
write.csv(si_beta(f_now, f_fut), "interaction_turnover.csv", row.names = FALSE)
```

Species names are matched after normalisation, so `Bombus terrestris`,
`Bombus_terrestris` and `Bombus.terrestris` are treated as the same species.
`si_metaweb()` warns if species in the rasters are missing from the interaction
table. The complete script for a 187 × 171 system is at
`system.file("scripts", "chile_case_study.R", package = "spinet")`.

## Interaction sign structures

The interaction vocabulary follows the Relations Ontology terms used by Global
Biotic Interactions (GloBI), so GloBI-standard tables can be read directly.

```r
si_interaction("pollinatedBy")   # ++  mutualism      both coextinct
si_interaction("eatenBy")        # +-  antagonism     consumer dies, resource released
si_interaction("+0")             # +0  commensalism   one-way cascade
si_interaction("-0")             # -0  amensalism     target released
si_interaction("competesWith")   # --  competition    unipartite, overlap_sign = -1
si_interaction("00")             # 00  neutralism     null benchmark
```

## Main functions

| Module | Functions |
|---|---|
| Objects | `si_stack()` `si_interaction()` `si_metaweb()` `si_grid()` |
| Engine | `si_overlap()` `si_local()` `si_threshold()` `si_prune()` |
| Forbidden links | `si_forbidden()` `si_rule_phenology()` `si_rule_morphology()` `si_rule_elevation()` `si_rule_taxonomy()` `si_rule_custom()` |
| Metrics | `si_metrics()` `si_degenerate()` `si_null()` `si_species_metrics()` `si_metric_map()` |
| Change | `si_beta()` `si_rewiring()` `si_novelty()` `si_network_extinction()` `si_compare()` |
| Robustness | `si_robustness()` `si_coextinction()` `si_robustness_map()` |
| Function | `si_function()` `si_service_change()` `si_redundancy()` `si_uncertainty()` |
| Simulation | `si_simulate_species()` `si_simulate_climate()` `si_simulate_metaweb()` `si_experiment()` `si_sensitivity()` |

## Bundled data

- `data(chile_pp)`: a 20 × 20 subset of the Chilean plant–pollinator system at
  10 arc-minute resolution, with present-day and 2070 suitability surfaces
  (continuous and binary), the documented metaweb, and phenology derived from
  dated interaction records.
- `si_example_rasters()`: a second 20 × 20 subset stored as GeoTIFF files at
  20 arc-minute resolution, with a study-area polygon and an interaction table.

## Getting help

- Help for a function: `?si_overlap`
- All functions: `help(package = "spinet")`
- Questions and bug reports: <https://github.com/ehsanrahimi666/spinet/issues>

## Citation

Rahimi, E. & Jung, C. (2026) spinet: an R package for projecting ecological
interaction networks through space and time from species distribution models.
R package version 0.1.2. <https://github.com/ehsanrahimi666/spinet>

Methodological background:

- Rahimi, E. & Jung, C. (2024) A new SDM-based approach for assessing climate
  change effects on plant–pollinator networks. *Insects* 15, 842.
  <https://doi.org/10.3390/insects15110842>
- Rahimi, E. & Jung, C. (2026) How can we incorporate species interactions into
  SDM-based climate change modelling? *Community Ecology*.
  <https://doi.org/10.1007/s42974-026-00330-4>

## Licence

GPL-3
