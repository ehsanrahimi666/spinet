# spinet

**Spatial Interaction Networks Under Environmental Change**

`spinet` projects ecological interaction networks through space and time by
coupling per-species habitat suitability surfaces with an interaction
constraint layer. It sits *downstream* of species distribution modelling: it
consumes SDM output and never fits distribution models itself.

## Why the constraint layer is not optional

If a local network is built from co-occurrence alone, the local interaction
matrix is the outer product of the two suitability vectors, `W = p ⊗ a`. That
matrix has rank one, which has three exact consequences:

| | consequence |
|---|---|
| binary pathway | every local network is the complete bipartite graph: connectance ≡ 1, NODF ≡ 0 |
| continuous pathway | H₂′ ≡ 0, because the null expectation of H₂′ *is* the outer product of the marginals |
| any static constraint layer | interaction rewiring β_OS ≡ 0, because species shared by two realisations interact identically |

`spinet` therefore treats the constraint layer as the primary scientific object
rather than a refinement, and `si_degenerate()` warns which metrics are
mathematically constant *before* they are computed.

## Installation

```r
# from a local build
install.packages("spinet_0.1.0.tar.gz", repos = NULL, type = "source")

# from source
# remotes::install_github("ehsanrahimi666/spinet")
```

## Quick start — runnable, no external files

The package ships a real 20 x 20 Chilean plant–pollinator system, so this runs
as-is:

```r
library(spinet)
data(chile_pp)

# 1. suitability surfaces (grid + cells make the result mappable)
plants <- si_stack(chile_pp$plants_current_bin, level = "plant",
                   grid = chile_pp$grid, cells = chile_pp$cells)
polls  <- si_stack(chile_pp$pollinators_current_bin, level = "pollinator",
                   grid = chile_pp$grid, cells = chile_pp$cells)

# 2. the constraint layer
mw <- si_metaweb(chile_pp$metaweb)          # 114 links, connectance 0.285

# 3. the spatial interaction field
f_current <- si_overlap(plants, polls, mw, method = "binary",
                        scenario = "current")

# 4. per-cell network metrics, mapped
m <- si_metrics(f_current, what = c("links", "connectance", "NODF", "H2prime"))
summary(m)
terra::plot(si_metric_map(m, c("links", "connectance")))

# 5. compare with the future projection
plants_f <- si_stack(chile_pp$plants_future_bin, level = "plant",
                     grid = chile_pp$grid, cells = chile_pp$cells)
polls_f  <- si_stack(chile_pp$pollinators_future_bin, level = "pollinator",
                     grid = chile_pp$grid, cells = chile_pp$cells)
f_future <- si_overlap(plants_f, polls_f, mw, method = "binary",
                       scenario = "future")

si_beta(f_current, f_future)                # beta_WN = beta_ST + beta_OS
table(si_network_extinction(f_current, f_future)$status)
```

Expected output of the last two lines:

```
 beta_S beta_WN beta_ST beta_OS
 0.4775  0.5837  0.4321  0.0000     <- rewiring is zero: the layer is static

 persists   extinct colonised    absent
      485       110        23      2182
```

## The same thing from GeoTIFFs

`chile_pp` is a matrix, so it cannot show the raster reading path. The package
also ships the same subset as GeoTIFFs:

```r
library(terra)
si_example_rasters("all")          # 8 rasters + study-area polygon + metaweb

region <- vect(si_example_rasters("region"))
P  <- si_stack(si_example_rasters("plants_current_binary"),
               mask = region, level = "plant")
A  <- si_stack(si_example_rasters("pollinators_current_binary"),
               mask = region, level = "pollinator")
mw <- si_metaweb(as.matrix(read.csv(si_example_rasters("metaweb"),
                                    row.names = 1)))

f <- si_overlap(P, A, mw, method = "binary")
plot(si_metric_map(si_metrics(f, what = c("links", "connectance"))))
```

`si_stack()` accepts either convention and tells them apart automatically: a
directory or several file paths means one file per species, with `band`
selecting the variant inside each; a single multi-layer raster means the layers
are the species.

## Using your own species distribution models

Replace step 1 with paths to your own rasters. This is a **template** — the
paths, the study-area polygon and the interaction table are yours to supply:

```r
study_area <- terra::vect("gis/study_area.shp")
obs        <- read.csv("data/observed_interactions.csv")   # edge list or matrix

# band 1 = continuous suitability, band 2 = thresholded binary, a common
# MaxEnt / flexsdm convention. Use the same `mask` for every stack so that
# cell i means the same place in each.
plants <- si_stack("sdm/plants_current",      band = 2, mask = study_area)
polls  <- si_stack("sdm/pollinators_current", band = 2, mask = study_area)

mw <- si_metaweb(obs, names_A = plants$names, names_B = polls$names)
f  <- si_overlap(plants, polls, mw, method = "binary")
```

A complete worked script for a 187 x 171 system is at
`system.file("scripts", "chile_case_study.R", package = "spinet")`.

## All six interaction sign structures

The interaction vocabulary follows the Relations Ontology terms used by Global
Biotic Interactions, so any GloBI-standard table can be read directly.

```r
si_interaction("pollinatedBy")   # ++  mutualism      both coextinct
si_interaction("eatenBy")        # +-  antagonism     consumer dies, resource released
si_interaction("+0")             # +0  commensalism   one-way cascade
si_interaction("-0")             # -0  amensalism     target released
si_interaction("competesWith")   # --  competition    unipartite, overlap_sign = -1
si_interaction("00")             # 00  neutralism     null benchmark
```

The sign structure controls the coextinction rule, whether increasing spatial
overlap is good news or bad news, and which metrics are interpretable.

## Modules

| Module | Key functions |
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

`data(chile_pp)` — a 20 × 20 subset of a Chilean plant–pollinator system with
present and 2070 suitability surfaces, the documented metaweb, and phenology
derived from dated interaction records.

## Citation

Rahimi, E. & Jung, C. (2026) spinet: spatial interaction networks under
environmental change. R package version 0.1.0.

## Licence

GPL-3.
