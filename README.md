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

## Quick start

```r
library(spinet)

# 1. read suitability surfaces (band 1 continuous, band 2 binary is a common
#    MaxEnt/flexsdm convention)
plants <- si_stack("sdm/plants_current",      band = 2, mask = study_area)
polls  <- si_stack("sdm/pollinators_current", band = 2, mask = study_area)

# 2. build a constraint layer
mw <- si_metaweb(observed_interactions, names_A = plants$names,
                                        names_B = polls$names)

# 3. the spatial interaction field
f <- si_overlap(plants, polls, mw, method = "binary")

# 4. per-cell network metrics, mapped
m <- si_metrics(f, what = c("links", "connectance", "NODF", "H2prime"))
terra::plot(si_metric_map(m))

# 5. compare with a future projection
b <- si_beta(f_current, f_future)     # β_WN = β_ST + β_OS
si_network_extinction(f_current, f_future)
```

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
