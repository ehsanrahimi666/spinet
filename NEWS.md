# spinet 0.1.0

First public release.

## Core

* Object model: `si_stack()`, `si_interaction()`, `si_metaweb()`, `si_traits()`,
  `si_grid()`. The interaction vocabulary follows Relations Ontology terms used
  by Global Biotic Interactions, so all six sign structures (`++`, `+-`, `+0`,
  `-0`, `--`, `00`) are supported with their own cascade rules.
* `si_overlap()` builds a lazy spatial interaction field with seven combination
  rules, storing the two suitability matrices and one constraint matrix rather
  than one interaction matrix per cell. At 187 x 171 species over 2,800 cells
  this is 8.7 MB instead of 783 MB.
* `si_stack()` reads either one raster per species (with `band` selecting the
  variant) or a single multi-layer raster whose layers are the species, and
  distinguishes the two automatically.

## Forbidden links

* `si_forbidden()` combines phenological, morphological, elevational,
  taxonomic, abundance-based and user-defined rules multiplicatively, with
  per-rule attribution through `si_forbidden_partition()`.
* Phenology at three levels of evidence: `si_phenology_simulate()`,
  `si_phenology_from_records()` (from dated Darwin Core / GloBI records) and
  user-supplied, with `si_phenology_check()` for internal validation.

## Metrics and diagnostics

* Native, dependency-free NODF, WNODF, H2', evenness, linkage density,
  generality, vulnerability and modularity.
* `si_degenerate()` detects rank-one, complete-bipartite and floorless-
  continuous fields and warns which metrics are constant by construction.
* `si_null()` gives standardised effect sizes against Patefield, curveball and
  shuffle null models, returning `NA` rather than a spurious z-score when the
  field is degenerate.

## Change, robustness and function

* `si_beta()` (Poisot partition), `si_rewiring()`, `si_novelty()`,
  `si_network_extinction()`, `si_compare()`.
* `si_robustness()` with and without rewiring, sign-aware through
  `si_cascade_rule()`; `si_coextinction()` supports climate-ordered extinction
  sequences.
* `si_function()`, `si_service_change()`, `si_redundancy()`,
  `si_uncertainty()`.

## Simulation

* `si_simulate_species()`, `si_simulate_climate()` (five non-spatial and three
  spatially explicit scenarios), `si_simulate_metaweb()` (connectance and
  architecture varied independently), `si_experiment()` and `si_sensitivity()`.

## Data

* `data(chile_pp)`: a 20 x 20 Chilean plant-pollinator subset with present and
  2070 suitability, the documented metaweb and derived phenology.
* `si_example_rasters()`: the same subset as GeoTIFFs plus the study-area
  polygon, so the raster reading path can be demonstrated without external
  files.
