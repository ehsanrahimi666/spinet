# spinet 0.1.2

## Bug fixes

* `si_beta()`: in cells where no links are shared between the species present at
  both times (including networks that are lost or newly formed), `beta_OS` is now 0
  and `beta_ST` equals `beta_WN`, so the additive partition holds in every cell.
  These cells previously returned `NA`, which biased averages of `beta_ST` downwards.

* `si_connectance()` and the simulation functions now accept the object
  returned by `si_forbidden()` as a constraint layer, as `si_overlap()` already did.

* Species names are now matched after treating spaces, dots and underscores as
  equivalent separators and removing no-break spaces, so names read with
  `read.csv()` (which turns spaces and invalid characters into dots) match
  raster file names. `si_metaweb()` now warns when species cannot be matched
  and reports how many documented links could not be placed, instead of
  silently giving those species no links.

* `si_stack()` without a `mask` kept only cells that were non-missing in the
  first layer. It now keeps every cell in which any species has a prediction.
  `si_overlap()` aligns two stacks on the same grid by raster cell, keeping the
  cells shared by both, rather than failing when their missing-value patterns
  differ.

* `si_threshold()` with `rule = "quantile"` or `"prevalence"` scored every cell
  as present for species with zero suitability in most cells, because the
  threshold resolved to zero. Presence now also requires positive suitability,
  and a warning names the species affected.

* `si_robustness()` and `si_coextinction()` counted species with no
  interactions as secondary extinctions at the first step, biasing robustness
  downward. Isolated species are now removed before the simulation.

* `si_forbidden()` combined rules by position even when they listed species in
  different orders. Rules carrying species names are now aligned by name, and
  rules with different species sets are rejected.

* `si_nodf()`, `si_wnodf()`, `si_evenness()` and `si_modularity()` failed on
  matrices containing `NA`; missing entries are now treated as absent links.

## Input validation

* `si_stack()` rejects a non-finite `na_value`, rejects duplicated species
  names and warns when suitability values fall outside [0, 1].
* `si_rule_elevation()` swaps limits supplied in the wrong order, with a
  warning. `si_simulate_climate()` checks that shift distances fit the grid.

# spinet 0.1.1

## Bug fixes

* `si_rule_phenology()` returned the overlap matrix transposed when given
  observed activity matrices (`list(A = , B = )`), so plants and animals were
  swapped. Passed on to `si_metaweb()` this silently forbade every link when
  the two levels happened to have the same number of species, and produced a
  name-matching failure otherwise. The simulated-phenology path was unaffected,
  which is why the original tests missed it. The function now also checks that
  both activity matrices have the same number of time periods, and sets
  dimnames on the result so that alignment by species name works.

* `si_forbidden()` documented `observed` as overriding trait-based inference
  but implemented it as an additional multiplicative rule, so a documented
  interaction that a trait rule disallowed was removed, and any pair with no
  record was forbidden however plausible its traits. `observed` is now applied
  after the rules are combined and binarised, and a new argument
  `observed_action` chooses between `"override"` (the default, matching the
  documentation: a documented pair is permitted whatever the rules say) and
  `"filter"` (the previous behaviour, appropriate only when the interaction
  record is close to complete). The returned object reports how many documented
  links the rules would otherwise have forbidden, as `rescued_links`.

## Tests

* Regression tests for both bugs, plus an orientation test covering every
  `si_rule_*` constructor.

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
