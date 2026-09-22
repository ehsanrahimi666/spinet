test_that("interaction vocabulary resolves labels, RO ids and sign shorthands", {
  expect_equal(si_interaction("pollinatedBy")$sign, "++")
  expect_equal(si_interaction("RO_0002456")$sign, "++")
  expect_equal(si_interaction("+-")$cascade, "consumer")
  expect_equal(si_interaction("-+")$cascade, "consumer")   # order-insensitive
  expect_equal(si_interaction("--")$mode, "unipartite")
  expect_equal(si_interaction("--")$overlap_sign, -1L)
  expect_equal(si_interaction("++")$overlap_sign, 1L)
  expect_error(si_interaction("not_a_type"), "Unknown interaction type")
  expect_error(si_interaction("++", overlap_sign = 5), "must be -1 or 1")
})

test_that("cascade rules differ by sign structure", {
  expect_equal(unname(si_cascade_rule(si_interaction("++"))),
               c("extinct", "extinct"))
  expect_equal(unname(si_cascade_rule(si_interaction("+-"))),
               c("released", "extinct"))
  expect_equal(unname(si_cascade_rule(si_interaction("+0"))),
               c("neutral", "extinct"))
  expect_equal(unname(si_cascade_rule(si_interaction("--"))),
               c("released", "released"))
})

test_that("si_stack validates and normalises", {
  m <- matrix(runif(50), 10, 5)
  st <- si_stack(m, level = "plant")
  expect_s3_class(st, "si_stack")
  expect_equal(nrow(st$values), 10L)
  expect_false(st$binary)
  expect_true(si_stack(matrix(c(0, 1), 10, 5))$binary)
  expect_equal(si_match_names("Villa", "Villa_")$map$b, "Villa_")
  expect_length(si_match_names("Villa", "Villa_")$only_a, 0L)
  expect_error(si_stack("/no/such/dir/at/all"), "not found|No rasters")
})

test_that("metaweb construction from edge lists and matrices agree", {
  el <- data.frame(p = c("p1", "p1", "p2"), a = c("a1", "a2", "a2"))
  mw <- si_metaweb(el, names_A = paste0("p", 1:3), names_B = paste0("a", 1:2))
  expect_equal(sum(mw$matrix), 3)
  expect_equal(dim(mw$matrix), c(3L, 2L))
  expect_true(si_metaweb(NULL, paste0("p", 1:3), paste0("a", 1:2))$unconstrained)
})

# ---------------------------------------------------------------------------
# The structural results the package exists to enforce
# ---------------------------------------------------------------------------

test_that("unconstrained continuous fields are rank one: H2' is exactly zero", {
  set.seed(1)
  p <- runif(40); a <- runif(30)
  h <- si_h2prime(outer(p, a), components = TRUE)
  expect_lt(abs(h[["H2max"]] - h[["H2"]]), 1e-9)
  expect_equal(unname(h[["H2prime"]]), 0)
})

test_that("unconstrained binary fields are complete bipartite graphs", {
  set.seed(2)
  P <- matrix(rbinom(200 * 10, 1, 0.5), 200, 10)
  A <- matrix(rbinom(200 * 8, 1, 0.5), 200, 8)
  f <- si_overlap(P, A, NULL, method = "binary", check = FALSE)
  mt <- suppressWarnings(si_metrics(f, what = c("connectance", "NODF")))
  expect_equal(stats::sd(mt$connectance, na.rm = TRUE), 0)
  expect_equal(stats::sd(mt$NODF, na.rm = TRUE), 0)
  expect_equal(mean(mt$connectance, na.rm = TRUE), 1)
  expect_equal(mean(mt$NODF, na.rm = TRUE), 0)
})

test_that("si_degenerate detects all three failure modes", {
  set.seed(3)
  P <- matrix(runif(100 * 6), 100, 6, dimnames = list(NULL, paste0("A", 1:6)))
  A <- matrix(runif(100 * 5), 100, 5, dimnames = list(NULL, paste0("B", 1:5)))
  expect_true("rank-one" %in%
                si_degenerate(si_overlap(P, A, NULL, "product", check = FALSE),
                              warn = FALSE)$issue)
  expect_true("complete-bipartite" %in%
                si_degenerate(si_overlap(P, A, NULL, "binary", check = FALSE),
                              warn = FALSE)$issue)
  mw <- si_simulate_metaweb(6, 5, 0.4, seed = 1)
  f <- si_overlap(P, A, mw, "product", floor = 0, check = FALSE)
  expect_true("floorless-continuous" %in% si_degenerate(f, warn = FALSE)$issue)
  expect_warning(si_overlap(P, A, NULL, "product"), "degenerate")
  expect_equal(nrow(si_degenerate(si_example_field(), warn = FALSE)), 0L)
})

test_that("a static constraint layer permits no rewiring gains", {
  set.seed(4)
  A1 <- si_simulate_species(10, 15, 15, prefix = "A", seed = 4)
  B1 <- si_simulate_species(8, 15, 15, prefix = "B", seed = 5)
  A2 <- si_simulate_climate(A1, "decline", severity = 0.3)
  B2 <- si_simulate_climate(B1, "decline", severity = 0.3)
  mw <- si_simulate_metaweb(10, 8, 0.3, seed = 4,
                            names_A = A1$names, names_B = B1$names)
  th <- function(s) si_threshold(s, "fixed", 0.3)
  f1 <- si_overlap(th(A1), th(B1), mw, method = "binary", check = FALSE)
  f2 <- si_overlap(th(A2), th(B2), mw, method = "binary", check = FALSE)
  b <- suppressWarnings(si_beta(f1, f2))
  expect_equal(mean(b$beta_OS, na.rm = TRUE), 0)
  expect_equal(sum(si_rewiring(f1, f2)$realised), 0)
  expect_warning(si_beta(f1, f2), "identically zero")
})

test_that("rewiring requires variance in the phenological shift", {
  set.seed(6)
  A1 <- si_simulate_species(12, 15, 15, prefix = "A", seed = 6)
  B1 <- si_simulate_species(10, 15, 15, prefix = "B", seed = 7)
  A2 <- si_simulate_climate(A1, "decline", severity = 0.25)
  B2 <- si_simulate_climate(B1, "decline", severity = 0.25)
  th <- function(s) si_threshold(s, "fixed", 0.3)
  gains <- function(sd_) {
    ph <- si_phenology_simulate(12, 10, shift_mean = -10, shift_sd = sd_,
                                seed = 8)
    m1 <- si_metaweb(si_rule_phenology(ph, "current", 0.3)$matrix,
                     names_A = A1$names, names_B = B1$names)
    m2 <- si_metaweb(si_rule_phenology(ph, "future", 0.3)$matrix,
                     names_A = A1$names, names_B = B1$names)
    sum(si_rewiring(si_overlap(th(A1), th(B1), m1, method = "binary", check = FALSE),
                    si_overlap(th(A2), th(B2), m2, method = "binary", check = FALSE))$realised)
  }
  expect_equal(gains(0), 0)                 # uniform shift: no rewiring
  expect_gt(gains(60), 0)                   # variable shift: rewiring
})

# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

test_that("NODF behaves on reference matrices", {
  perfect <- outer(1:6, 1:6, function(i, j) as.numeric(j <= 7 - i))
  expect_gt(si_nodf(perfect), 90)
  expect_equal(si_nodf(matrix(1, 6, 6)), 0)   # complete graph
  expect_true(is.na(si_nodf(matrix(1, 1, 5))))
  expect_lt(si_nodf(diag(6)), 1)             # perfectly modular
})

test_that("H2prime is bounded and zero for outer products", {
  set.seed(7)
  w <- matrix(rpois(200, 3), 20, 10)
  h <- si_h2prime(w)
  expect_gte(h, 0); expect_lte(h, 1)
  expect_equal(si_h2prime(outer(runif(10), runif(8))), 0)
})

test_that("evenness is one for a flat matrix", {
  expect_equal(si_evenness(matrix(1, 4, 4)), 1)
  expect_lt(si_evenness(matrix(c(100, 1, 1, 1), 2, 2)), 1)
})

test_that("simulated architectures differ as designed", {
  nest <- si_simulate_metaweb(30, 25, 0.15, architecture = "nested", seed = 1)
  modu <- si_simulate_metaweb(30, 25, 0.15, architecture = "modular", seed = 1)
  expect_gt(si_nodf(nest$matrix), si_nodf(modu$matrix))
  if (requireNamespace("igraph", quietly = TRUE))
    expect_gt(si_modularity(modu$matrix), si_modularity(nest$matrix))
  expect_lt(abs(si_connectance(nest) - 0.15), 0.06)
})

# ---------------------------------------------------------------------------
# Robustness
# ---------------------------------------------------------------------------

test_that("rewiring raises robustness, and constrained rewiring raises it less", {
  set.seed(9)
  mw <- si_simulate_metaweb(20, 16, 0.2, architecture = "nested", seed = 9)
  w <- mw$matrix * matrix(runif(320, 1, 10), 20, 16)
  r0 <- si_robustness(w, remove = "A", n_sim = 30, seed = 1)
  r1 <- si_robustness(w, remove = "A", rewiring = matrix(1, 20, 16),
                      n_sim = 30, seed = 1)
  r2 <- si_robustness(w, remove = "A",
                      rewiring = si_rule_morphology(runif(20), runif(16),
                                                    rule = "gower")$matrix,
                      n_sim = 30, seed = 1)
  expect_gt(r1$robustness, r0$robustness)
  expect_gte(r1$robustness, r2$robustness)
  expect_true(all(c(r0$robustness, r1$robustness, r2$robustness) <= 1))
})

test_that("cascade rule controls which level can coextinct", {
  set.seed(10)
  w <- si_simulate_metaweb(12, 10, 0.3, seed = 10)$matrix
  # removing the resource level of an antagonism releases the other level
  expect_equal(si_robustness(w, remove = "B",
                             interaction = si_interaction("+-"),
                             n_sim = 5)$robustness, 1)
  expect_lt(si_robustness(w, remove = "A",
                          interaction = si_interaction("++"),
                          n_sim = 20)$robustness, 1)
})

test_that("climate-ordered extinction requires a change vector", {
  w <- si_simulate_metaweb(10, 8, 0.3, seed = 11)$matrix
  expect_error(si_robustness(w, sequence = "climate"), "suitability_change")
  expect_silent(si_robustness(w, sequence = "climate",
                              suitability_change = runif(10), n_sim = 2))
})

# ---------------------------------------------------------------------------
# Forbidden links and experiments
# ---------------------------------------------------------------------------

test_that("forbidden rules combine and can be attributed", {
  set.seed(12)
  r1 <- matrix(rbinom(300, 1, 0.7), 20, 15)
  r2 <- matrix(rbinom(300, 1, 0.6), 20, 15)
  fl <- si_forbidden(a = r1, b = r2)
  expect_equal(sum(fl$matrix > 0), sum(r1 > 0 & r2 > 0))
  p <- si_forbidden_partition(fl)
  expect_equal(nrow(p), 2L)
  expect_error(si_forbidden(r1, matrix(1, 3, 3)), "same dimensions")
  expect_error(si_forbidden(), "at least one rule")
  expect_error(si_rule_custom(matrix(2, 3, 3)), "\\[0, 1\\]")
})

test_that("morphological barrier forbids links the consumer cannot reach", {
  corolla <- c(10, 30); proboscis <- c(5, 35)
  m <- si_rule_morphology(corolla, proboscis, rule = "barrier")$matrix
  expect_equal(m[2, 1], 0)   # short proboscis cannot reach deep corolla
  expect_equal(m[1, 2], 1)
})

test_that("elevation rule forbids non-overlapping ranges", {
  m <- si_rule_elevation(c(0, 2000), c(500, 2500), c(0), c(400))$matrix
  expect_equal(m[1, 1], 1)
  expect_equal(m[2, 1], 0)
})

test_that("experiments run and reproduce the connectance dependence", {
  ex <- suppressWarnings(si_experiment(
    list(richness = 12, connectance = c(1, 0.3), severity = 0.4, grid = 10),
    replicates = 1, cells = 20, seed = 1))
  expect_s3_class(ex, "si_experiment")
  expect_equal(nrow(ex), 2L)
  # unconstrained runs cannot change topology
  expect_equal(ex$d_connectance[ex$connectance == 1], 0)
  expect_equal(ex$d_NODF[ex$connectance == 1], 0)
})

test_that("sensitivity needs a varying factor", {
  ex <- suppressWarnings(si_experiment(list(richness = 10, grid = 8),
                                       replicates = 2, cells = 10, seed = 2))
  expect_error(si_sensitivity(ex, "beta_WN"), "No design factor varies")
})

test_that("edge cases return NA rather than failing", {
  expect_true(is.na(si_nodf(matrix(0, 3, 3))))
  expect_true(is.na(si_h2prime(matrix(0, 3, 3))))
  expect_true(is.na(si_evenness(matrix(0, 3, 3))))
  m <- si_web_metrics(matrix(0, 0, 0))
  expect_equal(unname(m[["links"]]), 0)
  f <- si_example_field(n_cells = 5)
  expect_error(si_local(f, 999), "out of range")
})

# ---------------------------------------------------------------------------
# Regression tests for the two bugs found by running the code on real
# plant-hummingbird data. The original suite covered only simulated phenology
# and never passed `observed` to si_forbidden(), so both slipped through.
# ---------------------------------------------------------------------------

test_that("si_rule_phenology keeps level-A by level-B orientation for observed data", {
  nA <- 7L; nB <- 4L                      # deliberately non-square
  A <- matrix(0L, 12, nA, dimnames = list(month.abb, paste0("plant", 1:nA)))
  B <- matrix(0L, 12, nB, dimnames = list(month.abb, paste0("bird", 1:nB)))
  A[1:6, ] <- 1L; B[4:9, ] <- 1L          # every pair shares Apr-Jun
  r <- si_rule_phenology(list(A = A, B = B), threshold = 0.1)

  expect_equal(dim(r$matrix), c(nA, nB))
  expect_equal(rownames(r$matrix), colnames(A))
  expect_equal(colnames(r$matrix), colnames(B))
  expect_true(all(r$matrix == 1))         # was 0 everywhere before the fix

  # and the silent downstream failure it caused
  n <- 6L
  A2 <- matrix(0L, 12, n, dimnames = list(month.abb, paste0("plant", 1:n)))
  B2 <- matrix(0L, 12, n, dimnames = list(month.abb, paste0("bird", 1:n)))
  A2[1:6, ] <- 1L; B2[4:9, ] <- 1L
  mw <- si_metaweb(si_rule_phenology(list(A = A2, B = B2), threshold = 0.1)$matrix,
                   names_A = colnames(A2), names_B = colnames(B2))
  expect_equal(sum(mw$matrix > 0), n * n)  # was 0 before the fix
})

test_that("si_rule_phenology grades partial overlap and rejects mismatched periods", {
  A <- matrix(0L, 12, 2, dimnames = list(month.abb, c("early", "late")))
  B <- matrix(0L, 12, 2, dimnames = list(month.abb, c("early", "late")))
  A[1:4, 1] <- 1L; A[7:12, 2] <- 1L
  B[1:4, 1] <- 1L; B[7:12, 2] <- 1L
  g <- si_rule_phenology(list(A = A, B = B), graded = TRUE)$matrix
  expect_equal(unname(diag(g)), c(1, 1))
  expect_equal(unname(g[1, 2]), 0)
  expect_error(si_rule_phenology(list(A = matrix(1, 12, 3), B = matrix(1, 6, 2))),
               "same number of rows")
})

test_that("si_forbidden(observed=) overrides by default and filters on request", {
  trait <- matrix(c(1,0,0,1, 0,0,1,1, 1,1,0,0, 0,0,0,1, 1,0,1,0), 5, 4,
                  byrow = TRUE,
                  dimnames = list(paste0("p", 1:5), paste0("a", 1:4)))
  obs <- matrix(0, 5, 4, dimnames = dimnames(trait))
  obs["p1", "a2"] <- 1        # documented, but the trait rule forbids it
  obs["p3", "a1"] <- 1        # documented and trait-permitted

  ov <- si_forbidden(traits = trait, observed = obs)
  expect_equal(ov$observed_action, "override")
  expect_equal(unname(ov$matrix["p1", "a2"]), 1)   # rescued by the record
  expect_equal(unname(ov$matrix["p1", "a1"]), 1)   # rule-permitted, undocumented
  expect_equal(ov$rescued_links, 1L)
  expect_equal(sum(ov$matrix > 0), sum(trait > 0) + 1)

  ft <- si_forbidden(traits = trait, observed = obs, observed_action = "filter")
  expect_equal(ft$observed_action, "filter")
  expect_equal(unname(ft$matrix["p1", "a2"]), 0)   # not trait-permitted
  expect_equal(unname(ft$matrix["p3", "a1"]), 1)
  expect_equal(sum(ft$matrix > 0), sum(trait > 0 & obs > 0))

  # override survives binarisation of a graded rule
  gr <- si_forbidden(traits = trait * 0.4, observed = obs,
                     binarise = TRUE, threshold = 0.5)
  expect_equal(unname(gr$matrix["p1", "a2"]), 1)
  expect_equal(sum(gr$matrix > 0), 2)
})

test_that("every si_rule_* returns a level-A by level-B matrix", {
  nA <- 7L; nB <- 4L
  d <- function(x) dim(x$matrix)
  expect_equal(d(si_rule_morphology(runif(nA), runif(nB))), c(nA, nB))
  expect_equal(d(si_rule_morphology(runif(nA), runif(nB), rule = "gower")), c(nA, nB))
  expect_equal(d(si_rule_elevation(runif(nA), runif(nA) + 1,
                                   runif(nB), runif(nB) + 1)), c(nA, nB))
  expect_equal(d(si_rule_abundance(runif(nA), runif(nB))), c(nA, nB))
  expect_equal(d(si_rule_taxonomy(rep(c("F1", "F2"), length.out = nA),
                                  list(a1 = "F1"), paste0("a", 1:nB))), c(nA, nB))
  expect_equal(d(si_rule_phenology(si_phenology_simulate(nA, nB, seed = 1))), c(nA, nB))
  expect_equal(dim(si_phenology_overlap(si_phenology_simulate(nA, nB, seed = 1))),
               c(nA, nB))
})

# ---------------------------------------------------------------------------
# Input handling: missing values, name matching and threshold edge cases
# ---------------------------------------------------------------------------

test_that("species names match across files, spreadsheets and read.csv", {
  nn <- spinet:::.norm_names
  expect_equal(nn("Toxomerus_calceolatus.."), "Toxomerus_calceolatus")  # check.names
  expect_equal(nn("Bombus terrestris"), "Bombus_terrestris")
  expect_equal(nn("Villa_"), "Villa")
  expect_equal(nn("Apis\u00a0mellifera"), "Apis_mellifera")               # no-break space
  m <- matrix(1, 2, 2, dimnames = list(c("p1", "p2"), c("a1", "zz")))
  expect_warning(si_metaweb(m, names_A = c("p1", "p2"), names_B = c("a1", "a2")),
                 "no row or column")
})

test_that("rasters with missing values are handled without silent cell loss", {
  r <- terra::rast(nrows = 10, ncols = 10, nlyrs = 2)
  terra::values(r) <- runif(200); names(r) <- c("s1", "s2")
  r[[1]][1:20] <- NA                       # first layer missing, second not
  expect_equal(nrow(si_stack(r, quiet = TRUE)$values), 100L)

  b <- terra::rast(nrows = 10, ncols = 10, nlyrs = 2)
  terra::values(b) <- runif(200); names(b) <- c("b1", "b2")
  b[[1]][1:30] <- NA; b[[2]][1:30] <- NA
  f <- suppressMessages(si_overlap(si_stack(r, quiet = TRUE),
                                   si_stack(b, quiet = TRUE), check = FALSE))
  expect_equal(f$n_cells, 70L)             # aligned on shared cells
  expect_error(si_stack(matrix(c(NA, .5, .2, .9), 2), na_value = NA), "na_value")
})

test_that("input validation catches common mistakes", {
  expect_warning(si_stack(matrix(c(0, 250, 1000, 40), 2)), "outside")
  expect_error(si_stack(matrix(runif(8), 4, dimnames = list(NULL, c("sp", "sp")))),
               "Duplicated")
  expect_error(si_simulate_climate(si_simulate_species(3, 10, 10, seed = 1),
                                   "shift", distance = 15), "distance")
  expect_warning(e <- si_rule_elevation(1000, 200, 300, 800)$matrix, "swapped")
  expect_equal(e[1, 1], 1)
})

test_that("quantile thresholds never make zero-suitability cells present", {
  v <- cbind(sparse = c(rep(0, 95), runif(5, .5, 1)), dense = runif(100))
  th <- suppressWarnings(si_threshold(si_stack(v), "quantile", 0.9))
  expect_lte(sum(th$values[, "sparse"]), 5)
  expect_warning(si_threshold(si_stack(v), "quantile", 0.9), "threshold of zero")
})

test_that("metrics, robustness and coextinction ignore missing and isolated species", {
  mm <- matrix(c(1, NA, 1, 1, 0, 1, 1, 1, 0), 3)
  expect_false(is.na(si_nodf(mm)))
  expect_false(is.na(si_evenness(mm)))
  W <- si_simulate_metaweb(10, 8, .3, seed = 2)$matrix
  W0 <- cbind(W, empty = 0)
  expect_equal(si_robustness(W, n_sim = 30, seed = 1)$robustness,
               si_robustness(W0, n_sim = 30, seed = 1)$robustness)
  cx <- si_coextinction(W0, rnorm(10), rnorm(9), n_sim = 5)
  expect_lt(cx$coextinct_B[2], 1)
})

test_that("si_forbidden aligns named rules by species, not by position", {
  r1 <- matrix(c(1, 0, 0, 1), 2, dimnames = list(c("p1", "p2"), c("a1", "a2")))
  fl <- si_forbidden(a = r1, b = r1[c(2, 1), ])
  expect_equal(unname(fl$matrix), unname(r1))
  bad <- r1; rownames(bad) <- c("p1", "p9")
  expect_error(si_forbidden(a = r1, b = bad), "same species")
})
