
test_that("si_forbidden objects are accepted wherever a constraint layer is", {
  ph <- si_phenology_simulate(8, 6, seed = 1)
  fl <- si_forbidden(phenology = si_rule_phenology(ph, threshold = 0.3),
                     binarise = TRUE)
  expect_equal(si_connectance(fl), mean(fl$metaweb$matrix > 0))
  expect_equal(si_connectance(fl), si_connectance(fl$metaweb))
})
