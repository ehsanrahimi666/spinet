
test_that("beta_WN = beta_ST + beta_OS wherever beta_WN is defined", {
  data(chile_pp, package = "spinet", envir = environment())
  st <- function(x, l) si_stack(x, level = l, grid = chile_pp$grid,
                                cells = chile_pp$cells, quiet = TRUE)
  mw <- si_metaweb(chile_pp$metaweb)
  f0 <- si_overlap(st(chile_pp$plants_current_bin, "plant"),
                   st(chile_pp$pollinators_current_bin, "pollinator"),
                   mw, method = "binary", check = FALSE)
  f1 <- si_overlap(st(chile_pp$plants_future_bin, "plant"),
                   st(chile_pp$pollinators_future_bin, "pollinator"),
                   mw, method = "binary", check = FALSE)
  b <- suppressWarnings(si_beta(f0, f1))
  d <- is.finite(b$beta_WN)
  expect_true(all(is.finite(b$beta_ST[d]) & is.finite(b$beta_OS[d])))
  expect_equal(b$beta_WN[d], b$beta_ST[d] + b$beta_OS[d])
  lost <- si_network_extinction(f0, f1)$status == "extinct"
  expect_true(all(b$beta_OS[lost] == 0))
  expect_equal(b$beta_ST[lost], b$beta_WN[lost])
})
