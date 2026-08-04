library(spinet)

## ---- 1. data ---------------------------------------------------------------
## The package ships a real 20 x 20 Chilean plant-pollinator system, so this
## runs with no external files.
data(chile_pp)

plants <- si_stack(chile_pp$plants_current_bin,      level = "plant")
polls  <- si_stack(chile_pp$pollinators_current_bin, level = "pollinator")

## attach the raster grid so results can be mapped
plants$template <- polls$template <- chile_pp$grid
plants$cells    <- polls$cells    <- chile_pp$cells

## ---- 2. constraint layer ---------------------------------------------------
mw <- si_metaweb(chile_pp$metaweb)
mw

## ---- 3. the spatial interaction field --------------------------------------
f_current <- si_overlap(plants, polls, mw, method = "binary",
                        scenario = "current")

## ---- 4. per-cell metrics, mapped -------------------------------------------
m <- si_metrics(f_current, what = c("links", "connectance", "NODF", "H2prime"))
summary(m)
terra::plot(si_metric_map(m, c("links", "connectance")))

## ---- 5. compare with the future projection ---------------------------------
plants_f <- si_stack(chile_pp$plants_future_bin,      level = "plant")
polls_f  <- si_stack(chile_pp$pollinators_future_bin, level = "pollinator")
plants_f$template <- polls_f$template <- chile_pp$grid
plants_f$cells    <- polls_f$cells    <- chile_pp$cells

f_future <- si_overlap(plants_f, polls_f, mw, method = "binary",
                       scenario = "future")

b <- si_beta(f_current, f_future)          # beta_WN = beta_ST + beta_OS
round(colMeans(b[, c("beta_S", "beta_WN", "beta_ST", "beta_OS")],
               na.rm = TRUE), 4)

table(si_network_extinction(f_current, f_future)$status)
