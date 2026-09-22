## ---------------------------------------------------------------------------
## Supplementary Code S2. Complete Chilean analysis (187 plants x 171
## pollinators, 2,800 cells, present day and 2070 under SSP5-8.5).
##
## Set CHILE_DATA to a directory containing:
##   plants_current/, plants_future/, pollinators_current/, pollinators_future/
##       one two-band GeoTIFF per species (band 1 continuous, band 2 binary),
##       named after the species
##   Reference_matrix.csv   plant-by-pollinator matrix of documented links
##   chile_polygon.shp      national boundary used as the analysis mask
##   pollination_catalogue.csv   Muschett and Fonturbel (2022), ";"-separated
## Results are written to ANALYSIS_OUT (default "analysis_output").
## ---------------------------------------------------------------------------
suppressMessages({library(spinet); library(terra)})
root <- Sys.getenv("CHILE_DATA", "chile_data")
out  <- Sys.getenv("ANALYSIS_OUT", "analysis_output")
if (!dir.exists(root)) stop("Set CHILE_DATA to the input directory (see header).")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
chile <- vect(file.path(root, "chile_polygon.shp"))

## Step 1: suitability stacks (binary band) with one mask for all four
rd <- function(d, lvl) si_stack(file.path(root, d), band = 2, mask = chile, level = lvl, quiet = TRUE)
P_now <- rd("plants_current", "plant");      P_fut <- rd("plants_future", "plant")
A_now <- rd("pollinators_current", "pollinator"); A_fut <- rd("pollinators_future", "pollinator")

## Step 2: constraint layer from the documented interactions
web <- as.matrix(read.csv(file.path(root, "Reference_matrix.csv"), row.names = 1))
mw  <- si_metaweb(web, names_A = P_now$names, names_B = A_now$names)
cat("documented links:", sum(mw$matrix > 0), "| connectance:", round(si_connectance(mw), 4), "\n")

## Step 3: present-day and 2070 fields, same constraint layer
pol   <- si_interaction("pollinatedBy", "plant", "pollinator")
f_now <- si_overlap(P_now, A_now, mw, method = "binary", interaction = pol, scenario = "present", check = FALSE)
f_fut <- si_overlap(P_fut, A_fut, mw, method = "binary", interaction = pol, scenario = "2070", check = FALSE)

## Step 4: per-cell metrics and maps
what  <- c("links", "connectance", "NODF")
m_now <- si_metrics(f_now, what = what); m_fut <- si_metrics(f_fut, what = what)
writeRaster(si_metric_map(m_now, what), file.path(out, "metrics_present.tif"), overwrite = TRUE)
writeRaster(si_metric_map(m_fut, what), file.path(out, "metrics_2070.tif"), overwrite = TRUE)

## Step 5: change between periods
b  <- si_beta(f_now, f_fut)
rw <- si_rewiring(f_now, f_fut)
ex <- si_network_extinction(f_now, f_fut)
write.csv(b, file.path(out, "beta_partition.csv"), row.names = FALSE)
write.csv(ex, file.path(out, "network_fate.csv"), row.names = FALSE)

summ <- data.frame(
  metric  = c("links per cell", "connectance", "NODF"),
  present = c(mean(m_now$links), mean(m_now$connectance, na.rm = TRUE), mean(m_now$NODF, na.rm = TRUE)),
  y2070   = c(mean(m_fut$links), mean(m_fut$connectance, na.rm = TRUE), mean(m_fut$NODF, na.rm = TRUE)))
print(summ); print(table(ex$status))

## Unconstrained comparison: co-occurrence treated as interaction
f0 <- si_overlap(P_now, A_now, NULL, method = "binary", interaction = pol, check = FALSE)
L0 <- vapply(seq_len(f0$n_cells), function(k) sum(si_local(f0, k, drop = FALSE) > 0), numeric(1))
cat(sprintf("unconstrained links per cell %.1f (%.1f-fold the constrained value)\n",
            mean(L0), mean(L0) / mean(m_now$links)))

## Phenology from dated records and its internal check
cat_file <- file.path(root, "pollination_catalogue.csv")
if (file.exists(cat_file)) {
  txt <- iconv(readLines(cat_file, warn = FALSE), "UTF-8", "UTF-8", sub = "")
  rec <- read.csv(textConnection(txt), sep = ";")
  ph  <- si_phenology_from_records(rec, "scientificNamePlants", "scientificNameAnimals",
                                   "eventDate", names_A = P_now$names, names_B = A_now$names)
  print(si_phenology_check(ph, mw))
}
