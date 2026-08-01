# ---------------------------------------------------------------------------
# spinet: Chile plant-pollinator case study
#
# Reproduces the empirical component of the analysis: 187 plants x 171
# pollinators, current and 2070 (SSP585), with and without the observed
# metaweb, on both the binary and continuous pathways.
#
# Inputs expected (edit `root`):
#   <root>/plants_current/*.tif        187 GeoTIFFs, band 1 = continuous,
#   <root>/plants_future/*.tif             band 2 = binary (max-sensitivity)
#   <root>/pollinators_current/*.tif   171 GeoTIFFs
#   <root>/pollinators_future/*.tif
#   <root>/Reference_matrix.csv        observed metaweb, 187 x 171
#   <root>/chile_polygon.shp           study-area mask
#   <root>/pollination_catalogue.csv   Muschett & Fonturbel catalogue (";" sep)
# ---------------------------------------------------------------------------

library(spinet)
library(terra)

root <- Sys.getenv("SPINET_CHILE", unset = "~/chile")
out  <- file.path(root, "spinet_output")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# --- 1. read -----------------------------------------------------------------
chile <- vect(file.path(root, "chile_polygon.shp"))

read_level <- function(folder, band)
  si_stack(file.path(root, folder), band = band, mask = chile,
           level = folder, quiet = TRUE)

# band 1 = continuous cloglog suitability; band 2 = binary (max sensitivity)
Pc_c <- read_level("plants_current", 1);      Pc_b <- read_level("plants_current", 2)
Pf_c <- read_level("plants_future", 1);       Pf_b <- read_level("plants_future", 2)
Ac_c <- read_level("pollinators_current", 1); Ac_b <- read_level("pollinators_current", 2)
Af_c <- read_level("pollinators_future", 1);  Af_b <- read_level("pollinators_future", 2)

# --- 2. the observed metaweb -------------------------------------------------
R <- as.matrix(read.csv(file.path(root, "Reference_matrix.csv"), row.names = 1))
mw <- si_metaweb(R, names_A = Pc_c$names, names_B = Ac_c$names)
message(sprintf("observed metaweb: %d links, connectance %.4f",
                sum(mw$matrix > 0), si_connectance(mw)))

# an unconstrained layer, for the contrast
mw0 <- si_metaweb(NULL, names_A = Pc_c$names, names_B = Ac_c$names)

# --- 3. four analytical pathways --------------------------------------------
pol <- si_interaction("pollinatedBy", "plant", "pollinator")
make <- function(P, A, M, meth, fl, sc)
  si_overlap(P, A, M, method = meth, floor = fl, interaction = pol,
             scenario = sc, check = FALSE)

fields <- list(
  # A: published binary framework, no constraint layer
  A_cur = make(Pc_b, Ac_b, mw0, "binary",  0,    "current"),
  A_fut = make(Pf_b, Af_b, mw0, "binary",  0,    "future"),
  # B: published continuous framework, no constraint layer
  B_cur = make(Pc_c, Ac_c, mw0, "product", 0,    "current"),
  B_fut = make(Pf_c, Af_c, mw0, "product", 0,    "future"),
  # C: binary x observed metaweb
  C_cur = make(Pc_b, Ac_b, mw,  "binary",  0,    "current"),
  C_fut = make(Pf_b, Af_b, mw,  "binary",  0,    "future"),
  # D: continuous x observed metaweb, pruned
  D_cur = make(Pc_c, Ac_c, mw,  "product", 0.05, "current"),
  D_fut = make(Pf_c, Af_c, mw,  "product", 0.05, "future"))

# --- 4. metrics, turnover and extinction ------------------------------------
metrics <- c("links", "connectance", "NODF", "H2prime", "H2_gap", "evenness")
res <- lapply(c("A", "B", "C", "D"), function(p) {
  f1 <- fields[[paste0(p, "_cur")]]; f2 <- fields[[paste0(p, "_fut")]]
  m1 <- si_metrics(f1, what = metrics); m2 <- si_metrics(f2, what = metrics)
  bb <- suppressWarnings(si_beta(f1, f2))
  rw <- si_rewiring(f1, f2)
  ex <- si_network_extinction(f1, f2)
  data.frame(pathway = p,
             t(colMeans(m1[, metrics], na.rm = TRUE)),
             t(colMeans(m2[, metrics], na.rm = TRUE)) -
               t(colMeans(m1[, metrics], na.rm = TRUE)),
             beta_WN = mean(bb$beta_WN, na.rm = TRUE),
             beta_ST = mean(bb$beta_ST, na.rm = TRUE),
             beta_OS = mean(bb$beta_OS, na.rm = TRUE),
             rewiring_gains = sum(rw$realised, na.rm = TRUE),
             networks_current = sum(ex$status %in% c("persists", "extinct")),
             networks_extinct = sum(ex$status == "extinct"))
})
tab <- do.call(rbind, res)
write.csv(tab, file.path(out, "chile_pathway_comparison.csv"), row.names = FALSE)
print(tab)

# --- 5. maps -----------------------------------------------------------------
mC <- si_metrics(fields$C_cur, what = metrics)
writeRaster(si_metric_map(mC), file.path(out, "chile_metrics_current.tif"),
            overwrite = TRUE)
ex <- si_network_extinction(fields$C_cur, fields$C_fut)
r <- rast(fields$C_cur$template); r[fields$C_cur$cells] <- as.integer(ex$status)
writeRaster(r, file.path(out, "chile_network_status.tif"), overwrite = TRUE)

# --- 6. derived phenology from the catalogue --------------------------------
cat_file <- file.path(root, "pollination_catalogue.csv")
if (file.exists(cat_file)) {
  ph <- si_phenology_from_records(read.csv(cat_file, sep = ";"),
                                  species_A = "scientificNamePlants",
                                  species_B = "scientificNameAnimals",
                                  date = "eventDate",
                                  names_A = Pc_c$names, names_B = Ac_c$names)
  saveRDS(ph, file.path(out, "chile_phenology.rds"))
  message(sprintf("derived phenology: %d plants, %d pollinators dated",
                  sum(rowSums(ph$A) > 0), sum(rowSums(ph$B) > 0)))
}

message("done -> ", out)
