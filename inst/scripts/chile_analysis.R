## ---------------------------------------------------------------------------
## Chile plant-pollinator case study
##
## Reproduces the empirical component of Rahimi & Jung: 187 plants x 171
## pollinators, present and 2070 (SSP5-8.5), with and without the observed
## metaweb, on both the binary and continuous pathways.
##
## Configure with two environment variables, or edit the defaults below:
##   SPINET_CHILE  directory holding the input data (default "chile_data")
##   SPINET_OUT    directory for results          (default "spinet_output")
##
## Expected contents of SPINET_CHILE:
##   plants_current/       187 GeoTIFFs, band 1 continuous, band 2 binary
##   plants_future/        187 GeoTIFFs
##   pollinators_current/  171 GeoTIFFs
##   pollinators_future/   171 GeoTIFFs
##   Reference_matrix.csv  observed metaweb, 187 x 171
##   chile_polygon.shp     study-area mask (+ .dbf .shx .prj)
##   pollination_catalogue.csv   Muschett & Fonturbel catalogue, ";" separated
##
## Note: some catalogue files contain invalid UTF-8. If read.csv() fails, run
##   iconv -f UTF-8 -t UTF-8 -c pollination_catalogue.csv > catalogue_clean.csv
## ---------------------------------------------------------------------------

suppressMessages({library(spinet); library(terra)})

root <- Sys.getenv("SPINET_CHILE", unset = "chile_data")
out  <- Sys.getenv("SPINET_OUT",   unset = "spinet_output")
if (!dir.exists(root))
  stop("Input directory not found: '", root, "'\n",
       "  Set it with Sys.setenv(SPINET_CHILE = \"/path/to/data\") ",
       "or edit `root` above.\n  Working directory is: ", getwd())
dir.create(out, showWarnings = FALSE, recursive = TRUE)

t0 <- Sys.time()
say <- function(...) cat(sprintf("[%5.0fs] ", as.numeric(Sys.time() - t0, units = "secs")),
                         ..., "\n", sep = "")

chile <- vect(file.path(root, "chile_polygon.shp"))
rd <- function(f, b) si_stack(file.path(root, f), band = b, mask = chile, level = f, quiet = TRUE)
Pc_c <- rd("plants_current",1); Pf_c <- rd("plants_future",1)
Ac_c <- rd("pollinators_current",1); Af_c <- rd("pollinators_future",1)
Pc_b <- rd("plants_current",2); Pf_b <- rd("plants_future",2)
Ac_b <- rd("pollinators_current",2); Af_b <- rd("pollinators_future",2)
say("read ", nrow(Pc_c$values), " cells x ", ncol(Pc_c$values), " plants / ", ncol(Ac_c$values), " pollinators")

R  <- as.matrix(read.csv(file.path(root, "Reference_matrix.csv"), row.names=1))
mw <- si_metaweb(R,  names_A = Pc_c$names, names_B = Ac_c$names)
mw0<- si_metaweb(NULL, names_A = Pc_c$names, names_B = Ac_c$names)
say("metaweb: ", sum(mw$matrix>0), " links, C = ", round(si_connectance(mw),4))

pol <- si_interaction("pollinatedBy","plant","pollinator")
mk <- function(P,A,M,meth,fl,sc) si_overlap(P,A,M,method=meth,floor=fl,interaction=pol,scenario=sc,check=FALSE)
F <- list(
  A_cur=mk(Pc_b,Ac_b,mw0,"binary",0,"current"),  A_fut=mk(Pf_b,Af_b,mw0,"binary",0,"future"),
  B_cur=mk(Pc_c,Ac_c,mw0,"product",0,"current"), B_fut=mk(Pf_c,Af_c,mw0,"product",0,"future"),
  C_cur=mk(Pc_b,Ac_b,mw ,"binary",0,"current"),  C_fut=mk(Pf_b,Af_b,mw ,"binary",0,"future"),
  D_cur=mk(Pc_c,Ac_c,mw ,"product",0.05,"current"), D_fut=mk(Pf_c,Af_c,mw,"product",0.05,"future"))

MET <- c("links","connectance","NODF","H2prime","H2_gap","evenness")
N <- F$A_cur$n_cells
set.seed(1); samp <- sort(sample.int(N, min(N, 400)))   # sample for the costly pathways

rows <- list()
for (p in c("A","B","C","D")) {
  f1 <- F[[paste0(p,"_cur")]]; f2 <- F[[paste0(p,"_fut")]]
  k <- if (p %in% c("A","B")) samp else seq_len(N)      # A/B are provably constant
  m1 <- suppressWarnings(si_metrics(f1, what=MET, cells=k))
  m2 <- suppressWarnings(si_metrics(f2, what=MET, cells=k))
  bb <- suppressWarnings(si_beta(f1,f2, cells=k))
  rw <- si_rewiring(f1,f2, cells=k)
  ex <- si_network_extinction(f1,f2)
  rows[[p]] <- data.frame(pathway=p, n_cells=length(k),
    cur_links=mean(m1$links,na.rm=TRUE), fut_links=mean(m2$links,na.rm=TRUE),
    cur_C=mean(m1$connectance,na.rm=TRUE), fut_C=mean(m2$connectance,na.rm=TRUE),
    sd_C=sd(m1$connectance,na.rm=TRUE),
    cur_NODF=mean(m1$NODF,na.rm=TRUE), fut_NODF=mean(m2$NODF,na.rm=TRUE), sd_NODF=sd(m1$NODF,na.rm=TRUE),
    cur_H2=mean(m1$H2prime,na.rm=TRUE), fut_H2=mean(m2$H2prime,na.rm=TRUE),
    cur_H2gap=mean(m1$H2_gap,na.rm=TRUE), sd_H2gap=sd(m1$H2_gap,na.rm=TRUE),
    beta_S=mean(bb$beta_S,na.rm=TRUE), beta_WN=mean(bb$beta_WN,na.rm=TRUE),
    beta_ST=mean(bb$beta_ST,na.rm=TRUE), beta_OS=mean(bb$beta_OS,na.rm=TRUE),
    gains=sum(rw$realised,na.rm=TRUE), losses=sum(rw$lost,na.rm=TRUE),
    nets_cur=sum(ex$status %in% c("persists","extinct")), nets_ext=sum(ex$status=="extinct"))
  say("pathway ", p, " done")
  if (p=="C") { saveRDS(list(m1=m1,m2=m2,bb=bb,ex=ex), file.path(out, "pathway_C.rds")) }
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path(out, "chile_pathways.csv"), row.names=FALSE)
print(t(tab))

## maps from pathway C
mC <- suppressWarnings(si_metrics(F$C_cur, what=MET))
saveRDS(list(metrics=mC, template=F$C_cur$template, cells=F$C_cur$cells), file.path(out, "chile_maps.rds"))
writeRaster(si_metric_map(mC), file.path(out, "chile_metrics_current.tif"), overwrite=TRUE)
ex <- si_network_extinction(F$C_cur, F$C_fut)
r <- rast(F$C_cur$template); r[F$C_cur$cells] <- as.integer(ex$status); names(r) <- "status"
writeRaster(r, file.path(out, "chile_network_status.tif"), overwrite=TRUE)
rl <- rast(F$C_cur$template); rl[F$C_cur$cells] <- ex$links_future - ex$links_current; names(rl) <- "link_change"
writeRaster(rl, file.path(out, "chile_link_change.tif"), overwrite=TRUE)
say("maps written")

## catalogue phenology
cta <- read.csv(file.path(root, "pollination_catalogue.csv"), sep=";")
ph <- si_phenology_from_records(cta, "scientificNamePlants","scientificNameAnimals","eventDate",
                                names_A = Pc_c$names, names_B = Ac_c$names)
chk <- si_phenology_check(ph, mw)
saveRDS(list(ph=ph, check=chk), file.path(out, "chile_phenology.rds"))
print(ph$coverage); print(unlist(chk))

## species-level range change (for the intro/results)
rng <- data.frame(
  level = c(rep("plant", ncol(Pc_b$values)), rep("pollinator", ncol(Ac_b$values))),
  species = c(Pc_b$names, Ac_b$names),
  cur = c(colSums(Pc_b$values), colSums(Ac_b$values)),
  fut = c(colSums(Pf_b$values), colSums(Af_b$values)))
rng$change <- (rng$fut - rng$cur) / pmax(rng$cur, 1)
write.csv(rng, file.path(out, "chile_range_change.csv"), row.names=FALSE)
say("species declining: ", round(100*mean(rng$change < 0),1), "% | mean change plants ",
    round(100*mean(rng$change[rng$level=="plant"]),1), "% pollinators ",
    round(100*mean(rng$change[rng$level=="pollinator"]),1), "%")

## bundled subset for the package: 20 x 20 best-sampled species
deg_p <- rowSums(mw$matrix); deg_a <- colSums(mw$matrix)
ip <- order(-deg_p)[1:20]; ia <- order(-deg_a)[1:20]
sub <- list(
  plants_current = Pc_c$values[, ip, drop=FALSE], plants_future = Pf_c$values[, ip, drop=FALSE],
  polls_current  = Ac_c$values[, ia, drop=FALSE], polls_future  = Af_c$values[, ia, drop=FALSE],
  plants_current_bin = Pc_b$values[, ip, drop=FALSE], plants_future_bin = Pf_b$values[, ip, drop=FALSE],
  polls_current_bin  = Ac_b$values[, ia, drop=FALSE], polls_future_bin  = Af_b$values[, ia, drop=FALSE],
  metaweb = mw$matrix[ip, ia, drop=FALSE],
  phenology_A = ph$A[, ip, drop=FALSE], phenology_B = ph$B[, ia, drop=FALSE],
  cells = F$C_cur$cells,
  template_wkt = as.character(crs(F$C_cur$template)),
  template_ext = as.vector(ext(F$C_cur$template)),
  template_dim = c(nrow(F$C_cur$template), ncol(F$C_cur$template)))
saveRDS(sub, file.path(out, "chile_subset.rds"), compress="xz")
say("subset saved: ", round(file.size(file.path(out, "chile_subset.rds"))/1e6,2), " MB")
say("ALL DONE")
