## code to prepare `moss` and `moss_hex` data sets

library(tidyverse)
library(phylospatial)
library(sf)
library(ape)
library(terra)


# data ===========================

# load data
tree_file = "~/documents/spatial_phylogenetics/ca_bryo/ca_bryo_sphy/results/chronograms/moss_chrono.tree"
comm_file = "~/documents/spatial_phylogenetics/ca_bryo/ca_bryo_sphy/results/comm/site_by_species.rds"
rast_file = "~/documents/spatial_phylogenetics/ca_bryo/ca_bryo_sphy/data/cpad_cced_raster_15km.tif"
comm <- readRDS(comm_file)
tree <- read.tree(file = tree_file)#[[2]]
template <- rast(rast_file)[[2]]

# clean and intersect species names
tree$tip.label <- str_remove_all(tree$tip.label, "-")
colnames(comm) <- str_remove_all(colnames(comm), "-")
xcom <- comm[, colnames(comm) %in% tree$tip.label]
tree <- drop.tip(tree, setdiff(tree$tip.label, colnames(comm)))
xcom <- xcom[, tree$tip.label]


# raster ======================

# convert to raster
rr <- r <- to_spatial(xcom, template)

# aggregate
r <- terra::aggregate(r, 2, na.rm = TRUE)

# construct spatial phylo object
moss <- phylospatial(r, tree)

writeRaster(moss$spatial, "inst/extdata/moss_raster.tif", overwrite = TRUE)
# moss$spatial <- rast("inst/extdata/moss_raster.tif")
writeRaster(ps_get_comm(moss), "inst/extdata/moss_comm.tif", overwrite = TRUE)

# usethis::use_data(moss, overwrite = TRUE)
# saveRDS(moss, "inst/extdata/moss.rds")


# polygons ==================

p <- st_as_sf(as.polygons(r, round = FALSE, aggregate = FALSE, na.rm = FALSE))
moss_poly <- phylospatial(p, tree)

testthat::expect_equal(as.vector(moss_poly$comm), as.vector(moss$comm))
saveRDS(moss_poly$spatial, "inst/extdata/moss_polygons.rds")



# current protection level ===================

ps <- moss()
comm <- ps_get_comm(ps)[[1]]
reserves <- rast("~/documents/spatial_phylogenetics/ca_bryo/ca_bryo_sphy/data/protection_status.tif")
protected <- resample(reserves, comm, method = "mean")
protected <- mask(protected, comm)
protected[protected > .95] <- 1
names(protected) <- "protection"
crs(protected) <- crs(comm)
writeRaster(protected, "inst/extdata/moss_protection.tif", overwrite = T)


# population density (cost proxy) ==================

# GPW v4 population density (persons per km2), 2020, 2.5 arc-minute (~5 km) resolution
pop <- geodata::population(year = 2020, res = 0.5, path = tempdir())

# crop the global layer to the study area (plus a buffer) before projecting
pop <- crop(pop, ext(project(comm, crs(pop))) + 1, snap = "out")

# aggregate onto the moss grid as mean density
popdens <- project(pop, comm, method = "average")
popdens <- mask(popdens, comm)
names(popdens) <- "popdens"

writeRaster(popdens, "inst/extdata/moss_popdens.tif", overwrite = TRUE)
