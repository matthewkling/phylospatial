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





