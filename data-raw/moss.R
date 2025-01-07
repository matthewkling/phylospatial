## code to prepare `moss` data set

library(tidyverse)
library(phylospatial)
library(sf)
library(ape)
library(terra)

# load data
tree_file = "~/Documents/spatial_phylogenetics/ca_bryo_sphy/data/moss_chrono.tree"
comm_file = "~/Documents/spatial_phylogenetics/ca_bryo_sphy/results/comm/site_by_species.rds"
rast_file = "~/Documents/spatial_phylogenetics/ca_bryo_sphy/data/cpad_cced_raster_15km.tif"
comm <- readRDS(comm_file)
tree <- read.tree(file = tree_file)#[[2]]
template <- rast(rast_file)[[2]]

# clean and intersect species names
tree$tip.label <- str_remove_all(tree$tip.label, "-")
colnames(comm) <- str_remove_all(colnames(comm), "-")
xcom <- comm[, colnames(comm) %in% tree$tip.label]
tree <- drop.tip(tree, setdiff(tree$tip.label, colnames(comm)))
xcom <- xcom[, tree$tip.label]
# xcom[is.na(xcom)] <- 0 # NA values not allowed in sphy functions

## raster ##

# convert to raster
rr <- r <- to_spatial(xcom, template)

# aggregate
r <- terra::aggregate(r, 4, na.rm = TRUE)

# construct spatial phylo object
moss <- phylospatial(r, tree)

# save version with raster spatial data (since this can't be exported as rda)
saveRDS(moss, "inst/extdata/moss.rds")
writeRaster(moss$spatial, "inst/extdata/moss_raster.tif", overwrite = TRUE)
writeRaster(ps_get_comm(moss), "inst/extdata/moss_comm.tif", overwrite = TRUE)
write.tree(moss$tree, "inst/extdata/moss_tree.nex")


## polygons ##

pts <- st_as_sf(as.data.frame(rr, xy = TRUE), coords = c("x", "y"), crs = crs(rr))
p <- st_make_grid(st_as_sf(as.polygons(r)), n = dim(r)[2:1], square = FALSE) %>%
      st_as_sf() %>%
      mutate(hexid = 1:nrow(.))
pts <- pts %>%
      mutate(hexid = st_within(pts, p, sparse = FALSE) %>% apply(1, which)) %>%
      st_drop_geometry() %>%
      group_by(hexid) %>%
      summarize_all(mean)
p <- left_join(p, pts) %>%
      select(-hexid)

# drop empty polygons
# fixme: adding this section introduced "terminal scrambling" test error
p <- p[!is.na(rowSums(st_drop_geometry(p))), ]
p <- p[rowSums(st_drop_geometry(p)) > 0, ]

moss <- phylospatial(p, tree)

# save data
usethis::use_data(moss, overwrite = TRUE)




