
# instantiate a new `phylospatial` object
new_phylospatial <- function(comm, tree, spatial, dissim = NULL, data_type, clade_fun){
      stopifnot(inherits(tree, "phylo"))
      stopifnot(is.character(data_type))
      structure(list(comm = comm,
                     tree = tree,
                     spatial = spatial,
                     data_type = data_type,
                     clade_fun = clade_fun,
                     dissim = NULL),
                class = "phylospatial")
}


#' Create a spatial phylogenetic object
#'
#' This function creates a \code{phylospatial} object. This is the core data type in the phylospatial library, and
#' is a required input to most other functions in the package. The two essential components of a spatial phylogenetic object
#' are a phylogenetic tree and an community data set.
#'
#' @param comm Community data representing the distribution of terminal taxa across sites. Can be a matrix with a column per terminal and
#' a row per site, a \link[terra]{SpatRaster} with one layer per terminal, or a `sf` data with a column per terminal. Taxa whose names do
#' not match between column/layer names in \code{comm} and tip labels in \code{tree} will be dropped with a warning (unless `build = FALSE`).
#' @param tree Phylogeny of class \link[ape]{phylo}. Terminals whose names do not match \code{comm} will be dropped with a warning (unless
#' `build = FALSE`). If this argument is not provided, terminals are assumed to follow a "star" tree with uniform branch lengths, which
#' will lead to non-phylogenetic versions of any analyses done with the resulting `phylospatial` object. Must be a rooted tree.
#' @param spatial An optional `SpatRaster` layer or `sf` object indicating site locations. The number of cells or rows must match \code{comm}.
#' Ignored if \code{comm} is a `SpatRaster` or `sf` object.
#' @param data_type Character giving the data type of \code{comm}. Must be "binary", "probability", "abundance", "auto" (the default), or "other".
#' This determines how community values for clades are calculated from the values for terminal taxa. If "binary" (presence-absence),
#' a clade is considered present in a site if any terminal in the clade is present. If "probability," clade probabilities are calculated
#' as the probability that at least one terminal is present in a site. If "abundance," clade abundances are calculated as the sum of
#' abundances for terminals in the clade in a site. If "auto," an attempt is made to guess which of these three data types was provided.
#' This argument is ignored if `clade_fun` is provided, or if `build = FALSE`. If "other", a custom `clade_fun` must be supplied.
#' @param clade_fun Function to calculate the local community weight for a clade based on community weights for tips found in a given location.
#' Must be either NULL (the default, in which case the default function for the selected \code{data_type} is used) or a summary function that
#' takes a numeric vector and returns a single numeric output. Ignored if \code{comm} already includes clade ranges.
#' @param build Logical indicating whether `comm` already includes clade ranges that should be used instead of building new ones.
#' Default is `TRUE`. If `FALSE`, `clade_fun` is ignored, no checks are performed to harmonize the tip labels and the community data, and
#' the columns of `comm` must exactly match the order of `tree` edges including tips and larger clades. If clade ranges are included in
#' `comm` but `build = TRUE`, they will be dropped and new clade ranges will be built.
#' @param check Logical indicating whether community data should be validated. Default is TRUE.
#' @param area_tol Numeric value giving tolerance for variation in the area of sites. Default is `0.01`. If the coefficient of variation in
#' the area or length of spatial units (e.g. grid cells) exceeds this value, an error will result. This check is performed because various
#' other functions in the library assume that sites are equal area. This argument is ignored if `check = FALSE` or if no spatial data is provided.
#'
#' @details
#' This function formats the input data as a `phylospatial` object. Beyond validating, cleaning, and restructing the data, the main operation
#' it performs is to compute community occurrence data for every internal clade on the tree. For a given clade and site, community data for
#' all the terminals in the clade are used to calculate the clade's occurrence value in the site. As described above, this calculation can
#' happen in various ways, depending on what type of community data you have (e.g. binary, probability, or abundance) and how you want to
#' summarize them. By default, the function tries to detect your `data_type` and use it to automatically select an appropriate summary
#' function as described above, but you can override this by providing your own function to `clade_fun`.
#'
#' You can also disable construction of the clade community matrix columns altogether by setting `build = FALSE`). This is atypical, but you
#' might want to use this option if you have your own distribution data data on all clades (e.g. from modeling occurrence probabilities for
#' clades in addition to terminal species), or if your community data comes from a previously-constructed `phylospatial` object.
#'
#' @return A `phylospatial` object, which is a list containing the following elements:
#' \describe{
#'  \item{"data_type":}{ Character indicating the community data type}
#'  \item{"tree":}{ Phylogeny of class `phylo`}
#'  \item{"comm":}{ Community matrix, including a column for every terminal taxon and every larger clade. Column order corresponds to tree edge order.}
#'  \item{"spatial":}{ A `SpatRaster` or `sf` providing spatial coordinates for the rows in `comm`. May be missing if no spatial data was supplied.}
#'  \item{"dissim":}{ A community dissimilary matrix of class `dist` indicating pairwise phylogenetic dissimilarity between sites. Missing unless
#'  \code{ps_dissim(..., add = TRUE)} is called.}
#' }
#'
#' @examples
#' \donttest{
#' # load species distribution data and phylogeny
#' comm <- terra::rast(system.file("extdata", "moss_comm.tif", package = "phylospatial"))
#' tree <- ape::read.tree(system.file("extdata", "moss_tree.nex", package = "phylospatial"))
#'
#' # construct `phylospatial` object
#' ps <- phylospatial(comm, tree)
#' ps
#' }
#' @export
phylospatial <- function(comm, tree = NULL, spatial = NULL,
                         data_type = c("auto", "probability", "binary", "abundance", "other"),
                         clade_fun = NULL, build = TRUE, check = TRUE, area_tol = 0.01){

      # checks
      data_type <- match.arg(data_type)
      stopifnot("Tree must be an object of class 'phylo'." = inherits(tree, c("phylo", "NULL")))
      stopifnot("Community data must be a `matrix`, `SpatRaster`, or `sf`." = inherits(comm, c("matrix", "SpatRaster", "sf")))
      stopifnot("Spatial reference must be a `SpatRaster` or `sf` object." = inherits(spatial, c("NULL", "SpatRaster", "sf")))
      if(inherits(spatial, "SpatRaster")) stopifnot("`spatial` must have the same number of grid cells as rows in `comm` data matrix." =
                                                          terra::ncell(spatial) == nrow(comm))
      if(inherits(comm, "matrix")) stopifnot("`comm` must have column names." = !is.null(colnames(comm)))


      # unpack community and spatial data
      if(inherits(comm, "SpatRaster")){
            spatial <- comm[[1]]
            spatial[] <- NA
            names(spatial) <- NULL
            comm <- terra::values(comm)
      }
      if(inherits(comm, "sf")){
            spatial <- sf::st_as_sf(sf::st_geometry(comm))
            comm <- as.matrix(sf::st_drop_geometry(comm))
      }

      if(is.null(tree)){
            warning("No phylogenetic tree was provided; any analyses using this `phylospatial` object will be non-phylogenetic.")
            tree <- star_tree(comm)
      }
      stopifnot("`tree` must be a rooted phylogeny, but is unrooted." = ape::is.rooted(tree))

      if(!build){
            check <- FALSE
            stopifnot("A `tree` must be provided if `clade_comm = TRUE`." =
                            inherits(tree, "phylo"))
            stopifnot("If `clade_comm = TRUE`, then the taxa in `comm` must match the edges in `tree`, both in number and in order." =
                            ncol(comm) == length(tree$edge.length))
      }

      # enforce equal area
      if(check & !is.null(spatial)){
            if(inherits(spatial, "SpatRaster") & terra::crs(spatial) == ""){
                  message("Raster data has no CRS, so can't check that cells are equal-area; it will be assumed that they are.")
            }else{
                  size <- area(spatial)
                  stopifnot("Spatial units are not equal area: the coefficient of varition in cell/polygon/line sizes exceeds 0.01." =
                                  (stats::sd(size$size) / mean(size$size)) < area_tol)
            }
      }

      # harmonize terminal taxa in tree and community
      if(check){
            tips <- intersect(tree$tip.label, colnames(comm))
            tips <- intersect(tips, colnames(comm)[colSums(comm, na.rm = TRUE) > 0])
            drop <- setdiff(tree$tip.label, tips)
            if(length(drop) > 0){
                  warning("Dropping ", length(drop), " tips from tree that are missing from community matrix or have no occurrences.")
                  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, tips))
            }
            drop <- setdiff(colnames(comm), tips)
            if(length(drop) > 0){
                  warning("Dropping ", length(drop), " taxa from community matrix that have no occurrences or are missing in tree.")
                  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, tips))
            }
            comm <- comm[, tree$tip.label]
      }

      # detect and validate data types
      is_binary <- function(x) identical(as.numeric(as.vector(x)), as.numeric(as.logical(x)))
      if(data_type == "auto"){
            if(any(na.omit(comm) < 0) | any(!is.finite(na.omit(comm)))) stop("Negative or infinite community values detected.")
            if(min(comm, na.rm = TRUE) >= 0 & max(comm, na.rm = TRUE) <= 1) data_type <- "probability"
            if(max(comm, na.rm = TRUE) > 1) data_type <- "abundance"
            if(is_binary(comm)) data_type <- "binary"
            if(build & check) message(paste("Community data type detected:", data_type))
            check <- FALSE
      }
      if(data_type == "binary"){
            if(check) if(! is_binary(comm)) stop("Binary data may only consist of 0 and 1, but other values were detected.")
      }
      if(data_type == "probability"){
            if(check) if(min(comm, na.rm = TRUE) < 0 | max(comm, na.rm = TRUE) > 1) stop("Probability data must be between 0 and 1, but other values were detected.")
      }
      if(data_type == "abundance"){
            if(check) if(min(comm, na.rm = TRUE) < 0) stop("Abundance data must be between nonnegative, but negative values were detected.")
      }

      # build clade ranges
      if(data_type == "other" & !inherits(clade_fun, "function")) stop("If `data_type = 'other'`, `clade_fun` must be a custom function.")
      if(is.null(clade_fun)) clade_fun <- get_clade_fun(data_type)
      if(build) comm <- build_tree_ranges(tree, comm, clade_fun)

      # scale branch lengths
      tree$edge.length <- tree$edge.length / sum(tree$edge.length)

      # create phylospatial object
      new_phylospatial(comm, tree, spatial, dissim = NULL, data_type, clade_fun)
}

