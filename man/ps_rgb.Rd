% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps_ordinate.R
\name{ps_rgb}
\alias{ps_rgb}
\title{Map phylospatial data onto RGB color bands}
\usage{
ps_rgb(ps, method = c("nmds", "cmds", "pca"), trans = identity, spatial = TRUE)
}
\arguments{
\item{ps}{A \code{phylospatial} object with a non-null \code{dissim} component, generated by \link{ps_add_dissim}.}

\item{method}{Ordination method, either "pca" (principal component analysis implemented via \code{stats::prcomp()}),
"cmds" (classical MDS, implemented via \code{stats::cmdscale()}), or "nmds" (the default, nonmetric MDS,
implemented via \code{vegan::metaMDS()}; this is slower but often preferred).}

\item{trans}{A function giving a transformation to apply to each dimension of the ordinated data.
The default is the identity function. Specifying \code{rank} generates a more uniform color distribution.}

\item{spatial}{Logical indicating whether a spatial object (inherited from \code{ps}) should be returned.
Default is TRUE.}
}
\value{
A matrix or spatial object with three variables containing RGB color values in the range 0-1.
}
\description{
Perform an ordination that reduces a spatial phylogenetic data set into three dimensions that can be
plotted as the RGB bands of color space to visualize spatial patterns of community phylogenetic composition.
This function is a wrapper around \code{ps_ordinate()}.
}
\examples{
ps <- ps_add_dissim(moss())
RGB <- ps_rgb(ps, method = "cmds")
terra::plotRGB(RGB * 255, smooth = FALSE)

}
