% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps_diversity.R
\name{clade_dist}
\alias{clade_dist}
\title{Pairwise distances among clades or nodes}
\usage{
clade_dist(tree, lineal = FALSE, edges = TRUE)
}
\arguments{
\item{tree}{A phylogeny of class \code{"phylo"}.}

\item{lineal}{Logical indicating whether to retain distances for pairs of nodes that are lineal ancestors/descendants.
If \code{FALSE} (the default), these are set to \code{NA}, retaining values only for node pairs that are collateral kin.}

\item{edges}{Logical indicating whether to return a distance matrix with a row for every edge in \code{tree}.
If \code{TRUE} (the default), rows/columns of the result correspond to \code{tree$edge}. If \code{FALSE}, rows/columns
correspond to nodes as in \code{ape::dist.nodes()}.}
}
\value{
A matrix of pairwise distances between nodes.
}
\description{
This function runs \code{ape::dist.nodes()} with some additional filtering and sorting. By default,
it returns distances between every pair of non-nested clades, i.e. every pair of collateral (non-lineal) nodes
including terminals and internal nodes.
}
\examples{
clade_dist(ape::rtree(10))
}
