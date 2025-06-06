% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps_prioritize.R
\name{plot_lambda}
\alias{plot_lambda}
\title{Plot alternative lambda values}
\usage{
plot_lambda(lambda = c(-1, -0.5, 0, 0.5, 2, 1))
}
\arguments{
\item{lambda}{A vector of lambda values to plot}
}
\value{
Plots a figure
}
\description{
Show a plot illustrating alternative values for the \code{lambda} parameter in \link{ps_prioritize}. Lambda determines the shape of
the "benefit" function that determines the conservation value of protecting a given proportion of the geographic range of a
species or clade. Positive values place a higher priority on protecting additional populations of largely unprotected taxa,
whereas negative values place a higher priority on protecting additional populations of relatively well-protected taxa. The
default value used by \link{ps_prioritize} is 1.
}
\examples{
plot_lambda()
plot_lambda(seq(0, 3, .1))

}
