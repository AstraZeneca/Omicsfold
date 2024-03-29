% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_loadings.R
\name{get.diablo.top.loadings.with.stability}
\alias{get.diablo.top.loadings.with.stability}
\title{Get DIABLO top loadings and stabilities}
\usage{
get.diablo.top.loadings.with.stability(
  trained.model,
  perf.result,
  block,
  feature.count = 20
)
}
\arguments{
\item{trained.model}{A multiomics DIABLO model.}

\item{perf.result}{The results from a perf test on the above model.}

\item{block}{The name or index of the block to get the data for.}

\item{feature.count}{The number of top loadings from each component to
include.}
}
\value{
A dataframe containing the top (by absolute magnitude) factor
loadings for each component and selection stability values for those
features.
}
\description{
Get a dataframe showing the top ranked loading factors for a trained mixOmics
multi-omics DIABLO model along with selection stability values.  Where a
feature appears in more than one component, the selection stability value
shown is for the lowest numbered component.
}
\examples{
\dontrun{
get.diablo.top.loadings.with.stability(diablo.model, diablo.perf.result, 'species', 100)
}
}
