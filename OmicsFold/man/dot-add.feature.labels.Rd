% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blockrank.R
\name{.add.feature.labels}
\alias{.add.feature.labels}
\title{Internal function to add feature names}
\usage{
.add.feature.labels(blockrank.i, x.data)
}
\arguments{
\item{blockrank.i}{A list of numeric vectors containing BlockRank scores. One list item per data block}

\item{x.data}{List containing the centered and standardized original predictor matrix (with feature names)}
}
\value{
A list of numeric vectors containing BlockRank scores, with named features
}
\description{
Internal function to add feature names
}
