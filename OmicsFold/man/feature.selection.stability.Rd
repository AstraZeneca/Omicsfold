% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_loadings.R
\name{feature.selection.stability}
\alias{feature.selection.stability}
\title{Extract sPLS-DA feature selection stability}
\usage{
feature.selection.stability(splsda.model, splsda.perf, comp)
}
\arguments{
\item{splsda.model}{Trained mixOmics sPLS-DA model.}

\item{splsda.perf}{Performance evaluation of the sPLS-DA model generated by
mixOmics perf.}

\item{comp}{Component number of which to retrieve feature selection
stabilities.}
}
\value{
A data frame containing the features selected for the specificed
component and their relative stability, as a proportion of trained models
0-1.
}
\description{
Extract feature selection stabilities for a given component from a
performance validated mixOmics sPLS-DA model Feature stabilities are an
important indicator of the confidence that a selected feature is predictive
for the outcome classes in the model, and hence (in combination with the
loading weight) is likely to be of biological significance.
}
\examples{
\dontrun{
feature.selection.stability(splsda.analysis, perf.splsda.analysis, comp=1)
}
}
