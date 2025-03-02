% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnalyzeSingle.R
\name{AnalyzeSingle}
\alias{AnalyzeSingle}
\title{Analysis for single samples or multiple samples without a DEG file provided}
\usage{
AnalyzeSingle(results_list, knowledgebase_version)
}
\arguments{
\item{results_list}{Results from CalculateInteractions()}
}
\value{
CSV Files
}
\description{
Analyzes calculated interaction scores by filtering by significant p-values
}
\details{
Converts the results from CalculateInteractions() into two CSV files:
\itemize{
\item 1: A long table of all interaction scores
\item 2: A filtered version with interaction scores that have a p-value less than 0.05
}
}
\keyword{internal}
