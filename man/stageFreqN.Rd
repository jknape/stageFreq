% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/functions.r
\name{stageFreqN}
\alias{stageFreqN}
\title{Compute stage probabilites by numeric integration.}
\usage{
stageFreqN(sampleTimes, dDur, mortRate, grid)
}
\arguments{
\item{sampleTimes}{Vector of times at which the cohort was sampled.}

\item{dDur}{A matrix of stage duration densities.}

\item{mortRate}{Mortality rates for each stage.}

\item{grid}{A grid over which numeric integration of stage duration probabilities is done.}
}
\value{
A matrix of stage counts. Each row gives the number of individuals in each stage for a given sampling time.
}
\description{
Compute stage probabilites by numeric integration.
}

