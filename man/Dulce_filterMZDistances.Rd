% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnnotateMS2.r
\name{Dulce_filterMZDistances}
\alias{Dulce_filterMZDistances}
\title{MS/MS Distances Filtration and Highlighting}
\usage{
Dulce_filterMZDistances(data, sigdiff = NULL, n = NULL, mzdiff = NULL)
}
\arguments{
\item{data}{\code{data.frame} object with calculated distances.}

\item{sigdiff}{numeric vector with significant m/z differences values.}

\item{n}{numeric value. Indicates the amount of different significant
differences that a spectra should have in order to be included (not filtered out).
If \code{NULL}, then \code{n=length(sigdiff)}}

\item{mzdiff}{numeric value. Indicates the allowed range of mz values by which
differences should be filtered.}
}
\value{
This function returns a \code{data.frame} with the filtered spectra based
on the significant differences that were present.
}
\description{
This function filters the MS/MS spectra by the amount of
different meaningful differences.
}
\seealso{

}
