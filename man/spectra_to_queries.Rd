% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spectra_to_queries.R
\name{spectra_to_queries}
\alias{spectra_to_queries}
\title{Spectra to queries}
\usage{
spectra_to_queries(
  spectra = NULL,
  export = "data/interim/queries.tsv",
  beta = 0.5,
  dalton = 0.01,
  decimals = 4L,
  intensity_min = 0L,
  ions_max = 10L,
  n_skel_min = 5L,
  n_spec_min = 3L,
  ppm = 25L,
  senspe_min = 0L,
  sensitivity_min = 0L,
  specificity_min = 0L,
  zero_val = 0L
)
}
\arguments{
\item{spectra}{Spectra path}

\item{export}{Export path}

\item{beta}{Beta parameter of the F-score calculation}

\item{dalton}{Tolerance in Dalton. Default to 0.01}

\item{decimals}{Number of decimals for rounding. Default to 4}

\item{intensity_min}{Minimal intensity. Default to 0}

\item{ions_max}{Maximal number of ions in the query. Default to 10}

\item{n_skel_min}{Minimal number of individuals per skeleton. Default to 5}

\item{n_spec_min}{Minimal number of individuals where a signal has to be found. Default to 3}

\item{ppm}{Tolerance in parts per million Default to 25}

\item{senspe_min}{Minimal product of inner and outer ratios. Default to 0}

\item{sensitivity_min}{Minimal sensitivity. Default to 0}

\item{specificity_min}{Minimal specificity. Default to 0}

\item{zero_val}{Zero value for intensity. Default to 0}
}
\value{
A file with diagnostic query ions
}
\description{
This function converts spectra to queries.
}
\examples{
NULL
}
