% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiFitR.R
\name{multiCDF_cont}
\alias{multiCDF_cont}
\title{Multiple Cumulative Distribution Functions for Continuous Variables}
\usage{
multiCDF_cont(x, seq_length, distributions)
}
\arguments{
\item{x}{The variable of which to get the CDF}

\item{seq_length}{The length of sequence to fit the distribution to}

\item{distributions}{The distributions to fit x against}
}
\value{
A dataframe with x, the real density, and the pdf of the desired
 distributions with length (nrows) equal to seq_length +1.
}
\description{
This function gets the cumulative distribution function for selected distributions
against a continuous, non-negative input variable. Possible distributions include "normal",
"lognormal", "gamma", "exponential", "cauchy", "t", "weibull", "logistic",
and "all".
}
