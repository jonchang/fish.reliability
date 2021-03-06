% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/utils.R
\name{par_get_rtt_matrix}
\alias{par_get_rtt_matrix}
\title{Parallel version of BAMMtools getRateThroughTimeMatrix}
\usage{
par_get_rtt_matrix(ephy, start.time = NULL, end.time = NULL,
  nslices = 100, node = NULL, nodetype = "include",
  mc.cores = getOption("mc.cores", 2L))
}
\arguments{
\item{mc.cores}{number of cores to initialize a cluster with}

\item{...}{other parameters, see \link[BAMMtools]{getRateThroughTimeMatrix}}
}
\description{
Parallel version of BAMMtools getRateThroughTimeMatrix
}

