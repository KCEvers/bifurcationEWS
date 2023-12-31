% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bifurcation.R
\name{get_bifurcation_range}
\alias{get_bifurcation_range}
\title{Generate bifurcation parameter sequence}
\usage{
get_bifurcation_range(
  bifpar_start,
  bifpar_end,
  pre_steps = 0,
  baseline_steps = 50,
  transition_steps = 50,
  post_steps = 0
)
}
\arguments{
\item{bifpar_start}{Bifurcation parameter value at start}

\item{bifpar_end}{Bifurcation parameter value at end; NA for null models}

\item{pre_steps}{Number of pre-steps, useful if transition did not start when planned}

\item{baseline_steps}{Number of baseline steps where the control parameter doesn't change}

\item{transition_steps}{Number of transition steps where the control parameter changes for every index in a step-wise manner}

\item{post_steps}{Number of post-transition steps}
}
\value{
List of bifurcation parameter values
}
\description{
Generate bifurcation parameter sequence
}
