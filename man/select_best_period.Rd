% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bifurcation.R
\name{select_best_period}
\alias{select_best_period}
\title{Select the best-fitting periodicity per bifurcation parameter}
\usage{
select_best_period(
  k_spread,
  factor_k,
  thresh_node,
  thresh_coord_spread,
  thresh_peak_idx_spread,
  nr_smooth,
  min_length_regime,
  X_names
)
}
\arguments{
\item{k_spread}{Dataframe specifying the fit of each period length k in terms of the distance in peak and trough coordinates and indices}

\item{factor_k}{Weighting of period length k; heavier weight means shorter k is preferred; factor_k = 0 means the optimal period length is chosen based solely on minimum spread}

\item{thresh_node}{Threshold under which timeseries is classified as node}

\item{thresh_coord_spread}{Threshold for distance in peak and trough coordinates that determines whether a cluster (i.e. period length) fits; if exceeded, denoted as chaotic or transitioning}

\item{thresh_peak_idx_spread}{Same as thresh_coord_spread but for peak and trough indices}

\item{nr_smooth}{Number of exceptions in a stable periodicity window to smooth over; nr_smooth = 0 means no smoothing}

\item{min_length_regime}{Minimum number of consecutive steps in the bifurcation parameter that have the same periodicity to qualify as a regime}

\item{X_names}{Names of variables in model}
}
\value{
Dataframe with periodicities
}
\description{
Select the best-fitting periodicity per bifurcation parameter
}
