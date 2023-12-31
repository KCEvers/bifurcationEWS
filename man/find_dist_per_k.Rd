% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bifurcation.R
\name{find_dist_per_k}
\alias{find_dist_per_k}
\title{Compute mean maximum distance per cluster}
\usage{
find_dist_per_k(ks, coord, peak_idx)
}
\arguments{
\item{ks}{Vector of possible period lengths}

\item{coord}{Peak and trough coordinates}

\item{peak_idx}{Peak and trough indices}
}
\value{
Distance in peak and trough coordinates and indices per cluster partitioning
}
\description{
Compute mean maximum distance per cluster
}
