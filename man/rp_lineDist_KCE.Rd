% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RQA.R
\name{rp_lineDist_KCE}
\alias{rp_lineDist_KCE}
\title{Extract line lengths from recurrence matrix}
\usage{
rp_lineDist_KCE(
  RM,
  DLmin = 2,
  VLmin = 2,
  HLmin = 2,
  WLmin = 1,
  DLmax = length(Matrix::diag(RM)) - 1,
  VLmax = length(Matrix::diag(RM)) - 1,
  HLmax = length(Matrix::diag(RM)) - 1,
  WLmax = length(Matrix::diag(RM)) - 1,
  d = NULL,
  theiler = NA,
  AUTO = NULL,
  matrices = FALSE
)
}
\arguments{
\item{RM}{A binarized distance matrix}

\item{DLmin}{Minimal diagonal line length (default = \code{2})}

\item{VLmin}{Minimal vertical line length (default = \code{2})}

\item{HLmin}{Minimal horizontal line length (default = \code{2})}

\item{WLmin}{Minimal white line length (default = \code{1}) # KCE}

\item{DLmax}{Maximal diagonal line length (default = length of diagonal -1)}

\item{VLmax}{Maximal vertical line length (default = length of diagonal -1)}

\item{HLmax}{Maximal horizontal line length (default = length of diagonal -1)}

\item{WLmax}{Maximal white line length (default = length of diagonal -1) # KCE}

\item{d}{Vector of diagonals to be extracted from matrix \code{RP} before line length distributions are calculated. A one element vector will be interpreted as a windowsize, e.g., \code{d = 50} will extract the diagonal band \code{-50:50}. A two element vector will be interpreted as a band, e.g. \code{d = c(-50,100)} will extract diagonals \code{-50:100}. If \code{length(d) > 2}, the numbers will be interpreted to refer to individual diagonals, \code{d = c(-50,50,100)} will extract diagonals \verb{-50,50,100}. If \code{length(d)} is \code{NULL}, 1 or 2, the theiler window is applied before diagonals are extracted. The theiler window is ignored if \code{length(d)>2}, or if it is larger than the matrix or band indicated by parameter \code{d}. A warning will be given is a theiler window was already applied to the matrix.}

\item{theiler}{Theiler window}

\item{AUTO}{Auto-recurrence?}

\item{matrices}{Return matrices? (default = \code{FALSE})}
}
\value{
A list object with distributions of line lengths. If \code{matrices = TRUE} dataframes are returned whose columns represent the nonzero diagonals, verticals, or, horizontals.
}
\description{
Extract lengths of diagonal, vertical and horizontal line segments from a recurrence matrix.
}
\details{
Based on the Matlab function \code{linedists} by Stefan Schinkel, Copyright (C) 2009 Stefan Schinkel, University of Potsdam, http://www.agnld.uni-potsdam.de

References:
S. Schinkel, N. Marwan, O. Dimigen & J. Kurths (2009):
"Confidence Bounds of recurrence-based complexity measures
Physics Letters A,  373(26), pp. 2245-2250

Copyright (C) 2009 Stefan Schinkel, University of Potsdam
\url{http://www.agnld.uni-potsdam.de}
}
\author{
Fred Hasselman
}
