% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EWS.R
\name{warnings_to_ROC}
\alias{warnings_to_ROC}
\title{Convert warnings to Receiver Operating Curve (ROC)}
\usage{
warnings_to_ROC(
  EWS_warnings,
  sigma_crit_step,
  thresh_max_sigma_crit,
  grouping_vars = c("sigma_crit", "metric")
)
}
\arguments{
\item{EWS_warnings}{Dataframe with warnings}

\item{sigma_crit_step}{Step size in sequence of critical values of sigma}

\item{thresh_max_sigma_crit}{Maximum critical cut-off value to look for}

\item{grouping_vars}{Names of grouping variables}
}
\value{
Dataframe with true positive rate, true negative rate, false positive rate, and false negative rate per critical value
}
\description{
Convert warnings to Receiver Operating Curve (ROC)
}
