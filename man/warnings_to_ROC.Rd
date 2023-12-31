% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EWS.R
\name{warnings_to_ROC}
\alias{warnings_to_ROC}
\title{Convert warnings to Receiver Operating Curve (ROC)}
\usage{
warnings_to_ROC(EWS_warnings, grouping_vars)
}
\arguments{
\item{EWS_warnings}{Dataframe with warnings}

\item{grouping_vars}{Names of grouping variables}
}
\value{
Dataframe with true positive rate, true negative rate, false positive rate, and false negative rate per critical value
}
\description{
Convert warnings to Receiver Operating Curve (ROC)
}
