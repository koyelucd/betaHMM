% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/threshold_identification.R
\name{threshold_identification}
\alias{threshold_identification}
\title{Function to identify the thresholds in a DNA sample
       having multiple replicates}
\usage{
threshold_identification(
  data,
  M,
  N,
  parameter_estimation_only = FALSE,
  seed = NULL
)
}
\arguments{
\item{data}{A dataframe of dimension \eqn{C \times N} containing
methylation values for \eqn{C} CpG sites from \eqn{R=1}
treatment group and having \eqn{N} replicates or DNA samples
collected from \eqn{N} patients.}

\item{M}{Number of methylation states to be identified in a DNA sample.}

\item{N}{Number of patients or DNA sample replicates collected for each
treatment group.}

\item{parameter_estimation_only}{If only model parameters are to be estimated
then value is TRUE else FALSE.}

\item{seed}{Seed to allow for reproducibility (default = NULL).}
}
\value{
A list containing:
\itemize{
   \item A - The transition matrix for the BHMM model.
   \item tau - The initial distribution for the BHMM model.
   \item phi - The shape parameters for the observation sequence data
               in the BHMM model.
   }
}
\description{
Function to identify the thresholds in a DNA sample
       having multiple replicates
}
\details{
Function to identify the thresholds in a DNA sample
       having multiple replicates
}
