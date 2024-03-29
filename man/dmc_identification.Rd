% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/methods.R
\name{dmc_identification}
\alias{dmc_identification}
\alias{dmc_identification,betaHMMResults-method}
\alias{dmc_identification-methods}
\title{DMC identification from estimated betaHMM model parameters}
\usage{
dmc_identification(betaHMM_object, ...)

\S4method{dmc_identification}{betaHMMResults}(
  betaHMM_object,
  AUC_threshold = 0.8,
  uncertainty_threshold = 0.2,
  ...
)
}
\arguments{
\item{betaHMM_object}{An S4 object of class
\code{\link[betaHMM:betaHMMResults]{betaHMMResults}}}

\item{...}{extra arguments}

\item{AUC_threshold}{The threshold for AUC metric for each chromosome.}

\item{uncertainty_threshold}{The threshold for uncertainty of belonging to
a particular hidden state, for each chromosome.}
}
\value{
An S4 object of class \code{\link[betaHMM:dmcResults]{dmcResults}}.
}
\description{
This is the primary user interface for the
\code{\link[betaHMM:dmc_identification]{dmc_identification}} function.
Generic S4 methods are implemented to identify the DMCs from the estimated
betaHMM model parameters for each chromosome. The supported class is a
 \code{\link[betaHMM:betaHMMResults]{betaHMMResults}} object. The output
 is an S4 object of class of \code{\link[betaHMM:dmcResults]{dmcResults}}.
}
\examples{
## Use simulated data for the betaHMM workflow example
set.seed(12345)

## read files
data(sample_methylation_file)
data(sample_annotation_file)
# Run betaHMM function
beta_out <- betaHMM(sample_methylation_file[1:50,],
                    sample_annotation_file[1:50,],
                    M = 3, N = 4, R = 2,iterations=2,
                    parallel_process = FALSE, seed = 12345,
                    treatment_group = c("Benign","Tumour"))

## Run dmc_identification function
dmc_out <- dmc_identification(beta_out)

# Run dmr_identification function
dmr_out <- dmr_identification(dmc_out, parallel_process = FALSE)

# Plot functions
# Get the AUC values calculated for each hidden state
AUC_chr <- AUC(dmc_out)

## plot the uncertainty for each hidden state
plot(beta_out, chromosome = "1", what = "uncertainty")

}
\seealso{
\code{\link{betaHMM}}
}
