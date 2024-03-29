% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{dmrResults-class}
\alias{dmrResults-class}
\alias{dmrResults}
\title{dmrResults object and constructor}
\usage{
dmrResults(SummarizedExperiment)
}
\arguments{
\item{SummarizedExperiment}{a \code{dmrResults}
results.}
}
\value{
a  \code{\link[betaHMM:dmrResults]{dmrResults}} object
}
\description{
\code{dmrResults} is a subclass of \code{RangedSummarizedExperiment},
used to store the DMRs identified.
}
\details{
This constructor function would not typically be used by "end users".
This simple class extends the \code{RangedSummarizedExperiment} class of the
SummarizedExperiment package
to allow other packages to write methods for results
objects from the
\code{\link[betaHMM:dmr_identification]{dmr_identification}} function.
It is used by
to wrap up the results table.
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
