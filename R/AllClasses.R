#' @rdname betaHMMResults
#' @export
#' @import SummarizedExperiment

setClass("betaHMMResults",
         contains = "RangedSummarizedExperiment",
         representation = representation(
           annotatedData="DataFrame"

         ))

## TODO:
## setValidity( ... )


#' betaHMMResults object and constructor
#'
#' \code{betaHMMResults} is a subclass of \code{RangedSummarizedExperiment},
#' used to store the betaHMM results as well as the annotated data useful
#'  for plotting.
#'
#' This constructor function would not typically be used by "end users".
#' This simple class extends the \code{RangedSummarizedExperiment} class of the
#' SummarizedExperiment package
#' to allow other packages to write methods for results
#' objects from the betaHMM package. It is used by
#' to wrap up the results table.
#'
#' @param SummarizedExperiment a \code{RangedSummarizedExperiment} of \code{betaHMM} results
#' @param annotatedData The annotated data passed as an input argument to the betaHMM package.
#'
#'
#' @return a betaHMMResults object
#' @docType class
#' @rdname betaHMMResults
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom S4Vectors DataFrame
#' @export
betaHMMResults <- function(SummarizedExperiment,
                         annotatedData) {
  se <- SummarizedExperiment
  if (!is(se, "RangedSummarizedExperiment")) {
    stop("'SummarizedExperiment' must be a RangedSummarizedExperiment object")
  }
  if(is.null(annotatedData)) annotatedData <- DataFrame(matrix(0, nrow=0, ncol=0))
  object <- new("betaHMMResults", se, annotatedData=annotatedData)
  return(object)
}

#' @rdname dmrResults
#' @export
#' @import SummarizedExperiment

setClass("dmrResults",
         contains = "RangedSummarizedExperiment")

## TODO:
## setValidity( ... )


#' dmrResults object and constructor
#'
#' \code{dmrResults} is a subclass of \code{RangedSummarizedExperiment},
#' used to store the dmrs identified.
#'
#' This constructor function would not typically be used by "end users".
#' This simple class extends the \code{RangedSummarizedExperiment} class of the
#' SummarizedExperiment package
#' to allow other packages to write methods for results
#' objects from the betaHMM package. It is used by
#' to wrap up the results table.
#'
#' @param SummarizedExperiment a \code{RangedSummarizedExperiment} of \code{dmr_identification} results.
#'
#'
#' @return a dmrResults object
#' @docType class
#' @rdname dmrResults
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom S4Vectors DataFrame
#' @export
dmrResults <- function(SummarizedExperiment) {
  se <- SummarizedExperiment
  if (!is(se, "RangedSummarizedExperiment")) {
    stop("'SummarizedExperiment' must be a RangedSummarizedExperiment object")
  }
  if(is.null(annotatedData)) annotatedData <- DataFrame(matrix(0, nrow=0, ncol=0))
  object <- new("dmrResults", se)
  return(object)
}

