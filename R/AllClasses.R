#' @rdname betaHMMResults
#' @export
#' @importFrom methods setClass is as new
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment



setClass("betaHMMResults",
            contains = "RangedSummarizedExperiment",
            representation = representation(annotatedData="DataFrame"))

setValidity("betaHMMResults",function(object){
    if(nrow(object@annotatedData)!=nrow(object@assays) &
       ncol(object@assays)!=object@metadata[[1]])
    {
        "Conditional posterior probability needs to be calculated for each
        of the CpG sites and K hidden states."
    }
})


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
#' @param SummarizedExperiment a \code{RangedSummarizedExperiment} of
#' \code{betaHMM} results
#' @param annotatedData The annotated data passed as an input argument to
#' the betaHMM package.
#'
#'
#' @return a betaHMMResults object
#' @docType class
#' @rdname betaHMMResults
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom S4Vectors DataFrame
#' @export
#' @example inst/examples/betaHMM_package.R
betaHMMResults <- function(SummarizedExperiment,annotatedData) {
    se_obj<-se(SummarizedExperiment)
    annotatedData<-anno_data(annotatedData)
    object <- new("betaHMMResults", se_obj, annotatedData=annotatedData)
    return(object)
    }

#' @rdname dmrResults
#' @export
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment

setClass("dmrResults",
            contains = "RangedSummarizedExperiment")


setValidity("dmrResults",function(object){
    x<-as.data.frame(assay(object))
    if(any(is.na(x[,seq_len(2)])))
    {
        "Starting and ending CpG sites of a DMR cannot be empty."
    }
})


#' dmrResults object and constructor
#'
#' \code{dmrResults} is a subclass of \code{RangedSummarizedExperiment},
#' used to store the DMRs identified.
#'
#' This constructor function would not typically be used by "end users".
#' This simple class extends the \code{RangedSummarizedExperiment} class of the
#' SummarizedExperiment package
#' to allow other packages to write methods for results
#' objects from the
#' \code{\link[betaHMM:dmr_identification]{dmr_identification}} function.
#' It is used by
#' to wrap up the results table.
#'
#' @param SummarizedExperiment a \code{dmrResults}
#' results.
#'
#'
#' @return a  \code{\link[betaHMM:dmrResults]{dmrResults}} object
#' @docType class
#' @rdname dmrResults
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom S4Vectors DataFrame
#' @export
#' @example inst/examples/betaHMM_package.R
dmrResults <- function(SummarizedExperiment) {
    se_obj<-se(SummarizedExperiment)
    object <- new("dmrResults", se_obj)
    return(object)
}




#' @rdname dmcResults
#' @export
#' @importFrom methods setClass
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment


setClass("dmcResults",
            contains = "RangedSummarizedExperiment")


setValidity("dmcResults",function(object){
    if(nrow(object@assays) !=length(unlist(object@metadata[[6]])))
    {
        "Uncertainty needs to be calculated for each CpG site analysed."
    }
})


#' dmcResults object and constructor
#'
#' \code{dmcResults} is a subclass of \code{RangedSummarizedExperiment},
#' used to store the DMCs identified.
#'
#' This constructor function would not typically be used by "end users".
#' This simple class extends the \code{RangedSummarizedExperiment} class of the
#' SummarizedExperiment package
#' to allow other packages to write methods for results
#' objects from the
#' \code{\link[betaHMM:dmc_identification]{dmc_identification}} function.
#' It is used by
#' to wrap up the results table.
#'
#' @param SummarizedExperiment a \code{RangedSummarizedExperiment} of
#' \code{dmcResults} results.
#'
#'
#' @return a \code{\link[betaHMM:dmcResults]{dmcResults}} object
#' @docType class
#' @rdname dmcResults
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom S4Vectors DataFrame
#' @export
#' @example inst/examples/betaHMM_package.R
dmcResults <- function(SummarizedExperiment) {
    se_obj<-se(SummarizedExperiment)
    object <- new("dmcResults", se_obj)
    return(object)
}


## Threshold definition

#' @rdname threshold_Results
#' @export
#' @importFrom methods setClass
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment


setClass("threshold_Results",
            contains = "RangedSummarizedExperiment",
            representation = representation(
            annotatedData="DataFrame"

            ))


setValidity("threshold_Results",function(object){
    if(nrow(object@annotatedData)!=nrow(object@assays) &
       ncol(object@assays)!=object@metadata[[2]])
    {
        "Conditional posterior probability needs to be calculated for each
        of the CpG sites and K hidden states."
    }
})


#' threshold_Results object and constructor
#'
#' \code{threshold_Results} is a subclass of \code{RangedSummarizedExperiment},
#' used to store the threshold_identification results as well as the annotated
#'  data useful for plotting.
#'
#' This constructor function would not typically be used by "end users".
#' This simple class extends the \code{RangedSummarizedExperiment} class of the
#' SummarizedExperiment package
#' to allow other packages to write methods for results
#' objects from the
#' \code{\link[betaHMM:threshold_identification]{threshold_identification}}
#' function. It is used by
#' to wrap up the results table.
#'
#' @param SummarizedExperiment a \code{RangedSummarizedExperiment} of
#' \code{threshold_Results} object.
#' @param annotatedData The annotated data passed as an input argument to
#' the \code{threshold_identification} function.
#'
#'
#' @return a \code{\link[betaHMM:threshold_Results]{threshold_Results}} object
#' @docType class
#' @rdname threshold_Results
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom S4Vectors DataFrame
#' @export
#' @examples
#' ## Use simulated data for the betaHMM workflow example
#' set.seed(12345)
#' library(betaHMM)
#'
#' ## read files
#' data(sample_methylation_file)
#' head(sample_methylation_file)
#' data(sample_annotation_file)
#' head(sample_annotation_file)
#' ##merge data
#' df=merge(sample_annotation_file[,c('IlmnID','CHR','MAPINFO')],
#' sample_methylation_file,by='IlmnID')
#'
#' ## sort data
#' df=df[order(df$CHR,df$MAPINFO),]
#' thr_out=threshold_identification(df[,c(1,4:7)],package_workflow=TRUE,M=3,4,
#' parameter_estimation_only=TRUE,seed=12345)
#'

threshold_Results <- function(SummarizedExperiment,
                            annotatedData) {
    se_obj<-se(SummarizedExperiment)
    annotatedData<-anno_data(annotatedData)
    object <- new("threshold_Results", se_obj, annotatedData=annotatedData)
    return(object)
}
