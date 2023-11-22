#' HMM for beta valued DNA data
#'
#' This is the primary user interface for the \code{betaHMM} function
#' Generic S4 methods are implemented to eatimate the parameters of a
#' homogeneous hidden Markov model for the beta valued DNA methylation data.
#' The supported classes are \code{matrix}, \code{data.frame},
#' \code{RangedSummarizedExperiment} and \code{GRanges}. The output of
#' \code{betaHMM} method is an
#' S4 object of class \code{betaHMMResults}.
#'
#' @inheritParams betaHMMrun
#' @param methylation_data A dataframe of dimension \eqn{(C \times (N \times R)
#' )+1} containing methylation values for \eqn{C} CpG sites from \eqn{R}
#' treatment groups each having \eqn{N} DNA samples and the IlmnID
#' for each CpG site. Maybe provided as a matrix or data.frame, GRanges
#' or RangedSummarizedExperiment object.
#' @param annotation_file A dataframe containing the EPIC methylation
#' annotation file. Maybe provided as a matrix or data.frame, GRanges
#' or RangedSummarizedExperiment object.
#'
#' @return An S4 object of class \code{betaHMMResults}, where conditional
#' probabilities of each CpG site belonging to a hidden state is stored as a
#' SimpleList of assay data, and the corresponding estimated model parameters,
#' log-likelihood values, and most probable hidden state sequence for each
#' chromosome are stored as metadata.
#'
#' @aliases
#' betaHMM
#' betaHMM-methods
#' betaHMM,matrix,matrix-method
#' betaHMM,data.frame,data.frame-method
#' betaHMM,RangedSummarizedExperiment,RangedSummarizedExperiment.-method
#' betaHMM,matrix,data.frame-method
#' betaHMM,matrix,RangedSummarizedExperiment.-method
#' betaHMM,data.frame,matrix-method
#' betaHMM,data.frame,RangedSummarizedExperiment.-method
#' betaHMM,RangedSummarizedExperiment,data.frame-method
#' betaHMM,RangedSummarizedExperiment,matrix-method
#' betaHMM,GRanges,data.frame-method
#' betaHMM,GRanges,matrix-method
#' betaHMM,GRanges,RangedSummarizedExperiment-method
#' betaHMM,GRanges,GRanges-method
#' betaHMM,data.frame,GRanges-method
#' betaHMM,RangedSummarizedExperiment,GRanges-method
#' betaHMM,matrix,GRanges-method
#' @author Koyel Majumdar
#' @export
#' @keywords methods
#' @rdname betaHMM
#' @docType methods
#' @importFrom methods setMethod is as new
#' @importFrom stats p.adjust
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom GenomicRanges mcols
#' @example inst/examples/betaHMM_package.R
#'

setMethod("betaHMM", signature = signature(methylation_data = "matrix",
annotation_file = "matrix"), definition = function(methylation_data,
annotation_file, M = 3, N = 4, R = 2, treatment_group = NULL,
parallel_process = FALSE, seed = NULL,iterations=100,...) {
data1 <- as.data.frame(methylation_data)
data2 <- as.data.frame(annotation_file)
run <- betaHMMrun(methylation_data = data1, annotation_file = data2,
                    M = M, N = N, R = R, treatment_group = treatment_group,
                    parallel_process = parallel_process,
                    seed = seed,iterations=iterations, ...)
return(run)
})



##############################################################################
#' @rdname betaHMM
#' @export
setMethod("betaHMM", signature = signature(methylation_data = "data.frame",
annotation_file = "data.frame"), definition = function(methylation_data,
annotation_file, M = 3, N = 4, R = 2, treatment_group = NULL,
parallel_process = FALSE,
seed = NULL,iterations=100,
...) {
data1 <- methylation_data
data2 <- annotation_file
run <- betaHMMrun(methylation_data = data1, annotation_file = data2,
                    M = M, N = N, R = R, treatment_group = treatment_group,
                    parallel_process = parallel_process,
                    seed = seed,iterations=iterations, ...)
return(run)
})


##############################################################################
#' @rdname betaHMM
#' @export
setMethod("betaHMM",
signature=signature(methylation_data="RangedSummarizedExperiment",
annotation_file = "RangedSummarizedExperiment"),
definition = function(methylation_data, annotation_file,
M = 3, N = 4, R = 2, treatment_group = NULL, parallel_process = FALSE,
seed = NULL,iterations=100, ...) {
data1 <- assay(methylation_data)
data2 <- assay(annotation_file)
run <- betaHMMrun(methylation_data = data1, annotation_file = data2,
                    M = M, N = N, R = R, treatment_group = treatment_group,
                    parallel_process = parallel_process,
                    seed = seed,iterations=iterations, ...)
return(run)
})

##############################################################################
#' @rdname betaHMM
#' @export
setMethod("betaHMM",signature=signature(methylation_data="GRanges",
annotation_file = "GRanges"),
definition = function(methylation_data, annotation_file,
M = 3, N = 4, R = 2, treatment_group = NULL, parallel_process = FALSE,
seed = NULL,iterations=100, ...) {
data1 <- as.data.frame(mcols(methylation_data))
data2 <- as.data.frame(mcols(annotation_file))
run <- betaHMMrun(methylation_data = data1, annotation_file = data2,
M = M, N = N, R = R, treatment_group = treatment_group,
parallel_process = parallel_process,
seed = seed,iterations=iterations, ...)
return(run)
})

##############################################################################
#' @rdname betaHMM
#' @export
setMethod("betaHMM",
signature=signature(methylation_data="RangedSummarizedExperiment",
annotation_file = "matrix"), definition = function(methylation_data,
annotation_file, M = 3, N = 4, R = 2, treatment_group = NULL,
parallel_process = FALSE,
seed = NULL,iterations=100,
...) {
data1 <- assay(methylation_data)
data2 <- as.data.frame(annotation_file)
run <- betaHMMrun(methylation_data = data1, annotation_file = data2,
                    M = M, N = N, R = R, treatment_group = treatment_group,
                    parallel_process = parallel_process,
                    seed = seed,iterations=iterations, ...)
return(run)
})
##############################################################################
#' @rdname betaHMM
#' @export
setMethod("betaHMM",
signature=signature(methylation_data="RangedSummarizedExperiment",
annotation_file = "data.frame"), definition = function(methylation_data,
annotation_file, M = 3, N = 4, R = 2, treatment_group = NULL,
parallel_process = FALSE,
seed = NULL,iterations=100,
...) {
data1 <- assay(methylation_data)
data2 <- annotation_file
run <- betaHMMrun(methylation_data = data1, annotation_file = data2,
                    M = M, N = N, R = R, treatment_group = treatment_group,
                    parallel_process = parallel_process,
                    seed = seed,iterations=iterations, ...)
return(run)
})
##############################################################################
#' @rdname betaHMM
#' @export
setMethod("betaHMM",
signature=signature(methylation_data="RangedSummarizedExperiment",
annotation_file = "GRanges"), definition = function(methylation_data,
annotation_file, M = 3, N = 4, R = 2, treatment_group = NULL,
parallel_process = FALSE,
seed = NULL,iterations=100,...) {
data1 <- assay(methylation_data)
data2 <- as.data.frame(mcols(annotation_file))
run <- betaHMMrun(methylation_data = data1, annotation_file = data2,
M = M, N = N, R = R, treatment_group = treatment_group,
parallel_process = parallel_process,
seed = seed,iterations=iterations, ...)
return(run)
})

##############################################################################
#' @rdname betaHMM
#' @export
setMethod("betaHMM", signature = signature(methylation_data = "matrix",
annotation_file = "RangedSummarizedExperiment"),
definition = function(methylation_data, annotation_file,
M = 3, N = 4, R = 2, treatment_group = NULL, parallel_process = FALSE,
seed = NULL,iterations=100, ...) {
data1 <- as.data.frame(methylation_data)
data2 <- assay(annotation_file)
run <- betaHMMrun(methylation_data = data1, annotation_file = data2,
                    M = M, N = N, R = R, treatment_group = treatment_group,
                    parallel_process = parallel_process,
                    seed = seed,iterations=iterations, ...)
return(run)
})

##############################################################################
#' @rdname betaHMM
#' @export
setMethod("betaHMM", signature = signature(methylation_data = "data.frame",
annotation_file = "RangedSummarizedExperiment"),
definition = function(methylation_data, annotation_file,
M = 3, N = 4, R = 2, treatment_group = NULL, parallel_process = FALSE,
seed = NULL,iterations=100,...) {
data1 <- methylation_data
data2 <- assay(annotation_file)
run <- betaHMMrun(methylation_data = data1, annotation_file = data2,
                    M = M, N = N, R = R, treatment_group = treatment_group,
                    parallel_process = parallel_process,
                    seed = seed,iterations=iterations, ...)
return(run)
})
##############################################################################
#' @rdname betaHMM
#' @export
setMethod("betaHMM", signature = signature(methylation_data = "data.frame",
annotation_file = "GRanges"),
definition = function(methylation_data, annotation_file,
M = 3, N = 4, R = 2, treatment_group = NULL, parallel_process = FALSE,
seed = NULL,iterations=100,...) {
data1 <- methylation_data
data2 <- as.data.frame(mcols(annotation_file))
run <- betaHMMrun(methylation_data = data1, annotation_file = data2,
M = M, N = N, R = R, treatment_group = treatment_group,
parallel_process = parallel_process,
seed = seed,iterations=iterations, ...)
return(run)
})
##############################################################################
#' @rdname betaHMM
#' @export
setMethod("betaHMM", signature = signature(methylation_data = "matrix",
annotation_file = "data.frame"), definition = function(methylation_data,
annotation_file, M = 3, N = 4, R = 2,treatment_group = NULL,
parallel_process = FALSE,
seed = NULL,iterations=100,
...) {
data1 <- as.data.frame(methylation_data)
data2 <- annotation_file
run <- betaHMMrun(methylation_data = data1, annotation_file = data2,
                    M = M, N = N, R = R, treatment_group = treatment_group,
                    parallel_process = parallel_process,
                    seed = seed,iterations=iterations, ...)
return(run)
})

##############################################################################
#' @rdname betaHMM
#' @export
setMethod("betaHMM", signature = signature(methylation_data = "matrix",
annotation_file = "GRanges"), definition = function(methylation_data,
annotation_file, M = 3, N = 4, R = 2,treatment_group = NULL,
parallel_process = FALSE,
seed = NULL,iterations=100,
...) {
data1 <- as.data.frame(methylation_data)
data2 <- as.data.frame(mcols(annotation_file))
run <- betaHMMrun(methylation_data = data1, annotation_file = data2,
M = M, N = N, R = R, treatment_group = treatment_group,
parallel_process = parallel_process,
seed = seed,iterations=iterations, ...)
return(run)
})

##############################################################################
#' @rdname betaHMM
#' @export
setMethod("betaHMM",
signature=signature(methylation_data="GRanges",
annotation_file = "matrix"),
definition = function(methylation_data, annotation_file,
M = 3, N = 4, R = 2, treatment_group = NULL, parallel_process = FALSE,
seed = NULL,iterations=100, ...) {
data1 <- as.data.frame(mcols(methylation_data))
data2 <- as.data.frame(annotation_file)
run <- betaHMMrun(methylation_data = data1, annotation_file = data2,
M = M, N = N, R = R, treatment_group = treatment_group,
parallel_process = parallel_process,
seed = seed,iterations=iterations, ...)
return(run)
})
##############################################################################
#' @rdname betaHMM
#' @export
setMethod("betaHMM",
signature=signature(methylation_data="GRanges",
annotation_file = "data.frame"),
definition = function(methylation_data, annotation_file,
M = 3, N = 4, R = 2, treatment_group = NULL, parallel_process = FALSE,
seed = NULL,iterations=100, ...) {
data1 <- as.data.frame(mcols(methylation_data))
data2 <- annotation_file
run <- betaHMMrun(methylation_data = data1, annotation_file = data2,
M = M, N = N, R = R, treatment_group = treatment_group,
parallel_process = parallel_process,
seed = seed,iterations=iterations, ...)
return(run)
})
##############################################################################
#' @rdname betaHMM
#' @export
setMethod("betaHMM",
signature=signature(methylation_data="GRanges",
annotation_file = "RangedSummarizedExperiment"),
definition = function(methylation_data, annotation_file,
M = 3, N = 4, R = 2, treatment_group = NULL, parallel_process = FALSE,
seed = NULL,iterations=100, ...) {
data1 <- as.data.frame(mcols(methylation_data))
data2 <- assay(annotation_file)
run <- betaHMMrun(methylation_data = data1, annotation_file = data2,
M = M, N = N, R = R, treatment_group = treatment_group,
parallel_process = parallel_process,
seed = seed,iterations=iterations, ...)
return(run)
})
##############################################################################
#' @rdname betaHMM
#' @export
setMethod("betaHMM", signature = signature(methylation_data = "data.frame",
annotation_file = "matrix"), definition = function(methylation_data,
annotation_file, M = 3, N = 4, R = 2, treatment_group = NULL,
parallel_process = FALSE,
seed = NULL,iterations=100,
...) {
data1 <- methylation_data
data2 <- as.data.frame(annotation_file)
run <- betaHMMrun(methylation_data = data1, annotation_file = data2,
M = M, N = N, R = R, treatment_group = treatment_group,
parallel_process = parallel_process,
seed = seed,iterations=iterations, ...)
return(run)
})

##############################################################################
#' Accessors for the betaHMM package.
#'
#' The accessor methods for accessing the betaHMMResults/ dmcResults/
#' dmrResults/ threshold_Results metadata.
#'
#' @docType methods
#' @rdname packageHelpers
#' @aliases
#' annotatedData
#' annotatedData,betaHMMResults-method
#' annotatedData,threshold_Results-method
#' K
#' K,RangedSummarizedExperiment-method
#' K,betaHMMResults-method
#' K,dmcResults-method
#' K,threshold_Results-method
#' K,NULL-method
#' N
#' N,RangedSummarizedExperiment-method
#' N,betaHMMResults-method
#' N,dmcResults-method
#' N,NULL-method
#' R
#' R,RangedSummarizedExperiment-method
#' R,betaHMMResults-method
#' R,dmcResults-method
#' R,NULL-method
#' A
#' A,RangedSummarizedExperiment-method
#' A,betaHMMResults-method
#' A,NULL-method
#' tau
#' tau,RangedSummarizedExperiment-method
#' tau,betaHMMResults-method
#' tau,NULL-method
#' phi
#' phi,RangedSummarizedExperiment-method
#' phi,betaHMMResults-method
#' phi,threshold_Results-method
#' phi,NULL-method
#' llk
#' llk,RangedSummarizedExperiment-method
#' llk,betaHMMResults-method
#' llk,NULL-method
#' treatment_group
#' treatment_group,RangedSummarizedExperiment-method
#' treatment_group,betaHMMResults-method
#' treatment_group,dmcResults-method
#' treatment_group,NULL-method
#' hidden_states
#' hidden_states,RangedSummarizedExperiment-method
#' hidden_states,betaHMMResults-method
#' hidden_states,threshold_Results-method
#' hidden_states,NULL-method
#' model_parameters
#' model_parameters,RangedSummarizedExperiment-method
#' model_parameters,threshold_Results-method
#' model_parameters,NULL-method
#' threshold
#' threshold,RangedSummarizedExperiment-method
#' threshold,threshold_Results-method
#' threshold,NULL-method
#' @param object a \code{betaHMMResults}/ \code{dmcResults}/
#' \code{threshold_Results} object.
#' @return Output varies depending on the method.
#' @export
#' @example inst/examples/betaHMM_package.R
setMethod("annotatedData", "betaHMMResults",
            function(object) object@annotatedData)

#' @rdname packageHelpers
#' @export
setMethod("annotatedData", "threshold_Results",
            function(object) object@annotatedData)

### K
#' @rdname packageHelpers
#' @export
setMethod("K", "RangedSummarizedExperiment",
            function(object) metadata(object)$K)

#' @rdname packageHelpers
#' @export
setMethod("K", "betaHMMResults", function(object) metadata(object)$K)
#' @rdname packageHelpers
#' @export
setMethod("K", "dmcResults", function(object) metadata(object)$K)

#' @rdname packageHelpers
#' @export
setMethod("K", "threshold_Results", function(object) metadata(object)$K)

#' @rdname packageHelpers
#' @export
setMethod("K", "NULL", function(object) NA)

### N
#' @rdname packageHelpers
#' @export
setMethod("N", "RangedSummarizedExperiment",
            function(object) metadata(object)$N)

#' @rdname packageHelpers
#' @export
setMethod("N", "betaHMMResults", function(object) metadata(object)$N)

#' @rdname packageHelpers
#' @export
setMethod("N", "dmcResults", function(object) metadata(object)$N)



#' @rdname packageHelpers
#' @export
setMethod("N", "NULL", function(object) NA)

### R
#' @rdname packageHelpers
#' @export
setMethod("R", "RangedSummarizedExperiment",
            function(object) metadata(object)$R)

#' @rdname packageHelpers
#' @export
setMethod("R", "betaHMMResults", function(object) metadata(object)$R)


#' @rdname packageHelpers
#' @export
setMethod("R", "dmcResults", function(object) metadata(object)$R)



#' @rdname packageHelpers
#' @export
setMethod("R", "NULL", function(object) NA)

### A
#' @rdname packageHelpers
#' @export
setMethod("A", "RangedSummarizedExperiment",
            function(object) metadata(object)$A)

#' @rdname packageHelpers
#' @export
setMethod("A", "betaHMMResults", function(object) metadata(object)$A)

#' @rdname packageHelpers
#' @export
setMethod("A", "NULL", function(object) NA)


### tau
#' @rdname packageHelpers
#' @export
setMethod("tau", "RangedSummarizedExperiment",
            function(object) metadata(object)$tau)

#' @rdname packageHelpers
#' @export
setMethod("tau", "betaHMMResults", function(object) metadata(object)$tau)

#' @rdname packageHelpers
#' @export
setMethod("tau", "NULL", function(object) NA)


### treatment_group
#' @rdname packageHelpers
#' @export
setMethod("treatment_group", "RangedSummarizedExperiment",
    function(object) metadata(object)$treatment_group)

#' @rdname packageHelpers
#' @export
setMethod("treatment_group", "betaHMMResults",
            function(object) metadata(object)$treatment_group)

#' @rdname packageHelpers
#' @export
setMethod("treatment_group", "dmcResults",
            function(object) metadata(object)$treatment_group)


#' @rdname packageHelpers
#' @export
setMethod("treatment_group", "NULL", function(object) NA)


### llk
#' @rdname packageHelpers
#' @export
setMethod("llk", "RangedSummarizedExperiment",
            function(object) metadata(object)$llk)

#' @rdname packageHelpers
#' @export
setMethod("llk", "betaHMMResults", function(object) metadata(object)$llk)

#' @rdname packageHelpers
#' @export
setMethod("llk", "NULL", function(object) NA)

### phi
#' @rdname packageHelpers
#' @export
setMethod("phi", "RangedSummarizedExperiment",
            function(object) metadata(object)$phi)

#' @rdname packageHelpers
#' @export
setMethod("phi", "betaHMMResults", function(object) metadata(object)$phi)

#' @rdname packageHelpers
#' @export
setMethod("phi", "threshold_Results", function(object) metadata(object)$phi)

#' @rdname packageHelpers
#' @export
setMethod("phi", "NULL", function(object) NA)


### hidden_states
#' @rdname packageHelpers
#' @export
setMethod("hidden_states", "RangedSummarizedExperiment",
    function(object) metadata(object)$hidden_states)

#' @rdname packageHelpers
#' @export
setMethod("hidden_states", "betaHMMResults",
            function(object) metadata(object)$hidden_states)

#' @rdname packageHelpers
#' @export
setMethod("hidden_states", "threshold_Results",
            function(object) metadata(object)$hidden_states)

#' @rdname packageHelpers
#' @export
setMethod("hidden_states", "NULL", function(object) NA)


### chromosome_number
#' @rdname packageHelpers
#' @export
setMethod("chromosome_number", "RangedSummarizedExperiment",
    function(object) metadata(object)$chromosome_number)

#' @rdname packageHelpers
#' @export
setMethod("chromosome_number", "betaHMMResults",
            function(object) metadata(object)$chromosome_number)

#' @rdname packageHelpers
#' @export
setMethod("chromosome_number", "dmrResults",
            function(object) metadata(object)$chromosome_number)

#' @rdname packageHelpers
#' @export
setMethod("chromosome_number", "NULL", function(object) NA)

#' @rdname packageHelpers
#' @export
setMethod("AUC", "RangedSummarizedExperiment",
            function(object) metadata(object)$AUC)

#' @rdname packageHelpers
#' @export
setMethod("AUC", "dmcResults", function(object) metadata(object)$AUC)

#' @rdname packageHelpers
#' @export
setMethod("AUC", "NULL", function(object) NA)

### conditional probability values
#' @rdname packageHelpers
#' @export
setMethod("uncertainty", "RangedSummarizedExperiment",
    function(object) metadata(object)$uncertainty)

#' @rdname packageHelpers
#' @export
setMethod("uncertainty", "dmcResults",
            function(object) metadata(object)$uncertainty)

#' @rdname packageHelpers
#' @export
setMethod("uncertainty", "NULL", function(object) NA)


#' @rdname packageHelpers
#' @export
setMethod("model_parameters", "RangedSummarizedExperiment",
    function(object) metadata(object)$model_parameters)

#' @rdname packageHelpers
#' @export
setMethod("model_parameters", "threshold_Results",
    function(object) metadata(object)$model_parameters)

#' @rdname packageHelpers
#' @export
setMethod("model_parameters", "NULL", function(object) NA)

### threshold
#' @rdname packageHelpers
#' @export
setMethod("threshold", "RangedSummarizedExperiment",
    function(object) metadata(object)$threshold)

#' @rdname packageHelpers
#' @export
setMethod("threshold", "threshold_Results",
            function(object) metadata(object)$threshold)

#' @rdname packageHelpers
#' @export
setMethod("threshold", "NULL", function(object) NA)


### DMR methods
#' @title DMR identification from DMCs identified
#' @description This is the primary user interface for the
#' \code{\link[betaHMM:dmr_identification]{dmr_identification}} function.
#' Generic S4 methods are implemented to identify the DMRs from the DMCs
#' identified in each chromosome. The supported classes are \code{data.frame}
#' and \code{\link[betaHMM:dmcResults]{dmcResults}} object. The output is an
#' S4 object of class \code{\link[betaHMM:dmrResults]{dmrResults}}.
#' @rdname dmr_identification
#' @aliases
#' dmr_identification
#' dmr_identification-methods
#' dmr_identification,dmcResults-method
#' @export
#' @seealso \code{\link{betaHMM}}
#' @inheritParams dmr_identification_run
#' @param dmc_identification_object a
#' \code{\link[betaHMM:dmcResults]{dmcResults}} object or the assay data
#' from the dmcResults.
#' @return An S4 object of class \code{\link[betaHMM:dmrResults]{dmrResults}}
#' where the CpG site information for each DMR is stored as a SimpleList
#' of assay data and the chromosomes analysed by the model is stored as
#' the metadata.
#' @importFrom stats complete.cases
#' @example inst/examples/betaHMM_package.R
#'
setMethod(f = "dmr_identification",
            signature(dmc_identification_object = "dmcResults"),
            definition = function(dmc_identification_object, DMC_count = 2,
                                ...) {
            object_df <- assay(dmc_identification_object)
            dmr_df <-dmr_identification_run(dmc_identification_object =
                                            object_df,
                                            DMC_count = DMC_count, ...)
        return(dmr_df)
    })



##############################################################################
#' @rdname dmr_identification
#' @export
setMethod(f="dmr_identification",
            signature = signature(dmc_identification_object = "matrix"),
            definition = function(dmc_identification_object, DMC_count = 2,
                                ...) {
            dmc_identification_object <-
                as.data.frame(dmc_identification_object)
            dmr_df <- dmr_identification_run(dmc_identification_object =
                                                dmc_identification_object,
                                                DMC_count = DMC_count, ...)
        return(dmr_df)
    })



##############################################################################
#' @rdname dmr_identification
#' @export
setMethod(f="dmr_identification",
            signature(dmc_identification_object = "data.frame"),
    definition = function(dmc_identification_object, DMC_count = 2,
                                ...) {
            dmr_df <- dmr_identification_run(dmc_identification_object =
                                                dmc_identification_object,
                                                DMC_count = DMC_count, ...)
        return(dmr_df)
    })

### DMC methods
#' @title DMC identification from estimated betaHMM model parameters
#' @description This is the primary user interface for the
#' \code{\link[betaHMM:dmc_identification]{dmc_identification}} function.
#' Generic S4 methods are implemented to identify the DMCs from the estimated
#' betaHMM model parameters for each chromosome. The supported class is a
#'  \code{\link[betaHMM:betaHMMResults]{betaHMMResults}} object. The output
#'  is an S4 object of class of \code{\link[betaHMM:dmcResults]{dmcResults}}.
#' @rdname dmc_identification
#' @aliases
#' dmc_identification
#' dmc_identification-methods
#' dmc_identification,betaHMMResults-method
#' @export
#' @seealso \code{\link{betaHMM}}
#' @inheritParams dmc_identification_run
#' @param betaHMM_object An S4 object of class
#'  \code{\link[betaHMM:betaHMMResults]{betaHMMResults}}
#' @return An S4 object of class \code{\link[betaHMM:dmcResults]{dmcResults}}.
#' @importFrom stats complete.cases
#' @example inst/examples/betaHMM_package.R
#'
setMethod(f = "dmc_identification",
            signature(betaHMM_object = "betaHMMResults"),
    definition = function(betaHMM_object, AUC_threshold = 0.8,
        uncertainty_threshold = 0.2, ...) {

        dmc_df <- dmc_identification_run(betaHMM_object = betaHMM_object,
            AUC_threshold = AUC_threshold,
            uncertainty_threshold = uncertainty_threshold,
            ...)
        return(dmc_df)
    })



## Threshold indentification function
#' HMM for beta valued DNA data for a single treatment condition
#'
#' The supported classes are \code{matrix} and \code{data.frame}.
#' The output of \code{threshold_identification} is an S4 object of
#' class \code{threshold_Results}.
#'
#' @inheritParams threshold_identification_run
#' @param object1 Methylation data and IlmnID. Maybe provided as a matrix
#' or dataframe.
#' @return
#' An S4 object of class \code{threshold_Results}.
#'
#' @aliases
#' threshold_identification
#' threshold_identification-methods
#' threshold_identification,matrix-method
#' threshold_identification,data.frame-method
#' @export
#' @rdname threshold_identification
#' @importFrom methods setMethod is as new
#' @importFrom stats p.adjust
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

setMethod("threshold_identification",
            signature = signature(object1 = "matrix"),
    definition = function(object1, package_workflow = TRUE,
        annotation_file = NULL, M = 3, N = 4,
        parameter_estimation_only = FALSE,
        seed = NULL, ...) {
        data1 <- as.data.frame(object1)
        run <- threshold_identification_run(data = data1,
            package_workflow = package_workflow,
            annotation_file = annotation_file,
            M = M, N = N,
            parameter_estimation_only = parameter_estimation_only,
            seed = seed, ...)
        return(run)
    })



##############################################################################
#' @rdname threshold_identification
#' @export
setMethod("threshold_identification",
            signature = signature(object1 = "data.frame"),
    definition = function(object1, package_workflow = TRUE,
        annotation_file = NULL, M = 3, N = 4,
        parameter_estimation_only = FALSE,
        seed = NULL, ...) {
        data1 <- object1
        run <- threshold_identification_run(data = data1,
            package_workflow = package_workflow,
            annotation_file = annotation_file,
            M = M, N = N,
            parameter_estimation_only = parameter_estimation_only,
            seed = seed, ...)
        return(run)
    })

##############################################################################

