#' @rdname betaHMM
#' @export
setGeneric("betaHMM", function(methylation_data, annotation_file, ...)
standardGeneric("betaHMM"),
signature = c("methylation_data", "annotation_file"))


#' @rdname packageHelpers
#' @export
setGeneric("annotatedData", function(object) standardGeneric("annotatedData"))

#' @rdname packageHelpers
#' @export
setGeneric("K", function(object) standardGeneric("K"))



#' @rdname packageHelpers
#' @export
setGeneric("N", function(object) standardGeneric("N"))

#' @rdname packageHelpers
#' @export
setGeneric("R", function(object) standardGeneric("R"))

#' @rdname packageHelpers
#' @export
setGeneric("A", function(object) standardGeneric("A"))

#' @rdname packageHelpers
#' @export
setGeneric("phi", function(object) standardGeneric("phi"))

#' @rdname packageHelpers
#' @export
setGeneric("treatment_group", function(object)
    standardGeneric("treatment_group"))



#' @rdname packageHelpers
#' @export
setGeneric("llk", function(object) standardGeneric("llk"))

#' @rdname packageHelpers
#' @export
setGeneric("tau", function(object) standardGeneric("tau"))

#' @rdname packageHelpers
#' @export
setGeneric("hidden_states", function(object) standardGeneric("hidden_states"))

#' @rdname packageHelpers
#' @export
setGeneric("chromosome_number", function(object)
    standardGeneric("chromosome_number"))


#' @rdname plot
#' @export
setGeneric("plot", function(x, ...) standardGeneric("plot"))



#' @rdname summary
#' @export
setGeneric("summary", function(object, ...) standardGeneric("summary"))

#' @rdname dmr_identification
#' @export
setGeneric("dmr_identification",
            function(dmc_identification_object, ...)
            standardGeneric("dmr_identification"),
            signature = c("dmc_identification_object"))

#' @rdname dmc_identification
#' @export
setGeneric("dmc_identification", function(betaHMM_object, ...)
    standardGeneric("dmc_identification"), signature = c("betaHMM_object"))

#' @rdname packageHelpers
#' @export
setGeneric("AUC", function(object) standardGeneric("AUC"))

#' @rdname packageHelpers
#' @export
setGeneric("uncertainty", function(object) standardGeneric("uncertainty"))


#' @rdname threshold_identification
#' @export
setGeneric("threshold_identification", function(object1, ...)
    standardGeneric("threshold_identification"), signature = c("object1"))


#' @rdname packageHelpers
#' @export
setGeneric("model_parameters", function(object)
    standardGeneric("model_parameters"))


#' @rdname packageHelpers
#' @export
setGeneric("threshold", function(object) standardGeneric("threshold"))


