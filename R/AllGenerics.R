#' @rdname betaHMM
#' @export
setGeneric("betaHMM", function(object,...) standardGeneric("betaHMM"),
           signature=c("object"))


#' @rdname betaHMMHelpers
#' @export
setGeneric("annotatedData", function(object) standardGeneric("annotatedData"))

#' @rdname betaHMMHelpers
#' @export
setGeneric("K", function(object) standardGeneric("K"))



#' @rdname betaHMMHelpers
#' @export
setGeneric("N", function(object) standardGeneric("N"))

#' @rdname betaHMMHelpers
#' @export
setGeneric("R", function(object) standardGeneric("R"))

#' @rdname betaHMMHelpers
#' @export
setGeneric("A", function(object) standardGeneric("A"))

#' @rdname betaHMMHelpers
#' @export
setGeneric("phi", function(object) standardGeneric("phi"))

#' @rdname betaHMMHelpers
#' @export
setGeneric("llk", function(object) standardGeneric("llk"))

#' @rdname betaHMMHelpers
#' @export
setGeneric("tau", function(object) standardGeneric("tau"))

#' @rdname betaHMMHelpers
#' @export
setGeneric("hidden_states", function(object) standardGeneric("hidden_states"))


#' @rdname plot
#' @export
setGeneric("plot", function(x, ...) standardGeneric("plot"))

#' @rdname dmr_identification
#' @export
setGeneric("dmr_identification", function(object,...) standardGeneric("dmr_identification"),
           signature=c("object"))
