#' HMM for beta valued DNA data
#'
#' This is the primary user interface for the \code{betaHMM} package.
#' Generic S4 methods are implemented to perform a homogeneous hidden Markov model for the beta valued DNA
#'              methylation data.
#' The supported classes are \code{matrix} and \code{data.frame}.
#' The output of \code{betaHMM} is an S4 object of class \code{betaHMMResults}.
#'
#' @inheritParams betaHMMrun
#' @param object Annotated data to be clustered. Maybe provided as a matrix or dataframe.
#'
#' @return
#' An S4 object of class \code{betaHMMResults}.
#'
#' @aliases
#' betaHMM
#' betaHMM-methods
#' betaHMM,matrix-method
#' betaHMM,data.frame-method
#'
#' @author Koyel Majumdar
#' @export
#' @keywords methods
#' @rdname betaHMM
#' @docType methods
#' @import methods
#' @importFrom stats p.adjust
#'

setMethod("betaHMM",
          signature=signature(object="matrix"),
          definition=function(object, M=3,N=4,R=2,seed=NULL,...)
          {
            data <- as.data.frame(object)
            #arg.user <- list(...)
              run <- betaHMMrun(data=data, M=M,N=N,R=R,seed=seed,...)
            return(run)
          })



#########################################################################################
#' @rdname betaHMM
#' @export
setMethod("betaHMM", signature=signature(object="data.frame"),
          definition=function(object,  M=3,N=4,R=2,seed=NULL,...)
          {
            data <- object
            run <- betaHMMrun(data=data, M=M,N=N,R=R,seed=seed,...)
            return(run)
          })

#########################################################################################
#' Accessors for the assigned cluster labels of a coseqResults object.
#'
#' The counts slot holds the count data as a matrix of non-negative integer
#' count values, one row for each observational unit (gene or the like), and one
#' column for each sample.
#'
#' @docType methods
#' @aliases
#' annotatedData
#' annotatedData,betaHMMResults-method
#' K
#' K,RangedSummarizedExperiment-method
#' K,betaHMMResults-method
#' K,NULL-method
#' N
#' N,RangedSummarizedExperiment-method
#' N,betaHMMResults-method
#' N,NULL-method
#' R
#' R,RangedSummarizedExperiment-method
#' R,betaHMMResults-method
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
#' phi,NULL-method
#' llk
#' llk,RangedSummarizedExperiment-method
#' llk,betaHMMResults-method
#' llk,NULL-method
#' hidden_states
#' hidden_states,RangedSummarizedExperiment-method
#' hidden_states,betaHMMResults-method
#' hidden_states,NULL-method
#' @rdname betaHMMHelpers
#' @param object a \code{betaHMMResults} object.
#' @export
setMethod("annotatedData", "betaHMMResults", function(object) object@annotatedData)

### K
#' @rdname betaHMMHelpers
#' @export
setMethod("K", "RangedSummarizedExperiment", function(object) metadata(object)$K)

#' @rdname betaHMMHelpers
#' @export
setMethod("K", "betaHMMResults", function(object) metadata(object)$K)

#' @rdname betaHMMHelpers
#' @export
setMethod("K", "NULL", function(object) NA)

### N
#' @rdname betaHMMHelpers
#' @export
setMethod("N", "RangedSummarizedExperiment", function(object) metadata(object)$N)

#' @rdname betaHMMHelpers
#' @export
setMethod("N", "betaHMMResults", function(object) metadata(object)$N)

#' @rdname betaHMMHelpers
#' @export
setMethod("N", "NULL", function(object) NA)

### R
#' @rdname betaHMMHelpers
#' @export
setMethod("R", "RangedSummarizedExperiment", function(object) metadata(object)$R)

#' @rdname betaHMMHelpers
#' @export
setMethod("R", "betaHMMResults", function(object) metadata(object)$R)

#' @rdname betaHMMHelpers
#' @export
setMethod("R", "NULL", function(object) NA)

### A
#' @rdname betaHMMHelpers
#' @export
setMethod("A", "RangedSummarizedExperiment", function(object) metadata(object)$A)

#' @rdname betaHMMHelpers
#' @export
setMethod("A", "betaHMMResults", function(object) metadata(object)$A)

#' @rdname betaHMMHelpers
#' @export
setMethod("A", "NULL", function(object) NA)


### tau
#' @rdname betaHMMHelpers
#' @export
setMethod("tau", "RangedSummarizedExperiment", function(object) metadata(object)$tau)

#' @rdname betaHMMHelpers
#' @export
setMethod("tau", "betaHMMResults", function(object) metadata(object)$tau)

#' @rdname betaHMMHelpers
#' @export
setMethod("tau", "NULL", function(object) NA)

### llk
#' @rdname betaHMMHelpers
#' @export
setMethod("llk", "RangedSummarizedExperiment", function(object) metadata(object)$llk)

#' @rdname betaHMMHelpers
#' @export
setMethod("llk", "betaHMMResults", function(object) metadata(object)$llk)

#' @rdname betaHMMHelpers
#' @export
setMethod("llk", "NULL", function(object) NA)

### phi
#' @rdname betaHMMHelpers
#' @export
setMethod("phi", "RangedSummarizedExperiment", function(object) metadata(object)$phi)

#' @rdname betaHMMHelpers
#' @export
setMethod("phi", "betaHMMResults", function(object) metadata(object)$phi)

#' @rdname betaHMMHelpers
#' @export
setMethod("phi", "NULL", function(object) NA)


### hidden_states
#' @rdname betaHMMHelpers
#' @export
setMethod("hidden_states", "RangedSummarizedExperiment", function(object) metadata(object)$hidden_states)

#' @rdname betaHMMHelpers
#' @export
setMethod("hidden_states", "betaHMMResults", function(object) metadata(object)$hidden_states)

#' @rdname betaHMMHelpers
#' @export
setMethod("hidden_states", "NULL", function(object) NA)

### DMR methods
#' @title DMR identification from DMCs
#' @description Function to identify the DMRs for the DMCs identified using BHMM
#' @details Function to identify the DMRs for the DMCs identified using BHMM
#' @rdname dmr_identification
#' @aliases
#' dmr_identification
#' dmr_identification-methods
#' dmr_identification,betaHMMResults-method
#' @export
#' @seealso \code{\link{betaHMM}}
#' @inheritParams dmr_identification_run
#' @param object a \code{betaHMMResults} object.
#'
#' @return A dataframe containing the following columns:
#' \itemize{
#' \item start_CpG - The starting CpG site in the particular DMR
#' \item end_CpG -  The ending CpG site in the particular DMR
#' \item DMR_size - Number of CPG sites identified in the DMR
#' \item chr_dmr - The chromosome corresponding to the CpG sites in the DMR.
#' \item map_start - MAPINFO of starting CpG site in the particular DMR
#' \item map_end -  MAPINFO of the ending CpG site in the particular DMR}
#' @importFrom stats complete.cases
#'
setMethod(f="dmr_identification", signature(object="betaHMMResults"),
          definition=function(object,diff_meth_cluster=1,...) {


            #object <- x
            dmr_df<-dmr_identification_run(object=object,diff_meth_cluster=diff_meth_cluster,...)
            return(dmr_df)
          })



