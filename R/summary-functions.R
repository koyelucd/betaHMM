#' Summarize results from the three betaHMM workflow functions
#'
#' A function to summarize the \code{betaHMMResults}, \code{dmcResults} or
#' \code{dmrResults} objects.
#'
#' @rdname summary
#' @aliases
#' summary
#' summary-methods
#' summary,betaHMMResults-method
#' summary,dmcResults-method
#' summary,dmrResults-method
#'
#' @param object An object of class \code{'betaHMMResults'} or
#' \code{'dmcResults'} or \code{'dmrResults'}.
#' @param ... Additional arguments
#' @author Koyel majumdar
#' @seealso \code{\link{betaHMM}}
#' @example inst/examples/betaHMM_package.R
#' @return Summary of the \code{'betaHMMResults'} or
#' \code{'dmcResults'} or \code{'dmrResults'} object.
#' @keywords methods
#' @export
#' @importFrom utils tail
setMethod(f = "summary", signature(object = "betaHMMResults"),
            function(object, ...) {
    chr_unique <- chromosome_number(object)
    anno <- annotatedData(object)
    llk_ch <- vector()
    prop_hs <- vector()
    chr_text <- vector()
    for (i in seq(1, length(chr_unique))) {
        label_chromosome <- paste("chr", chr_unique[i], sep = " ")
        llk_ch[i] <- tail(llk(object)[[label_chromosome]], 1)
        chr_text[i] <- label_chromosome
        prop_hs[i] <- paste(round((table(
            hidden_states(object)[[label_chromosome]])/
            length(hidden_states(object)[[label_chromosome]])), 3),
            collapse = ",")
    }
    patients <- paste(N(object), collapse = ",")
    cat("*************************************************\n")
    cat("betaHMM workflow function: Parameter estimation")
    cat("\n")
    cat(paste0("Number of hidden states estimated:", K(object)))
    cat("\n")
    cat(paste0("No. of total CpG sites:", nrow(anno)))
    cat("\n")
    cat(paste0("Number of treatment conditions :", R(object)))
    cat("\n")
    cat(paste0("No. of patients in each treatment group:", patients))
    cat("\n")
    cat("*************************************************\n")
    cat("\n")
    cat("Summary of each chromosome analysed\n")
    cat("\n")
    tab2 <- data.frame(`Chromosome Number` = chr_unique, `log-likelihood` =
                        llk_ch, `Prop. of CpG sites in each hidden state` =
                        prop_hs, row.names = rep("",
        length(chr_unique)), check.names = FALSE)
    tab2 <- as.matrix(tab2)
    print(tab2)
})

##############################################################
#' @export
#' @rdname summary
setMethod(f = "summary", signature(object = "dmcResults"),
            function(object, ...) {
    dmc_df <- assay(object)
    chr_unique <- unique(dmc_df$CHR)
    no_of_dmc <- vector()
    cpg_site <- vector()
    for (i in seq(1, length(chr_unique))) {
        df <- dmc_df[dmc_df$CHR == chr_unique[i], ]
        no_of_dmc[i] <- length(which(df$DMC == 1))
        cpg_site[i] <- nrow(df)
    }
    patients <- paste(N(object), collapse = ",")
    cat("*************************************************\n")
    cat("betaHMM workflow function: DMC identification")
    cat("\n")
    cat(paste0("No. of total CpG sites:", nrow(dmc_df)))
    cat("\n")
    cat(paste0("Number of treatment conditions :", R(object)))
    cat("\n")
    cat(paste0("No. of patients in each treatment group:", patients))
    cat("\n")
    cat(paste0("Total number of DMCs:", table(dmc_df$DMC)[2]))
    cat("\n")
    cat("*************************************************\n")
    cat("\n")
    cat("\nNo. of CpG sites and DMCs in each chromosome:\n")
    tab2 <- data.frame(`Chromosome Number` = chr_unique, `No. of CpG sites` =
                        cpg_site, `No. of DMCs` = no_of_dmc,
                        row.names = rep("", length(chr_unique)),
                        check.names = FALSE)
    tab2 <- as.matrix(tab2)
    print(tab2)
})
##############################################################
#' @export
#' @rdname summary
setMethod(f = "summary", signature(object = "dmrResults"),
            function(object, ...) {
    dmr_df <- assay(object)
    chr_unique <- chromosome_number(object)
    no_of_dmr <- vector()
    avg_length <- vector()
    for (i in seq(1, length(chr_unique))) {
        df <- dmr_df[dmr_df$CHR == chr_unique[i], ]
        no_of_dmr[i] <- nrow(df)
        avg_length[i] <- round(mean(df$DMR_size), 3)
    }
    cat("*************************************************\n")
    cat("betaHMM workflow function: DMR identification")
    cat("\n")
    cat(paste0("Total number of DMRs:", nrow(dmr_df)))
    cat("\n")
    cat("*************************************************\n")
    cat("\n")
    cat("\n")
    cat("\nSummary of DMRs identified in each chromosome:\n")
    tab2 <- data.frame(`Chromosome Number` = chr_unique, `No. of DMRs` =
                        no_of_dmr, `Average length of DMRs` = avg_length,
                        row.names = rep("",length(chr_unique)),
                        check.names = FALSE)
    tab2 <- as.matrix(tab2)
    print(tab2)
})

