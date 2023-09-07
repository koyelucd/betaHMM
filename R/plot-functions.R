globalVariables(c("Patient_Sample", "label", "Uncertainty", "M", "beta_value"))
#' Visualize results from \code{betaHMM}/ \code{dmc_identification}/
#' \code{threshold_identification} functions.
#'
#' @rdname plot
#' @aliases
#' plot
#' plot-methods
#' plot,betaHMMResults-method
#'
#' @param x An object of class
#' \code{\link[betaHMM:betaHMMResults]{betaHMMResults}}/
#' \code{\link[betaHMM:dmcResults]{dmcResults}}/
#' \code{\link[betaHMM:dmcResults]{threshold_Results}} object.
#' @param chromosome The chromosome number for which the plot is to be
#' displayed.
#' @param what The different plots that can be obtained are either
#' 'fitted density','kernel density' or
#' 'uncertainty' (default = 'fitted density').
#' @param treatment_group The names of the different treatment groups
#' to be displayed in the plot.If no value is passed then the sample names
#' estimated by the \code{betaHMM} function are used.
#' @param AUC The AUC values for that chromosome.
#' @param uncertainty_threshold The uncertainty threshold value used for DMC
#' identification.
#' @param title The title that the user wants to display.
#' If no title is to be displayed the default is 'NULL'.
#' @param ... Other graphics parameters.
#'
#' @return This function displays the following plots as requested by the user
#' when analysing the \code{\link[betaHMM:betaHMMResults]{betaHMMResults}}
#' output:
#' \itemize{
#' \item fitted density estimates - Plot showing the fitted density estimates
#' of the clustering solution under the optimal model selected.
#' \item kernel density estimates - Plot showing the kernel density estimates
#' of the clustering solution under the optimal model selected.
#' \item uncertainty -  A boxplot showing the uncertainties in the
#' hidden state estimation.
#' }
#' @author Koyel Majumdar
#'
#' @export
#' @example inst/examples/betaHMM_package.R
#' @importFrom ggplot2 ggplot aes
#' @importFrom stats C density
#' @importFrom scales seq_gradient_pal
#'
setMethod(f = "plot", signature(x = "betaHMMResults"),
            definition = function(x, chromosome = NULL,
                                what = c("fitted density", "kernel density",
    "uncertainty"), treatment_group = NULL, AUC = NULL,
    uncertainty_threshold = 0.2, title = NULL, ...) {
    graph_objects <- c()
    object <- x
    graph_objects <- betaHMMGlobalplots(object,
                                        chromosome = chromosome,
                                        what = what,
                                        treatment_group = treatment_group,
                                        AUC = AUC,
                                        uncertainty_threshold =
                                            uncertainty_threshold,
                                        title = title, ...)
    return(graph_objects)
})

betaHMMGlobalplots<-function(x,chromosome=NULL,what="fitted density",
                            treatment_group=NULL,AUC=NULL,
                            uncertainty_threshold=0.2, title = NULL, ...) {
    if (is.null(chromosome)) stop("Chromosome number cannot be empty")
    chromosome <- as.numeric(chromosome)
    object <- x; plotdata <- annotatedData(object)
    plotdata <- as.data.frame(plotdata)
    label_chromosome <- paste("chr", chromosome, sep = " ")
    R <- R(object); N <- N(object); C <- nrow(data); K <- K(object)
    phi_complete <- phi(object)
    phi <- list()
    phi[["sp_1"]] <- phi_complete[[1]][[label_chromosome]]
    phi[["sp_2"]] <- phi_complete[[2]][[label_chromosome]]
    data_comp <- plotdata[plotdata$CHR == chromosome, ]
    data <- data_comp[, -c(seq(1, 3))];  C <- nrow(data)
    hidden_states <- hidden_states(object)[[label_chromosome]]
    if (is.null(title)) {
        title_text <- ""
    } else {
        title_text <- title
    }
    if (is.null(treatment_group)) {
        treatment_group <- treatment_group(object)
    }
    if (!is.null(AUC)) {
        auc <- as.data.frame(AUC[[label_chromosome]])
        auc <- auc[order(auc$AUC, decreasing = TRUE), ]
        auc$label <- paste("State ", auc$State, ", AUC = ",
                            round(auc$AUC, 2), sep = "")
    } else {
        label <- seq(1, K)
        auc <- as.data.frame(label)
    }
    if (what == "kernel density") {
        plot_graph <- kernel_density_plot(data, hidden_states, auc, C,
                                            treatment_group, title_text)
    }
    if (what == "fitted density") {
        plot_graph <- fitted_density_plot(phi, hidden_states, auc, R, K, C,
                                            treatment_group, title_text)
    }
    if (what == "uncertainty") {
        z <- assay(object)
        plot_graph <- uncertainty_plots(data_comp, z, hidden_states, C,
                                        uncertainty_threshold, title_text)
    }
    return(plot_graph)
}

