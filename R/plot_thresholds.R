#'
#' @rdname plot
#' @aliases
#' plot
#' plot-methods
#' plot,threshold_Results-method
#'
#' @export
#' @param x An object of class
#' \code{\link[betaHMM:betaHMMResults]{betaHMMResults}}/
#' \code{\link[betaHMM:dmcResults]{dmcResults}}/
#' \code{\link[betaHMM:dmcResults]{threshold_Results}} object.
#' @param plot_threshold The "TRUE" option displays the threshold points in the
#' graph for the 3 state betaHMM model (default = "FALSE").
#' @param title The title that the user wants to display.
#' If no title is to be displayed the default is 'NULL'.
#' @param ... Other graphics parameters.
#' @return The function displays the plot for the estimated shape parameters
#' and threshold for the methylation states in a single DNA treatment condition
#' from the \code{\link[betaHMM:threshold_Results]{threshold_Results}} object.
#'
#' @importFrom ggplot2 ggplot aes
#' @importFrom scales seq_gradient_pal
setMethod(f = "plot", signature(x = "threshold_Results"),
            definition = function(x, plot_threshold=TRUE,title = NULL, ...) {
    # x <- object
    graph_objects <- c()

    object <- x
    graph_objects <- thresholdGlobalplots(object,plot_threshold=plot_threshold,
                                            title = title, ...)
    return(graph_objects)
})

thresholdGlobalplots <- function(x, plot_threshold=TRUE,title = NULL, ...) {
    txt <- ifelse(is.null(title), "", title)
    data <- as.data.frame(annotatedData(x))
    df_plot <- subset(data, select = -c(IlmnID))
    data_x <- sort(df_plot[, 1]); K <- K(x)
    prop <- as.numeric(table(hidden_states(x)))/nrow(data)
    data_th_plot <- matrix(data = NA, nrow = 1, ncol = 3)
    data_th_plot <- as.data.frame(data_th_plot)
    colnames(data_th_plot) <- c("beta_value", "density", "Cluster")
    phi_complete <- model_parameters(x)$phi
    alpha <- phi_complete[[1]]; delta <- phi_complete[[2]]
    for (i in seq(1, K)) {
        beta_value <- data_x
        Cluster <- rep(i, length(data_x))
        density <- prop[i] * stats::dbeta(data_x, alpha[i], delta[i])
        temp <- cbind(beta_value, density, Cluster)
        data_th_plot <- rbind(data_th_plot, temp)  }
    data_th_plot <- as.data.frame(data_th_plot)
    data_th_plot <- data_th_plot[-1, ]
    data_th_plot$Cluster <- as.factor(data_th_plot$Cluster)
    plot_graph <- ggplot(data_th_plot) + geom_line(aes(beta_value,
                                                        density,
                                                        color = Cluster),
                                                    linetype = "solid") +
        labs(x="Beta value",y="Density",title = txt, color = "Hidden States")
    if (K == 3) {
        colours <- c("chartreuse3", "magenta", "cyan3")
        plot_graph <- plot_graph + scale_color_manual(values = colours)}
    p.data <- ggplot_build(plot_graph)$data[[1]]
    p.text <- lapply(split(p.data, f = p.data$group), function(df) {
        df[which.max(df$y), ]})
    p.text <- do.call(rbind, p.text); p.text$prop <- prop
    plot_graph <- plot_graph + annotate("text", x = p.text$x, y = 0.2,
                                        label = sprintf("%.3f", p.text$prop),
                                        colour = p.text$colour, vjust = 0)
    if(plot_threshold==TRUE)
    {
        if (!is.null(plot_graph)) {
            th_plot <- as.vector(unlist(threshold(x)))
            ano_th <- th_plot - 0.02 }
        num <- sort(p.text$y)
        ano_y <- num[length(num)] - 0.1
        plot_graph <- plot_graph + geom_vline(xintercept = th_plot,
                                              linetype = "dotted") +
            annotate("text", x = ano_th, y = ano_y,
                     label = th_plot, angle = 90)
    }
    return(plot_graph)}
