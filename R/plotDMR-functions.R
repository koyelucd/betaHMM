globalVariables(c("Patient_Sample", "label", "M",
                    "beta_value", "name", "CP", "DMC", "group", "value"))
#' Visualize the DMCs and DMRs identified
#'
#' Plot a \code{\link[betaHMM:dmcResults]{dmcResults}} object.
#'
#' @rdname plot
#' @aliases
#' plot
#' plot-methods
#' plot,dmcResults-method
#'
#' @param x An object of class \code{\link[betaHMM:dmcResults]{dmcResults}}.
#' @param start_CpG The IlmnID of starting CpG site.
#' @param end_CpG The IlmnID of ending CpG site/ the total number
#' of CpGs to be plotted excluing the starting CpG site.
#' @param treatment_group The names of the different treatment groups to be
#' displayed in the plot.
#' If no value is passed then the treatment group names from the
#' \code{betaHMM} function are used.
#' @param N Number of replicates in each treatment group.
#' @param title The title that the user wants to display.
#' If no title is to be displayed the default is 'NULL'.
#' @param ... Other graphics parameters.
#'
#' @return This function displays the DMCs and DMRs plot from
#' the \code{\link[betaHMM:dmcResults]{dmcResults}} object.
#'
#' @author Koyel Majumdar
#'
#' @export
#'
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @import cowplot
#' @importFrom stats C
#' @importFrom scales seq_gradient_pal
#' @importFrom stringr str_split
#' @example inst/examples/betaHMM_package.R
setMethod(f = "plot", signature(x = "dmcResults"),
            definition = function(x, start_CpG = NULL, end_CpG = NULL,
                                treatment_group = NULL, N = NULL,
                                title = NULL, ...) {
    graph_objects <- c()
    object <- x
    graph_objects <- dmcGlobalplots(object, start_CpG = start_CpG,
                                    end_CpG = end_CpG,
                                    treatment_group = treatment_group,
                                    N = N, title = title, ...)
    return(graph_objects)
})

dmcGlobalplots<-function(x,start_CpG=NULL,end_CpG=NULL,treatment_group=NULL,
                            N = NULL, title = NULL, ...) {
    if (is.null(start_CpG)) stop("start_CpG cannot be empty.")
    if (is.null(end_CpG)) stop("end_CpG cannot be empty.")
    if (is.numeric(end_CpG)) {start_index <- which(dmc_df$IlmnID == start_CpG)
        new_index<-start_index+end_CpG; end_CpG<-dmc_df[new_index, "IlmnID"]}
    object <- x; dmc_df <- as.data.frame(assay(object))
    if (dmc_df[dmc_df$IlmnID == start_CpG, "CHR"] !=
        dmc_df[dmc_df$IlmnID == end_CpG, "CHR"])
        stop("Start and End CpGs cannot be from different Chromosomes.")
    chromosome <- dmc_df[dmc_df$IlmnID == start_CpG, "CHR"]
    ddcc<-dmc_df[dmc_df$CHR==chromosome,]
    label_chrom <- paste("chr", chromosome, sep = " ")
    ddcc$Uncertainty <- unlist(uncertainty(object)[[label_chrom]])
    si <- which(ddcc$IlmnID == start_CpG); ei <- which(ddcc$IlmnID == end_CpG)
    dmc_df_chr <- ddcc[si:ei, ]
    col_meth<-which(!(colnames(dmc_df_chr) %in% c("CHR","MAPINFO","IlmnID",
                                                    "hidden_state", "DMC",
                                                    "Uncertainty")))
    if (is.null(R)) { R <- R(object) }
    if (is.null(N)) { N <- N(object) }
    if (is.null(treatment_group)) { treatment_group<-treatment_group(object) }
    if (!is.null(N) & !is.null(treatment_group)) { group_name <- NA
        for (j in seq(1, length(treatment_group))) {
            group_name <- c(group_name, rep(treatment_group[j], N[j]))}
        group_name <- group_name[-1]
        group_name <- paste(group_name, seq(1, length(group_name)), sep = "_")
        colnames(dmc_df_chr)[col_meth] <- group_name}
    cols <- colnames(dmc_df_chr)
    long_dat <- pivot_longer(dmc_df_chr, cols = all_of(col_meth)) %>%
        mutate(group=str_split(name,"_",simplify=TRUE)[, 1])
    g1<-ggplot(long_dat)+geom_point(aes(x=factor(MAPINFO),y=1,fill=Uncertainty,
                                        colour=factor(DMC)),shape=21,size=5,
                                    stroke=1.2)+theme_void()+
        scale_shape_manual(values=c(`1`=19,`0`=1))+
        scale_colour_manual(values=c(`1`="red",`0`="black"),
                            labels=c("non-DMC","DMC"),name=NULL)+
        scale_fill_gradient(low = "black", high = "white",name="Uncertainty")+
        guides(shape = "none") + labs(x = NULL, y = NULL) +
        theme(legend.position = "top", legend.justification = "right")
    g2b<-ggplot(long_dat,aes(x=factor(MAPINFO),y=value,color=group,
                                group = name)) + geom_line() + theme_void() +
        ylab("beta value") +
        scale_x_discrete(labels = paste0(unique(long_dat$IlmnID), " ")) +
        theme(strip.background = element_blank(),
        strip.text.x = element_blank(), axis.text.y = element_blank(),
        axis.title.y=element_text(angle=90),axis.text.x=element_text(angle=90,
    hjust=.95,vjust=.2),legend.title=element_blank(),legend.position="bottom")
    plotG<-plot_grid(g1,g2b,ncol=1,rel_heights=c(0.3,0.7),align="v",axis="l")
    return(plotG)}


