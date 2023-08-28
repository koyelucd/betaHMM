globalVariables(c("Patient_Sample","label","Uncertainty","M","beta_value"))
#' Visualize results from betaHMM clustering
#'
#' Plot a \code{betaHMMResults} object.
#'
#' @rdname plot
#' @aliases
#' plot
#' plot-methods
#' plot,betaHMMResults-method
#'
#' @param x An object of class
#' \code{\link[betaHMM:betaHMMResults]{betaHMMResults}} object.
#' @param chromosome The chromosome number for which the plot is to be displayed
#' @param what The different plots that can be obtained are either
#'             "fitted density","kernel density" or
#'             "uncertainty" (default = "fitted density").
#' @param treatment_group The names of the different treatment groups
#' to be displayed in the plot.
#'  If no value is passed then the sample names estimated by the \code{betaHMM}
#'  function are used.
#' @param AUC The AUC values for that chromosome.
#' @param uncertainty_threshold The uncertainty threshold value used for DMC
#'                               identification.
#' @param title The title that the user wants to display.
#' If no title is to be displayed the default is "NULL".
#' @param ... Other graphics parameters.
#'
#' @return This function displays the following plots as requested by the user:
#' \itemize{
#' \item fitted density estimates - Plot showing the fitted density estimates
#' of the clustering solution under the optimal model selected.
#' \item kernel density estimates - Plot showing the kernel density estimates
#' of the clustering solution under the optimal model selected.
#' \item uncertainty -  A boxplot showing the uncertainties in the
#' hidden state estimation.
#' }
#'
#' @author Koyel Majumdar
#'
#' @export
#' @example inst/examples/betaHMM_package.R
#' @importFrom ggplot2 ggplot aes
#' @importFrom stats C
#' @importFrom scales seq_gradient_pal
#'
setMethod(f = "plot", signature(x = "betaHMMResults"),
          definition = function(x, chromosome = NULL,
                              what = c("fitted density", "kernel density",
                                           "uncertainty"),
                              treatment_group = NULL,
                              AUC = NULL,
                              uncertainty_threshold = 0.2,
                              title = NULL, ...) {
            # x <- object
            graph_objects <- c()

            object <- x
            graph_objects <- betaHMMGlobalplots(object, chromosome = chromosome,
                                                what = what,
                                                treatment_group =
                                                  treatment_group,
                                                AUC = AUC,
                                                uncertainty_threshold =
                                                uncertainty_threshold,
                                                title = title, ...)
            return(graph_objects)
          })

betaHMMGlobalplots <- function(x, chromosome = NULL,
                       what = "fitted density",
                       treatment_group = NULL,
                       AUC = NULL,
                       uncertainty_threshold = 0.2,
                       title = NULL, ...)
{
  #UseMethod("plot")

  if(is.null(chromosome))
    stop("Chromosome number cannot be empty")

  chromosome <- as.numeric(chromosome)
  object <- x
  plotdata <- annotatedData(object)
  plotdata <- as.data.frame(plotdata)
  label_chromosome <- paste("chr", chromosome, sep = " ")
  #data=plotdata[,-c(1:3)]
  R <- R(object)
  N <- N(object)
  C <- nrow(data)
  K <- K(object)
  phi_complete <- phi(object)
  # phi=phi(object)
  # hidden_states=hidden_states(object)

    phi <- list()
    phi[["sp_1"]] <- phi_complete[[1]][[label_chromosome]]
    phi[["sp_2"]] <- phi_complete[[2]][[label_chromosome]]
    data_comp <- plotdata[plotdata$CHR == chromosome,]
    data <- data_comp[ , -c(1:3)]
    C <- nrow(data)
    hidden_states <- hidden_states(object)[[label_chromosome]]
  if(is.null(title))
  {
    title_text <- ""
  }else{
    title_text <- title
  }
  if(is.null(treatment_group))
  {

    treatment_group <- treatment_group(object)
  }
    if(!is.null(AUC))
    {
      auc <- as.data.frame(AUC[[label_chromosome]])
      auc <- auc[order(auc$AUC,decreasing = T), ]
      auc$label <- paste("State ",auc$State, ", AUC = ", round(auc$AUC, 2), sep = "")
      #auc_label=as.factor(auc$label)
    }else
    {label <- seq(1,K)
    auc <- as.data.frame(label)
    }
   # par(ask=TRUE)
  if(what == "kernel density")
  {
    if(is.null(data)){
      plot_graph <- NULL
      warning("data argument cannot be NULL for generation of kernel
              density plot.", call. = FALSE)
    }else
    {
      column_len <- ncol(data)
      if(column_len == sum(N))
      {
        call_data <- data
      }else if(column_len > sum(N)){
        call_data <- data[ ,1:sum(N)]
      }else{
        call_data <- NULL
      }
      if(is.null(call_data))
      {
        plot_graph <- NULL
        warning("Differential methylation analysis to be done on multiple
                treatment groups.", call. = FALSE)
      }
      else{
        data_ggplot <- as.data.frame(call_data)
        data_ggplot$mem_final <- as.factor(hidden_states)
        colnames(data_ggplot)[length(data_ggplot)] <- "Cluster"
        cols <- ncol(data_ggplot)
        rows <- nrow(data_ggplot)
        data_matrix <- as.matrix(data_ggplot[ ,1:(cols-1)])
        data_new <- as.vector(data_matrix)
        col_names <- colnames(data_ggplot)
        col_len <- length(col_names)
        Cluster <- vector()
        Patient_sample <- vector()
        for(i in 1:(col_len-1))
        {
          temp <- gsub("_", " ", col_names[i])
          ps_names <- rep(temp, rows)
          Patient_sample <- c(Patient_sample, ps_names)
          Cluster <- c(Cluster, data_ggplot[ ,cols])
        }
        data_plot <- data.frame(data_new = data_new, Cluster = Cluster,
                                Patient_sample = Patient_sample)
        colnames(data_plot) <- c("beta_value", "Cluster", "Patient_Sample")
        #data_plot$beta_value <- as.numeric(data_plot$beta_value)
        data_plot$Cluster<-factor(data_plot$Cluster,levels=auc$State)
        data_plot$Cluster_full <- as.factor(data_plot$Cluster)
        levels(data_plot$Cluster_full) <- auc$label
        color_length <- col_len-1
        colours <- scales::seq_gradient_pal(low = "#FFC20A",
                                            high = "#0C7BDC",
                                            space = "Lab"
        )(1:color_length/color_length)

        plot_graph <- ggplot2::ggplot(data_plot) +
          ggplot2::geom_density(ggplot2::aes(x = beta_value,
                                             color = Patient_Sample))+
          ggplot2::xlab("Beta Value")+
          ggplot2::ylab("Density")+
          ggplot2::scale_color_manual("DNA Samples",values=colours)+
          ggplot2::facet_wrap(~factor(Cluster_full,levels=auc$label),
                              scales = "free_y",
          )+
          ggplot2::theme(axis.title.x = ggplot2::element_text(size = 10),
                         axis.title.y = ggplot2::element_text(size = 10)) +
          ggplot2::ggtitle(title_text)

        cluster_size <- table(hidden_states)
        f_labels <- data.frame(Cluster_full = levels(data_plot$Cluster_full),
                             label =
                               as.vector(round((cluster_size[auc$State]/C),3)))
        plot_graph <- plot_graph +
          ggplot2::geom_text(x = 0.2, y = 1, ggplot2::aes(label = label),
                             data = f_labels)
      }
    }

    # }


  }
  if(what == "fitted density")
  {
    vec_C <- 1001
    alpha <- t(phi$sp_1)
    delta <- t(phi$sp_2)
    tau <- round(as.vector(table(hidden_states) /
                           length(hidden_states)), 3)

    density_vec <- vector()
    cluster_vec <- vector()
    sample_vec <- vector()
    beta_vec <- vector()
    vec_x <- seq(0.001, 0.999, length = vec_C)
    #treatment_group<-c("Sample A","Sample B")

    for(i in 1:R)
    {
      for(j in 1:K)
      {
        #j=1
        tmp_vec <- sapply(vec_x, function(x) {tau[j]*
            (stats::dbeta(x, alpha[j,i], delta[j,i]))})
        beta_vec <- c(beta_vec, vec_x)
        density_vec <- c(density_vec, tmp_vec)
        cluster_vec <- c(cluster_vec, rep(j, times = length(tmp_vec)))
        sample_vec <- c(sample_vec, rep(treatment_group[i],
                                        times = length(tmp_vec)))

      }
    }

    df_new_tmp <- data.frame(beta_vec = beta_vec, density_vec = density_vec,
                             cluster_vec = cluster_vec, sample_vec = sample_vec)
    # df_new_tmp$sample_vec <- as.factor(df_new_tmp$sample_vec)
    # df_new_tmp$cluster_vec <- as.factor(df_new_tmp$cluster_vec)
    # df_new_tmp$beta_vec <- as.numeric(df_new_tmp$beta_vec)
    # df_new_tmp$density_vec <- as.numeric(df_new_tmp$density_vec)
    df_new_tmp$cluster_vec<-factor(df_new_tmp$cluster_vec,levels=auc$State)
    df_new_tmp$Cluster_full <- as.factor(df_new_tmp$cluster_vec)
    levels(df_new_tmp$Cluster_full) <- auc$label
    color_length <- R
    colours <- scales::seq_gradient_pal(low = "#FFC20A",
                                        high = "#0C7BDC",
                                        space =
                                        "Lab")(1:color_length/color_length)
    cluster_size <- table(hidden_states)

    plot_graph <- ggplot2::ggplot(df_new_tmp, ggplot2::aes(x = beta_vec,
                                                           y = density_vec,
                                                           color = sample_vec))+
      ggplot2::geom_line()+
      ggplot2::scale_color_manual(values = colours)+
      ggplot2::facet_wrap(~factor(Cluster_full,levels=auc$label),
                          scales = "free_y"
      )+ ggplot2::labs(color = "Treatment Groups", x = "Beta value",
                       y = "Density")+
      ggplot2::ggtitle(title_text)
    f_labels <- data.frame(Cluster_full = levels(df_new_tmp$Cluster_full),
                           label =
                             as.vector(round((cluster_size[auc$State] / C), 3)))
    plot_graph <- plot_graph +
      ggplot2::geom_text(data = f_labels, ggplot2::aes(x = 0.2, y = 0.1,
                                                       label = label,
                                                       color=NA),
                         show.legend = F, fontface="bold" )

  }


  if(what == "uncertainty")
  {
    labels <- c(max_uc = "Maximum uncertainty")
    z <- assay(object)
    z_chr <- z[rownames(z) %in% data_comp$IlmnID,  ]
    classification_final <- apply(z_chr, 1, max)
    uncertainty <- 1 - classification_final
    tau <- (table(hidden_states)) / C
    unc_df <- cbind(uncertainty, hidden_states)
    unc_df <- as.data.frame(unc_df)
    colnames(unc_df) <- c("Uncertainty", "Cluster")
    unc_df$Cluster <- as.factor(unc_df$Cluster)
    unc_df_sorted <- unc_df[order(unc_df$Cluster),  ]
    max_unc <- uncertainty_threshold
    h <- max_unc + 0.015
    max_uncertainty <- data.frame(yintercept = h, max_uncertainty = factor(h))
    plot_graph <- ggplot2::ggplot(unc_df_sorted, ggplot2::aes(x = Cluster,
                                                              y = Uncertainty))+
                                                           # color=Cluster)) +
      ggplot2::geom_boxplot()+
      # ggplot2::theme(axis.text.x=ggplot2::element_blank(),
      #                axis.ticks.x=ggplot2::element_blank()
      # )+ggplot2::labs(color="Hidden States")+
      ggplot2::xlab("Hidden States")+
      #ggplot2::theme(legend.position = "none")+
      ggplot2::ggtitle(title_text)+
      #ggplot2::ggtitle("Boxplot for uncertainties in clustering solution")+
      ggplot2::coord_cartesian(ylim = c(0, 1))+
      ggplot2::geom_hline(ggplot2::aes(yintercept = max_unc,
                                       linetype =
                                         paste("Uncertainty threshold = ",
                                                      max_unc)),
                          color = "black")+
      ggplot2::scale_linetype_manual(name = "",
                                     values = 2,guide =
                                       ggplot2::guide_legend(override.aes =
                                                               list(color =
                                                                      "black")))



  }

 return(plot_graph)


}


