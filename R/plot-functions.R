globalVariables(c("Patient_Sample","label","Uncertainty","M","beta_value"))
#' Visualize results from betaHMM clustering
#'
#' Plot a betaHMMResults object.
#'
#' @rdname plot
#' @aliases
#' plot
#' plot-methods
#' plot,betaHMMResults-method
#'
#' @param x An object of class \code{"betaHMMResults"}
#' @param what The different plots that can be obtained are either
#'             "fitted density","kernel density" or
#'             "uncertainty" (default = "fitted density").
#' @param treatment_group The names of the different treatment groups.
#'  If no value is passed then default values of sample names, e.g. Sample 1,
#'  Sample 2, etc are used as legend text (default = NULL).
#' @param title The title that the user wants to display.
#' If no title is to be displayed the default is "NULL".
#' @param ... Other graphics parameters.
#'
#' @return This function displays the following plots as requested by the user:
#' \itemize{
#' \item fitted density estimates - Plot showing the fitted density estimates of the clustering solution under the optimal model selected.
#' \item kernel density estimates - Plot showing the kernel density estimates of the clustering solution under the optimal model selected.
#' \item uncertainty -  A boxplot showing the uncertainties in the optimal clustering solution.
#' }
#'
#' @author Koyel Majumdar
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes
#' @importFrom stats C
#' @importFrom scales seq_gradient_pal
#'
setMethod(f="plot", signature(x="betaHMMResults"),
          definition=function(x,  what = c("fitted density","kernel density",
                                           "uncertainty"),
                              treatment_group = NULL,
                              title = NULL,...) {
            # x <- object
            graph_objects <- c()

            object <- x
            graph_objects<-betaHMMGlobalplots(object,what=what,
                                              treatment_group=treatment_group,
                                              title=title,...)
            return(graph_objects)
          })

betaHMMGlobalplots<-function(x,
                       what = "fitted density",
                       treatment_group = NULL,
                       title = NULL,...)
{
  #UseMethod("plot")
  object <- x
  plotdata=annotatedData(object)
  data=plotdata[,-c(1:3)]
  R=R(object)
  N=N(object)
  C=nrow(data)
  K=K(object)
  phi=phi(object)
  hidden_states=hidden_states(object)
  if(is.null(title))
  {
    txt=""
  }else{
    txt=title
  }
  if(is.null(treatment_group))
  {

    treatment_group=sapply(1:R, function(x) paste0("Sample ",x))
  }

  if(what == "kernel density")
  {
    if(is.null(data)){
      plot_graph=NULL
      warning("data argument cannot be NULL for generation of kernel density plot.", call. = FALSE)
    }else
    {
      column_len=ncol(data)
      if(column_len==(N*R))
      {
        call_data=data
      }else if(column_len>(N*R)){
        call_data=data[,1:(N*R)]
      }else{
        call_data=NULL
      }
      if(is.null(call_data))
      {
        plot_graph=NULL
        warning("Differential methylation analysis to be done on multiple
                treatment groups.", call. = FALSE)
      }
      else{
        data_ggplot<-as.data.frame(call_data)
        #data_ggplot$mem_final<-as.factor(data_ggplot$mem_final)
        data_ggplot$mem_final<-as.factor(hidden_states)
        colnames(data_ggplot)[length(data_ggplot)]<-"Cluster"
        cols=ncol(data_ggplot)
        rows=nrow(data_ggplot)
        data_matrix<-as.matrix(data_ggplot[,1:(cols-1)])
        data_new<-as.vector(data_matrix)
        col_names<-colnames(data_ggplot)
        col_len<-length(col_names)
        Cluster<-vector()
        Patient_sample=vector()
        for(i in 1:(col_len-1))
        {
          temp=gsub("_"," ",col_names[i])
          ps_names<-rep(temp,rows)
          Patient_sample<-c(Patient_sample,ps_names)
          Cluster<-c(Cluster,data_ggplot[,cols])
        }
        data_plot<-cbind(data_new,Cluster,Patient_sample)
        data_plot<-as.data.frame(data_plot)
        colnames(data_plot)<-c("beta_value","Cluster","Patient_Sample")
        data_plot$beta_value<-as.numeric(data_plot$beta_value)
        color_length<-col_len-1
        colours<-scales::seq_gradient_pal(low="#FFC20A",
                                          high="#0C7BDC",
                                          space = "Lab"
        )(1:color_length/color_length)

        plot_graph<-ggplot2::ggplot(data_plot)+
          ggplot2::geom_density(aes(x=beta_value,color=Patient_Sample))+
          ggplot2::xlab("Beta Value")+
          ggplot2::ylab("Density")+
          ggplot2::scale_color_manual("Treatment group",values=colours)+
          ggplot2::facet_wrap(~Cluster,scales = "free_y"
          )+
          ggplot2::theme(axis.title.x = ggplot2::element_text(size=10),
                         axis.title.y = ggplot2::element_text(size=10)) +
          ggplot2::ggtitle(txt)
        #ggplot2::ggtitle("Density estimates for K.R clustering solution")

        cluster_size=table(hidden_states)
        f_labels<-data.frame(Cluster=seq(1,length(cluster_size),by=1),
                             label=as.vector(round((cluster_size/C),3)))
        plot_graph<-plot_graph+
          ggplot2::geom_text(x = 0.2, y = 1, ggplot2::aes(label = label),
                             data = f_labels)

      }
    }

    # }


  }
  if(what=="fitted density")
  {
    vec_C=1001
    alpha<-t(phi$sp_1)
    delta<-t(phi$sp_2)
    tau<-round(as.vector(table(hidden_states)/
                           length(hidden_states)),3)

    density_vec<-vector()
    cluster_vec<-vector()
    sample_vec<-vector()
    beta_vec<-vector()
    vec_x<-seq(0.001, 0.999, length=vec_C)
    #treatment_group<-c("Sample A","Sample B")

    for(i in 1:R)
    {
      for(j in 1:K)
      {
        #j=1
        tmp_vec<-sapply(vec_x,function(x) {tau[j]*
            (stats::dbeta(x,alpha[j,i]
                          ,delta[j,i]))})
        beta_vec<-c(beta_vec,vec_x)
        density_vec<-c(density_vec,tmp_vec)
        cluster_vec<-c(cluster_vec,rep(j,times=length(tmp_vec)))
        sample_vec<-c(sample_vec,rep(treatment_group[i],times=length(tmp_vec)))

      }
    }

    df_new_tmp<-as.data.frame(cbind(beta_vec,density_vec,cluster_vec,
                                    sample_vec))
    df_new_tmp$sample_vec<-as.factor(df_new_tmp$sample_vec)
    df_new_tmp$cluster_vec<-as.factor(df_new_tmp$cluster_vec)
    df_new_tmp$beta_vec<-as.numeric(df_new_tmp$beta_vec)
    df_new_tmp$density_vec<-as.numeric(df_new_tmp$density_vec)
    color_length<-R
    colours<-scales::seq_gradient_pal(low="#FFC20A",
                                      high="#0C7BDC",space =
                                        "Lab")(1:color_length/color_length)
    cluster_size=table(hidden_states)
    plot_graph<-ggplot2::ggplot(df_new_tmp,ggplot2::aes(x=beta_vec,
                                                        y=density_vec,
                                                        color=sample_vec))+
      ggplot2::geom_line()+
      ggplot2::scale_color_manual(values=colours)+
      ggplot2::facet_wrap(~cluster_vec,scales = "free_y"
      )+ ggplot2::labs(color="DNA sample", x="Beta value", y="Density")+
      ggplot2::ggtitle(txt)
    f_labels<-data.frame(cluster_vec=as.factor(seq(1,K,by=1)),
                         label=as.vector(round((cluster_size/C),3)))
    plot_graph<-plot_graph+
      ggplot2::geom_text(data = f_labels, ggplot2::aes(x = 0.2, y = 0.1,
                                                       label = label,
                                                       color=NA),
                         show.legend = F,fontface="bold" )

  }


  if(what == "uncertainty")
  {
    labels<-c(max_uc="Maximum uncertainty")
    classification_final=apply(assay(object),1,max)
    uncertainty=1-classification_final
    tau=(table(hidden_states))/C
    unc_df<-cbind(uncertainty,hidden_states)
    unc_df<-as.data.frame(unc_df)
    colnames(unc_df)<-c("Uncertainty","Cluster")
    unc_df$Cluster<-as.factor(unc_df$Cluster)
    unc_df_sorted<-unc_df[order(unc_df$Cluster),]
    max_unc=1-1/(K)
    h=max_unc+0.015
    max_uncertainty <- data.frame(yintercept=h, max_uncertainty=factor(h))
    plot_graph<-ggplot2::ggplot(unc_df_sorted, ggplot2::aes(x=Cluster,
                                                            y=Uncertainty,
                                                            color=Cluster)) +
      ggplot2::geom_boxplot()+
      ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                     axis.ticks.x=ggplot2::element_blank()
      )+ggplot2::labs(color="Cluster number")+
      ggplot2::xlab("Cluster number")+
      #ggplot2::theme(legend.position = "none")+
      ggplot2::ggtitle(txt)+
      #ggplot2::ggtitle("Boxplot for uncertainties in clustering solution")+
      ggplot2::coord_cartesian(ylim = c(0, 1))+
      ggplot2::geom_hline(ggplot2::aes(yintercept=max_unc,
                                       linetype="Maximum uncertainty"),
                          color="black")+
      ggplot2::scale_linetype_manual(name="",
                                     values = 2,guide =
                                       ggplot2::guide_legend(override.aes =
                                                               list(color =
                                                                      "black")))



  }

 return(plot_graph)


}


