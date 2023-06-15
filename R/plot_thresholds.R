#' @title The plot for threshold indentification
#' @description The plot for threshold indentification
#' @details The plot for threshold indentification
#' @export
#' @param x A \code{\link[betaHMM:threshold_identification]{threshold_identification}} object.
#' @param data A dataframe of dimension \eqn{C \times N} containing methylation
#'             values for \eqn{C} CpG sites from one treatment group where
#'             DNA samples are either collected from \eqn{N} patients
#'             or has \eqn{N} replicates for the treatment group and this
#'             dataset was passed as an argument to the
#'             \code{\link[betaHMM:threshold_identification]{threshold_identification}} function.
#' @param title The title for the plot.
#' @return The plot for the estimated shape parameters and threshold for the methylation states.
#' @importFrom ggplot2 ggplot aes
#' @importFrom  plotly ggplotly
#' @importFrom scales seq_gradient_pal

plot_thresholds<-function(x,data,title=NULL)
{
  if(is.null(title))
  {
    txt=""
  }else{
    txt=title
  }
  data_x=sort(data[,1])
  K=3
  prop=as.numeric(table(x$states))/nrow(data)
  data_th_plot<-matrix(data=NA,nrow=1,ncol=3)
  data_th_plot<-as.data.frame(data_th_plot)
  colnames(data_th_plot)<-c("beta_value","density","Cluster")
  alpha=x$model_params$phi$sp_1
  delta=x$model_params$phi$sp_2

    for(i in 1:K)
    {
      beta_value=data_x
      density=prop[i]*stats::dbeta(data_x,alpha[i],delta[i])
      Cluster<-rep(i,length(data_x))
      temp<-cbind(beta_value,density,Cluster)
      data_th_plot<-rbind(data_th_plot,temp)
    }

  data_th_plot<-as.data.frame(data_th_plot)
  data_th_plot<-data_th_plot[-1,]
  data_th_plot$Cluster<-as.factor(data_th_plot$Cluster)
  #txt=""
  #colours<-c("chartreuse3","magenta","cyan3")
  plot_graph<-ggplot2::ggplot(data_th_plot)+
    ggplot2::geom_line(ggplot2::aes(beta_value,density,color=Cluster),
                       linetype = "solid")+
    ggplot2::labs(x="Beta value",y="Density",title=txt,
                  color ="Cluster")
  # +
  #   ggplot2::scale_color_manual(values=colours)
  if(K==3)
  {
    colours<-c("chartreuse3","magenta","cyan3")
    plot_graph<-plot_graph+
      ggplot2::scale_color_manual(values=colours)
  }
  p.data <- ggplot2::ggplot_build(plot_graph)$data[[1]]


  p.text <- lapply(split(p.data, f = p.data$group), function(df){
    df[which.max(df$y), ]
  })
  p.text <- do.call(rbind, p.text)
  p.text$prop=prop
  plot_graph<-plot_graph + ggplot2::annotate('text', x = p.text$x,
                                             #y = p.text$y,
                                             y=0.2,
                                             label = sprintf('%.3f',
                                                             p.text$prop),
                                             colour=p.text$colour,
                                             vjust = 0)


    if(!is.null(plot_graph)){

        th_plot<-x$threshold$threholds
        ano_th<-x$threshold$threholds-.02
           }
      num=sort(p.text$y)
      ano_y=num[length(num)]-0.1
      plot_graph<-plot_graph+ggplot2::geom_vline(xintercept = th_plot,linetype="dotted")+ggplot2::annotate("text",x=ano_th,y=ano_y,label=th_plot,angle=90)

      plot_graph
}
