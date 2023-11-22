globalVariables(c("x", "Cluster"))
#' @keywords internal
#' @importFrom ggplot2 ggplot geom_text aes geom_density ggtitle
#' @importFrom ggplot2 scale_colour_manual scale_fill_gradient facet_wrap theme
#' @importFrom ggplot2 geom_boxplot coord_cartesian geom_hline labs xlab ylab
#' @importFrom ggplot2 scale_linetype_manual guide_legend element_text
#' @importFrom stats kmeans sd dbeta density
#' @importFrom scales seq_gradient_pal

initialize_parameters_th <- function(data, M, N, R, seed = NULL) {
    K <- M
    A <- 0.6 * diag(K) + rep(0.05/K, K)  ## transmission probabilities
    tau <- rep(1/K, K)  ## initial mixing proportions
    data <- as.data.frame(data)
    for (i in seq(1, ncol(data))) {
        data[, i] <- as.numeric(data[, i])
    }
    data <- as.matrix(data)
    k_cluster <- kmeans(data, M)
    mem <- k_cluster$cluster
    data_clust <- cbind(data, mem)
    x <- as.matrix(data_clust)
    mem <- x[, ncol(x), drop=FALSE]
    data_full <- x
    x <- x[, -ncol(x), drop=FALSE]
    C <- nrow(x)
    N <- ncol(x)
    mu <- vapply(seq(1, K), function(k) mean(x[mem == k, ]),numeric(1))
    sigma <- vapply(seq(1, K), function(k) sd(x[mem == k, ]),numeric(1))
    term <- vapply(seq(1, K), function(k) (mu[k]*(1-mu[k])/(sigma[k]^2)) - 1,
                   numeric(1))
    alpha <- mu * term
    beta <- (1 - mu) * term
    alpha <- sort(alpha)
    delta <- sort(beta, decreasing = TRUE)
    alpha[2] <- 1
    delta[2] <- 1
    phi <- list(sp_1 = alpha, sp_2 = delta)
    return(list(A = A, tau = tau, phi = phi))
}



threshold_values <- function(data, tau, phi) {
    threshold_func <- function(data_x, alpha, delta, tau, i,upper_lower) {
        num <- tau[i] * dbeta(data_x, alpha[i], delta[i])
        deno_list <- lapply(seq_along(alpha)[-i], function(j) {
            tau[j] * dbeta(data_x, alpha[j], delta[j])
        })
        deno <- Reduce(`+`, deno_list)
        r <- num/deno
        index <- which(r >= 1)
        if(upper_lower==1){
        l <- data_x[min(index)]}else{
        l <- data_x[max(index)]
        }
        return(l)
    }
    data_x <- sort(data[, 1])
    mode <- as.vector((phi$sp_1 - 1)/(phi$sp_1 + phi$sp_2 - 2))
    cluster <- c(1, 2, 3)
    mode_vec <- cbind(mode, cluster)
    mode_vec <- as.data.frame(mode_vec)
    mode_vec <- mode_vec[order(mode_vec$mode),, drop=FALSE ]
    hypo <- as.numeric(mode_vec[1, 2])
    hyper <- as.numeric(mode_vec[3, 2])
    th_1 <- threshold_func(data_x, phi$sp_1, phi$sp_2, tau, hypo,2)
    th_2 <- threshold_func(data_x, phi$sp_1, phi$sp_2, tau, hyper,1)
    th_vec <- c(th_1, th_2)
    th_new1 <- unique(round(th_vec, 3))

    return(list(thresholds = th_new1))
}

kernel_density_plot <- function(data, hidden_states, auc, C,
                                treatment_group, title_text) {
    data_ggplot <- as.data.frame(data)
    data_ggplot$mem_final <- as.factor(hidden_states)
    colnames(data_ggplot)[length(data_ggplot)] <- "Cluster"
    cols <- ncol(data_ggplot)
    rows <- nrow(data_ggplot)
    data_matrix <- as.matrix(data_ggplot[, seq(1, (cols - 1)), drop=FALSE])
    data_new <- as.vector(data_matrix)
    col_names <- colnames(data_ggplot)
    col_len <- length(col_names)
    Cluster <- vector()
    Patient_sample <- vector()
    result_list <- lapply(seq(1, (col_len - 1)), function(i) {
        temp <- gsub("_", " ", col_names[i])
        ps_names <- rep(temp, rows)
        Cluster <- data_ggplot[, cols, drop=FALSE]
        return(list(Patient_sample = ps_names, Cluster = Cluster))
    })
    Patient_sample <- unlist(lapply(result_list, function(x) x$Patient_sample))
    Cluster <- do.call(rbind, lapply(result_list, function(x) x$Cluster))
    data_plot <- data.frame(data_new = data_new, Cluster = Cluster,
                            Patient_sample = Patient_sample)
    colnames(data_plot) <- c("beta_value", "Cluster", "Patient_Sample")
    data_plot$Cluster <- factor(data_plot$Cluster, levels = auc$State)
    data_plot$Cluster_full <- as.factor(data_plot$Cluster)
    levels(data_plot$Cluster_full) <- auc$label
    color_length <- col_len - 1
    colours_kd <- (seq_gradient_pal(low = "#FFC20A", high = "#0C7BDC",
                                 space = "Lab"))(seq(1, color_length)/
                                                     color_length)
    plot_graph <- ggplot(data_plot) + geom_density(aes(x = beta_value,
                                                       color=Patient_Sample))+
        xlab("Beta Value") +
        ylab("Density") + scale_color_manual("DNA Samples", values = colours_kd) +
        facet_wrap(~factor(Cluster_full, levels = auc$label),
                   scales = "free_y", ) +
        theme(axis.title.x = element_text(size = 10),
              axis.title.y = element_text(size = 10)) +
        ggtitle(title_text)
    cluster_size <- table(hidden_states)
    y <- vapply(auc$State, function(i) {
        max_density <- max(density(data_plot[data_plot$Cluster == i, 1])$y)
        return(max_density)
    }, numeric(1))
    f_labels <- data.frame(Cluster_full = levels(data_plot$Cluster_full),
                           label=as.vector(round((cluster_size[auc$State]/C),
                                                 3)), x = 0.8, y = y)
    plot_graph <- plot_graph + geom_text( aes(x = x, y = y,label = label),
                                          data = f_labels)
    return(plot_graph)
}

fitted_density_plot <- function(phi, hidden_states, auc, R, K, C,
                                treatment_group = NULL, title_text = NULL) {
    vec_C <- 1001
    alpha <- t(phi$sp_1)
    delta <- t(phi$sp_2)
    tau <- round(as.vector(table(hidden_states)/length(hidden_states)), 3)
     vec_x <- seq(0.001, 0.999, length = vec_C)
    combinations <- expand.grid(i = seq(1, R), j = seq(1, K))
    process_combination <- function(i, j, vec_x, tau, alpha, delta,
                                    treatment_group) {
        tmp_vec <- vapply(vec_x, function(x) {
            tau[j] * dbeta(x, alpha[j, i], delta[j, i])
        }, numeric(1))
        data.frame(
            beta = vec_x,
            density = tmp_vec,
            cluster = rep(j, times = length(tmp_vec)),
            sample = rep(treatment_group[i], times = length(tmp_vec))
        )
    }
    result_df <- do.call(rbind, apply(combinations, 1, function(row) {
        process_combination(row["i"], row["j"], vec_x, tau, alpha, delta,
                            treatment_group)
    }))
    beta_vec <- result_df$beta
    density_vec <- result_df$density
    cluster_vec <- result_df$cluster
    sample_vec <- result_df$sample

    df_new_tmp <- data.frame(beta_vec = beta_vec, density_vec = density_vec,
                             cluster_vec= cluster_vec,sample_vec = sample_vec)
    df_new_tmp$cluster_vec <-factor(df_new_tmp$cluster_vec, levels= auc$State)
    df_new_tmp$Cluster_full <- as.factor(df_new_tmp$cluster_vec)
    levels(df_new_tmp$Cluster_full) <- auc$label
    color_length_fd <- R
    colours_fd <- (seq_gradient_pal(low = "#FFC20A", high = "#0C7BDC",
                                 space = "Lab"))(seq(1, color_length_fd)/
                                                     color_length_fd)
    cluster_size <- table(hidden_states)
    plot_graph <- ggplot(df_new_tmp, aes(x = beta_vec, y = density_vec,
                                         color = sample_vec)) + geom_line() +
        scale_color_manual(values = colours_fd) +
        facet_wrap(~factor(Cluster_full, levels = auc$label),
                   scales = "free_y") + labs(color = "Treatment Groups",
                                             x = "Beta value",
                                             y = "Density")+ggtitle(title_text)
    y <- vapply(auc$State, function(i) {
        max_value <- max(df_new_tmp[df_new_tmp$cluster_vec == i, 2])
        return(max_value)
    }, numeric(1))
    f_labels <- data.frame(Cluster_full = levels(df_new_tmp$Cluster_full),
                           label=as.vector(round((cluster_size[auc$State]/C),
                                                 3)), x = 0.92, y = y)
    plot_graph <- plot_graph + geom_text(data = f_labels, aes(x = x, y = y,
                                                              label = label,
                                                              color = NA),
                                         show.legend = FALSE,fontface="bold")
    return(plot_graph)
}

uncertainty_plots <- function(data_comp, z, hidden_states, C,
                              uncertainty_threshold, title_text) {
    labels <- c(max_uc = "Maximum uncertainty")
    z_chr <- z[rownames(z) %in% data_comp$IlmnID,, drop=FALSE ]
    classification_final <- apply(z_chr, 1, max)
    uncertainty <- 1 - classification_final
    tau <- (table(hidden_states))/C
    unc_df <- cbind(uncertainty, hidden_states)
    unc_df <- as.data.frame(unc_df)
    colnames(unc_df) <- c("Uncertainty", "Cluster")
    unc_df$Cluster <- as.factor(unc_df$Cluster)
    unc_df_sorted <- unc_df[order(unc_df$Cluster), , drop=FALSE]
    max_unc <- uncertainty_threshold
    h <- max_unc + 0.015
    max_uncertainty <- data.frame(yintercept = h, max_uncertainty = factor(h))
    plot_graph <- ggplot(unc_df_sorted, aes(x = Cluster, y = Uncertainty)) +
        geom_boxplot() + xlab("Hidden States") +
        ggtitle(title_text) + coord_cartesian(ylim = c(0, 1)) +
        geom_hline(aes(yintercept = max_unc,
                       linetype = paste("Uncertainty threshold = ",
                                        max_unc)), color = "black") +
        scale_linetype_manual(name = "", values = 2,
                              guide = guide_legend(override.aes =
                                                       list(color = "black")))
    return(plot_graph)
}

se<-function(SummarizedExperiment)
{
    SummarizedExperiment <- SummarizedExperiment
    if (!is(SummarizedExperiment, "RangedSummarizedExperiment")) {
        stop("'SummarizedExperiment' must be a
             RangedSummarizedExperiment object")
    }
    return(SummarizedExperiment)
}
anno_data<-function(annotatedData)
{
    if(is.null(annotatedData))
        annotatedData <- DataFrame(matrix(0, nrow=0,ncol=0))
    return(annotatedData)
}

beta_value_numeric<-function(data,anno_file)
{
    final_subset <- subset(data, select = -IlmnID)
    row <- nrow(final_subset)
    col <- ncol(final_subset)
    numeric_vec <- as.numeric(as.matrix(final_subset))
    final_subset <- matrix(numeric_vec, ncol = col)
    if (abs(max(final_subset) - 1) < 1e-06) {
        max_f <- max(final_subset[final_subset != max(final_subset)])
    } else { max_f <- max(final_subset) }
    if (abs(min(final_subset) - 0) < 1e-06) {
        min_f <- min(final_subset[final_subset != min(final_subset)])
    } else { min_f <- min(final_subset)}
    final_subset[final_subset > max_f] <- max_f
    final_subset[final_subset < min_f] <- min_f
    cols <- which(colnames(data) != "IlmnID")
    data[, cols] <- final_subset
    data <- as.data.frame(data)
    data2<-merge(data,anno_file[,c("IlmnID","CHR","MAPINFO"), drop=FALSE],
                by="IlmnID")
    col_order <- which(colnames(data2) == "IlmnID" | colnames(data2) == "CHR"|
                           colnames(data2) == "MAPINFO")
    col_order2 <- which(colnames(data2) != "IlmnID" & colnames(data2)!="CHR" &
                            colnames(data2) != "MAPINFO")
    data3 <- data2[, c(col_order, col_order2), drop=FALSE]
    return(data3)
}
