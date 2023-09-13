globalVariables(c("x", "Cluster"))
#' @keywords internal
#' @import ggplot2
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
    mem <- x[, ncol(x)]
    data_full <- x
    x <- x[, -ncol(x)]
    C <- nrow(x)
    N <- ncol(x)
    mu <- vector("numeric", K)
    sigma <- vector("numeric", K)
    sigma_sq <- vector("numeric", K)
    alpha <- vector("numeric", K)
    beta <- vector("numeric", K)
    term <- vector("numeric", K)
    for (k in seq(1, K)) {
        mu[k] <- mean(x[mem == k, ])
        sigma[k] <- sd(x[mem == k, ])
        term[k] <- (mu[k] * (1 - mu[k])/(sigma[k]^2)) - 1
        alpha[k] <- mu[k] * term[k]
        beta[k] <- (1 - mu[k]) * term[k]
    }
    alpha <- sort(alpha)
    delta <- sort(beta, decreasing = TRUE)
    alpha[2] <- 1
    delta[2] <- 1
    phi <- list(sp_1 = alpha, sp_2 = delta)
    return(list(A = A, tau = tau, phi = phi))
}

BaumWelch_th <- function(data, M = 3, N, R, n.iter = 100, seed = NULL) {
    if (M != 3) {
        stop("M cannot be more or less than 3 as this function tries to
             identify the threshold between the 3 methylation states. ")
    } else if (any(N < 1)) { stop("N cannot be
    less than 1 as one or more DNA replicates need to be analysed.")
    } else if (M!=round(M)|N!=round(N)){
        stop("M and N has to be whole numbers.")}
    K <- M
    C <- nrow(data)
    trained_params <- initialize_parameters_th(data, M, N, R, seed)
    oldlogL <- 1
    BW_limit_accuracy <- 1e-05
    logL_vec <- vector()
    logL_vec <- oldlogL
    R <- 1
    for (i in seq(1, n.iter)) {
        probabilities <- list()
        for (n in seq(1, N)) {
            probabilities[[n]] <-matrix(data[, n], ncol = K, nrow =nrow(data))
            probabilities[[n]] <- t(apply(X = probabilities[[n]],
                                        MARGIN = 1, FUN = dbeta,
                                        shape1 = trained_params$phi$sp_1,
                                        shape2 = trained_params$phi$sp_2))}
        prob <- matrix(probabilities[[1]], nrow = C, ncol = K)
        for (k in 2:(N)) {
            temp <- as.matrix(probabilities[[k]])
            prob <- prob * temp}
        probabilities <- matrix(0, ncol = K, nrow = C)
        probabilities <- prob
        forward_alpha <- forward(probabilities, trained_params)
        backward_beta <- backward(probabilities, trained_params)
        log_alpha <- forward_alpha$log_alpha
        logL_alpha_T <- forward_alpha$scaled_logL
        log_beta <- backward_beta$log_beta
        logL_beta_1 <- backward_beta$scaled_logL
        logL_alpha_beta_T <- 0
        mid_t <- round(C/2)
        t <- max(log_alpha[mid_t, ])
        for (i in seq(1, K)) {
            logL_alpha_beta_T <- logL_alpha_beta_T +
                exp(log_alpha[mid_t, i] + log_beta[mid_t, i])}
        logL_alpha_beta_T <- log(logL_alpha_beta_T)
        logL <- logL_alpha_beta_T
        if (logL == -Inf | logL == Inf | is.na(logL)) {
            logL <- logL_alpha_T}
        if (logL == -Inf | logL == Inf | is.na(logL)) {
            logL <- logL_beta_1}
        eta <- matrix(c(0), ncol = K, nrow = C)
        eta <- exp(log_alpha + log_beta - logL)
        xi <- array(NA, dim = c((C - 1), K, K))
        alpha_c <- log_alpha[seq(1, C - 1), ]
        beta_c <- log_beta[2:C, ]
        proba <- log(probabilities[2:C, ])
        for (j in seq(1, K)) {
            for (i in seq(1, K)) {
                xi[, i, j] <- exp(alpha_c[, i] + log(trained_params$A[i, j]) +
                                    proba[, j] + beta_c[, j] - logL)}}
        tau <- eta[1, ]
        A <- matrix(0, K, K)
        A1 <- apply(xi, c(2, 3), sum)
        A2 <- apply(A1, FUN = sum, MARGIN = 1)
        A <- A1/A2
        y1 <- vector("logical", K)
        y2 <- vector("logical", K)
        al_new2 <- vector("logical", K)
        be_new2 <- vector("logical", K)
        term <- vector("logical", K)
        if (N == 1) {
            y1 <- apply((eta * (log(data[, seq(1,N)]))), FUN= sum,MARGIN = 2)/
                (N * apply(eta, FUN = sum, MARGIN = 2))
            y2 <- apply((eta * (log(1 - data[, seq(1,N)]))),FUN=sum,MARGIN=2)/
                (N * apply(eta, FUN = sum, MARGIN = 2))
        } else {
            y1<-apply((eta*rowSums(log(data[,seq(1,N)]))),FUN=sum, MARGIN = 2)/
                (N * apply(eta, FUN = sum, MARGIN = 2))
            y2<-apply((eta*rowSums(log(1-data[,seq(1,N)]))),FUN=sum,MARGIN=2)/
                (N * apply(eta, FUN = sum, MARGIN = 2))
        }
        term <- ((exp(-y1) - 1) * (exp(-y2) - 1)) - 1
        al_new2 <- 0.5 + (0.5 * exp(-y2)/term)
        be_new2 <- (0.5 * exp(-y2) * (exp(-y1) - 1))/term
        phi <- list(sp_1 = al_new2, sp_2 = be_new2)
        trained_params$A <- A
        trained_params$tau <- tau
        trained_params$phi <- phi
        logL_vec <- c(logL_vec, logL)
        diff_logL <- (logL - oldlogL)/oldlogL
        if (diff_logL > 0 & diff_logL < BW_limit_accuracy) {
            reached_limit_of_accuracy <- TRUE; break}
        oldlogL <- logL};  log_vec <- logL_vec[-1]
    return(list(A = A, tau = tau, phi = phi, log_vec = logL_vec, z = eta))}

Viterbi_th <- function(data, M = 3, N = 4, R = 1, tau, A, phi) {
    K <- M
    C <- nrow(data)
    probabilities <- list()
    for (n in seq(1, N)) {
        probabilities[[n]] <- matrix(data[, n], ncol = K, nrow = nrow(data))
        probabilities[[n]] <- t(apply(X = probabilities[[n]],
                                    MARGIN = 1, FUN = dbeta,
                                    shape1 = phi$sp_1, shape2 = phi$sp_2))
    }
    prob <- matrix(probabilities[[1]], nrow = C, ncol = K)
    for (k in 2:(N)) {
        temp <- as.matrix(probabilities[[k]])
        prob <- prob * temp
    }
    probabilities <- matrix(0, ncol = K, nrow = C)
    probabilities <- prob
    omega <- matrix(0, ncol = K, nrow = C)
    prob <- tau * probabilities[1, ]
    omega[1, ] <- prob/sum(prob)
    for (c in 2:C) {
        prob <- apply(omega[c - 1, ] * A, 2, max) * probabilities[c, ]
        omega[c, ] <- prob/sum(prob)
    }
    decoding <- numeric(C)
    decoding[C] <- which.max(omega[C, ])
    for (c in (C - 1):1) {
        decoding[c] <- which.max(A[, decoding[c + 1]] * omega[c, ])
    }
    return(decoding)
}

threshold_values <- function(data, tau, phi) {
    threshold_func_low <- function(data_x, alpha, delta, tau, i) {
        num <- tau[i] * dbeta(data_x, alpha[i], delta[i])
        deno <- 0
        for (j in seq(1, length(alpha))) {
            if (j != i) {
                deno <- deno + (tau[j] * dbeta(data_x,
                                                        alpha[j], delta[j]))
            }
        }
        r <- num/deno
        index <- which(r >= 1)
        l <- data_x[min(index)]
        return(l)
    }
    threshold_func_upper <- function(data_x, alpha, delta, tau, i) {
        num <- tau[i] * dbeta(data_x, alpha[i], delta[i])
        deno <- 0
        for (j in seq(1, length(alpha))) {
            if (j != i) {
                deno <- deno + (tau[j] * dbeta(data_x, alpha[j], delta[j]))
            }
        }
        r <- num/deno
        index <- which(r >= 1)
        u <- data_x[max(index)]
        return(u)
    }
    data_x <- sort(data[, 1])
    mode <- (phi$sp_1 - 1)/(phi$sp_1 + phi$sp_2 - 2)
    cluster <- c(1, 2, 3)
    mode_vec <- cbind(mode, cluster)
    mode_vec <- as.data.frame(mode_vec)
    mode_vec <- mode_vec[order(mode_vec$mode), ]
    hypo <- as.numeric(mode_vec[1, 2])
    hyper <- as.numeric(mode_vec[3, 2])
    th_1 <- threshold_func_upper(data_x, phi$sp_1, phi$sp_2, tau, hypo)
    th_2 <- threshold_func_low(data_x, phi$sp_1, phi$sp_2, tau, hyper)
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
    data_matrix <- as.matrix(data_ggplot[, seq(1, (cols - 1))])
    data_new <- as.vector(data_matrix)
    col_names <- colnames(data_ggplot)
    col_len <- length(col_names)
    Cluster <- vector()
    Patient_sample <- vector()
    for (i in seq(1, (col_len - 1))) {
        temp <- gsub("_", " ", col_names[i])
        ps_names <- rep(temp, rows)
        Patient_sample <- c(Patient_sample, ps_names)
        Cluster <- c(Cluster, data_ggplot[, cols])
    }
    data_plot <- data.frame(data_new = data_new, Cluster = Cluster,
                            Patient_sample = Patient_sample)
    colnames(data_plot) <- c("beta_value", "Cluster", "Patient_Sample")
    data_plot$Cluster <- factor(data_plot$Cluster, levels = auc$State)
    data_plot$Cluster_full <- as.factor(data_plot$Cluster)
    levels(data_plot$Cluster_full) <- auc$label
    color_length <- col_len - 1
    colours <- (seq_gradient_pal(low = "#FFC20A", high = "#0C7BDC",
                                        space = "Lab"))(seq(1, color_length)/
                                                            color_length)
    plot_graph <- ggplot(data_plot) + geom_density(aes(x = beta_value,
                                                        color=Patient_Sample))+
        xlab("Beta Value") +
        ylab("Density") + scale_color_manual("DNA Samples", values = colours) +
        facet_wrap(~factor(Cluster_full, levels = auc$label),
        scales = "free_y", ) + theme(axis.title.x = element_text(size = 10),
                                    axis.title.y = element_text(size = 10)) +
        ggtitle(title_text)
    cluster_size <- table(hidden_states)
    y <- vector()
    for (i in auc$State) {
        y[i]<-max(density(data_plot[data_plot$Cluster==i,1])[["y"]])}
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
    alpha <- t(phi$sp_1); delta <- t(phi$sp_2)
    tau <- round(as.vector(table(hidden_states)/length(hidden_states)), 3)
    density_vec <- vector(); cluster_vec <- vector()
    sample_vec <- vector(); beta_vec <- vector()
    vec_x <- seq(0.001, 0.999, length = vec_C)
    for (i in seq(1, R)) {
        for (j in seq(1, K)) {
            tmp_vec <- vapply(vec_x, function(x) {
                tau[j]*(dbeta(x,alpha[j,i],delta[j, i]))}, numeric(1))
            beta_vec <- c(beta_vec, vec_x)
            density_vec <- c(density_vec, tmp_vec)
            cluster_vec <- c(cluster_vec, rep(j, times = length(tmp_vec)))
            sample_vec<-c(sample_vec, rep(treatment_group[i],
                                            times = length(tmp_vec)))}}
    df_new_tmp <- data.frame(beta_vec = beta_vec, density_vec = density_vec,
                            cluster_vec= cluster_vec,sample_vec = sample_vec)
    df_new_tmp$cluster_vec <-factor(df_new_tmp$cluster_vec, levels= auc$State)
    df_new_tmp$Cluster_full <- as.factor(df_new_tmp$cluster_vec)
    levels(df_new_tmp$Cluster_full) <- auc$label
    color_length <- R
    colours <- (seq_gradient_pal(low = "#FFC20A", high = "#0C7BDC",
                                        space = "Lab"))(seq(1, color_length)/
                                                            color_length)
    cluster_size <- table(hidden_states)
    plot_graph <- ggplot(df_new_tmp, aes(x = beta_vec, y = density_vec,
                                        color = sample_vec)) + geom_line() +
        scale_color_manual(values = colours) +
        facet_wrap(~factor(Cluster_full, levels = auc$label),
                    scales = "free_y") + labs(color = "Treatment Groups",
                                            x = "Beta value",
        y = "Density") + ggtitle(title_text)
    y <- vector()
    for (i in auc$State) {
        y[i] <- max(df_new_tmp[df_new_tmp$cluster_vec == i, 2])
    }
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
    z_chr <- z[rownames(z) %in% data_comp$IlmnID, ]
    classification_final <- apply(z_chr, 1, max)
    uncertainty <- 1 - classification_final
    tau <- (table(hidden_states))/C
    unc_df <- cbind(uncertainty, hidden_states)
    unc_df <- as.data.frame(unc_df)
    colnames(unc_df) <- c("Uncertainty", "Cluster")
    unc_df$Cluster <- as.factor(unc_df$Cluster)
    unc_df_sorted <- unc_df[order(unc_df$Cluster), ]
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
