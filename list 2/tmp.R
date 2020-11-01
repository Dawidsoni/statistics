library(ggplot2)

count_estimators_statistics <- function(estimators_frame, estimated_paramter=NULL, excluded_col_names=c()) {
  statistics_frame <- data.frame(matrix(nrow = 1, ncol = 0))
  if (!is.null(estimated_paramter)) {
    statistics_frame["real_value"] <- estimated_paramter
  }
  for (col_name in colnames(estimators_frame)) {
    if (col_name %in% excluded_col_names) {
      next
    }
    if (grepl(pattern = "_real$", x = col_name)) {
      statistics_frame[col_name] <- mean(estimators_frame[[col_name]])
      next
    }
    real_value_col <- paste0(col_name, "_real")
    real_value <- if (is.null(estimated_paramter)) estimators_frame[[real_value_col]] else estimated_paramter
    statistics_frame[paste0(col_name, '_variance')] <- var(c(estimators_frame[[col_name]]), na.rm = TRUE)
    statistics_frame[paste0(col_name, '_mse')] <- mean(
      (c(estimators_frame[[col_name]]) - real_value) ^ 2, na.rm = TRUE
    )
    statistics_frame[paste0(col_name, '_bias')] <- mean(
      c(estimators_frame[[col_name]]) - real_value, na.rm = TRUE
    )
  }
  return(statistics_frame)
}

get_beta_fisher_information_samples <- function(samples_count, shape1, shape2, estimators_count=10000) {
  fisher_samples <- rep(NA, estimators_count)
  for (i in 1:estimators_count) {
    generated_samples <- rbeta(samples_count, shape1, shape2)
    fisher_samples[[i]] <- sum(log(generated_samples)) + samples_count / shape1
  }
  return(fisher_samples)
}

show_transformed_fisher_graphs <- function(fisher_information, samples_count, shape1, shape2, estimators_count=10000) {
  transformed_samples <- rep(NA, estimators_count)
  for (i in 1:estimators_count) {
    beta_samples <- rbeta(samples_count, shape1, shape2)
    estimated_shape1 <- -samples_count / sum(log(beta_samples))
    transformed_samples[[i]] <- sqrt(samples_count * fisher_information) * (estimated_shape1 - shape1)
  }
  samples_df <- data.frame(sample = transformed_samples)
  stat_mean <- mean(transformed_samples)
  stat_sd <- sd(transformed_samples)
  ggplot(samples_df, aes(x = sample)) +
    geom_histogram(aes(y = ..density..), fill = "blue", color = "black", binwidth = 2.5) +
    stat_function(fun = dnorm, args = list(mean = stat_mean, sd = stat_sd), size = 3.0, n = 1000) +
    ggtitle(paste0("Transformed samples of Fisher information of Beta(shape1 = ", shape1, ", shape2 = ", shape2, ")")) +
    xlab("Transformed Fisher information samples") +
    theme(plot.title = element_text(hjust = 0.5))
  ggplot(samples_df, aes(sample = sample)) +
    geom_qq(color = "darkorange") +
    geom_qq_line(color = "darkblue", size = 1.0) +
    ggtitle(paste0("Transformed samples of Fisher information of Beta(shape1 = ", shape1, ", shape2 = ", shape2, ")")) +
    ylab("Transformed Fisher information samples") +
    theme(plot.title = element_text(hjust = 0.5))
}

perform_shape1_of_beta_experiment <- function(shape1) {
  fisher_samples <- get_beta_fisher_information_samples(samples_count = 50, shape1 = shape1, shape2 = 1.0)
  fisher_information <- var(fisher_samples)
  paste0("Estimated Fisher information: ", round(fisher_information, digits = 3))
  show_transformed_fisher_graphs(fisher_information, samples_count = 50, shape1 = shape1, shape2 = 1.0)
}

perform_shape1_of_beta_experiment(shape1 = 0.5)
perform_shape1_of_beta_experiment(shape1 = 1.0)
perform_shape1_of_beta_experiment(shape1 = 2.0)
perform_shape1_of_beta_experiment(shape1 = 5.0)
