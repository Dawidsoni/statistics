# Statistics - List 2
library(LaplacesDemon)
library(ggplot2)

set.seed(1)

# Task 1

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

count_binomal_estimators_frame <- function(samples_count, sample_size, success_prob, estimators_count=10000) {
  estimators_frame <- data.frame(matrix(nrow = estimators_count, ncol = 2))
  colnames(estimators_frame) <- c("prob_ge3", "prob_ge3_real")
  for (i in 1:estimators_count) {
    generated_samples <- rbinom(samples_count, sample_size, success_prob)
    estimated_p <- mean(generated_samples) / sample_size
    estimators_frame[i, "prob_ge3"] <- pbinom(q = 2, sample_size, estimated_p, lower.tail = FALSE)
    estimators_frame[i, "prob_ge3_real"] <- pbinom(q = 2, sample_size, success_prob, lower.tail = FALSE)
  }
  return(estimators_frame)
}


binomal_frame1 <- count_binomal_estimators_frame(samples_count = 50, sample_size = 5, success_prob = 0.1)
round(count_estimators_statistics(binomal_frame1), digits = 6)

binomal_frame2 <- count_binomal_estimators_frame(samples_count = 50, sample_size = 5, success_prob = 0.3)
round(count_estimators_statistics(binomal_frame2), digits = 6)

binomal_frame3 <- count_binomal_estimators_frame(samples_count = 50, sample_size = 5, success_prob = 0.5)
round(count_estimators_statistics(binomal_frame3), digits = 6)

binomal_frame4 <- count_binomal_estimators_frame(samples_count = 50, sample_size = 5, success_prob = 0.7)
round(count_estimators_statistics(binomal_frame4), digits = 6)

binomal_frame5 <- count_binomal_estimators_frame(samples_count = 50, sample_size = 5, success_prob = 0.9)
round(count_estimators_statistics(binomal_frame5), digits = 6)

# Task 2

count_poisson_estimators_frame <- function(samples_count, lambda_param, estimators_count=10000) {
  estimators_frame <- data.frame(matrix(nrow = estimators_count, ncol = 22))
  for (col_index in 1:22) {
    density_value <- (col_index - 1) %% 11
    col_name <- if (col_index <= 11) paste0("prob_e", density_value) else paste0("prob_e", density_value, "_real")
    colnames(estimators_frame)[col_index] <- col_name
  }
  for (i in 1:estimators_count) {
    generated_samples <- rpois(samples_count, lambda_param)
    estimated_lambda <- mean(generated_samples)
    for (density_value in 0:10) {
      estimators_frame[i, paste0("prob_e", density_value)] <- dpois(density_value, lambda = estimated_lambda)
      estimators_frame[i, paste0("prob_e", density_value, "_real")] <- dpois(density_value, lambda = lambda_param)
    }
  }
  return(estimators_frame)
}

poisson_frame1 <- count_poisson_estimators_frame(samples_count = 50, lambda_param = 0.5)
round(count_estimators_statistics(poisson_frame1), digits = 6)

poisson_frame2 <- count_poisson_estimators_frame(samples_count = 50, lambda_param = 1.0)
round(count_estimators_statistics(poisson_frame2), digits = 6)

poisson_frame3 <- count_poisson_estimators_frame(samples_count = 50, lambda_param = 2.0)
round(count_estimators_statistics(poisson_frame3), digits = 6)

poisson_frame4 <- count_poisson_estimators_frame(samples_count = 50, lambda_param = 5.0)
round(count_estimators_statistics(poisson_frame4), digits = 6)


# Task 4

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

perform_shape1_of_beta_experiment <- function(shape1, samples_count) {
  fisher_samples <- get_beta_fisher_information_samples(samples_count, shape1 = shape1, shape2 = 1.0)
  fisher_information <- var(fisher_samples)
  print(paste0("Estimated Fisher information: ", round(fisher_information, digits = 3)))
  show_transformed_fisher_graphs(fisher_information, samples_count, shape1 = shape1, shape2 = 1.0)
}

perform_shape1_of_beta_experiment(shape1 = 0.5, samples_count = 50)
perform_shape1_of_beta_experiment(shape1 = 1.0, samples_count = 50)
perform_shape1_of_beta_experiment(shape1 = 2.0, samples_count = 50)
perform_shape1_of_beta_experiment(shape1 = 5.0, samples_count = 50)


# Task 5

count_estimators_frame <- function(samples_func, estimators_count=10000) {
  estimators_frame <- data.frame(matrix(nrow = estimators_count, ncol = 4))
  colnames(estimators_frame) <- c("mean", "median", "weighted_mean", "normal_mean")
  for (i in 1:estimators_count) {
    generated_samples <- samples_func()
    estimators_frame[i, "mean"] <- mean(generated_samples)
    estimators_frame[i, "median"] <- median(generated_samples)
    weights <- seq_along(generated_samples) / sum(seq_along(generated_samples))
    estimators_frame[i, "weighted_mean"] <- sum(generated_samples * weights)
    prev_sequences <- (seq_along(generated_samples) - 1) / length(generated_samples)
    current_sequences <- seq_along(generated_samples) / length(generated_samples)
    normal_weights <- dnorm(qnorm(prev_sequences)) - dnorm(qnorm(current_sequences))
    sorted_samples <- sort(generated_samples)
    estimators_frame[i, "normal_mean"] <- sum(sorted_samples * normal_weights)
  }
  return(estimators_frame)
}

laplace_frame1 <- count_estimators_frame(samples_func = function() rlaplace(n = 50, location = 1.0, scale = 1.0))
round(count_estimators_statistics(laplace_frame1, estimated_paramter = 1.0), digits = 6)

laplace_frame2 <- count_estimators_frame(samples_func = function() rlaplace(n = 50, location = 4.0, scale = 1.0))
round(count_estimators_statistics(laplace_frame2, estimated_paramter = 4.0), digits = 6)

laplace_frame3 <- count_estimators_frame(samples_func = function() rlaplace(n = 50, location = 1.0, scale = 2.0))
round(count_estimators_statistics(laplace_frame3, estimated_paramter = 1.0), digits = 6)


# Task 6

binomal_frame6 <- count_binomal_estimators_frame(samples_count = 20, sample_size = 5, success_prob = 0.1)
round(count_estimators_statistics(binomal_frame1), digits = 6)

binomal_frame7 <- count_binomal_estimators_frame(samples_count = 100, sample_size = 5, success_prob = 0.1)
round(count_estimators_statistics(binomal_frame1), digits = 6)

poisson_frame5 <- count_poisson_estimators_frame(samples_count = 20, lambda_param = 0.5)
round(count_estimators_statistics(poisson_frame5), digits = 6)

poisson_frame6 <- count_poisson_estimators_frame(samples_count = 100, lambda_param = 0.5)
round(count_estimators_statistics(poisson_frame6), digits = 6)

perform_shape1_of_beta_experiment(shape1 = 0.5, samples_count = 20)

perform_shape1_of_beta_experiment(shape1 = 0.5, samples_count = 100)

laplace_frame4 <- count_estimators_frame(samples_func = function() rlaplace(n = 20, location = 1.0, scale = 1.0))
round(count_estimators_statistics(laplace_frame4, estimated_paramter = 1.0), digits = 6)

laplace_frame5 <- count_estimators_frame(samples_func = function() rlaplace(n = 100, location = 1.0, scale = 1.0))
round(count_estimators_statistics(laplace_frame5, estimated_paramter = 1.0), digits = 6)