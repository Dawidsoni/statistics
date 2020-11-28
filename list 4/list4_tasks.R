# Statistics - List 4

set.seed(100)

get_distributions_colnames <- function () {
  return(c(
    "norm_mean0_sd1", "norm_mean1_sd1", "norm_mean0_sd2", "norm_mean1_sd2", "logis_mean0_scale1", "logis_mean1_scale1",
    "logis_mean0_scale2", "logis_mean1_scale2", "cauchy_mean0_scale1", "cauchy_mean1_scale1", "cauchy_mean0_scale2",
    "cauchy_mean1_scale2"
  ))
}

get_experiments_pairs <- function () {
  experiments_df <- data.frame(matrix(nrow = 12, ncol = 2))
  colnames(experiments_df) <- c("distribution1", "distribution2")
  experiments_df["distribution1"] <- c(
    "norm_mean0_sd1", "norm_mean0_sd1", "norm_mean0_sd1", "norm_mean0_sd1", "logis_mean0_scale1", "logis_mean0_scale1",
    "logis_mean0_scale1", "logis_mean0_scale1", "cauchy_mean0_scale1", "cauchy_mean0_scale1", "cauchy_mean0_scale1",
    "cauchy_mean0_scale1"
  )
  experiments_df["distribution2"] <- c(
    "norm_mean0_sd1", "norm_mean1_sd1", "norm_mean0_sd2", "norm_mean1_sd2", "logis_mean0_scale1", "logis_mean1_scale1",
    "logis_mean0_scale2", "logis_mean1_scale2", "cauchy_mean0_scale1", "cauchy_mean1_scale1", "cauchy_mean0_scale2",
    "cauchy_mean1_scale2"
  )
  return(experiments_df)
}

get_distributions_means <- function () {
  means_df <- data.frame(matrix(nrow = 1, ncol = 12))
  colnames(means_df) <- get_distributions_colnames()
  means_df["norm_mean0_sd1"] <- 0.0
  means_df["norm_mean1_sd1"] <- 1.0
  means_df["norm_mean0_sd2"] <- 0.0
  means_df["norm_mean1_sd2"] <- 1.0
  means_df["logis_mean0_scale1"] <- 0.0
  means_df["logis_mean1_scale1"] <- 1.0
  means_df["logis_mean0_scale2"] <- 0.0
  means_df["logis_mean1_scale2"] <- 1.0
  means_df["cauchy_mean0_scale1"] <- 0.0
  means_df["cauchy_mean1_scale1"] <- 1.0
  means_df["cauchy_mean0_scale2"] <- 0.0
  means_df["cauchy_mean1_scale2"] <- 1.0
  return(means_df)
}

get_distributions_variances <- function () {
  vars_df <- data.frame(matrix(nrow = 1, ncol = 12))
  colnames(vars_df) <- get_distributions_colnames()
  vars_df["norm_mean0_sd1"] <- 1.0
  vars_df["norm_mean1_sd1"] <- 1.0
  vars_df["norm_mean0_sd2"] <- 4.0
  vars_df["norm_mean1_sd2"] <- 4.0
  vars_df["logis_mean0_scale1"] <- pi ^ 2 / 3.0
  vars_df["logis_mean1_scale1"] <- pi ^ 2 / 3.0
  vars_df["logis_mean0_scale2"] <- 4.0 * pi ^ 2 / 3.0
  vars_df["logis_mean1_scale2"] <- 4.0 * pi ^ 2 / 3.0
  vars_df["cauchy_mean0_scale1"] <- 1.0
  vars_df["cauchy_mean1_scale1"] <- 1.0
  vars_df["cauchy_mean0_scale2"] <- 4.0
  vars_df["cauchy_mean1_scale2"] <- 4.0
  return(vars_df)
}

generate_distributions_samples <- function (samples_count) {
  samples_df <- data.frame(matrix(nrow = samples_count, ncol = 12))
  colnames(samples_df) <- get_distributions_colnames()
  samples_df["norm_mean0_sd1"] <- rnorm(n = samples_count, mean = 0, sd = 1)
  samples_df["norm_mean1_sd1"] <- rnorm(n = samples_count, mean = 1, sd = 1)
  samples_df["norm_mean0_sd2"] <- rnorm(n = samples_count, mean = 0, sd = 2)
  samples_df["norm_mean1_sd2"] <- rnorm(n = samples_count, mean = 1, sd = 2)
  samples_df["logis_mean0_scale1"] <- rlogis(n = samples_count, location = 0, scale = 1)
  samples_df["logis_mean1_scale1"] <- rlogis(n = samples_count, location = 1, scale = 1)
  samples_df["logis_mean0_scale2"] <- rlogis(n = samples_count, location = 0, scale = 2)
  samples_df["logis_mean1_scale2"] <- rlogis(n = samples_count, location = 1, scale = 2)
  samples_df["cauchy_mean0_scale1"] <- rcauchy(n = samples_count, location = 0, scale = 1)
  samples_df["cauchy_mean1_scale1"] <- rcauchy(n = samples_count, location = 1, scale = 1)
  samples_df["cauchy_mean0_scale2"] <- rcauchy(n = samples_count, location = 0, scale = 2)
  samples_df["cauchy_mean1_scale2"] <- rcauchy(n = samples_count, location = 1, scale = 2)
  return(samples_df)
}

aggregate_means_diff <- function (mean1, mean2) {
  return(mean1 - mean2)
}

aggreagate_vars_quotient <- function (variance1, variance2) {
  return(variance2 / variance1)
}

simulate_coverage_frequency <- function (
  parameters_df, parameters_aggregation_func, interval_func, samples_count = 50, simulations_count = 10000
) {
  experiments_df <- get_experiments_pairs()
  experiment_names <- with(experiments_df, paste(distribution1, distribution2, sep=" "))
  coverage_df <- data.frame(matrix(data = 0, nrow = 1, ncol = dim(experiments_df)[1]))
  colnames(coverage_df) <- experiment_names
  samples_df1 <- generate_distributions_samples(simulations_count * samples_count)
  samples_df2 <- generate_distributions_samples(simulations_count * samples_count)
  for (i in 1:simulations_count) {
    min_sample_index <- samples_count * (i - 1) + 1
    max_sample_index <- samples_count * i
    samples_slice1 <- samples_df1[min_sample_index:max_sample_index, ]
    samples_slice2 <- samples_df2[min_sample_index:max_sample_index, ]
    for (exp_index in 1:dim(experiments_df)[1]) {
      distribution1 <- experiments_df[exp_index, "distribution1"]
      distribution2 <- experiments_df[exp_index, "distribution2"]
      samples1 <- samples_slice1[[distribution1]]
      samples2 <- samples_slice2[[distribution2]]
      real_value <- parameters_aggregation_func(c(parameters_df[[distribution1]]), c(parameters_df[[distribution2]]))
      estimated_interval <- interval_func(distribution1, distribution2, samples1, samples2)
      experiment_name <- paste(distribution1, distribution2, sep=" ")
      if (is.na(estimated_interval[1])) {
        coverage_df[experiment_name] <- NA
      } else if (estimated_interval[1] <= real_value && estimated_interval[2] >= real_value) {
        coverage_df[experiment_name] <- coverage_df[experiment_name] + 1
      }
    }
  }
  coverage_df <- coverage_df / simulations_count
  return(coverage_df[, colSums(is.na(coverage_df)) == 0])
}

means_df <- get_distributions_means()
vars_df <- get_distributions_variances()


# Task 2

get_means_diff_ci <- function (distribution1, distribution2, samples1, samples2) {
  diff_mean <- mean(samples1) - mean(samples2)
  variance1 <- vars_df[[distribution1]]
  variance2 <- vars_df[[distribution2]]
  diff_sd <- sqrt(variance1 / length(samples1) + variance2 / length(samples2))
  quantile <- qnorm(0.975)
  return(c(diff_mean - quantile * diff_sd, diff_mean + quantile * diff_sd))
}

simulate_coverage_frequency(means_df, aggregate_means_diff, get_means_diff_ci)


# Task 4

get_means_diff_unknown_eq_sd_ci <- function (distribution1, distribution2, samples1, samples2) {
  diff_mean <- mean(samples1) - mean(samples2)
  df1 <- length(samples1) - 1
  df2 <- length(samples2) - 1
  joint_df <- df1 + df2
  joint_var <- (df1 * var(samples1) + df2 * var(samples2)) / joint_df
  quantile <- qt(0.975, df = joint_df)
  diff_sd <- sqrt(joint_var * (1 / length(samples1) + 1 / length(samples2)))
  return(c(diff_mean - quantile * diff_sd, diff_mean + quantile * diff_sd))
}

simulate_coverage_frequency(means_df, aggregate_means_diff, get_means_diff_unknown_eq_sd_ci)


# Task 6

get_means_diff_unknown_neq_sd_ci <- function (distribution1, distribution2, samples1, samples2) {
  diff_mean <- mean(samples1) - mean(samples2)
  n1 <- length(samples1)
  n2 <- length(samples2)
  diff_sd <- sqrt(var(samples1) / n1 + var(samples2) / n2)
  df_enumerator <- (var(samples1) / n1 + var(samples2) / n2) ^ 2
  df_denominator <- var(samples1) ^ 2 / (n1 ^ 2 * (n1 - 1)) + var(samples2) ^ 2 / (n2 ^ 2 * (n2 - 1))
  joint_df <- ceiling(df_enumerator / df_denominator)
  quantile <- qt(0.975, df = joint_df)
  return(c(diff_mean - quantile * diff_sd, diff_mean + quantile * diff_sd))
}

simulate_coverage_frequency(means_df, aggregate_means_diff, get_means_diff_unknown_neq_sd_ci)


# Task 8

get_vars_quotient_ci <- function (distribution1, distribution2, samples1, samples2) {
  var1 <- 1 / length(samples1) * sum((samples1 - mean(samples1)) ^ 2)
  var2 <- 1 / length(samples2) * sum((samples2 - mean(samples2)) ^ 2)
  vars_quotient <- var2 / var1
  quantile1 <- qf(0.025, df1 = length(samples1), df2 = length(samples2))
  quantile2 <- qf(0.975, df1 = length(samples1), df2 = length(samples2))
  return(c(quantile1 * vars_quotient, quantile2 * vars_quotient))
}

simulate_coverage_frequency(vars_df, aggreagate_vars_quotient, get_vars_quotient_ci)


# Task 10

get_vars_quotient_uknown_sd_ci <- function (distribution1, distribution2, samples1, samples2) {
  var1 <- var(samples1)
  var2 <- var(samples2)
  vars_quotient <- var2 / var1
  quantile1 <- qf(0.025, df1 = length(samples1) - 1, df2 = length(samples2) - 1)
  quantile2 <- qf(0.975, df1 = length(samples1) - 1, df2 = length(samples2) - 1)
  return(c(quantile1 * vars_quotient, quantile2 * vars_quotient))
}

simulate_coverage_frequency(vars_df, aggreagate_vars_quotient, get_vars_quotient_uknown_sd_ci)


# Task 11

simulate_coverage_frequency(means_df, aggregate_means_diff, get_means_diff_ci, samples_count = 20)
simulate_coverage_frequency(means_df, aggregate_means_diff, get_means_diff_unknown_eq_sd_ci, samples_count = 20)
simulate_coverage_frequency(means_df, aggregate_means_diff, get_means_diff_unknown_neq_sd_ci, samples_count = 20)
simulate_coverage_frequency(vars_df, aggreagate_vars_quotient, get_vars_quotient_ci, samples_count = 20)
simulate_coverage_frequency(vars_df, aggreagate_vars_quotient, get_vars_quotient_uknown_sd_ci, samples_count = 20)

simulate_coverage_frequency(means_df, aggregate_means_diff, get_means_diff_ci, samples_count = 100)
simulate_coverage_frequency(means_df, aggregate_means_diff, get_means_diff_unknown_eq_sd_ci, samples_count = 100)
simulate_coverage_frequency(means_df, aggregate_means_diff, get_means_diff_unknown_neq_sd_ci, samples_count = 100)
simulate_coverage_frequency(vars_df, aggreagate_vars_quotient, get_vars_quotient_ci, samples_count = 100)
simulate_coverage_frequency(vars_df, aggreagate_vars_quotient, get_vars_quotient_uknown_sd_ci, samples_count = 100)
