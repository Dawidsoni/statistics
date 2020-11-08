# Statistics - List 3

set.seed(100)

get_distributions_colnames <- function () {
  return(c(
    "norm_sd1", "norm_sd2", "norm_sd3", "logis_sd1", "logis_sd2", "logis_sd3", "cauchy_sd1", "cauchy_sd2",
    "cauchy_sd3", "exp_lambda1", "exp_lambda1_2", "exp_lambda_1_3", "chi_v1", "chi_v2", "chi_v3"
  ))
}

get_distributions_means <- function () {
  means_df <- data.frame(matrix(nrow = 1, ncol = 15))
  colnames(means_df) <- get_distributions_colnames()
  means_df["norm_sd1"] <- 0.0
  means_df["norm_sd2"] <- 0.0
  means_df["norm_sd3"] <- 0.0
  means_df["logis_sd1"] <- 0.0
  means_df["logis_sd2"] <- 0.0
  means_df["logis_sd3"] <- 0.0
  means_df["cauchy_sd1"] <- 0.0
  means_df["cauchy_sd2"] <- 0.0
  means_df["cauchy_sd3"] <- 0.0
  means_df["exp_lambda1"] <- 1.0
  means_df["exp_lambda1_2"] <- 2.0
  means_df["exp_lambda_1_3"] <- 3.0
  means_df["chi_v1"] <- 1.0
  means_df["chi_v2"] <- 2.0
  means_df["chi_v3"] <- 3.0
  return(means_df)
}

get_distributions_variances <- function () {
  variances_df <- data.frame(matrix(nrow = 1, ncol = 15))
  colnames(variances_df) <- get_distributions_colnames()
  variances_df["norm_sd1"] <- 1.0
  variances_df["norm_sd2"] <- 4.0
  variances_df["norm_sd3"] <- 9.0
  variances_df["logis_sd1"] <- pi ^ 2 / 3.0
  variances_df["logis_sd2"] <- 4.0 * pi ^ 2 / 3.0
  variances_df["logis_sd3"] <- 9.0 * pi ^ 2 / 3.0
  variances_df["cauchy_sd1"] <- 1.0
  variances_df["cauchy_sd2"] <- 4.0
  variances_df["cauchy_sd3"] <- 9.0
  variances_df["exp_lambda1"] <- 1.0
  variances_df["exp_lambda1_2"] <- 4.0
  variances_df["exp_lambda_1_3"] <- 9.0
  variances_df["chi_v1"] <- 2.0
  variances_df["chi_v2"] <- 4.0
  variances_df["chi_v3"] <- 6.0
  return(variances_df)
}

get_distributions_positive_rates <- function () {
  rates_df <- data.frame(matrix(nrow = 1, ncol = 9))
  colnames(rates_df) <- get_distributions_colnames()[1:9]
  rates_df["norm_sd1"] <- 0.5
  rates_df["norm_sd2"] <- 0.5
  rates_df["norm_sd3"] <- 0.5
  rates_df["logis_sd1"] <- 0.5
  rates_df["logis_sd2"] <- 0.5
  rates_df["logis_sd3"] <- 0.5
  rates_df["cauchy_sd1"] <- 0.5
  rates_df["cauchy_sd2"] <- 0.5
  rates_df["cauchy_sd3"] <- 0.5
  return(rates_df)
}

generate_distributions_samples <- function (samples_count) {
  samples_df <- data.frame(matrix(nrow = samples_count, ncol = 15))
  colnames(samples_df) <- get_distributions_colnames()
  samples_df["norm_sd1"] <- rnorm(n = samples_count, mean = 0, sd = 1)
  samples_df["norm_sd2"] <- rnorm(n = samples_count, mean = 0, sd = 2)
  samples_df["norm_sd3"] <- rnorm(n = samples_count, mean = 0, sd = 3)
  samples_df["logis_sd1"] <- rlogis(n = samples_count, location = 0, scale = 1)
  samples_df["logis_sd2"] <- rlogis(n = samples_count, location = 0, scale = 2)
  samples_df["logis_sd3"] <- rlogis(n = samples_count, location = 0, scale = 3)
  samples_df["cauchy_sd1"] <- rcauchy(n = samples_count, location = 0, scale = 1)
  samples_df["cauchy_sd2"] <- rcauchy(n = samples_count, location = 0, scale = 2)
  samples_df["cauchy_sd3"] <- rcauchy(n = samples_count, location = 0, scale = 3)
  samples_df["exp_lambda1"] <- rexp(n = samples_count, rate = 1.0)
  samples_df["exp_lambda1_2"] <- rexp(n = samples_count, rate = 0.5)
  samples_df["exp_lambda_1_3"] <- rexp(n = samples_count, rate = 1.0 / 3.0)
  samples_df["chi_v1"] <- rchisq(n = samples_count, df = 1)
  samples_df["chi_v2"] <- rchisq(n = samples_count, df = 2)
  samples_df["chi_v3"] <- rchisq(n = samples_count, df = 3)
  return(samples_df)
}

simulate_coverage_frequency <- function (parameters_df, interval_func, samples_count = 50, simulations_count = 10000) {
  coverage_df <- data.frame(matrix(data = 0, nrow = 1, ncol = length(parameters_df)))
  colnames(coverage_df) <- colnames(parameters_df)
  samples_df <- generate_distributions_samples(simulations_count * samples_count)
  for (i in 1:simulations_count) {
    min_sample_index <- samples_count * (i - 1) + 1
    max_sample_index <- samples_count * i
    samples_slice <- samples_df[min_sample_index:max_sample_index, ]
    for (column_name in colnames(samples_slice)) {
      if (!(column_name %in% colnames(parameters_df))) {
        next
      }
      real_value <- c(parameters_df[column_name])
      estimated_interval <- interval_func(column_name, c(samples_slice[[column_name]]))
      if (estimated_interval[1] <= real_value && estimated_interval[2] >= real_value) {
        coverage_df[column_name] <- coverage_df[column_name] + 1
      }
    }
  }
  return(coverage_df / simulations_count)
}

means_df <- get_distributions_means()
variances_df <- get_distributions_variances()
positive_rates_df <- get_distributions_positive_rates()


# Task 2

get_mean_confidence_interval <- function (column_name, samples, significance_level = 0.05) {
  middle_point <- mean(samples)
  variance <- c(variances_df[[column_name]])
  deviation <- qnorm(1 - significance_level / 2) * sqrt(variance / length(samples))
  return(c(middle_point - deviation, middle_point + deviation))
}

simulate_coverage_frequency(means_df, get_mean_confidence_interval)


# Task 4

get_mean_unknown_sd_confidence_interval <- function (unused_column_name, samples, significance_level = 0.05) {
  middle_point <- mean(samples)
  deviation <- qt(1 - significance_level / 2, df = length(samples) - 1) * sqrt(var(samples) / length(samples))
  return(c(middle_point - deviation, middle_point + deviation))
}

simulate_coverage_frequency(means_df, get_mean_unknown_sd_confidence_interval)


# Task 6

get_sd_confidence_interval_using_chisq <- function (column_name, samples, significance_level = 0.05) {
  mean <- c(means_df[[column_name]])
  estimated_var <- mean((samples - mean) ^ 2)
  df <- length(samples)
  max_point <- df * estimated_var / qchisq(significance_level / 2, df)
  min_point <-  df * estimated_var / qchisq(1 - significance_level / 2, df)
  return(c(min_point, max_point))
}

get_sd_confidence_interval_using_student <- function (column_name, samples, significance_level = 0.05) {
  mean <- c(means_df[[column_name]])
  var_samples <- (samples - mean) ^ 2
  middle_point <- mean(var_samples)
  df <- length(var_samples) - 1
  deviation <- qt(1 - significance_level / 2, df) * sqrt(var(var_samples) / length(var_samples))
  return(c(middle_point - deviation, middle_point + deviation))
}

simulate_coverage_frequency(variances_df, get_sd_confidence_interval_using_chisq)
simulate_coverage_frequency(variances_df, get_sd_confidence_interval_using_student)


# Task 8

get_sd_unknown_mean_confidence_interval <- function (unused_column_name, samples, significance_level = 0.05) {
  df <- length(samples) - 1
  max_point <- df * var(samples) / qchisq(significance_level / 2, df)
  min_point <-  df * var(samples) / qchisq(1 - significance_level / 2, df)
  return(c(min_point, max_point))
}

simulate_coverage_frequency(variances_df, get_sd_unknown_mean_confidence_interval)


# Task 10

get_positive_rate_confidence_interval <- function (unused_column_name, samples, significance_level = 0.05) {
  estimated_probability <- mean(samples >= 0)
  variance <- estimated_probability * (1 - estimated_probability)
  deviation <- qnorm(1 - significance_level / 2) * sqrt(variance / length(samples))
  return(c(estimated_probability - deviation, estimated_probability + deviation))
}

simulate_coverage_frequency(positive_rates_df, get_positive_rate_confidence_interval)


# Task 11

simulate_coverage_frequency(means_df, get_mean_confidence_interval, samples_count = 20)
simulate_coverage_frequency(means_df, get_mean_confidence_interval, samples_count = 100)

simulate_coverage_frequency(means_df, get_mean_unknown_sd_confidence_interval, samples_count = 20)
simulate_coverage_frequency(means_df, get_mean_unknown_sd_confidence_interval, samples_count = 100)

simulate_coverage_frequency(variances_df, get_sd_confidence_interval_using_chisq, samples_count = 20)
simulate_coverage_frequency(variances_df, get_sd_confidence_interval_using_student, samples_count = 20)
simulate_coverage_frequency(variances_df, get_sd_confidence_interval_using_chisq, samples_count = 100)
simulate_coverage_frequency(variances_df, get_sd_confidence_interval_using_student, samples_count = 100)

simulate_coverage_frequency(variances_df, get_sd_unknown_mean_confidence_interval, samples_count = 20)
simulate_coverage_frequency(variances_df, get_sd_unknown_mean_confidence_interval, samples_count = 100)

simulate_coverage_frequency(positive_rates_df, get_positive_rate_confidence_interval, samples_count = 20)
simulate_coverage_frequency(positive_rates_df, get_positive_rate_confidence_interval, samples_count = 100)
