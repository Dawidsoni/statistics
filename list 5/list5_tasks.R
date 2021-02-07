# Statistics - List 5

set.seed(100)
library("dplyr")
library("ggplot2")

calculate_ranks <- function (samples) {
  order_of_samples <- order(samples)
  ranks <- rep(NA, length(samples))
  ranks[order_of_samples] <- seq_along(samples)
  return(ranks)
}

calculate_rank_statistic <- function (samples1, samples2, statistic_func) {
  count1 <- length(samples1)
  count2 <- length(samples2)
  total_count <- count1 + count2
  joint_ranks <- calculate_ranks(c(samples1, samples2))
  ranks1 <- joint_ranks[1:count1]
  ranks2 <- joint_ranks[(count1 + 1):total_count]
  statistic1 <- sum(statistic_func((ranks1 - 0.5) / total_count)) / count1
  statistic2 <- sum(statistic_func((ranks2 - 0.5) / total_count)) / count2
  return(sqrt((count1 * count2) / total_count) * (statistic1 - statistic2))
}

calculate_wilcoxon_statistic <- function (samples1, samples2) {
  statistic_func <- function (x) sqrt(3) * (2 * x - 1)
  rank_statistic <- calculate_rank_statistic(samples1, samples2, statistic_func)
  return(rank_statistic ^ 2)
}

calculate_ab_statistic <- function (samples1, samples2) {
  statistic_func <- function (x) sqrt(48) * (0.25 - abs(x - 0.5))
  rank_statistic <- calculate_rank_statistic(samples1, samples2, statistic_func)
  return(rank_statistic ^ 2)
}

calculate_lepage_statistic <- function (samples1, samples2) {
  wilcoxon_statistic <- calculate_wilcoxon_statistic(samples1, samples2)
  ab_statistic <- calculate_ab_statistic(samples1, samples2)
  return(wilcoxon_statistic + ab_statistic)
}

calculate_ks_statistic <- function (samples1, samples2) {
  count1 <- length(samples1)
  count2 <- length(samples2)
  total_count <- count1 + count2
  max_cdf_diff <- 0.0
  for (i in seq_along(samples1)) {
    cdf1 <- mean(samples1 <= samples1[i])
    cdf2 <- mean(samples2 <= samples1[i])
    max_cdf_diff <- max(max_cdf_diff, abs(cdf1 - cdf2))
  }
  for (i in seq_along(samples2)) {
    cdf1 <- mean(samples1 <= samples2[i])
    cdf2 <- mean(samples2 <= samples2[i])
    max_cdf_diff <- max(max_cdf_diff, abs(cdf1 - cdf2))
  }
  sqrt((count1 * count2) / total_count) * max_cdf_diff
}

create_mean_scaled_distributions_df <- function (samples_count) {
  samples_df <- data.frame(matrix(nrow = samples_count, ncol = 23))
  colnames(samples_df) <- c(
    "norm_mean_0", "norm_mean_0_2", "norm_mean_0_4", "norm_mean_0_6", "norm_mean_0_8", "norm_mean_1", "norm_mean_1_2",
    "norm_mean_1_4", "logis_mean_0", "logis_mean_0_2", "logis_mean_0_4", "logis_mean_0_6", "logis_mean_0_8",
    "logis_mean_1", "logis_mean_1_2", "logis_mean_1_4", "cauchy_mean_0", "cauchy_mean_0_5", "cauchy_mean_1",
    "cauchy_mean_1_5", "cauchy_mean_2", "cauchy_mean_2_5", "cauchy_mean_3"
  )
  samples_df["norm_mean_0"] <- rnorm(samples_count, mean = 0, sd = 1)
  samples_df["norm_mean_0_2"] <- rnorm(samples_count, mean = 0.2, sd = 1)
  samples_df["norm_mean_0_4"] <- rnorm(samples_count, mean = 0.4, sd = 1)
  samples_df["norm_mean_0_6"] <- rnorm(samples_count, mean = 0.6, sd = 1)
  samples_df["norm_mean_0_8"] <- rnorm(samples_count, mean = 0.8, sd = 1)
  samples_df["norm_mean_1"] <- rnorm(samples_count, mean = 1.0, sd = 1)
  samples_df["norm_mean_1_2"] <- rnorm(samples_count, mean = 1.2, sd = 1)
  samples_df["norm_mean_1_4"] <- rnorm(samples_count, mean = 1.4, sd = 1)
  samples_df["logis_mean_0"] <- rlogis(samples_count, location = 0, scale = 1)
  samples_df["logis_mean_0_2"] <- rlogis(samples_count, location = 0.2, scale = 1)
  samples_df["logis_mean_0_4"] <- rlogis(samples_count, location = 0.4, scale = 1)
  samples_df["logis_mean_0_6"] <- rlogis(samples_count, location = 0.6, scale = 1)
  samples_df["logis_mean_0_8"] <- rlogis(samples_count, location = 0.8, scale = 1)
  samples_df["logis_mean_1"] <- rlogis(samples_count, location = 1.0, scale = 1)
  samples_df["logis_mean_1_2"] <- rlogis(samples_count, location = 1.2, scale = 1)
  samples_df["logis_mean_1_4"] <- rlogis(samples_count, location = 1.4, scale = 1)
  samples_df["cauchy_mean_0"] <- rcauchy(samples_count, location = 0.0, scale = 1)
  samples_df["cauchy_mean_0_5"] <- rcauchy(samples_count, location = 0.5, scale = 1)
  samples_df["cauchy_mean_1"] <- rcauchy(samples_count, location = 1.0, scale = 1)
  samples_df["cauchy_mean_1_5"] <- rcauchy(samples_count, location = 1.5, scale = 1)
  samples_df["cauchy_mean_2"] <- rcauchy(samples_count, location = 2.0, scale = 1)
  samples_df["cauchy_mean_2_5"] <- rcauchy(samples_count, location = 2.5, scale = 1)
  samples_df["cauchy_mean_3"] <- rcauchy(samples_count, location = 3.0, scale = 1)
  return(samples_df)
}

create_sd_scaled_distributions_df <- function (samples_count) {
  samples_df <- data.frame(matrix(nrow = samples_count, ncol = 21))
  colnames(samples_df) <- c(
    "norm_sd_1", "norm_sd_1_5", "norm_sd_2", "norm_sd_2_5", "norm_sd_3", "norm_sd_3_5", "norm_sd_4", "logis_sd_1",
    "logis_sd_1_5", "logis_sd_2", "logis_sd_2_5", "logis_sd_3", "logis_sd_3_5", "logis_sd_4", "cauchy_sd_1",
    "cauchy_sd_2", "cauchy_sd_3", "cauchy_sd_4", "cauchy_sd_5", "cauchy_sd_6", "cauchy_sd_7"
  )
  samples_df["norm_sd_1"] <- rnorm(samples_count, mean = 0, sd = 1.0)
  samples_df["norm_sd_1_5"] <- rnorm(samples_count, mean = 0, sd = 1.5)
  samples_df["norm_sd_2"] <- rnorm(samples_count, mean = 0, sd = 2.0)
  samples_df["norm_sd_2_5"] <- rnorm(samples_count, mean = 0, sd = 2.5)
  samples_df["norm_sd_3"] <- rnorm(samples_count, mean = 0, sd = 3.0)
  samples_df["norm_sd_3_5"] <- rnorm(samples_count, mean = 0, sd = 3.5)
  samples_df["norm_sd_4"] <- rnorm(samples_count, mean = 0, sd = 4.0)
  samples_df["logis_sd_1"] <- rlogis(samples_count, location = 0, scale = 1.0)
  samples_df["logis_sd_1_5"] <- rlogis(samples_count, location = 0, scale = 1.5)
  samples_df["logis_sd_2"] <- rlogis(samples_count, location = 0, scale = 2.0)
  samples_df["logis_sd_2_5"] <- rlogis(samples_count, location = 0, scale = 2.5)
  samples_df["logis_sd_3"] <- rlogis(samples_count, location = 0, scale = 3.0)
  samples_df["logis_sd_3_5"] <- rlogis(samples_count, location = 0, scale = 3.5)
  samples_df["logis_sd_4"] <- rlogis(samples_count, location = 0, scale = 4.0)
  samples_df["cauchy_sd_1"] <- rcauchy(samples_count, location = 0, scale = 1.0)
  samples_df["cauchy_sd_2"] <- rcauchy(samples_count, location = 0, scale = 2.0)
  samples_df["cauchy_sd_3"] <- rcauchy(samples_count, location = 0, scale = 3.0)
  samples_df["cauchy_sd_4"] <- rcauchy(samples_count, location = 0, scale = 4.0)
  samples_df["cauchy_sd_5"] <- rcauchy(samples_count, location = 0, scale = 5.0)
  samples_df["cauchy_sd_6"] <- rcauchy(samples_count, location = 0, scale = 6.0)
  samples_df["cauchy_sd_7"] <- rcauchy(samples_count, location = 0, scale = 7.0)
  return(samples_df)
}

create_mean_sd_scaled_distributions_df <- function (samples_count) {
  samples_df <- data.frame(matrix(nrow = samples_count, ncol = 23))
  colnames(samples_df) <- c(
    "norm_mean_0_sd_1", "norm_mean_0_2_sd_1", "norm_mean_0_4_sd_1_5", "norm_mean_0_6_sd_2", "norm_mean_0_8_sd_2_5",
    "norm_mean_1_sd_3", "norm_mean_1_2_sd_3_5", "norm_mean_1_4_sd_4", "logis_mean_0_sd_1", "logis_mean_0_2_sd_1",
    "logis_mean_0_4_sd_1_5", "logis_mean_0_6_sd_2", "logis_mean_0_8_sd_2_5", "logis_mean_1_sd_3",
    "logis_mean_1_2_sd_3_5", "logis_mean_1_4_sd_4", "cauchy_mean_0_sd_1", "cauchy_mean_0_5_sd_2", "cauchy_mean_1_sd_3",
    "cauchy_mean_1_5_sd_4", "cauchy_mean_2_sd_5", "cauchy_mean_2_5_sd_6", "cauchy_mean_3_sd_7"
  )
  samples_df["norm_mean_0_sd_1"] <- rnorm(samples_count, mean = 0, sd = 1.0)
  samples_df["norm_mean_0_2_sd_1"] <- rnorm(samples_count, mean = 0.2, sd = 1)
  samples_df["norm_mean_0_4_sd_1_5"] <- rnorm(samples_count, mean = 0.4, sd = 1.5)
  samples_df["norm_mean_0_6_sd_2"] <- rnorm(samples_count, mean = 0.6, sd = 2.0)
  samples_df["norm_mean_0_8_sd_2_5"] <- rnorm(samples_count, mean = 0.8, sd = 2.5)
  samples_df["norm_mean_1_sd_3"] <- rnorm(samples_count, mean = 1.0, sd = 3.0)
  samples_df["norm_mean_1_2_sd_3_5"] <- rnorm(samples_count, mean = 1.2, sd = 3.5)
  samples_df["norm_mean_1_4_sd_4"] <- rnorm(samples_count, mean = 1.4, sd = 4.0)
  samples_df["logis_mean_0_sd_1"] <- rlogis(samples_count, location = 0, scale = 1.0)
  samples_df["logis_mean_0_2_sd_1"] <- rlogis(samples_count, location = 0.2, scale = 1.0)
  samples_df["logis_mean_0_4_sd_1_5"] <- rlogis(samples_count, location = 0.4, scale = 1.5)
  samples_df["logis_mean_0_6_sd_2"] <- rlogis(samples_count, location = 0.6, scale = 2.0)
  samples_df["logis_mean_0_8_sd_2_5"] <- rlogis(samples_count, location = 0.8, scale = 2.5)
  samples_df["logis_mean_1_sd_3"] <- rlogis(samples_count, location = 1.0, scale = 3.0)
  samples_df["logis_mean_1_2_sd_3_5"] <- rlogis(samples_count, location = 1.2, scale = 3.5)
  samples_df["logis_mean_1_4_sd_4"] <- rlogis(samples_count, location = 1.4, scale = 4.0)
  samples_df["cauchy_mean_0_sd_1"] <- rcauchy(samples_count, location = 0.0, scale = 1.0)
  samples_df["cauchy_mean_0_5_sd_2"] <- rcauchy(samples_count, location = 0.5, scale = 1.0)
  samples_df["cauchy_mean_1_sd_3"] <- rcauchy(samples_count, location = 1.0, scale = 3.0)
  samples_df["cauchy_mean_1_5_sd_4"] <- rcauchy(samples_count, location = 1.5, scale = 4.0)
  samples_df["cauchy_mean_2_sd_5"] <- rcauchy(samples_count, location = 2.0, scale = 5.0)
  samples_df["cauchy_mean_2_5_sd_6"] <- rcauchy(samples_count, location = 2.5, scale = 6.0)
  samples_df["cauchy_mean_3_sd_7"] <- rcauchy(samples_count, location = 3.0, scale = 7.0)
  return(samples_df)
}

simulate_test_statistic_values <- function (
  statistic_func, experiments_df, samples_func, samples_count, simulations_count = 10000
) {
  samples_df1 <- samples_func(simulations_count * samples_count)
  samples_df2 <- samples_func(simulations_count * samples_count)
  statistics_df <- data.frame(matrix(data = NA, nrow = simulations_count, ncol = dim(experiments_df)[1]))
  colnames(statistics_df) <- experiments_df$experiment_name
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
      estimated_statistic <- statistic_func(samples1, samples2)
      experiment_name <- experiments_df[exp_index, "experiment_name"]
      statistics_df[i, experiment_name] <- estimated_statistic
    }
  }
  return(statistics_df)
}


# Task 1

get_normal_distribution_experiment_pairs <- function () {
  experiments_df <- data.frame(matrix(nrow = 1, ncol = 2))
  colnames(experiments_df) <- c("distribution1", "distribution2")
  experiments_df["distribution1"] <- c("norm_mean_0")
  experiments_df["distribution2"] <- c("norm_mean_0")
  experiments_df["experiment_name"] <- with(experiments_df, paste(distribution1, distribution2, sep="#"))
  return(experiments_df)
}

estimate_test_statistic_critical_value <- function (
  statistic_func, experiments_df, samples_func, samples_count, simulations_count = 10000, significance_level = 0.05
) {
  simulated_test_statistics <- simulate_test_statistic_values(
    statistic_func, experiments_df, samples_func, samples_count, simulations_count
  )
  estimated_critical_values <- data.frame(matrix(data = NA, nrow = 1, ncol = dim(simulated_test_statistics)[2]))
  colnames(estimated_critical_values) <- colnames(simulated_test_statistics)
  critical_value_index <- floor((1 - significance_level) * simulations_count)
  for (column in colnames(estimated_critical_values)) {
    ordered_values <- sort(simulated_test_statistics[[column]])
    estimated_critical_values[column] <- ordered_values[critical_value_index]
  }
  return(estimated_critical_values)
}

analyse_critical_values <- function (samples_count, significance_level = 0.05) {
  norm_experiments_df <- get_normal_distribution_experiment_pairs()
  estimated_wilcoxon_critical_value <- estimate_test_statistic_critical_value(
    calculate_wilcoxon_statistic, norm_experiments_df, create_mean_scaled_distributions_df, samples_count,
    significance_level = significance_level
  )
  print(paste0("Estimated Wilcoxon test critical value: ", round(estimated_wilcoxon_critical_value, digits = 3)))
  print(paste0("Real Wilcoxon test critical value: ", round(qchisq(1 - significance_level, df = 1), digits = 3)))
  estimated_ab_critical_value <- estimate_test_statistic_critical_value(
    calculate_ab_statistic, norm_experiments_df, create_mean_scaled_distributions_df, samples_count,
    significance_level = significance_level
  )
  print(paste0("Estimated AB test critical value: ", round(estimated_ab_critical_value, digits = 3)))
  print(paste0("Real AB test critical value: ", round(qchisq(1 - significance_level, df = 1), digits = 3)))
  estimated_lepage_critical_value <- estimate_test_statistic_critical_value(
    calculate_lepage_statistic, norm_experiments_df, create_mean_scaled_distributions_df, samples_count,
    significance_level = significance_level
  )
  print(paste0("Estimated Lepage test critical value: ", round(estimated_lepage_critical_value, digits = 3)))
  print(paste0("Real Lepage test critical value: ", round(qchisq(1 - significance_level, df = 2), digits = 3)))
  estimated_ks_critical_value <- estimate_test_statistic_critical_value(
    calculate_ks_statistic, norm_experiments_df, create_mean_scaled_distributions_df, samples_count,
    significance_level = significance_level
  )
  print(paste0("Estimated KS test critical value: ", round(estimated_ks_critical_value, digits = 3)))
  print(paste0("Real KS test critical value: ", round(sqrt(-log(significance_level / 2) * 0.5), digits = 3)))
}

analyse_critical_values(samples_count = 20, significance_level = 0.05)


# Task 2

get_mean_scaled_experiment_pairs <- function () {
  experiments_df <- data.frame(matrix(nrow = 21, ncol = 5))
  colnames(experiments_df) <- c("distribution1", "distribution2", "experiment_name", "type", "parameter_value")
  experiments_df["distribution1"] <- c(rep("norm_mean_0", 7), rep("logis_mean_0", 7), rep("cauchy_mean_0", 7))
  experiments_df["distribution2"] <- c(
    "norm_mean_0_2", "norm_mean_0_4", "norm_mean_0_6", "norm_mean_0_8", "norm_mean_1", "norm_mean_1_2", "norm_mean_1_4",
    "logis_mean_0_2", "logis_mean_0_4", "logis_mean_0_6", "logis_mean_0_8", "logis_mean_1", "logis_mean_1_2",
    "logis_mean_1_4", "cauchy_mean_0", "cauchy_mean_0_5", "cauchy_mean_1", "cauchy_mean_1_5", "cauchy_mean_2",
    "cauchy_mean_2_5", "cauchy_mean_3"
  )
  experiments_df["experiment_name"] <- with(experiments_df, paste(distribution1, distribution2, sep="#"))
  experiments_df["type"] <- c(rep("normal", 7), rep("logistic", 7), rep("cauchy", 7))
  experiments_df["parameter_value"] <- c(
    seq(from = 0.2, to = 1.4, by = 0.2), seq(from = 0.2, to = 1.4, by = 0.2), seq(from = 0.0, to = 3.0, by = 0.5)
  )
  return(experiments_df)
}

estimate_test_power <- function (
  statistic_func, experiments_df, samples_func, h0_rejection_func, samples_count, simulations_count = 10000
) {
  simulated_test_statistics <- simulate_test_statistic_values(
    statistic_func, experiments_df, samples_func, samples_count, simulations_count
  )
  estimated_test_powers <- data.frame(matrix(data = NA, nrow = 1, ncol = dim(simulated_test_statistics)[2]))
  colnames(estimated_test_powers) <- colnames(simulated_test_statistics)
  for (column in colnames(estimated_test_powers)) {
    mean(h0_rejection_func(simulated_test_statistics[[column]]))
    estimated_test_powers[column] <- mean(h0_rejection_func(simulated_test_statistics[[column]]))
  }
  return(estimated_test_powers)
}

generate_tests_powers <- function (
  experiments_df, samples_func, samples_count, simulations_count = 10000, significance_level = 0.05
) {
  left_critical_value <- significance_level / 2
  right_critical_value <- 1 - significance_level / 2
  chisq1_h0_rejection_func <- function (x) !between(pchisq(x, df = 1), left_critical_value, right_critical_value)
  chisq2_h0_rejection_func <- function (x) !between(pchisq(x, df = 2), left_critical_value, right_critical_value)
  ks_h0_rejection_func <- function(x) x >= sqrt(-log(significance_level) * 0.5)
  wilcoxon_test_power_df <- estimate_test_power(
    calculate_wilcoxon_statistic, experiments_df, samples_func, chisq1_h0_rejection_func, samples_count,
    simulations_count
  )
  ab_test_power_df <- estimate_test_power(
    calculate_ab_statistic, experiments_df, samples_func, chisq1_h0_rejection_func, samples_count, simulations_count
  )
  lepage_test_power_df <- estimate_test_power(
    calculate_lepage_statistic, experiments_df, samples_func, chisq2_h0_rejection_func, samples_count, simulations_count
  )
  ks_test_power_df <- estimate_test_power(
    calculate_ks_statistic, experiments_df, samples_func, ks_h0_rejection_func, samples_count, simulations_count
  )
  joint_test_power_df <- data.frame(matrix(nrow = 4, ncol = dim(wilcoxon_test_power_df)[2] + 1))
  colnames(joint_test_power_df) <- c("test_name", colnames(wilcoxon_test_power_df))
  joint_test_power_df["test_name"] <- c("Wilcoxon test", "AB test", "Lepage test", "KS test")
  joint_test_power_df[1,  colnames(wilcoxon_test_power_df)] <- wilcoxon_test_power_df[1, ]
  joint_test_power_df[2,  colnames(ab_test_power_df)] <- ab_test_power_df[1, ]
  joint_test_power_df[3,  colnames(lepage_test_power_df)] <- lepage_test_power_df[1, ]
  joint_test_power_df[4,  colnames(ks_test_power_df)] <- ks_test_power_df[1, ]
  return(joint_test_power_df)
}

extract_points_from_tests_powers_df <- function (tests_powers_df, column_names, y_points) {
  points_df <- data.frame(matrix(nrow = 0, ncol = 3))
  for (column_index in seq_along(column_names)) {
    extracted_df <- tests_powers_df[c("test_name", column_names[column_index])]
    colnames(extracted_df)[colnames(extracted_df) == column_names[column_index]] <- "y"
    extracted_df["x"] <- y_points[column_index]
    points_df <- rbind(points_df, extracted_df)
  }
  return(points_df)
}

plot_tests_powers <- function (tests_powers_df, title, x_label, y_label) {
  print(ggplot(tests_powers_df) +
    geom_point(aes(x = x, y = y, colour = factor(test_name)), shape = 19, size = 2.2) +
    xlab(x_label) +
    ylab(y_label) +
    ggtitle(title) +
    labs(colour = "Test type") +
    theme(plot.title = element_text(hjust = 0.5)))
}

analyse_mean_scaled_tests_powers <- function (samples_count, significance_level = 0.05) {
  mean_scaled_experiments_df <- get_mean_scaled_experiment_pairs()
  tests_powers_df <- generate_tests_powers(
    mean_scaled_experiments_df, create_mean_scaled_distributions_df, samples_count,
    significance_level = significance_level
  )
  normal_powers_df <- mean_scaled_experiments_df[mean_scaled_experiments_df$type == "normal", ]
  normal_distribution_df <- extract_points_from_tests_powers_df(
    tests_powers_df, normal_powers_df$experiment_name, normal_powers_df$parameter_value
  )
  plot_tests_powers(
    normal_distribution_df, title = "A comparison of tests powers for Normal(0, 1) vs. Normal(x, 1)",
    x_label = "Mean parameter of a normal distribution", y_label = "Test power"
  )
  logis_powers_df <- mean_scaled_experiments_df[mean_scaled_experiments_df$type == "logistic", ]
  logis_distribution_df <- extract_points_from_tests_powers_df(
    tests_powers_df, logis_powers_df$experiment_name, logis_powers_df$parameter_value
  )
  plot_tests_powers(
    logis_distribution_df, title = "A comparison of tests powers for Logistic(0, 1) vs. Logistic(x, 1)",
    x_label = "Location parameter of a Logistic distribution", y_label = "Test power"
  )
  cauchy_powers_df <- mean_scaled_experiments_df[mean_scaled_experiments_df$type == "cauchy", ]
  cauchy_distribution_df <- extract_points_from_tests_powers_df(
    tests_powers_df, cauchy_powers_df$experiment_name, cauchy_powers_df$parameter_value
  )
  plot_tests_powers(
    cauchy_distribution_df, title = "A comparison of tests powers for Cauchy(0, 1) vs. Cauchy(x, 1)",
    x_label = "Location parameter of a Cauchy distribution", y_label = "Test power"
  )
}

analyse_mean_scaled_tests_powers(samples_count = 20, significance_level = 0.05)


# Task 3

get_sd_scaled_experiment_pairs <- function () {
  experiments_df <- data.frame(matrix(nrow = 21, ncol = 5))
  colnames(experiments_df) <- c("distribution1", "distribution2", "experiment_name", "type", "paramter_value")
  experiments_df["distribution1"] <- c(rep("norm_sd_1", 7), rep("logis_sd_1", 7), rep("cauchy_sd_1", 7))
  experiments_df["distribution2"] <- c(
    "norm_sd_1", "norm_sd_1_5", "norm_sd_2", "norm_sd_2_5", "norm_sd_3", "norm_sd_3_5", "norm_sd_4", "logis_sd_1",
    "logis_sd_1_5", "logis_sd_2", "logis_sd_2_5", "logis_sd_3", "logis_sd_3_5", "logis_sd_4", "cauchy_sd_1",
    "cauchy_sd_2", "cauchy_sd_3", "cauchy_sd_4", "cauchy_sd_5", "cauchy_sd_6", "cauchy_sd_7"
  )
  experiments_df["experiment_name"] <- with(experiments_df, paste(distribution1, distribution2, sep="#"))
  experiments_df["type"] <- c(rep("normal", 7), rep("logistic", 7), rep("cauchy", 7))
  experiments_df["parameter_value"] <- c(
    seq(from = 1.0, to = 4.0, by = 0.5), seq(from = 1.0, to = 4.0, by = 0.5), seq(from = 1.0, to = 7.0, by = 1.0)
  )
  return(experiments_df)
}

analyse_sd_scaled_tests_powers <- function (samples_count, significance_level = 0.05) {
  sd_scaled_experiments_df <- get_sd_scaled_experiment_pairs()
  tests_powers_df <- generate_tests_powers(
    sd_scaled_experiments_df, create_sd_scaled_distributions_df, samples_count, significance_level = significance_level
  )
  normal_powers_df <- sd_scaled_experiments_df[sd_scaled_experiments_df$type == "normal", ]
  normal_distribution_df <- extract_points_from_tests_powers_df(
    tests_powers_df, normal_powers_df$experiment_name, normal_powers_df$parameter_value
  )
  plot_tests_powers(
    normal_distribution_df, title = "A comparison of tests powers for Normal(0, 1) vs. Normal(0, x)",
    x_label = "Standard deviation parameter of a normal distribution", y_label = "Test power"
  )
  logis_powers_df <- sd_scaled_experiments_df[sd_scaled_experiments_df$type == "logistic", ]
  logis_distribution_df <- extract_points_from_tests_powers_df(
    tests_powers_df, logis_powers_df$experiment_name, logis_powers_df$parameter_value
  )
  plot_tests_powers(
    logis_distribution_df, title = "A comparison of tests powers for Logistic(0, 1) vs. Logistic(0, x)",
    x_label = "Scale parameter of a Logistic distribution", y_label = "Test power"
  )
  cauchy_powers_df <- sd_scaled_experiments_df[sd_scaled_experiments_df$type == "cauchy", ]
  cauchy_distribution_df <- extract_points_from_tests_powers_df(
    tests_powers_df, cauchy_powers_df$experiment_name, cauchy_powers_df$parameter_value
  )
  plot_tests_powers(
    cauchy_distribution_df, title = "A comparison of tests powers for Cauchy(0, 1) vs. Cauchy(1, x)",
    x_label = "Scale parameter of a Cauchy distribution", y_label = "Test power"
  )
}

analyse_sd_scaled_tests_powers(samples_count = 20, significance_level = 0.05)


# Task 4

get_mean_sd_scaled_experiment_pairs <- function () {
  experiments_df <- data.frame(matrix(nrow = 21, ncol = 5))
  colnames(experiments_df) <- c("distribution1", "distribution2", "experiment_name", "type", "parameter_value")
  experiments_df["distribution1"] <- c(
    rep("norm_mean_0_sd_1", 7), rep("logis_mean_0_sd_1", 7), rep("cauchy_mean_0_sd_1", 7)
  )
  experiments_df["distribution2"] <- c(
    "norm_mean_0_2_sd_1", "norm_mean_0_4_sd_1_5", "norm_mean_0_6_sd_2", "norm_mean_0_8_sd_2_5", "norm_mean_1_sd_3",
    "norm_mean_1_2_sd_3_5", "norm_mean_1_4_sd_4", "logis_mean_0_2_sd_1", "logis_mean_0_4_sd_1_5", "logis_mean_0_6_sd_2",
    "logis_mean_0_8_sd_2_5", "logis_mean_1_sd_3", "logis_mean_1_2_sd_3_5", "logis_mean_1_4_sd_4", "cauchy_mean_0_sd_1",
    "cauchy_mean_0_5_sd_2", "cauchy_mean_1_sd_3", "cauchy_mean_1_5_sd_4", "cauchy_mean_2_sd_5", "cauchy_mean_2_5_sd_6",
    "cauchy_mean_3_sd_7"
  )
  experiments_df["experiment_name"] <- with(experiments_df, paste(distribution1, distribution2, sep="#"))
  experiments_df["type"] <- c(rep("normal", 7), rep("logistic", 7), rep("cauchy", 7))
  experiments_df["parameter_value"] <- c(
    seq(from = 0.2, to = 1.4, by = 0.2), seq(from = 0.2, to = 1.4, by = 0.2), seq(from = 0.0, to = 3.0, by = 0.5)
  )
  return(experiments_df)
}

analyse_mean_sd_scaled_tests_powers <- function (samples_count, significance_level = 0.05) {
  mean_sd_scaled_experiments_df <- get_mean_sd_scaled_experiment_pairs()
  tests_powers_df <- generate_tests_powers(
    mean_sd_scaled_experiments_df, create_mean_sd_scaled_distributions_df, samples_count,
    significance_level = significance_level
  )
  normal_powers_df <- mean_sd_scaled_experiments_df[mean_sd_scaled_experiments_df$type == "normal", ]
  normal_distribution_df <- extract_points_from_tests_powers_df(
    tests_powers_df, normal_powers_df$experiment_name, normal_powers_df$parameter_value
  )
  plot_tests_powers(
    normal_distribution_df, title = "A comparison of tests powers for Normal(0, 1) vs. Normal(x, 2.5x + 0.5))",
    x_label = "Mean parameter of a normal distribution",
    y_label = "Test power"
  )
  logis_powers_df <- mean_sd_scaled_experiments_df[mean_sd_scaled_experiments_df$type == "logistic", ]
  logis_distribution_df <- extract_points_from_tests_powers_df(
    tests_powers_df, logis_powers_df$experiment_name, logis_powers_df$parameter_value
  )
  plot_tests_powers(
    logis_distribution_df, title = "A comparison of tests powers for Logistic(0, 1) vs. Logistic(x, 2.5x + 0.5)",
    x_label = "Location parameter of a Logistic distribution",
    y_label = "Test power"
  )
  cauchy_powers_df <- mean_sd_scaled_experiments_df[mean_sd_scaled_experiments_df$type == "cauchy", ]
  cauchy_distribution_df <- extract_points_from_tests_powers_df(
    tests_powers_df, cauchy_powers_df$experiment_name, cauchy_powers_df$parameter_value
  )
  plot_tests_powers(
    cauchy_distribution_df, title = "A comparison of tests powers for Cauchy(0, 1) vs. Cauchy(x, 2x + 1)",
    x_label = "Location parameter of a Cauchy distribution",
    y_label = "Test power"
  )
}

analyse_mean_sd_scaled_tests_powers(samples_count = 20, significance_level = 0.05)


# Task 5

analyse_critical_values(samples_count = 50, significance_level = 0.05)


# Task 6

analyse_mean_scaled_tests_powers(samples_count = 50, significance_level = 0.05)
analyse_sd_scaled_tests_powers(samples_count = 50, significance_level = 0.05)
analyse_mean_sd_scaled_tests_powers(samples_count = 50, significance_level = 0.05)

