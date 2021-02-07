# Statistics - List 6

library("ggplot2")
library("dplyr")
library("pracma")
library("latex2exp")
set.seed(100)


calculate_pearson_statistic <- function (samples, sets_count) {
  set_indexes <- pmin(floor(sets_count * samples), sets_count - 1)
  set_cardinalities <- table(set_indexes)
  cardinalities <- unname(set_cardinalities[order(names(set_cardinalities))])
  expected_cardinality <- length(samples) / sets_count
  return(sum((cardinalities - expected_cardinality) ^ 2 / expected_cardinality))
}


calculate_neyman_statistic <- function (samples, polynomials_count) {
  statistic_value <- 0.0
  for (i in 1:polynomials_count) {
    legendre_values_sum <- sqrt(2 * i + 1) * sum(legendre(i, 2 * samples - 1)[1, ])
    statistic_value <- statistic_value + (1 / sqrt(length(samples)) * legendre_values_sum) ^ 2
  }
  return(statistic_value)
}


calculate_ks_statistic <- function (samples) {
  empirical_cdf <- 0:length(samples) / length(samples)
  max_diff <- 0.0
  sorted_samples <- sort(samples)
  for (i in seq_along(sorted_samples)) {
    max_diff <- max(max_diff, sorted_samples[i] - empirical_cdf[i])
    max_diff <- max(max_diff, sorted_samples[i] - empirical_cdf[i + 1])
  }
  return(sqrt(length(samples)) * max_diff)
}


# Task 1

generate_uniform_samples <- function (samples_count) {
  samples_df <- data.frame(matrix(nrow = samples_count, ncol = 1))
  colnames(samples_df) <- c("uniform_0_1")
  samples_df["uniform_0_1"] <- runif(samples_count)
  return(samples_df)
}


simulate_test_statistic_values <- function (statistic_func, samples_func, samples_count, simulations_count = 10000) {
  samples_df <- samples_func(simulations_count * samples_count)
  statistics_df <- data.frame(matrix(data = NA, nrow = simulations_count, ncol = dim(samples_df)[2]))
  colnames(statistics_df) <- colnames(samples_df)
  for (i in 1:simulations_count) {
    min_sample_index <- samples_count * (i - 1) + 1
    max_sample_index <- samples_count * i
    for (experiment_name in colnames(samples_df)) {
      samples_slice <- samples_df[[experiment_name]][min_sample_index:max_sample_index]
      statistics_df[i, experiment_name] <- statistic_func(samples_slice)
    }
  }
  return(statistics_df)
}


estimate_test_statistic_critical_value <- function (
  statistic_func, samples_func, samples_count, simulations_count = 10000, significance_level = 0.05
) {
  simulated_test_statistics <- simulate_test_statistic_values(
    statistic_func, samples_func, samples_count, simulations_count
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


analyse_critical_values <- function (samples_count, simulations_count = 10000, significance_level = 0.05) {
  critical_value_func <- function (statistic_func) estimate_test_statistic_critical_value(
    statistic_func, generate_uniform_samples, samples_count, simulations_count, significance_level
  )
  pearson4_statistic_func <- function(samples) calculate_pearson_statistic(samples, sets_count = 4)
  estimated_pearson4_critical_value <- critical_value_func(pearson4_statistic_func)
  pearson8_statistic_func <- function(samples) calculate_pearson_statistic(samples, sets_count = 8)
  estimated_pearson8_critical_value <- critical_value_func(pearson8_statistic_func)
  neyman1_statistic_func <- function(samples) calculate_neyman_statistic(samples, polynomials_count = 1)
  estimated_neyman1_critical_value <- critical_value_func(neyman1_statistic_func)
  neyman4_statistic_func <- function(samples) calculate_neyman_statistic(samples, polynomials_count = 4)
  estimated_neyman4_critical_value <- critical_value_func(neyman4_statistic_func)
  neyman8_statistic_func <- function(samples) calculate_neyman_statistic(samples, polynomials_count = 8)
  estimated_neyman8_critical_value <- critical_value_func(neyman8_statistic_func)
  estimated_ks_critical_value <- critical_value_func(calculate_ks_statistic)
  print(paste0("Estimations for samples_count = ", samples_count))
  print(paste0("Estimated P_4 test critical value: ", round(estimated_pearson4_critical_value, digits = 3)))
  print(paste0("Real P_4 test critical value: ", round(qchisq(1 - significance_level, df = 3), digits = 3)))
  print(paste0("Estimated P_8 test critical value: ", round(estimated_pearson8_critical_value, digits = 3)))
  print(paste0("Real P_8 test critical value: ", round(qchisq(1 - significance_level, df = 7), digits = 3)))
  print(paste0("Estimated N_1 test critical value: ", round(estimated_neyman1_critical_value, digits = 3)))
  print(paste0("Real N_1 test critical value: ", round(qchisq(1 - significance_level, df = 1), digits = 3)))
  print(paste0("Estimated N_4 test critical value: ", round(estimated_neyman4_critical_value, digits = 3)))
  print(paste0("Real N_4 test critical value: ", round(qchisq(1 - significance_level, df = 4), digits = 3)))
  print(paste0("Estimated N_8 test critical value: ", round(estimated_neyman8_critical_value, digits = 3)))
  print(paste0("Real N_8 test critical value: ", round(qchisq(1 - significance_level, df = 8), digits = 3)))
  print(paste0("Estimated KS test critical value: ", round(estimated_ks_critical_value, digits = 3)))
  print(paste0("Real KS test critical value: ", round(sqrt(-log(significance_level) * 0.5), digits = 3)))
}


analyse_critical_values(samples_count = 10)
analyse_critical_values(samples_count = 20)
analyse_critical_values(samples_count = 30)
analyse_critical_values(samples_count = 50)
analyse_critical_values(samples_count = 75)
analyse_critical_values(samples_count = 100)


# Task 3

generate_samples_from_pdf <- function (
  pdf_func, max_pdf_value, lower_support = 0, upper_support = 1, samples_count = 1
) {
  generated_samples <- rep(NA, samples_count)
  for (i in 1:samples_count) {
    while (TRUE) {
      x_coord <- runif(n = 1, min = lower_support, max = upper_support)
      y_coord <- max_pdf_value * runif(n = 1, min = 0, max = 1)
      if (y_coord <= pdf_func(x_coord)) {
        break
      }
    }
    generated_samples[[i]] <- x_coord
  }
  return(generated_samples)
}


generate_simple_cos_samples <- function (samples_count) {
  samples_df <- data.frame(matrix(nrow = samples_count, ncol = 1))
  colnames(samples_df) <- c("cos_offset_0_4_cycles_1")
  cos_pdf_func <- function (x) 1 + 0.4 * cos(pi * x)
  samples_df["cos_offset_0_4_cycles_1"] <- generate_samples_from_pdf(
    cos_pdf_func, max_pdf_value = 1.4, samples_count = samples_count
  )
  return(samples_df)
}


estimate_test_power <- function (
  statistic_func, samples_func, h0_rejection_func, samples_count, simulations_count = 10000
) {
  simulated_test_statistics <- simulate_test_statistic_values(
    statistic_func, samples_func, samples_count, simulations_count
  )
  estimated_test_powers <- data.frame(matrix(data = NA, nrow = 1, ncol = dim(simulated_test_statistics)[2]))
  colnames(estimated_test_powers) <- colnames(simulated_test_statistics)
  for (column in colnames(estimated_test_powers)) {
    estimated_test_powers[column] <- mean(h0_rejection_func(simulated_test_statistics[[column]]))
  }
  return(estimated_test_powers)
}


generate_tests_powers <- function (samples_func, samples_count, simulations_count = 10000, significance_level = 0.05) {
  left_critical_value <- significance_level / 2
  right_critical_value <- 1 - significance_level / 2
  chisq_df_rejection_func <- function (df, x) !between(pchisq(x, df = df), left_critical_value, right_critical_value)
  chisq1_rejection_func <- function (x) chisq_df_rejection_func(df = 1, x)
  chisq3_rejection_func <- function (x) chisq_df_rejection_func(df = 3, x)
  chisq4_rejection_func <- function (x) chisq_df_rejection_func(df = 4, x)
  chisq7_rejection_func <- function (x) chisq_df_rejection_func(df = 7, x)
  chisq8_rejection_func <- function (x) chisq_df_rejection_func(df = 8, x)
  ks_h0_rejection_func <- function (x) (
    x < sqrt(-log(right_critical_value) * 0.5) | x >  sqrt(-log(left_critical_value) * 0.5)
  )
  pearson4_statistic_func <- function(samples) calculate_pearson_statistic(samples, sets_count = 4)
  pearson8_statistic_func <- function(samples) calculate_pearson_statistic(samples, sets_count = 8)
  neyman1_statistic_func <- function(samples) calculate_neyman_statistic(samples, polynomials_count = 1)
  neyman4_statistic_func <- function(samples) calculate_neyman_statistic(samples, polynomials_count = 4)
  neyman8_statistic_func <- function(samples) calculate_neyman_statistic(samples, polynomials_count = 8)
  estimate_power_func <- function (statistic_func, h0_rejection_func) estimate_test_power(
    statistic_func, samples_func, h0_rejection_func, samples_count, simulations_count
  )
  pearson4_df <- estimate_power_func(pearson4_statistic_func, chisq3_rejection_func)
  pearson8_df <- estimate_power_func(pearson8_statistic_func, chisq7_rejection_func)
  neyman1_df <- estimate_power_func(neyman1_statistic_func, chisq1_rejection_func)
  neyman4_df <- estimate_power_func(neyman4_statistic_func, chisq4_rejection_func)
  neyman8_df <- estimate_power_func(neyman8_statistic_func, chisq8_rejection_func)
  ks_df <- estimate_power_func(calculate_ks_statistic, ks_h0_rejection_func)
  joint_test_power_df <- data.frame(matrix(nrow = 6, ncol = dim(pearson4_df)[2] + 1))
  colnames(joint_test_power_df) <- c("test_name", colnames(pearson4_df))
  joint_test_power_df["test_name"] <- c("P_4 test", "P_8 test", "N_1 test", "N_4 test", "N_8 test", "KS test")
  joint_test_power_df[1,  colnames(pearson4_df)] <- pearson4_df[1, ]
  joint_test_power_df[2,  colnames(pearson8_df)] <- pearson8_df[1, ]
  joint_test_power_df[3,  colnames(neyman1_df)] <- neyman1_df[1, ]
  joint_test_power_df[4,  colnames(neyman4_df)] <- neyman4_df[1, ]
  joint_test_power_df[5,  colnames(neyman8_df)] <- neyman8_df[1, ]
  joint_test_power_df[6,  colnames(ks_df)] <- ks_df[1, ]
  return(joint_test_power_df)
}


generate_multi_samples_tests_powers <- function (
  samples_func, samples_counts, simulations_count = 10000, significance_level = 0.05
) {
  joint_test_power_df <- data.frame()
  for (samples_count in samples_counts) {
    test_power_df <- generate_tests_powers(samples_func, samples_count, simulations_count, significance_level)
    test_power_df["samples_count"] <- samples_count
    joint_test_power_df <- rbind(joint_test_power_df, test_power_df)
  }
  return(joint_test_power_df)
}


plot_tests_powers <- function (test_power_df, power_column, title, x_label, y_label) {
  print(ggplot(test_power_df, aes_string(y = power_column)) +
    geom_point(aes(x = samples_count, colour = factor(x = test_name)), shape = 19, size = 2.2) +
    scale_shape_identity() +
    xlab(x_label) +
    ylab(y_label) +
    ggtitle(title) +
    labs(colour = "Test type") +
    theme(plot.title = element_text(hjust = 0.5)))
}


analyse_simple_cos_test_power <- function () {
  test_power_df <- generate_multi_samples_tests_powers(
    generate_simple_cos_samples, samples_counts = c(10, 20, 30, 50, 75, 100)
  )
  plot_tests_powers(
    test_power_df, power_column = "cos_offset_0_4_cycles_1",
    title = TeX("Test powers for a distribution with PDF: $1 + 0.4 * cos(pi * x)$"),
    x_label = "Number of samples", y_label = "Test power"
  )
}


analyse_simple_cos_test_power()


# Task 4

generate_custom_cos_samples <- function (samples_count) {
  samples_df <- data.frame(matrix(nrow = samples_count, ncol = 5))
  colnames(samples_df) <- c(
    "cos_offset_0_5_cycles_2", "cos_offset_0_5_cycles_3", "cos_offset_0_6_cycles_4", "cos_offset_0_7_cycles_5",
    "cos_offset_0_7_cycles_6"
  )
  cos_pdf_func <- function (offset, cycles, x) 1 + offset * cos(cycles * pi * x)
  samples_df["cos_offset_0_5_cycles_2"] <- generate_samples_from_pdf(
    function (x) cos_pdf_func(offset = 0.5, cycles = 2, x), max_pdf_value = 1.5, samples_count = samples_count
  )
  samples_df["cos_offset_0_5_cycles_3"] <- generate_samples_from_pdf(
    function (x) cos_pdf_func(offset = 0.5, cycles = 3, x), max_pdf_value = 1.5, samples_count = samples_count
  )
  samples_df["cos_offset_0_6_cycles_4"] <- generate_samples_from_pdf(
    function (x) cos_pdf_func(offset = 0.6, cycles = 4, x), max_pdf_value = 1.6, samples_count = samples_count
  )
  samples_df["cos_offset_0_7_cycles_5"] <- generate_samples_from_pdf(
    function (x) cos_pdf_func(offset = 0.7, cycles = 5, x), max_pdf_value = 1.7, samples_count = samples_count
  )
  samples_df["cos_offset_0_7_cycles_6"] <- generate_samples_from_pdf(
    function (x) cos_pdf_func(offset = 0.7, cycles = 6, x), max_pdf_value = 1.7, samples_count = samples_count
  )
  return(samples_df)
}


analyse_custom_cos_test_powers <- function () {
  test_power_df <- generate_multi_samples_tests_powers(
    generate_custom_cos_samples, samples_counts = c(10, 20, 30, 50, 75, 100)
  )
  plot_tests_powers(
    test_power_df, power_column = "cos_offset_0_5_cycles_2",
    title = TeX("Test powers for a distribution with PDF: $1 + 0.5 * cos(2 * pi * x)$"),
    x_label = "Number of samples", y_label = "Test power"
  )
  plot_tests_powers(
    test_power_df, power_column = "cos_offset_0_5_cycles_3",
    title = TeX("Test powers for a distribution with PDF: $1 + 0.5 * cos(3 * pi * x)$"),
    x_label = "Number of samples", y_label = "Test power"
  )
  plot_tests_powers(
    test_power_df, power_column = "cos_offset_0_6_cycles_4",
    title = TeX("Test powers for a distribution with PDF: $1 + 0.6 * cos(4 * pi * x)$"),
    x_label = "Number of samples", y_label = "Test power"
  )
  plot_tests_powers(
    test_power_df, power_column = "cos_offset_0_7_cycles_5",
    title = TeX("Test powers for a distribution with PDF: $1 + 0.7 * cos(5 * pi * x)$"),
    x_label = "Number of samples", y_label = "Test power"
  )
  plot_tests_powers(
    test_power_df, power_column = "cos_offset_0_7_cycles_6",
    title = TeX("Test powers for a distribution with PDF: $1 + 0.7 * cos(6 * pi * x)$"),
    x_label = "Number of samples", y_label = "Test power"
  )
}


analyse_custom_cos_test_powers()
