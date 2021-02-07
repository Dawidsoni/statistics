# Task 1

count_estimators_frame <- function(generate_samples_func, estimators_count=10000) {
  estimators_frame <- data.frame(matrix(nrow = estimators_count, ncol = 4))
  colnames(estimators_frame) <- c("mean", "median", "weighted_mean", "normal_mean")
  for (i in 1:estimators_count) {
    generated_samples <- generate_samples_func()
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


count_estimators_statistics <- function(estimators_frame, estimated_paramter, excluded_col_names=c()) {
  statistics_frame <- data.frame(matrix(nrow = 1, ncol = 0))
  for (col_name in colnames(estimators_frame)) {
    if (col_name %in% excluded_col_names) {
      next
    }
    statistics_frame[paste0(col_name, '_variance')] <- var(c(estimators_frame[[col_name]]), na.rm = TRUE)
    statistics_frame[paste0(col_name, '_mse')] <- mean(
      (c(estimators_frame[[col_name]]) - estimated_paramter) ^ 2, na.rm = TRUE
    )
    statistics_frame[paste0(col_name, '_bias')] <- mean(
      c(estimators_frame[[col_name]]) - estimated_paramter, na.rm = TRUE
    )
  }
  return(statistics_frame)
}


count_norm_estimators_statistics <- function (n, mean, sd) {
  estimators_frame <- count_estimators_frame(function() rnorm(n=n, mean=mean, sd=sd))
  return(count_estimators_statistics(estimators_frame, estimated_paramter=mean))
}

count_norm_estimators_statistics(n=50, mean=1, sd=1)
count_norm_estimators_statistics(n=50, mean=4, sd=1)
count_norm_estimators_statistics(n=50, mean=1, sd=2)


# Task 5

find_root_using_newton_method <- function(
  data_points, l_function, ld_function, initial_estimate, max_steps = 1000, epsilon = 1e-2
) {
  best_estimate <- initial_estimate
  for (step in 1:max_steps) {
    last_estimate <- best_estimate
    best_estimate <- best_estimate - l_function(data_points, best_estimate) / ld_function(data_points, best_estimate)
    if (is.infinite(best_estimate) || is.na(best_estimate) || is.nan(best_estimate)) {
      return(c(NaN, NaN))
    }
    if (abs(last_estimate - best_estimate) < epsilon) {
      break
    }
  }
  return(c(best_estimate, step))
}


get_mle_estimators_frame <- function (
  generate_samples_func, first_derrivative_func, second_derrivative_func, initial_estimate=NULL, samples_count=10000
) {
  statistics_frame <- data.frame(matrix(nrow = samples_count, ncol = 2))
  colnames(statistics_frame) <- c("mle", "steps")
  is_mean_initial_estimate <- is.null(initial_estimate)
  for (i in 1:samples_count) {
    generated_samples <- generate_samples_func()
    if (is_mean_initial_estimate) {
      initial_estimate <- mean(generated_samples)
    }
    root_info <- find_root_using_newton_method(
      generated_samples, first_derrivative_func, second_derrivative_func, initial_estimate
    )
    statistics_frame[i, "mle"] <- root_info[1]
    statistics_frame[i, "steps"] <- root_info[2]
  }
  return(statistics_frame)
}


count_logistic_mle_estimator_statistics <- function (n, location, scale, initial_estimate=NULL) {
  generate_samples_func <- function() rlogis(n, location, scale)
  lexpr <- function(x, theta) exp((theta - x) / scale)
  dexpr <- function(x, theta) (sqrt(scale) / scale * lexpr(x, theta)) / (sqrt(scale) + sqrt(scale) * lexpr(x, theta))
  first_derrivative <- function (x, theta) length(x) / scale - 2 * sum(dexpr(x, theta))
  ddexpr <- function (x, theta) (1 / scale * lexpr(x, theta)) / ((sqrt(scale) + sqrt(scale) * lexpr(x, theta)) ^ 2)
  second_derivative <- function(x, theta) -2 * sum(ddexpr(x, theta))
  mle_frame <- get_mle_estimators_frame(
    generate_samples_func, first_derrivative, second_derivative, initial_estimate
  )
  statistics_frame <- count_estimators_statistics(mle_frame, location, excluded_col_names = c("steps"))
  statistics_frame["steps_mean"] <- mean(c(mle_frame[["steps"]]), na.rm = TRUE)
  statistics_frame["convergence_rate"] <- 1 - mean(is.nan(c(mle_frame[["steps"]])))
  return(statistics_frame)
}


count_logistic_mle_estimator_statistics(n=50, location=1, scale=1)
count_logistic_mle_estimator_statistics(n=50, location=4, scale=1)
count_logistic_mle_estimator_statistics(n=50, location=1, scale=2)

count_logistic_mle_estimator_statistics(n=50, location=1, scale=1, initial_estimate=2.0)
count_logistic_mle_estimator_statistics(n=50, location=1, scale=1, initial_estimate=10.0)
count_logistic_mle_estimator_statistics(n=50, location=1, scale=1, initial_estimate=-5.0)
count_logistic_mle_estimator_statistics(n=50, location=1, scale=1, initial_estimate=-10.0)

# Task 6

count_cauchy_mle_estimator_statistics <- function (n, location, scale, initial_estimate=NULL) {
  generate_samples_func <- function() rcauchy(n, location, scale)
  lexpr <- function(x, theta) ((x - theta) / scale) ^ 2
  first_derivative <- function (x, theta) sum((2 / scale ^ 2 * (x - theta)) / (1 + lexpr(x, theta)))
  second_derivative <- function(x, theta) sum((2 / scale ^ 2 * (lexpr(x, theta) - 1)) / ((1 + lexpr(x, theta)) ^ 2))
  mle_frame <- get_mle_estimators_frame(
    generate_samples_func, first_derivative, second_derivative, initial_estimate
  )
  statistics_frame <- count_estimators_statistics(mle_frame, location, excluded_col_names = c("steps"))
  statistics_frame["steps_mean"] <- mean(c(mle_frame[["steps"]]), na.rm = TRUE)
  statistics_frame["convergence_rate"] <- 1 - mean(is.nan(c(mle_frame[["steps"]])))
  return(statistics_frame)
}

count_cauchy_mle_estimator_statistics(n=50, location=1, scale=1, initial_estimate=1.5)
count_cauchy_mle_estimator_statistics(n=50, location=4, scale=1, initial_estimate=4.5)
count_cauchy_mle_estimator_statistics(n=50, location=1, scale=2, initial_estimate=1.5)

count_cauchy_mle_estimator_statistics(n=50, location=1, scale=1, initial_estimate=1.75)
count_cauchy_mle_estimator_statistics(n=50, location=1, scale=1, initial_estimate=0.25)
count_cauchy_mle_estimator_statistics(n=50, location=1, scale=1, initial_estimate=3.0)
count_cauchy_mle_estimator_statistics(n=50, location=1, scale=1, initial_estimate=-1.0)

# Task 7

count_norm_estimators_statistics(n=20, mean=1, sd=1)
count_norm_estimators_statistics(n=20, mean=4, sd=1)
count_norm_estimators_statistics(n=20, mean=1, sd=2)

count_logistic_mle_estimator_statistics(n=20, location=1, scale=1)
count_logistic_mle_estimator_statistics(n=20, location=4, scale=1)
count_logistic_mle_estimator_statistics(n=20, location=1, scale=2)

count_cauchy_mle_estimator_statistics(n=20, location=1, scale=1, initial_estimate=1.5)
count_cauchy_mle_estimator_statistics(n=20, location=4, scale=1, initial_estimate=4.5)
count_cauchy_mle_estimator_statistics(n=20, location=1, scale=2, initial_estimate=1.5)

count_norm_estimators_statistics(n=100, mean=1, sd=1)
count_norm_estimators_statistics(n=100, mean=4, sd=1)
count_norm_estimators_statistics(n=100, mean=1, sd=2)

count_logistic_mle_estimator_statistics(n=100, location=1, scale=1)
count_logistic_mle_estimator_statistics(n=100, location=4, scale=1)
count_logistic_mle_estimator_statistics(n=100, location=1, scale=2)

count_cauchy_mle_estimator_statistics(n=100, location=1, scale=1, initial_estimate=1.5)
count_cauchy_mle_estimator_statistics(n=100, location=4, scale=1, initial_estimate=4.5)
count_cauchy_mle_estimator_statistics(n=100, location=1, scale=2, initial_estimate=1.5)
