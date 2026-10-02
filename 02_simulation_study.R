# Simulation study and exported tables

# ------------------------------------------------------------
# Settings

# Main simulation settings
nsim <- 100
n_schools <- 100
n_per_school <- 50
n_imputations <- 5
n_iterations <- 5
n_workers <- 16

# Reproducible random-number streams
simulation_seed <- 2026
missingness_seed <- 5000
runtime_seed <- 9000
seed_setting_gap <- 1000

age_missing_rate <- 0.20
study_hours_missing_rate <- 0.25
score_missing_rate <- 0.20

runtime_reps <- 20
icc_study_hours_value <- 5

true_beta_study_hours <- 1.8
true_random_intercept_sd <- 4
true_random_slope_sd <- 0.4
true_residual_sd <- 5
true_random_intercept_slope_cov <- 0

output_folder <- file.path(getwd(), "simulation_results")
graph_folder <- file.path(output_folder, "graphs")
table_folder <- file.path(output_folder, "tables")

if (dir.exists(output_folder)) {
  unlink(output_folder, recursive = TRUE, force = TRUE)
}

dir.create(graph_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(table_folder, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# True ICC

true_icc <- function(study_hours_value = icc_study_hours_value) {
  random_variance <-
    true_random_intercept_sd^2 +
    2 * study_hours_value * true_random_intercept_slope_cov +
    study_hours_value^2 * true_random_slope_sd^2

  residual_variance <- true_residual_sd^2
  random_variance / (random_variance + residual_variance)
}

# ------------------------------------------------------------
# Simulate one complete two-level dataset

sim_data <- function(n_schools,
                     n_per_school,
                     beta_study_hours = true_beta_study_hours) {
  n_total <- n_schools * n_per_school

  school_id <- as.integer(
    rep(seq_len(n_schools), each = n_per_school)
  )

  age <- rnorm(
    n_total,
    mean = 16,
    sd = 1
  )

  study_hours <- rnorm(
    n_total,
    mean = 5,
    sd = 2
  )

  school_intercept <- rnorm(
    n_schools,
    mean = 0,
    sd = true_random_intercept_sd
  )

  school_slope <- rnorm(
    n_schools,
    mean = 0,
    sd = true_random_slope_sd
  )

  score <-
    50 +
    2.5 * age +
    beta_study_hours * study_hours +
    school_intercept[school_id] +
    school_slope[school_id] * study_hours +
    rnorm(
      n_total,
      mean = 0,
      sd = true_residual_sd
    )

  data.frame(
    school_id = school_id,
    age = age,
    study_hours = study_hours,
    score = score
  )
}

# ------------------------------------------------------------
# Add MCAR missingness

add_missingness <- function(data, missingness_setting) {
  incomplete <- data
  n_total <- nrow(incomplete)

  # Predictors missing.
  if (missingness_setting == "Predictors missing") {
    age_rows <- sample(
      seq_len(n_total),
      size = round(age_missing_rate * n_total)
    )

    study_rows <- sample(
      seq_len(n_total),
      size = round(study_hours_missing_rate * n_total)
    )

    incomplete[age_rows, "age"] <- NA
    incomplete[study_rows, "study_hours"] <- NA
  }

  # Outcome missing.
  if (missingness_setting == "Outcome missing") {
    score_rows <- sample(
      seq_len(n_total),
      size = round(score_missing_rate * n_total)
    )

    incomplete[score_rows, "score"] <- NA
  }

  incomplete
}

# ------------------------------------------------------------
# Building the mice method vector and predictor matrix

make_mice_setup <- function(data, imputation_method, missingness_setting) {
  method <- c(
    school_id = "",
    age = "",
    study_hours = "",
    score = ""
  )

  if (missingness_setting == "Predictors missing") {
    method["age"] <- imputation_method
    method["study_hours"] <- imputation_method
  }

  if (missingness_setting == "Outcome missing") {
    method["score"] <- imputation_method
  }

  predictor_matrix <- mice::make.predictorMatrix(data)
  predictor_matrix[,] <- 0
  predictor_matrix[, "school_id"] <- -2

  if (missingness_setting == "Predictors missing") {
    predictor_matrix["age", "study_hours"] <- 1
    predictor_matrix["age", "score"] <- 1

    predictor_matrix["study_hours", "age"] <- 1
    predictor_matrix["study_hours", "score"] <- 1
  }

  if (missingness_setting == "Outcome missing") {
    predictor_matrix["score", "age"] <- 1

    # Random slope for study hours in the score model
    predictor_matrix["score", "study_hours"] <- 2
  }

  list(
    method = method,
    predictor_matrix = predictor_matrix
  )
}

# ------------------------------------------------------------
# Fitting the final analysis model

fit_analysis_model <- function(data) {
  fit <- try(
    glmmTMB::glmmTMB(
      score ~ age + study_hours + (1 + study_hours | school_id),
      data = data,
      family = gaussian(),
      REML = FALSE
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    stop("Final glmmTMB analysis model failed to fit.")
  }

  if (!is.null(fit$fit$convergence) && fit$fit$convergence != 0) {
    stop("Final glmmTMB analysis model did not converge.")
  }

  if (!is.null(fit$sdr$pdHess) && !isTRUE(fit$sdr$pdHess)) {
    stop("Final glmmTMB analysis model has a non-positive-definite Hessian.")
  }

  fit
}

# ------------------------------------------------------------
# Extract results from one fitted model

extract_model_results <- function(model,
                                  true_beta = true_beta_study_hours,
                                  study_hours_value = icc_study_hours_value) {
  coefficient_table <- summary(model)$coefficients$cond

  estimate <- coefficient_table["study_hours", "Estimate"]
  standard_error <- coefficient_table["study_hours", "Std. Error"]

  residual_df <- nobs(model) - length(glmmTMB::fixef(model)$cond)

  if (!is.finite(residual_df) || residual_df <= 0) {
    stop("Residual degrees of freedom are not valid for the analysis model.")
  }

  t_critical <- qt(0.975, df = residual_df)
  lower <- estimate - t_critical * standard_error
  upper <- estimate + t_critical * standard_error

  random_matrix <- as.matrix(VarCorr(model)$cond$school_id)

  random_intercept_variance <- random_matrix["(Intercept)", "(Intercept)"]
  random_slope_variance <- random_matrix["study_hours", "study_hours"]
  random_covariance <- random_matrix["(Intercept)", "study_hours"]
  residual_sd <- sigma(model)

  random_variance_at_value <-
    random_intercept_variance +
    2 * study_hours_value * random_covariance +
    study_hours_value^2 * random_slope_variance

  icc <- random_variance_at_value /
    (random_variance_at_value + residual_sd^2)

  c(
    estimate = estimate,
    se = standard_error,
    covered = as.numeric(lower <= true_beta && true_beta <= upper),
    ci_width = upper - lower,
    random_intercept_sd = sqrt(random_intercept_variance),
    random_slope_sd = sqrt(random_slope_variance),
    random_intercept_slope_cov = random_covariance,
    residual_sd = residual_sd,
    icc = icc,
    df = residual_df
  )
}

# ------------------------------------------------------------
# Run one MI analysis and pool the study-hours coefficient

run_mi_analysis <- function(incomplete_data,
                            imputation_method,
                            missingness_setting,
                            m = n_imputations,
                            maxit = n_iterations) {
  setup <- make_mice_setup(
    data = incomplete_data,
    imputation_method = imputation_method,
    missingness_setting = missingness_setting
  )

  if (imputation_method == "2l.lmer") {
    imputed <- mice::mice(
      incomplete_data,
      method = setup$method,
      predictorMatrix = setup$predictor_matrix,
      m = m,
      maxit = maxit,
      printFlag = FALSE,
      REML = FALSE
    )
  } else {
    imputed <- mice::mice(
      incomplete_data,
      method = setup$method,
      predictorMatrix = setup$predictor_matrix,
      m = m,
      maxit = maxit,
      printFlag = FALSE
    )
  }

  completed_data <- mice::complete(imputed, action = "all")
  fitted_models <- lapply(completed_data, fit_analysis_model)

  model_results <- vapply(
    fitted_models,
    extract_model_results,
    FUN.VALUE = c(
      estimate = 0,
      se = 0,
      covered = 0,
      ci_width = 0,
      random_intercept_sd = 0,
      random_slope_sd = 0,
      random_intercept_slope_cov = 0,
      residual_sd = 0,
      icc = 0,
      df = 0
    )
  )

  pooled <- mice::pool.scalar(
    Q = model_results["estimate", ],
    U = model_results["se", ]^2
  )

  pooled_t_critical <- qt(0.975, df = pooled$df)
  pooled_se <- sqrt(pooled$t)
  pooled_lower <- pooled$qbar - pooled_t_critical * pooled_se
  pooled_upper <- pooled$qbar + pooled_t_critical * pooled_se

  c(
    estimate = pooled$qbar,
    se = pooled_se,
    covered = as.numeric(
      pooled_lower <= true_beta_study_hours &&
        true_beta_study_hours <= pooled_upper
    ),
    ci_width = pooled_upper - pooled_lower,
    random_intercept_sd = mean(model_results["random_intercept_sd", ]),
    random_slope_sd = mean(model_results["random_slope_sd", ]),
    random_intercept_slope_cov = mean(model_results["random_intercept_slope_cov", ]),
    residual_sd = mean(model_results["residual_sd", ]),
    icc = mean(model_results["icc", ]),
    df = pooled$df,
    fmi = if (!is.null(pooled$fmi)) pooled$fmi else NA_real_
  )
}

# ------------------------------------------------------------
# We keep the simulation going if one analysis fails.
# This way, one bad fit does not stop the whole simulation,
# and we can still see how often each method works across all repetitions later.

run_safely <- function(code) {
  tryCatch(
    list(
      result = code,
      success = TRUE
    ),
    error = function(e) {
      list(
        result = c(
          estimate = NA_real_,
          se = NA_real_,
          covered = NA_real_,
          ci_width = NA_real_,
          random_intercept_sd = NA_real_,
          random_slope_sd = NA_real_,
          random_intercept_slope_cov = NA_real_,
          residual_sd = NA_real_,
          icc = NA_real_,
          df = NA_real_,
          fmi = NA_real_
        ),
        success = FALSE
      )
    }
  )
}

# ------------------------------------------------------------
# One simulation repetition for one missingness setting

run_one_simulation <- function(sim_id, missingness_setting) {
  setting_offset <- match(
    missingness_setting,
    missingness_settings
  )

    set.seed(simulation_seed + sim_id)

  complete_data <- sim_data(
    n_schools = n_schools,
    n_per_school = n_per_school
  )

  set.seed(missingness_seed + sim_id + setting_offset * seed_setting_gap)

  incomplete_data <- add_missingness(
    data = complete_data,
    missingness_setting = missingness_setting
  )

  complete_fit <- run_safely({
    model <- fit_analysis_model(complete_data)
    base <- extract_model_results(model)
    c(base, fmi = NA_real_)
  })

  complete_case_fit <- run_safely({
    complete_cases <- incomplete_data[complete.cases(incomplete_data), , drop = FALSE]
    model <- fit_analysis_model(complete_cases)
    base <- extract_model_results(model)
    c(base, fmi = NA_real_)
  })

  glmmtmb_fit <- run_safely({
    run_mi_analysis(
      incomplete_data = incomplete_data,
      imputation_method = "2l.glmmTMB",
      missingness_setting = missingness_setting
    )
  })

  lmer_fit <- run_safely({
    run_mi_analysis(
      incomplete_data = incomplete_data,
      imputation_method = "2l.lmer",
      missingness_setting = missingness_setting
    )
  })

  method_results <- list(
    "Complete data" = complete_fit,
    "Complete-case analysis" = complete_case_fit,
    "glmmTMB" = glmmtmb_fit,
    "lmer" = lmer_fit
  )

  output_rows <- lapply(names(method_results), function(method_name) {
    method_result <- method_results[[method_name]]

    row <- as.data.frame(as.list(method_result$result), check.names = FALSE)
    row$sim_id <- sim_id
    row$missingness_setting <- missingness_setting
    row$method <- method_name
    row$success <- method_result$success
    row
  })

  do.call(rbind, output_rows)
}

# ------------------------------------------------------------
# Parallel simulation

missingness_settings <- c(
  "Predictors missing",
  "Outcome missing"
)

simulation_grid <- expand.grid(
  sim_id = seq_len(nsim),
  missingness_setting = missingness_settings,
  stringsAsFactors = FALSE
)

cluster <- parallel::makeCluster(n_workers)

parallel::clusterEvalQ(cluster, {
  library(glmmTMB)
  library(mice)
  library(lme4)
  library(MASS)
})

parallel::clusterExport(
  cluster,
  varlist = c(
    "simulation_grid",
    "missingness_settings",
    "nsim",
    "n_schools",
    "n_per_school",
    "n_imputations",
    "n_iterations",
    "simulation_seed",
    "missingness_seed",
    "seed_setting_gap",
    "age_missing_rate",
    "study_hours_missing_rate",
    "score_missing_rate",
    "icc_study_hours_value",
    "true_beta_study_hours",
    "true_random_intercept_sd",
    "true_random_slope_sd",
    "true_residual_sd",
    "true_random_intercept_slope_cov",
    "mice.impute.2l.glmmTMB",
    "sim_data",
    "add_missingness",
    "make_mice_setup",
    "fit_analysis_model",
    "extract_model_results",
    "run_mi_analysis",
    "run_safely",
    "run_one_simulation"
  ),
  envir = environment()
)

simulation_list <- pblapply(
  seq_len(nrow(simulation_grid)),
  function(row_id) {
    run_one_simulation(
      sim_id = simulation_grid$sim_id[row_id],
      missingness_setting = simulation_grid$missingness_setting[row_id]
    )
  },
  cl = cluster
)

parallel::stopCluster(cluster)

simulation_results <- do.call(rbind, simulation_list)
rownames(simulation_results) <- NULL

simulation_results$bias <- simulation_results$estimate - true_beta_study_hours
simulation_results$percent_bias <- 100 * simulation_results$bias / true_beta_study_hours

saveRDS(
  simulation_results,
  file = file.path(output_folder, "simulation_raw_results.rds")
)

write.csv(
  simulation_results,
  file = file.path(output_folder, "simulation_raw_results.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Summaries

successful_results <- subset(simulation_results, success & is.finite(estimate))

split_groups <- split(
  successful_results,
  interaction(
    successful_results$missingness_setting,
    successful_results$method,
    drop = TRUE
  )
)

performance_summary <- do.call(
  rbind,
  lapply(split_groups, function(group_data) {
    empirical_se <- sd(group_data$estimate, na.rm = TRUE)
    mean_model_se <- mean(group_data$se, na.rm = TRUE)

    data.frame(
      missingness_setting = group_data$missingness_setting[1],
      method = group_data$method[1],
      n_success = nrow(group_data),
      estimate = mean(group_data$estimate, na.rm = TRUE),
      bias_percent = mean(group_data$percent_bias, na.rm = TRUE),
      empirical_se = empirical_se,
      model_se = mean_model_se,
      se_ratio = mean_model_se / empirical_se,
      rmse = sqrt(
        mean(
          (group_data$estimate - true_beta_study_hours)^2,
          na.rm = TRUE
        )
      ),
      coverage_percent = 100 * mean(group_data$covered, na.rm = TRUE),
      ci_width = mean(group_data$ci_width, na.rm = TRUE),
      fmi = mean(group_data$fmi, na.rm = TRUE)
    )
  })
)

variance_summary <- do.call(
  rbind,
  lapply(split_groups, function(group_data) {
    data.frame(
      missingness_setting = group_data$missingness_setting[1],
      method = group_data$method[1],
      n_success = nrow(group_data),
      random_intercept_sd = mean(group_data$random_intercept_sd, na.rm = TRUE),
      random_slope_sd = mean(group_data$random_slope_sd, na.rm = TRUE),
      residual_sd = mean(group_data$residual_sd, na.rm = TRUE),
      icc_at_study_hours_5 = mean(group_data$icc, na.rm = TRUE)
    )
  })
)

success_groups <- split(
  simulation_results,
  interaction(
    simulation_results$missingness_setting,
    simulation_results$method,
    drop = TRUE
  )
)

success_summary <- do.call(
  rbind,
  lapply(success_groups, function(group_data) {
    n_attempted <- nrow(group_data)
    n_success <- sum(group_data$success, na.rm = TRUE)

    data.frame(
      missingness_setting = group_data$missingness_setting[1],
      method = group_data$method[1],
      n_success = n_success,
      n_attempted = n_attempted,
      n_failed = n_attempted - n_success,
      success_percent = 100 * n_success / n_attempted
    )
  })
)

round_numeric_columns <- function(data, digits = 3) {
  numeric_columns <- vapply(data, is.numeric, logical(1))
  data[numeric_columns] <- lapply(data[numeric_columns], round, digits = digits)
  data
}

performance_summary_report <- round_numeric_columns(performance_summary, 3)
variance_summary_report <- round_numeric_columns(variance_summary, 3)
success_summary_report <- round_numeric_columns(success_summary, 3)

write.csv(
  performance_summary_report,
  file = file.path(table_folder, "table_1_statistical_performance.csv"),
  row.names = FALSE
)

write.csv(
  variance_summary_report,
  file = file.path(table_folder, "table_2_variance_components_and_icc.csv"),
  row.names = FALSE
)

write.csv(
  success_summary_report,
  file = file.path(table_folder, "table_4_success_rates.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# Runtime benchmark


runtime_rows <- list()
runtime_row_id <- 1

for (runtime_id in seq_len(runtime_reps)) {
  for (missingness_setting in missingness_settings) {

    setting_offset <- match(
      missingness_setting,
      missingness_settings
    )

    set.seed(
      runtime_seed +
        runtime_id +
        setting_offset * seed_setting_gap
    )

    runtime_complete <- sim_data(
      n_schools = n_schools,
      n_per_school = n_per_school
    )

    runtime_incomplete <- add_missingness(
      data = runtime_complete,
      missingness_setting = missingness_setting
    )

    method_order_runtime <- if (
      runtime_id %% 2 == 1
    ) {
      c(
        "glmmTMB",
        "lmer"
      )
    } else {
      c(
        "lmer",
        "glmmTMB"
      )
    }

    for (
      method_name in method_order_runtime
    ) {

      imputation_method <- if (
        method_name == "glmmTMB"
      ) {
        "2l.glmmTMB"
      } else {
        "2l.lmer"
      }

      elapsed <- system.time({
        runtime_result <- run_safely({
          run_mi_analysis(
            incomplete_data =
              runtime_incomplete,
            imputation_method =
              imputation_method,
            missingness_setting =
              missingness_setting
          )
        })
      })[["elapsed"]]

      runtime_rows[[runtime_row_id]] <-
        data.frame(
          runtime_id = runtime_id,
          missingness_setting =
            missingness_setting,
          method = method_name,
          elapsed_seconds = elapsed,
          success =
            runtime_result$success,
          stringsAsFactors = FALSE
        )

      runtime_row_id <-
        runtime_row_id + 1
    }
  }
}

runtime_results <- do.call(
  rbind,
  runtime_rows
)

runtime_successful <- subset(
  runtime_results,
  success
)

runtime_groups <- split(
  runtime_successful,
  interaction(
    runtime_successful$missingness_setting,
    runtime_successful$method,
    drop = TRUE
  )
)

runtime_summary <- do.call(
  rbind,
  lapply(
    runtime_groups,
    function(group_data) {
      data.frame(
        missingness_setting =
          group_data$missingness_setting[1],
        method =
          group_data$method[1],
        n_success =
          nrow(group_data),
        mean_seconds =
          mean(
            group_data$elapsed_seconds
          ),
        median_seconds =
          median(
            group_data$elapsed_seconds
          ),
        sd_seconds =
          sd(
            group_data$elapsed_seconds
          ),
        iqr_seconds =
          IQR(
            group_data$elapsed_seconds
          ),
        minimum_seconds =
          min(
            group_data$elapsed_seconds
          ),
        maximum_seconds =
          max(
            group_data$elapsed_seconds
          )
      )
    }
  )
)

runtime_summary_report <-
  round_numeric_columns(
    runtime_summary,
    3
  )

write.csv(
  runtime_results,
  file = file.path(
    output_folder,
    "serial_runtime_raw.csv"
  ),
  row.names = FALSE
)

write.csv(
  runtime_summary_report,
  file = file.path(
    table_folder,
    "table_3_serial_runtime_benchmark.csv"
  ),
  row.names = FALSE
)

cat("\nSimulation and runtime benchmark finished.\n")