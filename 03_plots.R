# ------------------------------------------------------------
# Plot helpers

method_order <- c(
  "Complete data",
  "Complete-case analysis",
  "glmmTMB",
  "lmer"
)

method_labels <- c(
  "Complete data" = "Complete data",
  "Complete-case analysis" = "Complete-case",
  "glmmTMB" = "MI: 2l.glmmTMB",
  "lmer" = "MI: 2l.lmer"
)

setting_labels <- c(
  "Predictors missing" = "Age + study hours missing",
  "Outcome missing" = "Score missing"
)

plot_data <- successful_results

plot_data$method <- factor(
  plot_data$method,
  levels = method_order
)
plot_data$missingness_setting <- factor(
  plot_data$missingness_setting,
  levels = missingness_settings
)

report_theme <- function() {
  theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 13,
        hjust = 0,
        margin = margin(b = 8)
      ),
      strip.background = element_rect(
        fill = "grey94",
        colour = "grey70",
        linewidth = 0.35
      ),
      strip.text = element_text(
        face = "bold",
        size = 10.5
      ),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(size = 8.5),
      axis.text.y = element_text(size = 8.5),
      axis.line = element_line(linewidth = 0.4),
      axis.ticks = element_line(linewidth = 0.35),
      legend.title = element_blank(),
      legend.position = "top",
      plot.margin = margin(8, 10, 8, 8)
    )
}

make_pretty_axis <- function(values,
                             include = NULL,
                             n_breaks = 8,
                             lower_bound = -Inf,
                             upper_bound = Inf) {
  values <- c(values[is.finite(values)], include)
  values <- values[is.finite(values)]

  if (length(values) == 0) {
    return(list(
      breaks = waiver(),
      limits = NULL
    ))
  }

  value_range <- range(values)

  if (diff(value_range) == 0) {
    padding <- if (value_range[1] == 0) {
      1
    } else {
      abs(value_range[1]) * 0.10
    }
    value_range <- value_range + c(-padding, padding)
  }

  breaks <- pretty(
    value_range,
    n = n_breaks
  )

  breaks <- breaks[
    breaks >= lower_bound &
      breaks <= upper_bound
  ]

  if (length(breaks) < 2) {
    breaks <- pretty(
      value_range,
      n = n_breaks
    )
  }

  limits <- c(
    max(lower_bound, min(breaks)),
    min(upper_bound, max(breaks))
  )

  list(
    breaks = breaks,
    limits = limits
  )
}

mc_mean_interval <- function(data, value_name) {
  grouped <- split(
    data,
    interaction(
      data$missingness_setting,
      data$method,
      drop = TRUE
    )
  )

  rows <- lapply(grouped, function(group_data) {
    values <- group_data[[value_name]]
    values <- values[is.finite(values)]
    n_values <- length(values)

    if (n_values == 0) {
      return(NULL)
    }

    mean_value <- mean(values)

    if (n_values < 2) {
      lower <- NA_real_
      upper <- NA_real_
    } else {
      t_value <- qt(
        0.975,
        df = n_values - 1
      )
      half_width <- t_value *
        sd(values) /
        sqrt(n_values)

      lower <- mean_value - half_width
      upper <- mean_value + half_width
    }

    data.frame(
      missingness_setting =
        as.character(group_data$missingness_setting[1]),
      method =
        as.character(group_data$method[1]),
      n = n_values,
      mean = mean_value,
      lower = lower,
      upper = upper
    )
  })

  output <- do.call(
    rbind,
    rows
  )

  output$method <- factor(
    output$method,
    levels = method_order
  )

  output$missingness_setting <- factor(
    output$missingness_setting,
    levels = missingness_settings
  )

  output
}

# ------------------------------------------------------------
# Main figures

# ----------------------------------------------------------
# Figure 1: percent bias

bias_plot_data <- mc_mean_interval(
  plot_data,
  "percent_bias"
)

bias_axis <- make_pretty_axis(
  c(
    bias_plot_data$lower,
    bias_plot_data$upper,
    bias_plot_data$mean
  ),
  include = 0,
  n_breaks = 8
)

figure_bias <- ggplot(
  bias_plot_data,
  aes(
    x = method,
    y = mean
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = lower,
      ymax = upper
    ),
    width = 0.12,
    linewidth = 0.55
  ) +
  geom_point(size = 2.8) +
  facet_wrap(
    ~ missingness_setting,
    nrow = 1,
    labeller = as_labeller(setting_labels)
  ) +
  scale_x_discrete(
    labels = method_labels
  ) +
  scale_y_continuous(
    breaks = bias_axis$breaks,
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  coord_cartesian(
    ylim = bias_axis$limits
  ) +
  labs(
    title = "Bias in the estimated study-hours effect",
    x = NULL,
    y = "Percent bias (%)"
  ) +
  report_theme()

ggsave(
  file.path(
    graph_folder,
    "figure_1_percent_bias.png"
  ),
  figure_bias,
  width = 13.5,
  height = 5.6,
  dpi = 300
)

# ----------------------------------------------------------
# Figure 2: coverage

coverage_groups <- split(
  plot_data,
  interaction(
    plot_data$missingness_setting,
    plot_data$method,
    drop = TRUE
  )
)

coverage_plot_data <- do.call(
  rbind,
  lapply(coverage_groups, function(group_data) {
    covered_values <- group_data$covered
    covered_values <- covered_values[!is.na(covered_values)]
    n_values <- length(covered_values)

    if (n_values == 0) {
      return(NULL)
    }

    coverage_value <- mean(covered_values)

    mc_limits <- qbinom(
      c(0.025, 0.975),
      size = n_values,
      prob = coverage_value
    ) / n_values

    data.frame(
      missingness_setting = as.character(
        group_data$missingness_setting[1]
      ),
      method = as.character(group_data$method[1]),
      n = n_values,
      coverage = coverage_value,
      coverage_lower = mc_limits[1],
      coverage_upper = mc_limits[2],
      coverage_label = paste0(
        sprintf('%.1f', 100 * coverage_value),
        '%'
      )
    )
  })
)

coverage_plot_data$method <- factor(
  coverage_plot_data$method,
  levels = rev(method_order)
)

coverage_plot_data$missingness_setting <- factor(
  coverage_plot_data$missingness_setting,
  levels = missingness_settings
)

coverage_min <- max(
  0.85,
  floor(
    100 * min(
      coverage_plot_data$coverage_lower,
      na.rm = TRUE
    )
  ) / 100 - 0.01
)

coverage_max <- min(
  1.005,
  ceiling(
    100 * max(
      coverage_plot_data$coverage_upper,
      na.rm = TRUE
    )
  ) / 100 + 0.005
)

coverage_breaks <- seq(
  floor(coverage_min * 100),
  100,
  by = 2
) / 100

figure_coverage <- ggplot(
  coverage_plot_data,
  aes(
    y = method,
    x = coverage
  )
) +
  geom_segment(
    aes(
      x = coverage_lower,
      xend = coverage_upper,
      yend = method
    ),
    linewidth = 0.6,
    colour = 'grey35'
  ) +
  geom_segment(
    aes(
      x = coverage_lower,
      xend = coverage_lower,
      y = as.numeric(method) - 0.10,
      yend = as.numeric(method) + 0.10
    ),
    linewidth = 0.6,
    colour = 'grey35',
    inherit.aes = FALSE
  ) +
  geom_segment(
    aes(
      x = coverage_upper,
      xend = coverage_upper,
      y = as.numeric(method) - 0.10,
      yend = as.numeric(method) + 0.10
    ),
    linewidth = 0.6,
    colour = 'grey35',
    inherit.aes = FALSE
  ) +
  geom_vline(
    xintercept = 0.95,
    linetype = 'dashed',
    linewidth = 0.55
  ) +
  geom_point(size = 2.9) +
  geom_label(
    aes(label = coverage_label),
    nudge_x = 0.004,
    hjust = 0,
    size = 3.0,
    label.size = 0,
    fill = 'white'
  ) +
  facet_wrap(
    ~ missingness_setting,
    nrow = 1,
    labeller = as_labeller(setting_labels)
  ) +
  scale_y_discrete(
    labels = method_labels
  ) +
  scale_x_continuous(
    limits = c(coverage_min, coverage_max),
    breaks = coverage_breaks,
    labels = scales::label_percent(
      accuracy = 1
    ),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  labs(
    title = 'Coverage of nominal 95% confidence intervals',
    x = 'Coverage',
    y = NULL
  ) +
  coord_cartesian(clip = 'off') +
  report_theme() +
  theme(
    axis.text.y = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(
    graph_folder,
    'figure_2_coverage.png'
  ),
  figure_coverage,
  width = 13.5,
  height = 5.4,
  dpi = 300
)

# ----------------------------------------------------------
# Figure 3: ICC recovery

true_icc_value <- true_icc(
  icc_study_hours_value
)

icc_plot_data <- variance_summary[
  variance_summary$missingness_setting %in% missingness_settings,
  c(
    "missingness_setting",
    "method",
    "icc_at_study_hours_5"
  )
]

icc_plot_data$method <- factor(
  icc_plot_data$method,
  levels = method_order
)

icc_plot_data$missingness_setting <- factor(
  icc_plot_data$missingness_setting,
  levels = missingness_settings
)

icc_values <- c(
  icc_plot_data$icc_at_study_hours_5,
  true_icc_value
)

icc_lower <- floor(
  (min(icc_values, na.rm = TRUE) - 0.010) / 0.005
) * 0.005

icc_upper <- ceiling(
  (max(icc_values, na.rm = TRUE) + 0.010) / 0.005
) * 0.005

icc_breaks <- seq(
  icc_lower,
  icc_upper,
  by = 0.005
)

figure_icc <- ggplot(
  icc_plot_data,
  aes(
    x = method,
    y = icc_at_study_hours_5
  )
) +
  geom_hline(
    yintercept = true_icc_value,
    linetype = "dashed",
    linewidth = 0.55
  ) +
  geom_point(size = 3.0) +
  geom_text(
    aes(
      label = sprintf(
        "%.3f",
        icc_at_study_hours_5
      )
    ),
    nudge_y = 0.0018,
    size = 3.0
  ) +
  facet_wrap(
    ~ missingness_setting,
    nrow = 1,
    labeller = as_labeller(setting_labels)
  ) +
  scale_x_discrete(
    labels = method_labels
  ) +
  scale_y_continuous(
    breaks = icc_breaks,
    labels = scales::label_number(
      accuracy = 0.001
    ),
    expand = expansion(mult = c(0.03, 0.08))
  ) +
  coord_cartesian(
    ylim = c(
      icc_lower,
      icc_upper
    )
  ) +
  labs(
    title = "Recovery of the school-level ICC",
    x = NULL,
    y = "Estimated ICC"
  ) +
  report_theme()

ggsave(
  file.path(
    graph_folder,
    "figure_3_icc_recovery.png"
  ),
  figure_icc,
  width = 13.5,
  height = 5.4,
  dpi = 300
)

# ----------------------------------------------------------
# Secondary Figure A1: estimate distributions

estimate_axis <- make_pretty_axis(
  plot_data$estimate,
  include = true_beta_study_hours,
  n_breaks = 9
)

figure_estimates <- ggplot(
  plot_data,
  aes(
    x = method,
    y = estimate
  )
) +
  geom_hline(
    yintercept = true_beta_study_hours,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_boxplot(
    width = 0.50,
    fill = "grey92",
    colour = "black",
    linewidth = 0.45,
    outlier.size = 1.5,
    outlier.alpha = 0.55
  ) +
  facet_wrap(
    ~ missingness_setting,
    nrow = 1,
    labeller = as_labeller(setting_labels)
  ) +
  scale_x_discrete(
    labels = method_labels
  ) +
  scale_y_continuous(
    breaks = estimate_axis$breaks,
    labels = scales::label_number(
      accuracy = 0.01
    ),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  coord_cartesian(
    ylim = estimate_axis$limits
  ) +
  labs(
    title = "Distribution of estimated study-hours coefficients",
    x = NULL,
    y = "Estimated coefficient"
  ) +
  report_theme()

ggsave(
  file.path(
    graph_folder,
    "appendix_estimate_distributions.png"
  ),
  figure_estimates,
  width = 13.5,
  height = 5.6,
  dpi = 300
)

# ----------------------------------------------------------
# Secondary Figure A2: empirical vs. model-based standard errors

se_dumbbell_data <- performance_summary[
  performance_summary$missingness_setting %in% missingness_settings,
  c(
    "missingness_setting",
    "method",
    "empirical_se",
    "model_se"
  )
]

se_dumbbell_data$method <- factor(
  se_dumbbell_data$method,
  levels = rev(method_order)
)

se_dumbbell_data$missingness_setting <- factor(
  se_dumbbell_data$missingness_setting,
  levels = missingness_settings
)

se_axis <- make_pretty_axis(
  c(
    se_dumbbell_data$empirical_se,
    se_dumbbell_data$model_se
  ),
  n_breaks = 8,
  lower_bound = 0
)

figure_se_calibration <- ggplot(
  se_dumbbell_data,
  aes(
    y = method
  )
) +
  geom_segment(
    aes(
      x = empirical_se,
      xend = model_se,
      yend = method
    ),
    linewidth = 0.65,
    colour = "grey45"
  ) +
  geom_point(
    aes(
      x = empirical_se,
      shape = "Empirical SE"
    ),
    size = 3.0,
    stroke = 0.7
  ) +
  geom_point(
    aes(
      x = model_se,
      shape = "Model-based SE"
    ),
    size = 3.0,
    stroke = 0.8
  ) +
  facet_wrap(
    ~ missingness_setting,
    nrow = 1,
    labeller = as_labeller(setting_labels)
  ) +
  scale_y_discrete(
    labels = method_labels
  ) +
  scale_x_continuous(
    breaks = se_axis$breaks,
    labels = scales::label_number(
      accuracy = 0.005
    ),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  coord_cartesian(
    xlim = se_axis$limits
  ) +
  scale_shape_manual(
    values = c(
      "Empirical SE" = 16,
      "Model-based SE" = 1
    )
  ) +
  labs(
    title = "Empirical and model-based standard errors",
    x = "Standard error",
    y = NULL,
    shape = NULL
  ) +
  report_theme() +
  theme(
    axis.text.y = element_text(size = 9)
  )

ggsave(
  file.path(
    graph_folder,
    "appendix_standard_errors.png"
  ),
  figure_se_calibration,
  width = 13.5,
  height = 5.5,
  dpi = 300
)

# ----------------------------------------------------------
# Secondary Figure A3: variance-component recovery

variance_component_plot_data <- rbind(
  data.frame(
    missingness_setting = variance_summary$missingness_setting,
    method = variance_summary$method,
    component = "Random-intercept SD",
    recovery_ratio = variance_summary$random_intercept_sd /
      true_random_intercept_sd
  ),
  data.frame(
    missingness_setting = variance_summary$missingness_setting,
    method = variance_summary$method,
    component = "Random-slope SD",
    recovery_ratio = variance_summary$random_slope_sd /
      true_random_slope_sd
  ),
  data.frame(
    missingness_setting = variance_summary$missingness_setting,
    method = variance_summary$method,
    component = "Residual SD",
    recovery_ratio = variance_summary$residual_sd /
      true_residual_sd
  )
)

variance_component_plot_data <- subset(
  variance_component_plot_data,
  missingness_setting %in% missingness_settings
)

variance_component_plot_data$method <- factor(
  variance_component_plot_data$method,
  levels = method_order
)

variance_component_plot_data$missingness_setting <- factor(
  variance_component_plot_data$missingness_setting,
  levels = missingness_settings
)

variance_component_plot_data$component <- factor(
  variance_component_plot_data$component,
  levels = c(
    "Random-intercept SD",
    "Random-slope SD",
    "Residual SD"
  )
)

variance_component_axis <- make_pretty_axis(
  c(
    variance_component_plot_data$recovery_ratio,
    1
  ),
  n_breaks = 8
)

figure_variance_components <- ggplot(
  variance_component_plot_data,
  aes(
    x = recovery_ratio,
    y = component,
    shape = method
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  # Separate variance-component rows.
  geom_hline(
    yintercept = c(1.5, 2.5),
    linewidth = 0.5,
    colour = "grey70"
  ) +
  geom_point(
    size = 3.0,
    stroke = 0.8,
    position = position_dodge(width = 0.55)
  ) +
  facet_wrap(
    ~ missingness_setting,
    nrow = 1,
    labeller = as_labeller(setting_labels)
  ) +
  scale_shape_manual(
    values = c(
      "Complete data" = 16,
      "Complete-case analysis" = 17,
      "glmmTMB" = 15,
      "lmer" = 18
    ),
    labels = method_labels
  ) +
  scale_x_continuous(
    breaks = variance_component_axis$breaks,
    labels = scales::label_number(
      accuracy = 0.05
    ),
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  coord_cartesian(
    xlim = variance_component_axis$limits
  ) +
  labs(
    title = "Recovery of variance components",
    x = "Estimated / true value",
    y = NULL,
    shape = NULL
  ) +
  report_theme() +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(
    graph_folder,
    "appendix_variance_component_recovery.png"
  ),
  figure_variance_components,
  width = 13.5,
  height = 5.8,
  dpi = 300
)

# ----------------------------------------------------------
# Secondary Figure A4: method success rate

success_plot_data <- aggregate(
  success ~ missingness_setting + method,
  data = simulation_results,
  FUN = mean
)

success_plot_data$success_percent <-
  100 * success_plot_data$success

success_plot_data$method <- factor(
  success_plot_data$method,
  levels = method_order
)

success_plot_data$missingness_setting <- factor(
  success_plot_data$missingness_setting,
  levels = missingness_settings
)

figure_success <- ggplot(
  success_plot_data,
  aes(
    x = method,
    y = success_percent
  )
) +
  geom_col(
    width = 0.62,
    fill = "grey80",
    colour = "black",
    linewidth = 0.35
  ) +
  geom_text(
    aes(
      label = paste0(
        round(success_percent, 1),
        "%"
      )
    ),
    vjust = -0.45,
    size = 3.0
  ) +
  facet_wrap(
    ~ missingness_setting,
    nrow = 1,
    labeller = as_labeller(setting_labels)
  ) +
  scale_x_discrete(
    labels = method_labels
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 10),
    limits = c(0, 105),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Successful analysis and imputation runs",
    x = NULL,
    y = "Successful runs (%)"
  ) +
  report_theme()

ggsave(
  file.path(
    graph_folder,
    "appendix_success_rate.png"
  ),
  figure_success,
  width = 13.5,
  height = 5.6,
  dpi = 300
)

# Boxplots with individual benchmark runs.
runtime_plot_data <- runtime_successful

runtime_plot_data$method <- factor(
  runtime_plot_data$method,
  levels = c("glmmTMB", "lmer")
)

runtime_plot_data$missingness_setting <- factor(
  runtime_plot_data$missingness_setting,
  levels = missingness_settings
)

runtime_max <- max(
  runtime_plot_data$elapsed_seconds,
  na.rm = TRUE
)

runtime_upper <- ceiling(runtime_max / 5) * 5

if (!is.finite(runtime_upper) || runtime_upper <= 0) {
  runtime_upper <- 10
}

# To use readable runtime-axis intervals.
runtime_step <- if (runtime_upper <= 60) 5 else 10

runtime_breaks <- seq(
  0,
  runtime_upper,
  by = runtime_step
)

figure_runtime <- ggplot(
  runtime_plot_data,
  aes(
    x = method,
    y = elapsed_seconds
  )
) +
  geom_boxplot(
    width = 0.42,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.65
  ) +
  geom_jitter(
    width = 0.075,
    height = 0,
    size = 2.0,
    alpha = 0.60,
    shape = 16
  ) +
  facet_wrap(
    ~ missingness_setting,
    nrow = 1,
    labeller = as_labeller(setting_labels)
  ) +
  scale_x_discrete(
    labels = method_labels
  ) +
  scale_y_continuous(
    limits = c(0, runtime_upper),
    breaks = runtime_breaks,
    labels = function(x) sprintf("%g", x),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = "Serial runtime benchmark",
    x = NULL,
    y = "Runtime per MI analysis (seconds)"
  ) +
  report_theme() +
  theme(
    axis.text.x = element_text(
      size = 9.5,
      margin = margin(t = 5)
    ),
    axis.text.y = element_text(size = 9.5),
    axis.title.y = element_text(
      face = "bold",
      size = 11
    ),
    strip.text = element_text(
      face = "bold",
      size = 10.5
    ),
    panel.spacing = grid::unit(0.65, "lines")
  )

ggsave(
  file.path(
    graph_folder,
    "figure_4_serial_runtime.png"
  ),
  figure_runtime,
  width = 9.5,
  height = 4.8,
  dpi = 300
)

cat("Graphs saved in: ", normalizePath(graph_folder, winslash = "/", mustWork = FALSE), "\n", sep = "")