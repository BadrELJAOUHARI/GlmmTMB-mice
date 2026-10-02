mice.impute.2l.glmmTMB <- function(y,
                                   ry,
                                   x,
                                   type,
                                   wy = NULL,
                                   ...) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("Package 'glmmTMB' is required.")
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required.")
  }

  if (is.null(wy)) {
    wy <- !ry
  }

  x <- as.data.frame(x, check.names = FALSE)

  if (is.null(names(type))) {
    names(type) <- colnames(x)
  }

  x <- cbind("(Intercept)" = 1, x)
  type <- c("(Intercept)" = 2, type)

  cluster_name <- names(type[type == -2])
  random_names <- names(type[type == 2])
  fixed_names <- names(type[type > 0])

  if (length(cluster_name) != 1) {
    stop("Exactly one cluster variable must be coded as -2.")
  }

  cluster_name <- cluster_name[[1]]

  cluster_values <- unique(x[[cluster_name]])

  x_fixed <- as.matrix(x[, fixed_names, drop = FALSE])
  z_random <- as.matrix(x[, random_names, drop = FALSE])

  x_observed <- x[ry, , drop = FALSE]
  y_observed <- y[ry]
  x_fixed_observed <- x_fixed[ry, , drop = FALSE]
  z_random_observed <- z_random[ry, , drop = FALSE]

  fixed_terms <- setdiff(fixed_names, "(Intercept)")
  random_terms <- setdiff(random_names, "(Intercept)")

  fixed_part <- if (length(fixed_terms) == 0) {
    "1"
  } else {
    paste(fixed_terms, collapse = " + ")
  }

  random_part <- if (length(random_terms) == 0) {
    "1"
  } else {
    paste("1 +", paste(random_terms, collapse = " + "))
  }

  imputation_formula <- as.formula(
    paste0(
      "y_observed ~ ",
      fixed_part,
      " + (",
      random_part,
      " | ",
      cluster_name,
      ")"
    )
  )

  fit_data <- data.frame(
    y_observed = y_observed,
    x_observed,
    check.names = FALSE
  )
  fit_data[[cluster_name]] <- factor(fit_data[[cluster_name]])

  fit <- try(
    glmmTMB::glmmTMB(
      formula = imputation_formula,
      data = fit_data,
      family = gaussian(),
      REML = FALSE,
      ...
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    stop("glmmTMB imputation model failed to fit.")
  }


  # Draw residual variance.

  sigma_hat <- sigma(fit)
  beta_hat <- glmmTMB::fixef(fit)$cond

  residual_df <- nrow(fit_data) - length(beta_hat)

  if (
    !is.finite(sigma_hat) ||
      sigma_hat <= 0 ||
      !is.finite(residual_df) ||
      residual_df <= 0
  ) {
    stop("Residual variance could not be drawn from the glmmTMB imputation model.")
  }

  sigma2_star <- residual_df * sigma_hat^2 / rchisq(1, residual_df)

  # Draw fixed effects.

  beta_covariance <- try(
    as.matrix(vcov(fit)$cond) * (sigma2_star / sigma_hat^2),
    silent = TRUE
  )

  if (
    inherits(beta_covariance, "try-error") ||
      any(!is.finite(beta_covariance))
  ) {
    stop("Fixed-effect covariance was unavailable for the glmmTMB parameter draw.")
  }

  beta_covariance <- (beta_covariance + t(beta_covariance)) / 2

  beta_chol <- try(
    chol(beta_covariance),
    silent = TRUE
  )

  if (inherits(beta_chol, "try-error")) {
    stop("Fixed-effect covariance was not positive definite, so beta could not be drawn.")
  }

  beta_star <- as.numeric(
    beta_hat +
      t(beta_chol) %*% rnorm(length(beta_hat))
  )
  names(beta_star) <- names(beta_hat)

  beta_star <- beta_star[fixed_names]

  if (anyNA(beta_star)) {
    stop("Could not align the fixed-effect draw with the imputation design matrix.")
  }

  # Draw random-effect covariance.

  random_effect_estimates <- try(
    as.matrix(glmmTMB::ranef(fit)$cond[[cluster_name]]),
    silent = TRUE
  )

  if (
    inherits(random_effect_estimates, "try-error") ||
      nrow(random_effect_estimates) == 0
  ) {
    stop("Could not extract random effects from the glmmTMB imputation model.")
  }

  random_effect_estimates <- random_effect_estimates[
    ,
    random_names,
    drop = FALSE
  ]

  lambda <- t(random_effect_estimates) %*% random_effect_estimates
  df_psi <- nrow(random_effect_estimates)

  temp_psi_star <- stats::rWishart(
    1,
    df_psi,
    diag(nrow(lambda))
  )[, , 1]

  temp <- MASS::ginv(lambda)
  eigen_temp <- eigen(temp)

  if (sum(eigen_temp$values > 0) == length(eigen_temp$values)) {
    decomposition <- eigen_temp$vectors %*%
      diag(
        sqrt(eigen_temp$values),
        nrow = length(eigen_temp$values)
      )

    psi_star <- MASS::ginv(
      decomposition %*% temp_psi_star %*% t(decomposition)
    )
  } else {
    temp_svd <- try(svd(lambda), silent = TRUE)

    if (!inherits(temp_svd, "try-error")) {
      decomposition <- temp_svd$u %*%
        diag(
          sqrt(temp_svd$d),
          nrow = length(temp_svd$d)
        )

      psi_star <- MASS::ginv(
        decomposition %*% temp_psi_star %*% t(decomposition)
      )
    } else {
      psi_star <- temp
      warning("psi fixed to estimate")
    }
  }


  # Draw cluster effects and impute missing values.

  imputed_y <- y
  observed_clusters <- unique(x_observed[[cluster_name]])

  for (cluster_value in cluster_values) {
    missing_in_cluster <- wy & x[[cluster_name]] == cluster_value

    if (!any(missing_in_cluster)) {
      next
    }

    if (cluster_value %in% observed_clusters) {
      observed_in_cluster <- x_observed[[cluster_name]] == cluster_value

      x_cluster <- x_fixed_observed[
        observed_in_cluster,
        ,
        drop = FALSE
      ]

      z_cluster <- z_random_observed[
        observed_in_cluster,
        ,
        drop = FALSE
      ]

      y_cluster <- y_observed[observed_in_cluster]

      sigma2_matrix <- diag(
        sigma2_star,
        nrow = nrow(z_cluster)
      )

      gain_matrix <- psi_star %*%
        t(z_cluster) %*%
        MASS::ginv(
          z_cluster %*%
            psi_star %*%
            t(z_cluster) +
            sigma2_matrix
        )

      random_mean <- gain_matrix %*%
        (y_cluster - x_cluster %*% beta_star)

      random_covariance <- psi_star -
        gain_matrix %*%
        z_cluster %*%
        psi_star
    } else {
      random_mean <- matrix(
        0,
        nrow = nrow(psi_star),
        ncol = 1
      )
      random_covariance <- psi_star
    }

    # Same random-effect draw used by mice::2l.lmer (eigen first then SVD).
    random_covariance <- random_covariance -
      upper.tri(random_covariance) * random_covariance +
      t(lower.tri(random_covariance) * random_covariance)

    random_eigen <- eigen(random_covariance)

    if (sum(random_eigen$values > 0) == length(random_eigen$values)) {
      transform <- random_eigen$vectors %*%
        sqrt(diag(
          random_eigen$values,
          nrow = length(random_eigen$values)
        ))

      random_draw <- as.numeric(
        random_mean +
          transform %*% rnorm(length(random_mean))
      )
    } else {
      random_svd <- try(svd(random_covariance), silent = TRUE)

      if (!inherits(random_svd, "try-error")) {
        transform <- random_svd$u %*%
          sqrt(diag(
            random_svd$d,
            nrow = length(random_svd$d)
          ))

        random_draw <- as.numeric(
          random_mean +
            transform %*% rnorm(length(random_mean))
        )
      } else {
        random_draw <- as.numeric(random_mean)
        warning("Random effect fixed to estimate")
      }
    }

    x_missing <- x_fixed[
      missing_in_cluster,
      ,
      drop = FALSE
    ]

    z_missing <- z_random[
      missing_in_cluster,
      ,
      drop = FALSE
    ]

    imputed_y[missing_in_cluster] <- as.numeric(
      x_missing %*% beta_star +
        z_missing %*% random_draw +
        rnorm(
          sum(missing_in_cluster),
          mean = 0,
          sd = sqrt(sigma2_star)
        )
    )
  }

  imputed_y[wy]
}
