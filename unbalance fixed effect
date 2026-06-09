############################################################
## Unbalanced Nested Fixed-Effect Design Simulation
## Revised parameters to avoid SD = 0
############################################################

rm(list = ls())
set.seed(12345)

############################################################
## 1. Design Parameters
############################################################

a <- 3
b <- 4

## Smaller and less extreme cluster sizes
n1 <- 8
n_vec <- c(n1, 2 * n1, 2 * n1)

N_plan <- b * sum(n_vec)

alpha_sig <- 0.05
target_power <- 0.95

## Larger residual variance to increase uncertainty
sigma2_true <- 2.50

## Smaller interim information fractions
theta_values <- c(0.25, 0.40, 0.55)

beta_scale_values <- c(0.25, 0.50, 1.00)

nsim <- 10000

############################################################
## 2. Planned Treatment Effects
############################################################

## Smaller treatment effects to avoid saturation
alpha_plan <- c(-0.25, 0, 0.25)
alpha_plan <- alpha_plan - mean(alpha_plan)

n_i_plan <- b * n_vec
w_plan <- n_i_plan / sum(n_i_plan)

delta2_plan <- sum(w_plan * alpha_plan^2)

############################################################
## 3. Generate Fixed Nested Effects beta_{j(i)}
############################################################

generate_beta_fixed <- function(a, b, beta_scale) {
  
  beta_mat <- matrix(0, nrow = a, ncol = b)
  
  for (i in 1:a) {
    temp <- seq(-1, 1, length.out = b)
    temp <- temp - mean(temp)
    beta_mat[i, ] <- beta_scale * temp
  }
  
  return(beta_mat)
}

############################################################
## 4. Generate Unbalanced Nested Fixed-Effect Data
############################################################

generate_nested_fixed_data_unbalanced <- function(a, b, n_vec, alpha, beta_mat, sigma2) {
  
  dat <- data.frame()
  cluster_id <- 1
  
  for (i in 1:a) {
    for (j in 1:b) {
      
      n_ij <- n_vec[i]
      
      y_ijk <- alpha[i] + beta_mat[i, j] +
        rnorm(n_ij, mean = 0, sd = sqrt(sigma2))
      
      temp <- data.frame(
        y = y_ijk,
        treatment = factor(i),
        cluster = factor(cluster_id),
        cluster_within_treatment = factor(j),
        n_ij_plan = n_ij
      )
      
      dat <- rbind(dat, temp)
      cluster_id <- cluster_id + 1
    }
  }
  
  return(dat)
}

############################################################
## 5. Take Interim Data
############################################################

take_interim_data_unbalanced <- function(dat, theta) {
  
  interim_dat <- do.call(
    rbind,
    lapply(split(dat, dat$cluster), function(x) {
      n_int <- max(2, floor(theta * nrow(x)))
      x[1:n_int, ]
    })
  )
  
  rownames(interim_dat) <- NULL
  return(interim_dat)
}

############################################################
## 6. ANOVA Components
############################################################

compute_components <- function(dat) {
  
  y <- dat$y
  N <- length(y)
  
  grand_mean <- mean(y)
  treatment_mean <- tapply(y, dat$treatment, mean)
  
  SSTO <- sum((y - grand_mean)^2)
  
  SSE <- sum(
    unlist(
      lapply(split(dat, dat$cluster), function(x) {
        sum((x$y - mean(x$y))^2)
      })
    )
  )
  
  SSB_A <- sum(
    unlist(
      lapply(split(dat, dat$cluster), function(x) {
        trt <- as.character(unique(x$treatment))
        n_j <- nrow(x)
        n_j * (mean(x$y) - treatment_mean[trt])^2
      })
    )
  )
  
  return(list(
    SSTO = SSTO,
    SSE = SSE,
    SSB_A = SSB_A,
    N = N
  ))
}

############################################################
## 7. Variance Estimators
############################################################

compute_variance_estimators_unbalanced_fixed <- function(dat, alpha_plan, beta_mat) {
  
  comp <- compute_components(dat)
  
  SSTO <- comp$SSTO
  SSE  <- comp$SSE
  N    <- comp$N
  
  J <- length(unique(dat$cluster))
  
  cluster_info <- unique(dat[, c(
    "treatment",
    "cluster",
    "cluster_within_treatment"
  )])
  
  fixed_mean <- numeric(nrow(cluster_info))
  n_ij_star <- numeric(nrow(cluster_info))
  
  for (r in 1:nrow(cluster_info)) {
    
    i <- as.numeric(as.character(cluster_info$treatment[r]))
    j <- as.numeric(as.character(cluster_info$cluster_within_treatment[r]))
    cl <- cluster_info$cluster[r]
    
    n_ij_star[r] <- sum(dat$cluster == cl)
    fixed_mean[r] <- alpha_plan[i] + beta_mat[i, j]
  }
  
  eta_bar_w_plan <- sum(n_ij_star * fixed_mean) / N
  
  fixed_adjustment <- sum(
    n_ij_star * (fixed_mean - eta_bar_w_plan)^2
  )
  
  sigma2_blind <- SSTO / (N - 1)
  sigma2_unblind <- SSE / (N - J)
  
  sigma2_adjusted <- (SSTO - fixed_adjustment) / (N - 1)
  sigma2_adjusted <- max(sigma2_adjusted, 1e-8)
  
  return(c(
    blind = sigma2_blind,
    unblind = sigma2_unblind,
    adjusted = sigma2_adjusted
  ))
}

############################################################
## 8. Estimate Fixed Nested Component Under Blinding
############################################################

estimate_fixed_nested_component_blind_unbalanced <- function(dat) {
  
  cluster_means <- tapply(dat$y, dat$cluster, mean)
  n_j <- as.numeric(table(dat$cluster))
  
  J <- length(cluster_means)
  N <- sum(n_j)
  
  grand_mean <- sum(n_j * cluster_means) / N
  
  MS_cluster <- sum(n_j * (cluster_means - grand_mean)^2) / (J - 1)
  
  MSE <- sum(
    unlist(
      lapply(split(dat, dat$cluster), function(x) {
        sum((x$y - mean(x$y))^2)
      })
    )
  ) / (N - J)
  
  n_eff <- (N - sum(n_j^2) / N) / (J - 1)
  
  beta_component_hat <- (MS_cluster - MSE) / n_eff
  beta_component_hat <- max(beta_component_hat, 0)
  
  return(beta_component_hat)
}

############################################################
## 9. Sample Size Re-estimation by Noncentral F Inversion
############################################################

reestimate_sample_size_unbalanced_fixed <- function(
    sigma2_int,
    delta2_plan,
    a,
    b,
    N_interim,
    alpha_sig = 0.05,
    target_power = 0.95,
    N_max = 5000
) {
  
  df1 <- a - 1
  J <- a * b
  
  N_min_valid <- J + df1 + 1
  candidate_N <- N_min_valid:N_max
  
  power_values <- sapply(candidate_N, function(N) {
    
    df2 <- N - J
    
    if (df2 <= 0) return(0)
    
    lambda_hat <- N * delta2_plan / sigma2_int
    
    fcrit <- qf(1 - alpha_sig, df1, df2)
    
    1 - pf(fcrit, df1, df2, ncp = lambda_hat)
  })
  
  ok <- which(power_values >= target_power)
  
  if (length(ok) == 0) {
    
    N_required <- NA
    N_final <- NA
    
  } else {
    
    N_required <- candidate_N[min(ok)]
    
    ## Keep final sample size at least as large as interim sample size
    ## but do not force rounding to multiples of a*b.
    N_final <- max(N_interim, ceiling(N_required))
  }
  
  return(list(
    N_required = N_required,
    N_final = N_final,
    sigmaA2_hat = sigma2_int,
    df1 = df1
  ))
}

############################################################
## 10. One Simulation Run
############################################################

one_simulation <- function(theta, beta_scale) {
  
  beta_mat <- generate_beta_fixed(
    a = a,
    b = b,
    beta_scale = beta_scale
  )
  
  dat_full <- generate_nested_fixed_data_unbalanced(
    a = a,
    b = b,
    n_vec = n_vec,
    alpha = alpha_plan,
    beta_mat = beta_mat,
    sigma2 = sigma2_true
  )
  
  dat_int <- take_interim_data_unbalanced(dat_full, theta)
  
  N_interim <- nrow(dat_int)
  
  var_est <- compute_variance_estimators_unbalanced_fixed(
    dat = dat_int,
    alpha_plan = alpha_plan,
    beta_mat = beta_mat
  )
  
  beta_component_hat <- estimate_fixed_nested_component_blind_unbalanced(dat_int)
  
  ssr_blind <- reestimate_sample_size_unbalanced_fixed(
    sigma2_int = var_est["blind"],
    delta2_plan = delta2_plan,
    a = a,
    b = b,
    N_interim = N_interim,
    alpha_sig = alpha_sig,
    target_power = target_power
  )
  
  ssr_unblind <- reestimate_sample_size_unbalanced_fixed(
    sigma2_int = var_est["unblind"],
    delta2_plan = delta2_plan,
    a = a,
    b = b,
    N_interim = N_interim,
    alpha_sig = alpha_sig,
    target_power = target_power
  )
  
  ssr_adjusted <- reestimate_sample_size_unbalanced_fixed(
    sigma2_int = var_est["adjusted"],
    delta2_plan = delta2_plan,
    a = a,
    b = b,
    N_interim = N_interim,
    alpha_sig = alpha_sig,
    target_power = target_power
  )
  
  out <- data.frame(
    theta = theta,
    beta_scale = beta_scale,
    N_interim = N_interim,
    
    sigma2_blind = var_est["blind"],
    sigma2_unblind = var_est["unblind"],
    sigma2_adjusted = var_est["adjusted"],
    
    beta_component_hat = beta_component_hat,
    
    N_required_blind = ssr_blind$N_required,
    N_required_unblind = ssr_unblind$N_required,
    N_required_adjusted = ssr_adjusted$N_required,
    
    N_final_blind = ssr_blind$N_final,
    N_final_unblind = ssr_unblind$N_final,
    N_final_adjusted = ssr_adjusted$N_final
  )
  
  return(out)
}

############################################################
## 11. Full Simulation Loop
############################################################

all_results <- data.frame()

for (theta in theta_values) {
  for (beta_scale in beta_scale_values) {
    
    cat("Running theta =", theta,
        " beta_scale =", beta_scale, "\n")
    
    temp_results <- do.call(
      rbind,
      replicate(
        nsim,
        one_simulation(theta = theta, beta_scale = beta_scale),
        simplify = FALSE
      )
    )
    
    all_results <- rbind(all_results, temp_results)
  }
}

############################################################
## 12. Summary of Variance Estimators
############################################################

variance_summary <- aggregate(
  cbind(
    sigma2_blind,
    sigma2_unblind,
    sigma2_adjusted,
    beta_component_hat
  ) ~ theta + beta_scale,
  data = all_results,
  FUN = mean
)

print(variance_summary)

############################################################
## 13. Summary Function
############################################################

summary_fun <- function(x) {
  c(
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    q025 = quantile(x, 0.025, na.rm = TRUE),
    q975 = quantile(x, 0.975, na.rm = TRUE)
  )
}

############################################################
## 14. Summary of Required Sample Size
############################################################

sample_size_summary_required <- aggregate(
  cbind(
    N_required_blind,
    N_required_unblind,
    N_required_adjusted
  ) ~ theta + beta_scale,
  data = all_results,
  FUN = summary_fun
)

clean_sample_size_summary_required <- do.call(
  data.frame,
  sample_size_summary_required
)

print(clean_sample_size_summary_required)

############################################################
## 15. Summary of Final Re-estimated Sample Size
############################################################

sample_size_summary_final <- aggregate(
  cbind(
    N_final_blind,
    N_final_unblind,
    N_final_adjusted
  ) ~ theta + beta_scale,
  data = all_results,
  FUN = summary_fun
)

clean_sample_size_summary_final <- do.call(
  data.frame,
  sample_size_summary_final
)

print(clean_sample_size_summary_final)

############################################################
## 16. Save Results
############################################################

write.csv(
  all_results,
  "unbalanced_nested_fixed_effect_full_simulation_results_no_zero_sd.csv",
  row.names = FALSE
)

write.csv(
  variance_summary,
  "unbalanced_nested_fixed_effect_variance_summary_no_zero_sd.csv",
  row.names = FALSE
)

write.csv(
  clean_sample_size_summary_required,
  "unbalanced_nested_fixed_effect_required_sample_size_summary_no_zero_sd.csv",
  row.names = FALSE
)

write.csv(
  clean_sample_size_summary_final,
  "unbalanced_nested_fixed_effect_final_sample_size_summary_no_zero_sd.csv",
  row.names = FALSE
)

############################################################
## End of Revised Simulation Code
############################################################
