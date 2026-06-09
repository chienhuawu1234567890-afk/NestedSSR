############################################################
## Balanced Nested Design Simulation
## Adjusted Blinded Variance Estimation and SSR
############################################################

rm(list = ls())
set.seed(12345)

############################################################
## 1. Design Parameters
############################################################

a  <- 3          # number of treatment groups
b  <- 4          # number of clusters per treatment
n0 <- 12         # planned subjects per cluster

N_plan <- a * b * n0

alpha_sig <- 0.05
target_power <- 0.90

sigma2_true  <- 1.00
sigmaB2_true <- 0.50

theta_values   <- c(0.3, 0.5, 0.7)
sigmaB2_values <- c(0.25, 0.50, 1.00)

nsim <- 10000

############################################################
## 2. Planned Treatment Effects
############################################################

alpha_plan <- c(-0.5, 0, 0.5)
alpha_plan <- alpha_plan - mean(alpha_plan)

delta2_plan <- mean(alpha_plan^2)

############################################################
## 3. Helper Function: Generate Balanced Nested Data
############################################################

generate_nested_data <- function(a, b, n0, alpha, sigma2, sigmaB2) {
  
  dat <- data.frame()
  cluster_id <- 1
  
  for (i in 1:a) {
    for (j in 1:b) {
      
      u_ij <- rnorm(1, mean = 0, sd = sqrt(sigmaB2))
      
      y_ijk <- alpha[i] + u_ij + rnorm(n0, mean = 0, sd = sqrt(sigma2))
      
      temp <- data.frame(
        y = y_ijk,
        treatment = factor(i),
        cluster = factor(cluster_id)
      )
      
      dat <- rbind(dat, temp)
      cluster_id <- cluster_id + 1
    }
  }
  
  return(dat)
}

############################################################
## 4. Helper Function: Take Interim Data
############################################################

take_interim_data <- function(dat, theta) {
  
  n_interim <- max(2, floor(theta * n0))
  
  interim_dat <- do.call(
    rbind,
    lapply(split(dat, dat$cluster), function(x) {
      x[1:n_interim, ]
    })
  )
  
  rownames(interim_dat) <- NULL
  return(interim_dat)
}

############################################################
## 5. Helper Function: ANOVA Components
############################################################

compute_components <- function(dat) {
  
  y <- dat$y
  N <- length(y)
  
  grand_mean <- mean(y)
  
  cluster_mean <- tapply(y, dat$cluster, mean)
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
## 6. Variance Estimators
############################################################

compute_variance_estimators <- function(dat, alpha_plan) {
  
  comp <- compute_components(dat)
  
  SSTO <- comp$SSTO
  SSE  <- comp$SSE
  N    <- comp$N
  
  a <- length(unique(dat$treatment))
  J <- length(unique(dat$cluster))
  
  n_i_star <- as.numeric(table(dat$treatment))
  
  alpha_bar_w_plan <- sum(n_i_star * alpha_plan) / N
  
  treatment_adjustment <- sum(
    n_i_star * (alpha_plan - alpha_bar_w_plan)^2
  )
  
  sigma2_blind <- SSTO / (N - 1)
  
  sigma2_unblind <- SSE / (N - J)
  
  sigma2_adjusted <- (SSTO - treatment_adjustment) / (N - 1)
  
  sigma2_adjusted <- max(sigma2_adjusted, 1e-8)
  
  return(c(
    blind = sigma2_blind,
    unblind = sigma2_unblind,
    adjusted = sigma2_adjusted
  ))
}

############################################################
## 7. Estimate Between-Cluster Variance Under Blinding
############################################################

estimate_sigmaB2_blind <- function(dat) {
  
  cluster_means <- tapply(dat$y, dat$cluster, mean)
  n0_star <- as.numeric(table(dat$cluster))[1]
  J <- length(cluster_means)
  
  MS_cluster <- n0_star * sum((cluster_means - mean(dat$y))^2) / (J - 1)
  
  MSE <- sum(
    unlist(
      lapply(split(dat, dat$cluster), function(x) {
        sum((x$y - mean(x$y))^2)
      })
    )
  ) / (J * (n0_star - 1))
  
  sigmaB2_hat <- (MS_cluster - MSE) / n0_star
  
  return(max(sigmaB2_hat, 0))
}

############################################################
## 8. Sample Size Re-estimation by Noncentral F Inversion
############################################################

reestimate_sample_size <- function(
    sigma2_int,
    sigmaB2_hat,
    delta2_plan,
    a,
    b,
    n0_star,
    N_star,
    alpha_sig = 0.05,
    target_power = 0.90,
    N_max = 2000
) {
  
  nbar_B_star <- n0_star
  
  sigmaA2_hat <- sigma2_int + nbar_B_star * sigmaB2_hat
  
  df1 <- a - 1
  df2 <- a * b - a
  
  fcrit <- qf(1 - alpha_sig, df1, df2)
  
  candidate_N <- N_star:N_max
  
  power_values <- sapply(candidate_N, function(N) {
    lambda_hat <- N * delta2_plan / sigmaA2_hat
    1 - pf(fcrit, df1, df2, ncp = lambda_hat)
  })
  
  ok <- which(power_values >= target_power)
  
  if (length(ok) == 0) {
    N_reest <- NA
  } else {
    N_reest <- candidate_N[min(ok)]
  }
  
  return(list(
    N_reest = N_reest,
    sigmaA2_hat = sigmaA2_hat,
    df1 = df1,
    df2 = df2
  ))
}

############################################################
## 9. One Simulation Run
############################################################

one_simulation <- function(theta, sigmaB2_true) {
  
  dat_full <- generate_nested_data(
    a = a,
    b = b,
    n0 = n0,
    alpha = alpha_plan,
    sigma2 = sigma2_true,
    sigmaB2 = sigmaB2_true
  )
  
  dat_int <- take_interim_data(dat_full, theta)
  
  N_star <- nrow(dat_int)
  n0_star <- as.numeric(table(dat_int$cluster))[1]
  
  var_est <- compute_variance_estimators(dat_int, alpha_plan)
  
  sigmaB2_hat <- estimate_sigmaB2_blind(dat_int)
  
  ssr_blind <- reestimate_sample_size(
    sigma2_int = var_est["blind"],
    sigmaB2_hat = sigmaB2_hat,
    delta2_plan = delta2_plan,
    a = a,
    b = b,
    n0_star = n0_star,
    N_star = N_star,
    alpha_sig = alpha_sig,
    target_power = target_power
  )
  
  ssr_unblind <- reestimate_sample_size(
    sigma2_int = var_est["unblind"],
    sigmaB2_hat = sigmaB2_hat,
    delta2_plan = delta2_plan,
    a = a,
    b = b,
    n0_star = n0_star,
    N_star = N_star,
    alpha_sig = alpha_sig,
    target_power = target_power
  )
  
  ssr_adjusted <- reestimate_sample_size(
    sigma2_int = var_est["adjusted"],
    sigmaB2_hat = sigmaB2_hat,
    delta2_plan = delta2_plan,
    a = a,
    b = b,
    n0_star = n0_star,
    N_star = N_star,
    alpha_sig = alpha_sig,
    target_power = target_power
  )
  
  out <- data.frame(
    theta = theta,
    sigmaB2_true = sigmaB2_true,
    N_star = N_star,
    
    sigma2_blind = var_est["blind"],
    sigma2_unblind = var_est["unblind"],
    sigma2_adjusted = var_est["adjusted"],
    
    sigmaB2_hat = sigmaB2_hat,
    
    sigmaA2_blind = ssr_blind$sigmaA2_hat,
    sigmaA2_unblind = ssr_unblind$sigmaA2_hat,
    sigmaA2_adjusted = ssr_adjusted$sigmaA2_hat,
    
    N_blind = ssr_blind$N_reest,
    N_unblind = ssr_unblind$N_reest,
    N_adjusted = ssr_adjusted$N_reest
  )
  
  return(out)
}

############################################################
## 10. Full Simulation Loop
############################################################

all_results <- data.frame()

for (theta in theta_values) {
  for (sigmaB2 in sigmaB2_values) {
    
    cat("Running theta =", theta, 
        " sigmaB2 =", sigmaB2, "\n")
    
    temp_results <- do.call(
      rbind,
      replicate(
        nsim,
        one_simulation(theta = theta, sigmaB2_true = sigmaB2),
        simplify = FALSE
      )
    )
    
    all_results <- rbind(all_results, temp_results)
  }
}

############################################################
## 11. Summary of Variance Estimators
############################################################

variance_summary <- aggregate(
  cbind(
    sigma2_blind,
    sigma2_unblind,
    sigma2_adjusted,
    sigmaB2_hat,
    sigmaA2_blind,
    sigmaA2_unblind,
    sigmaA2_adjusted
  ) ~ theta + sigmaB2_true,
  data = all_results,
  FUN = mean
)

print(variance_summary)

############################################################
## 12. Summary of Re-estimated Sample Size
############################################################

sample_size_summary <- aggregate(
  cbind(
    N_blind,
    N_unblind,
    N_adjusted
  ) ~ theta + sigmaB2_true,
  data = all_results,
  FUN = function(x) {
    c(
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      q025 = quantile(x, 0.025, na.rm = TRUE),
      q975 = quantile(x, 0.975, na.rm = TRUE)
    )
  }
)

print(sample_size_summary)

############################################################
## 13. Convert Sample Size Summary to Clean Data Frame
############################################################

clean_sample_size_summary <- do.call(
  data.frame,
  sample_size_summary
)

print(clean_sample_size_summary)

############################################################
## 14. Save Results
############################################################

write.csv(all_results,
          "balanced_nested_full_simulation_results.csv",
          row.names = FALSE)

write.csv(variance_summary,
          "balanced_nested_variance_summary.csv",
          row.names = FALSE)

write.csv(clean_sample_size_summary,
          "balanced_nested_sample_size_summary.csv",
          row.names = FALSE)

############################################################
## End of Simulation Code
############################################################
