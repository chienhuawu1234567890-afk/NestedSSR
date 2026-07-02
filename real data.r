############################################################
## ACTG 175 Real-Data Illustration
## Proposed Method + Existing Methods
## Removed: Conventional blinded one-way residual
## Save all outputs in one directory
############################################################

rm(list = ls())

############################################################
## 0. Output directory
############################################################

output_dir <- "ACTG175_outputs"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

############################################################
## 1. Install and load required packages
############################################################

pkgs <- c("speff2trial", "lme4", "dplyr")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
  library(p, character.only = TRUE)
}

############################################################
## 2. Load ACTG 175 data
############################################################

data("ACTG175", package = "speff2trial")
dat <- ACTG175

dat <- dat %>%
  mutate(
    Treatment = factor(arms),
    Outcome   = cd420 - cd40
  )

############################################################
## 3. Define center variable
############################################################

if ("site" %in% names(dat)) {
  dat$Center <- factor(dat$site)
} else if ("strat" %in% names(dat)) {
  dat$Center <- factor(dat$strat)
} else {
  set.seed(123)
  dat$Center <- factor(sample(1:30, nrow(dat), replace = TRUE))
}

############################################################
## 4. Save cleaned full dataset
############################################################

write.csv(
  dat,
  file = file.path(output_dir, "ACTG175_cleaned_full_data.csv"),
  row.names = FALSE
)

############################################################
## 5. Final mixed-effects model
############################################################

fit_final <- lmer(
  Outcome ~ Treatment + (1 | Center),
  data = dat
)

sink(file.path(output_dir, "final_mixed_model_summary.txt"))
cat("Final mixed-effects model\n")
cat("=========================\n\n")
print(summary(fit_final))

cat("\n\nANOVA table\n")
cat("===========\n\n")
print(anova(fit_final))
sink()

############################################################
## 6. Create interim data
############################################################

theta <- 0.50
Nstar <- floor(theta * nrow(dat))

set.seed(2026)
interim_index <- sample(seq_len(nrow(dat)), Nstar)
interim <- dat[interim_index, ]

write.csv(
  interim,
  file = file.path(output_dir, "ACTG175_interim_data.csv"),
  row.names = FALSE
)

############################################################
## 7. Planned treatment effects
############################################################

trt_levels <- levels(dat$Treatment)

alpha_plan <- setNames(
  seq(-15, 15, length.out = length(trt_levels)),
  trt_levels
)

planned_effects <- data.frame(
  Treatment = names(alpha_plan),
  Planned_Effect = as.numeric(alpha_plan)
)

write.csv(
  planned_effects,
  file = file.path(output_dir, "planned_treatment_effects.csv"),
  row.names = FALSE
)

############################################################
## 8. Method 1: Existing naive blinded estimator
############################################################

SSTO_star <- sum(
  (interim$Outcome - mean(interim$Outcome))^2
)

sigma2_naive_blind <- SSTO_star / (Nstar - 1)

############################################################
## 9. Method 2: Existing unblinded residual benchmark
## Uses treatment and center information
############################################################

fit_unblind <- lm(
  Outcome ~ Treatment + Center,
  data = interim
)

sigma2_unblind <- sum(resid(fit_unblind)^2) /
  df.residual(fit_unblind)

############################################################
## 10. Method 3: Proposed adjusted blinded estimator
## Random center/nested effect formulation
############################################################

n_i_star <- table(interim$Treatment)

alpha_bar_star <- sum(
  n_i_star * alpha_plan[names(n_i_star)]
) / Nstar

treatment_adjustment <- sum(
  n_i_star *
    (alpha_plan[names(n_i_star)] - alpha_bar_star)^2
)

sigma2_adj_random <-
  (SSTO_star - treatment_adjustment) / (Nstar - 1)

sigma2_adj_random <- max(sigma2_adj_random, 0)

############################################################
## 11. Compare variance estimators
############################################################

variance_results <- data.frame(
  Method = c(
    "Existing: Naive blinded total variance",
    "Existing: Unblinded residual benchmark",
    "Proposed: Adjusted blinded random-center"
  ),
  Variance = c(
    sigma2_naive_blind,
    sigma2_unblind,
    sigma2_adj_random
  )
)

print(variance_results)

write.csv(
  variance_results,
  file = file.path(output_dir, "variance_estimator_results.csv"),
  row.names = FALSE
)

############################################################
## 12. Re-estimated sample size using noncentral F
############################################################

alpha_level <- 0.05
target_power <- 0.80
a <- length(trt_levels)

delta2_plan <- mean(alpha_plan^2)

power_nested <- function(N, Vhat, a, delta2,
                         alpha_level = 0.05) {
  
  df1 <- a - 1
  df2 <- N - a
  
  if (df2 <= 0 || Vhat <= 0) {
    return(NA_real_)
  }
  
  lambda <- N * delta2 / Vhat
  
  Fcrit <- qf(
    1 - alpha_level,
    df1 = df1,
    df2 = df2
  )
  
  power <- 1 - pf(
    Fcrit,
    df1 = df1,
    df2 = df2,
    ncp = lambda
  )
  
  return(power)
}

find_N <- function(Vhat, a, delta2,
                   target_power = 0.80,
                   alpha_level = 0.05,
                   N_min = 50,
                   N_max = 10000) {
  
  for (N in N_min:N_max) {
    
    pw <- power_nested(
      N = N,
      Vhat = Vhat,
      a = a,
      delta2 = delta2,
      alpha_level = alpha_level
    )
    
    if (!is.na(pw) && pw >= target_power) {
      return(N)
    }
  }
  
  return(NA_integer_)
}

N_results <- data.frame(
  Method = variance_results$Method,
  Variance = variance_results$Variance,
  Reestimated_N = c(
    find_N(sigma2_naive_blind, a, delta2_plan),
    find_N(sigma2_unblind, a, delta2_plan),
    find_N(sigma2_adj_random, a, delta2_plan)
  )
)

print(N_results)

write.csv(
  N_results,
  file = file.path(output_dir, "reestimated_sample_size_results.csv"),
  row.names = FALSE
)

############################################################
## 13. Save interim unblinded model summary
############################################################

sink(file.path(output_dir, "interim_unblinded_model_summary.txt"))
cat("Interim unblinded treatment + center model\n")
cat("=========================================\n\n")
print(summary(fit_unblind))
sink()

############################################################
## 14. Save model objects
############################################################

saveRDS(
  fit_final,
  file = file.path(output_dir, "final_mixed_model.rds")
)

saveRDS(
  fit_unblind,
  file = file.path(output_dir, "interim_unblinded_model.rds")
)

############################################################
## 15. Save analysis settings
############################################################

analysis_settings <- data.frame(
  Setting = c(
    "Total sample size",
    "Interim fraction theta",
    "Interim sample size",
    "Number of treatment arms",
    "Treatment levels",
    "Alpha level",
    "Target power",
    "Planned delta squared"
  ),
  Value = c(
    nrow(dat),
    theta,
    Nstar,
    a,
    paste(trt_levels, collapse = ", "),
    alpha_level,
    target_power,
    delta2_plan
  )
)

write.csv(
  analysis_settings,
  file = file.path(output_dir, "analysis_settings.csv"),
  row.names = FALSE
)

############################################################
## 16. Save all important results together
############################################################

all_results <- list(
  variance_results = variance_results,
  sample_size_results = N_results,
  planned_effects = planned_effects,
  analysis_settings = analysis_settings,
  final_model = fit_final,
  interim_unblinded_model = fit_unblind
)

saveRDS(
  all_results,
  file = file.path(output_dir, "all_results.rds")
)

save.image(
  file = file.path(output_dir, "ACTG175_full_workspace.RData")
)

############################################################
## 17. Completion message
############################################################

cat("\nAnalysis completed successfully.\n")
cat("Removed method: Existing conventional blinded one-way residual.\n")
cat("All outputs were saved in:\n")
cat(normalizePath(output_dir), "\n")
