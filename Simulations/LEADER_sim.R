library(readr)
library(survival)
library(dplyr)
library(ggplot2)
library(tidyr)
library(WeightIt)
library(SuperLearner)
library(geepack) # Required for pseudo-observations

generate_data <- function(n_total = 1000, cuts, scenario = "baseline") {
  
  # --- Covariates ---
  x1 <- rnorm(n_total, mean = -0.5, sd = 1)
  x2 <- rbinom(n_total, 1, 0.5)
  u <- rnorm(n_total, 0.5, 1) # Unmeasured confounder for scenario 4
  
  # --- Selection Model Parameters (S=1 is source) ---
  alpha <- list(
    baseline = c(-0.4, 0.5, -0.5, 0),
    outcome_wrong = c(-0.4, 0.5, -0.5, 0),
    selection_wrong = c(-0.4, 0.5, -0.5, 0.6),
    both_wrong = c(-0.4, 0.5, -0.5, 0.6),
    positivity_practical = c(-0.4, 2, -2, 0), # Induce practical violation
    positivity_structural = c(-0.4, 0.5, -0.6, 0) # Rule applied later
  )
  
  params_alpha <- alpha[[scenario]]
  linear_pred_s <- params_alpha[1] + params_alpha[2] * x1 + params_alpha[3] * x2 + params_alpha[4] * u
  
  pi <- plogis(linear_pred_s)
  S <- rbinom(n_total, 1, pi)
  
  # Apply structural positivity violation
  if (scenario == "positivity_structural") {
    S[x1 > 0] <- 0
  }
  
  # --- Treatment Assignment (in source S=1 only) ---
  A <- rbinom(n_total, 1, 0.5)
  
  # --- Outcome Model Parameters (Weibull PH) ---
  beta <- list(
    baseline = c(0.3, -0.4, -0.5, 0, 0),
    outcome_wrong = c(0.3, -0.4, -0.5, 0.4, 0.7),
    selection_wrong = c(0.3, -0.4, -0.5, 0, 0),
    both_wrong = c(0.3, -0.4, -0.5, 0.4, 0.7),
    positivity_practical = c(0.3, -0.4, -0.5, 0, 0),
    positivity_structural = c(0.3, -0.4, -0.5, 0, 0)
  )
  
  params_beta <- beta[[scenario]]
  weibull_shape <- 1.5 # nu
  
  # --- Generate Potential Survival Times T ---
  v <- runif(n_total)
  
  linear_pred_t <- params_beta[1] * x1 +
    params_beta[2] * x2 + 
    params_beta[3] * A + 
    params_beta[4] * u + 
    params_beta[5] * A * u
  
  # Invert CDF: T = (-log(v) / (lambda0 * exp(eta)))^(1/nu)
  lambda0 <- 0.2
  T_event <- (-log(v) / (lambda0 * exp(linear_pred_t)))^(1 / weibull_shape)
  
  # --- Censoring Model ---
  gamma <- c(-3, 0.2)
  cens_rate <- exp(gamma[1])
  C_time <- rexp(n_total, rate = cens_rate)
  
  true_surv_func <- function(t, shape, lambda0,
                             A, S, x1, x2, u,
                             params_beta, scenario) {
    
    lp <- params_beta[1] * x1 +
      params_beta[2] * x2 + 
      params_beta[3] * A + 
      params_beta[4] * u + 
      params_beta[5] * A * u
    
    lambda <- lambda0*exp(lp)
    return(exp(-lambda * t^shape))
    
  }
  
  # Calculate survival probabilities for each person at each time point
  results <- lapply(cuts, function(t) {
    
    surv_0_vals <- true_surv_func(t, A = 0, S = S, x1 = x1, x2 = x2, u = u, shape = weibull_shape,
                                  lambda0 = lambda0, params_beta = params_beta, scenario = scenario)
    surv_1_vals <- true_surv_func(t, A = 1, S = S, x1 = x1, x2 = x2, u = u, lambda0 = lambda0,
                                  shape = weibull_shape, params_beta = params_beta, scenario = scenario)
    
    # Average over the large population to get the true marginal survival probabilities
    true_surv_0 <- mean(surv_0_vals)
    true_surv_1 <- mean(surv_1_vals)
    
    data.frame(
      time = t,
      surv0 = true_surv_0,
      surv1 = true_surv_1,
      diff = true_surv_1 - true_surv_0
    )
    
  })
  
  truth <- do.call(rbind, results)
  
  # --- Observed Data ---
  time <- pmin(T_event, C_time)
  status <- as.numeric(T_event <= C_time)
  
  # Final dataset
  df <- data.frame(id = 1:n_total, S, A, x1, x2, u, time, status)
  
  # For target population, A, time, status are unobserved
  df$A[df$S == 0] <- NA
  df$time[df$S == 0] <- NA
  df$status[df$S == 0] <- NA
  
  return(list(df = df, truth = truth))
  
}

dr_balance <- function(S, X, Z1, time1, event1, cuts,
                       method = c("ebal", "optweight", "super"),
                       sl.lib = c("SL.mean", "SL.glm")) {
  
  X1 <- subset(X, S == 1)
  
  n <- nrow(X)
  n0 <- sum(S == 0)
  n1 <- sum(S == 1)
  
  # pseudo obersvations for survival
  lfit <- survfit(Surv(time1, event1) ~ 1)
  psi1.mat <- as.matrix(pseudo(lfit, times = cuts, type = "pstate"))
  
  fmla <- as.formula(paste("s ~", paste(colnames(X), collapse= "+")))
  
  if (method == "super") {
    wfit <- weightit(fmla, method = method, data = data.frame(s = S, X), estimand = "ATC",
                     SL.library = sl.lib, SL.method = "method.balance", criterion = "ks.max")
  } else {
    wfit <- weightit(fmla, method = method, data = data.frame(s = S, X), estimand = "ATC")
  }
  
  samp <- wfit$weights
  
  pseudo.out <- sapply(1:ncol(psi1.mat), function(i, ...) {
    
    # outcome model
    psi1 <- psi1.mat[,i]
    out <- SuperLearner(Y = psi1, X = data.frame(X1, Z = Z1),
                        SL.library = sl.lib, family = gaussian())
    
    # prediction over every observation
    mu_tmp <- cbind(c(predict(out, newdata = data.frame(X, Z = 0))$pred),
                    c(predict(out, newdata = data.frame(X, Z = 1))$pred))
    
    # Expand
    psi <- rep(0, n)
    psi[S == 1] <- psi1
    Z <- rep(0, n)
    Z[S == 1] <- Z1
    
    # weights multipliers
    prZ <- mean(Z1)
    w0 <- (1 - Z) * S * samp / (1 - prZ)
    w1 <- Z * S * samp / prZ
    
    mu <- cbind(Z*mu_tmp[,2] + (1 - Z)*mu_tmp[,1], mu_tmp)
    
    aug_est <- sum((w1 - w0) * (psi - mu[,1]))/n1 + mean(mu[S == 0,3] - mu[S == 0,2])
    eic <- c((w1 - w0) * (psi - mu[,1]))/mean(I(S == 1)) + c(I(S == 0)*(mu[,3] - mu[,2] - aug_est))/mean(I(S == 0))
    
    return(list(estimate = aug_est, variance = var(eic) / n))
    
  }, simplify = FALSE)
  
  estimate <- do.call(c, lapply(pseudo.out, function(item) item$estimate))
  variance <- do.call(c, lapply(pseudo.out, function(item) item$variance))
  names(estimate) <- names(variance) <- cuts
  
  return(list(estimate = estimate, variance = variance, wfit = wfit))
  
}


# --- [Your generate_data and dr_balance functions here] ---
# NOTE: The user's original functions are assumed to be loaded.

# --- Simulation Wrapper ---

# Define the scenarios to iterate over
scenarios <- c("baseline", "outcome_wrong", "selection_wrong", "both_wrong", 
               "positivity_practical", "positivity_structural")
n_sim <- 1000 # Number of iterations
n_total <- 1000 # Sample size for each dataset
cuts <- c(1, 2, 3, 4, 5, 6) # Time points for survival estimation

results_by_scenario <- vector("list", length(scenarios))

set.seed(42)

# Use lapply to run the simulation for each scenario
for (j in 1:length(scenarios)) {
  
  scenario_name <- scenarios[[j]]
  
  # Lists to store results from each iteration
  estimates_list <- vector("list", n_sim)
  variances_list <- vector("list", n_sim)
  truths_list <- vector("list", n_sim)
  
  # Main simulation loop for the current scenario
  for (i in 1:n_sim) {
    
    print(i)
    
    # 1. Generate Data
    sim_data <- generate_data(n_total = n_total, scenario = scenario_name, cuts = cuts)
    df <- sim_data$df
    truth_diff <- sim_data$truth$diff
    
    # 2. Prepare data for the estimator
    S <- df$S
    W <- data.frame(x1 = df$x1, x2 = df$x2)
    Z1 <- df$A[df$S == 1]
    time1 <- df$time[df$S == 1]
    event1 <- df$status[df$S == 1]
    
    # 3. Run the doubly robust estimator
    fit <- tryCatch({
      dr_balance(S = S, X = W, Z1 = Z1, time1 = time1, event1 = event1,
                 cuts = cuts, method = "optweight", sl.lib = c("SL.mean", "SL.glm"))
    }, error = function(e) {
      # Handle potential errors in weightit or SuperLearner
      list(estimate = rep(NA, length(cuts)), variance = rep(NA, length(cuts)))
    })
    
    # 4. Store results
    estimates_list[[i]] <- fit$estimate
    variances_list[[i]] <- fit$variance
    truths_list[[i]] <- truth_diff
    
  }
  
  # --- Process and Evaluate Results ---
  
  # Convert lists of vectors into matrices for easier column-wise calculations
  estimates_mat <- do.call(rbind, estimates_list)
  variances_mat <- do.call(rbind, variances_list)
  truths_mat <- do.call(rbind, truths_list)
  
  # Flag ridiculous outliers as NA (e.g., |estimate| > 10 or variance > 100)
  outlier_est <- abs(estimates_mat) > 10
  outlier_var <- variances_mat > 100
  estimates_mat[outlier_est] <- NA
  variances_mat[outlier_var] <- NA
  
  # Report how many iterations were flagged
  n_flagged <- sum(rowSums(outlier_est | outlier_var, na.rm = TRUE) > 0)
  if (n_flagged > 0) {
    message(sprintf("Scenario '%s': %d iterations flagged with outliers", scenario_name, n_flagged))
  }
  
  # Calculate the "true" parameter value by averaging over all simulations
  avg_truth <- colMeans(truths_mat, na.rm = TRUE)
  
  # Calculate Bias: Average Estimate - Average Truth
  bias <- colMeans(estimates_mat, na.rm = TRUE) - avg_truth
  
  # Calculate RMSE: sqrt(mean((estimate - truth)^2))
  rmse <- sqrt(colMeans((estimates_mat - truths_mat)^2, na.rm = TRUE))
  
  # Calculate Coverage Probability
  # Lower and upper bounds of the 95% confidence interval for each iteration
  lower_ci <- estimates_mat - 1.96 * sqrt(variances_mat)
  upper_ci <- estimates_mat + 1.96 * sqrt(variances_mat)
  
  # Check if the true value for that iteration falls within the CI
  coverage_indicators <- (matrix(rep(avg_truth, n_sim), byrow = T, nrow = n_sim) >= lower_ci) &
    (matrix(rep(avg_truth, n_sim), byrow = T, nrow = n_sim) <= upper_ci)
  
  # The coverage probability is the proportion of times the CI contained the truth
  coverage <- colMeans(coverage_indicators, na.rm = TRUE)
  
  # Return a summary data frame for the current scenario
  results_by_scenario[[j]] <- data.frame(
    scenario = scenario_name,
    time = cuts,
    truth = avg_truth,
    bias = bias,
    rmse = rmse,
    coverage = coverage
  )
  
}

# Combine results from all scenarios into a single data frame
final_results <- do.call(rbind, results_by_scenario)

# --- Display Final Results ---
print(final_results)

# Re-factor the 'scenario' column for a more logical plot legend order
final_results$scenario <- factor(final_results$scenario, 
                                 levels = c("baseline", 
                                            "outcome_wrong", 
                                            "selection_wrong",
                                            "both_wrong",
                                            "positivity_practical",
                                            "positivity_structural"),
                                 labels = c("Baseline",
                                            "Incorrect Outcome Model", 
                                            "Incorrect Selection Model", 
                                            "Exchangeability Violation",
                                            "Practical Positivity Violation",
                                            "Structural Positivity Violation"))

# 1. Reshape the data from wide to long format
results_long <- final_results %>%
  select(scenario, time, bias, coverage) %>%
  pivot_longer(
    cols = c("bias", "coverage"),
    names_to = "metric",
    values_to = "value"
  ) %>%
  # Capitalize the metric names for better plot labels
  mutate(metric = case_when(
    metric == "bias" ~ "Bias",
    metric == "coverage" ~ "Coverage Probability"
  ))

write_csv(results_long, "~/Documents/LEADER/simulation_table.csv")

# 2. Create a helper data frame for the horizontal reference lines
ref_lines <- data.frame(
  metric = c("Bias", "Coverage Probability"),
  hline = c(0, 0.95)
)

# 3. Create the faceted plot
faceted_plot <- ggplot(results_long, aes(x = time, y = value, color = scenario, group = scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  # Add facet_wrap to create separate panels for each metric
  facet_wrap(~ metric, scales = "free_y") +
  # Add the appropriate horizontal reference line to each facet
  geom_hline(data = ref_lines, aes(yintercept = hline), linetype = "dashed", color = "black") +
  labs(
    title = "Estimator Performance by Scenario",
    subtitle = "Bias and 95% Coverage Probability Across Follow-up Times",
    x = "Follow-up Time (t)",
    y = "Value", # A generic y-axis label, as the meaning is in the facet title
    color = "Scenario"
  ) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6)) +
  theme_bw() +
  theme(legend.position = "bottom")

# Display the final plot
png("~/Documents/LEADER/simulation.png", width = 1000, height = 500)
print(faceted_plot)
dev.off()