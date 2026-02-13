library(geepack)
library(lme4)
library(dplyr)
library(swdpwr)
library(MASS)
library(tidyr)
library(readr)

###############################################################################
## POWER ANALYSIS: 2-ARM (25 mg vs. 5 mg Psilocybin)
## Outcome: MADRS (0-60)
## Effect estimates: Derived from EPiSODE (HAMD17) -> mapped to MADRS scale
##
## HAMD17 -> MADRS Mapping
## -----------------------------------------------------------------------
## Goodwin et al. (NEJM 2022) reports BOTH scales at baseline in the same
## TRD patients (N = 233):
##   Overall:  MADRS 32.5 (SD 6.0)  |  HAMD17 22.2 (SD 2.9)
##   25 mg:    MADRS 31.9 (SD 5.4)  |  HAMD17 21.8 (SD 3.0)
##   10 mg:    MADRS 33.0 (SD 6.3)  |  HAMD17 22.4 (SD 2.8)
##    1 mg:    MADRS 32.7 (SD 6.2)  |  HAMD17 22.2 (SD 2.9)
##
##   Mean ratio:  32.5 / 22.2 = 1.464
##
## Goodwin's own power calculation assumed a MADRS change-score SD of 11.0,
## which serves as an external anchor for the mapped variance parameters.
##
## EPiSODE (HAMD17) Treatment Phase 1, Week 6 (while-on-treatment; Table 2)
## -----------------------------------------------------------------------
##   25 mg (N=47): BL 22 (SD~5.1), change -6.06 (SD 7.93)
##    5 mg (N=48): BL 22 (SD~4.7), change -3.46 (SD 5.82)
##   Placebo (N=47): BL 22 (SD~3.6), change -1.51 (SD 5.36)
##   Model-estimated 25mg vs 5mg: -3.09 (95% CI: -5.50, -0.69)
##
## EPiSODE male subgroup (eTable 11, N=83; relevant for veteran cohort)
##   25 mg males (n=26): BL 21.96 (SD 5.54), Wk 6 change -7.85 (SD 8.48)
##    5 mg males (n=30): BL 20.93 (SD 3.96), Wk 6 change -4.03 (SD 5.80)
##
## Mapped MADRS parameters (x 1.464)
## -----------------------------------------------------------------------
##   Baseline MADRS mean:  22 x 1.464 ~ 32  (consistent with Goodwin)
##
##   5 mg reference change:       -3.46 x 1.464 ~ -5.1
##   25 mg change:                -6.06 x 1.464 ~ -8.9
##   Model DiD (25mg vs 5mg):     -3.09 x 1.464 ~ -4.5
##
##   Change-score SD (25 mg):  7.93 x 1.464 ~ 11.6  [cf. Goodwin's 11.0]
##   Change-score SD (5 mg):   5.82 x 1.464 ~  8.5
##
##   Male subgroup (x 1.464):
##     5 mg change:  -4.03 x 1.464 ~ -5.9
##     25 mg change: -7.85 x 1.464 ~ -11.5
##     Change-score SD (25 mg males): 8.48 x 1.464 ~ 12.4
##     Change-score SD (5 mg males):  5.80 x 1.464 ~ 8.5
##
## Baseline MADRS SDs (from Goodwin overall ~ 6.0; Raison: 5.7 and 4.5)
##   We use 6.0 as a common baseline SD across arms.
##
## Sensitivity parameters
## -----------------------------------------------------------------------
##   rho (within-subject correlation): 0.3, 0.5, 0.7
##     - No trial in the project reports the within-subject correlation
##     - Crude back-calculation from EPiSODE observed SDs yields rho ~ 0.1
##       (25mg arm) to 0.6 (5mg arm), supporting a broad sensitivity range
##   sig.diff.mult (change-score SD multiplier): 0.85, 1.0, 1.15
##     - Captures uncertainty in the HAMD17 -> MADRS variance mapping
##     - At mult = 1.0: 25mg SD = 11.6, 5mg SD = 8.5
##     - At mult = 0.85: consistent with Goodwin's assumed SD ~ 11.0
##
## Parametrization
## -----------------------------------------------------------------------
##   mu = 5 mg reference absolute change (MADRS)
##   lambda = scaling factor for 25 mg relative to 5 mg
##   DiD (interaction) = (lambda - 1) x mu
##   lambda = 1     -> null (no difference between arms)
##   lambda ~ 1.75  -> MADRS DiD ~ -3.8 (overall estimate)
##   lambda ~ 1.95  -> MADRS DiD ~ -4.8 (male subgroup estimate)
##
## Raison et al. (JAMA 2023) MADRS parameters (for context)
## -----------------------------------------------------------------------
##   Psilocybin (n=51): BL MADRS 35.5 (SD 5.7), Day 43 change -19.1
##   Niacin (n=53):     BL MADRS 35.0 (SD 4.5), Day 43 change  -6.8
##   Difference: -12.3 (95% CI: -17.5 to -7.2)
##   [Much larger effect -- MDD (not TRD), active placebo only]
###############################################################################

### FUNCTIONS FOR COVARIANCE

# Solve for the post-treatment marginal SD from change-score SD, baseline SD,
# and assumed within-subject correlation.
# Identity: Var(Y1 - Y0) = sig1^2 + sig0^2 - 2*rho*sig0*sig1 = sig.diff^2
# Rearranges to a quadratic in sig1.
marg_var <- function(sig.diff, sig0, rho) {
  
  a <- 1
  b <- -2*rho*sig0
  c <- sig0^2 - sig.diff^2
  
  disc <- b^2 - 4*a*c
  if (disc < 0) stop("No real solution: change-score SD too small given rho and sig0")
  
  roots <- c((-b - sqrt(disc))/(2*a),
             (-b + sqrt(disc))/(2*a))
  
  root <- roots[which(roots > 0)]
  if (length(root) == 0) stop("No positive root found")
  if (length(root) > 1) root <- min(root)  # take smaller positive root
  return(root)
  
}

# Exchangeable covariance pattern
xch <- function(sig2, rho, p){
  
  if(length(sig2) != 1)
    stop("length(sig2) != 1")
  
  if(length(rho) != 1)
    stop("length(rho) != 1")
  
  R <- matrix(rho, p, p)
  diag(R) <- 1
  D <- diag(sqrt(sig2), p, p)
  V <- D %*% R %*% D
  return(V)
  
}

# AR(1) covariance pattern
ar1 <- function(sig2, phi, p){
  
  if(length(sig2) != 1)
    stop("length(sig2) != 1")
  
  if(length(phi) != 1)
    stop("length(phi) != 1")
  
  if (length(p) == 1)
    times <- 1:p
  else
    times <- p
  
  H <- abs(outer(times, times, "-"))
  V <- sig2 * phi^H
  return(V)
  
}

# Unstructured covariance pattern
un <- function(sig2, R){
  
  if (!isSymmetric(R))
    stop("R must be a square/symmetric matrix")
  
  if( any(eigen(R)$values <= 0) )
    stop("R must be positive definite")
  
  if (any(diag(R) != 1))
    stop("R must have 1 along the diagonals")
  
  if (any(abs(R) > 1))
    stop("off diagonal entries of R must be between (-1,1)")
  
  if (length(sig2) == 1)
    sig2 <- rep(sig2, times = nrow(R))
  else if (length(sig2) != nrow(R))
    stop("length(sig2) != nrow(R)")
  
  if (any(sig2 <= 0))
    stop("sig2 must be strictly positive")
  
  p <- nrow(R)
  D <- diag(sqrt(sig2), p, p)
  V <- D %*% R %*% D
  
  return(V)
  
}

### SIMULATIONS

## Parameter values -- MADRS scale, EPiSODE-derived
alpha <- 0.05
n <- 78
days_seq <- c(0, 42) # repeated measure days (baseline, Week 6)
trt_seq <- c("Trt5", "Trt25") # treatment groups
p <- length(trt_seq)

# lambda: scaling factor for 25 mg relative to 5 mg
#   lambda = 1     -> no difference between arms (null)
#   lambda ~ 1.75  -> MADRS DiD ~ -3.8 (EPiSODE overall mapped estimate)
#   lambda ~ 1.95  -> MADRS DiD ~ -4.8 (EPiSODE male mapped estimate)
lambda_seq <- c(1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5)
rho_seq <- c(0.3, 0.5, 0.7) # within-subject correlation

# mu: reference group (5 mg) absolute change from baseline on MADRS
#   EPiSODE overall: -3.46 x 1.464 ~ -5.1
mu_seq <- c(-5.0)

# Change-score SD sensitivity
#   Base values (MADRS, mapped from EPiSODE x 1.464):
#     25 mg: 7.93 x 1.464 ~ 11.6
#      5 mg: 5.82 x 1.464 ~  8.5
#   Multiplier captures mapping uncertainty
sig.diff.mult_seq <- c(0.85, 1.0, 1.15)

# Baseline MADRS SD (common across arms; from Goodwin overall)
sig0 <- 6.0

# Base change-score SDs (MADRS, mapped from EPiSODE)
sig.diff.base <- c(Trt25 = 11.6, Trt5 = 8.5)

## Simulation scenarios
n.iter <- 1000
scenarios <- expand.grid(lambda = lambda_seq, rho = rho_seq, mu = mu_seq,
                         sig.diff.mult = sig.diff.mult_seq)

results <- data.frame()

# Start simulation
for (i in 1:nrow(scenarios)) {
  
  print(paste0("Scenario ", i, " / ", nrow(scenarios)))
  
  lambda <- scenarios$lambda[i]
  rho <- scenarios$rho[i]
  mu <- scenarios$mu[i]
  sig.diff.mult <- scenarios$sig.diff.mult[i]
  
  # Scale change-score SDs by multiplier
  sig.diff.25 <- sig.diff.base["Trt25"] * sig.diff.mult
  sig.diff.5  <- sig.diff.base["Trt5"]  * sig.diff.mult
  
  ## Mean matrices
  trt <- rep(trt_seq, each = floor(n/p))
  bl <- matrix(32, nrow = n, ncol = length(days_seq)) # MADRS baseline ~ 32
  
  trt.mult0 <- matrix(rep(c(0, mu), times = n),
                       nrow = n, byrow = TRUE)
  
  mu.mat0 <- bl + trt.mult0       # for 5 mg (reference)
  mu.mat1 <- bl + lambda*trt.mult0 # for 25 mg
  
  ## Covariance matrices -- MADRS scale
  # Derive post-treatment marginal SD from change-score SD, baseline SD,
  # and assumed within-subject correlation
  sig1_diag <- c(sig0, marg_var(sig.diff.25, rho = rho, sig0 = sig0)) # 25 mg
  sig0_diag <- c(sig0, marg_var(sig.diff.5,  rho = rho, sig0 = sig0)) # 5 mg
  
  R <- xch(sig2 = 1, rho = rho, p = length(days_seq))
  Sigma1 <- un(sig2 = sig1_diag^2, R = R)
  Sigma0 <- un(sig2 = sig0_diag^2, R = R)
  
  test.stat <- rep(NA, n.iter) # initialize rejection
  
  for (j in 1:n.iter) {
    
    # Generate Data, Heteroskedasticity with Treatment
    y1_wide <- mu.mat1 + mvrnorm(n, rep(0, length(days_seq)), Sigma1)
    y0_wide <- mu.mat0 + mvrnorm(n, rep(0, length(days_seq)), Sigma0)
    y_wide <- y1_wide*as.numeric(trt == "Trt25") + 
      y0_wide*as.numeric(trt == "Trt5")
    
    # Restrict outcomes to MADRS range [0, 60]
    y_wide[y_wide < 0] <- 0
    y_wide[y_wide > 60] <- 60
    colnames(y_wide) <- days_seq
    
    # Wide to Long
    data <- data.frame(y_wide) %>% 
      pivot_longer(
        cols = 1:length(days_seq), 
        names_to = "time",
        values_to = "y")
    
    data$id <- factor(rep(1:nrow(y_wide), each = length(days_seq)))
    data$time <- factor(data$time, levels = paste0("X", days_seq), label = days_seq)
    data$trt <- factor(rep(trt, each = length(days_seq)), levels = c("Trt5","Trt25")) 
    mat_data <- data.frame(y = data$y, id = data$id, time = data$time, model.matrix(~ trt*time, data = data)[,-1])
    
    # GEE
    mod_gee <- geeglm(y ~ . - time - id, data = mat_data, id = id,
                       waves = time, corstr = "exchangeable")
    
    # Covariance Estimates
    vcov_gee = mod_gee$geese$vbeta
    
    ## Setup contrast vectors
    contr_bl <- rep(0, length(coef(mod_gee)))
    names(contr_bl) <- names(coef(mod_gee))
    
    # DiD hypothesis: treatment x time42 interaction
    contr <- as.matrix(rbind(
      data.frame(t(contr_bl)) %>%
        mutate(trtTrt25.time42 = 1)))
    
    # Point estimates
    gee_est = c((contr) %*% coef(mod_gee))
    gee_se = sqrt(diag(contr %*% vcov_gee %*% t(contr)))
    
    # Test statistics
    z_vals <- c(Z_trt25 = gee_est[1]/gee_se[1])
    
    pvals <- 2*pnorm(abs(z_vals), lower.tail = FALSE) # two-sided
    test.stat[j] <- as.numeric(pvals < alpha) 
    
  }
  
  results <- rbind(results, data.frame(n = n,
                                       baseline_mean = 32,
                                       baseline_change = mu,
                                       correlation = rho,
                                       sig_diff_mult = sig.diff.mult,
                                       sig_diff_25 = as.numeric(sig.diff.25),
                                       sig_diff_5 = as.numeric(sig.diff.5),
                                       rel_diff = lambda,
                                       abs_diff = (lambda - 1)*mu, 
                                       hypothesis = c("TRT25 VS. TRT5"),
                                       pwr = mean(test.stat, na.rm = T)))
  
}

write_csv(results, "~/Documents/VA/Psilocybin/DiD_GEE_2arm_MADRS.csv")
