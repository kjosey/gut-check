# Dependencies
library(parallel)
library(swdpwr)
library(geepack)
library(sandwich)
library(dplyr)
options(dplyr.summarise.inform = FALSE)

set.seed(42) # for reproducibility purposes

# data dimensions (subject to change)
n.school <- c(40, 60, 80) # number of schools
n.iter <- 1000 # number of simulations
n.semester <- 3 # number of semesters
n.student <- 1000 # average number of students

# parameters
p0 <- 0.8 # baseline effect
p1 <- 0.75*p0 # relative risk

sig2.seq <- seq(0.01, 1, by = 0.01) # sequence of possible overdispersions

# calculate icc
# icc.seq <- c(0, 0.05, 0.1, 0.15, 0.2)
icc.seq <- c(0, 0.025, 0.05, 0.075, 0.1)
icc.sig2.seq <- sapply(sig2.seq, function(z, ...) {
  iccCounts:::r_Pois(log(p0), log(sqrt(z)))
})

out_list <- mclapply(n.school, function(n, ...) {
  
  # school characteristics
  school_id <- rep(1:n, times = n.semester) # school
  area_id <- rep(rep(1:4, times = n/4), n.semester) # region
  treat_id <- rep(rep(c(0,1), times = c(n/2, n/2)), n.semester) # treatment
  
  # need to vary school size
  if (n.student < 800) {
    school_size <- rep(ceiling(rnorm(n, mean = 600, sd = 50)), n.semester) 
  } else if (n.student < 1000) {
    school_size <- rep(ceiling(rnorm(n, mean = 800, sd = 75)), n.semester)
  } else {
    school_size <- rep(ceiling(rnorm(n, mean = 1000, sd = 100)), n.semester)
  }
  
  # vary semesters
  if (n.semester == 2) {
    semester_id <- rep(c(1:2), each = n) # 2 semesters
  } else {
    semester_id <- rep(c(1:3), each = n) # 3 semesters
  }
  
  # combine data and expand to individual level
  dat_id <- data.frame(school = school_id,
                       area = area_id, 
                       treat = treat_id,
                       semester = semester_id)
  dat <- dat_id[rep(1:nrow(dat_id), times = school_size),]
  dat$area <- factor(dat$area)
  dat$semester <- factor(dat$semester)
  
  # coefficients
  beta <- c(log(p0), log(p1) - log(p0), rnorm((2*(nlevels(dat$area) - 1)), 0 , sd = 0.1)) # random effects for region modeled as fixed effects
  output <- data.frame() # initialize data frame
  
  for (icc in icc.seq) {
    
    # identify random intercept variance associated with ICC
    sig2 <- sig2.seq[which.min(abs(icc.sig2.seq - icc))]
    
    # loop through model fits
    test_gee <- sapply(1:n.iter, function(i, ...) {
      
      print(i)
      
      # simulate outcome
      X <- model.matrix(~ treat*area, data = dat)
      alpha <- rnorm(n, 0, sqrt(sig2)) # random intercept
      beta[3:(2*(nlevels(dat$area) - 1) + 2)] <- rnorm((2*(nlevels(dat$area) - 1)), 0, sd = 0.1)
      mu <- exp(alpha[dat$school] + c(X %*% beta))
      dat$y <- rpois(nrow(X), lambda = mu)
      
      # aggregate to school level
      dat_reduce <- dat %>% group_by(school, area, treat, semester) %>% 
        summarize(y = sum(y), size = n()) %>%
        arrange(school, area, treat, semester)
      
      # Model Fit
      fit <- geeglm(y ~ treat, family = poisson(link = "log"),
                    data = dat_reduce, offset = log(size), id = school,
                    corstr = "exchangeable", waves = semester)
      
      # Construct Power Results
      est <- coef(fit)[2]
      se <- sqrt(vcov(fit)[2,2])
      cut <- est + qnorm(0.95)*se
      return(as.numeric(cut < 0))
      
    })
    
    output <- rbind(output, data.frame(icc = round(icc, 3), pwr = mean(test_gee))) # power
    
  } 
  
  return(output) 
  
}, mc.cores = 3)

output <- cbind(n = rep(n.school, each = length(icc.seq)), do.call(rbind, out_list))

write.csv(output, file = "~/Documents/breathe_power_results.csv")
