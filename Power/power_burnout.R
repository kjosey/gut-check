library(lme4)
library(geepack)
library(sandwich)
library(clubSandwich)
library(merDeriv)
library(truncnorm)
library(margins)
library(devtools) 
library(MCPanel)

# base data
set.seed(42) # this makes the simulation exactly reproducible
n = 700 # number of persons
m = 2 # observation periods
prob = 0.7 # 1 - attrition rate

# fixed effects
hosp.b = c(-0.025, -0.05, -0.075, -0.1, -0.125)  # decrease over time for VA = 5% decrease
mi.b = seq(0.0005, 0.004, by = 0.0005) # 1 unit increase in MI increases burnout by 100*0.XX%
er.b = seq(0.02, 0.11, by = 0.01) # 1 unit increase in RE increases burnout by 100*0.XX%

did2_pwr <- function(n, m, beta, gamma, exposure = c("mi", "er")) {

  id = factor(rep(1:n, each = m))
  
  # covariates
  t = rep(1:m, times = n) - 1 # time
  w = rep(c(1, 0), each = m*(n/2)) # hospital
  alpha = rep(rnorm(n, 0.25, 0.04), each = m) # random intercept
  
  if (exposure == "er") {
    x = as.vector(replicate(n, rtruncnorm(m, mean = 1.3, sd = 0.4, a = 1e-6)))
    x = c(scale(x, center = TRUE, scale = FALSE))
  } else if (exposure == "mi") {
    x = as.vector(replicate(n, rtruncnorm(m, mean = 37, sd = 13, a = 10, b = 100)))
    x = c(scale(x, center = TRUE, scale = FALSE))
  }
  
  # linear probabilities
  eta = alpha + 0.02*t + x*beta + w*I(t > 0)*gamma
  eta[eta < 0] <- 1e-6
  eta[eta > 1] <- 1 - 1e-6
  
  # outcomes
  y = rbinom(n*m, 1, prob = eta)
  
  # examine marginals
  # fit <- glm(y ~ er*t + hosp*t, data = data.frame(t = t, er = er, y = y),
  #              family = binomial(link = "identity"))
  # fit.2 <- glm(y.2 ~ mi*t + hosp*t, data = data.frame(t = t, mi = mi, y.2 = y.2),
  #              family = binomial(link = "identity"))
  
  # MCAR
  # z.mcar = matrix(y, ncol = m, nrow = n, byrow = TRUE)
  # 
  # for(i in 1:n){
  #   
  #   for(j in 2:m){
  #     
  #     if (is.na(z.mcar[i,j - 1])){ 
  #       
  #       z.mcar[i,j] <- NA
  #       
  #     } else {
  #       
  #       z.mcar[i,j] <- ifelse(runif(1) > prob, NA, z.mcar[i,j])
  #       
  #     }
  #   }
  # }
  # 
  # y.mcar = as.vector(t(z.mcar))

  # MAR
  z.mar = matrix(y, ncol = m, nrow = n, byrow = TRUE)

  for(i in 1:n){

    for(j in 2:m){

      if (is.na(z.mar[i,j - 1])){

        z.mar[i,j] <- NA

      } else if (z.mar[i,j - 1] > 0) {

        z.mar[i,j] <- ifelse(runif(1) > 0.8*prob, NA, z.mar[i,j])

      } else {
        
        z.mar[i,j] <- ifelse(runif(1) > prob, NA, z.mar[i,j])
        
      }
    }
  }
  
  y.mar = as.vector(t(z.mar))

  # NMAR
  # z.nmar = matrix(y, ncol=m, nrow=n, byrow=TRUE)
  # 
  # for(i in 1:n){
  #   
  #   for(j in 2:m){
  #     
  #     dif = y.nmar[i,j]-y.nmar[i,j-1]
  #     
  #     if (dif > 0) {  # if burnout stays high twice, drops out
  #
  #       z.nmar[i,j] = ifelse(runif(1) > prob, NA, z.nmar[i,j])
  #
  #     }
  #     
  #   }
  # }
  # 
  # y.nmar = as.vector(t(z.nmar))

  # Imputation
  data.full <- data.frame(t0 = as.numeric(t > 0), t = t, w = w, x = x, id = id, y.mar = y.mar)
  data <- data.full[complete.cases(data.full),]
  # m.data <- mice(data.full, method = "pmm", m = 1)
  # data <- complete(m.data, action = 1)
  
  # GEE
  # fit <- geeglm(y.mcar ~ x + t + w:t0, family = binomial(link = "identity"),
  #               data = data, id = id, waves = t, corstr = "exchangeable")
  # Sig <- vcov(fit)
  # se <- sqrt(diag(Sig))
  # coefs <- fit$coefficients
  # ci <- cbind(coefs - 1.96*se, coefs + 1.96*se) 
  
  # Mixed Model
  fit <- glmer(y.mar ~ x + t + w:t0 + (1|id), family = binomial(), data = data)
  # sand <- try(sandwich(fit, bread = bread(fit, full = TRUE), mean = meat(fit, level = 2)))
  
  # if (inherits(sand, "try-error"))
  #   Sig <- vcov(fit)
  # else
  #   Sig <- sand[1:4, 1:4]
  
  Sig <- vcov(fit)
  se <- sqrt(diag(Sig))
  coefs <- fit@beta
  ci <- cbind(coefs - 1.96*se, coefs + 1.96*se)
  pwr <- as.numeric(ci[,1] > 0 | ci[,2] < 0)
  
  return(pwr)
  
}

val.mi <- sapply(mi.b, function(b, ...) {
  
  print(b)
  pwrmat <- replicate(1000, did2_pwr(n = n, m = m, beta = b, gamma = hosp.b[2], exposure = "mi"))
  rowMeans(pwrmat)
  
})

val.er <- sapply(er.b, function(b, ...) {
  
  print(b)
  pwrmat <- replicate(1000, did2_pwr(n = n, m = m, beta = b, gamma = hosp.b[2], exposure = "er"))
  rowMeans(pwrmat)
  
})

# val.hosp <- sapply(hosp.b, function(c, ...) {
#   
#   print(c)
#   pwrmat <- replicate(1000, did2_pwr(n = n, m = m, gamma = c, beta = er.b[3],  exposure = "er"))
#   rowMeans(pwrmat)
#   
# })

colnames(val.mi) <- mi.b
colnames(val.er) <- er.b
# colnames(val.hosp) <- hosp.b

write_csv(data.frame(t(val.mi)), file = "~/Documents/mi.csv")
write_csv(data.frame(t(val.er)), file = "~/Documents/er.csv")
# val.hosp
