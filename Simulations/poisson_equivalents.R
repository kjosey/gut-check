library(dplyr)

n <- 10000
x1 <- rbinom(n, 1, 0.2)
x2 <- rbinom(n, 1, 0.4)
x3 <- rbinom(n, 1, 0.6)
x4 <- rbinom(n, 1, 0.8)
zip <- sample(1:500, n, replace = T)
ind_data <- data.frame(zip, x1, x2, x3, x4)

w1 <- rnorm(500)
w2 <- rnorm(500)
zip_data <- data.frame(zip = 1:500, w1 = w1, w2 = w2)
data <- merge(ind_data, zip_data, by = "zip")

lambda <- with(data, exp(-3 + 0.5*x1 - x2 + x3 - 0.5*x4 + 0.1*w1 - 0.1*w2))
data$y <- rbinom(n, 1, lambda)

agg_data <- data %>% 
  group_by(zip, x1, x2, x3, x4) %>% 
  summarise(y = sum(y), offset = n())

new_data <- merge(agg_data, zip_data, by = "zip")
new_data$ybar <- with(new_data, y/offset)

target <- glm(y ~ x1 + x2 + x3 + x4 + w1 + w2, data = data, family = quasipoisson(link = "log"))            
test1 <- glm(y ~ x1 + x2 + x3 + x4 + w1 + w2, offset = log(offset), data = new_data, family = quasipoisson(link = "log"))  
test2 <- glm(ybar ~ x1 + x2 + x3 + x4 + w1 + w2, weights = offset, data = new_data, family = quasipoisson(link = "log"))  
test3 <- glm(ybar ~ x1 + x2 + x3 + x4 + w1 + w2, data = new_data, family = quasipoisson(link = "log")) 

summary(target) # what we estimate
summary(test1) # with an offset
summary(test2) # with weights parameterized by rate
summary(test3) # rate without weights
