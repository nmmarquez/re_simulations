set.seed(123)
rm(list=ls())
require(lme4)
require(dplyr)

# This will go over a couple of different models using both the base lm function 
# as well as the lme4 random effects model


# the years that we model
Y <- 1980:2016 
# age group ids for model 
A <- 2:21
# locations taht we model
L <- 1:188 
# observation for each country age year
df <- expand.grid(location_id=L, age_group_id=A, year_id=Y)
N <- nrow(df)
# and well put a couple of covariates on there for good measure
K <- 2 # number of covariates
df[, paste0("x", 1:K)] <- sapply(1:K, function(x) rnorm(N))

# Example #1 simple linear effects model
# we will simulate a response with some error to show that models 
# will return the true value of beta 0 and 1
B0 <- 3 # this is the intercept term
B1 <- -6 # this is the slope on the variable x1
epsilon_sigma <- .4 # this is the standard deviation of the random error

df$y <- B0 + df$x1 * B1 + rnorm(N, sd=epsilon_sigma)

model1 <- lm(y ~ 1 + x1, data=df)
summary(model1)

# check out that the betas and the sd of the residual errror are pretty spot on

# Example #2 linear model with location specific fixed effects
loc_eff <- data.frame(location_id=L, loc_effect=runif(length(L), -10, 4))
df_loc <- dplyr::left_join(df, loc_eff)

df_loc$y <- B0 + df_loc$x1 * B1 + df_loc$loc_eff + rnorm(N, sd=epsilon_sigma)

model2 <- lm(y ~ 1 + x1 + as.factor(location_id), data=df_loc)
summary(model2)

# you will see the slope on x1 is good but know the model intercept shifted
head(loc_eff)

# example #3 linear model with random eff on location
# this will highlight how lme4 and codem treat and expect "random" effects 
loc_eff_sigma <- 1.2
loc_eff_rand <- data.frame(location_id=L, 
                           loc_effect=rnorm(length(L), sd=loc_eff_sigma))

df_rand <- dplyr::left_join(df, loc_eff_rand)

df_rand$y <- B0 + df_rand$x1 * B1 + df_rand$loc_eff + rnorm(N, sd=epsilon_sigma)

model3 <- lmer(y ~ 1 + x1 + (1|location_id), data=df_rand)
summary(model3)

# check out the random effects variance and where they come from