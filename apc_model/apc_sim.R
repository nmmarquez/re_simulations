set.seed(123)
rm(list=ls())
library(dplyr)

correlated_error <-function(N, sd=.1){
    # correlated errors froma random normal distribution ala a RW
    errors <- rep(0, N)
    errors[1] <- rnorm(1, sd=sd)
    for(i in 2:N){
        errors[i] <- rnorm(1, errors[i-1], sd=sd)
    }
    errors
}

# age specific mort rates
ages <- seq(35, 80)
A <- length(ages)
age_mort <- log(1:A) * 6 + 35
age_mort[ages > 47] <- age_mort[which(ages==48)]
plot(ages, age_mort)

# rate ratio by year setting the 12 year as the reference year
years <- 1960:2010
Y <- length(years)
year_rate_ratio <- (1:Y - 3) * -.02 + correlated_error(Y, sd=.03) + 1
year_rate_ratio <- year_rate_ratio / year_rate_ratio[12]

# rate ratio by cohort setting the 15 cohort as the reference
cohorts <- 1920:1975
C <- length(cohorts)

cohort_ratio <- (.0022 * (1:C - 15)**2 + .0022 * 1:C + 
                     correlated_error(C, .05) + 1)
cohort_ratio <- cohort_ratio / cohort_ratio[15]
plot(cohorts, cohort_ratio, type="l", xlim=c(min(cohorts), max(years)), 
     ylim=c(min(year_rate_ratio), max(cohort_ratio)), ylab="rate ratio", 
     xlab="year")
lines(years, year_rate_ratio, col=2)
legend("topright", legend=c("Cohort","Period"), lty=c(1,1), 
       lwd=c(.5,.5), col=c("black","red"))

# make sure ratios are above zero
all(c(cohort_ratio > 0, year_rate_ratio > 0))

# create the rates from the paramters
df <- data.frame(mort_rate=c(cohort_ratio %*% t(age_mort)))
df$cohort <- rep(cohorts, A)
df$age <- rep(ages, each=C)
df$year <- df$cohort + df$age

year_df <- data.frame(year=years, rate_ratio=year_rate_ratio)
df <- left_join(df, year_df)
df$mort_rate <- df$mort_rate * df$rate_ratio
df$rate_ratio <- NULL

# we can remove data where we dont have an obeservation for that APC combination
df <- df[!is.na(df$mort_rate),]
df$death_count <- rpois(nrow(df), df$mort_rate)

# make sure everything is a factor for the regression
for(variable in c("age", "year", "cohort")){
    df[,variable] <- as.factor(as.integer(df[,variable]))
}

# change the ref group
df <- within(df, year <- relevel(year, ref = 12))
df <- within(df, cohort <- relevel(cohort, ref = 15))

# run the model with all paramters
glm1 <- glm(death_count ~ age + year + cohort, data=df, family=poisson)
# the model does not converge because of classical APC model identifiebility
summary(glm1)

# lets implement the hack of setting two age groups the same like the paper
df_hack <- df
df_hack[df_hack$age == "50", "age"] <- "70"
glm2 <- glm(death_count ~ age + year + cohort, data=df_hack, family=poisson)
summary(glm2)

# oh crap it kinda works!!! still need to do some vetting though

# To do
# 1) extract the model age mort rates and compare them to actual mort rates
# 2) do the same for year and cohort rate ratios
# 3) change the true values to have parallel decreasing cohort and year effects
