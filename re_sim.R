rm(list=ls())

# probability of dying of diptheria
p_no_vac <- .47

# probability of dying of diptheria if you are vacinated 
p_yes_vac <- .28


# simulate data from 200 repondents 100 vacintaed 100 not.
N <- 100
pop_no_vac <- rbinom(N, size=1, prob=p_no_vac)
pop_yes_vac <- rbinom(N, size=1, prob=p_yes_vac)

df <- data.frame(death=c(pop_no_vac, pop_yes_vac), vac= c(rep(0, N), rep(1, N)))

model_int <- glm(death ~ 1, data=df, family="binomial")
exp(model_int$coefficients)
