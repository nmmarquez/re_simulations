rm(list=ls())
pacman::p_load(greta, DiagrammeR)

x <- as_data(iris$Petal.Length)
y <- as_data(iris$Sepal.Length)

# variables and priors
int = normal(0, 1)
coef = normal(0, 3)
sd = student(3, 0, 1, truncation = c(0, Inf))

# operations
mean <- int + coef * x

# likelihood
distribution(y) = normal(mean, sd)

# defining the model
m <- model(int, coef, sd)

# plotting
plot(m)

# sampling
draws <- mcmc(m, n_samples = 1000)
plot(draws)

lm1 <- lm(Sepal.Length ~ Petal.Length, data=iris)
mvtnorm::rmvnorm(10, lm1$coefficients, vcov(lm1))

