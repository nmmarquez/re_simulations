rm(list=ls())

# etsimate some noisy sin data
N <- 1000
x <- 1:N / 10
obs <- sin(x) + rnorm(N, sd=.1)
plot(x, obs, type="l")



# Lets now plot the diff of the data
plot(obs, type="l")
plot(diff(obs), type="l")
plot(diff(obs, differences=2), type="l")
plot(diff(obs, differences=3), type="l")

# and the ACF
acf(obs)
acf(diff(obs, differences=1))
acf(diff(obs, differences=2))
acf(diff(obs, differences=3))

# now when we take the seasonal lag of approx 2*pi*10
acf(diff(obs, differences=1, lag = 63))


# what does dicky fuller say about this??
tseries::adf.test(obs)

# but whats about Kwiatkowski-Phillips-Schmidt-Shin????
tseries::kpss.test(obs)

# seems like these are only good for catching the non seasonal stationarity

# but what about the ndiffs func?
forecast::ndiffs(obs)

# still a fail but check out nsdiffs for number of seasonal diffs

obs


def

