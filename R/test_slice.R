# Testing slice sampler

source("~/research/pitman_yor/battiston_sampler/R/b_sampler.R")

n <- 100
# Normal
for (i in 1:n) {
  s.samp[i] <- slice(50, dnorm, min=0, max=100, mean=50, sd=10, log=TRUE,
                     iter=25)
}

rnorm.samp <- rnorm(n, 50, 10)

qqplot(rnorm.samp, s.samp)
abline(0,1)

# Gamma
for (i in 1:n) {
  s.samp[i] <- slice(0.9, dgamma, min=0, max=10,  shape=1, rate=1, log=TRUE,
                     iter=25)
}

rgamma.samp <- rgamma(n, shape=1, rate=1)

qqplot(rgamma.samp, s.samp)
abline(0,1)
