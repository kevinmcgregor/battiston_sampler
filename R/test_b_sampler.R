# Testing Battiston sampler
source("~/research/pitman_yor/battiston_sampler/R/b_sampler.R")

# Testing sampler for concentration parameter ####
set.seed(24601)
J <- 25
conc <- 10
dsct <- 0.8
p.shape <- 1
p.scale <- 1
N <- rep(1000, J)
n.tab <- rep(50, J)
q <- rbeta(J, conc, N)
Q <- 1/p.scale - sum(log(q))

mc <- map_conc(conc, p.shape, Q, n.tab, dsct)
s.samp <- slice_conc(mc, p.shape, Q, n.tab, dsct)

xv <- seq(0.1, 5, by=0.1)
yv <- rep(0, length(xv))
for (i in 1:length(xv)) {
  yv[i] <- prob_conc(xv[i], p.shape, Q, n.tab, dsct)
}
plot(xv, yv, type="l")
abline(v=mc, col="red", lwd=2, lty=2)


