# Testing Battiston sampler
source("~/research/pitman_yor/battiston_sampler/R/b_sampler.R")

# Testing individual functions for sampler for concentration parameter ####
set.seed(24601)
J <- 100
conc <- 1
dsct <- 0.5
p.shape <- 2
p.scale <- 10
N <- rep(1000, J)
n.tab <- rep(400, J)
q <- rbeta(J, conc, N)
Q <- 1/p.scale - sum(log(q))

mc <- map_conc(conc, p.shape, Q, n.tab, dsct)
#s.samp <- slice_conc(mc, p.shape, Q, n.tab, dsct)

xv <- seq(max(0.001, mc-mc/2), mc+mc/2, by=0.01)
yv <- rep(0, length(xv))
for (i in 1:length(xv)) {
  yv[i] <- prob_conc(xv[i], p.shape, Q, n.tab, dsct)
}
#plot(xv, exp(yv-max(yv)), type="l")
#abline(v=mc, col="red", lwd=2, lty=2)

# Testing the final sampler for the concentration parameter
n.samp <- 1000
conc.samp <- rep(0, n.samp)
for (i in 1:n.samp) {
  conc.samp[i] <- samp_conc(conc, p.shape, p.scale, n.tab, N, dsct)
}

# Scaling curve to plot over histogram
h.vals <- hist(conc.samp, breaks=20)
max.hval <- max(h.vals$counts)
lines(xv, max.hval*exp(yv-max(yv)))
abline(v=mc, col="red", lwd=2, lty=2)



