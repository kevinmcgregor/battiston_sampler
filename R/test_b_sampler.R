# Testing Battiston sampler
source("~/research/pitman_yor/battiston_sampler/R/b_sampler.R")
source("~/research/pitman_yor/basic_py/R/hpy_functions.R")

# Sampling from HPY
N <- 1000
J <- 10
conc.top <- 1
dsct.top <- 0.2
conc.local <- rep(5, J)
dsct.local <- rep(0.3, J)
s.hpy <- sample_hpy(J, rep(N, J), conc.top, dsct.top, conc.local, dsct.local)
Y <- s.hpy$abund

# Running Battiston sampler
s <- b_sampler(Y, 1000, 100, p.shape = 1, p.scale = 1, mc.cores=10)


# Testing individual functions for sampler for concentration parameter ####
set.seed(24601)
J <- 10
conc <- 0.5
dsct <- 0.5
p.shape <- 2
p.scale <- 10
N <- rep(1000, J)
Y <- matrix(100, J, 10)
n.tab <- rep(100, J)
n.s.tab <- matrix(10, J, 10)
q <- rbeta(J, conc, N)
Q <- 1/p.scale - sum(log(q))

mc <- map_conc(conc, p.shape, Q, n.tab, dsct)
#s.samp <- slice_conc(mc, p.shape, Q, n.tab, dsct)

xv <- seq(max(0.001, mc-mc/2), mc+mc/2, by=0.001)
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


# Testing the sampler for the discount parameter
xv <- seq(0.01, 0.99, by=0.001)
yv <- rep(0, length(xv))
for (i in 1:length(xv)) {
  yv[i] <- prob_dsct(xv[i], conc, n.tab, n.s.tab, Y)
}

n.samp <- 1000
dsct.samp <- rep(0, n.samp)
for (i in 1:n.samp) {
  dsct.samp[i] <- samp_dsct(dsct, conc, n.tab, n.s.tab, Y)
}

# Scaling curve to plot over histogram
sc <- 1
h.vals <- hist(dsct.samp, breaks=20)
max.hval <- max(h.vals$counts)
lines(xv, sc*max.hval*exp(yv-max(yv)))
