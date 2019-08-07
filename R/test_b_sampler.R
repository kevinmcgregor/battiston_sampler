# Testing Battiston sampler
source("~/research/pitman_yor/battiston_sampler/R/b_sampler.R")
source("~/research/pitman_yor/basic_py/R/hpy_functions.R")

plot.dir <- "/mnt/GREENWOOD_BACKUP/home/kevin.mcgregor/research/pitman_yor/battiston_sampler/plots/fixed_params"

# Sampling from HPY
N <- 1000
J <- 25
conc.top <- 3
dsct.top <- 0.5
conc.local <- rep(5, J)
dsct.local <- rep(0.3, J)
s.hpy <- sample_hpy(J, rep(N, J), conc.top, dsct.top, conc.local, dsct.local)

Y <- s.hpy$abund

# Output csv to try on Chris's C code
output.csv <- TRUE
o.filename <- "/mnt/GREENWOOD_BACKUP/home/kevin.mcgregor/research/pitman_yor/chris_todd_sampler/code/test/in_data.csv"
if (output.csv) {
  Y.t <- t(Y)
  rownames(Y.t) <- paste0("otu", 1:NROW(Y.t))
  colnames(Y.t) <- paste0("samp", 1:NCOL(Y.t))
  write.csv(Y.t, file=o.filename, quote=FALSE)
}

# Y <- matrix(0, NROW(s.hpy$abund), NCOL(s.hpy$abund))
# for (i in 1:J) {
#   Y[i,] <- t(rmultinom(1, N, s.hpy$abund[i,]))
# }

# Prior hyper-params
p.shape = 1 #1.1
p.scale = NCOL(Y)

# Running Battiston sampler
s <- b_sampler(Y, 1000, 200, p.shape = p.shape, p.scale = p.scale)

# Running parallel chains
library(parallel)
n.chain <- 4
Y.list <- vector("list", n.chain)
for (i in 1:n.chain) Y.list[[i]] <- Y

s.p <- mclapply(Y.list, b_sampler, p.shape = p.shape, p.scale = p.scale,
                n.iter=1000, n.burn=200, mc.cores=1)#n.chain)

# Checking traceplot of tables
par(mfrow=c(5,2), mar=rep(1,4))
for (j in 1:J) {
  y.rge <- range(s$n.tab[,j], s.hpy$t.tab[j])
  plot(s$n.tab[,j], type="l", ylim=y.rge)
  abline(h=s.hpy$t.tab[j], col="red", lty=2, lwd=2)
}

est.tab <- colMeans(s$n.tab[100:NROW(s$n.tab),])
plot(s.hpy$t.tab, est.tab); abline(0,1)

# Checking traceplot for top-level concentration param
par(mfrow=c(1,1), mar=rep(2.5, 4))
y.rge <- range(s$gamma, conc.top)
plot(s$gamma, type="l", ylim=y.rge)
abline(h=conc.top, col="red", lty=2, lwd=2)

# Checking traceplot for top-level discount param
par(mfrow=c(1,1), mar=rep(2.5, 4))
y.rge <- range(s$alpha, dsct.top)
plot(s$alpha, type="l", ylim=y.rge)
abline(h=dsct.top, col="red", lty=2, lwd=2)

# Checking traceplots of local-level concentration params
par(mfrow=c(5,2), mar=rep(1,4))
for (j in 1:J) {
  y.rge <- range(s$theta[,j], conc.local[j])
  plot(s$theta[,j], type="l", ylim=y.rge)
  abline(h=conc.local[j], col="red", lty=2, lwd=2)
}

# Checking traceplots of local-level discount params
par(mfrow=c(5,2), mar=rep(1,4))
for (j in 1:J) {
  y.rge <- range(s$sigma[,j], dsct.local[j])
  plot(s$sigma[,j], type="l", ylim=y.rge)
  abline(h=dsct.local[j], col="red", lty=2, lwd=2)
}

# Posterior means
gamma.post <- mean(s$gamma)
alpha.post <- mean(s$alpha)
theta.post <- colMeans(s$theta)
sigma.post <- colMeans(s$sigma)
n.s.tab.post <- matrix(0, J, NCOL(Y))
for (i in 1:J) {
  for (j in 1:NCOL(Y)) {
    n.s.tab.post[i,j] <- sum(s$t.c[[i]][,1]==j)
  }
}

# True vs estimated Simpson's diversity index for population 1
getSimpson(s.hpy$tab, s.hpy$n.s.tab, 1, 1, 0.5, 5, 0.3)
getSimpson(s$t.c, n.s.tab.post, 1, gamma.post, alpha.post, theta.post[1], sigma.post[1])
# Sample Simpson's diversity index for population 1
p1 <- Y[1,]/sum(Y[1,])
p1 <- p1[p1!=0]
1-sum(p1^2)

# Credible intervals
p <- c(0.025, 0.975)
gamma.ci <- quantile(s$gamma, probs = p)
alpha.ci <- quantile(s$alpha, probs = p)
theta.ci <- t(apply(s$theta, 2, quantile, probs = p))
sigma.ci <- t(apply(s$sigma, 2, quantile, probs = p))

###########################
# Checking traceplots for parallel chains ####
col <- rainbow(n.chain)

# tables
pdf(paste0(plot.dir, "/table_traceplots_fixed.pdf"), height=8, width=12)
par(mfrow=c(5,2), mar=rep(1,4))
for (j in 1:J) {
  tmp <- unlist(lapply(s.p, function(x) {x$n.tab[,j]}))
  y.rge <- range(tmp, s.hpy$t.tab[j])
  plot(s.p[[1]]$n.tab[,j], type="l", ylim=y.rge, col=col[1])
  for (i in 1:n.chain) {
    lines(s.p[[i]]$n.tab[,j], ylim=y.rge, col=col[[i]])
    abline(h=s.hpy$t.tab[j], col="black", lty=2, lwd=2)
  }
}
dev.off()

# top-level concentration param
pdf(paste0(plot.dir, "/gamma_traceplot_fixed.pdf"), height=8, width=12)
par(mfrow=c(1,1), mar=rep(2,4))
tmp <- unlist(lapply(s.p, function(x) {x$gamma}))
y.rge <- range(tmp, conc.top)
plot(s.p[[1]]$gamma, type="l", ylim=y.rge, col=col[1])
for (i in 1:n.chain) {
  lines(s.p[[i]]$gamma, ylim=y.rge, col=col[[i]])
}
abline(h=conc.top, col="black", lty=2, lwd=2)
dev.off()

# top-level discount param
pdf(paste0(plot.dir, "/alpha_traceplot_fixed.pdf"), height=8, width=12)
par(mfrow=c(1,1), mar=rep(2,4))
tmp <- unlist(lapply(s.p, function(x) {x$alpha}))
y.rge <- range(tmp, dsct.top)
plot(s.p[[1]]$alpha, type="l", ylim=y.rge, col=col[1])
for (i in 1:n.chain) {
  lines(s.p[[i]]$alpha, ylim=y.rge, col=col[[i]])
}
abline(h=dsct.top, col="black", lty=2, lwd=2)
dev.off()

# local-level concentration params
pdf(paste0(plot.dir, "/theta_traceplots_fixed.pdf"), height=8, width=12)
par(mfrow=c(5,2), mar=rep(1,4))
for (j in 1:min(10,J)) {
  tmp <- unlist(lapply(s.p, function(x) {x$theta[,j]}))
  y.rge <- range(tmp, conc.local[j])
  plot(s.p[[1]]$theta[,j], type="l", ylim=y.rge, col=col[1])
  for (i in 1:n.chain) {
    lines(s.p[[i]]$theta[,j], ylim=y.rge, col=col[[i]])
    abline(h=conc.local[j], col="black", lty=2, lwd=2)
  }
}
dev.off()

# local-level discount params
pdf(paste0(plot.dir, "/sigma_traceplots_fixed.pdf"), height=8, width=12)
par(mfrow=c(5,2), mar=rep(1,4))
for (j in 1:min(10,J)) {
  tmp <- unlist(lapply(s.p, function(x) {x$sigma[,j]}))
  y.rge <- range(tmp, dsct.local[j])
  plot(s.p[[1]]$sigma[,j], type="l", ylim=y.rge, col=col[1])
  for (i in 1:n.chain) {
    lines(s.p[[i]]$sigma[,j], ylim=y.rge, col=col[[i]])
    abline(h=dsct.local[j], col="black", lty=2, lwd=2)
  }
}
dev.off()

# Testing individual functions for sampler for concentration parameter ####
# set.seed(24601)
# J <- 10
# conc <- 0.5
# dsct <- 0.5
# p.shape <- 10
# p.scale <- 10
# N <- rep(1000, J)
# Y <- matrix(100, J, 10)
# n.tab <- rep(100, J)
# n.s.tab <- matrix(10, J, 10)

# Checking top level
q <- rbeta(1, conc.top, sum(s.hpy$t.tab))
Q <- 1/p.scale - sum(log(q))

mc <- map_conc(conc.top, p.shape, Q, s.hpy$n.sp, dsct.top)
#s.samp <- slice_conc(mc, p.shape, Q, n.tab, dsct)

# Testing the final sampler for the concentration parameter
n.samp <- 1000
conc.samp <- rep(0, n.samp)
for (i in 1:n.samp) {
  conc.samp[i] <- samp_conc(conc.top, p.shape, 
                            p.scale, s.hpy$n.sp, 1, sum(s.hpy$t.tab), dsct.top, Q=Q)
}

xv <- seq(min(conc.samp), max(conc.samp), by=0.01)
yv <- rep(0, length(xv))
for (i in 1:length(xv)) {
  yv[i] <- prob_conc(xv[i], p.shape, Q, s.hpy$n.sp, dsct.top)
}

# Scaling curve to plot over histogram
h.vals <- hist(conc.samp, breaks=20)
max.hval <- max(h.vals$counts)
lines(xv, max.hval*exp(yv-max(yv)))
abline(v=mc, col="red", lwd=2, lty=2)


# Testing the sampler for the discount parameter
t.tots <- matrix(apply(s.hpy$n.s.tab, 2, sum), ncol=NCOL(Y))

xv <- seq(0.01, 0.99, by=0.001)
yv <- rep(0, length(xv))
for (i in 1:length(xv)) {
  yv[i] <- prob_dsct(xv[i], conc.top, NCOL(Y), t(rep(1, NCOL(Y))), t.tots)
}

n.samp <- 1000
dsct.samp <- rep(0, n.samp)
dsct.samp.buntine <- rep(0, n.samp)
for (i in 1:n.samp) {
  #print(i)
  dsct.samp[i] <- samp_dsct(0.5, conc.top, NCOL(Y), t(rep(1, NCOL(Y))), t.tots)
  dsct.samp.buntine[i] <- gStirling::sampDsct(0.9, NCOL(Y), NCOL(Y), 
                                    t.tots, t(rep(1, NCOL(Y))), conc.top)
}

# Scaling curve to plot over histogram
sc <- 1
h.vals <- hist(dsct.samp, breaks=20, col="blue",
               xlim=range(dsct.samp, dsct.samp.buntine))
max.hval <- max(h.vals$counts)
lines(xv, sc*max.hval*exp(yv-max(yv)))
hist(dsct.samp.buntine, breaks=20, add=TRUE, col="red")


