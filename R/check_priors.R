# Checking parameter estimates as a function of the prior hyperparams
# Make sure local-level params and table counts are fixed in sampler before running this
source("~/research/pitman_yor/battiston_sampler/R/b_sampler.R")
source("~/research/pitman_yor/basic_py/R/hpy_functions.R")

plot.dir <- "/mnt/GREENWOOD_BACKUP/home/kevin.mcgregor/research/pitman_yor/battiston_sampler/plots/check_priors/"

# Sampling from HPY
N <- 1000
J <- 25
conc.top <- 5 #15
dsct.top <- 0.9 #0.5
conc.local <- rep(5, J)
dsct.local <- rep(0.3, J)
s.hpy <- sample_hpy(J, rep(N, J), conc.top, dsct.top, conc.local, dsct.local)

Y <- s.hpy$abund

# Shape/scale values in gamma prior to loop over
shape.vals <- c(0.01, 0.1, 1, 5, 10, 25)
scale.vals <- c(0.01, 0.1, 1, 5, 10, 25)

# Containers for error of top-level concentration and discount params
conc.err <- dsct.err <- matrix(0, length(shape.vals), length(scale.vals))

for (i in 1:length(shape.vals)) {
  cat("i =", i, "\n")
  for (j in 1:length(scale.vals)) {
    cat(" j =", j, "\n")
    p.sh <- shape.vals[i]
    p.sc <- scale.vals[j]
    
    s <- b_sampler(Y, 2000, 500, p.shape = p.sh, p.scale = p.sc)
    conc.err[i,j] <- mean(s$gamma) - conc.top
    dsct.err[i,j] <- mean(s$alpha) - dsct.top
  }
}

# Plotting results
# Conc
pdf(paste0(plot.dir, "/conc_prior3.pdf"), height=10, width=12)
par(mfrow=c(3,2), mai=rep(1,4), mar=rep(4,4))
for (i in 1:length(shape.vals)) {
  plot(scale.vals, conc.err[i,], type="o", col="darkgreen", lwd=2, 
       ylim=c(min(0, min(conc.err[i,])), max(conc.err[i,])),
       main=paste0("Prior shape = ", shape.vals[i]),
       xlab="Prior scale", ylab="Error of estimate",
       cex.lab=1.4, cex.axis=1.4)
  abline(h=0, col="red", lty=2, lwd=2)
}
dev.off()

# Dsct
pdf(paste0(plot.dir, "/dsct_prior3.pdf"), height=10, width=12)
par(mfrow=c(3,2), mai=rep(1,4), mar=rep(4,4))
for (i in 1:length(shape.vals)) {
  plot(scale.vals, dsct.err[i,], type="o", col="darkgreen", lwd=2, 
       ylim=c(min(dsct.err[i,]), max(0, max(dsct.err[i,]))),
       main=paste0("Prior shape = ", shape.vals[i]),
       xlab="Prior scale", ylab="Error of estimate",
       cex.lab=1.4, cex.axis=1.4)
  abline(h=0, col="red", lty=2, lwd=2)
}
dev.off()

