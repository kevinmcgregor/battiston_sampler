# Testing Buntine's sampler on basic PY process (not hierarchical)
# Testing samplers for concentration and discount

source("~/research/pitman_yor/basic_py/R/py_functions.R")

# Small sampler for concentration/discount assuming table counts are known
cd_samp <- function(s, n.samp, n.burn) {
  n <- length(s$samp)
  p.shape <- 1
  p.scale <- 1
  
  gamma <- 1
  alpha <- 0.5
  
  gamma.s <- alpha.s <- rep(0, n.samp)
  
  for (i in 1:(n.samp+n.burn)) {
    idx <- i - n.burn
    if (i==1) cat("Beginning burn-in:", "\n")
    if (i==n.burn+1) cat("Beginning sampling:", "\n")
    if (i%%100==0) cat(" ", i, "\n")
    
    gamma <- gStirling::sampConc(gamma, n, s$t.s, p.shape, p.scale, alpha)
    alpha <- gStirling::sampDsct(alpha, s$t.s, s$t.s, s$t.c[,2], rep(1, s$t.s), gamma)

    if (i>n.burn) {
      gamma.s[idx] <- gamma
      alpha.s[idx] <- alpha
    }
  }
  
  return(list(gamma=gamma.s, alpha=alpha.s))
}


N <- 1000
conc.vals <- c(0.1, 1, 5, 10)
dsct.vals <- c(0.2, 0.4, 0.6, 0.8)
n.rep <- 4

err.conc <- err.dsct <- array(0, c(n.rep, length(conc.vals), length(dsct.vals)))
for (r in 1:n.rep) {
  for (i in 1:length(conc.vals)) {
    cat ("i =", i, ", conc =", conc.vals[i], "\n")
    for (j in 1:length(dsct.vals)) {
      cat (" j =", j, ", dsct =", dsct.vals[i], "\n")
      s.py <- sample_py(N, conc.vals[i], dsct.vals[j])
      s <- cd_samp(s.py, 1000, 200)
      err.conc[r,i,j] <- (mean(s$gamma) - conc.vals[i])/conc.vals[i]
      err.dsct[r,i,j] <- (mean(s$alpha) - dsct.vals[j])/dsct.vals[i]
    }
  }
}

# Averaging over replications
err.conc.avg <- apply(err.conc, 2:3, mean)
err.dsct.avg <- apply(err.dsct, 2:3, mean)

# Plotting results
# Conc
#pdf(paste0(plot.dir, "/conc_prior3.pdf"), height=10, width=12)
par(mfrow=c(2,2), mai=rep(1,4), mar=rep(4,4))
for (i in 1:length(conc.vals)) {
  plot(conc.vals, err.conc.avg[i,], type="o", col="darkgreen", lwd=2, 
       ylim=c(min(0, min(err.conc.avg[i,])), max(err.conc.avg[i,])),
       main=paste0("Discount = ", dsct.vals[i]),
       xlab="True concentration", ylab="Error of estimate",
       cex.lab=1.4, cex.axis=1.4)
  abline(h=0, col="red", lty=2, lwd=2)
}
#dev.off()



