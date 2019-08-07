# Checking top-level concentration param
# Simulating data under different top-level discount params and checking bias
# of top-level concentration param
library(parallel)

source("~/research/pitman_yor/battiston_sampler/R/b_sampler.R")
source("~/research/pitman_yor/basic_py/R/hpy_functions.R")

# Params
N <- 1000
J <- 25
conc.top <- 0.5
dsct.vals <- seq(0.1,0.9, by=0.05)
conc.local <- rep(5, J)
dsct.local <- rep(0.3, J)
p.shape = sqrt(10)
p.scale = sqrt(10)

# When running sampler, make sure to fix top-level discount in code
n.chain <- 4
Y.list <- vector("list", n.chain)
bias <- rep(0, length(dsct.vals))
for (i in 1:length(dsct.vals)) {
  cat("dsct = ", dsct.vals[i], "\n")
  
  dsct.top <- dsct.vals[i]
  
  # Simulate data
  s.hpy <- sample_hpy(J, rep(N, J), conc.top, dsct.top, conc.local, dsct.local)
  Y <- s.hpy$abund

  for (j in 1:n.chain) Y.list[[j]] <- Y
  
  s.p <- mclapply(Y.list, b_sampler, p.shape = p.shape, p.scale = p.scale,
                  n.iter=1000, n.burn=200, mc.cores=n.chain)
  
  # Bias of top-level concentration
  gamma.samp <- unlist(lapply(s.p, function(x){x$gamma}))
  bias[i] <- mean(gamma.samp) - conc.top
}

plot(dsct.vals, bias, type="o", lwd=2)
abline(h=0, lty=2, col="red", lwd=2)
