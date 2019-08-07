# Checking performance of Battiston sampler in simulated data
# (Ran simulation on Graham cluster)

in.dir <- ""
out.dir <- ""

n.sim <- 100
J <- 10
missing <- c(4,54,84) # Removing missing ones while they're re-run

# Containers for true values
theta <- sigma <- n.tab <- matrix(0, length(n.sim), J)
gamma <- alpha <- rep(0, length(n.sim))

# Containers for estimated values
theta.est <- sigma.est <- n.tab.est <- matrix(0, length(n.sim), J)
gamma.est <- alpha.est <- rep(0, length(n.sim))

# Function to extract post mean over all chains
postMean <- function(l, param) {
  if (!is.null(dim(l[[1]][[param]]))) {
  colMeans(do.call(rbind, lapply(l, function(x){
    apply(x[[param]], 2, mean)
  })))
  } else {
    mean(unlist(lapply(l, function(x){
      mean(x[[param]])
    })))
  }
}

for (s in 1:n.sim) {
  if (s %in% missing) next()
  cat("s =", s, "\n")
  
  # Read in file
  load(paste0(in.dir, "/b_sampler_sim", s, ".RData"))
  
  theta[s,] <- conc.local
  sigma[s,] <- dsct.local
  gamma[s] <- conc.top
  alpha[s] <- dsct.top
  n.tab[s,] <- s.hpy$t.tab
  
  theta.est[s,] <- postMean(s.p, "theta")
  sigma.est[s,] <- postMean(s.p, "sigma")
  gamma.est[s] <- postMean(s.p, "gamma")
  alpha.est[s] <- postMean(s.p, "alpha")
  n.tab.est[s,] <- postMean(s.p, "n.tab")
  
}

save(theta, sigma, gamma, alpha, n.tab, theta.est, 
     sigma.est, gamma.est, alpha.est, n.tab.est, 
     file=paste0(out.dir, "/sim_results.RData"))

