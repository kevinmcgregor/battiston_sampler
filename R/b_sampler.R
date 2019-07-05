# Functions for "multi-armed bandit" sampler from Battiston et al. (2018)


#' HPY sampler from Battiston et al. (2018)
#'
#' @param Y Matrix of taxa counts (rows are populations, columns are species)
#' @param n.iter Number of MCMC iterations (after burn-in)
#' @param quiet If TRUE, console output is suppressed
#' @param n.burn Number of burn-in iterations
#'
#' @return
#' @export
#' @importFrom parallel mclapply
#'
#' @examples
b_sampler <- function(Y, n.iter, n.burn, p.shape, p.scale, 
                      mc.cores=1, quiet=FALSE) {
  if (!is.numeric(Y) | !is.matrix(Y) |
      any(Y<0) | any(Y!=floor(Y))) stop("Y must be a numeric matrix of positive counts")
  
  J <- NROW(Y)
  K <- NCOL(Y) # Number of distinct species in joint sample
  
  # Initializing PY parameters
  gamma <- 1 # Top-level concentation
  alpha <- 0.2 # Top-level discount
  theta <- rep(5, J) # Population-level concentration
  sigma <- rep(0.3, J) # Population-level discount
  
  # Initializing table info
  sp.vec <- vector("list", length=J)
  n <- rowSums(Y) # The number of individuals in each population
  tab <- vector("list", length=J) # List to hold table indicators for each population
  t.c <- vector("list", length=J) # Table counts (and species corresponding to each table)
  n.tab <- rep(0, J) # Number of tables in each population
  n.s.tab <- matrix(1, J, K) # Number of tables for a given species in each population
  for (j in 1:J) {
    # Initially one table for each species
    sp.vec[[j]] <- rep(1:K, Y[j,])
    tab[[j]] <- rep(1:sum(Y[j,]>0), Y[j,][Y[j,]>0])
    t.c[[j]] <- cbind((1:K)[Y[j,]>0], Y[j,][Y[j,]>0])
    n.tab[j] <- NROW(t.c[[j]])
  }
  mc.list <- vector("list", length=J) # List for mclapply
  
  gamma.s <- rep(0, n.iter)
  alpha.s <- rep(0, n.iter)
  theta.s <- matrix(0, n.iter, J)
  sigma.s <- matrix(0, n.iter, J)
  n.tab.s <- matrix(0, n.iter, J)
  n.s.tab.s <- array(0, dim=c(n.iter, J, K))
  
  # Loop over MCMC iterations
  for (i in 1:(n.burn+n.iter)) {
    cat("i =", i, "\n")
    idx <- i - n.burn
    
    if (!quiet & i==1) cat("Beginning burn-in:", "\n")
    if (!quiet & i==n.burn+1) cat("Beginning sampling:", "\n")
    if (!quiet & i%%100==0) cat(" ", i, "\n")
    
    # Loop over populations
    for (j in 1:J) {
      #cat("j =", j, "\n")
      # Loop over individuals in a population
      for (p in 1:n[j]){
        sp.cur <- sp.vec[[j]][p]
        # Remove current individual from its table
        t.cur <- tab[[j]][p]
        t.c[[j]][t.cur,2] <- t.c[[j]][t.cur,2] - 1
        if (t.c[[j]][t.cur,2]==0) {
          # Drop table and make adjustments
          t.c[[j]] <- t.c[[j]][-t.cur,]
          tab[[j]][tab[[j]]>t.cur] <- tab[[j]][tab[[j]]>t.cur] - 1
          n.tab[j] <- n.tab[j] - 1
          n.s.tab[j, sp.cur] <- n.s.tab[j, sp.cur] - 1
        }
        
        # Reassign individual to new or existing table
        #new.t <- (theta[j]+n.tab[j]*sigma[j])*(sum(n.s.tab[,sp.cur])-alpha)/
        #          ((theta[j]+n[j]-1)*(gamma+sum(n.tab)))
        #prob.unsc <- c(freq.t, new.t)
        num <- (theta[j]+n.tab[j]*sigma[j])*(sum(n.s.tab[, sp.cur])-alpha)
        den <- (gamma+sum(n.tab))*(Y[j,sp.cur]-1-n.s.tab[j, sp.cur]*sigma[j]) +
                                (theta[j]+n.tab[j]*sigma[j])*(sum(n.s.tab[,sp.cur])-alpha)
        is.new <- rbinom(1, 1, num/den)
        #if (is.na(is.new)) browser()
        if (is.new) {
          # Allocate new table
          t.c[[j]] <- rbind(t.c[[j]], c(sp.cur, 1))
          n.tab[j] <- n.tab[j] + 1
          n.s.tab[j, sp.cur] <- n.s.tab[j, sp.cur] + 1
          tab[[j]][p] <- n.tab[j]
        } else {
          # Sample existing table (from proper species)
          freq.t <- ifelse(t.c[[j]][,1]==sp.cur, (t.c[[j]][,2]-sigma[j])/(theta[j]+n[j]-1), 0)
          wh.t <- sample(1:length(freq.t), 1, prob=freq.t)
          t.c[[j]][wh.t, 2] <- t.c[[j]][wh.t, 2] + 1
          tab[[j]][p] <- wh.t
        }
      }
      
      # Putting relevant quantities into list for use in mclapply
      mc.list[[j]] <- list(conc=theta[j],n.tab=n.tab[j],n=n[j],
                           dsct=sigma[j], n.s.tab=n.s.tab[j,,drop=FALSE],
                           Y=Y[j,,drop=FALSE], J=1, p.shape=p.shape,
                           p.scale=p.scale)
      
      # theta[j] <- samp_conc(theta[j], p.shape, p.scale, n.tab[j], 1, n[j], sigma[j])
      # sigma[j] <- samp_dsct(sigma[j], theta[j], n.tab[j], 
      #                       n.s.tab[j,,drop=FALSE], Y[j,,drop=FALSE]) 
    }
    
    # Sample local-level PY parameters
    theta <- unlist(parallel::mclapply(mc.list, samp_conc_list, mc.cores = mc.cores))
    sigma <- unlist(parallel::mclapply(mc.list, samp_dsct_list, mc.cores = mc.cores))
    
    for (j in 1:J) {
      if (i>n.burn) {
        theta.s[idx,j] <- theta[j]
        sigma.s[idx,j] <- sigma[j]
        n.tab.s[idx,j] <- n.tab[j]
      }
    }
    
    # Sample top-level PY parameters
    gamma <- samp_conc(gamma, p.shape, p.scale, n.tab, J, n, alpha)
    alpha <- samp_dsct(alpha, gamma, n.tab, n.s.tab, Y) 
    # 
    if (i>n.burn) {
      gamma.s[idx] <- gamma
      alpha.s[idx] <- alpha
      n.s.tab.s[idx,,] <- n.s.tab
    }
  }
  
  return(list(gamma=gamma.s, alpha=alpha.s, theta=theta.s, sigma=sigma.s,
              tab=tab, t.c=t.c, n.tab=n.tab.s, n.s.tab=n.s.tab.s,
              p.shape=p.shape, p.scale=p.scale))
}

#' Sample the concentration parameter hierarchical in Pitman-Yor process
#'
#' @param conc Current value of the concentration parameter
#' @param J Number of populations
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
samp_conc <- function(conc, p.shape, p.scale, n.tab, J, n, dsct) {
  max.conc <- 2000
  
  # Sample auxiliary Beta variables
  q <- rbeta(J, conc, n)
  Q <- 1/p.scale - sum(log(q))
  
  conc.map <- map_conc(conc, p.shape, Q, n.tab, dsct)
  conc.ret <- slice(conc.map, prob_conc, max=max.conc, p.shape=p.shape, 
                         Q=Q, n.tab=n.tab, dsct=dsct)

  return(conc.ret)
}

# Wrapper function for samp_conc for use in mclapply
samp_conc_list <- function(l) {
  samp_conc(l$conc, l$p.shape, l$p.scale, l$n.tab, l$J, l$n, l$dsct)
}

# Sample discount parameter
samp_dsct <- function(dsct, conc, n.tab, n.s.tab, Y) {
  dsct.ret <- slice(dsct, prob_dsct, min=0, max=1, 
                    conc=conc, n.tab=n.tab, n.s.tab=n.s.tab, Y=Y)
  return(dsct.ret)
}

# Wrapper function for samp_dsct for use in mclapply
samp_dsct_list <- function(l) {
  samp_dsct(l$dsct, l$conc, l$n.tab, l$n.s.tab, l$Y)
}

#' Inverse digamma function
#'
#' @param x The value at which to evaluate the inverse of the digamma function 
#'
#' @return The inverse of the digamma function evaluated at x
#' @export
#'
#' @examples
inv_digamma <- function(x) {
  if (x<(-2.22)) {
    g <- -1/(x-digamma(1))
  } else {
    g <- exp(x)+0.5
  }
  
  for(i in 1:5) {
    g <- g - (digamma(g)-x)/trigamma(g)
  }
  
  return(g)
}



map_conc <- function(conc, p.shape, Q, n.tab, dsct, err.tol=1e-4, max.iter=10) {
  J <- length(n.tab)
  i <- 1
  conc.old <- conc*1.1
  while (i < max.iter & abs((conc-conc.old)/conc) > err.tol) {
    #print(conc)
    conc.old <- conc
    p <- sum(digamma(n.tab+conc/dsct))
    p <- p + dsct*(p.shape-1)/(conc) - dsct*Q
    conc <- dsct*inv_digamma(p/J)
    i <- i + 1
  }
  return(conc)
}

prob_conc <- function(conc, p.shape, Q, n.tab, dsct) {
  log_prob <- -conc*Q+(p.shape-1)*log(conc)
  log_prob <- log_prob + sum(lgamma(n.tab+conc/dsct) - lgamma(conc/dsct))
  return(log_prob)
}


#' Title
#'
#' @param dsct 
#' @param conc 
#' @param n.tab 
#' @param n.s.tab 
#' @param Y 
#' @param s.table Optional: previously calculated Stirling number table
#' to reuse for faster computation
#'
#' @return
#' @importFrom gStirling gStirling
#' @export
#'
#' @examples
prob_dsct <- function(dsct, conc, n.tab, n.s.tab, Y, s.table=NULL) {
  # Find max Stirling number to calculate
  mt <- max(n.s.tab)
  my <- max(Y)
  if (is.null(s.table)) {
    s.table <- gStirling::gStirling(my, mt, dsct)
  }

  # Adding the log-Stirling numbers
  log_prob <- 0
  for (i in 1:NROW(Y)) {
    log_prob <- log_prob + sum(s.table[Y[i,]+my*(n.s.tab[i,]-1)])
  }
  
  log_prob <- log_prob + sum(n.tab*log(dsct)+
                    lgamma(n.tab+conc/dsct) - lgamma(conc/dsct))
  
  return(log_prob)
}

#' Slice sampler for a single parameter
#'
#' @param map Initial value of parameter, possibly an MAP estimate
#' @param l.prob Log-probability density function for parameter
#' @param iter Number of iterations to run slice sampler (can be small)
#' @param min Minimum bound for parameter
#' @param max Maximum bound for parameter
#' @param ... Arguments passed on to l.prob()
#'
#' @return A sample from the distribution characterized by l.prob()
#' @export
#'
#' @examples
slice <- function(map, l.prob, iter=5, min=0, max, ...) {
  bounds <- c(min, max)
  x <- map
  
  for (i in 1:iter) {
    #print(conc)
    y <- l.prob(x, ...)
    l <- bounds[1]
    r <- bounds[2]
    y <- y + log(runif(1))
    accept <- FALSE
    while (!accept) {
      x.try <- l + runif(1)*(r-l)
      if (is.na(l.prob(x.try, ...)) | is.na(y)) browser()
      if (l.prob(x.try, ...) > y) {
        x <- x.try
        accept <- TRUE
      } else {
        if (x.try < x) {
          l <- x.try
        } else {
          r <- x.try
        }
      }
    }
  }
  
  return(x)
}

# Find the bounds of a slice
# find_bounds <- function(conc.map, y, p.shape, Q, n.tab, dsct, sl.size, max.size) {
#   U <- runif(1)
#   V <- runif(1)
#   L <- conc.map - sl.size*U
#   R <- L + sl.size
#   J <- floor(max.size*V)
#   K <- (max.size - 1) - J
#   
#   # Finding left bound
#   while (J > 0 & y < prob_conc(L, p.shape, Q, n.tab, dsct)) {
#     L <- L - sl.size
#     J <- J - 1
#     if (L<=0) {
#       L <- 0
#       break()
#     }
#   }
#   
#   # Finding right bound
#   while (K > 0 & y < prob_conc(R, p.shape, Q, n.tab, dsct)) {
#     R <- R + sl.size
#     K <- K - 1
#   }
#   
#   return(c(L,R))
# }


# logStirling <- function(N, M, a) {
#   if (N <= 0 | floor(N)!=N | !is.numeric(N)) stop("N must be a positive integer")
#   if (M <= 0 | M > N | floor(M)!=M | !is.numeric(M)) stop("M must be a positive integer with M<=N")
#   if (!is.numeric(a) | a<=0 | a>=1) stop("a must be numeric and between 0 and 1")
#   
#   # At which value of N should we use the asymptotic form?
#   N.asymp <- 64
#   
#   if (N==M) {
#     r <- 0
#   } else if (M==1) {
#     r <- lgamma(N-a)-lgamma(1-a)
#   } else if (N>=N.asymp) {
#     r <- lgamma(N)-lgamma(1-a)-lgamma(M)-(M-1)*log(a)-a*log(N)
#   } else {
#     s1 <- sum(unlist(lapply(1:(N-M+1), logStirling, M=1, a=a)))
#     r <- s1
#   }
#   
#   return(r)
# }
