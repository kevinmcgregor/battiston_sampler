# Functions for "multi-armed bandit" sampler from Battiston et al. (2018)


#' HPY sampler from Battiston et al. (2018)
#'
#' @param Y Matrix of taxa counts (rows are populations, columns are species)
#' @param n.iter Number of MCMC iterations
#' @param quiet If TRUE, console output is suppressed
#'
#' @return
#' @export
#'
#' @examples
b_sampler <- function(Y, n.iter, quiet=FALSE) {
  if (!is.numeric(Y) | !is.matrix(Y) |
      any(Y<0) | any(Y!=floor(Y))) stop("Y must be a numeric matrix of positive counts")
  
  J <- NROW(Y)
  K <- NCOL(Y) # Number of distinct species in joint sample
  
  # Initializing parameters
  gamma <- 1 # Top-level concentation
  sigma <- 0.5 <- # Top-level discount
    
  # Containers
  n <- rowSums(Y) # The number of individuals in each population
  tab <- vector("list", length=J) # List to hold table indicators for each population
  for (j in 1:J) {
    tab[[j]] <- rep(0, n[j])
  }
  t.c <- vector("list", length=J) # Table counts (and species corresponding to each table)
  n.tab <- rep(0, J) # Number of tables in each population
  
  
  # Loop over MCMC iterations
  for (i in 1:n.iter) {
    if (!quiet & i%%100==0) cat("Sample: ", i, "\n")
    
    # Loop over populations
    for (j in 1:J) {
      
      # Loop over individuals in a population
      for (p in 1:n[j]){
        # Remove current individual from its table
        t.cur <- tab[[j]][p]
        t.c[[j]][t.cur,] <- t.c[[j]][t.cur,] - 1
        if (t.c[[j]][t.cur,]==0) {
          # Drop table and make adjustments
          t.c[[j]] <- t.c[[j]][-t.cur,]
          
        }
        
        p.new <- 0.5
        is.new <- rbinom(1, 1, p.new)
        if (is.new) {
          
        } else {
          
        }
      }
    }
  }
  
  return(1)
}

