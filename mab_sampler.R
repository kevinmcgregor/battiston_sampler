# Functions for "multi-armed bandit" sampler from Battiston et al. (2018)


#' Multi-armed bandit sampler
#'
#' @param Y Matrix of taxa counts (rows are populations, columns are species)
#' @param n.iter Number of MCMC iterations
#' @param quiet If TRUE, console output is suppressed
#'
#' @return
#' @export
#'
#' @examples
mab_sampler <- function(Y, n.iter, quiet=FALSE) {
  if (!is.numeric(Y) | !is.matrix(Y) |
      any(Y<0) | any(Y!=floor(Y))) stop("Y must be a numeric matrix of positive counts")
  
  J <- NROW(Y)
  K <- NCOL(Y) # Number of distinct species in joint sample
  
  # Initializing parameters
  beta0 <- 0
  gamma <- 1 # Top-level concentation
  sigma <- 0.5 <- # Top-level discount
    
  # Containers
  n <- rowSums(Y) # The number of individuals in each population
  
  # Loop over additional samples
  for (i in 1:n.iter) {
    if (!quiet & i%%100==0) cat("Sample: ", i, "\n")
    
    # Loop over populations
    for (j in 1:J) {
      
      # Loop over individuals in a population
      for (p in 1:n[j]){
        # Remove current individual from its table
        
      }
    }
  }
  
  return(1)
}

