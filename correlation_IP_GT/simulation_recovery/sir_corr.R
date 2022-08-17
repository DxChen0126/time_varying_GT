### MEDIUM GT SETTING

# IP mean 6.5 sd 3.5; GT mean 7.0 sd 4.0

# IP lnorm setting
# mean <- 6.5
# sigma <- 3.5
# 
# log(mean^2/sqrt(mean^2 + sigma^2)) # 1.74
# sqrt(log(sigma^2/mean^2 + 1)) # 0.50
# 
# # GT lnorm setting
# mean <- 7.0
# sigma <- 4.0
# 
# log(mean^2/sqrt(mean^2 + sigma^2)) # 1.80
# sqrt(log(sigma^2/mean^2 + 1)) # 0.53



##' @param x a vector of nodes
##' @param size number of nodes to pick at random
sample2 <- function(x, size) {
  if(length(x)==1) {
    rep(x, size)
  } else{
    sample(x, size, replace=TRUE)
  }
}

sir.full2 <- function(size,
                     R0=2.5,
                     meanlog=1.74,
                     sdlog=0.50,
                     meanlog2=1.80,
                     sdlog2=0.53,
                     rho=0,
                     I0,
                     seed = NULL,
                     imax,
                     keep.intrinsic=FALSE){
  
  if (!is.null(seed)) set.seed(seed)
  
  V <- 1:size
  
  if (missing(I0)) {
    if (missing(I0)) stop("specify the initial conditions")
    
  }
  initial_infected <- 1:I0
  
  if (missing(imax)) imax <- size
  
  queue_v <- queue_t <- queue_infector <- rep(NA, I0)
  
  queue_v[1:I0] <- initial_infected 
  queue_t[1:I0] <- 0
  
  t_infected <- t_symptomatic <- rep(NA, size)
  t_infected[initial_infected] <- 0
  
  t_gillespie <- NULL
  c_infected <- 0
  
  if (keep.intrinsic) {
    intrinsic_generation <- vector('list', length(V))
  } else {
    intrinsic_generation <- NULL
  }
  
  done <- rep(FALSE, size)
  infected_by <- rep(NA, size)
  
  stop <- FALSE
  
  while (!stop) {
    j.index <- which.min(queue_t)
    j <- queue_v[j.index]
    
    infected_by[j] <- queue_infector[j.index]
    t_infected[j] <- queue_t[j.index]
    
    t <- queue_t[j.index]; t_gillespie <- c(t_gillespie, t)
    
    incubation <- rlnorm(1, meanlog, sdlog)
    t_symptomatic[j] <- t+incubation
    
    c_infected <- c_infected +1
    
    ncontact <- rpois(1, R0)
    
    n <- V[V != j]
    
    if (ncontact > 0) {
      queue_v <- c(queue_v, sample2(n, ncontact))
      queue_infector <- c(queue_infector, rep(j, ncontact))
    }
    
    muloggen <- meanlog2 + sdlog2 * rho * (log(incubation) - meanlog)/sdlog
    sigmaloggen <- sqrt(sdlog2^2 * (1 - rho^2))
    
    generation <- rlnorm(ncontact, muloggen, sigmaloggen)
    
    if (keep.intrinsic) intrinsic_generation[[j]] <- generation
    
    if (ncontact > 0) {
      queue_t <- c(queue_t, t + generation)
    }
    
    done[j] <- TRUE
    
    filter2 <- !done[queue_v]
    queue_v <- queue_v[filter2]
    queue_infector <- queue_infector[filter2]
    queue_t <- queue_t[filter2]
    
    stop <- (c_infected == length(V) || all(done[queue_v]) || c_infected == imax)
  }
  
  return(
    list(
      data=data.frame(
        time=t_gillespie[(I0):c_infected],
        infected=(I0):c_infected
      ),
      intrinsic_generation=intrinsic_generation,
      t_infected=t_infected,
      t_symptomatic=t_symptomatic,
      infected_by=infected_by
    )
  )
}
