library(dplyr)

# with respect to Park SW's simulation code

##' @param x a vector of nodes
##' @param size number of nodes to pick at random
sample2 <- function(x, size) {
  if(length(x)==1) {
    rep(x, size)
  } else{
    sample(x, size, replace=TRUE) 
  }
}

makeconcrete <- function(x){
  if(x > (floor(x) + 0.5)){
    x <- floor(x) + 1
  }
  if(x < (floor(x) + 0.5)){
    x <- floor(x)
  }
  return(x)
}

mysir <- function(size,
                  R0=2.5,
                  alphaGT,
                  betaGT,
                  alphaIP,
                  betaIP,
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
    
    
    incubation <- rgamma(1, alphaIP, betaIP)
    t_symptomatic[j] <- t+incubation
    
    c_infected <- c_infected +1
    
    ncontact <- rpois(1, R0)
    # the number of cases cotacted by infector j
    
    n <- V[V != j]
    # V is the total population
    
    if (ncontact > 0) {
      queue_v <- c(queue_v, sample2(n, ncontact)) 
      queue_infector <- c(queue_infector, rep(j, ncontact))
    }
    
    generation <- rgamma(ncontact, alphaGT, betaGT)
    
    
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
      epidata = data.frame(
        time = t_gillespie[(I0):c_infected],
        infected = (I0):c_infected
      ),
      intrinsic_generation = intrinsic_generation,
      pairdata = data.frame(
        Infector.ID = infected_by[!is.na(infected_by)],
        Infectee.ID = which(!is.na(infected_by)),
        infector_Tinfect = t_infected[infected_by[!is.na(infected_by)]],
        infector_Tonset = t_symptomatic[infected_by[!is.na(infected_by)]],
        infectee_Tinfect = t_infected[which(!is.na(infected_by))],
        infectee_Tonset = t_symptomatic[which(!is.na(infected_by))],
        infector_TinfectD = unlist(lapply(t_infected[infected_by[!is.na(infected_by)]], 
                                          makeconcrete)),
        infector_onsetD = unlist(lapply(t_symptomatic[infected_by[!is.na(infected_by)]], 
                                        makeconcrete)),
        infectee_TinfectD = unlist(lapply(t_infected[which(!is.na(infected_by))], 
                                          makeconcrete)),
        infectee_onsetD = unlist(lapply(t_symptomatic[which(!is.na(infected_by))], 
                                        makeconcrete))
      )
    )
  )
}

# test <- mysir(size=1000, R0=2.5, alphaGT=(8/4)^2, betaGT= 8/(4^2), alphaIP=(6.5/3.5)^2, betaIP=6.5/(3.5^2), 
#              I0=10, seed=1, keep.intrinsic = FALSE)
# test$pairdata
# View(test$pairdata)
