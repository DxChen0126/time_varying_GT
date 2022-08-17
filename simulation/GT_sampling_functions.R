# we have to do truncated gamma given exposure windows
# with respect to R Programs for Truncated Distributions by Saralees Nadarajah and Samuel Kotz.

qtrunc <- function(p, spec, lb = -Inf, ub = Inf, ...)
{
  tt <- p
  G <- get(paste("p", spec, sep = ""), mode = "function")
  Gin <- get(paste("q", spec, sep = ""), mode = "function")
  tt <- Gin(G(lb, ...) + p*(G(ub, ...) - G(lb, ...)), ...)
  return(tt)
}

rtrunc <- function(n, spec, lb = -Inf, ub = Inf, ...)
{
  x <- u <- runif(n, min = 0, max = 1)
  x <- qtrunc(u, spec, lb = lb, ub = ub,...)
  return(x)
}


# for a given infector in the pairs, sample its possible IP values based on distribution

sampleIP.infector <- function(inputdata, allpair, IPdist, distfrom){
  # distfrom can be "gamma", "weibull", "lnorm"
  # for gamma input IPdist is shape and rate
  # for weibull input IPdist is shape and scale
  # for lognormal input IPdist is meanlog and sdlog
  ID.I <- inputdata$Infector.ID
  Tonset.I <- inputdata$Tonset.I
  Texpo.early.I <- inputdata$Texpo.early.I
  Texpo.late.I <- inputdata$Texpo.late.I
  
  refpair <- subset(allpair, Infector.ID == ID.I)
  
  # bound at infectee's earliest late exposure
  Texpo.S.early <- min(refpair$Texpo.late.S, na.rm = T)
  
  # in case all infectees without exposure
  if(is.infinite(Texpo.S.early)){
    Texpo.S.early <- NA
  }
  # the earliest onset time of all infectees
  Tonset.S.early <- min(refpair$Tonset.S)
  
  #sampling
  # 1) If only infector onset available
  if(is.na(Texpo.early.I) & is.na(Texpo.late.I)){
    tmp <- min(Tonset.S.early, Texpo.S.early, Tonset.I, na.rm = T)
    # there should be a time gap between infector's infection to infector onset 
    # or the earliest point that infectee could be infected 
    # (time for infector to become infectious)
    # that's why + runif(1, 0, 1)
    lowerb <- Tonset.I - tmp + runif(1, 0, 1)
    upperb <- ifelse(lowerb >= 21, lowerb + 7, 21)
    # some extremely long IP cases
  }
  
  # 2) If infector's expo.late and onset available but expo.early NA
  if(is.na(Texpo.early.I) & !is.na(Texpo.late.I)){
    tmp <- min(Tonset.S.early, Texpo.S.early, Texpo.late.I, na.rm = T)
    # if the lower bound is at infectee's earliest onset or earliest expo point
    # or lower bound equals to infector onset, the again runif(1, 0, 1) to
    # ensure from infection to being infectious/onset
    tmp <- ifelse(tmp == Texpo.late.I & tmp < Tonset.I, tmp, tmp - runif(1, 0, 1))
    lowerb <- Tonset.I - tmp 
    upperb <- ifelse(lowerb >= 21, lowerb + 7, 21)
    # some extremely long IP cases
  }
  
  # 3) If infector's expo.early and onset available but expo.late NA
  if(!is.na(Texpo.early.I) & is.na(Texpo.late.I)){
    # in case expo and onset same day
    upperb <- max(Tonset.I - Texpo.early.I, 1)
    tmp <- min(Texpo.S.early, Tonset.S.early, Tonset.I, na.rm = T) - runif(1, 0, 1)
    # if expo.S.early prior than expo.early.I, the expo.early.S should be limited 
    # to later than infector's earliest expo
    lowerb <- ifelse(tmp <= Texpo.early.I, upperb - runif(1, 0, 1), Tonset.I - tmp)
  }
  
  # 4) if all info available
  if(!is.na(Texpo.early.I) & !is.na(Texpo.late.I)){
    # in case expo and onset same day
    upperb <- max(Tonset.I - Texpo.early.I, 1)
    tmp <- min(Texpo.S.early, Tonset.S.early, Texpo.late.I, na.rm = T)
    tmp <- ifelse(tmp == Texpo.late.I & tmp < Tonset.I, tmp, tmp - runif(1, 0, 1))
    # if expo.s.early prior than expo.early.I, the expo.early.S should be limited 
    # to later than infector's earliest expo
    lowerb <- ifelse(tmp <= Texpo.early.I, upperb - runif(1, 0, 1), Tonset.I - tmp)
  }
  
  if(distfrom == "gamma"){
  IP.infector <- rtrunc(1, "gamma", lb = lowerb, ub = upperb,
                        shape = IPdist[1], rate = IPdist[2])
  }

  if(distfrom == "weibull"){
    IP.infector <- rtrunc(1, "weibull", lb = lowerb, ub = upperb,
                          shape = IPdist[1], scale = IPdist[2])
  }

  if(distfrom == "lnorm"){
    IP.infector <- rtrunc(1, "lnorm", lb = lowerb, ub = upperb,
                          meanlog = IPdist[1], sdlog = IPdist[2])
  }
  
  return(IP.infector)
}


# give each row of entire dataframe its corresponding infector's IP

matchIP.infector <- function(inputpair, allpair, IPdist, distfrom){
  allinfectorIP <- unlist(
    lapply(1:nrow(inputpair), 
           function(i){sampleIP.infector(inputpair[i,], allpair, IPdist, distfrom)}))
  refID <- unique(inputpair$Infector.ID)
  n.uniq <- length(refID)
  # make sure same infector has same sampled values
  for(i in 1:n.uniq){
    inds <- which(inputpair$Infector.ID == refID[i])
    allinfectorIP[inds] <- allinfectorIP[inds[1]]
  }
  return(allinfectorIP)
}

# sample infectee's IP based on previous infector's IP
sampleIP.infectee <- function(inputdata, IPdist, distfrom){
  
  Tonset.S <- inputdata$Tonset.S
  Texpo.early.S <- inputdata$Texpo.early.S
  Texpo.late.S <- inputdata$Texpo.late.S
  
  Tonset.I <- inputdata$Tonset.I
  # we cannot simultaneously sample IP infector and infectee
  # because one infector may have multiple infectee
  # before sample IP infectee, we have to fix infector's IP
  IP.I <- inputdata$IP.infector
  Tinfect.I <- Tonset.I - IP.I
  
  # 1) if only infectee onset available
  
  if(is.na(Texpo.early.S) & is.na(Texpo.late.S)){
    tmp <- Tonset.S - Tinfect.I 
    upperb <- ifelse(tmp <= 1, tmp, tmp - runif(1, 0, 1))
    lowerb <- ifelse(upperb <= 1, 0, runif(1, 0, 1))
  }
  
  # 2) if infectee latest exposure and onset available, but no earliest expo
  if(is.na(Texpo.early.S) & !is.na(Texpo.late.S)){
    tmp <- Tonset.S - Tinfect.I 
    upperb <- ifelse(tmp <= 1, tmp, tmp - runif(1, 0, 1))
    # infectee latest expo should be later than infector's infection as well
    # try to make ub > lb
    tmp2 <- ifelse(Tonset.S - Texpo.late.S >= upperb, upperb - runif(1, 0, 1),
                   Tonset.S - Texpo.late.S)
    lowerb <- ifelse(tmp2 <= 0, 0, tmp2)
  }
  
  # 3) if infectee earliest expo and onset available, but no latest expo
  if(!is.na(Texpo.early.S) & is.na(Texpo.late.S)){
    # earliest expo should be later than infector's infection
    tmp <- Tonset.S - Tinfect.I 
    tmp <- ifelse(tmp <= 1, tmp, tmp - runif(1, 0, 1))
    upperb <- ifelse(Texpo.early.S <= Tinfect.I, 
                     tmp, Tonset.S - Texpo.early.S)
    # in case onset and expo at same day
    upperb <- ifelse(upperb == 0, min(1, tmp), upperb)
    lowerb <- ifelse(upperb <= 1, 0, runif(1, 0, 1))
  }
  
  # 4) if all info available
  if(!is.na(Texpo.early.S) & !is.na(Texpo.late.S)){
    # earliest expo should be later than infector's infection
    tmp <- Tonset.S - Tinfect.I
    tmp <- ifelse(tmp <= 1, tmp, tmp - runif(1, 0, 1))
    upperb <- ifelse(Texpo.early.S <= Tinfect.I, 
                     tmp, Tonset.S - Texpo.early.S)
    # in case onset and expo at same day
    upperb <- ifelse(upperb == 0, min(1, tmp), upperb)
    # infectee latest expo should be later than infector's infection as well
    # infectee latest expo should be later than infector's infection as well
    # try to make ub > lb
    tmp2 <- ifelse(Tonset.S - Texpo.late.S >= upperb, upperb - runif(1, 0, 1),
                   Tonset.S - Texpo.late.S)
    lowerb <- ifelse(tmp2 <= 0, 0, tmp2)
  }
  
  if(distfrom == "gamma"){
    IP.infectee <- rtrunc(1, "gamma", lb = lowerb, ub = upperb,
                          shape = IPdist[1], rate = IPdist[2])
  }
  
  if(distfrom == "weibull"){
    IP.infectee <- rtrunc(1, "weibull", lb = lowerb, ub = upperb,
                          shape = IPdist[1], scale = IPdist[2])
  }
  
  if(distfrom == "lnorm"){
    IP.infectee <- rtrunc(1, "lnorm", lb = lowerb, ub = upperb,
                          meanlog = IPdist[1], sdlog = IPdist[2])
  }
  
  return(c(IP.infectee, lowerb, upperb, Tinfect.I))
}

# sample GT
sampleGT <- function(inputpair, allpair, IPdistInfector, InfectorIPdistfrom, IPdistInfectee, InfecteeIPdistfrom, iseed){
  #  n.total <- nrow(inputpair)
  #  GTsim <- numeric(n.total)
  set.seed(iseed)
  inputpair$IP.infector <- matchIP.infector(inputpair, allpair, IPdistInfector, InfectorIPdistfrom)
  extendinfo <- unlist(lapply(1:nrow(inputpair),
                                         function(i){sampleIP.infectee(inputpair[i,], IPdistInfectee, InfecteeIPdistfrom)}))
  refseq <- seq(1, length(extendinfo))
  indIP.infectee <- refseq[refseq %% 4 == 1]
  indlowerb <- refseq[refseq %% 4 == 2]
  indupperb <- refseq[refseq %% 4 == 3]
  indtinfect.I <- refseq[refseq %% 4 == 0]
  inputpair$IP.infectee <- extendinfo[indIP.infectee]
  infecteeIPlb <- extendinfo[indlowerb]
  infecteeIPub <- extendinfo[indupperb]
  Tinfect.I <- extendinfo[indtinfect.I]
  
  GTsim <- with(inputpair, Tonset.S - IP.infectee - Tonset.I + IP.infector)
  
  return(list(
    GTvalues = GTsim,
    IP.S.lb = infecteeIPlb,
    IP.S.ub = infecteeIPub,
    Tinfector.infect = Tinfect.I)
  )
}


### simulate infectee onset (check if GT matches data)
#first use the GTsim values to extend the pairs info

# simSonset <- function(inputpair, IPdistInfectee, InfecteeIPdistfrom, GTdist, GTdistfrom){
#   
#   inputpair$GTlb <- with(inputpair, Tonset.S - infecteeIPub - Tinfector.infect)
#   inputpair$GTub <- with(inputpair, Tonset.S - infecteeIPlb - Tinfector.infect)
#   
#   if(GTdistfrom == "gamma"){
#     Tinfectee.infect <- unlist(lapply(1:nrow(inputpair),
#     function(i){inputpair[i,]$Tinfector.infect + rtrunc(1, "gamma", lb = inputpair[i,]$GTlb, ub = inputpair[i,]$GTub,
#                 shape = GTdist[1], rate = GTdist[2])}))
#   }
#   
#   if(GTdistfrom == "weibull"){
#     Tinfectee.infect <- unlist(lapply(1:nrow(inputpair),
#                 function(i){inputpair[i,]$Tinfector.infect + rtrunc(1, "weibull", lb = inputpair[i,]$GTlb, ub = inputpair[i,]$GTub,
#                 shape = GTdist[1], scale = GTdist[2])}))
#   }
#   
#   if(GTdistfrom == "weibull"){
#     Tinfectee.infect <- unlist(lapply(1:nrow(inputpair),
#                 function(i){inputpair[i,]$Tinfector.infect + rtrunc(1, "lnorm", lb = inputpair[i,]$GTlb, ub = inputpair[i,]$GTub,
#                 meanlog = GTdist[1], sdlog = GTdist[2])}))
#   }
#   
#   if(InfecteeIPdistfrom == "gamma"){
#     IPinfectee.sim <- unlist(lapply(1:nrow(inputpair),
#                                     function(i){ rtrunc(1, "gamma", lb = inputpair[i,]$GTlb, 
#                                                ub = inputpair[i,]$GTub, shape = GTdist[1], rate = GTdist[2])}))
#       
#       
#       
#       
#       rgamma(nrow(inputpair), shape = IPdistInfectee[1], rate = IPdistInfectee[2])
#   }
#   
#   if(InfecteeIPdistfrom == "weibull"){
#     IPinfectee.sim <- rweibull(nrow(inputpair), shape = IPdistInfectee[1], scale = IPdistInfectee[2])
#   }
#   
#   if(InfecteeIPdistfrom == "lnorm"){
#     IPinfectee.sim <- rlnorm(nrow(inputpair), meanlog = IPdistInfectee[1], sdlog = IPdistInfectee[2])
#   }
#   
#   Tonset.S.sim <- Tinfectee.infect + IPinfectee.sim
#   
#   return(Tonset.S.sim)
# }











