
# function to add exposure window to cases
add_expowindow <- function(inputpair, windowlength, nseed){
  # inputpair should be the same formart as simpair
  
  ninfectee <- nrow(inputpair)
  set.seed(nseed)
  
  # assign expo window to infectees
  #windowlength <- 1
  #inputpair
  
  inputpair$Texpo.early.S <- inputpair$infectee_TinfectD - runif(ninfectee, 0, windowlength)
  inputpair$Texpo.late.S <- inputpair$infectee_TinfectD + runif(ninfectee, 0, windowlength)
  inputpair$Tonset.S <- inputpair$infectee_onsetD
  # make sure latest expo no later than onset
  ind.adj1 <- which(inputpair$Texpo.late.S >= inputpair$Tonset.S)
  inputpair[ind.adj1,]$Texpo.late.S <- inputpair[ind.adj1,]$Tonset.S
  inputpair[ind.adj1,]$Texpo.early.S <- inputpair[ind.adj1,]$Texpo.late.S - windowlength
  
  # assign expo window to infectors
  Infector.ID.uniq <- unique(inputpair$Infector.ID)
  ninfector <- length(Infector.ID.uniq)
  expo.early.I.ref <- expo.late.I.ref <- numeric(ninfector)
  expo.early.I <- expo.late.I <- numeric(ninfectee)
  for(i in 1:ninfector){
    ind.tmp <- which(inputpair$Infector.ID == Infector.ID.uniq[i])
    ind.ref <- which(inputpair$Infector.ID == Infector.ID.uniq[i])[1]
    expo.early.I.ref[i] <- inputpair[ind.ref,]$infector_TinfectD - runif(1, 0, windowlength)
    expo.late.I.ref[i] <- inputpair[ind.ref,]$infector_TinfectD + runif(1, 0, windowlength)
    expo.early.I[ind.tmp] <- rep(expo.early.I.ref[i], length(ind.tmp))
    expo.late.I[ind.tmp] <- rep(expo.late.I.ref[i], length(ind.tmp))
  }
  
  inputpair$Texpo.early.I <- expo.early.I
  inputpair$Texpo.late.I <- expo.late.I
  inputpair$Tonset.I <- inputpair$infector_onsetD
  
  # make sure latest expo no later than onset
  ind.adj2 <- which(inputpair$Texpo.late.I >= inputpair$Tonset.I)
  inputpair[ind.adj2,]$Texpo.late.I <- inputpair[ind.adj2,]$Tonset.I
  inputpair[ind.adj2,]$Texpo.early.I <- inputpair[ind.adj2,]$Texpo.late.I - windowlength

  return(inputpair)
}


### estimate back/for IP

est_IP <- function(inputpair, startseq, endseq){
  # requires fitdistrplus package
  # for simulation purpose we just fit gamma distribution
  
  nwindows <- length(startseq)
  
  dfIP.back <- dfIP.forward <- data.frame(
    alpha = numeric(nwindows),
    beta = numeric(nwindows)
  )
  
  for(i in 1:nwindows){
    subpair <- subset(inputpair, Tonset.I >= startseq[i] & Tonset.I <= endseq[i])
    infector.uniq <- unique(subpair$Infector.ID)
    ind.infector <- numeric(length(infector.uniq))
    for(j in 1:length(infector.uniq)){
      ind.infector[j] <- which(subpair$Infector.ID == infector.uniq[j])[1]
    }
    dfback <- data.frame(left = subpair[ind.infector,]$Tonset.I - subpair[ind.infector,]$Texpo.late.I,
                         right = subpair[ind.infector,]$Tonset.I - subpair[ind.infector,]$Texpo.early.I)
    dffor <- data.frame(left = subpair$Tonset.S - subpair$Texpo.late.S,
                        right = subpair$Tonset.S - subpair$Texpo.early.S)
    est.back <- fitdistcens(dfback, "gamma")
    est.for <- fitdistcens(dffor, "gamma")
    dfIP.back$alpha[i] <- est.back$estimate[1]
    dfIP.back$beta[i] <- est.back$estimate[2]
    dfIP.forward$alpha[i] <- est.for$estimate[1]
    dfIP.forward$beta[i] <- est.for$estimate[2]
  }
  
  return(list(
    IPinfector = dfIP.back,
    IPinfectee = dfIP.forward
  ))
  
}



est_IP_weight <- function(inputpair, startseq, endseq){
  # requires fitdistrplus package
  # for simulation purpose we just fit gamma distribution
  
  nwindows <- length(startseq)
  
  dfIP.back <- dfIP.forward <- data.frame(
    alpha = numeric(nwindows),
    beta = numeric(nwindows)
  )
  
  for(i in 1:nwindows){
    subpair <- subset(inputpair, Tonset.I >= startseq[i] & Tonset.I <= endseq[i])
#    infector.uniq <- unique(subpair$Infector.ID)
#    ind.infector <- numeric(length(infector.uniq))
#    for(j in 1:length(infector.uniq)){
#      ind.infector[j] <- which(subpair$Infector.ID == infector.uniq[j])[1]
#    }
    dfback <- data.frame(left = subpair$Tonset.I - subpair$Texpo.late.I,
                         right = subpair$Tonset.I - subpair$Texpo.early.I)
    dffor <- data.frame(left = subpair$Tonset.S - subpair$Texpo.late.S,
                        right = subpair$Tonset.S - subpair$Texpo.early.S)
    est.back <- fitdistcens(dfback, "gamma")
    est.for <- fitdistcens(dffor, "gamma")
    dfIP.back$alpha[i] <- est.back$estimate[1]
    dfIP.back$beta[i] <- est.back$estimate[2]
    dfIP.forward$alpha[i] <- est.for$estimate[1]
    dfIP.forward$beta[i] <- est.for$estimate[2]
  }
  
  return(list(
    IPinfector = dfIP.back,
    IPinfectee = dfIP.forward
  ))
  
}


# function to create GT samples

getGTsamples <- function(dfIP.backward, dfIP.forward, inputpair, startseq, endseq){
  # requires previous GT sampling function
  # requires foreach doparallel and fitdistrplus packages
  
  ncores <- detectCores()
  myCluster <- makeCluster(ncores - 2, # number of cores to use
                           type = "PSOCK") # type of cluster
  nwindow <- length(startseq)
  samplesize <- numeric(nwindow)
  reslist <- vector("list", nwindow)
  
 progressbar <- txtProgressBar(min = 0, max = nwindow, style = 3)
  
  for(j in 1:nwindow){
    nsim = 1000 
    IPdistInfector <- c(dfIP.backward$alpha[j], dfIP.backward$beta[j])
    IPdistInfectee <- c(dfIP.forward$alpha[j], dfIP.forward$beta[j])
    tmppair <- subset(inputpair, Tonset.I >= startseq[j] & Tonset.I <= endseq[j])
    samplesize[j] <- length(which(inputpair$Tonset.I >= startseq[j] & inputpair$Tonset.I <= endseq[j]))
    reslist[[j]] <- foreach(i = 1:nsim, .combine = 'append') %dopar% {
      simsamples <- sampleGT(tmppair, inputpair, IPdistInfector = IPdistInfector, InfectorIPdistfrom = "gamma", 
                             IPdistInfectee = IPdistInfectee, InfecteeIPdistfrom = "gamma", iseed = i)
      simsamples[[1]]
    }
    # rm(simsamples)
    setTxtProgressBar(progressbar, j)
  }
  stopCluster(myCluster)
  
  return(
    list(
      GTsamples = reslist,
      samplesize = samplesize
    )
  )
}

# function to estimate GT
estGT <- function(GTsamples, samplesize, startdate, enddate){

  ind.inf <- which(is.infinite(GTsamples))
  if(length(ind.inf) == 0){
    GTsamples <- GTsamples
    samplesize <- samplesize
  }
  if(length(ind.inf) != 0){
    GTsamples <- GTsamples[-ind.inf]
    samplesize <- samplesize - length(ind.inf)/1000
  }
  
  GTmat <- matrix(GTsamples, ncol = samplesize, byrow = T)
  GTmatlist <- split(GTmat, rep(1:ncol(GTmat), each = nrow(GTmat)))
  rm(GTmat)
  
  test <- optim(
    par = c(0,0),
    fn  = llh.fx,
    method = "Nelder-Mead",
    hessian = T,
    y.mc = GTmatlist
  )
  mu.est <- exp(test$par)
  UB.est <- qlnorm(0.975, test$par, sqrt(diag(solve(test$hessian))))
  LB.est <- qlnorm(0.025, test$par, sqrt(diag(solve(test$hessian))))
  
  data.frame(
    mu = unlist(mu.est)[1],
    muLB = unlist(LB.est)[1],
    muUB = unlist(UB.est)[1],
    sd = unlist(mu.est)[2],
    sdLB = unlist(LB.est)[2],
    sdUB = unlist(UB.est)[2],
    samplesize = samplesize,
    timewindow = paste0("D", startdate, " - ", "D", enddate)
  )
}

GTestimated <- function(GTsamplelist, startseq, endseq){
  nwindow <- length(startseq)
  ncores <- detectCores()
  myCluster <- makeCluster(ncores - 2, # number of cores to use
                           type = "PSOCK") # type of cluster
  reslist <- foreach(i = 1:nwindow, .combine = 'rbind') %dopar% {
#    rawGTsamples <- GTsamplelist$GTsamples[[i]]
#    rawGTsamples[is.infinite(rawGTsamples)] <- runif(1, 0, 1)
    estGT(GTsamplelist$GTsamples[[i]], GTsamplelist$samplesize[i],
          startseq[i], endseq[i])
  }
  return(reslist)
} 

#  empirical realized GT


empGTrealized <- function(inputpair, startseq, endseq){
  nwindow <- length(startseq)
  mu <- sigma <- samplesize <- numeric(nwindow)
  for(i in 1:nwindow){
    tmppair <- subset(inputpair, infector_onsetD >= startseq[i] & infector_onsetD <= endseq[i])
    tmpdata <- with(tmppair, infectee_Tinfect - infector_Tinfect)
    mu[i] <- mean(tmpdata)
    sigma[i] <- sd(tmpdata)
    samplesize[i] <- nrow(tmppair)
  }
  return(
    list(
      empmean = mu,
      empsd = sigma,
      samplesize = samplesize,
      timewindow = paste0("D", startseq, " - ", "D", endseq)
    )
  )
}
  

######################## backward est
#  estimate realized GT
empGTback.real <- function(inputpair, startdate, enddate){
  tmp <- subset(inputpair, infectee_onsetD >= startdate & infectee_onsetD <= enddate)
  tmpdata <- with(tmp, infectee_Tinfect - infector_Tinfect)
  samplesize <- nrow(tmp)
  data.frame(
    mu = mean(tmpdata),
    sd = sd(tmpdata),
    samplesize = samplesize,
    timewindow = paste0("D", startdate, " - ", "D", enddate)
  )
}


GTrealized.back <- function(inputpair, startseq, endseq){
  # View(simpair)
  nwindow <- length(startseq)
  ncores <- detectCores()
  myCluster <- makeCluster(ncores - 2, # number of cores to use
                           type = "PSOCK") # type of cluster
  reslist <- foreach(i = 1:nwindow, .combine = 'rbind') %dopar% {
    empGTback.real(inputpair, startseq[i], endseq[i])
  }
  return(reslist)
}

est_IP.back <- function(inputpair, startseq, endseq){
  # requires fitdistrplus package
  # for simulation purpose we just fit gamma distribution
  
  nwindows <- length(startseq)
  
  dfIP.back <- dfIP.forward <- data.frame(
    alpha = numeric(nwindows),
    beta = numeric(nwindows)
  )
  
  for(i in 1:nwindows){
    subpair <- subset(inputpair, Tonset.S >= startseq[i] & Tonset.S <= endseq[i])
    infector.uniq <- unique(subpair$Infector.ID)
    ind.infector <- numeric(length(infector.uniq))
    for(j in 1:length(infector.uniq)){
      ind.infector[j] <- which(subpair$Infector.ID == infector.uniq[j])[1]
    }
    dffor <- data.frame(left = subpair[ind.infector,]$Tonset.I - subpair[ind.infector,]$Texpo.late.I,
                         right = subpair[ind.infector,]$Tonset.I - subpair[ind.infector,]$Texpo.early.I)
    dfback <- data.frame(left = subpair$Tonset.S - subpair$Texpo.late.S,
                        right = subpair$Tonset.S - subpair$Texpo.early.S)
    est.back <- fitdistcens(dfback, "gamma")
    est.for <- fitdistcens(dffor, "gamma")
    dfIP.back$alpha[i] <- est.back$estimate[1]
    dfIP.back$beta[i] <- est.back$estimate[2]
    dfIP.forward$alpha[i] <- est.for$estimate[1]
    dfIP.forward$beta[i] <- est.for$estimate[2]
  }
  
  return(list(
    IPinfector = dfIP.forward,
    IPinfectee = dfIP.back
  ))
  
}


getGTsamples.back <- function(dfIP.backward, dfIP.forward, inputpair, startseq, endseq){
  # requires previous GT sampling function
  # requires foreach doparallel and fitdistrplus packages
  
  ncores <- detectCores()
  myCluster <- makeCluster(ncores - 2, # number of cores to use
                           type = "PSOCK") # type of cluster
  nwindow <- length(startseq)
  samplesize <- numeric(nwindow)
  reslist <- vector("list", nwindow)
  
  progressbar <- txtProgressBar(min = 0, max = nwindow, style = 3)
  
  for(j in 1:nwindow){
    nsim = 1000 
    IPdistInfector <- c(dfIP.forward$alpha[j], dfIP.forward$beta[j])
    IPdistInfectee <- c(dfIP.backward$alpha[j], dfIP.backward$beta[j])
    tmppair <- subset(inputpair, Tonset.S >= startseq[j] & Tonset.S <= endseq[j])
    samplesize[j] <- length(which(inputpair$Tonset.S >= startseq[j] & inputpair$Tonset.S <= endseq[j]))
    reslist[[j]] <- foreach(i = 1:nsim, .combine = 'append') %dopar% {
      simsamples <- sampleGT(tmppair, inputpair, IPdistInfector = IPdistInfector, InfectorIPdistfrom = "gamma", 
                             IPdistInfectee = IPdistInfectee, InfecteeIPdistfrom = "gamma", iseed = i)
      simsamples[[1]]
    }
    # rm(simsamples)
    setTxtProgressBar(progressbar, j)
  }
  stopCluster(myCluster)
  
  return(
    list(
      GTsamples = reslist,
      samplesize = samplesize
    )
  )
}




