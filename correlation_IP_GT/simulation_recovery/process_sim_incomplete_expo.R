# function to add exposure window to cases
add_expowindow_incomplete <- function(inputpair, windowlength, nseed){
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
  
  # after initial phase,
  # make 1/3 infectees completely missing expo, 1/3 infectees missing early expo
  # make 1/3 infectors completely missing expo, 1/3 infectors missing early expo
  # all miss and half miss can be intersection
  
  Infector.ID.uniq.2 <- unique(inputpair[inputpair$infector_onsetD >= 10,]$Infector.ID)
  Infectee.ID.uniq.2 <- unique(inputpair[inputpair$infector_onsetD >= 10,]$Infectee.ID)
  
  id.infector.allmiss <- sample(Infector.ID.uniq.2, round(ninfector/3), replace = F)
  id.infector.halfmiss <- sample(Infector.ID.uniq.2, round(ninfector/3), replace = F)
  for(ii in 1:round(ninfector/3)){
    inputpair[inputpair$Infector.ID == id.infector.allmiss[ii],]$Texpo.early.I <- NA
    inputpair[inputpair$Infector.ID == id.infector.allmiss[ii],]$Texpo.late.I <- NA
    inputpair[inputpair$Infector.ID == id.infector.halfmiss[ii],]$Texpo.early.I <- NA
  }
  
  id.infectee.allmiss <- sample(Infectee.ID.uniq.2, round(ninfectee/3), replace = F)
  id.infectee.halfmiss <- sample(Infectee.ID.uniq.2, round(ninfectee/3), replace = F)
  
  for(is in 1:round(ninfectee/3)){
    inputpair[inputpair$Infectee.ID == id.infectee.allmiss[is],]$Texpo.early.S <- NA
    inputpair[inputpair$Infectee.ID == id.infectee.allmiss[is],]$Texpo.late.S <- NA
    inputpair[inputpair$Infectee.ID == id.infectee.halfmiss[is],]$Texpo.early.S <- NA
  }
  
  return(inputpair)
}


### estimate back/for IP

est_IP_incomplete <- function(inputpair, startseq, endseq){
  # requires fitdistrplus package
  # for simulation purpose we just fit gamma distribution
  
  nwindows <- length(startseq)
  
  dfIP.back <- dfIP.forward <- data.frame(
    alpha = numeric(nwindows),
    beta = numeric(nwindows)
  )
  
  
  for(i in 1:nwindows){
    
    subpair <- subset(inputpair, Tonset.I >= startseq[i] & Tonset.I <= endseq[i])
    
    subpairinfector <- subset(subpair, !is.na(Texpo.early.I) & !is.na(Texpo.late.I))
    subpairinfectee <- subset(subpair, !is.na(Texpo.early.S) & !is.na(Texpo.early.S))
    
    infector.uniq <- unique(subpairinfector$Infector.ID)
    ind.infector <- numeric(length(infector.uniq))
    for(j in 1:length(infector.uniq)){
      ind.infector[j] <- which(subpairinfector$Infector.ID == infector.uniq[j])[1]
    }
    
    dfback <- data.frame(left = subpairinfector[ind.infector,]$Tonset.I - subpairinfector[ind.infector,]$Texpo.late.I,
                         right = subpairinfector[ind.infector,]$Tonset.I - subpairinfector[ind.infector,]$Texpo.early.I)
    dffor <- data.frame(left = subpairinfectee$Tonset.S - subpairinfectee$Texpo.late.S,
                        right = subpairinfectee$Tonset.S - subpairinfectee$Texpo.early.S)
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
















