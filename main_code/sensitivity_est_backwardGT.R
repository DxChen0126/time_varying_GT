rm(list = ls())

# read data
allpairs <- read.csv("allpairs0908.csv")

# format dates
colexpo <- grep(c("expo"), colnames(allpairs)) 
for(i in colexpo){
  allpairs[, i] <- as.Date(allpairs[, i], "%m/%d/%Y")
}

colonset <- grep(c("Onset"), colnames(allpairs)) 
for(i in colonset){
  allpairs[, i] <- as.Date(allpairs[, i], "%m/%d/%Y")
}

# transfer date into number
appendcol <- c()
for(i in colexpo){
  appendcol <- c(appendcol, as.numeric(allpairs[,i] - as.Date("2020-01-01")))
}
for(i in colonset){
  appendcol <- c(appendcol, as.numeric(allpairs[,i] - as.Date("2020-01-01")))
}

appendcol <- matrix(appendcol, ncol = 6, byrow = F)

allpairs[, 22:27] <- appendcol
colnames(allpairs)[22:27] <- c("Texpo.early.I", "Texpo.late.I", "Texpo.early.S", "Texpo.late.S",
                               "Tonset.I", "Tonset.S")


# best fit IP for backward infectee is weibull, for forward infector is gamma

library(mixdist)


dfIP.backward <- read.csv("dfIPinfectee.backward_fitdistWB.csv")

dfIP.forward <- read.csv("dfIPinfector.forward_fitdistG.csv")


####### estimate GT

source("GT_sampling_functions.R")
source("GT_estimation_function.R")


startseq <- c(as.Date("2020-01-01"), seq.Date(from = as.Date("2020-01-21"), to = as.Date("2020-02-05"), by = "1 day"))
endseq <- c(as.Date("2020-01-26"), seq.Date(from = as.Date("2020-01-27"), to = as.Date("2020-02-10"), by = "1 day"),
            as.Date("2020-02-29"))


library(foreach)
library(doParallel)
library(fitdistrplus)

ncores <- detectCores()

myCluster <- makeCluster(ncores - 1, # number of cores to use
                         type = "PSOCK") # type of cluster

registerDoParallel(myCluster)
reslist <- vector("list", 17)
nwindow <- 17

progressbar <- txtProgressBar(min = 0, max = nwindow, style = 3)

for(j in 1:nwindow){
  nsim = 1000 
  IPdistInfector <- c(dfIP.forward$alpha[j], dfIP.forward$beta[j])
  IPdistInfectee <- c(dfIP.backward$alpha[j], dfIP.backward$beta[j])
  tmppair <- subset(allpairs, Onset_Infectee >= startseq[j] & Onset_Infectee <= endseq[j])  
  reslist[[j]] <- foreach(i = 1:nsim, .combine = 'append') %dopar% {
    simsamples <- sampleGT(tmppair, allpairs, IPdistInfector = IPdistInfector, InfectorIPdistfrom = "gamma", 
                           IPdistInfectee = IPdistInfectee, InfecteeIPdistfrom = "weibull", iseed = i)
    list(
      GTsamples = simsamples[[1]],
      IP.S.lb = simsamples[[2]],
      IP.S.ub = simsamples[[3]],
      Tinfector.infect = simsamples[[4]]
    )
  }
  setTxtProgressBar(progressbar, j)
}
stopCluster(myCluster)

save(reslist, file = "backwardGT1000sim_samples_MAY06.rda")

load("backwardGT1000sim_samples_MAY06.rda")

samplesize <- numeric(17)
for(i in 1:17){
  samplesize[i] <- length(which(allpairs$Onset_Infectee >= startseq[i] & allpairs$Onset_Infectee <= endseq[i]))
}

refseq <- seq(1, 4000)
indGT <- refseq[refseq %% 4 == 1]
GTlist <- vector("list", 17)
for(i in 1:17){
  GTsamples <- reslist[[i]][indGT]
  GTmat <- matrix(unlist(GTsamples), ncol = samplesize[i], byrow = T)
  GTmatlist <- split(GTmat, rep(1:ncol(GTmat), each = nrow(GTmat)))
  GTlist[[i]] <- GTmatlist
  rm(GTsamples)
  rm(GTmat)
  rm(GTmatlist)
}
rm(reslist)

source("GT_estimation_function.R")

pars.est <- vector("list", 17)
negsumllk <- vector("list", 17)

# lapply(GTlist, function(x) sum(is.infinite(unlist(x))))

progressbar <- txtProgressBar(min = 0, max = nwindow, style = 3)

for(i in 1:17){
  
  test <- optim(
    par = c(0, 0),
    # fn  = llh.fx,
     fn = llh.fx.wb,
    # fn = llh.fx.ln,
    method = "Nelder-Mead",
    hessian = T,
    y.mc = GTlist[[i]]
  )
  
  # for wb only
  tmppar <- exp(test$par)
  pars.est[[i]] <- as.numeric(weibullparinv(shape = tmppar[1], scale = tmppar[2])[1:2])
  pars.est[[i]] <- exp(test$par)
  negsumllk[[i]] <- test$value
  
  setTxtProgressBar(progressbar, i)
}

indseq <- seq(1, 34)
indpar1 <- indseq[indseq %% 2 == 1]
indpar2 <- indseq[indseq %% 2 == 0]


df.estGT <- data.frame(
  mu = unlist(pars.est)[indpar1],
  sd = unlist(pars.est)[indpar2],
  samplesize = samplesize,
  negsumllk = unlist(negsumllk),
  timewindow = paste0(format(startseq, format =  "%b %d"), " - ",
                      format(endseq, format =  "%b %d"))
)

df.estGT.g <- df.estGT
df.estGT.wb <- df.estGT
df.estGT.ln <- df.estGT

write.csv(df.estGT.ln, "backwardGT_est_lnorm.csv")
write.csv(df.estGT.g, "backwardGT_est_gamma.csv")
write.csv(df.estGT.wb, "backwardGT_est_weibull.csv")

df.estGT.g <- read.csv("backwardGT_est_gamma.csv")

# gamma best fit

# bootCI

par1UB.boot <- par2UB.boot <- numeric(17)
par1LB.boot <- par2LB.boot <- numeric(17)
parsboot <- vector("list", 17)

ncores <- detectCores()
myCluster <- makeCluster(ncores - 2, # number of cores to use
                         type = "PSOCK") # type of cluster

registerDoParallel(myCluster)

start.t <- Sys.time()

progressbar <- txtProgressBar(min = 0, max = 17, style = 3)

for(i in 1:17){
  
  y.mc <- GTlist[[i]]
  
  TMP <- foreach(
    b = 1:1000
  ) %dopar% {
    
    set.seed(b)
    
    ind.b = sample(1:length(y.mc), length(y.mc), replace = T)
    y.mc.b = y.mc[ind.b]
    
    test = optim(
      par = c(0,0),
      # fn  = llh.fx.ln,
      # fn = llh.fx.wb,
      fn = llh.fx,
      method = "Nelder-Mead",
      hessian = T,
      y.mc = y.mc.b
    )
    return(test)
  }
  
   tmpboot <- sapply(TMP, function(b) { exp(b$par) })
  
  # for wb only
  # tmpboot <- sapply(TMP, function(b) {as.numeric(weibullparinv(shape = exp(b$par)[1], scale = exp(b$par)[2])[1:2])})
  
  parsboot[[i]] <- tmpboot
  
  par1UB.boot[i] <- quantile(tmpboot[1,], 0.975)
  par2UB.boot[i] <- quantile(tmpboot[2,], 0.975)
  
  par1LB.boot[i] <- quantile(tmpboot[1,], 0.025)
  par2LB.boot[i] <- quantile(tmpboot[2,], 0.025)
  
  rm(TMP)
  rm(tmpboot)
  rm(y.mc)
  
  setTxtProgressBar(progressbar, i)
}


stopCluster(myCluster)


end.t = Sys.time()
end.t - start.t 

save(parsboot, file = "backwardGT1000_bootpars_gamma_MAY06.rda")

rm(parsboot)

df.estGT.g$muLB.boot <- par1LB.boot
df.estGT.g$muUB.boot <- par1UB.boot
df.estGT.g$sdLB.boot <- par2LB.boot
df.estGT.g$sdUB.boot <- par2UB.boot

write.csv(df.estGT.g, "backwardGT_est_gamma.csv")



