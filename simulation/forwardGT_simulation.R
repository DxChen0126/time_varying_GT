rm(list = ls())

source("SIR simulation code.R")
source("process_sim.R")
source("GT_sampling_functions.R")
source("GT_estimation_function.R")


library(foreach)
library(doParallel)
library(fitdistrplus)
library(dplyr)

# save 1000 times raw simulation data

nsim = 1000
ncores <- detectCores()
myCluster <- makeCluster(ncores - 2, # number of cores to use
                         type = "PSOCK") # type of cluster

nset <- 3
mGT <- c(10, 7, 4)
sdGT <- c(6, 4, 2)
simlist <- vector("list", nset)

progressbar <- txtProgressBar(min = 0, max = nset, style = 3)

for(j in 1:nset){
simlist[[j]] <- foreach(i = 1:nsim, .combine = 'append') %dopar% {
  mysir(size=1000, R0=2.5, alphaGT=(mGT[j]/sdGT[j])^2, betaGT= mGT[j]/(sdGT[j]^2), alphaIP=(6.5/3.5)^2, betaIP=6.5/(3.5^2), 
        I0=10, seed=i, keep.intrinsic = FALSE)
  }
setTxtProgressBar(progressbar, j)
}

stopCluster(myCluster)

save(simlist, file = "SIR_sims_list_rawdata.rda")

load("SIR_sims_list_rawdata.rda")

# for different parameter setting we need different time windows

 startseqlist <- list(
   longGTstartseq = c(1, seq(10, 60, 10)),
   mediumGTstartseq = c(1, seq(5, 40, 7)),
   shortGTstartseq = c(1, seq(9, 29, 4))
 )
 
 endseqlist <- list(
   longGTstartseq = c(15, seq(19, 59, 10), 100),
   mediumGTstartseq = c(10, seq(11, 39, 7), 60),
   shortGTstartseq = c(10, seq(12, 28, 4), 40)
 )


# startseqlist <- list(
#   longGTstartseq = c(10, 35, 50),
#   mediumGTstartseq = c(15, 30, 40),
#   shortGTstartseq = c(10, 20, 25)
# )
# 
# endseqlist <- list(
#   longGTstartseq = c(25, 45, 65),
#   mediumGTstartseq = c(25, 35, 55),
#   shortGTstartseq = c(15, 25, 40)
# )

empGTlist <- vector("list", nset)
ind.pairdata <- seq(3, 3000, 3)

ncores <- detectCores()
myCluster <- makeCluster(ncores - 2, # number of cores to use
                         type = "PSOCK") # type of cluster

progressbar <- txtProgressBar(min = 0, max = nset, style = 3)

for(j in 1:nset){
  pairdatalist <- simlist[[j]][ind.pairdata]
  empGTlist[[j]] <- foreach(i = 1:nsim, .combine = 'append') %dopar% {
    empGTrealized(pairdatalist[[i]], startseqlist[[j]], endseqlist[[j]])
  }
  setTxtProgressBar(progressbar, j)
  rm(pairdatalist)
}

stopCluster(myCluster)

save(empGTlist, file = "empirical_GTrealized_list_rawdata_0202.rda")

rm(empGTlist)

realizedGTmean <- realizedGTsd <- sampleGT <- vector("list", nset)
ind.empmean <- seq(1, 4000, 4)
ind.empsd <- seq(2, 4000, 4)
ind.samplesize <- seq(3, 4000, 4)
nwindow <- length(startseqlist[[3]])

for(i in 1:nset){
  empmean <- empGTlist[[i]][ind.empmean]
  empmean <- matrix(unlist(empmean), ncol = nwindow, nrow = nsim, byrow = T)
  realizedGTmean[[i]] <- colMeans(empmean)
  empsd <- empGTlist[[i]][ind.empsd]
  empsd <- matrix(unlist(empsd), ncol = nwindow, nrow = nsim, byrow = T)
  realizedGTsd[[i]] <- colMeans(empsd)
  empsamplesize <- empGTlist[[i]][ind.samplesize]
  empsamplesize <- matrix(unlist(empsamplesize), ncol = nwindow, nrow = nsim, byrow = T)
  sampleGT[[i]] <- colMeans(empsamplesize)
}

realizedGT.approx.pars <- list(
  empmean = realizedGTmean,
  empsd = realizedGTsd,
  samplesize.ave = sampleGT,
  startseq = startseqlist,
  endseq = endseqlist
)

save(realizedGT.approx.pars, file = "GTrealized_approx_parameters_0202.rda")

rm(realizedGT.approx.pars)
rm(realizedGTmean)
rm(realizedGTsd)
rm(empmean)
rm(empsamplesize)
rm(empsd)
########################################################################################
load("SIR_sims_list_rawdata.rda")
######### randomly select 50 simulated SIR rawdata and estimated GT with fabricated exposure
# length of 1, 7 and 14 days; for long/medium and short GT SIR simulations

# make this into a function
 startseqlist <- list(
   longGTstartseq = c(1, seq(10, 60, 10)),
   mediumGTstartseq = c(1, seq(5, 40, 7)),
   shortGTstartseq = c(1, seq(9, 29, 4))
 )
 
 endseqlist <- list(
   longGTstartseq = c(15, seq(19, 59, 10), 100),
   mediumGTstartseq = c(10, seq(11, 39, 7), 60),
   shortGTstartseq = c(10, seq(12, 28, 4), 40)
 )

# startseqlist <- list(
#   longGTstartseq = c(10, 35, 50),
#   mediumGTstartseq = c(15, 30, 40),
#   shortGTstartseq = c(10, 20, 25)
# )
# 
# endseqlist <- list(
#   longGTstartseq = c(25, 45, 65),
#   mediumGTstartseq = c(25, 35, 55),
#   shortGTstartseq = c(15, 25, 40)
# )




GTfromsim <- function(nselect, settingchoice, simlist){
  # nselect is the numeber of simulated data to estimate
  # for example we randomly pict 50 from 1000 simulations
  
  GTlength <- settingchoice[1]
  expolength <- settingchoice[2]
  
  
  # GTlength is either 1, 2, 3, representing long/medium/short GT
  # expolength is exposure length, i.e. 1, 7, or 14 days
  
  estimatedGT <- vector("list", nselect)
  ind.pairdata <- seq(3, 3000, 3)
  
  for(n in 1:nselect){
    set.seed(n)
   
    ind.sample <- sample(seq(1, 1000, 1), 1, replace = F)
    simpair <- simlist[[GTlength]][ind.pairdata][ind.sample]
    
    startseq <- startseqlist[[GTlength]]
    endseq <- endseqlist[[GTlength]]
    nwindow <- length(startseq)
    
    # add exposure window
    simpair <- add_expowindow(simpair$pairdata, windowlength = expolength, nseed = n)
    
    # estimate IP
    IPestlist <- est_IP(simpair, startseq, endseq)
    
    # sample GT
    GTsamplelist <- getGTsamples(dfIP.backward = IPestlist$IPinfector,
                                 dfIP.forward = IPestlist$IPinfectee,
                                 inputpair = simpair,
                                 startseq, endseq)
    
    save(GTsamplelist, file = paste0("./RecordGTsamplingData_0202/",  
                                     "GTchoice", GTlength, "_","expochoice", expolength, "_",n, "th_","_GTsample.rda"))
    rm(simpair); rm(IPestlist)
    
    # estimate GT
    
    GT.est <- GTestimated(GTsamplelist, startseq, endseq)
    estimatedGT[[n]] <- GT.est
    
    save(GT.est, file = paste0("./RecordGTsamplingData_0202/",  
                               "GTchoice", GTlength, "_","expochoice", expolength, "_",n, "th_","estimatedGT.rda"))
    
    rm(GTsamplelist); rm(GT.est)
    
    print(n)
  }
  
  return(estimatedGT)
}
  
settings <- matrix(c(rep(c(1, 2, 3), each = 3), rep(c(14, 7, 1), 3)), 
                   ncol = 2, byrow = F)


GT_est_total <- vector("list", 9)

for(isetting in 1:9){
  GT_est_total[[isetting]] <- GTfromsim(nselect = 50, 
                                        settingchoice = settings[isetting,], 
                                        simlist)
  print(isetting)
}


# save(GT_est_total, file = "GT_est_allsims_0202.rda")

GT_est_total_alt <- GT_est_total

save(GT_est_total_alt, file = "GT_est_allsims.rda")


indicator.mu <- indicator.sd <- estmu.mean <- estsd.mean <- bias.mu <- bias.sd <- width.mu <- width.sd <- vector("list", 9)


for(isetting in 1:9){
  settingchoice = settings[isetting,]
  GTlength <- settingchoice[1]
  expolength <- settingchoice[2]
  startseq <- startseqlist[[GTlength]]
  endseq <- endseqlist[[GTlength]]
  
  estmuLB <- matrix(unlist(lapply(GT_est_total_alt[[isetting]], function(x) x$muLB)), ncol = 7, nrow = 50,
                    byrow = T)
  estmuUB <- matrix(unlist(lapply(GT_est_total_alt[[isetting]], function(x) x$muUB)), ncol = 7, nrow = 50,
                    byrow = T)
  estmu <- matrix(unlist(lapply(GT_est_total_alt[[isetting]], function(x) x$mu)), ncol = 7, nrow = 50,
                  byrow = T)
  mu.mean <- colMeans(estmu)  
  
  estsdLB <- matrix(unlist(lapply(GT_est_total_alt[[isetting]], function(x) x$sdLB)), ncol = 7, nrow = 50,
                    byrow = T)
  estsdUB <- matrix(unlist(lapply(GT_est_total_alt[[isetting]], function(x) x$sdUB)), ncol = 7, nrow = 50,
                    byrow = T)
  estsd <- matrix(unlist(lapply(GT_est_total_alt[[isetting]], function(x) x$sd)), ncol = 7, nrow = 50,
                  byrow = T)
  sd.mean <- colMeans(estsd)  
  
  ref.mu <- realizedGT.approx.pars$empmean[[GTlength]]
  ref.sd <- realizedGT.approx.pars$empsd[[GTlength]]
 
  ind.mu <- ind.sd <-  Bias.mu <- Bias.sd <- Width.mu <- Width.sd <- numeric(7)
  
  for(i in 1:7){
    ind.mu[i] <- sum(estmuLB[,i] <= ref.mu[i] & ref.mu[i] <= estmuUB[,i])
    ind.sd[i] <- sum(estsdLB[,i] <= ref.sd[i] & ref.sd[i] <= estsdUB[,i])
    Bias.mu[i] <- mu.mean[i]/ref.mu[i] -1
    Bias.sd[i] <- sd.mean[i]/ref.sd[i] -1
    Width.mu[i] <- mean(estmuUB[,i] - estmuLB[,i])
    Width.sd[i] <- mean(estsdUB[,i] - estsdLB[,i])
  }
  
    indicator.mu[[isetting]] <- ind.mu
    indicator.sd[[isetting]] <- ind.sd
    estmu.mean[[isetting]] <- mu.mean
    estsd.mean[[isetting]] <- sd.mean
    bias.mu[[isetting]] <- Bias.mu
    bias.sd[[isetting]] <- Bias.sd
    width.mu[[isetting]] <- Width.mu
    width.sd[[isetting]] <- Width.sd
  
  
}


df.longGTsim <- data.frame(
  empmean = realizedGT.approx.pars$empmean[[1]],
  empsd = realizedGT.approx.pars$empsd[[1]],
  empsamplesize.ave = realizedGT.approx.pars$samplesize.ave[[1]],
  timewindow = paste0("D",realizedGT.approx.pars$startseq[[1]], " - ", "D", realizedGT.approx.pars$endseq[[1]]),
  expo14.muCIrecover = indicator.mu[[1]],
  expo14.estmu = estmu.mean[[1]],
  expo14.biasmu = bias.mu[[1]],
  expo14.widthmu = width.mu[[1]],
  expo14.sdCIrecover = indicator.sd[[1]],
  expo14.estsd = estsd.mean[[1]],
  expo14.biassd = bias.sd[[1]],
  expo14.widthsd = width.sd[[1]],
  expo7.muCIrecover = indicator.mu[[2]],
  expo7.estmu = estmu.mean[[2]],
  expo7.biasmu = bias.mu[[2]],
  expo7.widthmu = width.mu[[2]],
  expo7.sdCIrecover = indicator.sd[[2]],
  expo7.estsd = estsd.mean[[2]],
  expo7.biassd = bias.sd[[2]],
  expo7.widthsd = width.sd[[2]],
  expo1.muCIrecover = indicator.mu[[3]],
  expo1.estmu = estmu.mean[[3]],
  expo1.biasmu = bias.mu[[3]],
  expo1.widthmu = width.mu[[3]],
  expo1.sdCIrecover = indicator.sd[[3]],
  expo1.estsd = estsd.mean[[3]],
  expo1.biassd = bias.sd[[3]],
  expo1.widthsd = width.sd[[3]]
)

df.longGTsim



df.medGTsim <- data.frame(
  empmean = realizedGT.approx.pars$empmean[[2]],
  empsd = realizedGT.approx.pars$empsd[[2]],
  empsamplesize.ave = realizedGT.approx.pars$samplesize.ave[[2]],
  timewindow = paste0("D",realizedGT.approx.pars$startseq[[2]], " - ", "D", realizedGT.approx.pars$endseq[[2]]),
  expo14.muCIrecover = indicator.mu[[4]],
  expo14.estmu = estmu.mean[[4]],
  expo14.biasmu = bias.mu[[4]],
  expo14.widthmu = width.mu[[4]],
  expo14.sdCIrecover = indicator.sd[[4]],
  expo14.estsd = estsd.mean[[4]],
  expo14.biassd = bias.sd[[4]],
  expo14.widthsd = width.sd[[4]],
  expo7.muCIrecover = indicator.mu[[5]],
  expo7.estmu = estmu.mean[[5]],
  expo7.biasmu = bias.mu[[5]],
  expo7.widthmu = width.mu[[5]],
  expo7.sdCIrecover = indicator.sd[[5]],
  expo7.estsd = estsd.mean[[5]],
  expo7.biassd = bias.sd[[5]],
  expo7.widthsd = width.sd[[5]],
  expo1.muCIrecover = indicator.mu[[6]],
  expo1.estmu = estmu.mean[[6]],
  expo1.biasmu = bias.mu[[6]],
  expo1.widthmu = width.mu[[6]],
  expo1.sdCIrecover = indicator.sd[[6]],
  expo1.estsd = estsd.mean[[6]],
  expo1.biassd = bias.sd[[6]],
  expo1.widthsd = width.sd[[6]]
)

df.medGTsim


df.shortGTsim <- data.frame(
  empmean = realizedGT.approx.pars$empmean[[3]],
  empsd = realizedGT.approx.pars$empsd[[3]],
  empsamplesize.ave = realizedGT.approx.pars$samplesize.ave[[3]],
  timewindow = paste0("D",realizedGT.approx.pars$startseq[[3]], " - ", "D", realizedGT.approx.pars$endseq[[3]]),
  expo14.muCIrecover = indicator.mu[[7]],
  expo14.estmu = estmu.mean[[7]],
  expo14.biasmu = bias.mu[[7]],
  expo14.widthmu = width.mu[[7]],
  expo14.sdCIrecover = indicator.sd[[7]],
  expo14.estsd = estsd.mean[[7]],
  expo14.biassd = bias.sd[[7]],
  expo14.widthsd = width.sd[[7]],
  expo7.muCIrecover = indicator.mu[[8]],
  expo7.estmu = estmu.mean[[8]],
  expo7.biasmu = bias.mu[[8]],
  expo7.widthmu = width.mu[[8]],
  expo7.sdCIrecover = indicator.sd[[8]],
  expo7.estsd = estsd.mean[[8]],
  expo7.biassd = bias.sd[[8]],
  expo7.widthsd = width.sd[[8]],
  expo1.muCIrecover = indicator.mu[[9]],
  expo1.estmu = estmu.mean[[9]],
  expo1.biasmu = bias.mu[[9]],
  expo1.widthmu = width.mu[[9]],
  expo1.sdCIrecover = indicator.sd[[9]],
  expo1.estsd = estsd.mean[[9]],
  expo1.biassd = bias.sd[[9]],
  expo1.widthsd = width.sd[[9]]
)

df.shortGTsim

df.7window.sim.sum <- rbind(df.longGTsim, df.medGTsim, df.shortGTsim)
df.7window.sim.sum$GTsetting <- c(rep("longGT", 7), rep("mediumGT", 7), rep("shortGT", 7))

write.csv(df.7window.sim.sum, "simulation_summary_7timewindows.csv")



