rm(list = ls())

source("SIR simulation code.R")
source("process_sim.R")
source("GT_sampling_functions.R")
source("GT_estimation_function.R")


library(foreach)
library(doParallel)
library(fitdistrplus)
library(dplyr)

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

nset <- 3
mGT <- c(10, 7, 4)
sdGT <- c(6, 4, 2)
nsim <- 1000

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

save(empGTlist, file = "empirical_GTrealized_list_rawdata_0408.rda")

# rm(empGTlist)

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

save(realizedGT.approx.pars, file = "GTrealized_approx_parameters_0408.rda")
load("GTrealized_approx_parameters_0408.rda")

rm(realizedGT.approx.pars)
rm(realizedGTmean)
rm(realizedGTsd)
rm(empmean)
rm(empsamplesize)
rm(empsd)

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

settings <- matrix(c(c(1, 2, 3), rep(4, 3)), 
                   ncol = 2, byrow = F)


GT_est_total <- vector("list", 3)

for(isetting in 1:3){
  GT_est_total[[isetting]] <- GTfromsim(nselect = 50, 
                                        settingchoice = settings[isetting,], 
                                        simlist)
  print(isetting)
}


save(GT_est_total, file = "GT_est_addsims_0408.rda")


estmuLB <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$muLB)), ncol = 7, nrow = 50,
                  byrow = T)
estmuUB <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$muUB)), ncol = 7, nrow = 50,
                  byrow = T)
estmu <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$mu)), ncol = 7, nrow = 50,
                byrow = T)
mu.mean <- colMeans(estmu)  

ref <- realizedGT.approx.pars$empmean[[3]]
realizedGT.approx.pars$samplesize.ave[[3]]

indicator <- numeric(7)
bias <- numeric(7)
width <- numeric(7)
for(i in 1:7){
 indicator[i] <- sum(estmuLB[,i] <= ref[i] & ref[i] <= estmuUB[,i])
  bias[i] <- mu.mean[i]/ref[i] -1
  width[i] <- mean(estmuUB[,i] - estmuLB[,i])
}

mu.mean
bias
width
indicator



indicator.mu <- indicator.sd <- estmu.mean <- estsd.mean <- bias.mu <- bias.sd <- width.mu <- width.sd <- vector("list", 3)


for(isetting in 1:3){
  settingchoice = settings[isetting,]
  GTlength <- settingchoice[1]
  expolength <- settingchoice[2]
  startseq <- startseqlist[[GTlength]]
  endseq <- endseqlist[[GTlength]]
  
  estmuLB <- matrix(unlist(lapply(GT_est_total[[isetting]], function(x) x$muLB)), ncol = 7, nrow = 50,
                    byrow = T)
  estmuUB <- matrix(unlist(lapply(GT_est_total[[isetting]], function(x) x$muUB)), ncol = 7, nrow = 50,
                    byrow = T)
  estmu <- matrix(unlist(lapply(GT_est_total[[isetting]], function(x) x$mu)), ncol = 7, nrow = 50,
                  byrow = T)
  mu.mean <- colMeans(estmu)  
  
  estsdLB <- matrix(unlist(lapply(GT_est_total[[isetting]], function(x) x$sdLB)), ncol = 7, nrow = 50,
                    byrow = T)
  estsdUB <- matrix(unlist(lapply(GT_est_total[[isetting]], function(x) x$sdUB)), ncol = 7, nrow = 50,
                    byrow = T)
  estsd <- matrix(unlist(lapply(GT_est_total[[isetting]], function(x) x$sd)), ncol = 7, nrow = 50,
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
  expo4.muCIrecover = indicator.mu[[1]],
  expo4.estmu = estmu.mean[[1]],
  expo4.biasmu = bias.mu[[1]],
  expo4.widthmu = width.mu[[1]],
  expo4.sdCIrecover = indicator.sd[[1]],
  expo4.estsd = estsd.mean[[1]],
  expo4.biassd = bias.sd[[1]],
  expo4.widthsd = width.sd[[1]]
)

df.longGTsim



df.medGTsim <- data.frame(
  empmean = realizedGT.approx.pars$empmean[[2]],
  empsd = realizedGT.approx.pars$empsd[[2]],
  empsamplesize.ave = realizedGT.approx.pars$samplesize.ave[[2]],
  timewindow = paste0("D",realizedGT.approx.pars$startseq[[2]], " - ", "D", realizedGT.approx.pars$endseq[[2]]),
  expo4.muCIrecover = indicator.mu[[2]],
  expo4.estmu = estmu.mean[[2]],
  expo4.biasmu = bias.mu[[2]],
  expo4.widthmu = width.mu[[2]],
  expo4.sdCIrecover = indicator.sd[[2]],
  expo4.estsd = estsd.mean[[2]],
  expo4.biassd = bias.sd[[2]],
  expo4.widthsd = width.sd[[2]]
)

df.medGTsim


df.shortGTsim <- data.frame(
  empmean = realizedGT.approx.pars$empmean[[3]],
  empsd = realizedGT.approx.pars$empsd[[3]],
  empsamplesize.ave = realizedGT.approx.pars$samplesize.ave[[3]],
  timewindow = paste0("D",realizedGT.approx.pars$startseq[[3]], " - ", "D", realizedGT.approx.pars$endseq[[3]]),
  expo4.muCIrecover = indicator.mu[[3]],
  expo4.estmu = estmu.mean[[3]],
  expo4.biasmu = bias.mu[[3]],
  expo4.widthmu = width.mu[[3]],
  expo4.sdCIrecover = indicator.sd[[3]],
  expo4.estsd = estsd.mean[[3]],
  expo4.biassd = bias.sd[[3]],
  expo4.widthsd = width.sd[[3]]
)

df.shortGTsim

df.expo4.sim.sum <- rbind(df.longGTsim, df.medGTsim, df.shortGTsim)
df.expo4.sim.sum$GTsetting <- c(rep("longGT", 7), rep("mediumGT", 7), rep("shortGT", 7))

write.csv(df.expo4.sim.sum, "simulation_summary_expo4d_0409.csv")







