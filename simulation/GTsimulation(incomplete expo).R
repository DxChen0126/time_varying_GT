rm(list = ls())

source("SIR simulation code.R")
source("process_sim_incomplete_expo.R")
source("GT_sampling_functions.R")
source("GT_estimation_function.R")

library(foreach)
library(doParallel)
library(fitdistrplus)
library(dplyr)

# test med GT with 7/4/1 day expo

startseq = c(1, seq(5, 40, 7))
endseq = c(10, seq(11, 39, 7), 60)


load("SIR_sims_list_rawdata.rda")

estimatedGT <- vector("list", 50)
ind.pairdata <- seq(3, 3000, 3)

for(n in 1:50){
  set.seed(n)
  ind.sample <- sample(seq(1, 1000, 1), 1, replace = F)
  simpair <- simlist[[2]][ind.pairdata][ind.sample]
  nwindow <- length(startseq)
 # expolength <- 7
 # expolength <- 1
  expolength <- 4
  # add exposure window
  simpair <- add_expowindow_incomplete(simpair$pairdata, windowlength = expolength, nseed = n)
  
  # estimate IP
  IPestlist <- est_IP_incomplete(simpair, startseq, endseq)
  
  # sample GT
  GTsamplelist <- getGTsamples(dfIP.backward = IPestlist$IPinfector,
                               dfIP.forward = IPestlist$IPinfectee,
                               inputpair = simpair,
                               startseq, endseq)
  
  save(GTsamplelist, file = paste0("./RecordGTsampling_expo_incomplete/", "expo_4day","_", n, "th_","GTsample.rda"))
  rm(simpair); rm(IPestlist)
  
  # estimate GT
  
  GT.est <- GTestimated(GTsamplelist, startseq, endseq)
  estimatedGT[[n]] <- GT.est
  
  save(GT.est, file = paste0("./RecordGTsampling_expo_incomplete/", "expo_4day","_", n, "th_","estimatedGT.rda"))
  
  rm(GTsamplelist); rm(GT.est)
  
  print(n)
}

estimatedGT[[2]]

# save(estimatedGT, file = "GT_est_sim_medGT_incomplete_expo7day.rda")
# save(estimatedGT, file = "GT_est_sim_medGT_incomplete_expo1day.rda")
save(estimatedGT, file = "GT_est_sim_medGT_incomplete_expo4day.rda")

estmuLB <- matrix(unlist(lapply(estimatedGT, function(x) x$muLB)), ncol = 7, nrow = 50,
                  byrow = T)
estmuUB <- matrix(unlist(lapply(estimatedGT, function(x) x$muUB)), ncol = 7, nrow = 50,
                  byrow = T)
estmu <- matrix(unlist(lapply(estimatedGT, function(x) x$mu)), ncol = 7, nrow = 50,
                byrow = T)
mu.mean <- colMeans(estmu)  

estsdLB <- matrix(unlist(lapply(estimatedGT, function(x) x$sdLB)), ncol = 7, nrow = 50,
                  byrow = T)
estsdUB <- matrix(unlist(lapply(estimatedGT, function(x) x$sdUB)), ncol = 7, nrow = 50,
                  byrow = T)
estsd <- matrix(unlist(lapply(estimatedGT, function(x) x$sd)), ncol = 7, nrow = 50,
                byrow = T)
sd.mean <- colMeans(estsd)  

load("GTrealized_approx_parameters_0408.rda")

ref.mu <- realizedGT.approx.pars$empmean[[2]]
ref.sd <- realizedGT.approx.pars$empsd[[2]]

ind.mu <- ind.sd <-  Bias.mu <- Bias.sd <- Width.mu <- Width.sd <- numeric(7)

for(i in 1:7){
  ind.mu[i] <- sum(estmuLB[,i] <= ref.mu[i] & ref.mu[i] <= estmuUB[,i])
  ind.sd[i] <- sum(estsdLB[,i] <= ref.sd[i] & ref.sd[i] <= estsdUB[,i])
  Bias.mu[i] <- mu.mean[i]/ref.mu[i] -1
  Bias.sd[i] <- sd.mean[i]/ref.sd[i] -1
  Width.mu[i] <- mean(estmuUB[,i] - estmuLB[,i])
  Width.sd[i] <- mean(estsdUB[,i] - estsdLB[,i])
}

ind.mu
ind.sd
Bias.mu
Bias.sd
Width.mu
Width.sd

mu.mean
sd.mean


df.medGTsim <- data.frame(
  empmean = realizedGT.approx.pars$empmean[[2]],
  empsd = realizedGT.approx.pars$empsd[[2]],
  empsamplesize.ave = realizedGT.approx.pars$samplesize.ave[[2]],
  timewindow = paste0("D",realizedGT.approx.pars$startseq[[2]], " - ", "D", realizedGT.approx.pars$endseq[[2]]),
  muCIrecover = ind.mu,
  estmu = mu.mean,
  biasmu = Bias.mu,
  widthmu = Width.mu,
  sdCIrecover = ind.sd,
  estsd = sd.mean,
  biassd = Bias.sd,
  widthsd = Width.sd
)

df.medGTsim

write.csv(df.medGTsim, "medGT_4day_incomplete_expo_simSummary.csv")




