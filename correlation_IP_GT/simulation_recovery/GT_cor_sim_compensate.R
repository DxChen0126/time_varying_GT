source("process_sim.R")
source("GT_sampling_corr_functions.R")
source("GT_estimation_corr_function.R")


library(foreach)
library(doParallel)
library(fitdistrplus)
library(dplyr)

load("sir_sim_corr_0628.rda")

startseq <- seq(5, 30, 5)
endseq <- seq(10, 35, 5)

######### medGT expo 7 # compensate 20 times sim

nselect = 20

tstart <- Sys.time()

GT_est_total <- vector("list", 4)

for(m in 1:4){
  
  estimatedGT <- vector("list", nselect)
  
  for(n in 51: (50 + nselect)){
    
    set.seed(n)
    ind.sample <- sample(seq(1, 1000, 1), 1, replace = F)
    
    simdata <- simlist_corr[[m]][ind.sample][[1]]
    pairdata <- get_pair_from_sim(simdata)
    
    
    # add exposure window
    simpair <- add_expowindow(pairdata, windowlength = 7, nseed = n)
    
    # estimate IP
    IPestlist <- tryCatch(expr = {est_IP_corr(simpair, startseq, endseq)},
                          error = function(e){
                            return("IP estimation failed")
                          })
    if(IPestlist != c("IP estimation failed")){
      
      # sample GT
      GTsamplelist <- getGTsamples(dfIP.backward = IPestlist$IPinfector,
                                   dfIP.forward = IPestlist$IPinfectee,
                                   inputpair = simpair,
                                   startseq, endseq)
      
      save(GTsamplelist, file = paste0("./record_sim/",  
                                       "settingchoice_", m, "_","meanexpo7","_",n, "th_","_GTsample.rda"))
      
      # estimate GT
      
      GT.est <- tryCatch(expr = {GTestimated(GTsamplelist, startseq, endseq, IPestlist$IPinfector)},
                         error = function(e){
                           return("estimation failed")
                         })
      
      estimatedGT[[n]] <- GT.est
      
      save(GT.est, file = paste0("./record_sim/",  
                                 "settingchoice_", m, "_","meanexpo7", "_", n, "th_","estimatedGT.rda"))
      
      rm(GTsamplelist); rm(GT.est)
      
    }
    print(n)
  }
  
  GT_est_total[[m]] <- estimatedGT
  
  rm(estimatedGT)
}

tend <- Sys.time()

tend-tstart

View(GT_est_total[[1]])

ind.GT1 <- c()

for(i in 51:70){
  if(is.null(GT_est_total[[1]][[i]]) | GT_est_total[[1]][i] == "estimation failed"){
    ind.GT1 <- c(ind.GT1, i)
  }
}

ind.GT1 # 2 times failed

ind.GT2 <- c()

for(i in 51:70){
  if(is.null(GT_est_total[[2]][[i]]) | GT_est_total[[2]][i] == "estimation failed"){
    ind.GT2 <- c(ind.GT2, i)
  }
}

ind.GT2 # 4 times failed

ind.GT3 <- c()

for(i in 51:70){
  if(is.null(GT_est_total[[3]][[i]]) | GT_est_total[[3]][i] == "estimation failed"){
    ind.GT3 <- c(ind.GT3, i)
  }
}

ind.GT3 # 4 times failed

ind.GT4 <- c()

for(i in 51:70){
  if(is.null(GT_est_total[[4]][[i]]) | GT_est_total[[4]][i] == "estimation failed"){
    ind.GT4 <- c(ind.GT4, i)
  }
}

ind.GT4 # 3 times failed

GT_est_total[[1]] <- GT_est_total[[1]][-ind.GT1]

GT_est_total[[2]] <- GT_est_total[[2]][-ind.GT2]

GT_est_total[[3]] <- GT_est_total[[3]][-ind.GT3]

GT_est_total[[4]] <- GT_est_total[[4]][-ind.GT4]


save(GT_est_total, file = "GT_est_corrsims_medGT_expo7_compensate_0801.rda")


load("GT_est_corrsims_medGT_expo7_compensate_0801.rda")


estsd.rho0 <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$sd)), ncol = 6, 
                     byrow = T)

estsdLB.rho0 <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$sdLB)), ncol = 6, 
                       byrow = T)

estsdUB.rho0 <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$sdUB)), ncol = 6, 
                       byrow = T)


# - 22, - 24


estmu.rho0 <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$mu)), ncol = 6,
                     byrow = T)


estmuLB.rho0 <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$muLB)), ncol = 6, 
                       byrow = T)
estmuUB.rho0 <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$muUB)), ncol = 6,
                       byrow = T)


estsd.rho25 <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$sd)), ncol = 6, 
                      byrow = T)


estmu.rho25 <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$mu)), ncol = 6,
                      byrow = T)


estmuLB.rho25 <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$muLB)), ncol = 6, 
                        byrow = T)
estmuUB.rho25 <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$muUB)), ncol = 6,
                        byrow = T)



estsd.rho50 <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$sd)), ncol = 6, 
                      byrow = T)

estsdLB.rho50 <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$sdLB)), ncol = 6, 
                        byrow = T)
estsdUB.rho50 <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$sdUB)), ncol = 6, 
                        byrow = T)




estmu.rho5 <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$mu)), ncol = 6,
                     byrow = T)


estmuLB.rho5 <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$muLB)), ncol = 6, 
                       byrow = T)
estmuUB.rho5 <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$muUB)), ncol = 6,
                       byrow = T)

# - 45

estsd.rho75 <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$sd)), ncol = 6, 
                      byrow = T)

estsdLB.rho75 <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$sdLB)), ncol = 6, 
                        byrow = T)
estsdUB.rho75 <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$sdUB)), ncol = 6, 
                        byrow = T)


estmu.rho75 <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$mu)), ncol = 6,
                      byrow = T)


estmuLB.rho75 <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$muLB)), ncol = 6, 
                        byrow = T)
estmuUB.rho75 <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$muUB)), ncol = 6,
                        byrow = T)

# - 15


load("GT_est_corrsims_medGT_expo7.rda")


estsd.rho0 <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$sd)), ncol = 6, 
                     byrow = T)

estsdLB.rho0 <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$sdLB)), ncol = 6, 
                       byrow = T)

estsdUB.rho0 <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$sdUB)), ncol = 6, 
                       byrow = T)


# - 22, - 24


estmu.rho0 <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$mu)), ncol = 6,
                     byrow = T)


estmuLB.rho0 <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$muLB)), ncol = 6, 
                       byrow = T)
estmuUB.rho0 <- matrix(unlist(lapply(GT_est_total[[1]], function(x) x$muUB)), ncol = 6,
                       byrow = T)


estsd.rho25 <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$sd)), ncol = 6, 
                      byrow = T)


estmu.rho25 <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$mu)), ncol = 6,
                      byrow = T)


estmuLB.rho25 <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$muLB)), ncol = 6, 
                        byrow = T)
estmuUB.rho25 <- matrix(unlist(lapply(GT_est_total[[2]], function(x) x$muUB)), ncol = 6,
                        byrow = T)



estsd.rho50 <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$sd)), ncol = 6, 
                      byrow = T)

estsdLB.rho50 <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$sdLB)), ncol = 6, 
                        byrow = T)
estsdUB.rho50 <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$sdUB)), ncol = 6, 
                        byrow = T)




estmu.rho5 <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$mu)), ncol = 6,
                     byrow = T)


estmuLB.rho5 <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$muLB)), ncol = 6, 
                       byrow = T)
estmuUB.rho5 <- matrix(unlist(lapply(GT_est_total[[3]], function(x) x$muUB)), ncol = 6,
                       byrow = T)

# - 45

estsd.rho75 <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$sd)), ncol = 6, 
                      byrow = T)

estsdLB.rho75 <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$sdLB)), ncol = 6, 
                        byrow = T)
estsdUB.rho75 <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$sdUB)), ncol = 6, 
                        byrow = T)


estmu.rho75 <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$mu)), ncol = 6,
                      byrow = T)


estmuLB.rho75 <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$muLB)), ncol = 6, 
                        byrow = T)
estmuUB.rho75 <- matrix(unlist(lapply(GT_est_total[[4]], function(x) x$muUB)), ncol = 6,
                        byrow = T)

load("empirical_corr_ref.rda")


# chechk the performance on mean expo1 medium setting of rho 0, 0.25, 0.5, 0.75

indicator.mu <- indicator.sd <- indicator.rho <- vector("list", 4)

estmu.mean <- estsd.mean <- bias.mu <- bias.sd <- width.mu <- width.sd <- vector("list", 4)
estrho.mean <- bias.rho <- width.rho <- vector("list", 4)

for(i in 1:4){
  estmuLB <- matrix(unlist(lapply(GT_est_total[[i]], function(x) x$muLB)), ncol = 6, 
                    byrow = T)
  estmuUB <- matrix(unlist(lapply(GT_est_total[[i]], function(x) x$muUB)), ncol = 6,
                    byrow = T)
  estmu <- matrix(unlist(lapply(GT_est_total[[i]], function(x) x$mu)), ncol = 6,
                  byrow = T)
  mu.mean <- colMeans(estmu)  
  
  estsdLB <- matrix(unlist(lapply(GT_est_total[[i]], function(x) x$sdLB)), ncol = 6,
                    byrow = T)
  estsdUB <- matrix(unlist(lapply(GT_est_total[[i]], function(x) x$sdUB)), ncol = 6, 
                    byrow = T)
  estsd <- matrix(unlist(lapply(GT_est_total[[i]], function(x) x$sd)), ncol = 6, 
                  byrow = T)
  sd.mean <- colMeans(estsd)  
  
  estrhoLB <- matrix(unlist(lapply(GT_est_total[[i]], function(x) x$rhoLB)), ncol = 6, 
                     byrow = T)
  estrhoUB <- matrix(unlist(lapply(GT_est_total[[i]], function(x) x$rhoUB)), ncol = 6, 
                     byrow = T)
  estrho <- matrix(unlist(lapply(GT_est_total[[i]], function(x) x$rho)), ncol = 6,
                   byrow = T)
  rho.mean <- colMeans(estrho)  
  
  ref.mu <- df.sum[[i]]$meanGTforonset
  ref.sd <- df.sum[[i]]$sdGTforonset
  ref.rho <- df.sum[[i]]$corres.GTIPback
  
  ind.mu <- ind.sd <-  Bias.mu <- Bias.sd <- Width.mu <- Width.sd <- numeric(6)
  ind.rho <- Bias.rho <- Width.rho <- numeric(6)
  
  for(j in 1:6){
    ind.mu[j] <- sum(estmuLB[, j] <= ref.mu[j] & ref.mu[j] <= estmuUB[, j])
    ind.sd[j] <- sum(estsdLB[, j] <= ref.sd[j] & ref.sd[j] <= estsdUB[, j])
    ind.rho[j] <- sum(estrhoLB[, j] <= ref.rho[j] & ref.rho[j] <= estrhoUB[, j])
    
    Bias.mu[j] <- mu.mean[j]/ref.mu[j] -1
    Bias.sd[j] <- sd.mean[j]/ref.sd[j] -1
    Bias.rho[j] <- rho.mean[j]/ref.rho[j] -1
    
    Width.mu[j] <- mean(estmuUB[, j] - estmuLB[, j])
    Width.sd[j] <- mean(estsdUB[, j] - estsdLB[, j])
    Width.rho[j] <- mean(estrhoUB[, j] - estrhoLB[, j])
  }
  
  indicator.mu[[i]] <- ind.mu
  indicator.sd[[i]] <- ind.sd
  indicator.rho[[i]] <- ind.rho
  
  estmu.mean[[i]] <- mu.mean
  estsd.mean[[i]] <- sd.mean
  estrho.mean[[i]] <- rho.mean
  
  bias.mu[[i]] <- Bias.mu
  bias.sd[[i]] <- Bias.sd
  bias.rho[[i]] <- Bias.rho
  
  width.mu[[i]] <- Width.mu
  width.sd[[i]] <- Width.sd
  width.rho[[i]] <- Width.rho
  
}

df.expo7GTsim_rho0 <- data.frame(
  empmean = df.sum[[1]]$meanGTforonset,
  empsd = df.sum[[1]]$sdGTforonset,
  emprho = df.sum[[1]]$corres.GTIPback,
  empsamplesize.ave = df.sum[[1]]$samplesize.GTIPback,
  timewindow = paste0("D", startseq, " - ", "D", endseq),
  muCIrecover = indicator.mu[[1]],
  estmu = estmu.mean[[1]],
  biasmu = bias.mu[[1]],
  widthmu = width.mu[[1]],
  sdCIrecover = indicator.sd[[1]],
  estsd = estsd.mean[[1]],
  biassd = bias.sd[[1]],
  widthsd = width.sd[[1]],
  rhoCIrecover = indicator.rho[[1]],
  estrho = estrho.mean[[1]],
  biasrho = bias.rho[[1]],
  widthrho = width.rho[[1]]
)


df.expo7GTsim_rho_0.25 <- data.frame(
  empmean = df.sum[[2]]$meanGTforonset,
  empsd = df.sum[[2]]$sdGTforonset,
  emprho = df.sum[[2]]$corres.GTIPback,
  empsamplesize.ave = df.sum[[2]]$samplesize.GTIPback,
  timewindow = paste0("D", startseq, " - ", "D", endseq),
  muCIrecover = indicator.mu[[2]],
  estmu = estmu.mean[[2]],
  biasmu = bias.mu[[2]],
  widthmu = width.mu[[2]],
  sdCIrecover = indicator.sd[[2]],
  estsd = estsd.mean[[2]],
  biassd = bias.sd[[2]],
  widthsd = width.sd[[2]],
  rhoCIrecover = indicator.rho[[2]],
  estrho = estrho.mean[[2]],
  biasrho = bias.rho[[2]],
  widthrho = width.rho[[2]]
) 


df.expo7GTsim_rho_0.5 <- data.frame(
  empmean = df.sum[[3]]$meanGTforonset,
  empsd = df.sum[[3]]$sdGTforonset,
  emprho = df.sum[[3]]$corres.GTIPback,
  empsamplesize.ave = df.sum[[3]]$samplesize.GTIPback,
  timewindow = paste0("D", startseq, " - ", "D", endseq),
  muCIrecover = indicator.mu[[3]],
  estmu = estmu.mean[[3]],
  biasmu = bias.mu[[3]],
  widthmu = width.mu[[3]],
  sdCIrecover = indicator.sd[[3]],
  estsd = estsd.mean[[3]],
  biassd = bias.sd[[3]],
  widthsd = width.sd[[3]],
  rhoCIrecover = indicator.rho[[3]],
  estrho = estrho.mean[[3]],
  biasrho = bias.rho[[3]],
  widthrho = width.rho[[3]]
) 


df.expo7GTsim_rho_0.75 <- data.frame(
  empmean = df.sum[[4]]$meanGTforonset,
  empsd = df.sum[[4]]$sdGTforonset,
  emprho = df.sum[[4]]$corres.GTIPback,
  empsamplesize.ave = df.sum[[4]]$samplesize.GTIPback,
  timewindow = paste0("D", startseq, " - ", "D", endseq),
  muCIrecover = indicator.mu[[4]],
  estmu = estmu.mean[[4]],
  biasmu = bias.mu[[4]],
  widthmu = width.mu[[4]],
  sdCIrecover = indicator.sd[[4]],
  estsd = estsd.mean[[4]],
  biassd = bias.sd[[4]],
  widthsd = width.sd[[4]],
  rhoCIrecover = indicator.rho[[4]],
  estrho = estrho.mean[[4]],
  biasrho = bias.rho[[4]],
  widthrho = width.rho[[4]]
) 


df.simsum.expo7.med <- rbind(df.expo7GTsim_rho0, df.expo7GTsim_rho_0.25, 
                             df.expo7GTsim_rho_0.5, df.expo7GTsim_rho_0.75)
df.simsum.expo7.med$rho_simset <- rep(c(0, 0.25, 0.5, 0.75), each = 6)

write.csv(df.simsum.expo7.med, file = "sim_performance_expo7_medsetting_compensate20sim.csv")

View(df.sum[[4]])  

write.csv(df.sum[[4]], file = "empirical_simdata_rho0.75.csv")  


df.sum.emp <- rbind(df.sum[[1]], df.sum[[2]], df.sum[[3]], df.sum[[4]])
df.sum.emp$rho_simset <- rep(c(0, 0.25, 0.5, 0.75), each = 6)

write.csv(df.sum.emp, file = "empirical_simdata.csv")

