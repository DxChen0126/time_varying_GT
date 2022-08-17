library(incidence)
library(EpiEstim)
library(foreach)
library(doParallel)
library(fitdistrplus)
library(dplyr)

### my saved SImulated raw data

load("SIR_sims_list_rawdata.rda")

## simlist 1, 2, 3: long medium, short setting of GT

### my saved estimated GT data

load("GT_est_allsims.rda")

View(GT_est_total_alt)

### GT_est_total_alt 3, 6, 9: long/medium/short setting of GT with 1 day expo

### first get "true" Re estimates over my 1000 sims for time window 5-11
### R0 * N_susceptible/N_total

GetRe_fun <- function(inputpair, setting){
  # input pair is the pairdata from my SImlist
  infector.uniq <- unique(inputpair$Infector.ID)
  ninfector.uniq <- length(infector.uniq)
  infector.induniq <- numeric(ninfector.uniq)
  # get unique infector index
  for(i in 1:ninfector.uniq){
    infector.induniq[i] <- which(inputpair$Infector.ID == infector.uniq[i])[1]
  }
  # make combined epi df
  combine.epi <- data.frame(
    ID = c(inputpair$Infector.ID[infector.induniq], inputpair$Infectee.ID),
    onset = c(inputpair$infector_onsetD[infector.induniq], inputpair$infectee_onsetD),
    infect = c(inputpair$infector_TinfectD[infector.induniq], inputpair$infectee_TinfectD)
  )
  # remove duplicates
  ID.uniq <- unique(combine.epi$ID)
  ncase <- length(ID.uniq)
  ind.unique <- numeric(ncase)
  for(i in 1:ncase){
    ind.unique[i] <- which(combine.epi$ID == ID.uniq[i])[1]
  }
  unique.epi <- combine.epi[ind.unique,]
  
  # setting: 1, 2, 3, stand for long/med/short
  if(setting > 1){
  
#  Re.emp <- numeric(10)
  R0 <- 2.5
  Ntotal <- 1000
  # for(i in 1:10){
  #   # only for time window 1 - 10 for short/med GT setting 
  #   # and 1-15 for long GT setting
  #   Nsusceptible <- 1000 - length(which(unique.epi$infect <= (i)))
  #   Re.emp[i] <- R0 * Nsusceptible/Ntotal
  # }
  Nsusceptible <- 1000 - length(which(unique.epi$infect <= 10))
  Re <- R0 * Nsusceptible/Ntotal
  }
  
  if(setting == 1){
    
#    Re.emp <- numeric(15)
    R0 <- 2.5
    Ntotal <- 1000
    # for(i in 1:15){
    #   # only for time window 1 - 10 for short/med GT setting 
    #   # and 1-15 for long GT setting
    #   Nsusceptible <- 1000 - length(which(unique.epi$infect <= (i)))
    #   Re.emp[i] <- R0 * Nsusceptible/Ntotal
    # }
    Nsusceptible <- 1000 - length(which(unique.epi$infect <= 15))
    Re <- R0 * Nsusceptible/Ntotal
  }
  
  
  # get geometric mean of Rt in this time window
#  output <- exp(mean(log(Re.emp)))
  output <- Re
  return(output)
}

# test
GetRe_fun(simlist[[3]][3]$pairdata, setting = 3) # works

### Get Re at first time window for each GT setting

Re.list <- vector("list", 3)

ind.pairdata <- seq(3, 3000, 3)

for(j in 1:3){

pairdatalist <- simlist[[j]][ind.pairdata]
ncores <- detectCores()
myCluster <- makeCluster(ncores - 2, # number of cores to use
                         type = "PSOCK") # type of cluster
nsim <- 1000
Re.list[[j]] <- foreach(i = 1:nsim, .combine = 'append') %dopar% {
  GetRe_fun(pairdatalist[[i]], setting = j)
}
stopCluster(myCluster)
}

Re.list 

save(Re.list, file = "Re_list_raw100sim_0803.rda")

Re.ref <- unlist(lapply(Re.list, mean))

### estimate SI in 50 sims in first time window 
### for each GT setting and save the results

# make a function to estimate SI

estSI <- function(inputpair, setting){
  muSI <- sdSI <- c()
  start <- end <- c()
  
  if(setting > 1){
    start = 1
    end = 10
  }
  
  if(setting == 1){
    start = 1
    end = 15
  }
  
    ind.tmp <- which(inputpair$infector_onsetD >= start & 
                       inputpair$infector_onsetD <= end)
    tmppair <- inputpair[ind.tmp,]
    SItmp <- tmppair$infectee_onsetD - tmppair$infector_onsetD
    estSItmp <- fitdist(SItmp, "norm")
    muSI <- estSItmp$estimate[1]
    sdSI <- estSItmp$estimate[2]
  
  return(data.frame(
    mean = muSI,
    sd = sdSI
  ))
}

# test
# estSI(test, setting = 2) # works

estSIlist <- vector("list", 3)

nselect <- 50

for(j in 1:3){

progressbar <- txtProgressBar(min = 0, max = 50, style = 3)

estimatedSI <- vector("list", nselect)

for(n in 1:nselect){
  set.seed(n)
  ind.sample <- sample(seq(1, 1000, 1), 1, replace = F)
  simpair <- simlist[[j]][ind.pairdata][ind.sample]$pairdata
  estimatedSI[[n]] <- estSI(simpair, setting = j)
  setTxtProgressBar(progressbar, n)
}

estSIlist[[j]] <- estimatedSI

}

save(estSIlist, file = "SI_est_allsims_allGTsetting_0803.rda")



###
# make a function to get incidence
GetIncidence_fun <- function(inputpair){
  # input pair is the pairdata from my simlist
  infector.uniq <- unique(inputpair$Infector.ID)
  ninfector.uniq <- length(infector.uniq)
  infector.induniq <- numeric(ninfector.uniq)
  # get unique infector index
  for(i in 1:ninfector.uniq){
    infector.induniq[i] <- which(inputpair$Infector.ID == infector.uniq[i])[1]
  }
  # make combined epi df
  combine.epi <- data.frame(
    ID = c(inputpair$Infector.ID[infector.induniq], inputpair$Infectee.ID),
    onset = c(inputpair$infector_onsetD[infector.induniq], inputpair$infectee_onsetD)
  )
  # remove duplicates
  ID.uniq <- unique(combine.epi$ID)
  ncase <- length(ID.uniq)
  ind.unique <- numeric(ncase)
  for(i in 1:ncase){
    ind.unique[i] <- which(combine.epi$ID == ID.uniq[i])[1]
  }
  unique.epi <- combine.epi[ind.unique,]
  epiincidence <- incidence(unique.epi$onset)
  return(epiincidence)
}

# test function
# GetIncidence_fun(test) # works

# just check Re obtained in time window 5-11

# R estimated based on SI

Rest.SI.list <- vector("list", 3)

for(j in 1:3){
  
Rest.SI <- RLB.SI <- RUB.SI <- numeric(nselect)

progressbar <- txtProgressBar(min = 0, max = 50, style = 3)

for(n in 1:nselect){
  set.seed(n)
  ind.sample <- sample(seq(1, 1000, 1), 1, replace = F)
  simpair <- simlist[[j]][ind.pairdata][ind.sample]$pairdata
  # get incidence
  tmpincidence <- GetIncidence_fun(simpair)
  estSI <- estSIlist[[j]][[n]]
  
  if(j > 1){

  R.SI <- wallinga_teunis(tmpincidence, method="parametric_si",
                          config = list(t_start = seq(1, 20), t_end = seq(10, 29),
                            mean_si = estSI$mean, std_si = estSI$sd,
                            n_sim = 100))
    Rest.SI[n] <- R.SI$R$`Mean(R)`[1]
    RLB.SI[n] <- R.SI$R$`Quantile.0.025(R)`[1]
    RUB.SI[n] <- R.SI$R$`Quantile.0.975(R)`[1]
  }
  
  if(j == 1){
    
    R.SI <- wallinga_teunis(tmpincidence, method="parametric_si",
                            config = list(t_start = seq(1, 20), t_end = seq(15, 34),
                                          mean_si = estSI$mean, std_si = estSI$sd,
                                          n_sim = 100))
    Rest.SI[n] <- R.SI$R$`Mean(R)`[1]
    RLB.SI[n] <- R.SI$R$`Quantile.0.025(R)`[1]
    RUB.SI[n] <- R.SI$R$`Quantile.0.975(R)`[1]
  }
  
  setTxtProgressBar(progressbar, n)
}

Rest.SI.df <- data.frame(
  est = Rest.SI,
  LB = RLB.SI,
  UB = RUB.SI
)

Rest.SI.list[[j]] <- Rest.SI.df
}

save(Rest.SI.list, file = "R_est_bySI_allGTsetting_0803.rda")

# R estimated based on GT

Rest.GT.list <- vector("list", 3)

for(j in 1:3){
  
  Rest.GT <- RLB.GT <- RUB.GT <- numeric(nselect)
  
  progressbar <- txtProgressBar(min = 0, max = 50, style = 3)
  
  for(n in 1:nselect){
    set.seed(n)
    ind.sample <- sample(seq(1, 1000, 1), 1, replace = F)
    simpair <- simlist[[j]][ind.pairdata][ind.sample]$pairdata
    # get incidence
    tmpincidence <- GetIncidence_fun(simpair)
    estGT <- GT_est_total_alt[[3 + (j-1)*3]][[n]]
    
    if(j > 1){
    
    R.GT <- wallinga_teunis(tmpincidence, method="parametric_si",
                            config = list(t_start = seq(1, 20), t_end = seq(10, 29),
                                          mean_si = estGT$mu[1], std_si = estGT$sd[1],
                                          n_sim = 100))
    Rest.GT[n] <- R.GT$R$`Mean(R)`[1]
    RLB.GT[n] <- R.GT$R$`Quantile.0.025(R)`[1]
    RUB.GT[n] <- R.GT$R$`Quantile.0.975(R)`[1]
    }
    
    if(j == 1){
      
      R.GT <- wallinga_teunis(tmpincidence, method="parametric_si",
                              config = list(t_start = seq(1, 20), t_end = seq(15, 34),
                                            mean_si = estGT$mu[1], std_si = estGT$sd[1],
                                            n_sim = 100))
      Rest.GT[n] <- R.GT$R$`Mean(R)`[1]
      RLB.GT[n] <- R.GT$R$`Quantile.0.025(R)`[1]
      RUB.GT[n] <- R.GT$R$`Quantile.0.975(R)`[1]
    }
    
    
    setTxtProgressBar(progressbar, n)
  }
  
  Rest.GT.df <- data.frame(
    est = Rest.GT,
    LB = RLB.GT,
    UB = RUB.GT
  )
  
  Rest.GT.list[[j]] <- Rest.GT.df
}

save(Rest.GT.list, file = "R_est_byGT_allGTsetting_0803.rda")

### check sim performance

Re.ref

Rest.SI.list

R.SI.simperform <- vector("list", 3)

for(i in 1:3){
  Rref <- Re.ref[i]
  RestSI <- Rest.SI.list[[i]]
  ind <- 0
  for(j in 1:50){
    if(Rref >= RestSI$LB[j] & Rref <= RestSI$UB[j]){
      ind <- ind + 1
    }
  }
  estmean <- mean(RestSI$est)
  bias <- (estmean - Rref)/Rref
  outdf <- data.frame(
    Rest = estmean,
    RLB = mean(RestSI$LB),
    RUB = mean(RestSI$UB),
    CIcover = ind,
    Bias = bias,
    Rref = Rref
  )
  R.SI.simperform[[i]] <- outdf
}

R.SI.simperform 

R.GT.simperform <- vector("list", 3)

for(i in 1:3){
  Rref <- Re.ref[i]
  RestGT <- Rest.GT.list[[i]]
  ind <- 0
  for(j in 1:50){
    if(Rref >= RestGT$LB[j] & Rref <= RestGT$UB[j]){
      ind <- ind + 1
    }
  }
  estmean <- mean(RestGT$est)
  bias <- (estmean - Rref)/Rref
  outdf <- data.frame(
    Rest = estmean,
    RLB = mean(RestGT$LB),
    RUB = mean(RestGT$UB),
    CIcover = ind,
    Bias = bias,
    Rref = Rref
  )
  R.GT.simperform[[i]] <- outdf
}

R.GT.simperform 

R.SI.simperform 

simperform.combine <- append(R.GT.simperform, R.SI.simperform)

write.csv(simperform.combine, "simperform.combine_R0_0803.csv")


######## get est SI and empirical SI summary for the first window

empSI <- function(inputpair, setting){
  muSI <- sdSI <- c()
  start <- end <- c()
  
  if(setting > 1){
    start = 1
    end = 10
  }
  
  if(setting == 1){
    start = 1
    end = 15
  }
  
  ind.tmp <- which(inputpair$infector_onsetD >= start & 
                     inputpair$infector_onsetD <= end)
  tmppair <- inputpair[ind.tmp,]
  SItmp <- tmppair$infectee_onsetD - tmppair$infector_onsetD
  muSI <- mean(SItmp)
  sdSI <- sd(SItmp)
  
  return(c(muSI, sdSI))
}

# test
# estSI(test, setting = 2) # works

empSIall <- matrix(nrow = 3, ncol = 2)


for(j in 1:3){
  
  empSImat <- matrix(nrow = 1000, ncol = 2)
  
  for(i in 1:1000){
    simpair <- simlist[[j]][ind.pairdata][i]$pairdata
    empSImat[i,] <- empSI(simpair, setting = j)
  }
  
  empSIall[j,] <- colMeans(empSImat)
  
}

save(empSIall, file = "SI_emp_allsims_allGTsetting_0803.rda")

empSIall <- as.data.frame(empSIall)
colnames(empSIall) <- c("mu", "sd")
write.csv(empSIall, "empSIsum.1stwindow.csv")
# summary over est SI
estSIlist

estSIsum <- vector("list", 3)

for(i in 1:3){
  SIest <- estSIlist[[i]]
  SIestmat <- matrix(unlist(SIest), nrow = 50, ncol = 2, byrow = T)
  SI.mu <- mean(SIestmat[,1])
  SI.sd <- mean(SIestmat[,2])
  estSIsum[[i]] <- c(SI.mu, SI.sd)
}

estSIsum <- matrix(unlist(estSIsum), nrow = 3, byrow = T)
estSIsum <- as.data.frame(estSIsum)
colnames(estSIsum) <- c("mu", "sd")
write.csv(estSIsum, "estSIsum.1stwindow.csv")
### R0 based on GT with different expo 4, 7, 10

# 4 day
load("GT_est_addsims_0408.rda")

# R estimated based on GT

Rest.GT.list.expo4 <- vector("list", 3)

for(j in 1:3){
  
  Rest.GT <- RLB.GT <- RUB.GT <- numeric(nselect)
  
  progressbar <- txtProgressBar(min = 0, max = 50, style = 3)
  
  for(n in 1:nselect){
    set.seed(n)
    ind.sample <- sample(seq(1, 1000, 1), 1, replace = F)
    simpair <- simlist[[j]][ind.pairdata][ind.sample]$pairdata
    # get incidence
    tmpincidence <- GetIncidence_fun(simpair)
    estGT <- GT_est_total[[j]][[n]]
    
    if(j > 1){
      
      R.GT <- wallinga_teunis(tmpincidence, method="parametric_si",
                              config = list(t_start = seq(1, 20), t_end = seq(10, 29),
                                            mean_si = estGT$mu[1], std_si = estGT$sd[1],
                                            n_sim = 100))
      Rest.GT[n] <- R.GT$R$`Mean(R)`[1]
      RLB.GT[n] <- R.GT$R$`Quantile.0.025(R)`[1]
      RUB.GT[n] <- R.GT$R$`Quantile.0.975(R)`[1]
    }
    
    if(j == 1){
      
      R.GT <- wallinga_teunis(tmpincidence, method="parametric_si",
                              config = list(t_start = seq(1, 20), t_end = seq(15, 34),
                                            mean_si = estGT$mu[1], std_si = estGT$sd[1],
                                            n_sim = 100))
      Rest.GT[n] <- R.GT$R$`Mean(R)`[1]
      RLB.GT[n] <- R.GT$R$`Quantile.0.025(R)`[1]
      RUB.GT[n] <- R.GT$R$`Quantile.0.975(R)`[1]
    }
    
    
    setTxtProgressBar(progressbar, n)
  }
  
  Rest.GT.df <- data.frame(
    est = Rest.GT,
    LB = RLB.GT,
    UB = RUB.GT
  )
  
  Rest.GT.list.expo4[[j]] <- Rest.GT.df
}

Rest.GT.list.expo4

save(Rest.GT.list.expo4, file = "R_est_byGTexpo4_allGTsetting_0803.rda")


R.GT.simperform.expo4 <- vector("list", 3)

for(i in 1:3){
  Rref <- Re.ref[i]
  RestGT <- Rest.GT.list.expo4[[i]]
  ind <- 0
  for(j in 1:50){
    if(Rref >= RestGT$LB[j] & Rref <= RestGT$UB[j]){
      ind <- ind + 1
    }
  }
  estmean <- mean(RestGT$est)
  bias <- (estmean - Rref)/Rref
  outdf <- data.frame(
    Rest = estmean,
    RLB = mean(RestGT$LB),
    RUB = mean(RestGT$UB),
    CIcover = ind,
    Bias = bias,
    Rref = Rref
  )
  R.GT.simperform.expo4[[i]] <- outdf
}

R.GT.simperform.expo4 

View(GT_est_total_alt[[5]][[33]])

# 7 day

Rest.GT.list.expo7 <- vector("list", 3)

for(j in 1:3){
  
  Rest.GT <- RLB.GT <- RUB.GT <- numeric(nselect)
  
  progressbar <- txtProgressBar(min = 0, max = 50, style = 3)
  
  for(n in 1:nselect){
    set.seed(n)
    ind.sample <- sample(seq(1, 1000, 1), 1, replace = F)
    simpair <- simlist[[j]][ind.pairdata][ind.sample]$pairdata
    # get incidence
    tmpincidence <- GetIncidence_fun(simpair)
    estGT <- GT_est_total_alt[[2 + (j-1)*3]][[n]]
    
    if(j > 1){
      
      if(estGT$mu[1] < 1){
        estGT$mu[1] <- estGT$mu[2]
        estGT$sd[1] <- estGT$sd[2]
      }
      
      R.GT <- wallinga_teunis(tmpincidence, method="parametric_si",
                              config = list(t_start = seq(1, 20), t_end = seq(10, 29),
                                            mean_si = estGT$mu[1], std_si = estGT$sd[1],
                                            n_sim = 100))
      Rest.GT[n] <- R.GT$R$`Mean(R)`[1]
      RLB.GT[n] <- R.GT$R$`Quantile.0.025(R)`[1]
      RUB.GT[n] <- R.GT$R$`Quantile.0.975(R)`[1]
    }
    
    if(j == 1){
      if(estGT$mu[1] < 1){
        estGT$mu[1] <- estGT$mu[2]
        estGT$sd[1] <- estGT$sd[2]
      }
      
      R.GT <- wallinga_teunis(tmpincidence, method="parametric_si",
                              config = list(t_start = seq(1, 20), t_end = seq(15, 34),
                                            mean_si = estGT$mu[1], std_si = estGT$sd[1],
                                            n_sim = 100))
      Rest.GT[n] <- R.GT$R$`Mean(R)`[1]
      RLB.GT[n] <- R.GT$R$`Quantile.0.025(R)`[1]
      RUB.GT[n] <- R.GT$R$`Quantile.0.975(R)`[1]
    }
    
    
    setTxtProgressBar(progressbar, n)
 #   print(j)
  }
  
  Rest.GT.df <- data.frame(
    est = Rest.GT,
    LB = RLB.GT,
    UB = RUB.GT
  )
  
  Rest.GT.list.expo7[[j]] <- Rest.GT.df
}

Rest.GT.list.expo7

save(Rest.GT.list.expo7, file = "R_est_byGTexpo7_allGTsetting_0803.rda")


R.GT.simperform.expo7 <- vector("list", 3)

for(i in 1:3){
  Rref <- Re.ref[i]
  RestGT <- Rest.GT.list.expo7[[i]]
  ind <- 0
  for(j in 1:50){
    if(Rref >= RestGT$LB[j] & Rref <= RestGT$UB[j]){
      ind <- ind + 1
    }
  }
  estmean <- mean(RestGT$est)
  bias <- (estmean - Rref)/Rref
  outdf <- data.frame(
    Rest = estmean,
    RLB = mean(RestGT$LB),
    RUB = mean(RestGT$UB),
    CIcover = ind,
    Bias = bias,
    Rref = Rref
  )
  R.GT.simperform.expo7[[i]] <- outdf
}

R.GT.simperform.expo7 

# 10 day

load("GT_est_addsims_expo10_0411.rda")

View(GT_est_total[[1]][[1]])

Rest.GT.list.expo10 <- vector("list", 3)

for(j in 1:3){
  
  Rest.GT <- RLB.GT <- RUB.GT <- numeric(nselect)
  
  progressbar <- txtProgressBar(min = 0, max = 50, style = 3)
  
  for(n in 1:nselect){
    set.seed(n)
    ind.sample <- sample(seq(1, 1000, 1), 1, replace = F)
    simpair <- simlist[[j]][ind.pairdata][ind.sample]$pairdata
    # get incidence
    tmpincidence <- GetIncidence_fun(simpair)
    estGT <- GT_est_total[[j]][[n]]
    
    if(j > 1){
      
        if(estGT$mu[1] < 1){
          estGT$mu[1] <- estGT$mu[2]
          estGT$sd[1] <- estGT$sd[2]
        }
      
      R.GT <- wallinga_teunis(tmpincidence, method="parametric_si",
                              config = list(t_start = seq(1, 20), t_end = seq(10, 29),
                                            mean_si = estGT$mu[1], std_si = estGT$sd[1],
                                            n_sim = 100))
      Rest.GT[n] <- R.GT$R$`Mean(R)`[1]
      RLB.GT[n] <- R.GT$R$`Quantile.0.025(R)`[1]
      RUB.GT[n] <- R.GT$R$`Quantile.0.975(R)`[1]
    }
    
    if(j == 1){
      
      if(estGT$mu[1] < 1){
        estGT$mu[1] <- estGT$mu[2]
        estGT$sd[1] <- estGT$sd[2]
      }
      
      R.GT <- wallinga_teunis(tmpincidence, method="parametric_si",
                              config = list(t_start = seq(1, 20), t_end = seq(15, 34),
                                            mean_si = estGT$mu[1], std_si = estGT$sd[1],
                                            n_sim = 100))
      Rest.GT[n] <- R.GT$R$`Mean(R)`[1]
      RLB.GT[n] <- R.GT$R$`Quantile.0.025(R)`[1]
      RUB.GT[n] <- R.GT$R$`Quantile.0.975(R)`[1]
    }
    
    
    setTxtProgressBar(progressbar, n)
  }
  
  Rest.GT.df <- data.frame(
    est = Rest.GT,
    LB = RLB.GT,
    UB = RUB.GT
  )
  
  Rest.GT.list.expo10[[j]] <- Rest.GT.df
}

Rest.GT.list.expo10

save(Rest.GT.list.expo10, file = "R_est_byGTexpo10_allGTsetting_0803.rda")


R.GT.simperform.expo10 <- vector("list", 3)

for(i in 1:3){
  Rref <- Re.ref[i]
  RestGT <- Rest.GT.list.expo10[[i]]
  ind <- 0
  for(j in 1:50){
    if(Rref >= RestGT$LB[j] & Rref <= RestGT$UB[j]){
      ind <- ind + 1
    }
  }
  estmean <- mean(RestGT$est)
  bias <- (estmean - Rref)/Rref
  outdf <- data.frame(
    Rest = estmean,
    RLB = mean(RestGT$LB),
    RUB = mean(RestGT$UB),
    CIcover = ind,
    Bias = bias,
    Rref = Rref
  )
  R.GT.simperform.expo10[[i]] <- outdf
}

R.GT.simperform.expo10 

R.SI.mtx <- matrix(unlist(R.SI.simperform), nrow = 3, byrow = T)
R.GT.mtx.expo1 <- matrix(unlist(R.GT.simperform), nrow = 3, byrow = T)
R.GT.mtx.expo4 <- matrix(unlist(R.GT.simperform.expo4), nrow = 3, byrow = T)
R.GT.mtx.expo7 <- matrix(unlist(R.GT.simperform.expo7), nrow = 3, byrow = T)
R.GT.mtx.expo10 <- matrix(unlist(R.GT.simperform.expo10), nrow = 3, byrow = T)

simperform.sum <- rbind(R.SI.mtx,
                        R.GT.mtx.expo1,
                        R.GT.mtx.expo4,
                        R.GT.mtx.expo7,
                        R.GT.mtx.expo10)

simperform.sum <- as.data.frame(simperform.sum)

colnames(simperform.sum) <- c("R.est", "R.LB", "R.UB", "CIcover", "Bias", "R.ref")
simperform.sum$based_on <- c(rep("SI", 3), rep("GT_expo1", 3), rep("GT_expo4", 3),
                             rep("GT_expo7", 3), rep("GT_expo10", 3))
simperform.sum$GT_setting <- rep(c("long", "med", "short"), 5)
simperform.sum$timewindow <- rep(c("D1-15", "D1-10", "D1-10"), 5)


write.csv(simperform.sum, "simR0.perform_allexpo_0803.csv")












