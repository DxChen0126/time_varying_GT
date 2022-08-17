source("sir_corr.R")
library(dplyr)

nsim <- 1000

corr <- c(0, 0.25, 0.5, 0.75)

simlist_corr <- vector('list', length(corr))

for (j in 1:length(corr)) {
  print(j)
  simlist <- vector('list', nsim)
  
  i <- 1
  
  while (i <= nsim) {
    print(i)
    sir_sim <- sir.full2(size=1000, I0=10, seed=i, rho=corr[j], keep.intrinsic = FALSE)
    
    if (nrow(sir_sim$data) > 100) {
      simlist[[i]] <- sir_sim
      i <- i +1
    }
  }
  
  simlist_corr[[j]] <- simlist
}

save(simlist_corr, file="sir_sim_corr_0628.rda")

load("sir_sim_corr_0628.rda")

startseq <- seq(5, 30, 5)
endseq <- seq(10, 35, 5)

df.sum <- vector("list", 4)

for(j in 1:4){
  corlist_GTIPback <- vector("list", 1000)
  corlist_GTIPfor <- vector("list", 1000)
  samplesizelist.1 <- vector("list", 1000)
  samplesizelist.2 <- vector("list", 1000)
  GTlistfor.m <- GTlistforonset.m <- GTlistfor.sd <- GTlistforonset.sd <- vector("list", 1000)
  IPlistfor.m <- IPlistback.m <- IPlistfor.sd <- IPlistback.sd <- vector("list", 1000)
  
  simrho_j <- simlist_corr[[j]]
  
  pb <- txtProgressBar(min = 0,      
                       max = 1000, 
                       style = 3)  
  
  for(i in 1:1000){
    simdata <- simrho_j[[i]]
    
    Infector.ID = simdata$infected_by[!is.na(simdata$infected_by)]
    Infectee.ID = which(!is.na(simdata$infected_by))
    infector_Tinfect = simdata$t_infected[simdata$infected_by[!is.na(simdata$infected_by)]]
    infector_Tonset = simdata$t_symptomatic[simdata$infected_by[!is.na(simdata$infected_by)]]
    infectee_Tinfect = simdata$t_infected[which(!is.na(simdata$infected_by))]
    infectee_Tonset = simdata$t_symptomatic[which(!is.na(simdata$infected_by))]
    
    
    df.pair <- data.frame(
      infectorID = Infector.ID,
      infectorTinfect = infector_Tinfect,
      infectorTonset = infector_Tonset,
      infecteeID = Infectee.ID,
      infecteeTinfect = infectee_Tinfect,
      infecteeTonset = infectee_Tonset
    )
    
    df.pair <- df.pair %>% mutate(
      GT = infecteeTinfect - infectorTinfect,
      infectorIP = infectorTonset - infectorTinfect
    )
    
    corlist_GTIPback[[i]] <- corlist_GTIPfor [[i]] <- samplesizelist.1[[i]] <- samplesizelist.2[[i]] <- numeric(6)
    GTlistfor.m[[i]] <- GTlistforonset.m[[i]] <- GTlistfor.sd[[i]] <- GTlistforonset.sd[[i]] <- numeric(6)
    IPlistfor.m[[i]] <- IPlistback.m[[i]] <- IPlistfor.sd[[i]] <- IPlistback.sd[[i]] <- numeric(6)
    
    
    for(m in 1:6){
      tmppair1 <- subset(df.pair, startseq[m] <= infectorTinfect & infectorTinfect <= endseq[m])
      corlist_GTIPfor[[i]][m] <- cor(tmppair1$GT, tmppair1$infectorIP) 
      samplesizelist.1[[i]][m] <- nrow(tmppair1)
      GTlistfor.m[[i]][m] <- mean(tmppair1$GT)
      GTlistfor.sd[[i]][m] <- sd(tmppair1$GT)
      IPlistfor.m[[i]][m] <- mean(tmppair1$infectorIP)
      IPlistfor.sd[[i]][m] <- sd(tmppair1$infectorIP)
      
      tmppair2 <- subset(df.pair, startseq[m] <= infectorTonset & infectorTonset <= endseq[m])
      corlist_GTIPback[[i]][m] <- cor(tmppair2$GT, tmppair2$infectorIP) 
      samplesizelist.2[[i]][m] <- nrow(tmppair2)
      GTlistforonset.m[[i]][m] <- mean(tmppair2$GT)
      GTlistforonset.sd[[i]][m] <- sd(tmppair2$GT)
      IPlistback.m[[i]][m] <- mean(tmppair2$infectorIP)
      IPlistback.sd[[i]][m] <- sd(tmppair2$infectorIP)
    }
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  corres.GTIPfor <- colMeans(matrix(unlist(corlist_GTIPfor), 
                                    ncol = 6, nrow = 1000, byrow = T))
  corres.GTIPback <- colMeans(matrix(unlist(corlist_GTIPback), 
                                     ncol = 6, nrow = 1000, byrow = T))
  samplesize.1 <- colMeans(matrix(unlist(samplesizelist.1), 
                                  ncol = 6, nrow = 1000, byrow = T))
  samplesize.2 <- colMeans(matrix(unlist(samplesizelist.2), 
                                  ncol = 6, nrow = 1000, byrow = T))
  meanGTfor <- colMeans(matrix(unlist(GTlistfor.m), ncol = 6, nrow = 1000, byrow = T))
  sdGTfor <- colMeans(matrix(unlist(GTlistfor.sd), ncol = 6, nrow = 1000, byrow = T))
  meanIPfor <- colMeans(matrix(unlist(IPlistfor.m), ncol = 6, nrow = 1000, byrow = T))
  sdIPfor <- colMeans(matrix(unlist(IPlistfor.sd), ncol = 6, nrow = 1000, byrow = T))
  meanGTforonset <- colMeans(matrix(unlist(GTlistforonset.m), ncol = 6, nrow = 1000, byrow = T))
  sdGTforonset <- colMeans(matrix(unlist(GTlistforonset.sd), ncol = 6, nrow = 1000, byrow = T))
  meanIPback <- colMeans(matrix(unlist(IPlistback.m), ncol = 6, nrow = 1000, byrow = T))
  sdIPback <- colMeans(matrix(unlist(IPlistback.sd), ncol = 6, nrow = 1000, byrow = T))
  
  df.sum[[j]] <- data.frame(
    corres.GTIPfor = corres.GTIPfor,
    samplesize.GTIPfor = samplesize.1,
    meanGTfor = meanGTfor,
    sdGTfor = sdGTfor,
    meanIPfor = meanIPfor,
    sdIPfor = sdIPfor,
    corres.GTIPback = corres.GTIPback, 
    samplesize.GTIPback = samplesize.2,
    meanGTforonset = meanGTforonset,
    sdGTforonset = sdGTforonset,
    meanIPback = meanIPback,
    sdIPback = sdIPback
  )
  
}


save(df.sum, file = "empirical_corr_ref.rda")

### process sim

