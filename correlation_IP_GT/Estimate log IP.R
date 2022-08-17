rm(list = ls())

allpairs <- read.csv("allpairs0908.csv")

library(fitdistrplus)
library(tidyverse)

# format date

colexpo <- grep(c("expo"), colnames(allpairs)) 
for(i in colexpo){
  allpairs[, i] <- as.Date(allpairs[, i], "%m/%d/%Y")
}

colonset <- grep(c("Onset"), colnames(allpairs)) 
for(i in colonset){
  allpairs[, i] <- as.Date(allpairs[, i], "%m/%d/%Y")
}

#################(weighted) Incubation period of infector and infectee (referenced by infector onset)

##### infector
df.infector.full <- allpairs %>% drop_na(expo.early.I, expo.late.I)

# check if there is any infector who has earliest possible exposure and onset at same date
which(df.infector.full$expo.early.I == df.infector.full$Onset_Infector) # 111, 112
df.infector.full[c(111, 112),]$expo.early.I <- df.infector.full[c(111, 112),]$Onset_Infector - 1

startseq <- c(as.Date("2020-01-01"), seq.Date(from = as.Date("2020-01-15"), to = as.Date("2020-01-30"), by = "1 day"))
endseq <- c(as.Date("2020-01-20"), seq.Date(from = as.Date("2020-01-21"), to = as.Date("2020-02-04"), by = "1 day"),
            as.Date("2020-02-29"))

# temporal estimates of IP (weighted) infector
n <- length(startseq)

par1 <- par2 <- numeric(n)
par1.boot <- par2.boot <- vector("list", n)
samplesize <- numeric(n)

progressbar <- txtProgressBar(min = 0, max = n, style = 3)

for(i in 1:n){
  # i <- 1
  subdata <- subset(df.infector.full, Onset_Infector >= startseq[i] & Onset_Infector <= endseq[i])
  samplesize[i] <- nrow(subdata)
  dftmp <- data.frame(left = with(subdata, as.numeric(Onset_Infector - expo.late.I)),
                      right = with(subdata, as.numeric(Onset_Infector - expo.early.I)))
  
  fit.ln <- fitdistcens(dftmp, "lnorm")
  fit.ln.boot <- bootdistcens(fit.ln, niter = 1000)
  
  par1[i] <- fit.ln$estimate[1] # meanlog
  par2[i] <- fit.ln$estimate[2] # sdlog
  
  par1.boot[[i]] <- fit.ln.boot$estim[, 1]
  par2.boot[[i]] <- fit.ln.boot$estim[, 2]
  
  setTxtProgressBar(progressbar, i)
}

close(progressbar)

df.IPlnormpars.weightinfector <- data.frame(
  meanlog = par1,
  sdlog = par2,
  timewindow = paste(startseq, endseq),
  samplesize = samplesize
)

save(df.IPlnormpars.weightinfector, file = "IPlnormpars_weightinfector.rda")

muIP <- sdIP <- numeric(n)
muIP.lb <- muIP.ub <- sdIP.lb <- sdIP.ub <- numeric(n)

for(i in 1:n){
  muIP[i] <- exp(par1[i] + (1/2)*par2[i]^2)
  sdIP[i] <- exp(par1[i] + (1/2)*par2[i]^2)*sqrt(exp(par2[i]^2) - 1)
  muIP.lb[i] <- quantile(exp(par1.boot[[i]] + (1/2)*par2.boot[[i]]^2), 0.025)
  muIP.ub[i] <- quantile(exp(par1.boot[[i]] + (1/2)*par2.boot[[i]]^2), 0.975)
  sdIP.lb[i] <- quantile(exp(par1.boot[[i]] + (1/2)*par2.boot[[i]]^2)*sqrt(exp(par2.boot[[i]]^2) - 1), 0.025)
  sdIP.ub[i] <- quantile(exp(par1.boot[[i]] + (1/2)*par2.boot[[i]]^2)*sqrt(exp(par2.boot[[i]]^2) - 1), 0.975)
}

df.IPlnormest.weightinfector <- data.frame(
  mean = muIP,
  mean.lb = muIP.lb,
  mean.ub = muIP.ub,
  sd = sdIP,
  sd.lb = sdIP.lb,
  sd.ub = sdIP.ub,
  timewindow = paste(startseq, endseq),
  samplesize = samplesize
)

write.csv(df.IPlnormest.weightinfector, "IPlnormest_weightinfector.csv")

########### infectee

df.infectee.full <- allpairs %>% drop_na(expo.early.S, expo.late.S)

# check if there is any infector who has earliest possible exposure and onset at same date
which(df.infectee.full$expo.early.S == df.infectee.full$Onset_Infectee) # none

n <- length(startseq)

par1 <- par2 <- numeric(n)
par1.boot <- par2.boot <- vector("list", n)
samplesize <- numeric(n)

progressbar <- txtProgressBar(min = 0, max = n, style = 3)

for(i in 1:n){
  # i <- 1
  subdata <- subset(df.infectee.full, Onset_Infector >= startseq[i] & Onset_Infector <= endseq[i])
  samplesize[i] <- nrow(subdata)
  dftmp <- data.frame(left = with(subdata, as.numeric(Onset_Infectee - expo.late.S)),
                      right = with(subdata, as.numeric(Onset_Infectee - expo.early.S)))
  
  fit.ln <- fitdistcens(dftmp, "lnorm")
  fit.ln.boot <- bootdistcens(fit.ln, niter = 1000)
  
  par1[i] <- fit.ln$estimate[1] # meanlog
  par2[i] <- fit.ln$estimate[2] # sdlog
  
  par1.boot[[i]] <- fit.ln.boot$estim[, 1]
  par2.boot[[i]] <- fit.ln.boot$estim[, 2]
  
  setTxtProgressBar(progressbar, i)
}

close(progressbar)

df.IPlnormpars.infectee <- data.frame(
  meanlog = par1,
  sdlog = par2,
  timewindow = paste(startseq, endseq),
  samplesize = samplesize
)

save(df.IPlnormpars.infectee, file = "IPlnormpars_infectee.rda")

muIP <- sdIP <- numeric(n)
muIP.lb <- muIP.ub <- sdIP.lb <- sdIP.ub <- numeric(n)

for(i in 1:n){
  muIP[i] <- exp(par1[i] + (1/2)*par2[i]^2)
  sdIP[i] <- exp(par1[i] + (1/2)*par2[i]^2)*sqrt(exp(par2[i]^2) - 1)
  muIP.lb[i] <- quantile(exp(par1.boot[[i]] + (1/2)*par2.boot[[i]]^2), 0.025)
  muIP.ub[i] <- quantile(exp(par1.boot[[i]] + (1/2)*par2.boot[[i]]^2), 0.975)
  sdIP.lb[i] <- quantile(exp(par1.boot[[i]] + (1/2)*par2.boot[[i]]^2)*sqrt(exp(par2.boot[[i]]^2) - 1), 0.025)
  sdIP.ub[i] <- quantile(exp(par1.boot[[i]] + (1/2)*par2.boot[[i]]^2)*sqrt(exp(par2.boot[[i]]^2) - 1), 0.975)
}

df.IPlnormest.infectee <- data.frame(
  mean = muIP,
  mean.lb = muIP.lb,
  mean.ub = muIP.ub,
  sd = sdIP,
  sd.lb = sdIP.lb,
  sd.ub = sdIP.ub,
  timewindow = paste(startseq, endseq),
  samplesize = samplesize
)

write.csv(df.IPlnormest.infectee, "IPlnormest_infectee.csv")


