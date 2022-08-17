rm(list = ls())
library(fitdistrplus)
library(foreach)
library(doParallel)
library(mixdist)
source("GT_sampling_functions.R")


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


###### estimate IP in sliding windows without distinguishing infector-infectee
# first make sure all cases earliest expo should be at least one day earlier than onset
ind.notice <- which(allpairs$expo.early.I == allpairs$Onset_Infector)
allpairs[ind.notice,]$expo.early.I <- allpairs[ind.notice,]$Onset_Infector - 1

ind.notice2 <- which(allpairs$expo.early.S == allpairs$Onset_Infectee)
# this is empty

startseq <- c(as.Date("2020-01-01"), seq.Date(from = as.Date("2020-01-15"), to = as.Date("2020-01-30"), by = "1 day"))
endseq <- c(as.Date("2020-01-20"), seq.Date(from = as.Date("2020-01-21"), to = as.Date("2020-02-04"), by = "1 day"),
            as.Date("2020-02-29"))

# temporal estimates of IP
n <- length(startseq)

par1 <- par2 <- vector("list", n)
par1.boot <- par2.boot <- vector("list", n)
AIC <- vector("list", n)
samplesize <- numeric(n)

progressbar <- txtProgressBar(min = 0, max = n, style = 3)

for(i in 1:n){
  subdata1 <- subset(allpairs, Onset_Infector >= startseq[i] & Onset_Infector <= endseq[i])
  subdata2 <- subset(subdata1, !is.na(expo.early.I) & !is.na(expo.late.I))
  subdata3 <- subset(subdata1, !is.na(expo.early.S) & !is.na(expo.late.S))
  # get unique infector
  uniqueinfector <- unique(subdata2$Infector.ID)
  ind.infector.uniq <- c()
  for(a in 1:length(uniqueinfector)){
    ind.infector.uniq <- c(ind.infector.uniq, 
                           which(subdata2$Infector.ID == uniqueinfector[a])[1])
  }
  subdata2 <- subdata2[ind.infector.uniq,]
  # get unique infecteee (if a case is both infector and infectee at same time window)
  ind.infector <- subdata2$Infector.ID
  ind.infectee <- subdata3$Infectee.ID
  ind.subtract <- c()
  for(b in 1:length(ind.infector)){
    ind.subtract <- c(ind.subtract, which(ind.infectee == ind.infector[b]))
  }
  if(length(ind.subtract) >= 1){
    subdata3 <- subdata3[-ind.subtract,]
    }
  # now this is the final sample size
  samplesize[i] <- nrow(subdata2) + nrow(subdata3)
  
  dftmp <- data.frame(left = c(with(subdata2, as.numeric(Onset_Infector - expo.late.I)),
                               with(subdata3, as.numeric(Onset_Infectee - expo.late.S))),
                      right = c(with(subdata2, as.numeric(Onset_Infector - expo.early.I)),
                                with(subdata3, as.numeric(Onset_Infectee - expo.early.S))))
  fit.g <- fitdistcens(dftmp, "gamma")
  fit.wb <- fitdistcens(dftmp, "weibull")
  fit.ln <- fitdistcens(dftmp, "lnorm")
  
  fit.g.boot <- bootdistcens(fit.g, niter = 1000)
  fit.wb.boot <- bootdistcens(fit.wb, niter = 1000)
  fit.ln.boot <- bootdistcens(fit.ln, niter = 1000)
  
  par1[[i]] <- c(fit.g$estimate[1], fit.wb$estimate[1], fit.ln$estimate[1])
  par2[[i]] <- c(fit.g$estimate[2], fit.wb$estimate[2], fit.ln$estimate[2])
  par1.boot[[i]] <- list(fit.g.boot$estim[, 1], fit.wb.boot$estim[, 1], fit.ln.boot$estim[, 1])
  par2.boot[[i]] <- list(fit.g.boot$estim[, 2], fit.wb.boot$estim[, 2], fit.ln.boot$estim[, 2])
  AIC[[i]] <- c(fit.g$aic, fit.wb$aic,fit.ln$aic)
  
  setTxtProgressBar(progressbar, i)
}

IPaggregate.reslist <- list(
  par1 = par1,
  par2 = par2,
  par1boot = par1.boot,
  par2boot = par2.boot,
  samplesize = samplesize,
  timewindow = paste0(format(startseq, format =  "%b %d"), " - ",
                      format(endseq, format =  "%b %d")),
  AIC = AIC,
  dists = c("Gamma", "Weibull", "Log-Normal")
)

save(IPaggregate.reslist, file = "IPaggregate.reslist.rda")

load("IPaggregate.reslist.rda")

IPaggregate.reslist$par1
IPaggregate.reslist$par2
IPaggregate.reslist$AIC
IPaggregate.reslist$timewindow

AICs <- unlist(IPaggregate.reslist$AIC)

df.AIC.IPaggregate <- data.frame(
  timewindow = rep(IPaggregate.reslist$timewindow, 3),
  AIC = c(AICs[seq(1, 49, 3)], AICs[seq(2, 50, 3)], AICs[seq(3, 51, 3)]),
  dist = c(rep("Gamma", 17), rep("Weibull", 17), rep("Log-Normal", 17)),
  samplesize = rep(IPaggregate.reslist$samplesize, 3)
)

write.csv(df.AIC.IPaggregate, "df.AIC.IPaggregate.csv")

# best fitted by gamma

df.IPaggregate.gamma <- data.frame(
  shape = unlist(IPaggregate.reslist$par1)[seq(1, 49, 3)],
  rate = unlist(IPaggregate.reslist$par2)[seq(1, 49, 3)]
)

write.csv(dfIPaggregate.gamma, "pars_IPaggregate_gamma.csv")

IPaggregate.mu.LB <- IPaggregate.mu.UB <- IPaggregate.sd.LB <- IPaggregate.sd.UB <-numeric(17)

for(i in 1:17){
  par1boot <- unlist(IPaggregate.reslist$par1boot[[i]][1])
  par2boot <- unlist(IPaggregate.reslist$par2boot[[i]][1])
  IPaggregate.mu.LB[i] <- quantile(par1boot/par2boot, 0.025)
  IPaggregate.mu.UB[i] <- quantile(par1boot/par2boot, 0.975)
  IPaggregate.sd.LB[i] <- quantile(sqrt(par1boot/par2boot^2), 0.025)
  IPaggregate.sd.UB[i] <- quantile(sqrt(par1boot/par2boot^2), 0.975)
}

IPaggregate.gamma <- data.frame(
  mu = df.IPaggregate.gamma$shape/df.IPaggregate.gamma$rate,
  muLB = IPaggregate.mu.LB,
  muUB = IPaggregate.mu.UB,
  sd = sqrt(df.IPaggregate.gamma$shape/df.IPaggregate.gamma$rate^2),
  sdLB = IPaggregate.sd.LB,
  sdUB = IPaggregate.sd.UB
)

write.csv(IPaggregate.gamma, "IPaggregate_fitgamma.csv")


#################################### estimate GT
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

startseq <- c(as.Date("2020-01-01"), seq.Date(from = as.Date("2020-01-15"), to = as.Date("2020-01-30"), by = "1 day"))
endseq <- c(as.Date("2020-01-20"), seq.Date(from = as.Date("2020-01-21"), to = as.Date("2020-02-04"), by = "1 day"),
            as.Date("2020-02-29"))



################## update sensitivity results

IPdistInfector <- IPdistInfectee <- df.IPaggregate.gamma

as.numeric(IPdistInfector[1,])

ncores <- detectCores()

myCluster <- makeCluster(ncores - 1, # number of cores to use
                         type = "PSOCK") # type of cluster

registerDoParallel(myCluster)
reslist <- vector("list", 17)
nwindow <- 17

progressbar <- txtProgressBar(min = 0, max = nwindow, style = 3)

for(j in 1:nwindow){
  nsim = 1000 # 
  tmppair <- subset(allpairs, Onset_Infector >= startseq[j] & Onset_Infector <= endseq[j])  
  reslist[[j]] <- foreach(i = 1:nsim, .combine = 'append') %dopar% {
    simsamples <- sampleGT(tmppair, allpairs, IPdistInfector = as.numeric(IPdistInfector[j,]), InfectorIPdistfrom = "gamma", 
                          IPdistInfectee = as.numeric(IPdistInfectee[j,]), InfecteeIPdistfrom = "gamma", iseed = i)
    GTsamples <- simsamples[[1]]
  }
  setTxtProgressBar(progressbar, j)
}
stopCluster(myCluster)
# ABOUT HALF AN HOUR TO FINISH THIS

# View(reslist[[17]])

save(reslist, file = "GTsamples_IP_nodistinguish.rda")

load("GTsamples_IP_nodistinguish.rda")

source("GT_estimation_function.R")

GTlist <- vector("list", 17)
for(i in 1:17){
  GTsamples <- reslist[[i]]
  GTmat <- matrix(unlist(GTsamples), nrow = 1000, byrow = T)
  GTmatlist <- split(GTmat, rep(1:ncol(GTmat), each = nrow(GTmat)))
  GTlist[[i]] <- GTmatlist
  rm(GTsamples)
  rm(GTmat)
  rm(GTmatlist)
}

lapply(GTlist, function(x) sum(is.infinite(unlist(x))))
# no infinite samples

pars.est <- vector("list", 17)
negsumllk <- vector("list", 17)

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
#  pars.est[[i]] <- exp(test$par)
  negsumllk[[i]] <- test$value
  
  setTxtProgressBar(progressbar, i)
}

indseq <- seq(1, 34)
indpar1 <- indseq[indseq %% 2 == 1]
indpar2 <- indseq[indseq %% 2 == 0]


df.estGT <- data.frame(
  mu = unlist(pars.est)[indpar1],
  sd = unlist(pars.est)[indpar2],
  negsumllk = unlist(negsumllk),
  timewindow = paste0(format(startseq, format =  "%b %d"), " - ",
                      format(endseq, format =  "%b %d"))
)

df.estGT.g <- df.estGT
df.estGT.wb <- df.estGT
df.estGT.ln <- df.estGT


write.csv(df.estGT.g, "GTfor_IPnostrat_estgamma.csv")
write.csv(df.estGT.wb, "GTfor_IPnostrat_estweibull.csv")
write.csv(df.estGT.ln, "GTfor_IPnostrat_estlnorm.csv")

df.estGT.ln <- read.csv("GTfor_IPnostrat_estlnorm.csv")

# best fit by lnorm
# get bootCI

par1UB.boot <- par2UB.boot <- numeric(17)
par1LB.boot <- par2LB.boot <- numeric(17)
parsboot <- vector("list", 17)

ncores <- detectCores()
myCluster <- makeCluster(ncores - 1, # number of cores to use
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
      fn  = llh.fx.ln,
      # fn = llh.fx.wb,
      # fn = llh.fx,
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

save(parsboot, file = "forwardGT_sensitivity_bootpars_lnorm_MAY06.rda")

df.estGT.ln$muLB.boot <- par1LB.boot
df.estGT.ln$muUB.boot <- par1UB.boot
df.estGT.ln$sdLB.boot <- par2LB.boot
df.estGT.ln$sdUB.boot <- par2UB.boot

write.csv(df.estGT.ln, "GTfor_IPnostrat_estlnorm.csv")

df.estGT.ln <- read.csv("GTfor_IPnostrat_estlnorm.csv")

################## compare GT based on main analysis, fixed IP, semi fixed IP

GTmain <- read.csv("forwardGT_lnorm_May06.csv")
GTmain[,4]

GTcomb <- rbind(GTmain[, -c(4, 11:15)], df.estGT.ln[,-1])

GTcomb$type <- c(rep("main", 17), rep("IP without stratification", 17))
GTcomb$X <- rep(seq(1, 17), 2)


p.GT.comb.m <- ggplot(GTcomb, aes(x = X, y = mu, ymin = muLB.boot, ymax = muUB.boot, group = type, color = type)) +
  geom_point(position = position_dodge(0.7), size = 1.0) + 
  geom_errorbar(position = position_dodge(0.7), width = 0, size = 0.5) +
  geom_text(mapping = aes(y = muUB.boot+0.30, label = sprintf("%0.2f", round(mu, digits = 2))), size = 3.25, position = position_dodge(1.0)) +
  labs(x="Mid point of time window", y="Mean generation time (days)") + 
  scale_x_continuous(limits = c(0, 18), breaks = seq(1, 17, 1), labels =
                       c("Jan 10",
                         "Jan 18", "Jan 19", "Jan 20", "Jan 21", "Jan 22",
                         "Jan 23", "Jan 24", "Jan 25", "Jan 26", "Jan 27", "Jan 28",
                         "Jan 29", "Jan 30", "Jan 31", "Feb 1",
                         "Feb 14")) +
  scale_x_cut(breaks = c(1.5, 16.5), space = 0.5, which = c(1, 3), scales = c(0.125, 0.125)) + 
  theme(
    axis.ticks.x.top = element_blank(),
    axis.text.x.top = element_blank(),
    axis.line.x.top = element_blank(),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    axis.text.x = element_text(),
    axis.text = element_text(size=8.5),
    axis.line = element_line(color="black", size=0.75),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  ) + scale_y_continuous(breaks = seq(3, 10, 1), limits = c(2, 10)) 
#  ggtitle("Comparison of estimated mean of forward generation time")

p.GT.comb.m  
  
  
  
p.GT.comb.sd <- ggplot(GTcomb, aes(x = X, y = sd, ymin = sdLB.boot, ymax = sdUB.boot, group = type, color = type)) +
  geom_point(position = position_dodge(0.7), size = 1.0) +  
  geom_errorbar(position = position_dodge(0.7), width = 0, size = 0.5) +
  geom_text(mapping = aes(y = sdUB.boot+0.30, label = sprintf("%0.2f", round(sd, digits = 2))), size = 3.25, position = position_dodge(1.0)) +
  labs(x="Mid point of time window", y="Standard deviation of generation time (days)") +
  scale_x_continuous(limits = c(0, 18), breaks = seq(1, 17, 1), labels =
                       c("Jan 10",
                         "Jan 18", "Jan 19", "Jan 20", "Jan 21", "Jan 22",
                         "Jan 23", "Jan 24", "Jan 25", "Jan 26", "Jan 27", "Jan 28",
                         "Jan 29", "Jan 30", "Jan 31", "Feb 1",
                         "Feb 14")) +
  scale_x_cut(breaks = c(1.5, 16.5), space = 0.5, which = c(1, 3), scales = c(0.125, 0.125)) + 
  theme(
    axis.ticks.x.top = element_blank(),
    axis.text.x.top = element_blank(),
    axis.line.x.top = element_blank(),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    axis.text.x = element_text(),
    axis.text = element_text(size=8.5),
    axis.line = element_line(color="black", size=0.75),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  ) + scale_y_continuous(breaks = seq(2, 6, 1), limits = c(1, 6.5)) 
p.GT.comb.sd   

pGT.sensitivity <- plot_grid(print(p.GT.comb.m), print(p.GT.comb.sd), nrow = 2,
                             labels = c("a", "b"))


ggsave("Fig_sensitivity.pdf", width = 12.5, height = 8, units = "in", plot = pGT.sensitivity)






