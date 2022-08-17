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


library(foreach)
library(doParallel)

####### estimate GT

source("GT_sampling_corr_functions.R")

startseq <- c(as.Date("2020-01-01"), seq.Date(from = as.Date("2020-01-15"), to = as.Date("2020-01-30"), by = "1 day"))
endseq <- c(as.Date("2020-01-20"), seq.Date(from = as.Date("2020-01-21"), to = as.Date("2020-02-04"), by = "1 day"),
            as.Date("2020-02-29"))

ncores <- detectCores()
myCluster <- makeCluster(ncores - 2, # number of cores to use
                         type = "PSOCK") # type of cluster

registerDoParallel(myCluster)
reslist <- vector("list", 17)
nwindow <- 17

load("IPlnormpars_weightinfector.rda")
load("IPlnormpars_infectee.rda")

df.IPlnormpars.weightinfector
df.IPlnormpars.infectee

progressbar <- txtProgressBar(min = 0, max = nwindow, style = 3)

for(j in 1:nwindow){
  nsim = 1000 
  IPdistInfector <- c(df.IPlnormpars.weightinfector[j,1], df.IPlnormpars.weightinfector[j,2])
  IPdistInfectee <- c(df.IPlnormpars.infectee[j,1], df.IPlnormpars.infectee[j,2])
  tmppair <- subset(allpairs, Onset_Infector >= startseq[j] & Onset_Infector <= endseq[j])  
  reslist[[j]] <- foreach(i = 1:nsim, .combine = 'append') %dopar% {
    simsamples <- sampleGT(tmppair, allpairs, IPdistInfector = IPdistInfector, InfectorIPdistfrom = "lnorm", 
                           IPdistInfectee = IPdistInfectee, InfecteeIPdistfrom = "lnorm", iseed = i)
    list(
      GTsamples = simsamples[[1]],
      logGTsamples = simsamples[[2]],
      logIPsamples = simsamples[[3]]
    )
  }
  setTxtProgressBar(progressbar, j)
}
stopCluster(myCluster)

save(reslist, file = "forwardGT1000sim_corr_JUNE20.RData")

load("forwardGT1000sim_corr_JUNE20.RData")

samplesize <- numeric(17)
for(i in 1:17){
  samplesize[i] <- length(which(allpairs$Onset_Infector >= startseq[i] & allpairs$Onset_Infector <= endseq[i]))
}

refseq <- seq(1, 3000)
indGTlog <- refseq[refseq %% 3 == 2]
indIPlog <- refseq[refseq %% 3 == 0]

datalist <- vector("list", 17)
for(i in 1:17){
  GTlogsamples <- reslist[[i]][indGTlog]
  IPlogsamples <- reslist[[i]][indIPlog]
  GTlogmat <- matrix(unlist(GTlogsamples), ncol = samplesize[i], byrow = T)
  IPlogmat <- matrix(unlist(IPlogsamples), ncol = samplesize[i], byrow = T)
  datalist[[i]] <- vector("list", samplesize[i])
  for(j in 1:samplesize[i]){
    datalist[[i]][[j]] <- cbind(GTlogmat[,j], IPlogmat[,j])
  }
  
  rm(GTlogsamples)
  rm(IPlogsamples)
  rm(GTlogmat)
  rm(IPlogmat)
}

source("GT_estimation_corr_function.R")

par.est <- vector("list", 17)
# llk.est <- numeric(17)
nwindow <- 17

progressbar <- txtProgressBar(min = 0, max = nwindow, style = 3)

for(i in 1:17){
  
  
  test <- optim(
    par = c(0.5, 0.5, 0),
    fn = llh.fx.ln,
    method = "L-BFGS-B",
    hessian = T,
    pars.fix = c(df.IPlnormpars.weightinfector[i,1],
                 df.IPlnormpars.weightinfector[i,2]),
    data.mc = datalist[[i]]
  )
  
  par.est[[i]] <- test$par
#  llk.est[i] <- test$value
  print(i)
  
  setTxtProgressBar(progressbar, i)
}

par.est
# llk.est


indmu <- seq(1, 51, 3)
indsd <- seq(2, 51, 3)
ind.rho <- seq(3, 51, 3)

GT.est_corr.df <- data.frame(
  mean = unlist(par.est)[indmu],
  sd = unlist(par.est)[indsd],
  rho = tanh(unlist(par.est)[ind.rho]),
  llk = llk.est
)

save(GT.est_corr.df, file = "estGTforward_corr_6July.rda")

load("estGTforward_corr_6July.rda")

# if we set rho to be fixed

par.est <- vector("list", 17)
# llk.est <- numeric(17)
nwindow <- 17

progressbar <- txtProgressBar(min = 0, max = nwindow, style = 3)

for(i in 1:17){
  
  
  test <- optim(
    par = c(0.5, 0.5),
    fn = llh.fx.ln.adj,
    method = "L-BFGS-B",
    hessian = T,
    pars.fix = c(0.75, # 0.5, 0.25
      df.IPlnormpars.weightinfector[i,1],
      df.IPlnormpars.weightinfector[i,2]),
    data.mc = datalist[[i]]
  )
  
  par.est[[i]] <- test$par
  #  llk.est[i] <- test$value
  print(i)
  
  setTxtProgressBar(progressbar, i)
}

par.est
# llk.est


indmu <- seq(1, 34, 2)
indsd <- seq(2, 34, 2)

GT.est_corr.df_rho0.25 <- data.frame(
  mean = unlist(par.est)[indmu],
  sd = unlist(par.est)[indsd]
)

GT.est_corr.df_rho0.5 <- data.frame(
  mean = unlist(par.est)[indmu],
  sd = unlist(par.est)[indsd]
)

GT.est_corr.df_rho0.75 <- data.frame(
  mean = unlist(par.est)[indmu],
  sd = unlist(par.est)[indsd]
)

GT.est_corr_rhofix.df <- rbind(GT.est_corr.df_rho0.25, 
                               GT.est_corr.df_rho0.5, 
                               GT.est_corr.df_rho0.75)

save(GT.est_corr.df, file = "estGTforward_corr_6July.rda")

load("estGTforward_corr_6July.rda")

save(GT.est_corr_rhofix.df, file = "estGTforward_corr_rhofix_7July.rda")

load("estGTforward_corr_rhofix_7July.rda")

# BOOTSTRAP CI
library(foreach)
library(doParallel)


meanUB.boot <- meanLB.boot <- numeric(17)
sdUB.boot <- sdLB.boot <- numeric(17)
rhoLB.boot <- rhoUB.boot <- numeric(17)

parsboot <- vector("list", 17)

ncores <- detectCores()

myCluster <- makeCluster(ncores - 2, # number of cores to use
                         type = "PSOCK") # type of cluster

registerDoParallel(myCluster)

start.t <- Sys.time()

progressbar <- txtProgressBar(min = 0, max = 17, style = 3)

for(i in 1:17){
  
  data.mc <- datalist[[i]]
  
TMP <- foreach(
  b = 1:1000,
  .packages = "mvtnorm"
) %dopar% {
  
  set.seed(b)
  
  ind.b = sample(1:length(data.mc), length(data.mc), replace = T)
  data.mc.b = data.mc[ind.b]
  
  test <- optim(
    par = c(0.5, 0.5, 0),
    fn = llh.fx.ln,
    method = "L-BFGS-B",
    hessian = T,
    pars.fix = c(df.IPlnormpars.weightinfector[i,1],
                 df.IPlnormpars.weightinfector[i,2]),
    data.mc = data.mc.b
  )
  
  return(test)
}
  
  tmpboot <- sapply(TMP, function(b) { b$par})
  
  parsboot[[i]] <- tmpboot
  indmean <- seq(1, 3000, 3)
  indsd <- seq(2, 3000, 3)
  indrho <- seq(3, 3000, 3)
  
  meanboot <- unlist(tmpboot)[indmean]
  sdboot <- unlist(tmpboot)[indsd]
  rhoboot <- tanh(unlist(tmpboot)[indrho])
  
  meanLB.boot[i] <- quantile(meanboot, 0.025)
  meanUB.boot[i] <- quantile(meanboot, 0.975)
  
  sdLB.boot[i] <- quantile(sdboot, 0.025)
  sdUB.boot[i] <- quantile(sdboot, 0.975)
  
  rhoLB.boot[i] <- quantile(rhoboot, 0.025)
  rhoUB.boot[i] <- quantile(rhoboot, 0.975)
  
  setTxtProgressBar(progressbar, i)
}


stopCluster(myCluster)


end.t = Sys.time()
end.t - start.t 

save(parsboot, file = "forwardGT1000_bootpars_corr_JUNE29.RData")


########### boot CI for fix rho
setrho <- c(0.25, 0.5, 0.75)

muUBlist <- muLBlist <- sdUBlist <- sdLBlist <- vector("list", 3)
parsbootlist <- vector("list", 3)

start.t <- Sys.time()

for(m in 1:3){

meanUB.boot <- meanLB.boot <- numeric(17)
sdUB.boot <- sdLB.boot <- numeric(17)
rhoLB.boot <- rhoUB.boot <- numeric(17)

parsboot <- vector("list", 17)

ncores <- detectCores()

myCluster <- makeCluster(ncores - 2, # number of cores to use
                         type = "PSOCK") # type of cluster

registerDoParallel(myCluster)

progressbar <- txtProgressBar(min = 0, max = 17, style = 3)

for(i in 1:17){
  
  data.mc <- datalist[[i]]
  
  TMP <- foreach(
    b = 1:1000,
    .packages = "mvtnorm"
  ) %dopar% {
    
    set.seed(b)
    
    ind.b = sample(1:length(data.mc), length(data.mc), replace = T)
    data.mc.b = data.mc[ind.b]
    
    test <- optim(
      par = c(0.5, 0.5),
      fn = llh.fx.ln.adj,
      method = "L-BFGS-B",
      hessian = T,
      pars.fix = c(setrho[m],
                   df.IPlnormpars.weightinfector[i,1],
                   df.IPlnormpars.weightinfector[i,2]),
      data.mc = data.mc.b
    )
    return(test)
  }
  
  tmpboot <- sapply(TMP, function(b) { b$par})
  
  parsboot[[i]] <- tmpboot
  indmean <- seq(1, 2000, 2)
  indsd <- seq(2, 2000, 2)
  
  meanboot <- unlist(tmpboot)[indmean]
  sdboot <- unlist(tmpboot)[indsd]
  
  meanLB.boot[i] <- quantile(meanboot, 0.025)
  meanUB.boot[i] <- quantile(meanboot, 0.975)
  
  sdLB.boot[i] <- quantile(sdboot, 0.025)
  sdUB.boot[i] <- quantile(sdboot, 0.975)
  
  setTxtProgressBar(progressbar, i)
}

stopCluster(myCluster)

muUBlist[[m]] <- meanUB.boot
muLBlist[[m]] <- meanLB.boot

sdUBlist[[m]] <- sdUB.boot
sdLBlist[[m]] <- sdLB.boot

parsbootlist[[m]] <- parsboot

print(m)

}


end.t = Sys.time()
end.t - start.t 



GTfor_corr_bootCI <- data.frame(
  muUB = meanUB.boot,
  muLB = meanLB.boot,
  sdUB = sdUB.boot,
  sdLB = sdLB.boot,
  rhoUB = rhoUB.boot,
  rhoLB = rhoLB.boot
)


GTfor_corr_bootCI.rhofix <- data.frame(
  muUB = unlist(muUBlist),
  muLB = unlist(muLBlist),
  sdUB = unlist(sdUBlist),
  sdLB = unlist(sdLBlist),
  rhofix = c(rep(0.25, 17), rep(0.5, 17), rep(0.75, 17))
)


GT.est_corr_CI.df <- cbind(GT.est_corr.df, GTfor_corr_bootCI)

GT.est_corr_CI_rhofix.df <- cbind(GT.est_corr_rhofix.df, GTfor_corr_bootCI.rhofix)

GT.est_corr_CI.df$timewindow <- paste(startseq, endseq)

GT.est_corr_CI_rhofix.df$timewindow <- rep(paste(startseq, endseq), 3)

write.csv(GT.est_corr_CI.df, file = "estGT_forward_corr.csv")

write.csv(GT.est_corr_CI_rhofix.df, file = "estGT_forward_corr_rhofix.csv")

GT.est_corr_CI.df <- read.csv("estGT_forward_corr.csv")

GT.est_corr_CI_rhofix.df <- read.csv("estGT_forward_corr_rhofix.csv")
#################
library(ggplot2)
library(ggbreak)

df.estGT.ln <- read.csv("forwardGT_lnorm_May06.csv")

df.estGT.ln$X <- seq(1, 17, 1)

GTcompare <- data.frame(
  X = rep(seq(1, 17, 1), 5),
  type = c(rep("independent", 17), rep("correlated, rho to be estimated", 17), rep("correlated, rho set at 0.25", 17),
           rep("correlated, rho set at 0.5", 17), rep("correlated, rho set at 0.75", 17)),
  mu = c(df.estGT.ln$mu, GT.est_corr_CI.df$mean, GT.est_corr_CI_rhofix.df$mean),
  muLB.boot = c(df.estGT.ln$muLB.boot, GT.est_corr_CI.df$muLB, GT.est_corr_CI_rhofix.df$muLB),
  muUB.boot = c(df.estGT.ln$muUB.boot, GT.est_corr_CI.df$muUB, GT.est_corr_CI_rhofix.df$muUB),
  sd = c(df.estGT.ln$sd, GT.est_corr_CI.df$sd, GT.est_corr_CI_rhofix.df$sd),
  sdLB.boot = c(df.estGT.ln$sdLB.boot, GT.est_corr_CI.df$sdLB, GT.est_corr_CI_rhofix.df$sdLB),
  sdUB.boot = c(df.estGT.ln$sdUB.boot, GT.est_corr_CI.df$sdUB, GT.est_corr_CI_rhofix.df$sdUB)
)

p.GT.forward.main.bootCI <- ggplot(GTcompare, aes(x = X, y = mu, ymin = muLB.boot, ymax = muUB.boot, group = type, color = type)) +
  geom_point(size = 2.0, position = position_dodge(0.8), aes(shape = type)) + 
  geom_errorbar(width = 0, size = 0.5, position = position_dodge(0.8)) +
#  geom_text(mapping = aes(y = muUB.boot+0.30, label = sprintf("%0.2f", round(mu, digits = 2))), position = position_dodge(0.8), size = 3.0)+
  labs(x="Mid point of time window", y="Mean generation time (days)") +
  scale_x_continuous(limits = c(0, 18), breaks = seq(1, 17, 1), labels =
                       c("Jan 10",
                         "Jan 18", "Jan 19", "Jan 20", "Jan 21", "Jan 22",
                         "Jan 23", "Jan 24", "Jan 25", "Jan 26", "Jan 27", "Jan 28",
                         "Jan 29", "Jan 30", "Jan 31", "Feb 1",
                         "Feb 14")) +
  scale_x_cut(breaks = c(1.5, 16.5), space = 0.5, which = c(1, 3), scales = c(0.125, 0.125)) + 
  #  scale_x_continuous(breaks = seq(1, 17, 1), labels =
  #                       c("Jan 10", "Jan 18", "Jan 19", "Jan 20", "Jan 21", "Jan 22",
  #                         "Jan 23", "Jan 24", "Jan 25", "Jan 26", "Jan 27", "Jan 28",
  #                         "Jan 29", "Jan 30", "Jan 31", "Feb 1", "Feb 14")) +
  geom_segment(aes(x = 7, y = 8.5, xend = 7, yend = 7.5), 
               arrow = arrow(length = unit(0.05, "inches")), lwd=0.5, color = "red") +
  geom_segment(aes(x = 7, y = 8.5, xend = 7.5, yend = 8.5), 
               lwd=0.5, color = "red") +
  annotate("text", x = 9 + 0.5, y = 8.5, label = "Lockdown in Wuhan", size = 3.5) + 
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.text.x.top = element_blank(),
    axis.line.x.top = element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    axis.text.x = element_text(),
    axis.text = element_text(size=8.5),
    axis.line = element_line(size=0.75),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  ) + 
  coord_cartesian(ylim=c(2.5, 9.5)) + scale_y_continuous(breaks = seq(3, 9, 1)) 
#  ggtitle("Estimated mean of forward generation time (Mainland China, 2020)")

p.GT.forward.sd.bootCI <- ggplot(GTcompare, aes(x = X, y = sd, ymin = sdLB.boot, ymax = sdUB.boot, group = type, color = type)) +
  geom_point(size = 2.0, position = position_dodge(0.8), aes(shape = type)) + 
  scale_x_continuous(limits = c(0, 18), breaks = seq(1, 17, 1), labels =
                       c("Jan 10",
                         "Jan 18", "Jan 19", "Jan 20", "Jan 21", "Jan 22",
                         "Jan 23", "Jan 24", "Jan 25", "Jan 26", "Jan 27", "Jan 28",
                         "Jan 29", "Jan 30", "Jan 31", "Feb 1",
                         "Feb 14")) + 
  geom_errorbar(width = 0, size = 0.5, position = position_dodge(0.8)) +
#  geom_text(mapping = aes(y = sdUB.boot+0.30, label = sprintf("%0.2f", round(sd, digits = 2))), size = 3.5, position = position_dodge(0.7)) +
  #  geom_text(mapping = aes(y = 4, label = n_sample), size = 2.5) +
  labs(x="Mid point of time window", y="Standard deviation of generation time (days)") +
  scale_x_cut(breaks = c(1.5, 16.5), space = 0.5, which = c(1, 3), scales = c(0.125, 0.125)) +
  geom_segment(aes(x = 7, y = 6.5, xend = 7, yend = 5), 
               arrow = arrow(length = unit(0.05, "inches")), lwd=0.5, color = "red") +
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
    axis.line = element_line(size=0.75),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  ) + scale_y_continuous(breaks = seq(2, 7, 1), limits = c(1, 8)) 
#  ggtitle("Estimated standard deviation of forward generation time")

# library(aplot)

# pGT.main.bootCI <- plot_list(gglist = list(p.GT.forward.main.bootCI, p.GT.forward.sd.bootCI),
#           ncol = 1, nrow = 2, labels = c("A", "B"), tag_size = 14)

library(cowplot)

legend_common <- get_legend(p.GT.forward.sd.bootCI)

pGT.main.bootCI <- plot_grid(print(p.GT.forward.main.bootCI), print(p.GT.forward.sd.bootCI), nrow = 2,
                             labels = c("a", "b"))



ggsave("forwardGT_compare.pdf", width = 12.5, height = 8, units = "in")

prho <- ggplot(GT.est_corr_CI.df, aes(x = seq(1, 17, 1), y = rho, ymin = rhoLB, ymax = rhoUB)) +
  geom_point(size = 1.0) + 
  geom_errorbar(width = 0, size = 0.5) +
  geom_text(mapping = aes(y = rhoUB+0.030, label = sprintf("%0.2f", round(rho, digits = 2))))+
  labs(x="Mid point of time window", y="Correlation coefficient (estimated)") +
  scale_x_continuous(limits = c(0, 18), breaks = seq(1, 17, 1), labels =
                       c("Jan 10",
                         "Jan 18", "Jan 19", "Jan 20", "Jan 21", "Jan 22",
                         "Jan 23", "Jan 24", "Jan 25", "Jan 26", "Jan 27", "Jan 28",
                         "Jan 29", "Jan 30", "Jan 31", "Feb 1",
                         "Feb 14")) +
  scale_x_cut(breaks = c(1.5, 16.5), space = 0.5, which = c(1, 3), scales = c(0.125, 0.125)) + 
  #  scale_x_continuous(breaks = seq(1, 17, 1), labels =
  #                       c("Jan 10", "Jan 18", "Jan 19", "Jan 20", "Jan 21", "Jan 22",
  #                         "Jan 23", "Jan 24", "Jan 25", "Jan 26", "Jan 27", "Jan 28",
  #                         "Jan 29", "Jan 30", "Jan 31", "Feb 1", "Feb 14")) +
  geom_segment(aes(x = 7, y = 0.85, xend = 7, yend = 0.7), 
               arrow = arrow(length = unit(0.05, "inches")), lwd=0.5, color = "red") +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.text.x.top = element_blank(),
    axis.line.x.top = element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    axis.text.x = element_text(),
    axis.text = element_text(size=8.5),
    axis.line = element_line(size=0.75),
    panel.grid = element_blank(),
    legend.position = "top"
  ) + 
  coord_cartesian(ylim=c(-0.05, 0.85)) + scale_y_continuous(breaks = seq(0, 0.85, 0.05)) 

ggsave("time_moving_rho.pdf", width = 9, height = 6, units = "in")





















