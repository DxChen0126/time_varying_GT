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


# prepare IP estimates (best fit is weibull)

# 

dfIP.backward <- read.csv("dfIPinfector.backward_fitdistWB.csv")

dfIP.forward <- read.csv("dfIPinfectee.forward_fitdistWB.csv")


library(foreach)
library(doParallel)
library(mixdist)
####### estimate GT

source("GT_sampling_functions.R")

startseq <- c(as.Date("2020-01-01"), seq.Date(from = as.Date("2020-01-15"), to = as.Date("2020-01-30"), by = "1 day"))
endseq <- c(as.Date("2020-01-20"), seq.Date(from = as.Date("2020-01-21"), to = as.Date("2020-02-04"), by = "1 day"),
            as.Date("2020-02-29"))

ncores <- detectCores()
myCluster <- makeCluster(ncores - 2, # number of cores to use
                         type = "PSOCK") # type of cluster

registerDoParallel(myCluster)
reslist <- vector("list", 17)
nwindow <- 17

progressbar <- txtProgressBar(min = 0, max = nwindow, style = 3)

for(j in 1:nwindow){
  nsim = 1000 
  IPdistInfector <- c(dfIP.backward$alpha[j], dfIP.backward$beta[j])
  IPdistInfectee <- c(dfIP.forward$alpha[j], dfIP.forward$beta[j])
  tmppair <- subset(allpairs, Onset_Infector >= startseq[j] & Onset_Infector <= endseq[j])  
  reslist[[j]] <- foreach(i = 1:nsim, .combine = 'append') %dopar% {
    simsamples <- sampleGT(tmppair, allpairs, IPdistInfector = IPdistInfector, InfectorIPdistfrom = "weibull", 
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

save(reslist, file = "forwardGT1000sim_samples_May06.rda")

samplesize <- numeric(17)
for(i in 1:17){
  samplesize[i] <- length(which(allpairs$Onset_Infector >= startseq[i] & allpairs$Onset_Infector <= endseq[i]))
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

source("GT_estimation_function.R")

pars.est <- vector("list", 17)
negsumllk <- vector("list", 17)

lapply(GTlist, function(x) sum(is.infinite(unlist(x))))

progressbar <- txtProgressBar(min = 0, max = nwindow, style = 3)

for(i in 1:17){
  
  test <- optim(
    par = c(0, 0),
    # fn  = llh.fx,
    # fn = llh.fx.wb,
     fn = llh.fx.ln,
    method = "Nelder-Mead",
    hessian = T,
    y.mc = GTlist[[i]]
  )
  
  # for wb only
#  tmppar <- exp(test$par)
#  pars.est[[i]] <- as.numeric(weibullparinv(shape = tmppar[1], scale = tmppar[2])[1:2])
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


# BOOTSTRAP CI
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
    # fn  = llh.fx.ln,
    fn = llh.fx.wb,
   # fn = llh.fx,
    method = "Nelder-Mead",
    hessian = T,
    y.mc = y.mc.b
  )
  return(test)
}
  
 # tmpboot <- sapply(TMP, function(b) { exp(b$par) })
  
  # for wb only
   tmpboot <- sapply(TMP, function(b) {as.numeric(weibullparinv(shape = exp(b$par)[1], scale = exp(b$par)[2])[1:2])})
  
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

# save(parsboot, file = "forwardGT1000_bootpars_gamma_MAY06.rda")
save(parsboot, file = "forwardGT1000_bootpars_wb_MAY06.rda")
# save(parsboot, file = "forwardGT1000_bootpars_lnorm_MAY06.rda")

df.estGT.g$muLB.boot <- par1LB.boot
df.estGT.g$muUB.boot <- par1UB.boot
df.estGT.g$sdLB.boot <- par2LB.boot
df.estGT.g$sdUB.boot <- par2UB.boot

df.estGT.wb$muLB.boot <- par1LB.boot
df.estGT.wb$muUB.boot <- par1UB.boot
df.estGT.wb$sdLB.boot <- par2LB.boot
df.estGT.wb$sdUB.boot <- par2UB.boot


df.estGT.ln$muLB.boot <- par1LB.boot
df.estGT.ln$muUB.boot <- par1UB.boot
df.estGT.ln$sdLB.boot <- par2LB.boot
df.estGT.ln$sdUB.boot <- par2UB.boot

write.csv(df.estGT.g, "forwardGT_gamma_May06.csv")

write.csv(df.estGT.wb, "forwardGT_weibull_May06.csv")

write.csv(df.estGT.ln, "forwardGT_lnorm_May06.csv")

df.estGT.ln <- read.csv("forwardGT_lnorm_May06.csv")


library(ggplot2)
library(ggbreak)

df.estGT.ln$X <- seq(1, 17, 1)

p.GT.forward.main.bootCI <- ggplot(df.estGT.ln, aes(x = X, y = mu, ymin = muLB.boot, ymax = muUB.boot)) +
  geom_point(size = 1.0) + 
  geom_errorbar(width = 0, size = 0.5) +
  geom_text(mapping = aes(y = muUB.boot+0.30, label = sprintf("%0.2f", round(mu, digits = 2))), size = 3.5)+
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
  geom_segment(aes(x = 7, y = 8.5, xend = 7, yend = 7), 
              arrow = arrow(length = unit(0.05, "inches")), lwd=0.5, color = "red") +
  geom_segment(aes(x = 7, y = 8.5, xend = 7.5, yend = 8.5), lwd=0.5, color = "red") +
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
    panel.grid = element_blank()
  ) + 
  coord_cartesian(ylim=c(2.5, 8.5)) + scale_y_continuous(breaks = seq(4, 8, 1)) 
#  ggtitle("Estimated mean of forward generation time (Mainland China, 2020)")

p.GT.forward.sd.bootCI <- ggplot(df.estGT.ln, aes(x = X, y = sd, ymin = sdLB.boot, ymax = sdUB.boot)) +
  geom_point(size = 1.0) + 
  scale_x_continuous(limits = c(0, 18), breaks = seq(1, 17, 1), labels =
                       c("Jan 10",
                         "Jan 18", "Jan 19", "Jan 20", "Jan 21", "Jan 22",
                         "Jan 23", "Jan 24", "Jan 25", "Jan 26", "Jan 27", "Jan 28",
                         "Jan 29", "Jan 30", "Jan 31", "Feb 1",
                         "Feb 14")) + 
  geom_errorbar(width = 0, size = 0.5) +
  geom_text(mapping = aes(y = sdUB.boot+0.30, label = sprintf("%0.2f", round(sd, digits = 2))), size = 3.5) +
  #  geom_text(mapping = aes(y = 4, label = n_sample), size = 2.5) +
  labs(x="Mid point of time window", y="Standard deviation of generation time (days)") +
  scale_x_cut(breaks = c(1.5, 16.5), space = 0.5, which = c(1, 3), scales = c(0.125, 0.125)) +
  geom_segment(aes(x = 7, y = 5.5, xend = 7, yend = 4), 
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
    panel.grid = element_blank()
  ) + scale_y_continuous(breaks = seq(2, 5, 1), limits = c(1, 5.5)) 
#  ggtitle("Estimated standard deviation of forward generation time")

# library(aplot)

# pGT.main.bootCI <- plot_list(gglist = list(p.GT.forward.main.bootCI, p.GT.forward.sd.bootCI),
#           ncol = 1, nrow = 2, labels = c("A", "B"), tag_size = 14)

library(cowplot)

pGT.main.bootCI <- plot_grid(print(p.GT.forward.main.bootCI), print(p.GT.forward.sd.bootCI), nrow = 2,
          labels = c("a", "b"))

ggsave("Fig_forwardGTmain.pdf", width = 12.5, height = 8, units = "in", plot = pGT.main.bootCI)

