rm(list = ls())

allpairs <- read.csv("allpairs0908.csv")

library(ggplot2)
library(ggpubr)
library(fitdistrplus)

# format date

colexpo <- grep(c("expo"), colnames(allpairs)) 
for(i in colexpo){
  allpairs[, i] <- as.Date(allpairs[, i], "%m/%d/%Y")
}

colonset <- grep(c("Onset"), colnames(allpairs)) 
for(i in colonset){
  allpairs[, i] <- as.Date(allpairs[, i], "%m/%d/%Y")
}

length(unique(allpairs$Infector.ID))
length(unique(allpairs$Infectee.ID))


ind.ISboth <- which(allpairs$Infectee.ID %in% allpairs$Infector.ID)


# note in this file some infectors are repeated if they infected more than one case

ind.infector <- numeric(length(unique(allpairs$Infector.ID)))

for(i in 1:length(unique(allpairs$Infector.ID))){
  ind.infector[i] <- which(allpairs$Infector.ID == unique(allpairs$Infector.ID)[i])[1]
}

df.onset <- data.frame(
  dates = c(allpairs$Onset_Infector[ind.infector], allpairs$Onset_Infectee),
  group = c(rep("infector", length(unique(allpairs$Infector.ID))), rep("infectee", nrow(allpairs)))
)

dfinfector <- df.onset[df.onset$group == "infector",]

library(incidence)


dfinfector.incidence <- incidence(dfinfector$dates)
dfinfector.incidence$dates
dfinfector <- data.frame(
  dates = dfinfector.incidence$dates,
  counts = dfinfector.incidence$counts
)

ticks <- data.frame(
  x = seq(from = as.Date("2019-12-30"), to = as.Date("2020-02-29"), by = "1 day") - 0.5,
  y = c(rep(-1.25, 2), -5.75, rep(-1.25, 30), -5.75, rep(-1.25, 28))
)

# visualisation

p.infector.onset <- ggplot() + 
  geom_rect(data = dfinfector, 
            aes(xmin = dates - 0.5, xmax = dates + 0.5, ymin = 0, ymax = counts), fill = "#ADD8E6", color = "white") +
  scale_x_date(breaks = seq.Date(as.Date("2019-12-30"), as.Date("2020-02-29"), "day"), limits = c(as.Date("2019-12-30")-0.5, as.Date("2020-02-29")+0.5),
               labels = NULL, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0,40,10), limits = c(-10, 45), expand = c(0, 0)) + 
  ylab("Number of cases") + 
  labs(color = "") + 
  
  theme_bw() + theme(panel.background = element_blank(),
                     panel.grid = element_blank(),
                     panel.border = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.line.x = element_blank(),
                     axis.line.y = element_blank(),
                     legend.position = "top") +
  geom_hline(yintercept = 0, size = 0.75) + 
  geom_segment(data = ticks,
               aes(x = x, xend = x,
                   y = 0, yend = y), size = 0.5) +
  xlab(" ") + 
  geom_line(data = data.frame(x = c(as.Date("2020-01-23"), as.Date("2020-01-23")),
                              y = c(36, 40)), aes(x = x, y = y), size = 0.75, color = "red") +
  geom_line(data = data.frame(x = c(as.Date("2020-01-23"), as.Date("2020-01-24")),
                              y = c(40, 40)), aes(x = x, y = y), size = 0.75, color = "red") +
  geom_line(data = data.frame(x = c(as.Date("2020-01-23"), as.Date("2020-01-23")),
                              y = c(0, 35)), aes(x = x, y = y), size = 0.75, color = "red", linetype = "dotted") +
  geom_segment(aes(x = as.Date("2020-01-23"), y = 36.5, xend = as.Date("2020-01-23"), yend = 36), 
               arrow = arrow(length = unit(0.05, "inches")), lwd=0.5, color = "red") +
  annotate("text", x = as.Date("2020-01-30") + 0.5, y = 37.5, label = "Lockdown in Wuhan; \n Nation-wide control measures \n started to be implemented", size = 3.5) + 
  geom_line(data = data.frame(x = c(as.Date("2019-12-30")-0.5, 
                                    as.Date("2019-12-30")-0.5), 
                              y = c(0, Inf)), aes(x = x, y = y), size = 0.75) +
  annotate("text", x = c(as.Date("2019-12-30"),as.Date("2020-01-06"), as.Date("2020-01-13"), as.Date("2020-01-20"),
                         as.Date("2020-01-27"), as.Date("2020-02-03"), as.Date("2020-02-10"),
                         as.Date("2020-02-17"), as.Date("2020-02-24")), 
           y = rep(-2, 9), 
           label = c("30","6", "13", "20", "27", "3", "10", 
                     "17","24"), size = 3.5) +
  annotate("text", x = c(as.Date("2020-01-16"), as.Date("2020-02-15")), 
           y = rep(-5.75, 2), 
           label = c("January 2020", "February 2020"), size = 3.5) +
  ggtitle("Infector onset")
p.infector.onset


dfinfectee <- df.onset[df.onset$group == "infectee",]
dfinfectee.incidence <- incidence(dfinfectee$dates)
dfinfectee.incidence$dates
dfinfectee <- data.frame(
  dates = dfinfectee.incidence$dates,
  counts = dfinfectee.incidence$counts
)

ticks <- data.frame(
  x = seq(from = as.Date("2019-12-30"), to = as.Date("2020-02-29"), by = "1 day") - 0.5,
  y = c(rep(-1.3, 2), -6, rep(-1.3, 30), -6, rep(-1.3, 28))
)

p.infectee.onset <- ggplot() + 
  geom_rect(data = dfinfectee, 
            aes(xmin = dates - 0.5, xmax = dates + 0.5, ymin = 0, ymax = counts), fill = "#E9967A", color = "white") +
  scale_x_date(breaks = seq.Date(as.Date("2019-12-30"), as.Date("2020-02-29"), "day"), limits = c(as.Date("2019-12-30")-0.5, as.Date("2020-02-29")+0.5),
               labels = NULL, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0,50,10), limits = c(-10, 50), expand = c(0, 0)) + 
  ylab("Number of cases") + 
  labs(color = "") + 
  geom_line(data = data.frame(x = c(as.Date("2020-01-23"), as.Date("2020-01-23")),
                              y = c(0, 50)), aes(x = x, y = y), size = 0.75, color = "red", linetype = "dotted") +
  theme_bw() + theme(panel.background = element_blank(),
                     panel.grid = element_blank(),
                     panel.border = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.line.x = element_blank(),
                     axis.line.y = element_blank(),
                     axis.title = element_text(),
                     legend.position = "top") +
  geom_hline(yintercept = 0, size = 0.75) + 
  geom_segment(data = ticks,
               aes(x = x, xend = x,
                   y = 0, yend = y), size = 0.5) +
  xlab(" ") + 
  geom_line(data = data.frame(x = c(as.Date("2019-12-30")-0.5, 
                                    as.Date("2019-12-30")-0.5), 
                              y = c(0, Inf)), aes(x = x, y = y), size = 0.75) +
  annotate("text", x = c(as.Date("2019-12-30"),as.Date("2020-01-06"), as.Date("2020-01-13"), as.Date("2020-01-20"),
                         as.Date("2020-01-27"), as.Date("2020-02-03"), as.Date("2020-02-10"),
                         as.Date("2020-02-17"), as.Date("2020-02-24")), 
           y = rep(-2.25, 9), 
           label = c("30","6", "13", "20", "27", "3", "10", 
                     "17","24"), size = 3.5) +
  annotate("text", x = c(as.Date("2020-01-16"), as.Date("2020-02-15")), 
           y = rep(-6, 2), 
           label = c("January 2020", "February 2020"), size = 3.5) +
  ggtitle("Infectee onset")
  
p.infectee.onset  
  
  
ponset.strtf <- ggarrange(p.infector.onset, p.infectee.onset, 
                      nrow = 2, ncol = 1, align = "v", labels = c("a", "b"))

ggsave("Fig_onset_incidence.pdf", ponset.strtf, width = 12.5, height = 8, units = "in")

length(which(df.onset$dates[df.onset$group == "infector"] >= as.Date("2020-01-23") &
               df.onset$dates[df.onset$group == "infector"] <= as.Date("2020-01-29"))) # 199

length(which(df.onset$dates[df.onset$group == "infector"] < as.Date("2020-01-23"))) # 107

length(which(df.onset$dates[df.onset$group == "infector"] > as.Date("2020-01-29"))) # 122

length(which(df.onset$dates[df.onset$group == "infectee"] >= as.Date("2020-01-29") & 
               df.onset$dates[df.onset$group == "infectee"] <= as.Date("2020-02-03"))) # 231

##############################################################################################################

########################## Incubation period of infector and infectee (referenced by infector onset)

##### infector

length(which(is.na(allpairs[ind.infector,]$expo.early.I))) # 299
length(which(is.na(allpairs[ind.infector,]$expo.late.I))) # 167
length(which(is.na(allpairs[ind.infector,]$expo.early.I) & is.na(allpairs[ind.infector,]$expo.late.I))) # 164
length(ind.infector) # 428

ind.fullinfo.infector <- which((!is.na(allpairs[ind.infector,]$expo.early.I) & 
                                  !is.na(allpairs[ind.infector,]$expo.late.I)))



# only 126 infectors have full expo info

fullinfo.infector <- allpairs[ind.infector,][ind.fullinfo.infector,]



# check if there is any infector who has earliest possible exposure and onset at same date
which(fullinfo.infector$expo.early.I == fullinfo.infector$Onset_Infector)

# make this special case's earliest possible exposure point one day earlier (to ensure IP > 0)
fullinfo.infector[77,]$expo.early.I <- fullinfo.infector[77,]$Onset_Infector - 1

mean(fullinfo.infector$expo.late.I - fullinfo.infector$expo.early.I) # 3.4 days

startseq <- c(as.Date("2020-01-01"), seq.Date(from = as.Date("2020-01-15"), to = as.Date("2020-01-30"), by = "1 day"))
endseq <- c(as.Date("2020-01-20"), seq.Date(from = as.Date("2020-01-21"), to = as.Date("2020-02-04"), by = "1 day"),
            as.Date("2020-02-29"))

# temporal estimates of IP
n <- length(startseq)


for(i in 1:n){
  
  # i <- 1
  subdata <- subset(fullinfo.infector, Onset_Infector >= startseq[i] & Onset_Infector <= endseq[i])
  print(mean(subdata$expo.late.I - subdata$expo.early.I))
}


par1 <- par2 <- vector("list", n)
par1.boot <- par2.boot <- vector("list", n)
AIC <- vector("list", n)
samplesize <- numeric(n)

progressbar <- txtProgressBar(min = 0, max = n, style = 3)

for(i in 1:n){
  
  # i <- 1
  subdata <- subset(fullinfo.infector, Onset_Infector >= startseq[i] & Onset_Infector <= endseq[i])
  samplesize[i] <- nrow(subdata)
  dftmp <- data.frame(left = with(subdata, as.numeric(Onset_Infector - expo.late.I)),
                      right = with(subdata, as.numeric(Onset_Infector - expo.early.I)))
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

IPinfector.reslist.back <- list(
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

save(IPinfector.reslist.back, file = "IPinfector.reslist.back.RData")

shape = unlist(IPinfector.reslist.back$par1)[seq(2, 50, 3)]
scale = unlist(IPinfector.reslist.back$par2)[seq(2, 50, 3)]
shapeboot = lapply(IPinfector.reslist.back$par1boot, function(x) unlist(x[[2]]))
scaleboot = lapply(IPinfector.reslist.back$par2boot, function(x) unlist(x[[2]]))

muUB <- muLB <- sigmaUB <- sigmaLB <- numeric(17)

for(i in 1:17){
  shapeboot <- unlist(IPinfector.reslist.back$par1boot[[i]][[2]])
  scaleboot <- unlist(IPinfector.reslist.back$par2boot[[i]][[2]])
  estboot <- weibullparinv(shapeboot, scaleboot)
  muUB[i] <- quantile(estboot$mu, 0.975)
  muLB[i] <- quantile(estboot$mu, 0.025)
  sigmaUB[i] <- quantile(estboot$sigma, 0.975)
  sigmaLB[i] <- quantile(estboot$sigma, 0.025)
}

ests <- weibullparinv(shape, scale)

dfIP.backward <- data.frame(
  alpha = shape,
  beta = scale,
  mu = ests$mu,
  sigma = ests$sigma,
  samplesize = IPinfector.reslist.back$samplesize,
  timewindow = IPinfector.reslist.back$timewindow,
  muUB = muUB,
  muLB = muLB,
  sigmaUB = sigmaUB,
  sigmaLB = sigmaLB
)

write.csv(dfIP.backward, "dfIPinfector.backward_fitdistWB.csv")


IPinfector.reslist.back$par1
IPinfector.reslist.back$par2
IPinfector.reslist.back$AIC
IPinfector.reslist.back$timewindow

AICs <- unlist(IPinfector.reslist.back$AIC)


df.AIC.backinfector <- data.frame(
  timewindow = rep(IPinfector.reslist.back$timewindow, 3),
  AIC = c(AICs[seq(1, 49, 3)], AICs[seq(2, 50, 3)], AICs[seq(3, 51, 3)]),
  dist = c(rep("Gamma", 17), rep("Weibull", 17), rep("Log-Normal", 17)),
  samplesize = rep(IPinfector.reslist.back$samplesize, 3)
)

write.csv(df.AIC.backinfector, "df.AIC.backinfector.csv")


##### infectee
length(which(is.na(allpairs$expo.late.S) & is.na(allpairs$expo.early.S))) # 259


ind.fullinfo.infectee <- which((!is.na(allpairs$expo.early.S) & 
                                  !is.na(allpairs$expo.late.S)))

fullinfo.infectee <- allpairs[ind.fullinfo.infectee,]

mean(fullinfo.infectee$expo.late.S - fullinfo.infectee$expo.early.S) # 5.9 days

# temporal estimates of IP
n <- length(startseq)

for(i in 1:n){
  
  subdata <- subset(fullinfo.infectee, Onset_Infector >= startseq[i] & Onset_Infector <= endseq[i])
  print(mean(subdata$expo.late.S - subdata$expo.early.S))
}



par1 <- par2 <- vector("list", n)
par1.boot <- par2.boot <- vector("list", n)
AIC <- vector("list", n)
samplesize <- numeric(n)

progressbar <- txtProgressBar(min = 0, max = n, style = 3)

for(i in 1:n){
  
  subdata <- subset(fullinfo.infectee, Onset_Infector >= startseq[i] & Onset_Infector <= endseq[i])
  samplesize[i] <- nrow(subdata)
  dftmp <- data.frame(left = with(subdata, as.numeric(Onset_Infectee - expo.late.S)),
                      right = with(subdata, as.numeric(Onset_Infectee - expo.early.S)))
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

IPinfectee.reslist.for <- list(
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

save(IPinfectee.reslist.for, file = "IPinfectee.reslist.for.RData")

AICs <- unlist(IPinfectee.reslist.for$AIC)


df.AIC.forwardinfectee <- data.frame(
  timewindow = rep(IPinfectee.reslist.for$timewindow, 3),
  AIC = c(AICs[seq(1, 49, 3)], AICs[seq(2, 50, 3)], AICs[seq(3, 51, 3)]),
  dist = c(rep("Gamma", 17), rep("Weibull", 17), rep("Log-Normal", 17)),
  samplesize = rep(IPinfectee.reslist.for$samplesize, 3)
)

write.csv(df.AIC.forwardinfectee, "df.AIC.forwardinfectee.csv")


shape = unlist(IPinfectee.reslist.for$par1)[seq(2, 50, 3)]
scale = unlist(IPinfectee.reslist.for$par2)[seq(2, 50, 3)]
shapeboot = lapply(IPinfectee.reslist.for$par1boot, function(x) unlist(x[[2]]))
scaleboot = lapply(IPinfectee.reslist.for$par2boot, function(x) unlist(x[[2]]))

muUB <- muLB <- sigmaUB <- sigmaLB <- numeric(17)

for(i in 1:17){
  shapeboot <- unlist(IPinfectee.reslist.for$par1boot[[i]][[2]])
  scaleboot <- unlist(IPinfectee.reslist.for$par2boot[[i]][[2]])
  estboot <- weibullparinv(shapeboot, scaleboot)
  muUB[i] <- quantile(estboot$mu, 0.975)
  muLB[i] <- quantile(estboot$mu, 0.025)
  sigmaUB[i] <- quantile(estboot$sigma, 0.975)
  sigmaLB[i] <- quantile(estboot$sigma, 0.025)
}

ests <- weibullparinv(shape, scale)

dfIP.forward <- data.frame(
  alpha = shape,
  beta = scale,
  mu = ests$mu,
  sigma = ests$sigma,
  samplesize = IPinfectee.reslist.for$samplesize,
  timewindow = IPinfectee.reslist.for$timewindow,
  muUB = muUB,
  muLB = muLB,
  sigmaUB = sigmaUB,
  sigmaLB = sigmaLB
)

write.csv(dfIP.forward, "dfIPinfectee.forward_fitdistWB.csv")


##############################################################################################################

########################## Incubation period of infector and infectee (referenced by infectee onset)

startseq <- c(as.Date("2020-01-01"), seq.Date(from = as.Date("2020-01-21"), to = as.Date("2020-02-05"), by = "1 day"))
endseq <- c(as.Date("2020-01-26"), seq.Date(from = as.Date("2020-01-27"), to = as.Date("2020-02-10"), by = "1 day"),
            as.Date("2020-02-29"))

# temporal estimates of IP
n <- length(startseq)

par1 <- par2 <- vector("list", n)
par1.boot <- par2.boot <- vector("list", n)
AIC <- vector("list", n)
samplesize <- numeric(n)

progressbar <- txtProgressBar(min = 0, max = n, style = 3)

for(i in 1:n){
  
  subdata <- subset(fullinfo.infector, Onset_Infectee >= startseq[i] & Onset_Infectee <= endseq[i])
  samplesize[i] <- nrow(subdata)
  dftmp <- data.frame(left = with(subdata, as.numeric(Onset_Infector - expo.late.I)),
                      right = with(subdata, as.numeric(Onset_Infector - expo.early.I)))
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

IPinfector.reslist.for <- list(
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

IPinfector.reslist.for$samplesize

save(IPinfector.reslist.for, file = "IPinfector.reslist.for.RData")

IPinfector.reslist.for$par1
IPinfector.reslist.for$par2
IPinfector.reslist.for$AIC
IPinfector.reslist.for$timewindow

AICs <- unlist(IPinfector.reslist.for$AIC)


df.AIC.forinfector <- data.frame(
  timewindow = rep(IPinfector.reslist.for$timewindow, 3),
  AIC = c(AICs[seq(1, 49, 3)], AICs[seq(2, 50, 3)], AICs[seq(3, 51, 3)]),
  dist = c(rep("Gamma", 17), rep("Weibull", 17), rep("Log-Normal", 17)),
  samplesize = rep(IPinfector.reslist.for$samplesize, 3)
)

write.csv(df.AIC.forinfector, "df.AIC.forwardinfector.csv")

##### infectee

# temporal estimates of IP
n <- length(startseq)

par1 <- par2 <- vector("list", n)
par1.boot <- par2.boot <- vector("list", n)
AIC <- vector("list", n)
samplesize <- numeric(n)

progressbar <- txtProgressBar(min = 0, max = n, style = 3)

for(i in 1:n){
  
  subdata <- subset(fullinfo.infectee, Onset_Infectee >= startseq[i] & Onset_Infectee <= endseq[i])
  samplesize[i] <- nrow(subdata)
  dftmp <- data.frame(left = with(subdata, as.numeric(Onset_Infectee - expo.late.S)),
                      right = with(subdata, as.numeric(Onset_Infectee - expo.early.S)))
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

IPinfectee.reslist.back <- list(
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

IPinfectee.reslist.back$samplesize

save(IPinfectee.reslist.back, file = "IPinfectee.reslist.back.RData")

IPinfectee.reslist.back$par1
IPinfectee.reslist.back$par2
IPinfectee.reslist.back$AIC
IPinfectee.reslist.back$timewindow

AICs <- unlist(IPinfectee.reslist.back$AIC)


df.AIC.backinfectee <- data.frame(
  timewindow = rep(IPinfectee.reslist.back$timewindow, 3),
  AIC = c(AICs[seq(1, 49, 3)], AICs[seq(2, 50, 3)], AICs[seq(3, 51, 3)]),
  dist = c(rep("Gamma", 17), rep("Weibull", 17), rep("Log-Normal", 17)),
  samplesize = rep(IPinfectee.reslist.back$samplesize, 3)
)

write.csv(df.AIC.backinfectee, "df.AIC.backwardinfectee.csv")

shape = unlist(IPinfectee.reslist.back$par1)[seq(2, 50, 3)]
scale = unlist(IPinfectee.reslist.back$par2)[seq(2, 50, 3)]
shapeboot = lapply(IPinfectee.reslist.back$par1boot, function(x) unlist(x[[2]]))
scaleboot = lapply(IPinfectee.reslist.back$par2boot, function(x) unlist(x[[2]]))

muUB <- muLB <- sigmaUB <- sigmaLB <- numeric(17)

for(i in 1:17){
  shapeboot <- unlist(IPinfectee.reslist.back$par1boot[[i]][[2]])
  scaleboot <- unlist(IPinfectee.reslist.back$par2boot[[i]][[2]])
  estboot <- weibullparinv(shapeboot, scaleboot)
  muUB[i] <- quantile(estboot$mu, 0.975)
  muLB[i] <- quantile(estboot$mu, 0.025)
  sigmaUB[i] <- quantile(estboot$sigma, 0.975)
  sigmaLB[i] <- quantile(estboot$sigma, 0.025)
}

ests <- weibullparinv(shape, scale)

dfIP.backward <- data.frame(
  alpha = shape,
  beta = scale,
  mu = ests$mu,
  sigma = ests$sigma,
  samplesize = IPinfectee.reslist.back$samplesize,
  timewindow = IPinfectee.reslist.back$timewindow,
  muUB = muUB,
  muLB = muLB,
  sigmaUB = sigmaUB,
  sigmaLB = sigmaLB
)

write.csv(dfIP.backward, "dfIPinfectee.backward_fitdistWB.csv")



shape = unlist(IPinfector.reslist.for$par1)[seq(1, 49, 3)]
rate = unlist(IPinfector.reslist.for$par2)[seq(1, 49, 3)]
shapeboot = lapply(IPinfector.reslist.for$par1boot, function(x) unlist(x[[1]]))
rateboot = lapply(IPinfector.reslist.for$par2boot, function(x) unlist(x[[1]]))

muUB <- muLB <- sigmaUB <- sigmaLB <- numeric(17)

for(i in 1:17){
  shapeboot <- unlist(IPinfector.reslist.for$par1boot[[i]][[1]])
  rateboot <- unlist(IPinfector.reslist.for$par2boot[[i]][[1]])
  muUB[i] <- quantile(shapeboot/rateboot, 0.975)
  muLB[i] <- quantile(shapeboot/rateboot, 0.025)
  sigmaUB[i] <- quantile(sqrt(shapeboot/rateboot^2), 0.975)
  sigmaLB[i] <- quantile(sqrt(shapeboot/rateboot^2), 0.025)
}


dfIP.forward <- data.frame(
  alpha = shape,
  beta = rate,
  mu = shape/rate,
  sigma = sqrt(shape/rate^2),
  samplesize = IPinfector.reslist.for$samplesize,
  timewindow = IPinfector.reslist.for$timewindow,
  muUB = muUB,
  muLB = muLB,
  sigmaUB = sigmaUB,
  sigmaLB = sigmaLB
)

write.csv(dfIP.forward, "dfIPinfector.forward_fitdistG.csv")

