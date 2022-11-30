rm(list = ls())

allpairs <- read.csv("allpairs0908.csv")


library(ggplot2)
library(ggpubr)
library(fitdistrplus)
library(incidence)
library(EpiEstim)

# format date

colexpo <- grep(c("expo"), colnames(allpairs)) 
for(i in colexpo){
  allpairs[, i] <- as.Date(allpairs[, i], "%m/%d/%Y")
}

colonset <- grep(c("Onset"), colnames(allpairs)) 
for(i in colonset){
  allpairs[, i] <- as.Date(allpairs[, i], "%m/%d/%Y")
}

####### temporal SI

startseq <- c(as.Date("2020-01-01"), seq.Date(from = as.Date("2020-01-15"), to = as.Date("2020-01-30"), by = "1 day"))
endseq <- c(as.Date("2020-01-20"), seq.Date(from = as.Date("2020-01-21"), to = as.Date("2020-02-04"), by = "1 day"),
            as.Date("2020-02-29"))

muSI <- sdSI <- numeric(17)
muSI.LB <- muSI.UB <- sdSI.LB <- sdSI.UB <- numeric(17)
samplesize <- numeric(17)

for(i in 1:17){
  tmppair <- subset(allpairs, Onset_Infector >= startseq[i] & Onset_Infector <= endseq[i])
  samplesize[i] <- nrow(tmppair)
  SIdata <- as.numeric(tmppair$Onset_Infectee - tmppair$Onset_Infector)
  tmpfit <- fitdist(SIdata, "norm")
  muSI[i] <- tmpfit$estimate[1]
  sdSI[i] <- tmpfit$estimate[2]
  tmpboot <- bootdist(tmpfit)
  muSI.LB[i] <- quantile(tmpboot$estim[, 1], 0.025)
  muSI.UB[i] <- quantile(tmpboot$estim[, 1], 0.975)
  sdSI.LB[i] <- quantile(tmpboot$estim[, 2], 0.025)
  sdSI.UB[i] <- quantile(tmpboot$estim[, 2], 0.975)
}

slidSI.df <- data.frame(
  timewindow = paste0(startseq, " - ", endseq),
  X = seq(1, 17, 1),
  mu = muSI,
  muLB = muSI.LB,
  muUB = muSI.UB,
  sd = sdSI,
  sdLB = sdSI.LB,
  sdUB = sdSI.UB,
  samplesize = samplesize
)

write.csv(slidSI.df, "slidSI_7dslid_20220506.csv")

# plot
library(ggbreak)
library(cowplot)

slidSI.df <- read.csv("slidSI_7dslid_20220506.csv")

p.SI.forward <- ggplot(slidSI.df, aes(x = X, y = mu, ymin = muLB, ymax = muUB)) +
  geom_point(size = 1.0) + 
  geom_errorbar(width = 0, size = 0.5) +
  geom_text(mapping = aes(y = muUB+0.30, label = sprintf("%0.2f", round(mu, digits = 2))), size = 3.0)+
  labs(x="Mid point of time window", y="Mean serial interval (days)") +
  scale_x_continuous(limits = c(0, 18), breaks = seq(1, 17, 1), labels =
                       c("Jan 10",
                         "Jan 18", "Jan 19", "Jan 20", "Jan 21", "Jan 22",
                         "Jan 23", "Jan 24", "Jan 25", "Jan 26", "Jan 27", "Jan 28",
                         "Jan 29", "Jan 30", "Jan 31", "Feb 1",
                         "Feb 14")) +
  geom_segment(aes(x = 7, y = 10, xend = 7, yend = 8), 
               arrow = arrow(length = unit(0.05, "inches")), lwd=0.5, color = "red") +
  geom_segment(aes(x = 7, y = 10, xend = 7.5, yend = 10), 
               lwd=0.5, color = "red") +
  annotate("text", x = 9 + 0.5, y = 10, label = "Lockdown in Wuhan", size = 3.0) + 
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
    axis.text = element_text(size=7.5),
    axis.line = element_line(color="black", size=0.5),
    panel.grid = element_blank()
  ) + scale_y_continuous(breaks = seq(2, 10, 1),limits = c(1.5, 10.5))

p.SI.forward.sd <- ggplot(slidSI.df, aes(x = X, y = sd, ymin = sdLB, ymax = sdUB)) +
  geom_point(size = 1.0) + 
  geom_errorbar(width = 0, size = 0.5) +
  geom_segment(aes(x = 7, y = 6.5, xend = 7, yend = 5.8), 
               arrow = arrow(length = unit(0.05, "inches")), lwd=0.5, color = "red") +
  geom_text(mapping = aes(y = sdUB+0.30, label = sprintf("%0.2f", round(sd, digits = 2))), size = 3.0) +
  labs(x="Mid point of time window", y="Standard deviation of serial interval (days)") +
  scale_x_continuous(limits = c(0, 18), breaks = seq(1, 17, 1), labels =
                       c("Jan 10",
                         "Jan 18", "Jan 19", "Jan 20", "Jan 21", "Jan 22",
                         "Jan 23", "Jan 24", "Jan 25", "Jan 26", "Jan 27", "Jan 28",
                         "Jan 29", "Jan 30", "Jan 31", "Feb 1",
                         "Feb 14")) +
  
  scale_x_cut(breaks = c(1.5, 16.5), space = 0.5, which = c(1, 3), scales = c(0.125, 0.125)) + 
  theme(
    axis.line.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.text.x.top = element_blank(),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    axis.text.x = element_text(),
    axis.text = element_text(size=7.5),
    axis.line = element_line(color="black", size=0.5),
    panel.grid = element_blank()
  ) + scale_y_continuous(breaks = seq(4, 6, 0.5), limits = c(3.5, 6.5)) 

pSI <- plot_grid(print(p.SI.forward), print(p.SI.forward.sd), nrow = 2,
                 labels = c("a", "b"))

ggsave("Fig_forwardSI_est.pdf", width = 12.5, height = 8, units = "in", plot = pSI)

ggsave("Fig_S6_final_1128.pdf", width = 8.25, height = 8.5, units = "in", plot = pSI)


########### estimate Rt
GTest <- read.csv("forwardGT_lnorm_May06.csv")
slidSI <- read.csv("slidSI_7dslid_20220506.csv")


### onset time window

# for onset, time windows are
# Jan 08 - 14; Jan 09 - 15; ..., Jan 14 - 20 (use the first SI/GT estimate)
# start = seq(8, 14); end = seq(14, 20)
# Jan 15 - 21; Jan 16 - 22; ..., Jan 29 - Feb 04 (use the corresponding sliding window)
# start = seq(15, 29); end = seq(21, 35)
# Jan 30 - Feb 05; Jan 31 - Feb 06; ..., Feb 11 - Feb 17 (use the last SI/GT estimate)
# start = seq(30, 42); end = seq(36, 48)

# a total of 7 + 15 + 13 = 35 time windows
infector.Tonset.df <- data.frame(
  ID = allpairs$Infector.ID,
  onset = allpairs$Onset_Infector
)

infectee.Tonset.df <- data.frame(
  ID = allpairs$Infectee.ID,
  onset = allpairs$Onset_Infectee
)

all.Tonset.df <- rbind(infector.Tonset.df,
                       infectee.Tonset.df)

ID.uniq <- unique(all.Tonset.df$ID)
Tonset.all <- numeric(length(ID.uniq))

for(i in 1:length(ID.uniq)){
  subdata <- subset(all.Tonset.df, ID == ID.uniq[i])
  Tonset.all[i] <- as.character(subdata$onset[1])
}

Tonset.all <- as.Date(Tonset.all)

all.onset.inc <- incidence(Tonset.all)

all.onset.inc$dates







Rt.GT <- Rt.GT.LB <- Rt.GT.UB <- numeric(35)

all.onset.inc$dates

all.onset.inc$counts

res1 <- wallinga_teunis(
  incid = all.onset.inc,
  method = "parametric_si",
  config = 
    list(
      t_start = seq(8, 14),
      t_end =  seq(14, 20), 
      mean_si = GTest$mu[1],
      std_si = GTest$sd[1],
      n_sim = 100
    )
)

Rt.GT[1:7] <- res1$R$`Mean(R)`[1:7]  
Rt.GT.LB[1:7] <- res1$R$`Quantile.0.025(R)`[1:7] 
Rt.GT.UB[1:7] <- res1$R$`Quantile.0.975(R)`[1:7]  

for(i in 1:15){
  restmp <- wallinga_teunis(
    incid = all.onset.inc,
    method = "parametric_si",
    config = 
      list(
        t_start = seq(15, 29),
        t_end =  seq(21, 35), 
        mean_si = GTest$mu[i+1],
        std_si = GTest$sd[i+1],
        n_sim = 100
      )
  )
  
  Rt.GT[i+7] <- restmp$R$`Mean(R)`[i]  
  Rt.GT.LB[i+7] <- restmp$R$`Quantile.0.025(R)`[i] 
  Rt.GT.UB[i+7] <- restmp$R$`Quantile.0.975(R)`[i]  
  
  rm(restmp)
}

Rt.GT[8:22]


res2 <- wallinga_teunis(
  incid = all.onset.inc,
  method = "parametric_si",
  config = 
    list(
      t_start = seq(30, 42),
      t_end =  seq(36, 48), 
      mean_si = GTest$mu[17],
      std_si = GTest$sd[17],
      n_sim = 100
    )
)


Rt.GT[23:35] <- res2$R$`Mean(R)`[1:13]  
Rt.GT.LB[23:35] <- res2$R$`Quantile.0.025(R)`[1:13] 
Rt.GT.UB[23:35] <- res2$R$`Quantile.0.975(R)`[1:13]  


### onset time window

# for onset, time windows are
# Jan 08 - 14; Jan 09 - 15; ..., Jan 14 - 20 (use the first SI estimate)
# start = seq(8, 14); end = seq(14, 20)
# Jan 15 - 21; Jan 16 - 22; ..., Jan 29 - Feb 04 (use the corresponding sliding window)
# start = seq(15, 29); end = seq(21, 35)
# Jan 30 - Feb 05; Jan 31 - Feb 06; ..., Feb 11 - Feb 17 (use the last SI estimate)
# start = seq(30, 42); end = seq(36, 48)

# a total of 7 + 15 + 13 = 35 time windows

Rt.SI <- Rt.SI.LB <- Rt.SI.UB <- numeric(35)

all.onset.inc$dates

res3 <- wallinga_teunis(
  incid = all.onset.inc,
  method = "parametric_si",
  config = 
    list(
      t_start = seq(8, 14),
      t_end =  seq(14, 20), 
      mean_si = slidSI$mu[1],
      std_si = slidSI$sd[1],
      n_sim = 100
    )
)


Rt.SI[1:7] <- res3$R$`Mean(R)`[1:7]  
Rt.SI.LB[1:7] <- res3$R$`Quantile.0.025(R)`[1:7] 
Rt.SI.UB[1:7] <- res3$R$`Quantile.0.975(R)`[1:7]  

for(i in 1:15){
  restmp <- wallinga_teunis(
    incid = all.onset.inc,
    method = "parametric_si",
    config = 
      list(
        t_start = seq(15, 29),
        t_end =  seq(21, 35), 
        mean_si = slidSI$mu[i+1],
        std_si = slidSI$sd[i+1],
        n_sim = 100
      )
  )
  
  
  Rt.SI[i+7] <- restmp$R$`Mean(R)`[i]  
  Rt.SI.LB[i+7] <- restmp$R$`Quantile.0.025(R)`[i] 
  Rt.SI.UB[i+7] <- restmp$R$`Quantile.0.975(R)`[i]  
  
  rm(restmp)
}

Rt.SI[8:22]


res4 <- wallinga_teunis(
  incid = all.onset.inc,
  method = "parametric_si",
  config = 
    list(
      t_start = seq(30, 42),
      t_end =  seq(36, 48), 
      mean_si = slidSI$mu[17],
      std_si = slidSI$sd[17],
      n_sim = 100
    )
)


Rt.SI[23:35] <- res4$R$`Mean(R)`[1:13]  
Rt.SI.LB[23:35] <- res4$R$`Quantile.0.025(R)`[1:13] 
Rt.SI.UB[23:35] <- res4$R$`Quantile.0.975(R)`[1:13]  

Rt.sum <- data.frame(
  Rt.GT =  Rt.GT,
  Rt.GT.LB = Rt.GT.LB,
  Rt.GT.UB = Rt.GT.UB,
  Rt.SI =  Rt.SI,
  Rt.SI.LB = Rt.SI.LB,
  Rt.SI.UB = Rt.SI.UB,
  Dates.end = seq.Date(from = as.Date("2020-01-14"), to = as.Date("2020-02-17"), by = "1 day")
)


Rt.sum$Dates.end <- as.Date(Rt.sum$Dates.end)

ticks <- data.frame(
  x = seq(from = as.Date("2019-12-30"), to = as.Date("2020-02-29"), by = "1 day") - 0.5,
  y = c(rep(-0.1, 2), -0.3, rep(-0.1, 30), -0.3, rep(-0.1, 28))
)

p.Rt <- ggplot() + 
  geom_line(data = Rt.sum[5:35,], 
            aes(x = Dates.end, y = Rt.GT, color = "GT")) +
  geom_ribbon(data = Rt.sum[5:35,], aes(ymin = Rt.GT.LB, ymax = Rt.GT.UB, x = Dates.end), fill = "#F8766D", alpha = 0.3)+
  geom_line(data = Rt.sum[5:35,], 
            aes(x = Dates.end, y = Rt.SI, color = "SI")) +
  geom_ribbon(data = Rt.sum[5:35,], aes(ymin = Rt.SI.LB, ymax = Rt.SI.UB, x = Dates.end), fill = "#00BFC4", alpha = 0.3)+
  scale_x_date(breaks = seq.Date(as.Date("2019-12-30"), as.Date("2020-02-29"), "day"), limits = c(as.Date("2019-12-30")-0.5, as.Date("2020-02-29")+0.5),
               labels = NULL, expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.35, 4.5), breaks = seq(0, 4.0, 0.5)) + 
  theme_bw() + theme(panel.background = element_blank(),
                     panel.grid = element_blank(),
                     panel.border = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.line.x = element_blank(),
                     axis.line.y = element_blank(),
                     legend.position = c(0.8, 0.8),
                     legend.title = element_blank()) +
  geom_hline(yintercept = 0, size = 0.75) + 
  geom_hline(yintercept = 1, size = 0.5, color = "darkgrey", linetype = "dashed") + 
  geom_segment(data = ticks,
               aes(x = x, xend = x,
                   y = 0, yend = y), size = 0.5) +
  xlab(" ") + ylab("Rt") +
  geom_line(data = data.frame(x = c(as.Date("2019-12-30")-0.5, 
                                    as.Date("2019-12-30")-0.5), 
                              y = c(0, Inf)), aes(x = x, y = y), size = 1.25) +
  geom_line(data = data.frame(x = c(as.Date("2020-01-23"), 
                                    as.Date("2020-01-23")), 
                              y = c(0, Inf)), aes(x = x, y = y), size = 0.75, linetype = "dotted",
            color = "red") +
  annotate("text", x = c(as.Date("2019-12-30"),as.Date("2020-01-06"), as.Date("2020-01-13"), as.Date("2020-01-20"),
                         as.Date("2020-01-27"), as.Date("2020-02-03"), as.Date("2020-02-10"),
                         as.Date("2020-02-17"), as.Date("2020-02-24")), 
           y = rep(-0.15, 9), 
           label = c("30","6", "13", "20", "27", "3", "10", 
                     "17","24"), size = 3.5) +
  annotate("text", x = c(as.Date("2020-01-16"), as.Date("2020-02-15")), 
           y = rep(-0.35, 2), 
           label = c("January 2020", "February 2020"), size = 3.5) 
p.Rt

ticks <- data.frame(
  x = seq(from = as.Date("2019-12-30"), to = as.Date("2020-02-29"), by = "1 day") - 0.5,
  y = c(rep(-1.25, 2), -5.75, rep(-1.25, 30), -5.75, rep(-1.25, 28))
)

# visualisation

p.epi <- ggplot() + 
  geom_rect(data = dfTonset, 
            aes(xmin = dates - 0.5, xmax = dates + 0.5, ymin = 0, ymax = counts), color = "white", alpha = 0.5) +
  scale_x_date(breaks = seq.Date(as.Date("2019-12-30"), as.Date("2020-02-29"), "day"), limits = c(as.Date("2019-12-30") - 0.5, as.Date("2020-02-29")+0.5),
               labels = NULL, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 90, 10), limits = c(-10, 90), expand = c(0, 0)) + 
  ylab("Number of cases") + 
  theme_bw() + theme(panel.background = element_blank(),
                     panel.grid = element_blank(),
                     panel.border = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.line.x = element_blank(),
                     axis.line.y = element_blank(),
                     legend.position = c(0.8, 0.8),
                     legend.title = element_blank()) +
  geom_hline(yintercept = 0, size = 0.75) + 
  geom_segment(data = ticks,
               aes(x = x, xend = x,
                   y = 0, yend = y), size = 0.5) +
  xlab(" ") + 
  geom_line(data = data.frame(x = c(as.Date("2020-01-23"), as.Date("2020-01-23")),
                              y = c(80, 85)), aes(x = x, y = y), size = 0.75, color = "red") +
  geom_line(data = data.frame(x = c(as.Date("2020-01-23"), as.Date("2020-01-25")),
                              y = c(85, 85)), aes(x = x, y = y), size = 0.75, color = "red") +
  geom_line(data = data.frame(x = c(as.Date("2020-01-23"), as.Date("2020-01-23")),
                              y = c(0, 80)), aes(x = x, y = y), size = 0.75, color = "red", linetype = "dotted") +
  geom_segment(aes(x = as.Date("2020-01-23"), y = 85, xend = as.Date("2020-01-23"), yend = 80), 
               arrow = arrow(length = unit(0.05, "inches")), lwd=0.5, color = "red") +
  annotate("text", x = as.Date("2020-02-01") + 0.5, y = 82, label = "Lockdown in Wuhan; \n Nation-wide control measures \n started to be implemented", size = 3.5) + 
  geom_line(data = data.frame(x = c(as.Date("2019-12-30") - 0.5, 
                                    as.Date("2019-12-30") - 0.5),
                              y = c(0, Inf)), aes(x = x, y = y), size = 1.25) +
  annotate("text", x = c(as.Date("2019-12-30"),as.Date("2020-01-06"), as.Date("2020-01-13"), as.Date("2020-01-20"),
                         as.Date("2020-01-27"), as.Date("2020-02-03"), as.Date("2020-02-10"),
                         as.Date("2020-02-17"), as.Date("2020-02-24")), 
           y = rep(-2.75, 9), 
           label = c("30","6", "13", "20", "27", "3", "10", 
                     "17","24"), size = 3.5) +
  annotate("text", x = c(as.Date("2020-01-16"), as.Date("2020-02-15")), 
           y = rep(-7, 2), 
           label = c("January 2020", "February 2020"), size = 3.5) 
p.epi

p.fig5_v1 <- ggarrange(p.epi, p.Rt, nrow = 2, ncol = 1, align = "v", heights = c(1.25, 1),
                       labels = c("a", "b"))
p.fig5_v1

ggsave("Fig_epicurve_Rt.pdf", plot = p.fig5_v1, width = 12.5, height = 8, unit = "in")


