library(ggplot2)
library(ggpubr)
library(ggbreak)
library(cowplot)

dfempSI <- read.csv("forwardempSI.csv")
tw.uniq <- unique(dfempSI$timewindow)

empmean <- numeric(17)
for(i in 1:17){
  tmp <- subset(dfempSI, timewindow == tw.uniq[i])
  empmean[i] <- mean(tmp$empSI)
}
empmean

empq1 <- numeric(17)
empq3 <- numeric(17)


for(i in 1:17){
  tmp <- subset(dfempSI, timewindow == tw.uniq[i])
  empq1[i] <- quantile(tmp$empSI, 0.25)
  empq3[i] <- quantile(tmp$empSI, 0.75)
}

empq1
empq3


df.empSI.forward <- data.frame(
  X = seq(1, 17),
  empmean = empmean,
  empq1 = empq1,
  empq3 = empq3
)


p.empSI.forward <- ggplot(df.empSI.forward, aes(x = X, y = empmean, ymin = empq1, ymax = empq3)) +
  geom_point(size = 1.0) + 
  geom_errorbar(width = 0, size = 0.5) +
  annotate("text", y = df.empSI.forward$empmean + 0.5, x = df.empSI.forward$X, 
           label = sprintf("%0.2f", round(df.empSI.forward$empmean, digits = 2)), 
           size = 3.5, hjust = 1.1) +
  labs(x="Mid point of time window", y="Serial interval (days)") +
  scale_x_continuous(limits = c(0, 18), breaks = seq(1, 17, 1), labels =
                       c("Jan 10",
                         "Jan 18", "Jan 19", "Jan 20", "Jan 21", "Jan 22",
                         "Jan 23", "Jan 24", "Jan 25", "Jan 26", "Jan 27", "Jan 28",
                         "Jan 29", "Jan 30", "Jan 31", "Feb 1",
                         "Feb 14")) +
  scale_x_cut(breaks = c(1.5, 16.5), space = 0.5, which = c(1, 3), scales = c(0.125, 0.125)) +
  geom_segment(aes(x = 7, y = 11, xend = 7, yend = 9), 
               arrow = arrow(length = unit(0.05, "inches")), lwd=0.5, color = "red") +
  geom_segment(aes(x = 7, y = 11, xend = 7.5, yend = 11), 
               lwd=0.5, color = "red") +
  annotate("text", x = 9 + 0.5, y = 11, label = "Lockdown in Wuhan", size = 3.5) + 
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
  ) + scale_y_continuous(breaks = seq(1, 12, 2)) 


dfIP.infector.back <- read.csv("dfIPinfector.backward_fitdistWB.csv")
dfIP.infectee.for <- read.csv("dfIPinfectee.forward_fitdistWB.csv")


dfIP.comb.refbyI <- rbind(dfIP.infector.back, dfIP.infectee.for)
dfIP.comb.refbyI$strtfy <- c(rep("Infector", 17), rep("Infectee", 17))

p.IP.comb.refbyI <- ggplot(dfIP.comb.refbyI, aes(x = X, y = mu, ymin = muLB, 
                                                 ymax = muUB, group = strtfy, color = strtfy)) + 
  geom_point(position = position_dodge(0.7), size = 2) + 
  geom_errorbar(width = 0, size = 0.5, position = position_dodge(0.7)) +
  geom_text(mapping = aes(y = muUB + 0.50, label = sprintf("%0.2f", round(mu, digits = 2))), 
            size = 3.5, position = position_dodge(1.1)) +
  labs(x="Mid point of time window", y="Incubation period (days)") +
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
    legend.title = element_blank(),
    legend.position = "none"
  ) + scale_y_continuous(breaks = seq(1, 11, 2), limits = c(0, 12)) 



pforward.SIIP.comb <- plot_grid(print(p.empSI.forward), print(p.IP.comb.refbyI), nrow = 2,
                             labels = c("a", "b"))


ggsave("Fig_forward_empSI_estIP.pdf", width = 12.5, height = 8, units = "in", plot = pforward.SIIP.comb)

# fig S1

dfempSI <- read.csv("backwardempSI.csv")

dfempSI$timewindow <- factor(dfempSI$timewindow, levels = unique(dfempSI$timewindow)[1:17])

tw.uniq <- unique(dfempSI$timewindow)

empmean <- numeric(17)
empq1 <- numeric(17)
empq3 <- numeric(17)
for(i in 1:17){
  tmp <- subset(dfempSI, timewindow == tw.uniq[i])
  empmean[i] <- mean(tmp$empSI)
  empq1[i] <- quantile(tmp$empSI, 0.25)
  empq3[i] <- quantile(tmp$empSI, 0.75)
  
}
empmean
empq1
empq3

df.empSI.backward <- data.frame(
  X = seq(1, 17),
  empmean = empmean,
  empq1 = empq1,
  empq3 = empq3
)

p.empSI.backward <- ggplot(df.empSI.backward, aes(x = X, y = empmean, ymin = empq1, ymax = empq3)) +
  geom_point(size = 1.0) + 
  geom_errorbar(width = 0, size = 0.5) +
  annotate("text", y = df.empSI.backward$empmean + 0.5, x = df.empSI.backward$X, 
           label = sprintf("%0.2f", round(df.empSI.backward$empmean, digits = 2)), 
           size = 3.5, hjust = 1.1) +
  labs(x="Mid point of time window", y="Serial interval (days)") +
  scale_x_continuous(limits = c(0, 18), breaks = seq(1, 17, 1), labels =
                       c("Jan 13", "Jan 24", "Jan 25", "Jan 26", "Jan 27", "Jan 28",
                         "Jan 29", "Jan 30", "Jan 31", "Feb 1", "Feb 2", "Feb 3",
                         "Feb 4", "Feb 5", "Feb 6", "Feb 7", "Feb 17")) + 
  scale_x_cut(breaks = c(1.5, 16.5), space = 0.5, which = c(1, 3), scales = c(0.125, 0.125)) +
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
  ) + scale_y_continuous(breaks = seq(1, 12, 2)) 


dfIP.infectee.back <- read.csv("dfIPinfectee.backward_fitdistWB.csv")
dfIP.infector.for <- read.csv("dfIPinfector.forward_fitdistG.csv")

dfIP.comb.refbyS <- rbind(dfIP.infectee.back, dfIP.infector.for)
dfIP.comb.refbyS$strtfy <- c(rep("Infectee", 17), rep("Infector", 17))

p.IP.comb.refbyS <- ggplot(dfIP.comb.refbyS, aes(x = X, y = mu, ymin = muLB, 
                                                 ymax = muUB, group = strtfy, color = strtfy)) + 
  geom_point(position = position_dodge(0.7), size = 2) + 
  geom_errorbar(width = 0, size = 0.5, position = position_dodge(0.7)) +
  geom_text(mapping = aes(y = muUB + 0.30, label = sprintf("%0.2f", round(mu, digits = 2))), 
            size = 3.5, position = position_dodge(0.9)) +
  labs(x="Mid point of time window", y="Incubation period (days)") +
  scale_x_continuous(limits = c(0, 18), breaks = seq(1, 17, 1), labels =
                       c("Jan 13", "Jan 24", "Jan 25", "Jan 26", "Jan 27", "Jan 28",
                         "Jan 29", "Jan 30", "Jan 31", "Feb 1", "Feb 2", "Feb 3",
                         "Feb 4", "Feb 5", "Feb 6", "Feb 7", "Feb 17")) +
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
    legend.title = element_blank(),
    legend.position = "none"
  ) + scale_y_continuous(breaks = seq(2, 12, 2), limits = c(0, 12.5)) 



pbackward.SIIP.comb <- plot_grid(print(p.empSI.backward), print(p.IP.comb.refbyS), nrow = 2,
                                labels = c("a", "b"))


ggsave("Fig_backward_empSI_estIP.pdf", width = 12.5, height = 8, units = "in", plot = pbackward.SIIP.comb)











