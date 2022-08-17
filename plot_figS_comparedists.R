
rm(list = ls())

library(ggplot2)
library(ggpubr)


IPinfectee <- read.csv("df.IPfor.infectee.summary.csv")
IPinfectee$X <- rep(seq(1, 17), 3)
IPinfector <- read.csv("df.IPback.infector.summary.csv")
IPinfector$X <- rep(seq(1, 17), 3)
GTmain <- read.csv("GTmain1_summary.csv")
GTmain$X <- rep(seq(1, 17), 3)


p.infectee.mu <- ggplot(IPinfectee, aes(x = X, y = mu, ymin = muLB, 
                             ymax = muUB, group = dist, color = dist)) + 
  geom_point(position = position_dodge(0.7), size = 1.5) + 
  geom_errorbar(width = 0, size = 0.5, position = position_dodge(0.7)) +
  labs(x="Mid point of time window", y=" ") +
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
    axis.text = element_text(size=8),
    axis.line = element_line(color="black", size=0.75), 
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  ) + scale_y_continuous(breaks = seq(1, 11, 2), limits = c(0, 11.5)) +
  ggtitle("Mean of forward infectee's IP")


p.infector.mu <- ggplot(IPinfector, aes(x = X, y = mu, ymin = muLB, 
                                        ymax = muUB, group = dist, color = dist)) + 
  geom_point(position = position_dodge(0.7), size = 1.5) + 
  geom_errorbar(width = 0, size = 0.5, position = position_dodge(0.7)) +
  labs(x="Mid point of time window", y=" ") +
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
    axis.text = element_text(size=8),
    axis.line = element_line(color="black", size=0.75), 
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  ) + scale_y_continuous(breaks = seq(1, 11, 2), limits = c(0, 11.5)) +
  ggtitle("Mean of backward infectee's IP")


p.GT.mu <- ggplot(GTmain, aes(x = X, y = mu, ymin = muLB, 
                                        ymax = muUB, group = dist, color = dist)) + 
  geom_point(position = position_dodge(0.7), size = 1.5) + 
  geom_errorbar(width = 0, size = 0.5, position = position_dodge(0.7)) +
  labs(x="Mid point of time window", y=" ") +
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
    axis.text = element_text(size=8),
    axis.line = element_line(color="black", size=0.75), 
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  ) + scale_y_continuous(breaks = seq(1, 11, 2), limits = c(0, 11.5)) +
  ggtitle("Mean of forward GT\n(GT samples based on IPs fitted by Weibull)")


p.infectee.sd <- ggplot(IPinfectee, aes(x = X, y = sd, ymin = sdLB, 
                                        ymax = sdUB, group = dist, color = dist)) + 
  geom_point(position = position_dodge(0.7), size = 1.5) + 
  geom_errorbar(width = 0, size = 0.5, position = position_dodge(0.7)) +
  labs(x="Mid point of time window", y=" ") +
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
    axis.text = element_text(size=8),
    axis.line = element_line(color="black", size=0.75), 
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  ) + scale_y_continuous(breaks = seq(1, 7, 1), limits = c(0, 7)) +
  ggtitle("SD of forward infectee's IP")

p.infector.sd <- ggplot(IPinfector, aes(x = X, y = sd, ymin = sdLB, 
                                        ymax = sdUB, group = dist, color = dist)) + 
  geom_point(position = position_dodge(0.7), size = 1.5) + 
  geom_errorbar(width = 0, size = 0.5, position = position_dodge(0.7)) +
  labs(x="Mid point of time window", y=" ") +
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
    axis.text = element_text(size=8),
    axis.line = element_line(color="black", size=0.75), 
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  ) + scale_y_continuous(breaks = seq(1, 7, 1), limits = c(0, 7)) +
  ggtitle("SD of backward infector's IP")


p.GT.sd <- ggplot(GTmain, aes(x = X, y = sd, ymin = sdLB, 
                                        ymax = sdUB, group = dist, color = dist)) + 
  geom_point(position = position_dodge(0.7), size = 1.5) + 
  geom_errorbar(width = 0, size = 0.5, position = position_dodge(0.7)) +
  labs(x="Mid point of time window", y=" ") +
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
    axis.text = element_text(size=8),
    axis.line = element_line(color="black", size=0.75), 
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  ) + scale_y_continuous(breaks = seq(1, 7, 1), limits = c(0, 7)) +
  ggtitle("SD of forward GT\n(GT samples based on IPs fitted by Weibull)")


p.S3 <- plot_grid(print(p.infectee.mu), print(p.infectee.sd), print(p.infector.mu), 
                  print(p.infector.sd), print(p.GT.mu), print(p.GT.sd),
                  nrow = 3, ncol = 2, labels = c("a", "d", "b", "e", "c", "f"))

ggsave("Fig_comparedists.pdf", width = 16, height = 8, units = "in", plot = p.S3)












