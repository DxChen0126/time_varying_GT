library(ggplot2)
library(ggbreak)
library(cowplot)

df.estGT.g <- read.csv("backwardGT_est_gamma.csv")
df.estGT.g$X <- seq(1, 17, 1)


p.GT.backward.main.bootCI <- ggplot(df.estGT.g, aes(x = X, y = mu, ymin = muLB.boot, ymax = muUB.boot)) +
  geom_point(size = 1.0) + 
  geom_errorbar(width = 0, size = 0.5) +
  geom_text(mapping = aes(y = muUB.boot+0.30, label = sprintf("%0.2f", round(mu, digits = 2))), size = 3.5)+
  labs(x="Mid point of time window", y="Mean generation time (days)") +
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
  ) + 
 scale_y_continuous(breaks = seq(4, 7, 1), limits = c(3.5, 7.5)) 

p.GT.backward.sd.bootCI <- ggplot(df.estGT.g, aes(x = X, y = sd, ymin = sdLB.boot, ymax = sdUB.boot)) +
  geom_point(size = 1.0) + 
  scale_x_continuous(limits = c(0, 18), breaks = seq(1, 17, 1), labels =
                       c("Jan 13", "Jan 24", "Jan 25", "Jan 26", "Jan 27", "Jan 28",
                         "Jan 29", "Jan 30", "Jan 31", "Feb 1", "Feb 2", "Feb 3",
                         "Feb 4", "Feb 5", "Feb 6", "Feb 7", "Feb 17") ) + 
  geom_errorbar(width = 0, size = 0.5) +
  geom_text(mapping = aes(y = sdUB.boot+0.30, label = sprintf("%0.2f", round(sd, digits = 2))), size = 3.5) +
  #  geom_text(mapping = aes(y = 4, label = n_sample), size = 2.5) +
  labs(x="Mid point of time window", y="Standard deviation of generation time (days)") +
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
    axis.line = element_line(size=0.75),
    panel.grid = element_blank()
  ) + scale_y_continuous(breaks = seq(2, 5, 1), limits = c(1, 5.5)) 


pGT.back.bootCI <- plot_grid(print(p.GT.backward.main.bootCI), print(p.GT.backward.sd.bootCI), nrow = 2,
                             labels = c("a", "b"))


ggsave("Fig_backwardGT.pdf", width = 12.5, height = 8, units = "in", plot = pGT.back.bootCI)


