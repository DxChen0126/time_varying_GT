
library(latex2exp)
library(ggplot2)



arrow <- data.frame(
  x=0,
  y=0,
  x1=4,
  y1=0
)


pGT.diag.forward <- ggplot(arrow) +
  geom_segment(aes(0.5, 0, xend=0.5, yend=0.2), lty=2) +
  geom_segment(aes(1.5, 0, xend=1.5, yend=0.2), lty=2) +
  geom_segment(aes(2.5, 0, xend=2.5, yend=0.2), lty=2) +
  geom_segment(aes(3.5, 0, xend=3.5, yend=0.15), lty=2) +
  geom_segment(aes(0.25, 0, xend = 0.75, yend = 0), size=3, lty = 1, color = "#008080") +
  geom_point(aes(1.5, 0), size=3, col = "red") +
  geom_segment(aes(2.25, 0, xend = 2.75, yend = 0), size=3, lty = 1, color = "#008080") +
  geom_point(aes(3.5, 0), size=3) +
  geom_text(x=4.0, y=-0.025, label="Time", fontface = "bold") +
  annotate("text", x=1.47, y=-0.05, label="Infector onset") +
  annotate("text", x=0.5, y=-0.05, label="Infector's exposure window") +
  annotate("text", x=1.5, y=-0.025, label="Start of follow up", col="red") +
  geom_segment(aes(1.5, -0.024, xend=1.5, yend=-0.005), lty=1, arrow = arrow(length = unit(0.05, "inches")), lwd=0.5) +
  annotate("text", x=2.5, y=-0.05, label="Infectee's exposure window") +
  annotate("text", x=3.5, y=-0.05, label="Infectee onset") +
  geom_segment(aes(x, y, xend=x1, yend=y1), arrow = arrow(length = unit(0.2, "inches")), lwd=0.5) +
  geom_segment(aes(1.5, 0.1, xend=0.5, yend=0.1), arrow = arrow(length = unit(0.1, "inches")), lwd=1) +
  geom_segment(aes(0.5, 0.2, xend=2.5, yend=0.2), arrow = arrow(length = unit(0.1, "inches")), lwd=1.5, col="red") +
  geom_segment(aes(2.5, 0.05, xend=3.5, yend=0.05), arrow = arrow(length = unit(0.1, "inches")), lwd=1) +
  geom_segment(aes(1.5, 0.15, xend=3.5, yend=0.15), arrow = arrow(length = unit(0.1, "inches")), lwd=1) +
  geom_segment(aes(1.5, 0.05, xend=2.5, yend=0.05), lty=2) +
  annotate("text", x=2.0, y=0.165, label="Forward serial interval") +
  annotate("text", x=2.75, y=0.165, label=TeX("$SI^{forward}$")) +
  annotate("text", x=0.9, y=0.115, label="Backward incubation period") +
  annotate("text", x=1.5, y=0.115, label=TeX("$IP_{infector}^{backward}$")) +
  annotate("text", x=0.9, y=0.215, label="Forward generation time") +
  annotate("text", x=2.0, y=0.215, label=TeX("$GT^{forward} = SI^{forward} - $IP_{infectee}^{forward} + $IP_{infector}^{backward}$$")) +
  annotate("text", x=2.8, y=0.08, label="Forward incubation period") +
  annotate("text", x=2.8, y=0.065, label="(Infector onset referenced)") +
  annotate("text", x=3.4, y=0.07, label=TeX("$IP_{infectee}^{forward}$")) +
  annotate("text", x=2, y=0.025, label=TeX("$SI^{forward} - $IP_{infectee}^{forward}$")) +
#  ggtitle("Composition of forward generation time") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.background=element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

pGT.diag.forward


pGT.diag.backward <- ggplot(arrow) +
  geom_segment(aes(0.5, 0, xend=0.5, yend=0.2), lty=2) +
  geom_segment(aes(1.5, 0, xend=1.5, yend=0.2), lty=2) +
  geom_segment(aes(2.5, 0, xend=2.5, yend=0.2), lty=2) +
  geom_segment(aes(3.5, 0, xend=3.5, yend=0.15), lty=2) +
  geom_segment(aes(0.25, 0, xend = 0.75, yend = 0), size=3, lty = 1, color = "#008080") +
  geom_point(aes(1.5, 0), size=3) +
  geom_segment(aes(2.25, 0, xend = 2.75, yend = 0), size=3, lty = 1, color = "#008080") +
  geom_point(aes(3.5, 0), size=3, col = "red") +
  geom_text(x=4.0, y=-0.025, label="Time", fontface = "bold") +
  annotate("text", x=1.5, y=-0.05, label="Infector onset") +
  annotate("text", x=0.5, y=-0.05, label="Infector's exposure window") +
  annotate("text", x=3.5, y=-0.025, label="End of follow up", col="red") +
  geom_segment(aes(3.5, -0.024, xend=3.5, yend=-0.005), lty=1, arrow = arrow(length = unit(0.05, "inches")), lwd=0.5) +
  annotate("text", x=2.5, y=-0.05, label="Infectee's exposure window") +
  annotate("text", x=3.47, y=-0.05, label="Infectee onset") +
  geom_segment(aes(x, y, xend=x1, yend=y1), arrow = arrow(length = unit(0.2, "inches")), lwd=0.5) +
  geom_segment(aes(0.5, 0.075, xend=1.5, yend=0.075), arrow = arrow(length = unit(0.1, "inches")), lwd=1) +
  geom_segment(aes(2.5, 0.2, xend=0.5, yend=0.2), arrow = arrow(length = unit(0.1, "inches")), lwd=1.5, col="red") +
  geom_segment(aes(3.5, 0.05, xend=2.5, yend=0.05), arrow = arrow(length = unit(0.1, "inches")), lwd=1) +
  geom_segment(aes(3.5, 0.15, xend=1.5, yend=0.15), arrow = arrow(length = unit(0.1, "inches")), lwd=1) +
  geom_segment(aes(2.5, 0.05, xend=1.5, yend=0.05), lty=2) +
  annotate("text", x=2.0, y=0.165, label="Backward serial interval") +
  annotate("text", x=2.75, y=0.165, label=TeX("$SI^{backward}$")) +
  annotate("text", x=0.9, y=0.11, label="Incubation period") +
  annotate("text", x=0.9, y=0.09, label="(Infectee onset referenced)") +
#  annotate("text", x=1, y=0.095, label="approximately forward)") +
  annotate("text", x=1.6, y=0.1, label=TeX("$IP_{infector}^{infectee $ $ referenced}$")) +
  annotate("text", x=0.9, y=0.215, label="Backward generation time") +
  annotate("text", x=2.25, y=0.215, label=TeX("$GT^{backward} = SI^{backward} - $IP_{infectee}^{backward} + $IP_{infector}^{infectee $ $ referenced}$$")) +
  annotate("text", x=2.8, y=0.065, label="Backward incubation period") +
  annotate("text", x=3.4, y=0.065, label=TeX("$IP_{infectee}^{backward}$")) +
  annotate("text", x=2, y=0.025, label=TeX("$SI^{backward} - $IP_{infectee}^{backward}$")) +
#  ggtitle("Composition of backward generation time") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.background=element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

pGT.diag.backward

cols <- c("Infector" = "#48D1CC","Infectee" = "#FF8C00")
point_shapes <- c("Infector" = 21,"Infectee" = 24)




p.SI.diag.leftcens <- ggplot(arrow) +
  geom_point(aes(0.5, 0.25), shape = 21, fill = cols[1],size = 2) +
  geom_segment(aes(0.5, 0.25, xend = 1.2, yend = 0.25),lty = 1) +
  geom_point(aes(1.2, 0.25), shape = 24, fill = cols[2],) +
  geom_point(aes(0.75, 0.45), shape = 21, fill = cols[1], size = 2) +
  geom_segment(aes(0.75, 0.45, xend = 1.5, yend = 0.45),lty = 1) +
  geom_point(aes(1.5, 0.45), shape = 24, fill = cols[2]) +
  geom_point(aes(0.75, 0.35), shape = 21, fill = cols[1], size = 2) +
  geom_segment(aes(0.75, 0.35, xend = 1.75, yend = 0.35),lty = 1) +
  geom_point(aes(1.75, 0.35), shape = 24, fill = cols[2]) +
  geom_point(aes(1.0, 0.55), shape = 21, fill = cols[1], size = 2) +
  geom_segment(aes(1.0, 0.55, xend = 2.25, yend = 0.55),lty = 1) +
  geom_point(aes(2.25, 0.55), shape = 24, fill = cols[2]) +
  geom_point(aes(1.0, 0.65), shape = 21, fill = cols[1], size = 2) +
  geom_segment(aes(1.0, 0.65, xend = 1.85, yend = 0.65),lty = 1) +
  geom_point(aes(1.85, 0.65), shape = 24, fill = cols[2]) +
  geom_point(aes(1.0, 0.65), shape = 21, fill = cols[1], size = 2) +
  geom_segment(aes(1.0, 0.65, xend = 1.85, yend = 0.65),lty = 1) +
  geom_point(aes(1.85, 0.65), shape = 24, fill = cols[2]) +
  geom_point(aes(1.0, 0.75), shape = 21, fill = cols[1], size = 2) +
  geom_segment(aes(1.0, 0.75, xend = 2.0, yend = 0.75),lty = 1) +
  geom_point(aes(2.0, 0.75), shape = 24, fill = cols[2]) +
  geom_point(aes(1.0, 0.85), shape = 21, fill = cols[1], size = 2) +
  geom_segment(aes(1.0, 0.85, xend = 2.5, yend = 0.85),lty = 1) +
  geom_point(aes(2.5, 0.85), shape = 24, fill = cols[2]) +
  geom_point(aes(1.0, 0.95), shape = 21, fill = cols[1], size = 2) +
  geom_segment(aes(1.0, 0.95, xend = 2.25, yend = 0.95),lty = 1) +
  geom_point(aes(2.25, 0.95), shape = 24, fill = cols[2]) +
  geom_segment(aes(1.0, 0.15, xend = 1.0, yend = 1.0),lty = 2) +
  geom_segment(aes(0.25, 0.15, xend = 2.75, yend = 0.15), arrow = arrow(length = unit(0.2, "inches")), lwd=0.5) +
  geom_segment(aes(1, 0.11, xend = 1, yend = 0.14), arrow = arrow(length = unit(0.05, "inches")), lwd=0.5) +
  annotate("text", x = 1, y= 0.1, label="Start of follow up", color = "red") + 
  geom_rect(aes(xmin = 1, xmax = 2.5, ymin = 0.5, ymax = 1), alpha = 0.1, fill = "green") + 
  annotate("text", x = 1.5, y= 1, label="Observed pairs") + 
  geom_rect(aes(xmin = 0.25, xmax = 1, ymin = 0.15, ymax = 0.5), alpha = 0.1, fill = "red") + 
  annotate("text", x = 0.5, y= 0.5, label="Unobserved pairs") + 
  geom_text(x=2.75, y=0.1, label="Time", fontface = "bold") +
  scale_colour_manual(name="Onset", values=cols) + 
  geom_point(aes(0.1, 0.8), shape = 21, fill = cols[1],size = 2) +
  annotate("text", x = 0.5, y= 0.8, label="Infector's symptom onset")+
  geom_point(aes(0.1, 0.75), shape = 24, fill = cols[2],size = 2) +
  annotate("text", x = 0.5, y= 0.75, label="Infectee's symptom onset")+
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.background=element_blank(),
    plot.title = element_text(hjust = 0.5)
  )  
#  ggtitle("Illustration of left censoring in serial interval")
  
  
p.SI.diag.leftcens


p.SI.diag.rightcens <- ggplot(arrow) +
  geom_point(aes(0.25, 0.25), shape = 21, fill = cols[1],size = 2) +
  geom_segment(aes(0.25, 0.25, xend = 2, yend = 0.25),lty = 1) +
  geom_point(aes(2, 0.25), shape = 24, fill = cols[2],) +
  geom_point(aes(0.5, 0.35), shape = 21, fill = cols[1], size = 2) +
  geom_segment(aes(0.5, 0.35, xend = 2, yend = 0.35),lty = 1) +
  geom_point(aes(2, 0.35), shape = 24, fill = cols[2]) +
  geom_point(aes(0.5, 0.95), shape = 21, fill = cols[1], size = 2) +
  geom_segment(aes(0.5, 0.95, xend = 2.1, yend = 0.95),lty = 1) +
  geom_point(aes(2.1, 0.95), shape = 24, fill = cols[2]) +
  geom_point(aes(0.75, 0.85), shape = 21, fill = cols[1], size = 2) +
  geom_segment(aes(0.75, 0.85, xend = 2.25, yend = 0.85),lty = 1) +
  geom_point(aes(2.25, 0.85), shape = 24, fill = cols[2]) +
  geom_point(aes(0.75, 0.45), shape = 21, fill = cols[1], size = 2) +
  geom_segment(aes(0.75, 0.45, xend = 2, yend = 0.45),lty = 1) +
  geom_point(aes(2, 0.45), shape = 24, fill = cols[2]) +
  
  geom_point(aes(1.0, 0.75), shape = 21, fill = cols[1], size = 2) +
  geom_segment(aes(1.0, 0.75, xend = 2.5, yend = 0.75),lty = 1) +
  geom_point(aes(2.5, 0.75), shape = 24, fill = cols[2]) +
  
  geom_point(aes(1.25, 0.55), shape = 21, fill = cols[1], size = 2) +
  geom_segment(aes(1.25, 0.55, xend = 2.0, yend = 0.55),lty = 1) +
  geom_point(aes(2.0, 0.55), shape = 24, fill = cols[2]) +
  
  geom_point(aes(1.5, 0.65), shape = 21, fill = cols[1], size = 2) +
  geom_segment(aes(1.5, 0.65, xend = 2, yend = 0.65),lty = 1) +
  geom_point(aes(2, 0.65), shape = 24, fill = cols[2]) +
  
  geom_segment(aes(2.0, 0.15, xend = 2.0, yend = 1.0),lty = 2) +
  geom_segment(aes(0.25, 0.15, xend = 2.75, yend = 0.15), arrow = arrow(length = unit(0.2, "inches")), lwd=0.5) +
  geom_segment(aes(2, 0.11, xend = 2, yend = 0.14), arrow = arrow(length = unit(0.05, "inches")), lwd=0.5) +
  annotate("text", x = 2, y= 0.1, label="End of follow up", color = "red") + 
  geom_rect(aes(xmin = 0.25, xmax = 2, ymin = 0.15, ymax = 0.7), alpha = 0.1, fill = "green") + 
  annotate("text", x = 0.5, y= 0.7, label="Observed pairs") + 
  geom_rect(aes(xmin = 2, xmax = 2.75, ymin = 0.7, ymax = 1), alpha = 0.1, fill = "red") + 
  annotate("text", x = 2.25, y= 1, label="Unobserved pairs") + 
  geom_text(x=2.75, y=0.1, label="Time", fontface = "bold") +
  geom_point(aes(2.15, 0.5), shape = 21, fill = cols[1],size = 2) +
  annotate("text", x = 2.5, y= 0.5, label="Infector's symptom onset")+
  geom_point(aes(2.15, 0.45), shape = 24, fill = cols[2],size = 2) +
  annotate("text", x = 2.5, y= 0.45, label="Infectee's symptom onset")+
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.background=element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) 
#  ggtitle("Illustration of right censoring in serial interval")


p.SI.diag.rightcens


library(ggpubr)


p1 <- ggarrange(
  p.SI.diag.leftcens, NULL,  pGT.diag.forward, 
  nrow = 1, widths = c(1, 0.05, 1.25),
  labels = c("a", "", "c")
)

p2 <- ggarrange(
  p.SI.diag.rightcens, NULL, pGT.diag.backward, 
  nrow = 1, widths = c(1, 0.05, 1.25),
  labels = c("b", "", "d")
)




p.diag_combine <- ggarrange(p1, p2, nrow = 2, heights = c(1, 1))
  
  
# ggsave("Fig1_20220204.pdf", width = 18, height = 12, units = "in", 
#        plot = p.diag_combine)

ggsave("Fig_diagram_20220506.pdf", width = 18, height = 10, units = "in", 
               plot = p.diag_combine)











































