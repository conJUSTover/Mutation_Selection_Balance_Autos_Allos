library(ggplot2)
library(viridis)
library(patchwork)

auto <- read.csv("fuller_bifurcation.txt", header = T)

auto$q <- auto$G4 + 0.75*auto$G3 + 0.5*auto$G2 + 0.25*auto$G1

p <- ggplot(auto[auto$Generation < 500,], aes(x=Generation, y=q)) + 
  geom_hline(yintercept = 0.3608, linetype = 'dotted') + 
  geom_line(aes(color=as.factor(init_q))) + 
  scale_color_viridis(discrete=T) + theme_bw() + ggtitle("Autotetraploid") + 
  theme(legend.position = "none") + ylab("")


dip <- read.csv("dip_fuller_bifurcation.txt", header = T)

dip$q <- dip$G2 + 0.5*dip$G1

q <- ggplot(dip[dip$Generation < 500,], aes(x=Generation, y=q)) + 
  geom_hline(yintercept = 0.9201, linetype = 'dotted') + 
  geom_line(aes(color=as.factor(init_q))) + ggtitle("Diploid") + 
  scale_color_viridis(discrete=T) + theme_bw() + 
  theme(legend.position = "none") + ylab("Frequency of Selected Allele")

png("Bistability_Trajectories.png", width = 6, height = 3, res = 300, units = "in")
q | p
dev.off()



