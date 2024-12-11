library(ggplot2)
library(viridis)
library(patchwork)

setwd("~/GitHub/Mutation_Selection_Balance_Autos_Allos/python_simulations/bistability/")

auto_dom <- read.csv("fuller_bifurcation.txt", header = T)
auto_dom <- auto_dom[auto_dom$init_q %% 10 == 0,]

auto_add <- read.csv("additive/fuller_bifurcation.txt", header = T)
auto_rec <- read.csv("recessive/fuller_bifurcation.txt", header = T)

auto_dom$Dominance <- rep("Dominant")
auto_add$Dominance <- rep("Additive")
auto_rec$Dominance <- rep("Recessive")

auto <- rbind(auto_dom, auto_add, auto_rec)
auto$Dominance <- factor(auto$Dominance, levels = c("Dominant", 'Additive', 'Recessive'))


auto$q <- auto$G4 + 0.75*auto$G3 + 0.5*auto$G2 + 0.25*auto$G1

p <- ggplot(auto[auto$Generation < 2000,], aes(x=Generation, y=q)) + 
  #geom_hline(yintercept = 0.3608, linetype = 'dotted') + 
  geom_line(aes(color=as.factor(init_q))) + 
  scale_color_viridis(discrete=T) + theme_bw() + ggtitle("Autotetraploid") + 
  theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank()) + ylab("") + facet_grid(Dominance~.)


dip_dom <- read.csv("dip_fuller_bifurcation.txt", header = T)
dip_dom <- dip_dom[dip_dom$init_q %% 10 == 0,]

dip_add <- read.csv("additive/dip_fuller_bifurcation.txt", header = T)
dip_rec <- read.csv("recessive/dip_fuller_bifurcation.txt", header = T)

dip_dom$Dominance <- rep("Dominant")
dip_add$Dominance <- rep("Additive")
dip_rec$Dominance <- rep("Recessive")

dip <- rbind(dip_dom, dip_add, dip_rec)

dip$q <- dip$G2 + 0.5*dip$G1
dip$Dominance <- factor(dip$Dominance, levels = c("Dominant", 'Additive', 'Recessive'))

q <- ggplot(dip[dip$Generation < 2000,], aes(x=Generation, y=q)) + 
  #geom_hline(yintercept = 0.9201, linetype = 'dotted') + 
  geom_line(aes(color=as.factor(init_q))) + ggtitle("Diploid") + 
  scale_color_viridis(discrete=T) + theme_bw() + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank()) + ylab("Frequency of Selected Allele") + 
  facet_grid(Dominance~.)

png("Bistability_Trajectories.png", width = 6, height = 6, res = 300, units = "in")
q | p
dev.off()



