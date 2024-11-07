library(ggplot2)
library(viridis)

setwd("~/GitHub/Mutation_Selection_Balance_Autos_Allos/python_simulations/bistability/")

auto <- read.csv("fuller_bifurcation.txt", header = T)

auto$q <- auto$G4 + 0.75*auto$G3 + 0.5*auto$G2 + 0.25*auto$G1

ggplot(auto[auto$Generation < 2000,], aes(x=Generation, y=q)) + 
  geom_line(aes(color=as.factor(init_p))) + 
  geom_hline(yintercept = 0.3608, linetype = 'dotted') + 
  scale_color_viridis(discrete=T) + theme_bw() + 
  theme(legend.position = "none")


dip <- read.csv("dip_fuller_bifurcation.txt", header = T)

dip$q <- dip$G4 + 0.75*dip$G3 + 0.5*dip$G2 + 0.25*dip$G1

ggplot(dip[dip$Generation < 2000,], aes(x=Generation, y=q)) + 
  geom_line(aes(color=as.factor(init_q))) + 
  geom_hline(yintercept = 0.9201, linetype = 'dotted') + 
  scale_color_viridis(discrete=T) + theme_bw() + 
  theme(legend.position = "none")

