library(ggplot2)
library(viridis)
library(patchwork)

setwd("/Users/srgib/Desktop/Univeristy of Arizona/Gutenkunst Group/Mutation_Selection_Balance_Autos_Allos/python_simulations/bistability/")

auto_dom <- read.csv("dominant/fuller_bifurcation.txt", header = T)
auto_dom <- auto_dom[auto_dom$init_q %% 10 == 0,]

auto_add <- read.csv("additive/fuller_bifurcation.txt", header = T)
auto_rec <- read.csv("recessive/fuller_bifurcation.txt", header = T)

auto_dom$Dominance <- rep("Dominant")
auto_add$Dominance <- rep("Additive")
auto_rec$Dominance <- rep("Recessive")

auto <- rbind(auto_dom, auto_add, auto_rec)
auto$Dominance <- factor(auto$Dominance, levels = c("Dominant", 'Additive', 'Recessive'))

auto$q <- auto$G4 + 0.75*auto$G3 + 0.5*auto$G2 + 0.25*auto$G1

auto$w <- auto$G4*(1-auto$s) + auto$G3*(1-auto$s*auto$h3) + auto$G2*(1-auto$s*auto$h2) + auto$G1*(1-auto$s*auto$h1) + auto$G0

auto$w_var <- auto$G4*(auto$w - (1-auto$s))^2 + auto$G3*(auto$w -(1-auto$s*auto$h3))^2 + auto$G2*(auto$w - (1-auto$s*auto$h2))^2 + auto$G1*(auto$w -(1-auto$s*auto$h1))^2 + auto$G0*(auto$w - 1)^2

p <- ggplot(auto[auto$Generation < 2000,], aes(x=Generation, y=q)) + 
  #geom_hline(yintercept = 0.3608, linetype = 'dotted') + 
  geom_line(aes(color=as.factor(init_q))) + 
  scale_color_viridis(discrete=T) + theme_bw() + ggtitle("Autotetraploid") + 
  theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank()) + ylab("") + facet_grid(Dominance~.)


dip_dom <- read.csv("dominant/dip_fuller_bifurcation.txt", header = T)
dip_dom <- dip_dom[dip_dom$init_q %% 10 == 0,]

dip_add <- read.csv("additive/dip_fuller_bifurcation.txt", header = T)
dip_rec <- read.csv("recessive/dip_fuller_bifurcation.txt", header = T)

dip_dom$Dominance <- rep("Dominant")
dip_add$Dominance <- rep("Additive")
dip_rec$Dominance <- rep("Recessive")

dip <- rbind(dip_dom, dip_add, dip_rec)

dip$q <- dip$G2 + 0.5*dip$G1
dip$Dominance <- factor(dip$Dominance, levels = c("Dominant", 'Additive', 'Recessive'))

dip$w <- dip$G2*(1-dip$s) + dip$G1*(1-dip$s*dip$h) + dip$G0

dip$w_var <- dip$G2*(dip$w - (1-dip$s))^2 + dip$G1*(dip$w -(1-dip$s*dip$h))^2 + dip$G0*(dip$w - 1)^2


q <- ggplot(dip[dip$Generation < 2000,], aes(x=Generation, y=q)) + 
  #geom_hline(yintercept = 0.9201, linetype = 'dotted') + 
  geom_line(aes(color=as.factor(init_q))) + ggtitle("Diploid") + 
  scale_color_viridis(discrete=T) + theme_bw() + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank()) + ylab("Frequency of Selected Allele") + 
  facet_grid(Dominance~.)

png("bistability_allele_trajectories.png", width = 6, height = 6, res = 300, units = "in")
q | p
dev.off()

p1 <- ggplot(auto[auto$Generation < 2000,], aes(x=Generation, y=w)) + 
  #geom_hline(yintercept = 0.3608, linetype = 'dotted') + 
  geom_line(aes(color=as.factor(init_q))) + 
  scale_color_viridis(discrete=T) + theme_bw() + ggtitle("Autotetraploid") + 
  theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank()) + ylab("") + facet_grid(Dominance~.)

q1 <- ggplot(dip[dip$Generation < 2000,], aes(x=Generation, y=w)) + 
  #geom_hline(yintercept = 0.9201, linetype = 'dotted') + 
  geom_line(aes(color=as.factor(init_q))) + ggtitle("Diploid") + 
  scale_color_viridis(discrete=T) + theme_bw() + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank()) + ylab("Population Fitness") + 
  facet_grid(Dominance~.)

png("bistability_fitness_trajectories.png", width = 6, height = 6, res = 300, units = "in")
q1 | p1
dev.off()


p2 <- ggplot(auto[auto$Generation < 2000,], aes(x=Generation, y=1-q)) + 
  #geom_hline(yintercept = 0.3608, linetype = 'dotted') + 
  geom_line(aes(color=as.factor(init_q))) + 
  scale_color_viridis(discrete=T) + theme_bw() + ggtitle("Autotetraploid") + 
  theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank()) + ylab("") + facet_grid(Dominance~.)

q2 <- ggplot(dip[dip$Generation < 2000,], aes(x=Generation, y=1-q)) + 
  #geom_hline(yintercept = 0.9201, linetype = 'dotted') + 
  geom_line(aes(color=as.factor(init_q))) + ggtitle("Diploid") + 
  scale_color_viridis(discrete=T) + theme_bw() + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank()) + ylab("Frequency of Neutral Allele") + 
  facet_grid(Dominance~.)

png("bistability_neutral_allele_trajectories.png", width = 6, height = 6, res = 300, units = "in")
q2 | p2
dev.off()

p3 <- ggplot(auto[auto$Generation < 2000,], aes(x=Generation, y=w_var)) + 
  #geom_hline(yintercept = 0.3608, linetype = 'dotted') + 
  geom_line(alpha = 0.5, aes(color=as.factor(init_q))) + 
  scale_color_viridis(discrete=T) + theme_bw() + ggtitle("Autotetraploid") + 
  theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank()) + ylab("") + facet_grid(Dominance~.)

q3 <- ggplot(dip[dip$Generation < 2000,], aes(x=Generation, y=w_var)) + 
  #geom_hline(yintercept = 0.9201, linetype = 'dotted') + 
  geom_line(alpha = 0.5, aes(color=as.factor(init_q))) + ggtitle("Diploid") + 
  scale_color_viridis(discrete=T) + theme_bw() + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank()) + ylab("Fitness Variance") + 
  facet_grid(Dominance~.)

png("bistability_fit_var_trajectories.png", width = 6, height = 6, res = 300, units = "in")
q3 | p3
dev.off()
