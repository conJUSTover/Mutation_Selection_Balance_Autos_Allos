library(ggplot2)
library(viridis)
library(patchwork)

setwd("~/GitHub/Mutation_Selection_Balance_Autos_Allos/python_simulations/weak_mutation")

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
  geom_line(aes(color=as.factor(init_q))) + 
  scale_color_viridis(discrete=T) + theme_bw() + ggtitle("Autotetraploid") + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
  ylab("") + facet_grid(Dominance~.)


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
  geom_line(aes(color=as.factor(init_q))) + ggtitle("Diploid") + 
  scale_color_viridis(discrete=T) + theme_bw() + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank()) + ylab("Frequency of Derived Allele") + 
  facet_grid(Dominance~.) + xlab("")

allo_dom <- read.csv("dominant/allo_fuller_bifurcation.txt", header = T)

allo_add <- read.csv("additive/allo_fuller_bifurcation.txt", header = T)
allo_rec <- read.csv("recessive/allo_fuller_bifurcation.txt", header = T)

allo_dom$Dominance <- rep("Dominant")
allo_add$Dominance <- rep("Additive")
allo_rec$Dominance <- rep("Recessive")

allo <- rbind(allo_dom, allo_add, allo_rec)

allo$q <- allo$g11 + allo$g01 # this is fine, because the populations start in HWE
allo$Dominance <- factor(allo$Dominance, levels = c("Dominant", 'Additive', 'Recessive'))

allo$w <- allo$G22*(1-allo$s) + (allo$G21+allo$G12)*(1-allo$s*allo$h3) + (allo$G20 + allo$G02 + allo$G11)*(1-allo$s*allo$h2) + (allo$G10+allo$G01)*(1-allo$s*allo$h1) + allo$G00

allo$w_var <- allo$G22*(allo$w - (1-allo$s))^2 + (allo$G21+allo$G12)*(allo$w -(1-allo$s*allo$h3))^2 + (allo$G20 + allo$G02 + allo$G11)*(allo$w - (1-allo$s*allo$h2))^2 + (allo$G10 + allo$G01)*(allo$w -(1-allo$s*allo$h1))^2 + allo$G00*(allo$w - 1)^2

r <- ggplot(allo[allo$Generation < 2000,], aes(x=Generation, y=q)) + 
  geom_line(aes(color=as.factor(init_q))) + ggtitle("Allotetraploid") + 
  scale_color_viridis(discrete=T) + theme_bw() + 
  theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank()) + ylab("") + 
  facet_grid(Dominance~.) + xlab("")

allele_traj = q | p | r
ggsave("weak_mut_allele_trajectories.pdf", allele_traj, width = 6.5, height = 5, units = "in")


p1 <- ggplot(auto[auto$Generation < 2000,], aes(x=Generation, y=w)) + 
  geom_line(aes(color=as.factor(init_q))) + 
  scale_color_viridis(discrete=T) + theme_bw() + ggtitle("Autotetraploid") + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
  ylab("") + facet_grid(Dominance~.)

q1 <- ggplot(dip[dip$Generation < 2000,], aes(x=Generation, y=w)) + 
  geom_line(aes(color=as.factor(init_q))) + ggtitle("Diploid") + 
  scale_color_viridis(discrete=T) + theme_bw() + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank()) + ylab("Population Fitness") + 
  facet_grid(Dominance~.) + xlab("")

r1 <- ggplot(allo[allo$Generation < 2000,], aes(x=Generation, y=w)) + 
  geom_line(aes(color=as.factor(init_q))) + ggtitle("Allotetraploid") + 
  scale_color_viridis(discrete=T) + theme_bw() + 
  theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank()) + ylab("") + 
  facet_grid(Dominance~.) + xlab("")


fitness_traj = q1 | p1 | r1
ggsave("weak_mut_fitness_trajectories.pdf", fitness_traj, width = 6.5, height = 5, units = "in")


p2 <- ggplot(auto[auto$Generation < 2000,], aes(x=Generation, y=1-q)) + 
  geom_line(aes(color=as.factor(init_q))) + 
  scale_color_viridis(discrete=T) + theme_bw() + ggtitle("Autotetraploid") + 
  theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), strip.background = element_blank(), strip.text.y = element_blank()) +
  ylab("") + facet_grid(Dominance~.)

q2 <- ggplot(dip[dip$Generation < 2000,], aes(x=Generation, y=1-q)) + 
  geom_line(aes(color=as.factor(init_q))) + ggtitle("Diploid") + 
  scale_color_viridis(discrete=T) + theme_bw() + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank()) + 
  ylab("Frequency of Ancestral Allele") + 
  facet_grid(Dominance~.) + xlab("")

r2 <- ggplot(allo[allo$Generation < 2000,], aes(x=Generation, y=1-q)) + 
  geom_line(aes(color=as.factor(init_q))) + 
  scale_color_viridis(discrete=T) + theme_bw() + ggtitle("Allotetraploid") + 
  theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
  ylab("") + facet_grid(Dominance~.) + xlab("")


ancestral_allele_traj = q2 | p2 | r2
ggsave("weak_mut_ancestral_allele_trajectories.pdf", ancestral_allele_traj, width = 6.5, height = 5, units = "in")

p3 <- ggplot(auto[auto$Generation < 2000,], aes(x=Generation, y=w_var)) + 
  geom_line(alpha = 0.5, aes(color=as.factor(init_q))) + 
  scale_color_viridis(discrete=T) + theme_bw() + ggtitle("Autotetraploid") + 
  theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(), strip.background = element_blank(), strip.text.y = element_blank()) + 
  ylab("") + facet_grid(Dominance~.)

q3 <- ggplot(dip[dip$Generation < 2000,], aes(x=Generation, y=w_var)) + 
  geom_line(alpha = 0.5, aes(color=as.factor(init_q))) + ggtitle("Diploid") + 
  scale_color_viridis(discrete=T) + theme_bw() + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank()) + ylab("Fitness Variance") + 
  facet_grid(Dominance~.) + xlab("")

r3 <- ggplot(allo[allo$Generation < 2000,], aes(x=Generation, y=w_var)) + 
  geom_line(alpha = 0.5, aes(color=as.factor(init_q))) + 
  scale_color_viridis(discrete=T) + theme_bw() + ggtitle("Allotetraploid") + 
  theme(legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
  ylab("") + facet_grid(Dominance~.) + xlab("")


fit_var_traj = q3 | p3 | r3
ggsave("weak_mut_fitness_variance_trajectories.pdf", fit_var_traj, width = 6.5, height = 5, units = "in")
