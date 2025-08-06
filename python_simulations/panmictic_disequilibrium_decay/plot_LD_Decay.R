
library(ggplot2)
library(reshape2)

setwd("~/GitHub/Mutation_Selection_Balance_Autos_Allos/python_simulations/panmictic_disequilibrium_decay/")
auto_infile <- read.csv("auto_LD_decay.txt", header = T)
allo_infile <- read.csv("allo_LD_decay.txt", header = T)
dip_infile <- read.csv("dip_LD_decay.txt", header = T)

#Calculate tetraploid genotypes in allos
allo_infile$G0 <- allo_infile$G00
allo_infile$G1 <- allo_infile$G01 + allo_infile$G10 
allo_infile$G2 <- allo_infile$G02 + allo_infile$G11 + allo_infile$G20
allo_infile$G3 <- allo_infile$G12 + allo_infile$G21
allo_infile$G4 <- allo_infile$G22

#calculate ancestral allele frequencies 
auto_infile$p <- auto_infile$g0 + 0.5*auto_infile$g1
allo_infile$p <- allo_infile$g00 + 0.5*(allo_infile$g10 + allo_infile$g10)
dip_infile$p <- dip_infile$g0 + 0.5*dip_infile$g1

auto_infile$Ploidy <- rep("Autotetraploid")
allo_infile$Ploidy <- rep("Allotetraploid")
dip_infile$Ploidy <- rep("Diploid")

#Using the Panmictic Disequilibrium Measure from Gallais 2003 
auto_infile$D <- auto_infile$g0 - auto_infile$p^2
allo_infile$D <- allo_infile$g00 - allo_infile$p^2
dip_infile$D <- dip_infile$g0 - dip_infile$p

#Let's look at deviations
auto_subset <- auto_infile[,c('X..Generation', 'G0', 'G1', 'G2', 'G3', 'G4', 'Ploidy', 'D')]
allo_subset <- allo_infile[,c('X..Generation', 'G0', 'G1', 'G2', 'G3', 'G4', 'Ploidy', 'D')]
dip_subset <- dip_infile[,c('X..Generation', 'G0', 'G1', 'G2', 'G3', 'G4', 'Ploidy', 'D')]

df <- rbind(auto_subset, allo_subset, dip_subset)
names(df) <- c('Generation', 'G0', 'G1', 'G2', 'G3', 'G4', 'Model', 'D')


pdf('Panmictic_Disequilibrium_Decay.pdf', width = 3, height = 3.25)
ggplot(df[df$Generation < 11,], aes(Generation, D)) + 
  geom_line(aes(linetype = Model)) + 
  ylab(bquote(Delta)) + 
  theme_bw() + 
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) + 
  theme(
    legend.key.size = unit(1, "lines"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt"),
    legend.spacing.y = unit(0.1, "lines"),
    legend.position = "bottom"
  ) + 
  guides(linetype = guide_legend(nrow = 2))
dev.off()

