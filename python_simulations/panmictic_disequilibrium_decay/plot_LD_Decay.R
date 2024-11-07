library(ggplot2)
library(reshape2)

setwd("~/GitHub/Mutation_Selection_Balance_Autos_Allos/python_simulations/")
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

auto_infile$Ploidy <- rep("Auto")
allo_infile$Ploidy <- rep("Allo")
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
names(df) <- c('Generation', 'G0', 'G1', 'G2', 'G3', 'G4', 'Ploidy', 'D')



png('Panmictic_Disequilibrium_Decat.png', width = 5, height = 3, res = 300, units = 'in')
ggplot(df[df$Generation < 11,], aes(Generation, D)) + 
  geom_line(aes(linetype=Ploidy)) + ylab(bquote("\u0394"[0])) + 
  theme_bw() + scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10))
dev.off()
