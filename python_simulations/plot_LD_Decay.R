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

auto_infile$p
allo_infile$p
dip_infile$p

auto_infile$ploidy <- rep("Auto")
allo_infile$ploidy <- rep("Allo")
dip_infile$ploidy <- rep("Diploid")

#Calculate Expected Genotype Frequencies at HWE
auto_infile$G0HWE <- auto_infile$p^4
auto_infile$G1HWE <- 4*auto_infile$p^3*(1-auto_infile$p)
auto_infile$G2HWE <- 6*auto_infile$p^2*(1-auto_infile$p)^2
auto_infile$G3HWE <- 4*auto_infile$p*(1-auto_infile$p)^3
auto_infile$G4HWE <- (1-auto_infile$p)^4
auto_infile$SSG <- (auto_infile$G0 - auto_infile$G0HWE)^2 + 
  (auto_infile$G1 - auto_infile$G1HWE)^2 + (auto_infile$G2 - auto_infile$G2HWE)^2 + 
  (auto_infile$G3 - auto_infile$G3HWE)^2 + (auto_infile$G4 - auto_infile$G4HWE)^2

allo_infile$G0HWE <- allo_infile$p^4
allo_infile$G1HWE <- 4*allo_infile$p^3*(1-allo_infile$p)
allo_infile$G2HWE <- 6*allo_infile$p^2*(1-allo_infile$p)^2
allo_infile$G3HWE <- 4*allo_infile$p*(1-allo_infile$p)^3
allo_infile$G4HWE <- (1-allo_infile$p)^4
allo_infile$SSG <- (allo_infile$G0 - allo_infile$G0HWE)^2 + 
  (allo_infile$G1 - allo_infile$G1HWE)^2 + (allo_infile$G2 - allo_infile$G2HWE)^2 + 
  (allo_infile$G3 - allo_infile$G3HWE)^2 + (allo_infile$G4 - allo_infile$G4HWE)^2


dip_infile$G0HWE <- dip_infile$p^2
dip_infile$G2HWE <- 2*dip_infile$p*(1-dip_infile$p)
dip_infile$G4HWE <- (1-dip_infile$p)^2
dip_infile$SSG <- (dip_infile$G0 - dip_infile$G0HWE)^2 + 
  (dip_infile$G2 - dip_infile$G2HWE)^2 + (dip_infile$G4 - dip_infile$G4HWE)^2

#Calculate Expected Gamete Frequencies at HWE
auto_infile$g0HWE <- auto_infile$p^2
auto_infile$g1HWE <- 2*auto_infile$p*(1-auto_infile$p)
auto_infile$g2HWE <- (1-auto_infile$p)^2
auto_infile$SSg <- (auto_infile$g0 - auto_infile$g0HWE)^2 + 
  (auto_infile$g1 - auto_infile$g1HWE)^2 + (auto_infile$g2 - auto_infile$g2HWE)^2

allo_infile$g00HWE <- allo_infile$p^2
allo_infile$g01HWE <- allo_infile$p*(1-allo_infile$p)
allo_infile$g10HWE <- (1-allo_infile$p)*allo_infile$p
allo_infile$g11HWE <- (1-allo_infile$p)^2
allo_infile$SSg <- (allo_infile$g00 - allo_infile$g00HWE)^2 + 
  (allo_infile$g01 - allo_infile$g01HWE)^2 + (allo_infile$g10 - allo_infile$g10HWE)^2 + 
  (allo_infile$g11 - allo_infile$g11HWE)^2

dip_infile$g0HWE <- dip_infile$p
dip_infile$g2HWE <- (1-dip_infile$p)
dip_infile$SSg <- (dip_infile$g0 - dip_infile$g0HWE)^2 + (dip_infile$g2 - dip_infile$g2HWE)^2

#Let's look at deviations
auto_subset <- auto_infile[,c('X..Generation', 'G0', 'G1', 'G2', 'G3', 'G4', 'ploidy', 'SSG', 'SSg')]
allo_subset <- allo_infile[,c('X..Generation', 'G0', 'G1', 'G2', 'G3', 'G4', 'ploidy', 'SSG', 'SSg')]
dip_subset <- dip_infile[,c('X..Generation', 'G0', 'G1', 'G2', 'G3', 'G4', 'ploidy', 'SSG', 'SSg')]

df <- rbind(auto_subset, allo_subset, dip_subset)
names(df) <- c('Generation', 'G0', 'G1', 'G2', 'G3', 'G4', 'ploidy', 'SSG', 'SSg')
df_melt <- melt(df, id.vars = c('Generation', 'ploidy', 'SSG', 'SSg'))
names(df_melt) <- c('Generation', 'ploidy', 'SSG', 'SSg', 'Genotype', 'Frequency')

ggplot(df_melt[df_melt$Generation < 10,], aes(Generation, Frequency)) + 
  geom_line(aes(linetype=ploidy, color = Genotype)) + 
  theme_bw()

ggplot(df[df$Generation < 10,], aes(Generation, sqrt(SSG))) + 
  geom_line(aes(linetype=ploidy)) + 
  theme_bw()

ggplot(df[df$Generation < 10,], aes(Generation, sqrt(SSg))) + 
  geom_line(aes(linetype=ploidy)) + 
  theme_bw()