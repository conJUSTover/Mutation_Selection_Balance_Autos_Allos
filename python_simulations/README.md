# Python Simulations

This directory contains both simulation scripts and code to generate all data and figures associated with these simulation scripts. More details on each subdirectory are outlined below. 

## simulation_scripts

This subdirectory contains the discrete time, generation-based Python simulations for our model. There is a simulation script for each model (diploids, autotetraploids, and allotetraploids). 

This simulation avoids the use of a limit to obtain continuous time ODEs and may be more suitable for simulations with very strong mutation or selection. 

## panmictic_disequilibrium_decay

This subdirectory contains Bash commands, data, and an R plotting script for Figure 1 in the manuscript. 

Specifically, we simulate a short period (20 generations) and plot the decay of PD over the first 10 of those generations. To recreate the data, run the Bash commands included in cmd in the terminal. To recreate the figure, run the plot_LD_decay.R script.

## bistability
This subdirectory contains Bash commands, data, and an R plotting script for Figure 5 in the manuscript and Figures S8, S9 in the supplement. 

Specifically, we simulate 3000 generations and plot the fitness, fitness variance, and (derived) allele frequency for the first 2000 of those generations. 

This subdirectory is further split by dominance coefficients in the recessive, additive, and dominant sub-subdirectories. 

To recreate the data, within the folder for each dominance case, run the Bash commands included in cmd, cmd_dip, and cmd_allo in the terminal. Then, run the cmd2 and cmd2_allo Bash commands. To recreate the figure, run the bistability_approach.R script. 

## weak_mutation
This subdirectory contains Bash commands, data, and an R plotting script for Figure 4 in the manuscript and Figures S6, S7 in the supplement. 

Specifically, we simulate 3000 generations and plot the fitness, fitness variance, and (derived) allele frequency for the first 2000 of those generations. 

This subdirectory is further split by dominance coefficients in the recessive, additive, and dominant sub-subdirectories. 

To recreate the data, within the folder for each dominance case, run the Bash commands included in cmd, cmd_dip, and cmd_allo in the terminal. Then, run the cmd2 and cmd2_allo Bash commands. To recreate the figure, run the bistability_approach.R script. 
