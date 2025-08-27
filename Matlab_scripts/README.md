# Matlab Scripts

This directory contains both continuous time scripts for analyzing the system of ODEs we consider and code to recreate all data and figures associated with this analysis. More details on each subdirectory are outlined below. 

## base_symbolic_models

This subdirectory contains base models for the continuous time analysis. This uses Matlab's symbolic toolbox and associated vpasolve and eig functions to numerically estimate fixed points of mutation-selection balance and classify their stability. This analysis relies on a continuous itme limit and may not be suitable for very large mutation and selection. 

## equal_mut_bifn

This subdirectory contains additional Matlab functions and scripts and an R plotting script to recreate the data and plots for Figure 2 in the manuscript and Figures S1, S2, and S3 in the supplemental material. The data and scripts directly held within this directory have mu = 1e-8 and are used in all four figures. This subdirectory contains two sub-subdirectories (mu_1e-6 and mu_1e-7) which have part of the data used in Figures S1 and S2. 

Specifically, we numerically estimate the equilibrium solutions of the ODEs and classify their stability under an equal mutation rates regime. We plot allele frequency and mutation load as summary statistics. 

To recreate the data, run equal_mut_comparison_data.m in this subdirectory and within mu_1e-6 and mu_1e-7. To recreate all four figures, run plot_equal_mut.R.

## deleterious_bifn

This subdirectory contains additional Matlab functions and scripts and an R plotting script to recreate the data and plots for Figure 3 in the manuscript.

Specifically, we numerically estimate the equilibrium solutions of the ODEs and classify their stability under a biased mutation rates regime. This admits a pair of saddle-node bifurcations for the fully dominant case. We plot allele frequency and mutation load as summary statistics. 

To recreate the data, run generate_bifn_data.m. To recreate the figure, run plot_bifn.R. 

## Phase_planes

This subdirectory contains additional Matlab functions and scripts to recreate the data and plots for Figures S4 and S5 in the supplemental material.

Specifically, we plot the differential equations decomposed into selection and mutation fields against allele frequency for the diploids. For the autotetraploids, we plot three phase planes with the ODEs represented as vector fields. This provides a visualization of what the dynamics look like in the phase plane as the strength of selection changes. 

To recreate the figures, simply run diploids_phase_plane_fields.m and phase_plane_autos.m. 

## ODE15s_scripts

This subdirectory contains additional scripts for numerically solving the systems of ODEs using the Matlab function ode15s. 

Because random mating and selection/mutation affect the dynamics on drastically different time scales, the system of ODEs is stiff and requires use of a stiff solver. We found that the built-in ode15s routine in Matlab (variable-step, variable-order) works well for our problem. Additionally, to maintain numerical stability and ensure positivity of the resulting solutions, we logarithmically transformed the gamete frequencies to ensure their positivity during solving. These functions plot the resulting allele, gamete, and genotype frequencies over time for arbitrary selection and mutation parameters. The key advantage here is that ode15s is much faster to run for weak mutation and selection compared to our discrete time, generation-based Python simulations. This allows us to track the trajectories of the ODEs as they approach an equilibrium (whereas the symbolic scripts only determine the end equilibrium state and the Python simulation is prohibitively slow and memory consuming to run for such long periods). 

However, do note that for very large mutation or selection values, this approach may not be suitable. Please use the Python simulations in this case. 
