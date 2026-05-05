**This repo contains publicly available code for the paper** "Improving Minority Population Sampling with BISG Probabilities: Evidence from a Survey of Jewish Americans" by Kyla Chasalow, Eitan Hirsh, Kosuke Imai, and Laura Royden (2026)

Due to the sensitive nature of developing surname-based membership probabilities for a minority population, we are not publicly releasing the obituary data or probabilities learned from the obituary data. However, this repo contains a fully simulated version of our pipeline and methods constructed using data from the American Jewish Population Project to minimic some of the features of our actual application.

---------

**Guide to overall file system**

----------------------

`data` 
  
  * `data/AJPP_data/` data downloaded from Brandeis University’s American Jewish Population Project (https://ajpp.brandeis.edu/) on August 10, 2024 and our processed version
  
  * `data/AJPP_objects` key probability distributions we calculated using the AJPP data
  
  * `data/simulated` directory where data created by simulation are written

`src` 
  
  * `dataprep` - scripts that do AJPP data processing and some general processing helpers
  
  * `sampler` - scripts that implement sampler for heirarchical Bayesian model
  
  * `simulation`  - scripts to run simulated pipeline. 

`figures` is for any figures created by scripts in src

`results` is for any outputs from simulation that are not figures

---

**Guide to `src/simulation` scripts**

---------

**Helper files**

1. `data_generators.R`
2. `sampling_tools.R`
3. `simulation_tools.R`

**Main run files (in order)**

1. `1_simulate_surname_lists.R` - build stable lists of fake surnames
2. `2_simulate_data_generation_plots.R` - generates plots showing features of data simulation approach
3. `3_simulate_population_sampling` - simulate process of sampling from a population
4. `4_simulate_surname_MCMCprob_estimation` - script for playing around with MCMC sampler to estimate surname distribution on simulated data. Meant to be interactive and only for a single case
5. `5_simulate_surname_MCMCprob_estimation_grid` - same as 4. only now doing it systematically to produce final simulation results, while varying various properties of the set-up. Optionally, set `save_raw_output` = TRUE to save the individual chains for each setting run
6. `6_analyze_grid_results` - calculate some metrics and creates various plots based on simulation results files. 
7. `7_analyze_grid_results2` - if in 5_ script, `save_raw_output` was set to TRUE, then can use this script to get things like trace plots for each individual setting run. WARNING: this generates a lot of plots if the simulation grid is large.
