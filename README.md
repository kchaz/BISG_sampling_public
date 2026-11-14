**This repo contains publicly available code for the paper** "Improving Minority Population Sampling with BISG Probabilities: Evidence from a Survey of Jewish Americans" by Kyla Chasalow, Eitan Hirsh, Kosuke Imai, and Laura Royden (2026) available at https://arxiv.org/abs/2605.05384

Due to the sensitive nature of developing surname-based membership probabilities for a minority population, we are not publicly releasing the obituary data or probabilities learned from the obituary data. However, this repo contains a fully simulated version of our pipeline and methods constructed using data from the American Jewish Population Project to mimic some of the features of our actual application.


# Overall repo structure and a comment on file paths

Much of the code (particularly within the src/3_ folder below) was run on Harvard's cluster, where the voter file data are located, so some filepaths are relative to that cluster. Because we are not able to release the voter file data itself (so that code is not runnable anyway), we have not changed these paths. The rest of the filepaths are relative to a Github repo (link to public version here https://github.com/kchaz/BISG_sampling_public) which assumes that all code is located within src/ and all data in sub-directories of a directory called data/


# Guide to Code

Given all data files, the main analysis of the paper could be replicated by running the following numbered sub-directories of `src` in order.
In addition, there is a simulation directory containing simulations.

---------------------
`1_data_processing` 
---------------------
The scripts in this directory run our initial data cleaning processes on each of the data sources, including cleaning the obituary data and calculating P(Jewish|State) and P(State|Jewish) from the AJPP estimated counts. There is also a script that we used for scraping names from Wikipedia.

-------------
`2_sampler` 
------------
The scripts in this directory are for running an collapsed MCMC sampler for fitting the Bayesian model described in the paper. /HBayesSampler contains the general code for implementing the method. /sampler_run_code contains our calls for fitting the model on the obituary data.

* /HBayesSampler
	- `sampler.cpp` and `sampler_Rversion.R` - main scripts for running the sampler. The .cpp version should contain all the same functions as the R version, only implemented in Rcpp, which substantially increases speed. The main functions are runMCMC (for .cpp file) and R_runMCMC (for .R file). The names of the arguments are aligned with the variable names from the paper with a for alpha.
	- `sampler_analysis_helpers` - functions for analyzing output of the sampler functions
	- `data_prep_helpers.R` - a function for getting a dataframe with (State, surname) entries and converting it to a $n_{s,g}$ count matrix.

* /sampler_run_code
	- `1_run_sampler_on_obit_data.R` - run sampler on obituary data
	- `2_analyze_sampler_output.R` - get some statistics, build some plots
	- `3_reshape_output_to_dataframe.R` - take sampler output, calculate the actual posterior means of the \theta_{s,g} = P(S=s|G=g,R=1) probabilities, and then reshape output from sampler into a long format with columns surname, state, probability. Then join these with the P(R=1|G) estimates to get the P(S=s|G=g,R=1)P(R=1|G) numerator needed for calculating our final probabilities
	- `4_get_bounds.R` - script for calculating bounds on possible values of P(S=s|G=g,R=1) implied by P(R=1|G) and P(S|G) from voter file. These are stored so that they can later be used to deal with issue of incoherence between data sources that leads to estimated BISG probabilities outside of [0,1]. Note: there is a 'conservative' bound option here that we did not end up actually using (it doesn't guarantee probs in [0,1]) but explored. Note: this script uses data found only on the cluster.
		* Bound 1: non-conservative - this ensures that downstream $P(R=1|S,G)$ calculated is not $> 1$  - uses point estimate of $\Pr(R=1|G)$
		* Bound 2: conservative - this bound uses upper/lower confidence interval bounds on $\Pr(R=1|G)$ from the AJPP data calculation in recognition that point estimate could be wrong. The problem is that this does not gaurantee the downstream estimate of P(R=1|S,G) is in [0,1]. We did not end up using this.
	- `5_impose_bounds.R` - script for taking bounds from previous script and actually imposing them.

Unordered scripts:
	- `paths.R` - defines some paths where different data and functions are located. These are paths relative to the Github repo by default but when running this code using a Bash script, one can optionally feed custom path arguments and this script will then use those. 
	- `analyze_bounds.R` - quick script for exploring the bounds calculated in 4_
	- `examine_output.R` - load stats on proposal acceptance rate calculated and saved from running 1_ script

---------------------------
`3_sample_from_voterfile` 
-------------------------
These are the scripts for actually sampling from the voter file. The outputs from the sampler from 2_sampler/sampler_run_code/ are one
input into this

* Ordered scripts
	- `1_full_pipeline_code.R` -  Implements the full pipeline from calculating BISG probabilities to allocating by state to taking the final sample. Assumes that sampler from 2_sampler/ file has been run. Paths are specific to the data on the cluster
	- `2_run_checks.R` - runs some checks on some output created by previous script
	- `3_merge_in_unnormalized_probabilities.R` - a post-hoc script for merging the unnormalized (i.e. original BISG) estimates to the final sample. We had originally not included these but then realized we wanted them so we could analyze the distribution of the BISG probabilities of those sampled.
	
---------------------
`4_analyze_survey` 
---------------------
These are the scripts for analyzing the survey data and calculating survey weights

`1_counts_for_weights.R`
	Used to get gender, party, and age information from voter file merged with the BISG probabilities estimated previously

`2_estimate_weights_from_voterfile.R`
	Uses voter file and BISG probabilities to estimate P(X=x|State, Jewish) for X = gender, party, age

`3_aggregate_weights_across_states.R`
	Uses outputs of 1_ script to calculate P(X=x|Jewish) and then creates assorted plots to visualize these and compare them to Pew

`4_analyze_survey.Rmd`
	Script used to produce all the final survey estimates reported in the paper. 

`5_process_survey_for_release`
	Script used to produce released survey.
	NOTE: that in the released survey, the name of the BISG probabilities column was changed from a not so informative internal name
	(total_probs_bounded) to a more interpretable one (BISG)

-----------------------------
`5_calculations_for_paper` 
-----------------------------
These are assorted scripts for calculating statistics and creating plots reported in the paper. 
They are ordered roughly be theme (e.g., data validation, allocation calculations, post-survey analysis) but
aside from the 1-1_ and 1-2_ scripts, there is no required order for running these scripts.

------------
simulation
------------
This is for running a fully simulated version of our analysis process.

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


----
EDA
----
These are a few assorted scripts we used early on in the project to explore some properties of the datasets.
