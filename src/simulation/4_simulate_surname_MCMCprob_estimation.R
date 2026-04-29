#'----------------------------------------------------------------------------
#' This script contains main calls for generating simulated obituary data
#' and running sampler on a single simulated dataset
#' (in contrast to the script 5_ for running a loop over different settings)
#'
#' This script is meant to be exploratory and was used in development
#'
#'-----------------------------------------------------------------------
# If TRUE, will save the generated dataframe and other objects for inspection
SAVE_SIMULATION_OBJECTS = FALSE

if(SAVE_SIMULATION_OBJECTS){
  save_path <- paste0("../../data/simulated/")
}

# Libraries
library(readr)
library(gtools) # for rdirichlet
library(Rcpp)
library(abind)
library(glue)

set.seed(21)
source("data_generators.R") #functions for simulating data
source("simulation_error_calculators.R")

# Load relevant code scripts - set working directory to source first
source("../dataprep/data_processing_helpers.R")
source("../sampler/sampler_analysis_helpers.R")

#NOTE: for initial testing and development, used R version
#source("../sampler/HBayesSampler/sampler_Rversion.R")
sourceCpp("../sampler/sampler.cpp")

################################################################################
# Set-up for data generation
################################################################################

# State information
df_states <- read_csv("../../data/states.csv")
df_states = df_states[order(df_states$Abbreviation),] #make sure alphabetical order by abbreviation, not state name!
n_states <- dim(df_states)[1]

# AJPP info
pR1 = readRDS("../../data/AJPP_objects/pR1.RDS")
p_state_gvn_R1 = readRDS("../../data/AJPP_objects/p_state_gvn_R1.RDS")
p_state_gvn_R0 = readRDS("../../data/AJPP_objects/p_state_gvn_R0.RDS")
p_R1_gvn_state = readRDS("../../data/AJPP_objects/p_R1_gvn_state.RDS")
p_R0_gvn_state = readRDS("../../data/AJPP_objects/p_R0_gvn_state.RDS")

# Number of surnames to use, can go up to 10000
num_sur <- 1000

# Surnames
surnames <- readRDS(paste0(glue("../../data/simulated/surnames10000.RDS")))
surnames = surnames[1:num_sur]

# Generate true surname distributions for R=1 group
# These are the the estimands we want to recover
true_prob_mat = get_pmat(df_states$Abbreviation, surnames, rate = 1)
true_prob_marg = get_p_surname_gvn_R1(true_prob_mat, p_state_gvn_R1)


################################################################################
# Simulate Obituary Data
################################################################################
#in real case, have ~200,000 obs of 50,000 surnames so here, create a 1/4 ratio

n <- 100000
df <- create_df(true_prob_mat, n = n, props = p_state_gvn_R1, fips_mode = F)

# Reshape into count matrix
state_names <- sort(df_states$Abbreviation)
D <- get_count_mat(df, state_names, surnames)
ns <- apply(D, 1, sum)
ng <- apply(D, 2, sum)

# Checks
stopifnot(sum(D) == n)
stopifnot(length(ns) == num_sur)
stopifnot(length(ng) == nrow(df_states))

################################################################################
# Saves if requested
################################################################################
if(SAVE_SIMULATION_OBJECTS){
  write.csv(df, glue("{save_path}obit_data.csv"), row.names = F)
  saveRDS(true_prob_mat, glue("{save_path}true_prob_mat.RDS"))
  saveRDS(true_prob_marg, glue("{save_path}true_prob_marg.RDS"))
}

################################################################################
# Run sampler
################################################################################

#' ------------------------------------
# Initialize parameters
#' ------------------------------------
#' Gamma is a global hyper-parameter. I use inflated, unnormalized counts
#' because better behavior of Dirichlet when param are > 1. Had some issues
#' when did not have the +1 because then had some = 1
#' --------------------------------
gam <- ns + 1
init_eta <- 1
init_a <- rdirichlet(1,rbinom(length(gam),rbinom(1,100,.5),.5))

#' ------------------------------------
#' Run Sampler 1
#' ------------------------------------
niter <- 50
start = Sys.time()
out <- runMCMC1(D, #for R version, add R_ in front
  ng = ng,
  niter = niter,
  init_a, init_eta,
  gam = gam,
  lambda1 = 1, #gamma prior parameter. 1 makes it exponential
  lambda2 = 1/100, #gamma prior parameter
  testmode = T,
  verbose = T,
  increment = 1,
  means_burnin = 1 #don't start calculating theta posterior means until after this value
)
end = Sys.time()
print(end-start)

#add some dimension naming for later
rownames(out$a_chain) = rownames(D)
colnames(out$theta_sums) = colnames(D)
rownames(out$theta_sums) = rownames(D)
dimnames(out$a_chain) <- list("Parameter" = rownames(out$a_chain),
                              "Iteration" = NULL)


##########################################################################
# Quick functions for plotting and calculating some stats
# + looking at posterior means and comparing to truth
##########################################################################
upper_cutoffs = NULL#c(50,100)

b = 1 #burn-in

# accept rates
get_sampler1_accept_stats(out, burnin = b)

# plots
plot_eta_chain(out, burnin = b)
plot_a_chain(out, indices = 1:9, burnin = b) #can plot 6 at a time

# posterior means
pmeans = get_alpha_and_eta_postmeans(out, burnin = b, upper_cutoffs)
print(pmeans$eta_mean)

# posterior means of alpha's - MSE improve?
a_means = pmeans$a_means
a_output = get_alpha_error_analysis(D, a_means, true_prob_marg)
a_output

# posterior means of thetas - MSE improve compared to MLE
theta_output = analyze_theta_for_multiple_cutoffs(D = D,
                                   ng = ng,
                                   out = out,
                                   true_prob_mat = true_prob_mat,
                                   burnin = b,
                                   upper_cutoffs = upper_cutoffs)

print(theta_output)

rho_output = analyze_rho(out, ng, burnin = b, upper = NULL)
rho_output$rho_postmean
