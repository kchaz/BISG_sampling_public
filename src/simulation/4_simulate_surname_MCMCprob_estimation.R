#
# Author: Kyla Chasalow
#
#-------------------------------------------------------------------------------
# This script contains main calls to sampler for
# running sampler on a single simulated dataset
# (in contrast to the script for running a whole loop of samplers over
# different settings)
#-------------------------------------------------------------------------------

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
source("../2_sampler/HBayesSampler/data_prep_helpers.R")
source("../2_sampler/HBayesSampler/sampler_analysis_helpers.R")
source("../2_sampler/HBayesSampler/sampler_Rversion.R")
sourceCpp("../2_sampler/HBayesSampler/sampler.cpp")

################################################################################
# Set-up for data generation
################################################################################

# State information
df_states <- read_csv("../../data/states.csv")
df_states = df_states[order(df_states$Abbreviation),] #make sure alphabetical order by abbreviation, not state name!
n_states <- dim(df_states)[1]

# AJPP info
pR1 = readRDS("../../data/DataObjects/pR1.RDS")
p_state_gvn_R1 = readRDS("../../data/DataObjects/p_state_gvn_R1.RDS")
p_state_gvn_R0 = readRDS("../../data/DataObjects/p_state_gvn_R0.RDS")
p_R1_gvn_state = readRDS("../../data/DataObjects/p_R1_gvn_state.RDS")
p_R0_gvn_state = readRDS("../../data/DataObjects/p_R0_gvn_state.RDS")

# Number of surnames to use
num_sur <- 1000

# Directory for saving attributes of data generating process
formatted_time <- format(Sys.time(), "%Y%m%d%H%M%S")
folder_path <- paste0("simulated_data/",formatted_time,"_", num_sur, "names_obitdata")
if (!dir.exists(folder_path)){
  dir.create(folder_path)
}

# Surnames
surnames <- readRDS(paste0(glue("simulated_data/surnames{num_sur}.RDS")))
surnames = surnames[1:num_sur]
saveRDS(surnames, paste0(folder_path,"/possible_surnames.RDS"))

# Surname distribution for R=1 group
pmat = get_pmat(df_states$Abbreviation, surnames, rate = 1)
saveRDS(pmat, paste0(folder_path, "/p_S_gvn_R1_state.RDS"))

pmarg = get_p_surname_gvn_R1(pmat, p_state_gvn_R1)
saveRDS(pmarg, paste0(folder_path, "/p_S_gvn_R1.RDS"))


################################################################################
# Simulate Obituary Data and save it
################################################################################
#in real case, have a bit more than 200,000 obs of 50,000 surnames
#so here, create a 1/4 ratio

n <- 100000
df <- create_df(pmat, n = n, props = p_state_gvn_R1, fips_mode = F)
write_csv(df, paste0(folder_path, "/ObitData.csv"))


################################################################################
# Run sampler
################################################################################

# Estimands we want to recover
true_prob_mat <- readRDS(paste0(folder_path, "/p_S_gvn_R1_state.RDS"))
true_prob_marg <- readRDS(paste0(folder_path, "/p_S_gvn_R1.RDS"))

# load main dataset and reshape it into count matrix
loc_names <- sort(df_states$Abbreviation)
df <- read_csv(paste0(folder_path, "/ObitData.csv"), show_col_types = F)
n <- nrow(df)
D <- get_count_mat(df, loc_names, surnames)
ns <- apply(D, 1, sum)
ng <- apply(D, 2, sum)

stopifnot(sum(D) == n)
stopifnot(length(ns) == num_sur)
stopifnot(length(ng) == nrow(df_states))


#' ------------------------------------
# Initialize parameters
#' ------------------------------------
#' Gamma is a global hyper-parameter. I use inflated, unnormalized counts
#' because better behavior of Dirichlet when param are > 1. Had some issues
#' when did not have the +1 so that some = 1
#' --------------------------------
gam <- ns + 1
init_eta <- 1

#init_a <- gam / sum(gam) # very good initial guess
#init_a = rep(1/num_sur, num_sur) #uniform initial guess-still converges after burn-in
init_a <- rdirichlet(1,rbinom(length(gam),rbinom(1,100,.5),.5))

#' ------------------------------------
#' Run Sampler 1
#' ------------------------------------
niter <- 200
start = Sys.time()
out <- runMCMC(D, #for R version, add R_ in front
  ng = ng,
  niter = niter,
  init_a, init_eta,
  gam = gam,
  lambda1 = 1,    #gamma prior parameter. 1 makes it exponential
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


# SAVE OUTPUT
#formatted_time <- format(Sys.time(), "%Y%m%d%H%M%S")
#saveRDS(out, paste0("sampler_outputs/", formatted_time, "_sim", num_sur, "_iter", niter, ".RDS"))
#out = readRDS("sim_outputs/20241010131910_sim2000_5000.RDS")


##########################################################################
# Hand written plot func and stats
# and code for looking at posterior means and comparing to truth
##########################################################################
upper_cutoffs = NULL#c(50,100)

b = 10 #burn-in

# accept rates
get_sampler1_accept_stats(out, burnin = b)

# plots
plot_eta_chain(out, burnin = b)
plot_a_chain(out, indices = 1:9, burnin = b) #can plot 6 at a time

# posterior means
pmeans = get_alpha_and_eta_postmeans(out, burnin = b, upper_cutoffs)
print(pmeans$c_mean)

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
