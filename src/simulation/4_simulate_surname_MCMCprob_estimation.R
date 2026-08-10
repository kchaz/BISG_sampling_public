#'------------------------------------------------
#' This script contains main calls to sampler for
#' running sampler on a single simulated dataset
#' (in contrast to the script for running a whole loop)
#'------------------------------------------------

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
source("../sampler/HBayesSampler/sampler_analysis_helpers.R")

#TODO: use this instead
sourceCpp("../sampler/HBayesSampler/sampler.cpp")

#NOTE: for inital testing and development
source("../sampler/HBayesSampler/sampler_Rversion.R")

#TODO: temp
# par(mfrow= c(1,2))
# nrep = 10000
# holder1 = rep(NA, nrep)
# holder2 = rep(NA, nrep)
# for(i in 1:nrep){
#   out = R_sample_bounded_line(.2,.3,0,0,1,1)
#
#   #holder1[i] = out$lratio
#   holder1[i] = out$new[1]
#
#   out = sample_bounded_line(.2,.3,0,0,1,1)
#   holder2[i] = out$new[1]
# }
# hist(holder1)
# hist(holder2)


################################################################################
# Set-up for data generation
################################################################################

# State information
df_states <- read_csv("../../Data/states.csv")
df_states = df_states[order(df_states$Abbreviation),] #make sure alphabetical order by abbreviation, not state name!
n_states <- dim(df_states)[1]

# AJPP info
pR1 = readRDS("../../Data/DataObjects/pR1.RDS")
p_state_gvn_R1 = readRDS("../../Data/DataObjects/p_state_gvn_R1.RDS")
p_state_gvn_R0 = readRDS("../../Data/DataObjects/p_state_gvn_R0.RDS")
p_R1_gvn_state = readRDS("../../Data/DataObjects/p_R1_gvn_state.RDS")
p_R0_gvn_state = readRDS("../../Data/DataObjects/p_R0_gvn_state.RDS")

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

# or loading if only ran this part
#path = paste0("simulated_data/20251222114716_10000names_obitdata/")
#temp = stringr::str_extract(path, "(?<=_)\\d+(?=names)names")
#temp = gsub("names", "", temp)
#num_sur = as.numeric(temp)
#surnames <- readRDS(glue("simulated_data/surnames{num_sur}.RDS"))

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

#TODO: loop over all three options
#init_a <- gam / sum(gam) # very good initial guess
#init_a = rep(1/num_sur, num_sur) #uniform initial guess-still converges after burn-in
init_a <- rdirichlet(1,rbinom(length(gam),rbinom(1,100,.5),.5))

#' ------------------------------------
#' Run Sampler 1
#' ------------------------------------
niter <- 200
start = Sys.time()
out <- runMCMC1(D, #for R version, add R_ in front
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


# TEMP: dev for rho chain plot
# rho_output = analyze_rho(out, ng, burnin = b, upper = NULL)
# rho_output
#
#
# out[["rho_mat"]] = rho_output$rho_mat
#
#
# output_list = list(out)
#
# rho_mat_list <- lapply(output_list, function(x){x$rho_mat})
# rho_chains <- abind::abind(rho_mat_list, along = 3)
# dimnames(rho_chains) <- list("Iteration" = NULL,
#                              "State" = colnames(output_list[[1]]$rho_mat),
#                              "Chain" = 1)
#
# rho_chains = aperm(rho_chains, perm = c(1,3,2))  # Reshape to fit required order
#
# state_subset = c("NY", "FL","TX","WY")
# p = mcmc_trace(rho_chains, pars = state_subset)
# p




# what if I set c = 1 for all? Answer: MSE is larger... the larger c is, the
# more regularize towards n_s proportions
# out.copy = out
# out.copy$c_chain = rep(1, length(out$c_chain))
# Theta_means1 = get_thetameans1(D, ng, out.copy, burnin = b)
# m1 = get_theta_MSEs(Theta_means1, D, true_prob_mat, burnin = b)
# mean(m1$mse_postmean)
# mean(m1$mse_mle)

####################################################################
# SAMPLER 2
# note: sampler 2 may have a bug but I am not fixing it right now
####################################################################
# init_Theta <- t(rdirichlet(dim(D)[2], init_a)) #option 1
# init_Theta <- matrix(abs(rnorm(prod(dim(D)))), ncol = dim(D)[2], nrow = dim(D)[1])
# init_Theta <- D + 1 # 0 count states initialized to uniform
# init_Theta <- apply(init_Theta, 2, function(v) {
#   return(v / sum(v))
# })
# stopifnot(all(round(apply(init_Theta, 2, sum), 10) == 1))
# stopifnot(all(init_Theta > 0))
# stopifnot(all(init_Theta < 1))
#' #####################################
#' Run Sampler 2
#' #####################################
# niter <- 5000
# start = Sys.time()
# out2 <- runMCMC2(
#   D = D,
#   gam = gam,
#   niter = niter,
#   init_a = init_a,
#   init_c = init_c,
#   init_Theta = init_Theta,
#   lambda1 = 1,
#   lambda2 = 1/100,
#   u = .01,
#   testmode = T,
#   verbose = T
# )
#' #####################################
#' Analyse Sampler 2 Output
#' #####################################
# # accept rates
# get_sampler2_accept_stats(out2, b)
#
# # plots
# plot_c_chain(out2, burnin = b)
# plot_a_chain(out2, indices = 1:9, burnin = b) #can plot 9 at a time
#
# par(mfrow = c(2,3))
# plot_theta_chain(out2, 1, 1, burnin = b)
# plot_theta_chain(out2, 1, 2, burnin = b)
# plot_theta_chain(out2, 1, 3, burnin = b)
# plot_theta_chain(out2, 2, 1, burnin = b)
# plot_theta_chain(out2, 2, 2, burnin = b)
# plot_theta_chain(out2, 2, 3, burnin = b)
#
# # posterior means
# pmeans = get_alphac_postmeans(out2, burnin = b)
# print(pmeans$c_mean)
#
# a_means = pmeans$a_means
# get_a_error_analysis(D, a_means, true_prob_marg, plot = b)
#
# Theta_means2 = get_thetameans2(out2, burnin = b)

# # Look at how Theta's compare between samplers 1 and 2
# par(mfrow = c(1,1))
# mean(abs(Theta_means1 - Theta_means2))
# max(abs(Theta_means1 - Theta_means2))
# hist(Theta_means1 - Theta_means2, breaks = 50)
