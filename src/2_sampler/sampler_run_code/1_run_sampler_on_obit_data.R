#'
#' Apply MCMC sampler based on hierarchical Bayesian model to data
#'
#' Author: Kyla Chasalow
#'
# Libraries
library(readr)
library(gtools) # for rdirichlet
library(Rcpp)

# Scripts
source("paths.R")
source(fs::path(src_path,"2_sampler/HBayesSampler/data_prep_helpers.R"))
source(fs::path(src_path,"2_sampler/HBayesSampler/sampler_analysis_helpers.R"))
sourceCpp(fs::path(src_path,"2_sampler/HBayesSampler/sampler.cpp"))

###############################
# Load Data and Key Components
###############################
df = read_csv(paste0(obit_data_path,"obit_list_092424-fulldehyph-check.csv"))
n = nrow(df)

# get list of *all* 51 state names (+DC)
states <- read_csv(paste0(data_path,"states.csv"))
loc_names <- sort(states$Abbreviation)


# get list of surnames
surnames = unique(df$last)
num_sur = length(surnames)

# get count matrix
D <- get_count_mat(df, loc_names,  surnames,surname_colname = "last")
ns <- apply(D, 1, sum)
ng <- apply(D, 2, sum)


#############################
# Init param and hyperparam
#############################
gam <- ns + 1
# very good initial guess + some noise
noise = rbinom(length(gam),100,.5)
init_a <- (gam + noise) / sum(gam + noise)
init_eta <- rbinom(1,100,.5)


####################
# Run Sampler
####################
niter <- 45000
start = Sys.time()
out <- runMCMC(
   D,
   ng = ng,
   niter = niter,
   init_a,
   init_eta,
   gam = gam,
   lambda1 = 1,    #gamma prior parameter. 1 makes it exponential
   lambda2 = 1/100, #gamma prior parameter
   testmode = T,
   verbose = T,
   increment = 5, #save every 5th
   means_burnin = 15000 #don't start accumulating means until here
)
end = Sys.time()
print(end-start)

# Add back some dimension names
rownames(out$a_chain) = rownames(D)
colnames(out$theta_sums) = colnames(D)
rownames(out$theta_sums) = rownames(D)
dimnames(out$a_chain) <- list("Parameter" = rownames(out$a_chain),
                              "Iteration" = NULL)



#SAVE OUTPUT
formatted_time <- format(Sys.time(), "%Y%m%d%H%M%S")
saveRDS(out, paste0(save_sampler_output_path, formatted_time, "_", "realdf_", niter, ".RDS"))
