#'
#' Simulate performance of sampler with repeat generation of data,
#' and varying some sampler settings
#' ----------------------------------------------------------------
# Libraries
library(readr)
library(gtools) # for rdirichlet
library(Rcpp)
library(bayesplot)
library(coda) # gelman-rubin diagnostics
library(abind)
library(glue)
library(future.apply)

source("simulation_tools.R")
source("simulation_error_calculators.R")
source("data_generators.R") # functions for simulating data
source("../dataprep/data_processing_helpers.R")
source("../sampler/HBayesSampler/sampler_analysis_helpers.R")
sourceCpp("../sampler/HBayesSampler/sampler.cpp")
source("../sampler/HBayesSampler/sampler_Rversion.R")

#-------------------------------------------------------------------------------
# Things that are fixed throughout
#-------------------------------------------------------------------------------

# State information
df_states <- read_csv("../../Data/states.csv")
df_states <- df_states[order(df_states$Abbreviation), ] # make sure alphabetical order by abbreviation, not state name!
n_states <- dim(df_states)[1]

# AJPP info
pR1 <- readRDS("../../Data/DataObjects/pR1.RDS")
p_state_gvn_R1 <- readRDS("../../Data/DataObjects/p_state_gvn_R1.RDS")
p_state_gvn_R0 <- readRDS("../../Data/DataObjects/p_state_gvn_R0.RDS")
p_R1_gvn_state <- readRDS("../../Data/DataObjects/p_R1_gvn_state.RDS")
p_R0_gvn_state <- readRDS("../../Data/DataObjects/p_R0_gvn_state.RDS")

#-------------------------------------------------------------------------------
# Location and type of output to save, path for results df
#-------------------------------------------------------------------------------

save_path <- "results/"
fs::dir_create(save_path)

true_val_path <- fs::path(glue("{save_path}/true_probabilities"))
fs::dir_create(true_val_path)

per_setting_results_path <- fs::path(glue("{save_path}/per_setting_results"))
fs::dir_create(per_setting_results_path)

# grid_filepath <- fs::path(glue("{save_path}/grid_results.csv"))

tempfile_path <- fs::path(glue("{save_path}/tempfiles"))
fs::dir_create(tempfile_path)


#-------------------------------------------------------------------------------
# Setting up grid to vary over - THINGS THAT CAN VARY
#-------------------------------------------------------------------------------
# NOTE: number of states implicitly fixed at 51
# Real problem has about 230,000 obs x 50000 surnames x 51 states
# Setting that most mimics this is 50,000 obs and 10000 surnames
# though it's still different because this is spread over 51 states
# (worse than our actual then...)

num2run <- 5 # number of replicates per setting

setting_list <- list(
  "num_sur" = c(100, 1000, 10000),
  "n" = c(25000, 50000, 100000),
  "alpha_init_option" = c("uniform", "best_guess")
)

# TODO: currently can't easily restart this and add more here
# TODO: make that possible via changes to tempfile consolidator?
# will calculate posterior means at each of these``
# will run for max value iterations
niter_cutoffs <- c(1, 10,
                   100, 250, 500,
                   1000, 1500, 2000)
#-------------------------------------------------------------------------------
# Controls for what grid should save for each individual setting
#-------------------------------------------------------------------------------

save_raw_output <- T

#-------------------------------------------------------------------------------
# Step 1: generate and save a single true_prob_mat matrix for each num_sur requested
#-------------------------------------------------------------------------------
set.seed(21) # so creates same each time - this is fast, so ok if re-does it when re-start

for (num_sur in setting_list$num_sur) {
  surnames <- readRDS(glue("simulated_data/surnames{num_sur}.RDS"))
  true_prob_mat <- get_pmat(df_states$Abbreviation, surnames, rate = 1)
  true_probs_marg <- get_p_surname_gvn_R1(true_prob_mat, p_state_gvn_R1)
  saveRDS(true_prob_mat, fs::path(glue("{true_val_path}/true_prob_mat-{num_sur}names.RDS")))
  saveRDS(true_probs_marg, fs::path(glue("{true_val_path}/true_probs_marg-{num_sur}names.RDS")))
}

#-------------------------------------------------------------------------------
# Step 2: Define what should happen in a single step
# for each setting, create dataset, fit sampler, save output metrics,
# optionally create plots
#-------------------------------------------------------------------------------

single_run <- function(settings,
                       per_setting_results_path,
                       true_val_path,
                       df_states,
                       niter_cutoffs) {
  num_sur <- as.numeric(settings[["num_sur"]])
  n <- as.numeric(settings[["n"]])
  alpha_init_option <- settings[["alpha_init_option"]]
  save_path <- settings[["filename"]]

  niter <- max(niter_cutoffs)

  num_sur_str <- format(as.numeric(num_sur), scientific = FALSE)
  n_str <- format(as.numeric(n), scientific = FALSE)
  niter_str <- format(as.numeric(niter), scientific = FALSE)

  print(paste0(
    "Running setting: num_sur = ", num_sur_str,
    ", n = ", n_str,
    ", alpha_init_option = ", alpha_init_option,
    ", niter = ", niter_str
  ))

  # Get true parameters of DGP
  true_prob_mat <- readRDS(fs::path(glue("{true_val_path}/true_prob_mat-{num_sur}names.RDS")))
  true_probs_marg <- readRDS(fs::path(glue("{true_val_path}/true_probs_marg-{num_sur}names.RDS")))

  # Generate data under these
  surnames <- readRDS(glue("simulated_data/surnames{num_sur}.RDS"))
  df <- create_df(true_prob_mat, n = n, props = p_state_gvn_R1, fips_mode = F)

  # Convert data to counts
  loc_names <- sort(df_states$Abbreviation)
  D <- get_count_mat(df, loc_names, surnames)
  ns <- apply(D, 1, sum)
  ng <- apply(D, 2, sum)
  stopifnot(sum(D) == n)
  stopifnot(length(ns) == num_sur)
  stopifnot(length(ng) == nrow(df_states))

  # Set global hyperparameter
  gam <- ns + 1

  # Initialize parameters
  init_eta <- 1
  if (alpha_init_option == "uniform") {
    init_a <- rep(1 / num_sur, num_sur)
  } else if (alpha_init_option == "best_guess") {
    init_a <- gam / sum(gam) # very good initial guess
  } else {
    stop(glue("Invalid alpha_init_option = {alpha_init_option} given"))
  }

  # Fit sampler
  start <- Sys.time()
  out <- runMCMC1(D, # for R version, add R_ in front
    ng = ng,
    niter = niter,
    init_a, init_eta,
    gam = gam,
    lambda1 = 1, # gamma prior parameter. 1 makes it exponential
    lambda2 = 1 / 100, # gamma prior parameter
    testmode = T,
    verbose = T,
    increment = 1,
    means_burnin = 10 # don't start calculating theta posterior means until after this value
  )
  end <- Sys.time()

  # add some dimension naming for later
  rownames(out$a_chain) <- rownames(D)
  colnames(out$theta_sums) <- colnames(D)
  rownames(out$theta_sums) <- rownames(D)
  dimnames(out$a_chain) <- list(
    "Parameter" = rownames(out$a_chain),
    "Iteration" = NULL
  )

  # Calculate Metrics
  #--------------------------------------------
  metrics <- vector(mode = "list")

  # first add settings
  metrics[["num_sur"]] <- num_sur
  metrics[["n"]] <- n
  metrics[["alpha_init_option"]] <- alpha_init_option
  metrics[["niter"]] <- niter

  metrics[["sampler_time"]] <- end - start

  b <- 1 # burn-in - just fixing this to essentially not having it for now

  # accept rates
  temp <- get_sampler1_accept_stats(out, burnin = b)
  metrics[["c_rate"]] <- temp$c_rate
  metrics[["a_mean_prop_accept"]] <- temp$a_mean_prop_accept
  metrics[["a_max_accept_prop"]] <- temp$a_max_accept_prop
  metrics[["a_min_accept_prop"]] <- temp$a_min_accept_prop


  # posterior means
  pmeans <- get_alpha_and_eta_postmeans(out, burnin = b, upper_cutoffs = niter_cutoffs)

  # error analysis for alpha - will to separatley for each niter cutoff
  a_output <- get_alpha_error_analysis(D, pmeans$a_means, true_probs_marg)
  metrics <- c(metrics, a_output)

  # error analysis for theta
  theta_output <- analyze_theta_for_multiple_cutoffs(
    D = D,
    ng = ng,
    out = out,
    true_prob_mat = true_prob_mat,
    burnin = b,
    upper_cutoffs = niter_cutoffs
  )
  metrics <- c(metrics, theta_output)

  # analyze eta
  metrics[[glue("eta_postmean_{max(niter_cutoffs)}")]] = pmeans$eta_means[1,length(niter_cutoffs)] #get the max iter one

  # function to calculate rho
  rho_output = analyze_rho(out = out, ng = ng, burnin = b, upper = NULL) #just do full , haven't written it to allow multiple cutoffs yet
  metrics[["postmean_rho_mean"]] = rho_output$postmean_rho_mean
  metrics[["postmean_rho_sd"]] = rho_output$postmean_rho_sd
  metrics[["postmean_rho_min"]] = rho_output$postmean_rho_min
  metrics[["postmean_rho_max"]] = rho_output$postmean_rho_max
  metrics[["rho_postmean_over_individuals"]] = rho_output$rho_postmean_over_individuals

  # save the rho mat of values for each state for all iterations (excluding burn-in)
  out[["rho_mat"]] = rho_output[["rho_mat"]]

  # Add filename to for record keeping and restart capability
  metrics["filename"] <- save_path

  # Write a row to a dataframe
  # NOTE: warning - this assumes same order of columns as existing file
  # NOTE: update on previous warning - does not matter actually because writing a new file for each run.
  # and my merge file can deal with different col in different orders
  row <- as.data.frame(metrics)
  if (file.exists(save_path)) { # add to file if it exists
    write.table(row,
      file = save_path,
      sep = ",",
      row.names = F,
      col.names = F,
      append = T
    )
  } else {
    write.csv(row, save_path, row.names = F)
  }

  # Create a directory for saving output things specific to this setting
  setting_path <- fs::path(glue("{per_setting_results_path}/names{num_sur_str}/n{n_str}/alpha_{alpha_init_option}"))
  fs::dir_create(setting_path)

  # if (save_plots) {
  #   pdf(fs::path(glue("{setting_path}/c_chain_plot.pdf"))) # not a replicate-specific name
  #   plot_c_chain(out, burnin = b)
  #   dev.off()
  #
  #   pdf(fs::path(glue("{setting_path}/alpha_chain_plot_9surnames.pdf")))
  #   plot_a_chain(out, indices = 1:9, burnin = b)
  #   dev.off()
  # }

  # Save raw output if requested
  if (save_raw_output) {
    formatted_time <- format(Sys.time(), "%Y%m%d%H%M%S")
    saveRDS(out, fs::path(glue("{setting_path}/{formatted_time}_raw_output.RDS")))
  }
}

#------------------------------------------------------------------------
# Step 3: Define grid of values to loop over
#------------------------------------------------------------------------

# Consolidate files saved in tempfile folder for settings already run
# (not saving anything for now)
# Expect these to already have a temptfile_path column
df_already_run <- merge_tf_csvs(tempfile_path,
  pattern = "^tf_.*\\.csv$",
  output_file = NULL
)


df_settings <- filename_based_setting_updater(
  setting_list = setting_list,
  num2run = num2run,
  df_already_run = df_already_run,
  prefix = "tf_"
)


# # Update settings - old approach not robust to parallelizing stuff
# df_settings <- setting_updater(setting_list,
#   num2run = num2run,
#   existing_df = df_settings_already_run,
#   path_to_existing_df = grid_filepath,
#   run_grouping = 1, # 1 row is createdper setting run
#   return_num_runs = T
# )$settings
# df_settings$run_ID <- NULL # don't use this
#
# # Expand so that a row per setting to run, repeats for num2run >1
# df_expanded <- df_settings %>%
#   uncount(num2run)
# df_expanded <- df_expanded[order(df_expanded$alpha_init_option), ]

# # Create run ID that is unique for expanded version
# df_expanded <- df_expanded %>%
#   group_by(num_sur, n, alpha_init_option) %>%
#   mutate(row_num = row_number())
#
# df_expanded$run_ID <- df_expanded$num_runs + df_expanded$row_num
# df_expanded$row_num <- NULL


# This is to determine what kind of parallel processing to do
# depending on the system used
get_strategy <- function() {
  if (parallelly::supportsMulticore()) {
    return("future::multicore")
  } else {
    return("future::multisession")
  }
}

use_parallel <- F #TODO: haven't fully gotten the parallel to work

if (nrow(df_settings) > 0) {
  if (use_parallel) {
    # Set-up parallelization
    future::plan(
      strategy = get_strategy(),
      workers = 8
    )
  }

  # Now it will run the following with parallelization (and any other use of future package unless re-set it)
  future_apply(df_settings, 1, function(settings) {
    single_run(
      settings,
      per_setting_results_path = per_setting_results_path,
      true_val_path = true_val_path,
      df_states = df_states,
      niter_cutoffs = niter_cutoffs
    )
  },
  future.seed = T
  )

  # reset back so doesn't interfere
  future::plan(strategy = "future::sequential")
} else {
  print("No settings left to run")
}
