#'
#' Create per-setting plots
#' This can take some time to run
#'

library(tidyverse)
library(glue)
library(bayesplot)
library(coda) #gelman-rubin diagnostics

source("simulation_tools.R")

# Load simulation results
results_savepath = "../../results/simulation"
true_val_path = fs::path(glue("{results_savepath}/true_probabilities"))
per_setting_results_path = fs::path(glue("{results_savepath}/per_setting_results"))
grid_filepath = fs::path(glue("{results_savepath}/grid_results.csv"))
tempfile_path = fs::path(glue("{results_savepath}/tempfiles"))
df = merge_tf_csvs(tempfile_path, pattern = "^tf_.*\\.csv$")


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Processing raw outputs
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Get the raw output files grouped by their directory
get_leaf_raw_outputs <- function(base, pattern = "raw_output\\.RDS$") {

  # all dirs including base
  dirs <- c(base, list.dirs(base, recursive = TRUE, full.names = TRUE))

  # for each dir, list raw_output.RDS (if any)
  raw_by_dir <- lapply(dirs, function(d) {
    list.files(d, pattern = pattern, full.names = TRUE)
  })

  # keep only dirs that actually have files ending in raw_output.RDS
  keep <- lengths(raw_by_dir) > 0
  raw_by_dir <- raw_by_dir[keep]
  names(raw_by_dir) <- dirs[keep]  # name each element by its leaf directory path

  raw_by_dir
}

#----------------------------------------------------------------------------
theme_for_all_plots <- function(p, base_size = 14, grid = F) {
  p +
    theme_bw(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size + 2),
      axis.text  = element_text(size = base_size),
      strip.text = element_text(size = base_size),
      axis.ticks.length = unit(0.25, "cm"),
      axis.ticks = element_line(linewidth = 0.8),
      if(!grid){
        panel.grid = element_blank()
      }
    )
}

# Analysis for a single setting
single_analysis <- function(files_to_process, results_savepath,
                            create_density_plots = F,
                            create_trace_plots = T){

  #-------------------------------------------------------------------------------
  # Load all the files and put them in a 3-D array
  #-------------------------------------------------------------------------------
  output_list <- lapply(files_to_process, readRDS)

  # Get rid of old name
  for(i in 1:length(output_list)){
    lst = output_list[[i]]
    names(lst)[names(lst) == "c_chain"] <- "eta_chain"
    names(lst)[names(lst) == "num_accept_c"] <- "num_accept_eta"
    output_list[[i]] = lst
  }

  n_chains = length(output_list)
  n_iter = length(output_list[[1]]$eta_chain) #TODO: should this be + the init value? already have form loop?


  #-------------------------------------------------------------------------------
  # Create Nice MCMC Plots using bayesplot
  #-------------------------------------------------------------------------------

  # Get collections of chains
  eta_chain_list <- lapply(output_list, function(x){x$eta_chain})
  a_chain_list <- lapply(output_list, function(x){x$a_chain})

  # Turn them into single arrays
  eta_chains <- abind::abind(eta_chain_list, along = 2)
  eta_chains = array(eta_chains, dim = c(dim(eta_chains),1)) #make it 3D array
  a_chains <- abind::abind(a_chain_list, along = 3)

  # Reformat in order expected by bayesplot library
  dimnames(eta_chains) <- list("Iteration" = NULL,
                             "Chain" = 1:n_chains,
                             "Parameter" = "eta_scale")
  dimnames(a_chains) <- list("Parameter" = rownames(output_list[[1]]$a_chain),
                             "Iteration" = NULL,
                             "Chain" = 1:n_chains)
  a_chains = aperm(a_chains, perm = c(2,3,1))  # Reshape to fit required order


  bayesplot::color_scheme_set("mix-blue-pink")


  # Plot of density and trace plots across repeat chains - ideally see similar
  rowname_subset = rownames(output_list[[1]]$a_chain)[sample(1:100,12)]
  if(create_density_plots){
    # Density plots
    p = mcmc_dens_overlay(a_chains, pars = rowname_subset)
    p <- theme_for_all_plots(p, 14)

    ggplot2::ggsave(filename = fs::path(glue("{results_savepath}/a_chains_density_plot.pdf")),
                    plot = p,
                    width = 15,
                    height = 8,
                    device = "pdf")

    p = mcmc_dens_overlay(eta_chains, pars = "eta_scale")
    p <- theme_for_all_plots(p,14)
    ggplot2::ggsave(filename = fs::path(glue("{results_savepath}/eta_chains_density_plot.pdf")),
                    plot = p,
                    width = 8,
                    height = 4,
                    device = "pdf")
  }

  if(create_trace_plots){
    # Trace Plots
    p = mcmc_trace(a_chains, pars = rowname_subset)
    p <- theme_for_all_plots(p, 14) +
      labs(
        x = "Iteration",
        y = expression(alpha[s] ~ "parameter value")
      ) +
      theme(
        axis.title.x = element_text(size = 20, margin = margin(t = 8)),
        axis.title.y = element_text(size = 20, margin = margin(t = 8))
      )
    ggplot2::ggsave(filename = fs::path(glue("{results_savepath}/a_chains_trace_plot.pdf")),
                    plot = p,
                    width = 15,
                    height = 8,
                    device = "pdf")
    p = mcmc_trace(eta_chains, pars = "eta_scale")
    p <- theme_for_all_plots(p, 14) +
      labs(x = "Iteration",
           y =expression(eta ~ "parameter value")
           ) +
      theme(
        axis.title.x = element_text(size = 20, margin = margin(t = 8)),
        axis.title.y = element_text(size = 20, margin = margin(t = 8))
      )
    ggplot2::ggsave(filename = fs::path(glue("{results_savepath}/eta_chains_trace_plot.pdf")),
                    plot = p,
                    width = 10,
                    height = 4,
                    device = "pdf")
  }
}


#----------------------------------------------------------------------------
raw_lists <- get_leaf_raw_outputs(per_setting_results_path)
#----------------------------------------------------------------------------

# Add multi-chain trace plots to each directory
for (results_savepath in names(raw_lists)){
  print(results_savepath)
  files_to_process = raw_lists[[results_savepath]]
  single_analysis(files_to_process = files_to_process,
                  results_savepath = results_savepath,
                  create_density_plots = F,
                  create_trace_plots = T)
}


#-------------------------------------------------------------------------------
# GELMAN-RUBIN Diagnostics
#-------------------------------------------------------------------------------
#TODO: revisit this.
# code below continues from the code that is now in single_analysis() above
# but may be out of date.
#
# # Convert to format needed for Gelman-Rubin
# a_mcmc_list = vector("list", n_chains)
# eta_mcmc_list = vector("list", n_chains)
# for (i in 1:n_chains){
#   a_mcmc_list[[i]] = mcmc(a_chains[,i,]) #each object of list is an mcmc object
#   eta_mcmc_list[[i]] = mcmc(eta_chains[,i,])
# }
# a_mcmc_list = mcmc.list(a_mcmc_list)
# eta_mcmc_list = mcmc.list(eta_mcmc_list)
#
# #gelman rubin diagnostic
# eta_diagnostics = gelman.diag(eta_mcmc_list)
# a_diagnostics = gelman.diag(a_mcmc_list)
# a_diagnostics
# eta_diagnostics
#
# # 100 parameters took .85 seconds -- so 70000 would take 595 seconds...
# # assuming linear increase (which may not be).
# start = Sys.time()
# out = gelman.diag(a_mcmc_list)
# end = Sys.time()
# print(end-start)
#
#
# #calculate diagnostic at different points along chain
# cutoffs = seq(100,n_iter, by = 500)
#
# point.max = rep(NA, length(cutoffs))
# upper.CI.max = rep(NA, length(cutoffs))
# point.mean = rep(NA, length(cutoffs))
# upper.CI.mean = rep(NA, length(cutoffs))
#
# for (i in 1:length(cutoffs)){
#   print(paste0("Cutoff",i))
#   c = cutoffs[i]
#   chainsub = mcmc.list(lapply(a_mcmc_list, function(x){mcmc(x[1:c,])} ))
#   diagnostics = gelman.diag(chainsub, transform = T)#use transform?
#   point.max[i] = max(diagnostics$psrf[,"Point est."])
#   upper.CI.max[i] = max(diagnostics$psrf[,"Upper C.I."])
#   point.mean[i] = mean(diagnostics$psrf[,"Point est."])
#   upper.CI.mean[i] = mean(diagnostics$psrf[,"Upper C.I."])
# }
#
# par(mfrow = c(1,2))
# plot(cutoffs, point.max, "l", ylim = c(1, max(point.max)),
#      main = "Point Estimate", col = "red",
#      ylab = "GR Diagnostic")
# lines(cutoffs, point.mean, col = "blue")
# abline(h=1.1, col = "black", lty = 2)
# abline(h=1, col = "black", lty = 2)
#
# plot(cutoffs, upper.CI.max, "l", ylim = c(1, max(point.max)),
#      main = "Upper CI", col = "red",
#      ylab = "GR Diagnostic")
# lines(cutoffs, upper.CI.mean, col = "blue")
# abline(h=1.1, col = "black", lty = 2)
# abline(h=1, col = "black", lty = 2)
# legend("topright", legend = c("Max","Mean"), col = c("red","blue"), lwd =2)
#


