#
# Author: Kyla Chasalow
#
# Analyze files created by 5_ script
#
#
# Create per-setting plots
# This can take some time to run.
#

library(tidyverse)
library(glue)
library(bayesplot)
library(coda) #gelman-rubin diagnostics

source("simulation_tools.R")

save_path = "results/"
true_val_path = fs::path(glue("{save_path}/true_probabilities"))
per_setting_results_path = fs::path(glue("{save_path}/per_setting_results"))
grid_filepath = fs::path(glue("{save_path}/grid_results.csv"))

#df = read.csv(grid_filepath)

tempfile_path = fs::path(glue("{save_path}/tempfiles"))
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

  # keep only dirs that actually have raw_output.RDS
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
single_analysis <- function(files_to_process,
                            save_path,
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


  # Same process with the rho mat if exists (added this to code later so some
  # earlier outputs might not have it)
  if("rho_mat" %in% names(output_list[[1]])) {
    rho_mat_list <- lapply(output_list, function(x){x$rho_mat})
    rho_chains <- abind::abind(rho_mat_list, along = 3)
    dimnames(rho_chains) <- list("Iteration" = NULL,
                                 "State" = colnames(output_list[[1]]$rho_mat),
                                 "Chain" = 1)
    rho_chains = aperm(rho_chains, perm = c(1,3,2))  # Reshape to fit required order
  } else{
    rho_chains = NULL
  }

  bayesplot::color_scheme_set("mix-blue-pink")


  # Plot of density and trace plots across repeat chains - ideally see similar
  rowname_subset = rownames(output_list[[1]]$a_chain)[sample(1:100,12)]
  if(create_density_plots){
    # Density plots
    p = mcmc_dens_overlay(a_chains, pars = rowname_subset)
    p <- theme_for_all_plots(p, 14)

    ggplot2::ggsave(filename = fs::path(glue("{save_path}/a_chains_density_plot.pdf")),
                    plot = p,
                    width = 15,
                    height = 8,
                    device = "pdf")

    p = mcmc_dens_overlay(eta_chains, pars = "eta_scale")
    p <- theme_for_all_plots(p,14)
    ggplot2::ggsave(filename = fs::path(glue("{save_path}/eta_chains_density_plot.pdf")),
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
    ggplot2::ggsave(filename = fs::path(glue("{save_path}/a_chains_trace_plot.pdf")),
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
    ggplot2::ggsave(filename = fs::path(glue("{save_path}/eta_chains_trace_plot.pdf")),
                    plot = p,
                    width = 10,
                    height = 4,
                    device = "pdf")

    if(!is.null(rho_chains)){
      state_subset = c("NY", "FL","TX","WY")
      p = mcmc_trace(rho_chains, pars = state_subset)
      p <- theme_for_all_plots(p, 14) +
        labs(x = "Iteration",
             y = expression(rho[g] ~ "parameter value")
        ) +
        theme(
          axis.title.x = element_text(size = 20, margin = margin(t = 8)),
          axis.title.y = element_text(size = 20, margin = margin(t = 8))
        )
      ggplot2::ggsave(filename = fs::path(glue("{save_path}/rho_chains_trace_plot.pdf")),
                      plot = p,
                      width = 8,
                      height = 4,
                      device = "pdf")
    }
  }
}

#----------------------------------------------------------------------------
raw_lists <- get_leaf_raw_outputs(per_setting_results_path)
#----------------------------------------------------------------------------

# Add multi-chain trace plots to each directory
for (save_path in names(raw_lists)){
  print(save_path)
  files_to_process = raw_lists[[save_path]]
  single_analysis(files_to_process = files_to_process,
                  save_path = save_path,
                  create_density_plots = F,
                  create_trace_plots = T)
}
