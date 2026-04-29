library(tidyverse)
library(glue)
library(ggplot2)
library(ggExtra)
library(scales)
library(forcats)

set.seed(543)
source("sampling_tools.R")
source("data_generators.R")

#----------------------
# Set-up
#----------------------

# AJPP info
pR1 = .02 # rounded version of actual
p_R1_gvn_state = readRDS("../../data/AJPP_objects/p_R1_gvn_state.RDS")
p_R0_gvn_state = readRDS("../../data/AJPP_objects/p_R0_gvn_state.RDS")
p_state_gvn_R1 = readRDS("../../data/AJPP_objects/p_state_gvn_R1.RDS")
p_state_gvn_R0 = readRDS("../../data/AJPP_objects/p_state_gvn_R0.RDS")

# State information
df_states <- read_csv("../../data/states.csv")
df_states = df_states[order(df_states$Abbreviation),] #make sure alphabetical order by abbreviation, not state name!

# Surname info
surnames <- readRDS(paste0("../../data/simulated/surnames10000.RDS"))
n_names = 500
surnames = surnames[1:n_names]

# Target sample size and sampling frame sizes to mimic true rate of
# rare population in our application
Target = 1000
sf_size = 1200000
n1 = ceiling(sf_size/(1/.02+1)) #number of rare population members

# Set-up Simulation
nsim = 100
betavals = c(10, 100, 1000, 10000)


get_sampling_prob_suite <- function(sampling_frame, p_R1_gvn_state, Target){

  # Get per-state targets
  #-----------------------
  group_targets = get_per_state_targets_pois(Target,
                                              sampling_frame,
                                              p_R1_gvn_state)


  # Baseline - Random sample, still with targets from each state
  #---------------------------------------------------------
  sampling_frame$random_sample_per_state = 1
  sampling_frame = get_rescaled_sampling_probs(sampling_frame,
                                              group_column = "state",
                                              group_targets = group_targets,
                                              probability_column =  "random_sample_per_state",
                                              new_prob_column_name = "random_sample_per_state",
                                              fix_greater_than_1_probs = T
                                              )

  # create scaled sampling prob based on real Pr(R=1|s,g)
  #--------------------------------------------------------
  sampling_frame = get_rescaled_sampling_probs(sampling_frame,
                                               group_column = "state",
                                               group_targets = group_targets,
                                               probability_column = "p_R1_gvn_surname_state",
                                               new_prob_column_name = "targeted_prob",
                                               fix_greater_than_1_probs = T
                                               )


  prob_options = c("random_sample_per_state",
                   "targeted_prob")

  nice_labels = c("Random sample in each state",
                  "Surname-state probabilities")

  run_checks(sampling_frame, Target,
             prob_columns = prob_options,
             tol = 10
  )

  return(list("sampling_frame" = sampling_frame,
              "prob_options" = prob_options,
              "nice_labels" = nice_labels))

}

# Note that with fix_greater_than_1_probs = T
# I did an ad hoc fix of setting prob > 1 to 1
# This only happening for the very large beta...
# TODO: come back to this


holder = data.frame()

# Repeated sampling simulation
for (i in 1:nsim){
  print(paste0("Simulation: ", i))
  for(j in seq_along(betavals)){
    print(paste0(j,"th Beta value"))
    beta = betavals[j]

    # Generate sampling frame
    pmat = get_pmat(df_states$Abbreviation, surnames, rate = 1) # P(S|state,R=1)  values
    pmat_flip = flip_probabilities(pmat, sd = 1, beta = beta) # P(S|state,R=0) values

    sampling_frame = generate_sampling_frame(pmat, pmat_flip,
                                             n1, pR1, p_state_gvn_R1,
                                             p_R1_gvn_state, p_R0_gvn_state)

    out = get_sampling_prob_suite(sampling_frame, p_R1_gvn_state, Target)
    sampling_frame = out$sampling_frame
    prob_options = out$prob_options
    nice_labels = out$nice_labels

    temp_holder = rep(NA, length(prob_options))

    #for each sampling method
    for(k in seq_along(prob_options)){

      sample = rbinom(dim(sampling_frame)[1], 1, sampling_frame[[prob_options[k]]])

      nsampled = sum(sample)
      nraresampled = sum(sample == 1 & sampling_frame$rare_status == 1)

      yield = nraresampled/nsampled
      temp_holder[k] = yield
    }

    if (nrow(holder) == 0) {
      cols <- c("nsim", "beta", nice_labels)
      holder <- data.frame(matrix(ncol = length(cols), nrow = 0))
      colnames(holder) <- cols
    }

    # Make temp_holder a *data frame* with the same columns
    temp_holder <- data.frame(t(c(i, beta, temp_holder)))  # transpose to make it a 1-row df
    colnames(temp_holder) <- colnames(holder)              # assign matching column names

    # Bind
    holder <- rbind(holder, temp_holder)


  }
}

# Long/tidy format: one column per method becomes rows
holder_long <- holder %>%
  select(-nsim) %>%                               # ignore nsim
  pivot_longer(-beta, names_to = "method", values_to = "yield") %>%
  mutate(
    yield  = as.numeric(yield)
  )


holder_long$single_method = paste0(holder_long$method, ", β=", holder_long$beta)
holder_long$single_method[grepl("Random sample in each state", holder_long$single_method)] = "Random sample in each state"


p <- ggplot(holder_long, aes(x = yield, y = single_method)) +
  geom_boxplot(outlier.alpha = 0.6, width = 0.6) +
    labs(
      title = glue("Comparison of Success Rates in Sampling Rare Population"),
      subtitle = glue("{nsim} repeat simulations of sampling frame and sampling method"),
      x = "Proportion of Rare Population Sampled",
      y = NULL
    ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", size = 18, hjust = 0.5),
    strip.text   = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  ) + theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
      axis.title.x = element_text(size = 15, face = "bold"),
      axis.text.x  = element_text(size = 15),
      axis.text.y  = element_text(size = 15),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.2)
    )
print(p)

fig_save_path = "../../figures/simulation"
ggsave(
  filename = file.path(fig_save_path, "sampling_yield.pdf"),
  plot = p,
  width = 12,       # inches
  height = 6,      # inches
  device = cairo_pdf,  # high-quality vector text rendering
  dpi = 600
)
