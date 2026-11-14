#'
#' Author: Kyla Chasalow
#'
#' AI use: minimal, only some help with tweaking plot code
#'
#' Must set working directory to current script for this to run properly


library(tidyverse)
library(glue)
library(ggplot2)
library(ggExtra)
library(scales)
library(forcats)


set.seed(543)  # for reproducibility

source("sampling_tools.R")
source("data_generators.R")

# AJPP info
p_R1_gvn_state = readRDS("../../data/DataObjects/p_R1_gvn_state.RDS")
p_R0_gvn_state = readRDS("../../data/DataObjects/p_R0_gvn_state.RDS")
p_state_gvn_R1 = readRDS("../../data/DataObjects/p_state_gvn_R1.RDS")
p_state_gvn_R0 = readRDS("../../data/DataObjects/p_state_gvn_R0.RDS")
pR1 = .02

# State information
df_states <- read_csv("../../data/states.csv")
df_states = df_states[order(df_states$Abbreviation),] #make sure alphabetical order by abbreviation, not state name!

# Surname info
n_names = 500
surnames <- readRDS(paste0("simulated_data/surnames10000.RDS"))
surnames = surnames[1:n_names]

Target = 1000 #target sample size
sf_size = 1200000
n1 = ceiling(sf_size * .02)

# Set-up Simulation
nsim = 50
betavals = c(1, 10, 50, 100, 500, 1000, 5000, 10000)


get_sampling_prob_suite <- function(sampling_frame, p_R1_gvn_state, Target){

  # Get per-state targets
  #-----------------------

  # group_targets = get_per_state_targets_strat(Target,
  #                                       sampling_frame,
  #                                       p_R1_gvn_state)

  group_targets = get_per_state_targets_pois(Target,
                                              sampling_frame,
                                              p_R1_gvn_state)

  # Baselines
  #---------------------------------------------------------
  # Random sample, still with targets from each state
  sampling_frame$random_sample_per_state = 1
  sampling_frame = get_rescaled_sampling_probs(sampling_frame,
                                              group_column = "state",
                                              group_targets = group_targets,
                                              probability_column =  "random_sample_per_state",
                                              new_prob_column_name = "random_sample_per_state",
                                              fix_greater_than_1_probs = T #TODO: think about what to do with this
                                              )

  # create scaled sampling prob based on real Pr(R=1|s,g)
  #--------------------------------------------------------
  sampling_frame = get_rescaled_sampling_probs(
     sampling_frame,
     group_column = "state",
     group_targets = group_targets,
     probability_column = "p_R1_gvn_surname_state",
     new_prob_column_name = "targeted_prob",
     fix_greater_than_1_probs = T #TODO: think about what to do with this
     )


  prob_options = c("random_sample_per_state",
                   "targeted_prob")
                   #qcolnames
                   #prec_colnames)

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

#TODO: Note that I did a temp fix of setting prob > 1 to 1
#This only happening for the very large beta...

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

    temp_holder = rep(NA, length(prob_options) + 2) #+2 for the two BISG

    # Mean BISG prob for each group
    mean_BISG_rare = mean(sampling_frame$p_R1_gvn_surname_state[sampling_frame$rare_status==1])
    mean_BISG_notrare = mean(sampling_frame$p_R1_gvn_surname_state[sampling_frame$rare_status==0])
    temp_holder[length(prob_options)+1] = mean_BISG_rare
    temp_holder[length(prob_options)+2] = mean_BISG_notrare
    BISG_precision_labels = c("mean_BISG_rare", "mean_BISG_notrare")


    #for each sampling method
    for(k in seq_along(prob_options)){

      sample = rbinom(dim(sampling_frame)[1], 1, sampling_frame[[prob_options[k]]])

      nsampled = sum(sample)
      nraresampled = sum(sample == 1 & sampling_frame$rare_status == 1)

      yield = nraresampled/nsampled
      temp_holder[k] = yield
    }

    if (nrow(holder) == 0) {
      cols <- c("nsim", "beta", nice_labels, BISG_precision_labels)
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
  pivot_longer(-c(beta,
                  mean_BISG_rare,
                  mean_BISG_notrare)
               , names_to = "method", values_to = "yield") %>%
  mutate(
    yield  = as.numeric(yield)
  )

holder_long$single_method = paste0(holder_long$method, ", β=", holder_long$beta)
holder_long$single_method[grepl("Random sample in each state", holder_long$single_method)] = "Random sample in each state"
holder_long$dif_in_mean_BISG = holder_long$mean_BISG_rare - holder_long$mean_BISG_notrare


# save simulation output in case want to tweak plots below without re-running
save_path = "results/sampling_yield_results"
fs::dir_create(save_path)
write.csv(holder_long, fs::path(
  save_path, "yield_results.csv"
))


#---------------------------------------------------------------------------
# Old plot from earlier version of paper - plotting against beta is not so
# interpretable
#---------------------------------------------------------------------------
# p <- ggplot(holder_long, aes(x = yield, y = single_method)) +
#   geom_boxplot(outlier.alpha = 0.6, width = 0.6) +
#     labs(
#       title = glue("Comparison of Success Rates in Sampling Rare Population"),
#       subtitle = glue("{nsim} repeat simulations of sampling frame and sampling method"),
#       x = "Proportion of Rare Population Sampled",
#       y = NULL
#     ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title   = element_text(face = "bold", size = 18, hjust = 0.5),
#     strip.text   = element_text(face = "bold", size = 14),
#     panel.grid.minor = element_blank()
#   ) + theme_minimal(base_size = 16) +
#     theme(
#       plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
#       axis.title.x = element_text(size = 15, face = "bold"),
#       axis.text.x  = element_text(size = 15),
#       axis.text.y  = element_text(size = 15),
#       panel.grid.minor = element_blank(),
#       legend.position = "none",
#       panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.2)
#     )
# print(p)
#
# fig_save_path = "figures"
# ggsave(
#   filename = file.path(fig_save_path, "sampling_yield.pdf"),
#   plot = p,
#   width = 12,       # inches
#   height = 6,      # inches
#   device = cairo_pdf,  # high-quality vector text rendering
#   dpi = 600
# )
#---------------------------------------------------------------------------


yield_summaries <- holder_long %>%
  group_by(single_method) %>%
  summarise(
    mean_rare_BISG_value = mean(mean_BISG_rare),
    mean_dif = mean(dif_in_mean_BISG),
    mean_yield = mean(yield),
    rare_lower = quantile(mean_BISG_rare, 0.025), #show middle 95% of simulations
    rare_upper = quantile(mean_BISG_rare, 0.975),
    dif_lower = quantile(dif_in_mean_BISG, 0.025),
    dif_upper = quantile(dif_in_mean_BISG, 0.975),
    yield_lower = quantile(yield, 0.025),
    yield_upper = quantile(yield, 0.975),
    n_sim = n(),
    .groups = "drop"
  )

# separate the random sampling version -- BISG dif doesn't affect these
random_sample_mean = yield_summaries %>%
  filter(single_method == "Random sample in each state") %>%
  pull(mean_yield)

random_sample_lower = yield_summaries %>%
  filter(single_method == "Random sample in each state") %>%
  pull(yield_lower)

random_sample_upper = yield_summaries %>%
  filter(single_method == "Random sample in each state") %>%
  pull(yield_upper)

yield_summaries = yield_summaries %>% filter(
  single_method != "Random sample in each state"
)

# Main plot
q <- ggplot(
  yield_summaries,
  aes(x = mean_dif, y = mean_yield)
) +
  # Vertical confidence interval for mean yield
  geom_errorbar(
    aes(ymin = yield_lower, ymax = yield_upper),
    width = 0,
    show.legend = FALSE
  ) +
  # Horizontal confidence interval for mean difference
  geom_segment(
    aes(x = dif_lower, xend = dif_upper, yend = mean_yield),
    show.legend = FALSE
  ) +
  geom_point(size = 4) +
  labs(
    title = "BISG Probability Precision vs Rare Population Sampling Success Rates",
    #subtitle = glue("Means across {nsim} repeat simulations of sampling frame and sampling method"),
    x = "Mean difference in mean BISG probabilities between rare and non-rare group",
    y = "Mean proportion of sample that is rare population",
  ) +
  coord_cartesian(xlim = c(0, 1.05),
                  ylim = c(0, 1.05), expand = FALSE) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 26, hjust = 0.5),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text = element_text(size = 22),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "gray", fill = NA, linewidth = 0.2
    ),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  geom_hline(
    yintercept = random_sample_mean,
    colour = "blue",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = 0.40,
    y = random_sample_mean + 0.04,
    label = "Mean yield of taking a random sample in each state",
    colour = "blue",
    hjust = 0,
    size = 8
  )

fig_save_path = "figures"
ggsave(
  filename = file.path(fig_save_path, "sampling_yield2.pdf"),
  plot = q,
  width = 16,       # inches
  height = 10,      # inches
  device = cairo_pdf,  # high-quality vector text rendering
  dpi = 600
)



