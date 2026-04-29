#------------------------------------------------------------
# Small run of simulation for creating illustrative plots
# Not saving anything but the plots
#------------------------------------------------------------

library(readr)
library(gtools) # for rdirichlet
library(tidyverse)
library(glue)
library(ggplot2)
library(ggExtra)
library(scales)

set.seed(21)

source("data_generators.R") #functions for simulating data


################################################################################
# Set-up
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

# Surname Settings - toggle this to reproduce different simulations in the paper
surnames <- readRDS(paste0("../../data/simulated/surnames10000.RDS"))
n_names <- 500 # can set up to 10000
surnames = surnames[1:n_names]

# surname distribution for R=1 group
pmat = get_pmat(df_states$Abbreviation, surnames, rate = 1)
pmarg = get_p_surname_gvn_R1(pmat, p_state_gvn_R1)

# sampling frame size settings
sf_size = 1200000 #target sampling frame size
n1 = ceiling(sf_size/(1/.02+1)) #number of rare population members

# obituary data size settings
n <- 100000

# filepath to save figures
fig_save_path = "../../figures/simulation"


################################################################################
# Simulate Obituary Data
################################################################################

df <- create_df(pmat, n = n,  props = p_state_gvn_R1, fips_mode = F)

# Plot characteristics
#---------------------------
pdf(file.path(fig_save_path, "surname_probabilities.pdf"), width = 14, height = 5)
par(
  mfrow = c(1, 2),
  mar = c(6, 6, 1.5, 1),
  mgp = c(4, 1.0, 0),   # move label slightly farther from axis
  las = 1
)
# Panel (a): Histogram of P(s | R=1)
hist(
  pmarg,
  breaks = 25,
  freq   = FALSE,
  col    = "gray80",
  border = "white",
  main   = "",
  xlab   = "Surname probabilities marginalized over state",
  ylab   = "Density",
  cex.axis = 1.2,
  cex.lab  = 1.3
)
grid(col = "gray90", lty = "dotted")
# Panel (b): Pairwise simulated probabilities
i <- 9; j <- 10
plot(
  pmat[, i],
  pmat[, j],
  pch = 19,
  cex = 1.2,
  col = rgb(0, 0, 0, 0.6),
  xlab = glue("Surname probabilities for {colnames(pmat)[i]}"),
  ylab = glue("Surname probabilities for {colnames(pmat)[j]}"),
  cex.axis = 1.2,
  cex.lab  = 1.3
)
abline(0, 1, lty = "dashed", lwd = 1.2, col = "gray40")
grid(col = "gray90", lty = "dotted")
dev.off()


#################################################################################
# SAMPLING FRAME
# Generate a single fake sampling frame with ground truth labels
# For simulating sampling purposes. Make rare population realistically rare
#################################################################################

# Set how different R=0 group is from R=1 group
beta = 1000

# Get R=0 group P(s|g,R=0) probabilities
pmat_flip = flip_probabilities(pmat, sd = 1, beta = beta)

# Generate frame
sampling_frame = sampling_frame = generate_sampling_frame(pmat,
                                                          pmat_flip,
                                                          n1,
                                                          pR1,
                                                          p_state_gvn_R1,
                                                          p_R1_gvn_state,
                                                          p_R0_gvn_state)

# --- Count and normalize ---
counts <- table(sampling_frame$surname, sampling_frame$rare_status)
probs  <- prop.table(counts, margin = 2)  # normalize within each rare_status group

plotdf <- data.frame(
  surname = rownames(counts),
  prob_0  = probs[, "0"],
  prob_1  = probs[, "1"]
)

# Shared upper bound (2% padding)
upper <- max(plotdf$prob_1, plotdf$prob_0) * 1.02

# Base scatter
p_scatter <- ggplot(plotdf, aes(x = prob_1, y = prob_0, label = surname)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 2.4, alpha = 0.9) +
  coord_equal(expand = FALSE) +
  scale_x_continuous(limits = c(-0.002, upper),
                     labels = number_format(accuracy = 0.001)) +
  scale_y_continuous(limits = c(-0.002, upper),
                     labels = number_format(accuracy = 0.001)) +
  labs(
    x = "Fraction of R = 1 observations with surname s",
    y = "Fraction of R = 0 observations with surname s",
    title = "Comparison of Surname Observation Frequencies"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 16))

# Add marginal histograms (top & right)
p_with_marginals <- ggMarginal(
  p_scatter,
  type = "histogram",
  margins = "both",          # top & right
  bins = 25,
  fill = "gray75",
  color = "white",
  size = 6,                  # relative size of marginals vs main plot
  xparams = list(alpha = 0.9),
  yparams = list(alpha = 0.9)
)

print(p_with_marginals)

ggsave(
  filename = file.path(fig_save_path, "surname_probabilities_R0_vs_R1.pdf"),
  plot = p_with_marginals,
  width = 8,       # inches
  height = 8,      # inches
  device = cairo_pdf,
  dpi = 600
)


#################################################################################
# SAMPLING SIMULATION - varying beta
#################################################################################

betavals = c(10, 100, 1000, 10000)
frames = list()

for (beta in betavals){

  pmat_flip = flip_probabilities(pmat, sd = 1, beta = beta)

  sampling_frame = generate_sampling_frame(pmat, pmat_flip,
                                           n1, pR1, p_state_gvn_R1,
                                           p_R1_gvn_state, p_R0_gvn_state)

  sampling_frame$beta = beta
  frames[[beta]] = sampling_frame

}
sampling_frame <- do.call(rbind, frames)

# Before plotting, convert beta column to character labels formatted as expressions
sampling_frame$beta <- factor(sampling_frame$beta,
                              levels = betavals,
                              labels = paste0("beta==", betavals))

q <- ggplot(sampling_frame, aes(x = p_R1_gvn_surname_state)) +
  geom_histogram(
    bins = 45,
    color = "white",
    fill = "black",
    alpha = 0.9
  ) +
  facet_grid(
    rows = vars(rare_status),
    cols = vars(beta),
    scales = "free_y",
    labeller = labeller(
      rare_status = c(`0` = "Non-Rare Population", `1` = "Rare Population"),
      beta = label_parsed  # parse math expressions for beta labels
    )
  ) +
  labs(
    title = "Probability of Being in Rare Population by Rare Population Status",
    x = "Probability of being in rare population given surname and state",
    y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 18),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.text.x  = element_text(size = 15),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.y  = element_text(size = 15),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.6)  # light gray borders
  )

ggsave(
  filename = file.path(fig_save_path, "R1_probabilities_R0_vs_R1.pdf"),
  plot = q,
  width = 15,       # inches
  height = 8,      # inches
  device = cairo_pdf,  # high-quality vector text rendering
  dpi = 600
)
