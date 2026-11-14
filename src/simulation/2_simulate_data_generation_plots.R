#
# Author: Kyla Chasalow
#
# Small run of simulation for creating illustrative plots that
# appear in paper
#

library(readr) # for loading data
library(gtools) # for rdirichlet
library(tidyverse)
library(glue)
library(ggplot2)
library(patchwork)
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
pR1 = readRDS("../../data/DataObjects/pR1.RDS")
p_state_gvn_R1 = readRDS("../../data/DataObjects/p_state_gvn_R1.RDS")
p_state_gvn_R0 = readRDS("../../data/DataObjects/p_state_gvn_R0.RDS")
p_R1_gvn_state = readRDS("../../data/DataObjects/p_R1_gvn_state.RDS")
p_R0_gvn_state = readRDS("../../data/DataObjects/p_R0_gvn_state.RDS")


# Surname Settings - toggle this to reproduce different simulations in the paper
n_names <- 500 #
surnames <- readRDS(paste0("simulated_data/surnames10000.RDS"))
surnames = surnames[1:n_names]


# Surname distribution for R=1 group
pmat = get_pmat(df_states$Abbreviation, surnames, rate = 1)
pmarg = get_p_surname_gvn_R1(pmat, p_state_gvn_R1)

# Size settings
sf_size = 1200000 #target sampling frame size
n1 = ceiling(sf_size/(1/.02+1)) #number of rare population members



################################################################################
# Simulate Obituary Data
################################################################################
n <- 100000
df <- create_df(pmat, n = n,  props = p_state_gvn_R1, fips_mode = F)

# Surname prob plot
h <- hist(pmarg, breaks = 25, plot = FALSE)
hist_df <- data.frame(
  xmin = head(h$breaks, -1),
  xmax = tail(h$breaks, -1),
  density = h$density
)
i <- 9
j <- 10

scatter_df <- data.frame(
  x = pmat[, i],
  y = pmat[, j]
)

# Approximate base R's automatic axis ticks, including its 4% padding
base_breaks <- function(x) {
  limits <- range(x, finite = TRUE)
  if (diff(limits) == 0) return(pretty(limits))
  limits <- limits + c(-1, 1) * 0.04 * diff(limits)
  grDevices::axisTicks(limits, log = FALSE)
}


base_theme <- theme_classic(base_size = 18, base_family = "sans") +
  theme(
    text = element_text(colour = "black"),

    plot.title = element_text(
      size = 26,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 16)
    ),
    plot.title.position = "panel",

    plot.tag = element_text(
      size = 26,
      face = "bold",
      hjust = 0,
      vjust = 1
    ),
    plot.tag.position = c(0, 1),
    axis.title = element_text(size = 22),
    axis.title.x = element_text(margin = margin(t = 14)),
    axis.title.y = element_text(margin = margin(r = 14)),
    axis.text = element_text(size = 18, colour = "black"),
    axis.text.x = element_text(margin = margin(t = 6)),
    axis.text.y = element_text(
      angle = 0,
      margin = margin(r = 6)
    ),

    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    axis.ticks.length = grid::unit(4, "pt"),

    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    panel.ontop = FALSE,

    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.key = element_blank(),

    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(t = 16, r = 20, b = 16, l = 16)
  )

p_a <- ggplot(hist_df) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = density),
    fill = "black", colour = "white", linewidth = 0.4
  ) +
  scale_x_continuous(
    limits = range(h$breaks),
    breaks = base_breaks(h$breaks),
    expand = expansion(mult = 0.04)
  ) +
  scale_y_continuous(
    limits = c(0, max(h$density)),
    breaks = base_breaks(c(0, h$density)),
    expand = expansion(mult = 0.04)
  ) +
  labs(
    title = "Distribution of Surname Probabilities",
    x = "Surname probabilities, marginalized over state",
    y = "Density",
    tag = "A"
  ) +
  base_theme

p_b <- ggplot(scatter_df, aes(x, y)) +
  geom_point(shape = 16, size = 2.5, alpha = 0.6) +
  geom_abline(
    intercept = 0, slope = 1,
    linetype = "dashed", linewidth = 0.45, colour = "gray40"
  ) +
  scale_x_continuous(
    breaks = base_breaks(scatter_df$x),
    expand = expansion(mult = 0.04)
  ) +
  scale_y_continuous(
    breaks = base_breaks(scatter_df$y),
    expand = expansion(mult = 0.04)
  ) +
  # Prevent the reference line from affecting the plotted ranges
  coord_cartesian(
    xlim = range(scatter_df$x, finite = TRUE),
    ylim = range(scatter_df$y, finite = TRUE)
  ) +
  labs(
    title = "Comparing Surname Probabilities of Two States",
    x = sprintf("Surname probabilities for %s", colnames(pmat)[i]),
    y = sprintf("Surname probabilities for %s", colnames(pmat)[j]),
    tag = "B"
  ) +
  base_theme

surname_plot <- p_a + p_b + plot_layout(ncol = 2, widths = c(1, 1))

fig_save_path <- "figures"
fs::dir_create(fig_save_path)

ggsave(
  filename = file.path(fig_save_path, "surname_probabilities.pdf"),
  plot = surname_plot,
  device = "pdf",
  width = 20,
  height = 10,
  units = "in",
  bg = "white"
)


#################################################################################
# SAMPLING FRAME - single version
#
# Generate a fake sampling frame situation with ground truth labels
# For simulating sampling purposes. Make rare population realistically rare

#################################################################################
beta = 1000
pmat_flip = flip_probabilities(pmat, sd = 1, beta = beta)

sampling_frame = sampling_frame = generate_sampling_frame(pmat, pmat_flip,
                                                          n1, pR1, p_state_gvn_R1,
                                                          p_R1_gvn_state, p_R0_gvn_state)

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
    title = "A. Comparison of Surname Frequencies"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 26,
      margin = margin(b = 15)
    ),
    axis.title.x = element_text(size = 20, margin = margin(t = 12)),
    axis.title.y = element_text(size = 20, margin = margin(r = 12)),
    axis.text = element_text(size = 16, color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.4),
    plot.margin = margin(15, 20, 15, 15),
  )


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
  width = 12,       # inches
  height = 10,      # inches
  device = cairo_pdf,  # high-quality vector text rendering
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
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "B. Probability of Being in Rare Population by Rare Population Status",
    x = "Probability of being in rare population given surname and state",
    y = "Count"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 22,
                              hjust = 0.5,
                              margin = margin(b = 15)),
    strip.text = element_text(face = "bold", size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.4),
    axis.title.x = element_text(size = 20, margin = margin(t = 12)),
    axis.title.y = element_text(size = 20, margin = margin(r = 12)),
    axis.text.x  = element_text(size = 15),
    axis.text.y  = element_text(size = 15),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.6),  # light gray borders
    plot.margin = margin(15, 20, 15, 15),
  )
ggsave(
  filename = file.path(fig_save_path, "R1_probabilities_R0_vs_R1.pdf"),
  plot = q,
  width = 15,       # inches
  height = 8,      # inches
  device = cairo_pdf,  # high-quality vector text rendering
  dpi = 600
)
