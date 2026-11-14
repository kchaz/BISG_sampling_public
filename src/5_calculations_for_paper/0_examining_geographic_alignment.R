#
# Author: Kyla Chasalow
#
# AI use: used ChatGPT for help with generating and refining plot code.
#

# Packages
library(dplyr)
library(ggplot2)
library(forcats)


pR1 <- readRDS("../../data/DataObjects/pR1.RDS")
p_state_gvn_R1 <- readRDS("../../data/DataObjects/p_state_gvn_R1.RDS")
p_state_gvn_R0 <- readRDS("../../data/DataObjects/p_state_gvn_R0.RDS")
p_R1_gvn_state <- readRDS("../../data/DataObjects/p_R1_gvn_state.RDS")f
p_R0_gvn_state <- readRDS("../../data/DataObjects/p_R0_gvn_state.RDS")
state_lookup <- readRDS("../../data/DataObjects/fips2state.RDS")
states = read.csv("../../data/states.csv")
p_R1_gvn_state_upper <- readRDS("../../data/DataObjects/p_R1_gvn_state_upper.RDS")
p_R1_gvn_state_lower <- readRDS("../../data/DataObjects/p_R1_gvn_state_lower.RDS")
p_R0_gvn_state_upper <- readRDS("../../data/DataObjects/p_R0_gvn_state_upper.RDS")
p_R0_gvn_state_lower <- readRDS("../../data/DataObjects/p_R0_gvn_state_lower.RDS")



###################################################################################
# Goal: generate a plot showing that the upper vs lower bounds don't lead to
# very different proportions
###################################################################################


# Assumes p_R1_gvn_state, p_R1_gvn_state_lower, p_R1_gvn_state_upper are
# named numeric vectors with state identifiers in names.

# Create df with state ID and estimate columns
df_R1 <- tibble(
  state_id = names(p_R1_gvn_state),
  estimate = as.numeric(p_R1_gvn_state)
) %>%
  mutate(
    lower = as.numeric(p_R1_gvn_state_lower[state_id]),
    upper = as.numeric(p_R1_gvn_state_upper[state_id])
  )

# arrange so that ordered by estimate value
df_R1 <- df_R1 %>%
  arrange(estimate) %>%
  mutate(state_id = fct_inorder(state_id))

# Create plot
ggplot(df_R1, aes(x = state_id, y = estimate)) +
  geom_linerange(aes(ymin = lower, ymax = upper), linewidth = 0.6, alpha = 0.7) +
  geom_point(size = 1.8) +
  coord_flip() +
  labs(
    title = "Proportion of Adults that are Jewish by State",
    x = NULL,
    y = "Estimated Pr(R=1 | g)"
  ) +
  theme_minimal(base_size = 14) + theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.text.x = element_text(size = 14),   # <--- make tick labels larger
    axis.text.y = element_text(size = 14),   # <--- make tick labels larger
    axis.title.x = element_text(size = 15),  # <--- make x-axis title larger
    axis.title.y = element_text(size = 15),  # <--- make x-axis title larger
    axis.ticks.x = element_line(linewidth = 0.8),  # <--- make tick marks thicker
    axis.ticks.y = element_line(linewidth = 0.8),  # <--- make tick marks thicker
  )

# Save the plot to a PDF file
ggsave(
  filename = "../../paper/Figures/data_validation/pR1_state_bounds.pdf",  # save as PDF
  plot = last_plot(),        # saves the most recently created ggplot
  width = 8,                 # width in inches
  height = 12,               # height in inches (taller to fit 51 states)
  units = "in"               # specify inches as the measurement unit
)


#################################################################################
# Goal: compare obit data proportions to Brandeis ones
# and compare voter file proportionsof names from obit data to Brandeis ones
#################################################################################

# Obit data
df = read.csv("../../data/ObitData/obit_list_092424-fulldehyph-check.csv")
counts = table(df$state)

# Voter-file counts for surnames in obit data
df_vf = read.csv("../../data/cluster_output/dist_by_state_of_obit_surnames.csv")


# How prevalence in obit data relates to Brandeis data
#--------------------------------------------------------
# excluded states
excluded = states$Abbreviation[!states$Abbreviation %in% unique(df$state)]

zero_counts = rep(0, length(excluded))
names(zero_counts) = excluded

obit_fracs = c(counts/sum(counts), zero_counts)
obit_fracs = obit_fracs[names(p_state_gvn_R1)]
show_state = p_state_gvn_R1 > .01 | obit_fracs > 0.005
u = max(obit_fracs, p_state_gvn_R1)*1.05
c1 = cor(obit_fracs, p_state_gvn_R1)
print(c1)
c2 = cor(obit_fracs[obit_fracs>0], p_state_gvn_R1[obit_fracs > 0])
print(c2)

# Fraction coming from each state among people with surnames in obit data
vf_fracs = df_vf$Count_SurnameObs_in_Obit/sum(df_vf$Count_SurnameObs_in_Obit)
names(vf_fracs) = df_vf$State
c3 = cor(vf_fracs, p_state_gvn_R1)
print(c3)


# Make sure same order
obit_fracs = obit_fracs[states$Abbreviation]
vf_fracs = vf_fracs[states$Abbreviation]
p_state_gvn_R1 = p_state_gvn_R1[states$Abbreviation]

# ChatGPT generated version of plot above
plotdf <- data.frame(
  state = names(obit_fracs),
  ajpp  = as.numeric(p_state_gvn_R1),
  vf = as.numeric(vf_fracs),
  obit  = as.numeric(obit_fracs),
  row.names = NULL
)
p <- ggplot(plotdf, aes(y = ajpp, x = obit)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 2.6, alpha = 0.9) +
  geom_text(
    data = plotdf[show_state, ],
    aes(label = state),
    size = 3.5,
    hjust = -0.5,   # shift text slightly to the right of the point
    vjust = 0.4,    # small vertical adjustment
    color = "black"
  ) +
  coord_equal(xlim = c(0, u), ylim = c(0, u), expand = FALSE) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  labs(
    title = "Obituary Data vs AJPP Estimate",
    y = "AJPP estimate of P(State | Jewish)",
    x = "Fraction coming from each state\n among observed obituaries",
    tag = "A"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 26, hjust = 0.5),
    axis.text.x = element_text(size = 14),   # <--- make tick labels larger
    axis.text.y = element_text(size = 14),   # <--- make tick labels larger
    axis.title.x = element_text(size = 15),  # <--- make x-axis title larger
    axis.title.y = element_text(size = 15),  # <--- make x-axis title larger
    axis.ticks.x = element_line(linewidth = 0.8),  # <--- make tick marks thicker
    axis.ticks.y = element_line(linewidth = 0.8),  # <--- make tick marks thicker
    plot.tag = element_text(face = "bold", size = 26,
                            hjust = 0, vjust = 1),
    plot.tag.position = c(0, 1),
  )
print(p)

q <- ggplot(plotdf, aes(y = ajpp, x = vf)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 2.6, alpha = 0.9) +
  geom_text(
    data = plotdf[show_state, ],
    aes(label = state),
    size = 3.5,
    hjust = -0.5,   # shift text slightly to the right of the point
    vjust = 0.4,    # small vertical adjustment
    color = "black"
  ) +
  coord_equal(xlim = c(0, u), ylim = c(0, u), expand = FALSE) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  labs(
    title = "Voter File vs AJPP Estimate",
    y = "AJPP estimate of P(State | Jewish)",
    x = "Fraction coming from each state\n among voter file entries with surnames from the obituary data",
    tag = "B"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 26, hjust = 0.5),
    axis.text.x = element_text(size = 14),   # <--- make tick labels larger
    axis.text.y = element_text(size = 14),   # <--- make tick labels larger
    axis.title.x = element_text(size = 15),  # <--- make x-axis title larger
    axis.title.y = element_text(size = 15),  # <--- make x-axis title larger
    axis.ticks.x = element_line(linewidth = 0.8),  # <--- make tick marks thicker
    axis.ticks.y = element_line(linewidth = 0.8),  # <--- make tick marks thicker
    plot.tag = element_text(face = "bold", size = 26,
                            hjust = 0, vjust = 1),
    plot.tag.position = c(0, 1),
  )
print(q)

library(patchwork)

# Combine side by side
combined_plot <- p + q + plot_layout(guides = "collect")
print(combined_plot)
ggsave("../../paper/Figures/data_validation/obits_vs_ajpp_prop.pdf", combined_plot, width = 15, height = 10, dpi = 350)
