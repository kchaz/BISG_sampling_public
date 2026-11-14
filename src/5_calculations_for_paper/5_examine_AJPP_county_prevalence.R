
# state distributions for reference
p_state = readRDS("../../data/DataObjects/p_R1_gvn_state.RDS")

# Load the dataset that still has per-county
filename <- "../../data/JewishCountyCounts/brandeis_fips_pop_counts-proc.csv"
AJPP_data <- read.csv(filename)
stopifnot(sum(AJPP_data$county_given_jewish_adult) == 1)

# NOTE: I've called this 'county' but it is technically sometimes a
# collection of counties (multiple fips)
p_county = AJPP_data$jewish_adult_given_county

#-------------------------------------------------------------------------------
# Histogram of probabilities
# NOTE: plot code created with ChatGPT 6.1 Sol Medium and then tweaked
#-------------------------------------------------------------------------------
breaks <- hist(c(p_county, p_state), breaks = 20, plot = FALSE)$breaks

h1 <- hist(p_county, breaks = breaks, plot = FALSE)
h2 <- hist(p_state,  breaks = breaks, plot = FALSE)

# Normalize each histogram to proportions
h1$counts <- h1$counts / sum(h1$counts)
h2$counts <- h2$counts / sum(h2$counts)


draw_histogram <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(
    mar = c(5.5, 5.5, 5.5, 1.5),
    family = "sans",
    cex.lab = 1.3,
    cex.axis = 1.1,
    las = 1,
    bty = "l",
    mgp = c(3.2, 0.8, 0),
    tcl = -0.25,
    lwd = 2
  )

  county_color <- "#0072B2"
  state_color  <- "#D55E00"
  fills <- adjustcolor(c(state_color, county_color), alpha.f = 0.6)

  ymax <- min(1, max(h1$counts, h2$counts) * 1.2)

  plot(
    h2,
    freq = TRUE,
    col = fills[1],
    border = NA,#state_color,
    ylim = c(0, 1), yaxs = "i",
    main = "Distribution of Jewish prevalence by state versus county",
    xlab = "AJPP estimated probability",
    ylab = "Proportion of geographic units",
    axes = FALSE,
    cex.main = 2,
    cex.lab = 2
  )

  plot(
    h1,
    freq = TRUE,
    col = fills[2],
    border = NA,#county_color,
    add = TRUE
  )

  # Draw nothing if count is 0 zero
  draw_bar_outlines <- function(h, color, width = 3) {
    keep <- h$counts > 0

    rect(
      xleft   = head(h$breaks, -1)[keep],
      ybottom = 0,
      xright  = tail(h$breaks, -1)[keep],
      ytop    = h$counts[keep],
      col     = NA,
      border  = color,
      lwd     = width
    )
  }
  # After both histogram plots, before drawing the axes:
  draw_bar_outlines(h2, state_color)
  draw_bar_outlines(h1, county_color)

  # reset lwd before drawing axes
  axis(1, lwd = 0, lwd.ticks = 1)
  axis(2, lwd = 0, lwd.ticks = 1)

  legend(
    "topright",
    legend = c("P(Jewish | State)", "P(Jewish | County)"),
    fill = fills,
    border = c(state_color, county_color),
    bty = "n",
    cex = 2,
    inset = 0.02
  )
}
draw_histogram()

# Save as vector PDF (ideal for publications)
pdf(
  "../../paper/Figures/allocation/county_vs_state_histogram.pdf",
  width = 16, height = 12,
  pointsize = 16,
  useDingbats = FALSE
)
draw_histogram()
dev.off()

summary(AJPP_data$jewish_adult_given_county)

#-----------------------------------------------------------------
# What fraction of Jews live in the counties with top prevalence
#-----------------------------------------------------------------
sum(AJPP_data$jewish_adult_given_county < .02)
sum(AJPP_data$county_given_jewish_adult[AJPP_data$jewish_adult_given_county < .02])
sum(AJPP_data$jewish_adult_given_county > .02 & AJPP_data$jewish_adult_given_county < .05)
sum(AJPP_data$county_given_jewish_adult[AJPP_data$jewish_adult_given_county >= .02 & AJPP_data$jewish_adult_given_county < .05])
sum(AJPP_data$jewish_adult_given_county >= .05)
sum(AJPP_data$county_given_jewish_adult[AJPP_data$jewish_adult_given_county >= .05])

#-----------------------------------------------------------------
# Within-state variation
#-----------------------------------------------------------------
library(dplyr)
library(ggplot2)

# Summarize county prevalence within each state
state_summary <- AJPP_data %>%
  filter(!is.na(state), is.finite(jewish_adult_given_county)) %>%
  group_by(state) %>%
  summarise(
    minimum = min(jewish_adult_given_county),
    maximum = max(jewish_adult_given_county),
    median  = median(jewish_adult_given_county),
    .groups = "drop"
  ) %>%
  mutate(state = reorder(state, median))

# Horizontal range lines, endpoint dots, and median squares
state_plot <- ggplot(state_summary, aes(y = state)) +
  geom_segment(
    aes(x = minimum, xend = maximum, yend = state),
    linewidth = 0.8
  ) +
  geom_point(aes(x = minimum), size = 2.5) +
  geom_point(aes(x = maximum), size = 2.5) +
  geom_point(
    aes(x = median),
    shape = 22, size = 3.5, fill = "white", stroke = 1
  ) +
  scale_x_continuous(labels = scales::label_percent()) +
  labs(
    x = "AJPP estimated probability of Jewish given County",
    y = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    plot.caption = element_text(hjust = 0, size = 12)
  )

print(state_plot)

# Save with enough vertical space for all states
ggsave(
  "../../paper/Figures/allocation/county_prevalence_by_state.pdf",
  plot = state_plot,
  width = 10,
  height = max(6, 0.24 * nrow(state_summary) + 1),
  units = "in"
)

