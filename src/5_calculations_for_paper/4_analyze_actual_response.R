library(dplyr)
library(scales)
library(ggplot2)
library(patchwork)


################################################################################
# Data on response rates
################################################################################
df = read.csv("../../data/SurveyPrivate/state_table.csv")

# Calculate means and totals
apply(df[,2:6], 2, mean, na.rm = T)
apply(df[,2:6], 2, sum, na.rm = T)

# zeros
df$Responded[df$state == "ND"] = 0
df$Jewish[df$state == "ND"] = 0
df$Response.rate[df$state == "ND"] = 0

df$Jewish.rate[is.na(df$Jewish.rate) & df$state != "ND"] = 0

summary(df$Response.rate)
sd(df$Response.rate, na.rm = T)
df$state[which.min(df$Response.rate)]
df$state[which.max(df$Response.rate)]

max(df$Response.rate[df$Response.rate < .1])
min(df$Responded[df$Responded > 0])
mean(df$Responded)
median(df$Responded)


# Panel A: Response rates by state
p_response <- ggplot(df, aes(x = Response.rate)) +
  geom_histogram(
    bins = 25,
    boundary = 0,
    color = "white",
    fill = "gray20"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "State response rate",
    y = "Number of states",
    title = "Distribution of Response Rates",
    tag = "A"
  )

# Panel B: Jewish share of respondents by state
p_jewish <- ggplot(df, aes(x = Jewish.rate)) +
  geom_histogram(
    bins = 20,
    boundary = 0,
    color = "white",
    fill = "gray20"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "State percentage of respondents who are Jewish",
    y = "Number of states",
    title = "Jewish Share of Respondents",
    tag = "B"
  )

# Stack panels vertically and center their titles
p_combined <- (p_response / p_jewish) &
  theme_classic(base_size = 18) &
  theme(
    plot.title = element_text(
      face = "bold", hjust = 0.5, size = 26,
      margin = margin(b = 12)
    ),
    axis.title = element_text(size = 18),
    axis.title.x = element_text(margin = margin(t = 12)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text = element_text(size = 14, color = "black"),
    plot.margin = margin(12, 20, 12, 15),
    plot.tag = element_text(face = "bold", size = 26,
                            hjust = 0, vjust = 1),
    plot.tag.position = c(0, 1),
  )

p_combined

ggsave(
  filename = "../../paper/Figures/response_analysis/response_rates_combined.pdf",
  plot = p_combined,
  width = 10,
  height = 10
)


################################################################################
# survey data - calibration plot
################################################################################
library(readxl)
df <- read_excel("../../data/SurveyPrivate/survey_for_analysis_100125_with_total_prob.xlsx")
# Note: has +6 rows for the people not linked...those people correspond to 6 of the non-respondent rows

df_obs = df[!is.na(df$jewish),]

# Note: 6 people not linkable
nrow(df_obs)
sum(is.na(df_obs$total_probs_bounded))
nrow(df_obs) - sum(is.na(df_obs$total_probs_bounded))

# Filter out those 6
df_obs = df_obs[!is.na(df_obs$total_probs_bounded),]
df_obs$total_probs_bounded = pmin(df_obs$total_probs_bounded, 1)


get_calibration_data <- function(df, prob_variable, jewish_variable, n_bins = 10) {

  # Create quantile-based bins
  breaks <- quantile(
    df[[prob_variable]],
    probs = seq(0, 1, length.out = n_bins + 1),
    na.rm = TRUE,
    type = 7
  )
  breaks <- unique(breaks)
  # Calculate bin metrics
  df_cal <- df %>%
    mutate(
      prob_bin = cut(
        .data[[prob_variable]],
        breaks = breaks,
        include.lowest = TRUE,
        right = TRUE
      )
    ) %>%
    group_by(prob_bin) %>%
    summarise(
      n = n(),
      mean_prob = mean(.data[[prob_variable]], na.rm = TRUE),
      frac_jewish = mean(.data[[jewish_variable]], na.rm = TRUE),
      .groups = "drop"
    )
  return(df_cal)
}

# Variable 1: the unnormalized probabilties (total_probs_bounded) for each survey respondents
# Variable 2: whether or not that respondent was Jewish (is_jewish)
n_bins = 20
df_cal = get_calibration_data(df_obs,
                  prob_variable = "total_probs_bounded",
                  jewish_variable = "jewish",
                  n_bins = n_bins
                  )
#----------------------------------------------------
# Plot
#----------------------------------------------------
w <- ggplot(df_cal, aes(x = mean_prob, y = frac_jewish)) +
  geom_point(size = 3) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dotted",
    color = "black",
    alpha = 0.7
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = paste0("Mean predicted probability (", n_bins, " bins)"),
    y = "Fraction Jewish",
    title = "Calibration Plot"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
print(w)


#---------------------------------------------------------
# Summary stats on total prob bounded and sampling prob
# For all respondents and Jewish respondents
#---------------------------------------------------------
summary(df_obs$total_probs_bounded)
summary(df_obs$sampling_prob, na.rm = T)

summary(df_obs$total_probs_bounded[df_obs$jewish == 1])
summary(df_obs$sampling_prob[df_obs$jewish == 1], na.rm = T)

counts = table(df_obs$total_probs_bounded)
sum(counts==1)
