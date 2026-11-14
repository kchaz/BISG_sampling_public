library(tidyverse)
library(readr)
library(scales)

df = read_csv("../../data/ObitData/obit_list_092424-fulldehyph-check.csv")

#----------------------------------
# obituary observations per state

length(unique(df$state))
counts = table(df$state)

counts[counts<=50]
counts[counts>50 & counts <= 100]
counts[counts>101 & counts <= 500]
counts[counts>501 & counts <= 1000]
counts[counts>1001 & counts <= 5000]
counts[counts>5001 & counts <= 10000]
counts[counts>10001 & counts <= 20000]
counts[counts>20001 & counts <= 30000]
counts[counts>30000]

length(counts[counts>1000])

#-------------------------------------------------------------------------------
# Explore Age
#-------------------------------------------------------------------------------

# Extract 4-digit years starting with 19 or 20 to get years of birth and death
df$yob <- as.numeric(stringr::str_extract(df$dob, "\\b(19|20)\\d{2}\\b"))
df$yod <- as.numeric(stringr::str_extract(df$dod, "\\b(19|20)\\d{2}\\b"))

# age distribution
summary(df$age, na.rm = T)
sd(df$age, na.rm = T)

median(df$yob, na.rm = T)
max(df$yod, na.rm = T)
min(df$yod, na.rm = T)


# Turn this into a 1 x 3 plot, publication ready in terms of aesthetics
# On each plot, mark in the title or wherever else makes sense the % NA
# for that variable
age_na = mean(is.na(df$age))
yob_na = mean(is.na(df$yob))
yod_na = mean(is.na(df$yod))
par(mfrow = c(1,3))
hist(df$age)
hist(df$yob)
hist(df$yod)


# ChatGPT generated version of plot
#----------------------------------------------
# Variables to plot
vars <- c("age", "yob", "yod")

# NA rates
age_na <- mean(is.na(df$age))
yob_na <- mean(is.na(df$yob))
yod_na <- mean(is.na(df$yod))

# Labels for facet strips with % NA
facet_labels <- c(
  age = paste0("Age at Death (NA: ", percent(age_na, accuracy = 0.1), ")"),
  yob = paste0("Year of Birth (NA: ", percent(yob_na, accuracy = 0.1), ")"),
  yod = paste0("Year of Death (NA: ", percent(yod_na, accuracy = 0.1), ")")
)

# Long data for faceting; drop NAs per variable for clean histograms
df_long <- df %>%
  select(all_of(vars)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  filter(!is.na(value))

# Plot
p <- ggplot(df_long, aes(x = value)) +
  geom_histogram(bins = 30, color = "white") +
  facet_wrap(~ variable, nrow = 1, scales = "free_x",
             labeller = labeller(variable = facet_labels)) +
  labs(
    title = "Distributions of Age at Death, Year of Birth, and Year of Death of People in Jewish Obituaries 1999–2024",
    x = NULL,
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 16),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.3),
    panel.grid.major.y = element_line(linewidth = 0.3),
    plot.title = element_text(face = "bold", size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.ticks.x = element_line(linewidth = 0.8),
    axis.ticks.y = element_line(linewidth = 0.8),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),  # <-- adds light gray borders
    panel.spacing.x = unit(0.8, "lines"),
    panel.background = element_rect(fill = "white", color = "gray85")
)
print(p)
ggsave("../../paper/Figures/age_yob_yod.pdf", p, width = 15, height = 6, dpi = 320)


