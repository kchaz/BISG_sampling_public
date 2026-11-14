#https://rfortherestofus.com/2024/06/us-maps
library(tidyverse)
library(usmap)

p_R1_gvn_state <- readRDS("../../data/DataObjects/p_R1_gvn_state.RDS")


# Probabilities vector with state abbreviations
probs <- p_R1_gvn_state

# Convert to data frame and join with state shapefile
df_probs <- data.frame(
  state = names(p_R1_gvn_state),
  "Proportion Jewish" = as.numeric(p_R1_gvn_state),
  check.names = F #don't convert space to .
)

map = plot_usmap(
  color = "white",
  linewidth = 0.2,
  regions = "states",
  data = df_probs,
  values = "Proportion Jewish"
) +
  scale_fill_gradient(low = "white", high = "black", na.value = "blue") +
  theme(legend.position = 'top',
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 20, face = "bold")
        )+
  guides(
    fill = guide_colorbar(
      barwidth = unit(10, 'cm')
    )
  )

ggsave(
  filename = "../../paper/Figures/prevalence_map.pdf",
  plot = map,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600
)
