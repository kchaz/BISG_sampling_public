
library(tidyverse)
library(cowplot)
library(glue)

# For decoding bin counts below
decode_vec <- function(s) {
  parts <- strsplit(s, "\\|")[[1]]
  kv <- strsplit(parts, "=")
  vals <- as.numeric(sapply(kv, `[`, 2))
  names(vals) <- sapply(kv, `[`, 1)
  vals
}

df_tot = read.csv("../../data/cluster_output/probability_summaries_by_state.csv")
df_tot$num_zero = df_tot$state_counts - df_tot$state_nonzeroprob_counts


# How often bound method applied
df_tot$num_with_bounds_applied = NULL #decided to do this another way - rounding error issues here


# How often >1 trunc applied
range(round(df_tot$num_with_trunc_at_1/df_tot$state_nonzeroprob_counts,3))
sum(df_tot$num_with_trunc_at_1)/sum(df_tot$state_nonzeroprob_counts)


# HISTOGRAMS
#--------------------
create_plots <- function(counts, num_zero, title_label){

    # extract bin edges from names
    bins <- t(sapply(
      names(counts),
      function(s) as.numeric(
        unlist(regmatches(s, gregexpr("[0-9.]+", s)))
      )
    ))
    left  <- bins[,1]
    right <- bins[,2]

    df_plot <- data.frame(
      bin = names(counts),
      count = as.numeric(counts),
      left = left,
      right = right
    ) %>% filter(count > 0) #remove 0 count bins

    # Remove the zero counts from the first bin
    df_plot[1,"count"] = df_plot[1,"count"] - num_zero

    # Calculate proportion
    df_plot$prop = df_plot$count/sum(df_plot$count)

    p = ggplot(df_plot) +
      geom_rect(
        aes(xmin = left, xmax = right, ymin = 0, ymax = prop),
        fill = "black",
        color = "white",
        linewidth = .8
      ) +
      scale_y_continuous(
        limits = c(0, 1),
        expand = c(0, 0)
      ) +
      labs(
        x = "Value",
        y = "Proportion",
        title = glue("{title_label} > 0")
      ) +
      theme_minimal(base_size = 18) +
      theme(
        plot.title = element_text(
          size = 22,
          face = "bold",
          hjust = 0.5
        ),
        axis.title = element_text(size = 18),
        axis.text  = element_text(size = 18),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),

        # Grid tuning
        panel.grid.major = element_line(color = "gray85",
                                        linewidth = 0.4),
        panel.grid.minor = element_blank(),

        # Remove panel border explicitly
        panel.border = element_blank()
      )


    q = ggplot(df_plot[3:nrow(df_plot),]) +
      geom_rect(
        aes(xmin = left, xmax = right, ymin = 0, ymax = prop),
        fill = "black",
        color = "white",
        linewidth = .3
      ) +
      labs(
        x = "Probability",
        y = "Proportion",
        title = glue("{title_label} > 0.02")
      ) +
      theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(
          size = 20,
          face = "bold",
          hjust = 0.5
        ),
        axis.title = element_text(size = 18),
        axis.text  = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),

        # Grid tuning
        panel.grid.major = element_line(color = "gray85", linewidth = 0.4),
        panel.grid.minor = element_blank(),

        # Remove panel border explicitly
        panel.border = element_blank()
      )

    return(list("p"=p,
                "q"=q)
    )

}



num_zero = sum(df_tot$num_zero)

counts = Reduce(`+`, lapply(df_tot$total_probs_bounded_bins, decode_vec))

out = create_plots(counts, num_zero, title_label = "Unnormalized probabilities")

fp1 = ggdraw(out$p) +
  draw_plot(
    out$q + theme(
      plot.background = element_rect(
        fill = "white",
        color = "black",
        linewidth = 0.5
      )
    ),
    x = 0.35, y = 0.35,
    width = 0.55, height = 0.55
  )
print(fp1)
ggsave(
  filename = "../../paper/Figures/unnormalized_prob_histogram.pdf",
  plot = fp1,
  width = 14,
  height = 9,
  units = "in",
  dpi = 300
)



counts = Reduce(`+`, lapply(df_tot$sampling_prob_bins, decode_vec))
out = create_plots(counts, num_zero, title_label = "Sampling probabilities")

fp2 = ggdraw(out$p) +
  draw_plot(
    out$q + theme(
      plot.background = element_rect(
        fill = "white",
        color = "black",
        linewidth = 0.5
      )
    ),
    x = 0.35, y = 0.35,
    width = 0.55, height = 0.55
  )
print(fp2)
ggsave(
  filename = "../../paper/Figures/normalized_prob_histogram.pdf",
  plot = fp2,
  width = 14,
  height = 9,
  units = "in",
  dpi = 300
)



library(patchwork)
fp = fp1 + fp2
ggsave(
  filename = "../../paper/Figures/sampling_prob_histogram.pdf",
  plot = fp,
  width =25,
  height = 10,
  units = "in",
  dpi = 300
)

