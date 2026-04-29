library(tidyverse)
library(glue)
library(bayesplot)

source("simulation_tools.R")

results_savepath = "../../results/simulation/"
figure_savepath = "../../figures/simulation/"

true_val_path = fs::path(glue("{results_savepath}/true_probabilities"))
per_setting_results_path = fs::path(glue("{results_savepath}/per_setting_results"))

tempfile_path = fs::path(glue("{results_savepath}/tempfiles"))
df = merge_tf_csvs(tempfile_path, pattern = "^tf_.*\\.csv$")
df$n_label <- format(df$n, scientific = FALSE, big.mark = ",") #for plotting

# Get stats
max(df$a_max_accept_prop)
min(df$a_min_accept_prop)
mean(df$a_mean_prop_accept)
mean(df$a_min_accept_prop == 0)

mean(df$eta_rate)
max(df$eta_rate)
min(df$eta_rate)

# Get stats by num_sur - see goes up with higher num_sur
mean(df$a_mean_prop_accept[df$num_sur == 100])
mean(df$a_mean_prop_accept[df$num_sur == 1000])
mean(df$a_mean_prop_accept[df$num_sur == 10000])

mean(df$eta_rate[df$num_sur == 100])
mean(df$eta_rate[df$num_sur == 1000])
mean(df$eta_rate[df$num_sur == 10000])

# Get stats for the specific example given in paper
df %>% filter(num_sur == 10000 &
              n == 50000 &
              alpha_init_option == "uniform"
              ) %>%
  summarize(amean = mean(a_mean_prop_accept),
            cmean = mean(eta_rate))


#--------------------------------------
#--------------------------------------
# Posterior mean rho value plot
#--------------------------------------
#--------------------------------------
df_rho <- df[!is.na(df$postmean_rho_mean), ]

# replace with averages
df_rho <- df_rho %>% group_by(num_sur, n, n_label) %>%
  summarize(postmean_rho_mean = mean(postmean_rho_mean),
            postmean_rho_sd = mean(postmean_rho_sd),
            postmean_rho_max = mean(postmean_rho_max),
            postmean_rho_min = mean(postmean_rho_min)
            )

library(dplyr)
library(ggplot2)

df_rho <- df_rho %>%
  group_by(num_sur, n, n_label) %>%
  summarize(
    postmean_rho_mean = mean(postmean_rho_mean),
    postmean_rho_sd   = mean(postmean_rho_sd),
    postmean_rho_max  = mean(postmean_rho_max),
    postmean_rho_min  = mean(postmean_rho_min),
    .groups = "drop"
  )

p_rho = ggplot(
  df_rho,
  aes(
    x     = log10(num_sur),
    group = n_label
  )
) +
  # shaded min–max band
  geom_ribbon(
    aes(
      ymin = postmean_rho_min,
      ymax = postmean_rho_max,
      fill = n_label
    ),
    alpha = 0.05,
    color = NA
  ) +
  # faded min/max boundary lines (same color)
  geom_line(
    aes(y = postmean_rho_min, color = n_label),
    linewidth = 0.9,
    alpha = 0.5,
    linetype = "dashed"
  ) +
  geom_line(
    aes(y = postmean_rho_max, color = n_label),
    linewidth = 0.9,
    alpha = 0.5,
    linetype = "dashed"
  ) +
  # main mean line + points
  geom_line(
    aes(y = postmean_rho_mean, color = n_label),
    linewidth = 1.2
  ) +
  geom_point(
    aes(y = postmean_rho_mean, color = n_label),
    size = 2.8,
    alpha = 0.85
  ) +

  labs(
    title = expression("Relationship of " * rho[g] * " to Sample Size and Number of Surnames"),
    x     = expression("Number of surnames (" * log[10] * " scale)"),
    y     = expression("Posterior mean of mean " * rho[g] * ""),
    color = "Sample Size"
  ) +

  # keep one legend (color); drop fill legend
  guides(fill = "none") +

  scale_x_continuous(
    breaks = seq(floor(min(log10(df_rho$num_sur))), ceiling(max(log10(df_rho$num_sur))), by = 1),
    labels = seq(floor(min(log10(df_rho$num_sur))), ceiling(max(log10(df_rho$num_sur))), by = 1)
  ) +
  scale_y_continuous(limits = c(0, 1)) +

  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 12)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 10)),
    axis.text.x  = element_text(size = 15),
    axis.text.y  = element_text(size = 15),
    axis.ticks   = element_line(linewidth = 0.9),
    axis.ticks.length = unit(0.28, "cm"),
    legend.position = "right",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 15),
    legend.key.size = unit(1.3, "lines"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
    panel.grid.minor = element_blank()
  )
ggplot2::ggsave(filename = fs::path(glue("{figure_savepath}/rho_plot.pdf")),
                plot = p_rho,
                width = 16,
                height = 8,
                device = cairo_pdf # so that alpha shows up
)


# Note: MLE is unaffected by alpha init option so pool together
#
# Creates two plots
#
# Plot 1: metric by number of iterations, faceted by num surnames, and
# either filtered to a specific n or also faceted by n. Colored by alpha init option
# (or MLE)
#
#
# Plot 2: metric by n, faceted by num surnames, colored by alpha init option and
# MLE


plotter <- function(df,
                    metric = c("TV", "rmse"),
                    parameter = c("theta", "alpha"),
                    plot1_n_value = NULL
                    ) {

  required <- c("alpha_init_option", "num_sur", "n", "niter", "n_label")
  stopifnot(all(required %in% names(df)))

  # Get relevant columns
  cols = colnames(df)
  cols = cols[grepl(parameter, cols)]
  cols = cols[grepl(metric, cols)]

  postmean_cols = cols[grepl("postmean", cols)]
  mle_col = cols[grepl("mle", cols)]

  df_sub = df[,c(required, postmean_cols, mle_col)]

  #-----------------------------------------------------------------------------
  df_long <- df_sub %>%
    pivot_longer(
      cols = all_of(c(postmean_cols, mle_col)),
      names_to = "estimator_type",
      values_to = "estimator_value"
    )

  df_long$estimator_iter = suppressWarnings(
    as.numeric(sub(".*postmean([0-9]+)_.*", "\\1", df_long$estimator_type))
  )

  # Define groups alpha_init_options + mle
  df_long$group = df_long$alpha_init_option
  df_long$group[df_long$alpha_init_option == "uniform"] = "Posterior mean with uniform initial \u03B1"
  df_long$group[df_long$alpha_init_option == "best_guess"] = "Posterior mean with smoothed MLE initial \u03B1"
  df_long$group[df_long$estimator_type == mle_col] = "MLE" #this must be last

  # Get averages, sd
  df_plot <- df_long %>%
    group_by(group, num_sur, n, n_label, estimator_type, estimator_iter) %>%
    summarize(
      estimator_mean = mean(estimator_value),
      sd   = sd(estimator_value),
      nrep = n(),
      se   = sd / sqrt(nrep),
      ymin = estimator_mean - 2 * se,
      ymax = estimator_mean + 2 * se,
      .groups = "drop"
    )

  if(parameter == "alpha"){
    param_label = "\u03B1"
  }
  if(parameter == "theta"){
    param_label <- "\u03B8"
  }

  # Avoid chaing df_plot since also need it for 2nd plot
  df_iter    <- df_plot %>% filter(!is.na(estimator_iter))
  df_no_iter <- df_plot %>% filter(is.na(estimator_iter))

  if(!is.null(plot1_n_value )){
    df_iter <- df_iter %>% filter(n == plot1_n_value)
    df_no_iter <- df_no_iter %>% filter(n == plot1_n_value)
  }



  p <- ggplot() +
    # cases with an iteration (points/lines across x)
    geom_line(
      data = df_iter,
      aes(x = estimator_iter, y = estimator_mean, color = group,
          group = interaction(group, num_sur)
          ),
      linewidth = 1.2,
      alpha = 0.7
    ) +
    geom_point(
      data = df_iter,
      aes(x = estimator_iter, y = estimator_mean, color = group),
      size = 1.5,
      alpha = 0.9
    ) +
    # geom_errorbar(data = df_iter,
    #   aes(x = estimator_iter,
    #       y = estimator_mean,
    #       ymin = ymin,
    #       ymax = ymax,
    #       color = group),
    #   width = 0.08) +
    # cases with no iteration: horizontal line at estimator_value
    geom_hline(
      data = df_no_iter,
      aes(yintercept = estimator_mean, color = group),
      linewidth = 1.2,
      alpha = 0.7
    ) +
    facet_grid(
      rows = vars(n_label),
      cols = vars(num_sur),
      labeller = labeller(
        num_sur = function(x) paste("Number of surnames =", x),
        n_label = function(x) paste("Sample size =", x)
      )
    ) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = "Iteration",
         y = glue("Mean {metric}({param_label})"),
         color = "Estimator:",
         title = "") +
    theme_bw() +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text  = element_text(size = 20),

      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12),

      axis.ticks = element_line(linewidth = 0.8),
      axis.ticks.length = unit(0.25, "cm"),

      legend.key.size = unit(1.2, "lines"),
      strip.text = element_text(size = 14, face = "bold") #facet label size
    ) + scale_color_manual(
      values =c(
        "#0072B2",  # blue
        "#D55E00",  # orange/red
        "#009E73"   # bluish green
      )
    )

  #--------------------------------------------------------------
  # Plot 2 - subset only to max iter
  #--------------------------------------------------------------
  m = max(df_iter$estimator_iter)
  df_iter_max <- df_plot %>% filter(is.na(estimator_iter) | estimator_iter == m)

  q <- ggplot() +
    # cases with an iteration (points/lines across x)
    geom_line(
      data = df_iter_max,
      aes(x = n, y = estimator_mean, color = group,
          group = interaction(group, num_sur)
      ),
      linewidth = 1.2,
      alpha = 0.7
    ) +
    geom_point(
      data = df_iter_max,
      aes(x = n, y = estimator_mean, color = group),
      size = 1.5,
      alpha = 0.9
    ) +
    facet_grid(
      cols = vars(num_sur),
      labeller = labeller(
        num_sur = function(x) paste("Number of surnames =", x)
      )
    )  +
    scale_y_continuous(limits = c(0, NA)) +
    scale_x_continuous(limits = c(0, NA)) +
    labs(x = "Sample size",
         y = glue("Mean {metric}({param_label})"),
         color = "Estimator:",
         title = "") +
    theme_bw() +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text  = element_text(size = 20),

      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12),

      axis.ticks = element_line(linewidth = 0.8),
      axis.ticks.length = unit(0.25, "cm"),

      legend.key.size = unit(1.2, "lines"),
      strip.text = element_text(size = 14, face = "bold") #facet label size
    ) + scale_color_manual(
      values =c(
        "#0072B2",  # blue
        "#D55E00",  # orange/red
        "#009E73"   # bluish green
      )
    )
  #TODO: make it so x axis shows up with only values actually
  # have data for

  return(list("plot1" = p, "plot2" = q))

}

# NOTE: if do metric = rmse, see it get smaller as num_sur increases purely because
# the parameter values are of smaller magnitude because have to sum to 1 and there are
# more of them. TV does not have this problem and more intuitively reflects idea that
# in relative terms, esitmation is worse when have more parameters
out <- plotter(df, parameter = "alpha", metric = "TV",
               plot1_n_value = 50000)
p = out$plot1
q = out$plot2
print(p)
print(q)
ggplot2::ggsave(filename = fs::path(glue("{figure_savepath}/alpha_grid_results.pdf")),
                plot = p,
                width = 16,
                height = 8,
                device = cairo_pdf # so that alpha shows up
                )
ggplot2::ggsave(filename = fs::path(glue("{figure_savepath}/alpha_grid_results2.pdf")),
                plot = q,
                width = 16,
                height = 8,
                device = cairo_pdf # so that alpha shows up
)

out <- plotter(df, parameter = "theta", metric = "TV",
               plot1_n_value = 50000)
p = out$plot1
q = out$plot2
print(p)
print(q)

ggplot2::ggsave(filename = fs::path(glue("{figure_savepath}/theta_grid_results.pdf")),
                plot = p,
                width = 16,
                height = 8,
                device = cairo_pdf # so that alpha shows up
)
ggplot2::ggsave(filename = fs::path(glue("{figure_savepath}/theta_grid_results2.pdf")),
                plot = q,
                width = 16,
                height = 8,
                device = cairo_pdf # so that alpha shows up
)
