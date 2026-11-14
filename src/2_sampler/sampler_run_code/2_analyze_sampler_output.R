#'
#' Analyze output of MCMC sampler based on heirarchical Bayesian model
#'
#' Author: Kyla Chasalow
#'
# Libraries
library(readr)

# Scripts
source("paths.R")
source(fs::path(src_path,"2_sampler/HBayesSampler/data_prep_helpers.R"))
source(fs::path(src_path,"2_sampler/HBayesSampler/sampler_analysis_helpers.R"))


#################################################

save_path = paste0(save_sampler_output_path, "/chain_analysis/")
if (!dir.exists(save_path)) {
  # Create the directory (including any necessary parent directories)
  dir.create(save_path, recursive = TRUE)
}

paper_figure_save_path = "../../../paper/Figures/sampler_output/"

##########################################################################
#' LOAD SAMPLER OUTPUT
##########################################################################

# get filename (should only be one that matches this description)
file_list <- list.files(save_sampler_output_path,
                        pattern = "_realdf_45000\\.RDS$",
                        full.names = TRUE)


filename = file_list[1]

#filename = "20241005205601_realdf_100.RDS" #TEMP
out = readRDS(filename)

# for backwards compatibility
if ("c_chain" %in% names(out)){
  out$eta_chain = out$c_chain
  out$num_accept_eta = out$num_accept_c
}

###############################
# Load Main Data and Key Components
###############################
df = read_csv(paste0(obit_data_path,"obit_list_092424-fulldehyph-check.csv"))
n = nrow(df)

# get list of *all* 51 state names (+DC)
states <- read_csv(paste0(data_path,"states.csv"))
loc_names <- sort(states$Abbreviation)


# get list of surnames
surnames = unique(df$last)
num_sur = length(surnames)

# get count matrix
D <- get_count_mat(df, loc_names,  surnames,surname_colname = "last")
ns <- apply(D, 1, sum)
ng <- apply(D, 2, sum)

#rownames(out$a_chain) = rownames(D) #TEMP
#colnames(out$a_chain) = colnames(D) #TEMP


##########################################################################
#' LOCATION-JEWISH DATA (used for exploring outputs)
##########################################################################
# see generate_fips_state_proportions.R for code to generate these
pR1 = readRDS(paste0(data_object_path,"pR1.RDS"))
p_R1_gvn_fips = readRDS(paste0(data_object_path,"p_R1_gvn_fips.RDS"))
p_R0_gvn_fips = readRDS(paste0(data_object_path,"p_R0_gvn_fips.RDS"))
p_fips_gvn_R1 = readRDS(paste0(data_object_path,"p_fips_gvn_R1.RDS"))
p_fips_gvn_R0 = readRDS(paste0(data_object_path,"p_fips_gvn_R0.RDS"))

p_state_gvn_R1 = readRDS(paste0(data_object_path,"p_state_gvn_R1.RDS"))
p_state_gvn_R0 = readRDS(paste0(data_object_path,"p_state_gvn_R0.RDS"))
p_R1_gvn_state = readRDS(paste0(data_object_path,"p_R1_gvn_state.RDS"))
p_R0_gvn_state = readRDS(paste0(data_object_path,"p_R0_gvn_state.RDS"))


##########################################################################
#' Analysis
##########################################################################

#Selection of names of different frequencies to look at
NAMES = c("COHEN", #3274 COUNT
          "BERKOWITZ", #286 COUNT
          "SNYDER", #263 COUNT
          "BERGMAN", #100 count
          "HIRSH", #68 COUNT
          "ZUCKERBERG" #10 COUNT
          )  # 1 COUNT
num.obs = table(df$last)[NAMES]
labels = sapply(NAMES, function(x){paste0(x," (", num.obs[x], ")")})
ids = sapply(NAMES, function(x){which(rownames(out$a_chain) == x)}) #to get right order

################################################################
# Trace plots and accept rates
################################################################

# Get Stats
stats = get_sampler1_accept_stats(out)
saveRDS(stats, paste0(save_path,"main_sampler_stats.RDS"))

# Plot Trace Plots for subset of names
pdf(
  file = paste0(paper_figure_save_path, "tr_plot_eta.pdf"),
  width = 9,    # inches
  height = 5
)
plot_eta_chain(out, burnin = 10)
dev.off()


pdf(
  file = paste0(paper_figure_save_path,"tr_plot_a.pdf"),
  width = 9,    # inches
  height = 6
  )
plot_a_chain(out, indices = ids, labels = labels) #can plot 6 at a time
dev.off()

################################################################
# Rho
################################################################
# # with less burnin for trace plot
# rho_out = analyze_rho(out, ng, burnin = 100, upper = NULL) #note 100 corresponds to first 500
# pdf(
#   file = paste0(paper_figure_save_path, "tr_plot_rho_NJ.pdf"),
#   width = 9,    # inches
#   height = 5
# )
# plot_rho_chain(rho_out$rho_mat, "NJ")
# dev.off()
#Note: the above plot is basically same pattern as eta plot...unsurprisingly,
#because it's same transformation of eta at each step...so not using it

rho_out = analyze_rho(out, ng, burnin = 3001, upper = NULL)

# build/ensure plot_df has the right types
plot_df <- data.frame(
  state = as.character(names(rho_out$rho_postmean)),
  rho_g = as.numeric(rho_out$rho_postmean),
  stringsAsFactors = FALSE
)
plot_df = plot_df[order(plot_df$rho_g, decreasing = T),]

# keep rows with rho_g < 1
plot_df <- plot_df[plot_df$rho_g < 1, , drop = FALSE]

# append the "All Others" row as a data.frame (keeps types)
all_others <- data.frame(state = "All Others", rho_g = 1, stringsAsFactors = FALSE)
plot_df <- rbind(all_others, plot_df)

# set to factor with factor levels in the current row order (useful for plotting)
plot_df$state <- factor(plot_df$state, levels = plot_df$state)

library(ggplot2)
p <- ggplot(plot_df, aes(x = state, y = rho_g)) +
  geom_point(
    size = 2.2,
    color = "gray20"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "gray50"
  ) +
  coord_flip() +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Posterior Mean of Shrinkage Parameter by State",
    y = expression(rho[g]),
    x = "",
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold",
                              hjust = 0.5,
                              size = 26),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_text(size = 22)
  )

ggsave(
  filename = paste0(paper_figure_save_path,"posterior_mean_rho_g.pdf"),
  plot = p,
  width = 12,
  height = 10,
  units = "in",
  dpi = 1200
)


################################################################
# Calculate final posterior means
Theta_means = out$theta_sums/out$sums_denom

################################################################
#Plotting some example names - observed vs state proportions
################################################################
#normalize D by column to get per-state name proportions
Dnorm = apply(D,2, function(x){
  if(sum(x)>0){
    return(x/sum(x))
  }
  return(rep(0,length(x)))
}
)
rownames(Dnorm) = rownames(D)
apply(Dnorm, 2, sum)


namescatter <- function(name, Theta_means, plot_state_comparison = T){
  par(mfrow = c(1,1))

  u = max(Dnorm[name,], Theta_means[name,])*1.05

  plot(Dnorm[name,], Theta_means[name,], "n",
       main = paste0("Estimate vs State Proportion for ",name),
       ylab = expression(theta),
       xlab = "Observed State Proportion",
       xlim = c(0,u),
       ylim = c(0,u)
  )
  text(Dnorm[name,], Theta_means[name,],
       colnames(D), col = "black", cex = .80)
  ps = ns/n
  abline(v = ps[name], lty = 2, col = "blue")
  abline(h = ps[name], lty = 2, col = "blue")
  abline(a=0,b=1,lty=2,col = "gray")
  legend("topleft", legend = c("Overall Observed Proportion"),
         col = c("blue"), lty = 2,
         bty = "n")

  if (plot_state_comparison){
    par(mfrow = c(1,2))
    #plot(Theta_means1[name,], )
    plot(p_R1_gvn_state, Theta_means1[name,],
         ylab = expression(theta),
         xlab = "P(R=1|State)")
    abline(lm(Theta_means[name,] ~ p_R1_gvn_state), lty = 2, col = "gray")

    #plot(Theta_means1[name,], )
    plot(p_state_gvn_R1, Theta_means1[name,],
         ylab = expression(theta),
         xlab = "P(State|R=1)")
    abline(lm(Theta_means[name,] ~ p_state_gvn_R1), lty = 2, col = "gray")
  }
}


for (name in NAMES){
  pdf(paste0(paper_figure_save_path, name, "_compare_plot.pdf"))
  namescatter(name, Theta_means, plot_state_comparison = F)
  dev.off()
}



###############################################################################
# Some summary stats
###############################################################################

hist(Theta_means)
mean(Theta_means)
range(Theta_means)
