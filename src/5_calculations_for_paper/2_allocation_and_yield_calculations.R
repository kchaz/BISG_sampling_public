library(tidyverse)
library(ggplot2)
library(patchwork) #for creating plot grids
library(glue)

old_method_color = "black"
new_method_color = "blue"
results_color = "red"


################################################################################
#' Uses output of summarize_probabilities_by_state.R
#' (1) Compares allocation methods and plots scatterplot of them against each other
#' (2) Compares expected yields for basic methods and our method
#'
################################################################################

# Load File
#df = read.csv("../../data/cluster_output/normalizing_totals_by_state-original.csv") #one with bug targets
df = read.csv("../../data/cluster_output/probability_summaries_by_state.csv") #updated one

# Set overall target
Target = 50000

# State distributions from AJPP
p_R1_gvn_state = readRDS(fs::path("../../data/DataObjects/p_R1_gvn_state.RDS"))
p_state_gvn_R1 = readRDS(fs::path("../../data/DataObjects/p_state_gvn_R1.RDS"))

df$p_R1_gvn_state = p_R1_gvn_state[df$state]
df$p_state_gvn_R1 = p_state_gvn_R1[df$state]


get_Tg <- function(numerators, Target, tol = 5){
    sample_fracs = numerators/sum(numerators)
    Tg = round(Target*sample_fracs)
    stopifnot(abs(sum(Tg) - Target) <= tol)
    return(Tg)
}


##########################################
# what actually did (bug)
##########################################
numerators = sqrt(df$p_R1_gvn_state * df$prob_totals)
Tg = get_Tg(numerators, Target)
df$target_bug = Tg

# stat related to impact of bug
#cor(sqrt(df$prob_totals), df$state_nonzeroprob_counts)


##########################################
# strat with filter (ad-hoc - only replace the one)
##########################################
# numerators = sqrt(df$p_R1_gvn_state) * df$state_nonzeroprob_counts #frame count is filtered count
# sample_fracs = numerators/sum(numerators)
# Tg = round(Target*sample_fracs)
# stopifnot(abs(sum(Tg) - Target) <= 5)
#
# # Add targets to data frame
# df$target_strat_filter_ad_hoc = Tg

##########################################
# strat with no filter adjustment
##########################################
numerators = sqrt(df$p_R1_gvn_state) * df$state_counts
Tg = get_Tg(numerators, Target)
df$target_strat_non_filter = Tg


##########################################
# strat with filtering only using best case
# P(R=1|g,H=1)
##########################################
numerators = sqrt(df$p_R1_gvn_state * df$state_counts/df$state_nonzeroprob_counts) * df$state_nonzeroprob_counts
Tg = get_Tg(numerators, Target)
df$target_strat_filter = Tg

##########################################
# strat with filter-updated P(R=1|g,H=1)
# Problem with this is if doing only filtering
# You wouldn't have the P(R=1|s,g) to calculate this
##########################################
# numerators = sqrt(df$prob_squared_totals * df$state_nonzeroprob_counts)#frame count is filtered count
# sample_fracs = numerators/sum(numerators)
# Tg = round(Target*sample_fracs)
# stopifnot(abs(sum(Tg) - Target) <= 5)
# df$target_strat_filter2 = Tg

##########################################
# just do proportionate to state proportions P(R=1|g)
##########################################
numerators = df$p_R1_gvn_state
Tg = get_Tg(numerators, Target)
df$target_pR1_gvn_g = Tg

##########################################
# do proportionate to P(g | R=1)
##########################################
sample_fracs = df$p_state_gvn_R1
Tg = get_Tg(numerators, Target)
df$target_pg_gvn_R1 = Tg

###########################################
# T_opt calculation for Pois - ad hoc still using p(g|R=1) original
###########################################
# numerators = df$p_state_gvn_R1  * sqrt(df$prob_totals) / sqrt(df$prob_squared_totals)
# Tg = get_Tg(numerators, Target)
# df$target_pois = Tg

###########################################
# T_opt calculation for Pois - with NO filtering
# and a tiny positive probability assigned to all units with surname not observed in obit data
###########################################
df$epsilon = 10^-7 #based on  df$prob_min

df$num_zero = df$state_counts - df$state_nonzeroprob_counts

# Update the prob totals as if P(R=1|s,g) = \epsilon for these units
df$prob_totals_without_filter = df$prob_totals + df$num_zero * df$epsilon
df$prob_squared_totals_without_filter = df$prob_squared_totals + (df$num_zero * df$epsilon^2)
df$prob_4thpower_totals_without_filter = df$prob_4thpower_totals + (df$num_zero * df$epsilon^4)

numerators = df$p_state_gvn_R1 / sqrt(df$prob_squared_totals_without_filter / df$prob_totals_without_filter)
Tg = get_Tg(numerators, Target)
df$target_pois_without_filter = Tg

########################################################
# Sensitivity Analysis
# T_opt calculation for Pois - with NO filtering
# Suppose P(R=1|s,g) is REALLY small but when calculating our
# sampling probs, we use an estimate that is large
########################################################
eps_true = 10^-7
eps_used = 10^-2

# Update the prob totals
df$prob_totals_without_filter2 = df$prob_totals + df$num_zero * eps_used#This is only based on sampling probs
df$prob_squared_totals_without_filter2 = df$prob_squared_totals + (df$num_zero * eps_true * eps_used) #square comes from est * true being same in ideal calc
df$prob_4thpower_totals_without_filter2 = df$prob_4thpower_totals + (df$num_zero * (eps_true * eps_used)^2)

numerators = df$p_state_gvn_R1 / sqrt(df$prob_squared_totals_without_filter2 / df$prob_totals_without_filter2)
Tg = get_Tg(numerators, Target)
df$target_pois_without_filter2 = Tg


###########################################
# T_opt calculation for Pois - with full filter adjustment
###########################################
numerators = df$prob_totals  * sqrt(df$prob_totals) / sqrt(df$prob_squared_totals)
Tg = get_Tg(numerators, Target)
df$target_pois_filter = Tg


##########################################
##########################################
# Compare
##########################################
##########################################

label_list = list("target_bug" = "Actually used",
                 "target_pois_filter" = "Poisson sampling \n optimal with \n filter adjustment",
                 #"target_pois_without_filter" = "Poisson sampling\n no filtering",
                 "target_strat_non_filter" = "Disproportionate \n stratification with \n original counts",
                 "target_strat_filter" = "Disproproportionate \n stratification with \n filter adjustment",
                 "target_pg_gvn_R1" = "Proportional \nto Pr(g|R=1)",
                 "target_pR1_gvn_g" = "Proportional \nto Pr(R=1|g)"
)



# #################################

# columns to compare
cols <- names(label_list)
X <- df[, cols, drop = FALSE]
pretty_labs <- unname(unlist(label_list))

upper <- max(X, na.rm = TRUE)*1.1
lower <- 0

# --- lower: scatter + x=y + labels ---
panel_scatter <- function(x, y, ...) {
        points(x, y, pch = 19, ...)
        abline(a = 0, b = 1, lty = 2)
        text(x, y, labels = df$state, cex = 0.9, pos = 4)
}


# --- upper: shaded cell + correlation text ---
panel_cor_heat <- function(x, y, digits = 4, ...) {
        r <- cor(x, y, use = "pairwise.complete.obs")

        # map r in [-1,1] to a diverging palette
        pal <- colorRampPalette(c("navy", "white", "firebrick"))(201)
        idx <- 1 + round((r - (-1)) / 2 * (length(pal) - 1))
        idx <- max(1, min(length(pal), idx))

        # draw a full cell background
        usr <- par("usr"); on.exit(par(usr = usr))
        par(usr = c(0, 1, 0, 1))
        rect(0, 0, 1, 1, col = pal[idx], border = "grey85")

        # print correlation
        text(0.5, 0.5, sprintf(paste0("%.", digits, "f"), r), cex = 2.5)
}


pdf("../../paper/Figures/allocation_comparison_plot.pdf", width = 10, height = 10)
pairs(
        X,
        labels = pretty_labs,
        lower.panel = panel_scatter,
        upper.panel = panel_cor_heat,
        diag.panel  = NULL,              # optional: remove diagonal
        xlim = c(lower, upper),
        ylim = c(lower, upper)
)
dev.off()


################################################################################
################################################################################
################################################################################
# Expected sample yield table
################################################################################
################################################################################
################################################################################
get_expected_yield_from_stratified_sampling <- function(Target, Tgvec, pvec){
        return(sum(Tgvec * pvec)/Target)
}

get_SRS_yield_frac_variance <- function(N, n, p){
        1/n * p * (1-p) * (N - n)/(N-1)
}


################################################################################
# Simple random sample
################################################################################
N = sum(df$state_counts)
n = Target

# Note to self: proportionate allocation based on just state size will be equivalent
# in yield but different in variance
yield_frac1 = sum(df$state_counts * df$p_R1_gvn_state)/N
SD1 = sqrt(get_SRS_yield_frac_variance(N=N,
                                  n=n,
                                  p=yield_frac1))
print(c(yield_frac1, SD1)*100)


################################################################################
# Proportionate allocation based on just state size
# without regard to size of Jewish population
################################################################################
num = df$state_counts
sample_frac = num / sum(num)
Tg = round(Target * sample_frac)
stopifnot(round(sum(Tg)) == Target)
yield_frac2 = get_expected_yield_from_stratified_sampling(Target, Tg, df$p_R1_gvn_state)

per_state_var = get_SRS_yield_frac_variance(N = df$state_counts,
                                       n = Tg,
                                       p = df$p_R1_gvn_state
                       )
SD2 = sqrt(sum(per_state_var * Tg^2)/n^2)
print(c(yield_frac2, SD2)*100)
# Variance is only marginally smaller - not worth a separate row


################################################################################
# Disproportionate allocation
# stratified sample
################################################################################
yield_frac3 = get_expected_yield_from_stratified_sampling(Target,
                                            df$target_strat_non_filter,
                                            df$p_R1_gvn_state)

per_state_var = get_SRS_yield_frac_variance(N = df$state_counts,
                                       n = df$target_strat_non_filter,
                                       p = df$p_R1_gvn_state
)
SD3 = sqrt(sum(per_state_var * df$target_strat_non_filter^2)/n^2)
print(c(yield_frac3, SD3)*100)



################################################################################
# Disproportionate allocation
# Stratified sample
# Only now with filter adjustment
# i.e. sampling only from the filtered names
# and with updated probability P(R=1|g)
################################################################################

# Filter adjusted P(R=1|G=g,H=1)
df$p_R1_gvn_state_filtered_true = df$prob_totals / df$state_nonzeroprob_counts
df$p_R1_gvn_state_filtered_best_case = df$p_R1_gvn_state * df$state_counts / df$state_nonzeroprob_counts

# NOTE: want the allocation to be as if did not have P(R=1|s,g) estimates
# So for allocation, we will update Pr(Jewish | g) under assumption that
# filtering only remove non-Jewish. This is something can do without estimating P(R=1|s,g)
num = sqrt(df$p_R1_gvn_state_filtered_best_case) * df$state_nonzeroprob_counts
sample_frac = num / sum(num)
Tg = Target * sample_frac
stopifnot(round(sum(Tg)) == Target)

# BUT then we will calculate expected yield using the actual
yield_frac4 = get_expected_yield_from_stratified_sampling(Target, Tg, df$p_R1_gvn_state_filtered_true)


per_state_strat_yield_frac_var = get_SRS_yield_frac_variance(N = df$state_nonzeroprob_counts,
                                       n = Tg,
                                       p = df$p_R1_gvn_state_filtered_true
)
df$per_state_strat_yield_frac_sd = sqrt(per_state_strat_yield_frac_var)

SD4 = sqrt(sum(per_state_strat_yield_frac_var *  Tg^2)/n^2)
print(c(yield_frac4, SD4)*100)


################################################################################
# Poisson method - with filtering
################################################################################

get_pois_yield_analysis <- function(
                state_labels, Tg, prob_total, prob_squared_total, prob_4thpower_totals){

        # Expectation of frac of rare population sampled - by state
        per_state_pois_yield_count = (Tg / prob_total  * prob_squared_total)
        per_state_pois_yield_frac = per_state_pois_yield_count/Tg

        # Expectation of frac of rare population sampled - overall
        yield_count = sum(Tg/prob_total * prob_squared_total)
        yield_frac = yield_count/Target

        # Variance of sample size n_g  by state
        per_state_var_term2 = (Tg^2 / prob_total^2) * prob_squared_total
        per_state_pois_ng_vars = Tg - per_state_var_term2
        per_state_pois_ng_sd = sqrt(per_state_pois_ng_vars)
        range_ng_sd = range(per_state_pois_ng_sd)
        max_ng_sd = which.max(per_state_pois_ng_sd)
        min_ng_sd = which.min(per_state_pois_ng_sd)

        # Variance of overall sample size n
        SD_n = sqrt(Target - sum(per_state_var_term2)) #Overall

        # Variance of sample size n_{g1} - by state rare population
        p1 = (Tg/prob_total) * prob_squared_total
        p2 = (Tg^2/prob_total^2) * prob_4thpower_totals
        per_state_pois_nrare_vars = p1-p2
        per_state_pois_nrare_sd = sqrt(per_state_pois_nrare_vars)

        range_ngrare_sd = range(per_state_pois_nrare_sd)
        max_ngrare_sd = which.max(per_state_pois_nrare_sd)
        min_ngrare_sd = which.min(per_state_pois_nrare_sd)

        # Variance of sample size n_{1} - overall rare population
        overall_var_nrare = sum(per_state_pois_nrare_vars)
        overall_sd_nrare = sqrt(overall_var_nrare)

        # Var of per-state sampling fraction (with Taylor approximation)
        per_state_pois_yield_frac_var = per_state_pois_nrare_vars / Tg^2 + per_state_pois_yield_count^2/Tg^4 * per_state_pois_ng_vars
        per_state_pois_yield_frac_sd = sqrt(per_state_pois_yield_frac_var)

        # Var of overall sampling fraction (with Taylor approximation)
        term1 = overall_var_nrare/Target^2
        term2 = yield_count^2/Target^4 * SD_n^2
        pois_yield_frac_sd = sqrt(term1 + term2)

        list("range_ng_sd" = range_ng_sd,
             "max_ng_sd" = state_labels[max_ng_sd],
             "min_ng_sd" = state_labels[min_ng_sd],
             "SD_n" = SD_n,
             "range_ngrare_sd" = range_ngrare_sd,
             "max_ngrare_sd" = state_labels[max_ngrare_sd],
             "min_ngrare_sd" = state_labels[min_ngrare_sd],
             "overall_sd_nrare" = overall_sd_nrare,
             "pois_yield_frac" = yield_frac,
             "pois_yield_frac_sd" = pois_yield_frac_sd, #-----------------------
             "per_state_pois_ng_vars" = per_state_pois_ng_vars,
             "per_state_pois_nrare_vars" = per_state_pois_nrare_vars,
             "per_state_pois_yield_count" = per_state_pois_yield_count,
             "per_state_pois_yield_frac" = per_state_pois_yield_frac,
             "per_state_pois_yield_frac_sd" = per_state_pois_yield_frac_sd
        )
}


# Sensitivity calculation 1
# Poisson WITHOUT filtering, for unobs. surnames P(R=1|s,g) set to be and actually is small
out = get_pois_yield_analysis(
        state_labels = df$state,
        Tg = df$target_pois_without_filter,
        prob_total = df$prob_totals_without_filter,
        prob_squared_total = df$prob_squared_totals_without_filter,
        prob_4thpower_totals = df$prob_4thpower_totals_without_filter
)
print(out)


# Sensitivity calculation 2
# Poisson WITHOUT filtering, for unobs. surnames P(R=1|s,g) set to be large BUT actually is small
out = get_pois_yield_analysis(
    state_labels = df$state,
    Tg = df$target_pois_without_filter2,
    prob_total = df$prob_totals_without_filter2,
    prob_squared_total = df$prob_squared_totals_without_filter2,
    prob_4thpower_totals = df$prob_4thpower_totals_without_filter2
)
print(out)


# Our actual calculation
# Poisson WITH filtering
out = get_pois_yield_analysis(
        state_labels = df$state,
        Tg = df$target_pois_filter,
        prob_total = df$prob_totals,
        prob_squared_total = df$prob_squared_totals,
        prob_4thpower_totals = df$prob_4thpower_totals
)
print(out)


# Add by-state elements to dataset for use in plotting below
df$per_state_pois_ng_vars = out$per_state_pois_ng_vars
df$per_state_pois_nrare_vars = out$per_state_pois_nrare_vars
df$per_state_pois_yield_frac = out$per_state_pois_yield_frac
df$per_state_pois_yield_frac_sd = out$per_state_pois_yield_frac_sd
df$per_state_pois_yield_count = out$per_state_pois_yield_count


# Add some elements to dataset for use in plotting below


# Tg = df$target_pois_filter
#
# # Expectation of frac of rare population sampled - by state
# df$per_state_pois_yield_count = (Tg / df$prob_totals  * df$prob_squared_totals)
# df$per_state_pois_yield_frac = df$per_state_pois_yield_count/Tg
#
# # Expectation of frac of rare population sampled - overall
# yield_count5 = sum(Tg/df$prob_totals * df$prob_squared_totals)
# yield_frac5 = yield_count5/Target
#
# # Variance of sample size n_g  by state
# per_state_var_term2 = (Tg^2 / df$prob_totals^2) * df$prob_squared_totals
# df$per_state_pois_ng_vars = Tg - per_state_var_term2
# df$per_state_pois_ng_sd = sqrt(df$per_state_pois_ng_vars)
# range(df$per_state_pois_ng_sd)
# which.max(df$per_state_pois_ng_sd)
# which.min(df$per_state_pois_ng_sd)
#
# # Variance of overall sample size n
# SD_n = sqrt(Target - sum(per_state_var_term2)) #Overall
#
#
# # Variance of sample size n_{g1} - by state rare population
# p1 = (Tg/df$prob_totals) * df$prob_squared_totals
# p2 = (Tg^2/df$prob_totals^2) * df$prob_4thpower_totals
# df$per_state_pois_nrare_vars = p1-p2
# df$per_state_pois_nrare_sd = sqrt(df$per_state_pois_nrare_vars)
#
# range(df$per_state_pois_nrare_sd)
# which.max(df$per_state_pois_nrare_sd)
# which.min(df$per_state_pois_nrare_sd)
#
# # Variance of sample size n_{1} - overall rare population
# overall_var_n_rare = sqrt(sum(p1-p2)) #Overall
#
# # Var of per-state sampling fraction (with Taylor approximation)
# per_state_vars_frac = df$per_state_pois_nrare_vars / Tg^2 + df$per_state_pois_yield_count^2/Tg^4 * df$per_state_pois_ng_vars
# df$per_state_pois_yield_frac_sd = sqrt(per_state_vars_frac)
#
# # Var of overall sampling fraction (with Taylor approximation)
# term1 = overall_var_nrare/Target^2
# term2 = yield_count5^2/Target^4 * SD_n^2
# SD5 = sqrt(term1 + term2)
#
# print(c(yield_frac5, SD5)*100)




################################################################################
################################################################################
# Plot
################################################################################
################################################################################
library(ggrepel) # for better text labels

################################################################################
# Plot comparing allocations for strat vs pois
################################################################################

df$prop_strat <- df$target_strat_filter / Target
df$prop_pois  <- df$target_pois_filter / Target
upper = max(df$prop_strat, df$prop_pois)
c_cut <- 0.018
q = ggplot(df, aes(x = prop_strat, y = prop_pois)) +
        geom_point(size = 3) +

        # 45-degree reference line
        geom_abline(intercept = 0, slope = 1,
                    linetype = "dashed",
                    color = "gray60") +
        # Labels only where x > c_cut
        geom_text_repel(
                data = subset(df, prop_strat > c_cut | prop_pois > c_cut),
                aes(label = state),
                size = 5,
                box.padding = 0.5, #how tight the labels pack
                point.padding = 0.7,
                force = 1.5,
                max.overlaps = Inf,
                direction = "both",
                hjust = 0,
                nudge_x = 0.003, #how far from point
                nudge_y = 0.003,
                segment.color = "gray" #remove lines by setting to NA
        ) +
        scale_x_continuous(limits = c(0, upper*1.05)) +
        scale_y_continuous(limits = c(0, upper*1.05)) +
        labs(
                title = "State Allocations",
                x = "Disproportionate stratified sampling",
                y = "Poisson sampling",
                tag = "B"
        ) +
        theme_minimal(base_size = 20) +
        theme(
                plot.title = element_text(face = "bold",
                                          hjust = 0.5,
                                          size = 28),
                axis.text.y = element_text(size = 24),
                axis.text.x = element_text(size = 24),
                axis.title.y = element_text(size = 24, margin = margin(t = 20)),
                axis.title.x = element_text(size = 24, margin = margin(t = 16)),
                legend.position = "bottom",
                plot.tag = element_text(face = "bold", size = 28,
                                        hjust = 0, vjust = 1),
                plot.tag.position = c(0, 1)
        )


################################################################################
# Select data needed for plot of yield calculations and create a separate df
# work on
################################################################################

plot_df = df[,c("state",
                    "p_R1_gvn_state",
                    "target_strat_filter",
                    "target_pois_filter",
                    "per_state_pois_yield_frac",
                    "per_state_pois_yield_frac_sd",
                    "per_state_strat_yield_frac_sd")]


################################################################################
# Create "All other" category and do calculations for it
################################################################################
# Create an "All other" category for 25 states with smallest P(R=1|g) values
mask = df$state %in% df$state[order(df$p_R1_gvn_state)[1:25]]

# Total allocation for all other category with stratified sampling
all_other_T_strat = sum(df$target_strat_filter[mask])
p_R1_AO = sum(df$p_R1_gvn_state[mask]*df$target_strat_filter[mask]/all_other_T_strat) ## P(R=1|all other)
all_other_strat_yield_frac_sd = sqrt(sum(df$per_state_strat_yield_frac_sd[mask]^2 * df$target_strat_filter[mask]^2)/all_other_T_strat^2)


# Total allocation for all other category with Poisson sampling
all_other_T_pois = sum(df$target_pois_filter[mask])
p_R1_A0_pois = sum(df$per_state_pois_yield_count[mask])/all_other_T_pois

term1 = sum(df$per_state_pois_nrare_vars[mask])/all_other_T_pois^2
term2 = sum(df$per_state_pois_yield_count[mask])^2/all_other_T_pois^4 * sum(df$per_state_pois_ng_vars[mask])
all_other_pois_yield_frac_sd = sqrt(term1+term2)

row = data.frame("state" = "All Others",
        "p_R1_gvn_state" = p_R1_AO,
        "target_strat_filter" = all_other_T_strat,
        "target_pois_filter" = all_other_T_pois,
        "per_state_pois_yield_frac" = p_R1_A0_pois,
        "per_state_pois_yield_frac_sd" = all_other_pois_yield_frac_sd,
        "per_state_strat_yield_frac_sd" = all_other_strat_yield_frac_sd,
        stringsAsFactors = FALSE
)

plot_df = rbind(plot_df[!mask, , drop= F], row)

################################################################################
# Bring in ACTUAL response rates
################################################################################
df_response = read.csv("../../data/SurveyResults/state_table.csv")

plot_df = left_join(plot_df, df_response, by = "state")

# Fill in the actual all others rate
others = df_response$state[!df_response$state %in% plot_df$state]
num_responded = sum(df_response$Responded[df_response$state %in% others], na.rm = T)
num_Jewish = sum(df_response$Jewish[df_response$state %in% others], na.rm = T)
plot_df$Jewish.rate[plot_df$state == "All Others"] = num_Jewish/num_responded



################################################################################


# Enforce order by difference in yield frac
plot_df = plot_df[order(plot_df$per_state_pois_yield_frac -  plot_df$p_R1_gvn_state, decreasing = F),]
plot_df$state <- factor(plot_df$state, levels = plot_df$state) #necessary to enforce order

top_state <- tail(levels(plot_df$state), 1)
top_df <- subset(plot_df, state == top_state)

p = ggplot(plot_df, aes(y = state)) +
        # x values (red circles)
        geom_point(aes(x = p_R1_gvn_state,
                       color = "Stratified",
                       shape = "Stratified"),
                   size = 3) +

        # y values (blue squares)
        geom_point(aes(x = per_state_pois_yield_frac,
                       color = "Poisson",
                       shape = "Poisson"),
                   size = 3) +

        # Axis labels
        labs(x = "Expected proportion Jewish",
             y = NULL,
             color = "Sampling Method",
             shape = "Sampling Method",
             title = "Estimated Success Rates",
             tag = "A"
             ) +

        # Set x-axis limits
        coord_cartesian(xlim = c(0, 1)) +

        # Manual color and shape scales
        scale_color_manual(
                breaks = c("Stratified", "Poisson"),
                values = c("Stratified" = old_method_color,
                           "Poisson" = new_method_color)
        ) +
        scale_shape_manual(
                breaks = c("Stratified", "Poisson"),
                values = c("Stratified" = 16,
                           "Poisson" = 15)
        ) +
        geom_errorbarh(
                aes(
                        xmin = pmax(per_state_pois_yield_frac - 2 * per_state_pois_yield_frac_sd, 0),
                        xmax = per_state_pois_yield_frac + 2 * per_state_pois_yield_frac_sd
                ),
                height = 0.25,        # controls vertical size of the little end caps
                color = new_method_color,
                alpha = 0.3,          # faded
                size = 1
        ) +
        geom_errorbarh(
                aes(
                        xmin = pmax(p_R1_gvn_state - 2 * per_state_strat_yield_frac_sd, 0),
                        xmax = p_R1_gvn_state + 2 * per_state_strat_yield_frac_sd
                ),
                height = 0.25,        # controls vertical size of the little end caps
                color = old_method_color,
                alpha = 0.3,          # faded
                size = 1
        ) +
        theme_minimal(base_size = 20) +
        theme(
                plot.title = element_text(face = "bold",
                                          hjust = 0.5,
                                          size = 28),
                axis.text.y = element_text(size = 20),
                axis.text.x = element_text(size = 24),
                axis.title.x = element_text(size = 24, margin = margin(t = 16)),
                legend.position = "none", #remove
                legend.text = element_text(size = 24),
                legend.title = element_text(size = 24),
                plot.tag = element_text(face = "bold", size = 28,
                                        hjust = 0, vjust = 1),
                plot.tag.position = c(0, 1)
        ) + geom_text(
            data = top_df,
            aes(x = p_R1_gvn_state,
                y = state,
                label = "Stratified"),
            color = old_method_color,
            hjust = -0.3,   # push text to the right
            size = 8
        ) +
            geom_text(
                data = top_df,
                aes(x = per_state_pois_yield_frac,
                    y = state,
                    label = "Poisson"),
                color = new_method_color,
                hjust = -0.4,
                size = 8
            )

print(p)


f1 = p + q + plot_layout(ncol = 2) +
        plot_annotation(
                #title = "Sampling Allocations and Success Rates by State",
                theme = theme(
                        plot.tag = element_text(size = 24, face = "bold"),
                        plot.title = element_text(
                                size = 28,
                                face = "bold",
                                hjust = 0.5
                        ),
                        plot.margin = margin(t = 15)   # add space ABOVE
                ),
        )

print(f1)

ggsave(
        filename = "../../paper/Figures/response_analysis/yield_and_allocation_comparison_plot.pdf",
        plot = f1,
        width = 20,
        height = 10,
        units = "in",
        dpi = 1200
)



# Save results
#---------------
write.csv(df, "../../data/cluster_output/probability_summaries_by_state_with_yield_calculations.csv")



################################################################################
################################################################################
################################################################################
################################################################################
# Version of plot that has the actual rates of success from survey +
# a calibration plot created from survey results
################################################################################
################################################################################
################################################################################
################################################################################

df_MA = plot_df[plot_df$state == "MA",] #3rd state with good locationr a label


g = p + geom_point(
    aes(x = Jewish.rate),
    shape = 4,        # X mark
    color = results_color,
    size = 3,
    stroke = 1.2      # makes the X thicker (important for visibility)
) + geom_text(
    data = df_MA,
    aes(x = Jewish.rate,
        y = state,
        label = "Observed"),
    color = results_color,
    hjust = -.2,   # push text to the right
    size = 8
)

print(g)


################################################################################
library(readxl)
df <- read_excel("../../data/SurveyResults/survey_for_analysis_100125_with_total_prob.xlsx")
# Note: has +6 rows for the people not linked...those people correspond to 6 of the non-respondent rows

df_obs = df[!is.na(df$jewish),]

# Note: 6 people not linkable
nrow(df_obs)
sum(is.na(df_obs$total_probs_bounded))
nrow(df_obs) - sum(is.na(df_obs$total_probs_bounded))

# Filter out those 6
df_obs = df_obs[!is.na(df_obs$total_probs_bounded),]
df_obs$total_probs_bounded = pmin(df_obs$total_probs_bounded, 1)



get_calibration_data <- function(
        df,
        prob_variable,
        jewish_variable,
        n_bins = 10,
        binning = c("quantile", "equal_width")
) {
    binning <- match.arg(binning)

    if (binning == "quantile") {
        breaks <- quantile(
            df[[prob_variable]],
            probs = seq(0, 1, length.out = n_bins + 1),
            na.rm = TRUE,
            type = 7
        )
        breaks <- unique(breaks)
    } else if (binning == "equal_width") {
        breaks <- seq(0, 1, length.out = n_bins + 1)
    }

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
                              n_bins = n_bins,
                              binning = "equal_width"
)

#TODO: check calculation
df_cal$se = sqrt(df_cal$frac_jewish * (1-df_cal$frac_jewish))/sqrt(df_cal$n)

#----------------------------------------------------
# Plot
#----------------------------------------------------
w <- ggplot(df_cal, aes(x = mean_prob, y = frac_jewish)) +
    geom_point(size = 4,
               shape = 17) +
    geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        color = "black",
        alpha = 0.7
    ) +
    geom_errorbar(
        aes(
            ymin = frac_jewish - 2 * se,
            ymax = frac_jewish + 2 * se
        ),
        width = 0.02,   # controls horizontal cap width
        color = "black",
        linewidth = 0.8
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
        x = glue("Mean predicted probability in bin"),
        y = "Fraction Jewish respondents",
        title = "Calibration Plot",
        tag = "B"
    ) +
    scale_x_continuous(
        limits = c(0, 1),
        expand = expansion(mult = c(0.02, 0.02)) # bit of padding but less than default
    ) +
    scale_y_continuous(
        limits = c(0, 1),
        expand = expansion(mult = c(0.02, 0.02)) # bit of padding but less than default
    ) +
    theme_minimal(base_size = 16) +
    theme(
        plot.title = element_text(face = "bold",
                                  hjust = 0.5,
                                  size = 28),

        axis.text.y = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 24, margin = margin(t = 16)),
        axis.title.y = element_text(
            size = 24,
            margin = margin(r = 20)  # increase space between ticks and label
        ),
        axis.line = element_blank(),   # remove axes
        legend.position = "none",
        plot.tag = element_text(face = "bold", size = 28,
                                hjust = 0, vjust = 1),
        plot.tag.position = c(0, 1),
    )



print(w)

f2 = g + w + plot_layout(ncol = 2) +
    plot_annotation(
        #title = "Sampling Allocations and Success Rates by State",
        theme = theme(
            plot.tag = element_text(size = 24, face = "bold"),
            plot.title = element_text(
                size = 28,
                face = "bold",
                hjust = 0.5
            ),
            plot.margin = margin(t = 15)   # add space ABOVE
        ),
    )
print(f2)



ggsave(
    filename = "../../paper/Figures/response_analysis/yield_comparison_plot_with_actual.pdf",
    plot = f2,
    width = 20,
    height = 10,
    units = "in",
    dpi = 1200
)
