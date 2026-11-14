# developing the process for reshaping the theta matrix into a table
# and merging the P(R=1|G) probs with it
library(tidyverse)
source("paths.R")

#for saving some info about process
stats = list()

##########################################################################
# Load the final file
###########################################################################
#get filename (should only be one that matches this description)
file_list <- list.files(save_sampler_output_path,
                        pattern = "_realdf_45000\\.RDS$",
                        full.names = TRUE)
filename = file_list[1]
out = readRDS(filename)
timestamp <- gsub(".*/(\\d{14})_.*", "\\1", filename) #for saving later


##########################################################################
#' Load AJPP data derivatives
##########################################################################
# see generate_fips_state_proportions.R for code to generate these
pR1 = readRDS(paste0(data_object_path,"pR1.RDS"))
p_R1_gvn_state = readRDS(paste0(data_object_path,"p_R1_gvn_state.RDS"))
p_R0_gvn_state = readRDS(paste0(data_object_path,"p_R0_gvn_state.RDS"))


# Put num and denom together to get posterior means and reshape output matrix to data frame form
df = as.data.frame(out$theta_sums / out$sums_denom)

check = all(round(apply(df, 2, sum), 8) == 1)
if (!check){stop("Error: Not all columns of theta matrix sum to 1")}

df$surname = rownames(df)
long_df = tidyr::pivot_longer(df,
                              cols = -surname,
                              names_to = "state",
                              values_to = "p_surname_gvn_state_R1")

# Turn P(R=1|G) into a data frame
df_geography = data.frame("state" = names(p_R1_gvn_state),
                          "p_R1_gvn_state" = p_R1_gvn_state)

# Join the dataframe and calculate product for numerator
long_df = left_join(long_df, df_geography, by = "state")
long_df$numerator = long_df$p_surname_gvn_state_R1 * long_df$p_R1_gvn_state

# Some summary info for quick access
stats[["min_p_surname_gvn_state_R1"]] = min(long_df$p_surname_gvn_state_R1)
stats[["max_p_surname_gvn_state_R1"]] = max(long_df$p_surname_gvn_state_R1)
stats[["which.min_p_surname_gvn_state_R1"]] = long_df[which.min(long_df$p_surname_gvn_state_R1), ]
stats[["which.max_p_surname_gvn_state_R1"]] = long_df[which.max(long_df$p_surname_gvn_state_R1), ]

stats[["min_numerator"]] = min(long_df$numerator)
stats[["max_numerator"]] = max(long_df$numerator)
stats[["which.min_numerator"]] = long_df[which.min(long_df$numerator), ]
stats[["which.max_numerator"]] = long_df[which.max(long_df$numerator), ]
print(stats) #so that gets added to log
save_path = file.path(save_sampler_output_path, paste0(timestamp,"_basic_stats.RDS"))
saveRDS(stats, file = save_path)


# Sort by surname and then state
long_df <- arrange(long_df, surname)
long_df <- arrange(long_df, state)


save_dir = file.path(save_sampler_output_path, "state_obitdata_surname_probs")
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
  message("Directory created: ", save_dir)
}

# Save a file for each state
states = unique(long_df$state)
for (s in states){
  sub = filter(long_df, state == s) %>%
    select(surname, state, p_surname_gvn_state_R1, p_R1_gvn_state, numerator) #remove other quantities to save space
  save_path = file.path(save_dir, paste0(s,"_",timestamp,".csv"))
  write.csv(sub, save_path)
}
