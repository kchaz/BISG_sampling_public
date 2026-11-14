library(tidyverse)
#'
#' Want to get the count of observations with surnames in the obit data
#' for each state
#'
#' NOTE: later realized this script redundant with
#'
#'
finalprob_dir <- fs::path("/n/imai_lab/Lab/lroyden/files/vf_fullnames_with_jprobs")
save_dir <- fs::path("/n/imai_lab/Lab/lroyden/paper")

# Load each file
state_csv_files <- list.files(finalprob_dir,
                              pattern = "^[A-Za-z]{2}\\.csv$",
                              full.names = TRUE)
print(state_csv_files)

nf = length(state_csv_files)
r = 10

names = rep(NA, nf)
info = matrix(NA, nrow = nf, ncol = 3)

for (i in 1:nf){
  print(i)

  file = state_csv_files[i]
  state <- stringr::str_extract(file, "[A-Za-z]{2}(?=\\.csv$)")

  df = read.csv(file)

  df <- df %>%
    group_by(Last) %>%
    summarise(
      Count = sum(Count, na.rm = TRUE),
      last_probs = max(last_probs, na.rm = TRUE), #max should not make a difference
      total_probs = max(total_probs, na.rm = TRUE),
      total_probs_bounded = max(total_probs_bounded, na.rm = TRUE)
    )

  info[i,1] = sum(df$Count[df$last_probs > 0]) #should only be non-zero if surname in obit data
  info[i,2] = sum(df$Count) #note: double checked this num is same as get from the raw full name count files
  info[i,3] = sum(round(df$total_probs,r) != round(df$total_probs_bounded,10))
  names[i] = state
}

out_df = data.frame(info)
out_df$state = names
colnames(out_df) = c("Count_SurnameObs_in_Obit",
                     "Count_SurnameObs_Total",
                     "Count_Surnames_with_Bounds_Applied", "State")
save_path = file.path(save_dir, "dist_by_state_of_obit_surnames.csv")
print(dim(out_df))
write.csv(out_df, save_path, row.names = F)
print("saved")

# See comparing_Brandeis_proportions.R for analysis
