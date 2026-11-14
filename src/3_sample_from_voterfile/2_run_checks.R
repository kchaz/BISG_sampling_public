#'
#' Run various kinds of data quality checks and do some basic EDA
#'
#' Relies on some absolute filepaths to locations on cluster
#'
#'
library(tidyverse)
library(glue)

surname_dir <- fs::path("/n/imai_lab/Lab/lroyden/files/jewish_probs/state_vf_with_surname_only_jprobs")
finalprob_dir <- fs::path("/n/imai_lab/Lab/lroyden/files/vf_fullnames_with_jprobs")

#surname_dir = "test/surnames"
#finalprob_dir = "test/final"


##############################################
# Check files that contain only the surname
# part of the final probability
##############################################

# Get all CSV files with exactly two-letter names
state_csv_files <- list.files(surname_dir,
                              pattern = "^vf_surname_jprobs_[A-Za-z]{2}\\.csv$",
                              full.names = TRUE)

print(glue("{length(state_csv_files)} files to check"))
problem_log = list()
EDA_log = list()

for (file in state_csv_files){

  # prep
  log = list()
  state <- stringr::str_extract(file, "[A-Za-z]{2}(?=\\.csv$)")
  df = read.csv(file)

  # check for NA and invalid probabilities
  for (col in c("num", "denom", "total", "total_bounded", "total_bounded_cnsv")){
    num_errors = sum(is.na(df[, col]))
    if (num_errors > 0) {
      log[glue("{col} has NA")] <- num_errors
    }
    num_errors = sum(df[, col] > 1 + 1e-10)
    if (num_errors > 0) {
      log[glue("{col}: has >1 values")] <- num_errors
      log[glue("{col}: max >1 magnitude")] <- max(df[df[, col] > 1e-10, col])
    }
    num_errors = sum(df[, col] < 0)
    if (num_errors > 0) {
      log[glue("{col} has <0 values")] <- num_errors
    }
  }

  # check that all records have a state
  num_errors = sum(is.na(df$state))
  if (num_errors > 0) {
    log[glue("Rows with state entry = NA")] <- num_errors
  }
  # Check that state is same as state stated in file name
  num_errors = sum(df$state != state, na.rm = T)
  if (num_errors > 0) {
    log[glue("Rows with state entry != file name state")] <- num_errors
  }

  # Add to overall problem_log log
  if (length(log)>0){
    problem_log[[state]] = log
  }

  ####################################
  # EDA
  ####################################
  log = list()
  # Calculate how often the dif total versions differ
  log["total_!=_total_bounded"] = mean(df$total != df$total_bounded)
  log["total_!=_total_bounded_cnsv"] = mean(df$total != df$total_bounded_cnsv)
  log["total_bounded_!=_total_bounded_cnsv"] = mean(df$total_bounded_cnsv != df$total_bounded)
  EDA_log[[state]] = log

}
saveRDS(EDA_log, file.path(surname_dir, "vf_surname_jprobs-EDA_log.RDS"))
print("saved EDA log")
if(length(problem_log) > 0){
  saveRDS(problem_log, file.path(surname_dir, "vf_surname_jprobs-problem_log.RDS"))
  print("saved problem log")
}



##############################################
# Check files that contain the final prob
##############################################
print("Checking join files")

# Get all CSV files with exactly two-letter names
state_csv_files <- list.files(finalprob_dir,
                              pattern = "^[A-Za-z]{2}\\.csv$",
                              full.names = TRUE)

print(glue("{length(state_csv_files)} files to check"))

problem_log = list()
EDA_log = list()


for (file in state_csv_files){

  # prep
  log = c()
  state <- stringr::str_extract(file, "[A-Za-z]{2}(?=\\.csv$)")
  df = read.csv(file)

  # Check that no names have a NA count or a 0 count (sanity check)
  if (sum(is.na(df$Count)) > 0){
    log["countsNA"] = "NA count rows present"
  }
  if (sum(df$Count == 0, na.rm = T) > 0){
    log["counts0"] = "0 count rows present"
  }

  # check common names are present with non-NA first name ratios
  for(name in c("EMILY","DANIEL","WILLIAM","SARAH", "JOHN")){
    if(!name %in% df$First){
      log[name] = "absent"
    } else{
      if(all(is.na(df[df$First == name, "first_ratio"])))
      log[name] = "present but first name ratio is NA"
    }
  }

  # Check for final probabilities > 1, <0
  # These could happen - not entirely errors...but want to know how prevalent
  for(col in c("total_probs", "total_probs_bounded", "total_probs_bounded_cnsv")){

    num_errors <- sum(df[,col] > 1+1e-10, na.rm = T)
    if (num_errors > 0) {
      log[[glue("{col}: probabilties > 1")]] = num_errors
      log[glue("{col}: max >1 magnitude")] <- max(df[df[, col] > 1e-10, col], na.rm =T)
    }
    num_errors <- sum(df$last_probs < 0, na.rm = T)
    if (num_errors > 0) {
      log[["probabilties < 0"]] = num_errors
    }
  }

  # add if there are problems
  if (length(log)>0){
    problem_log[[state]] = log
  }

  ########################
  # EDA
  ########################
  log = list()

  # how many ratios did not use a first name
  log["frac_NA_first_ratio"] = mean(is.na(df$first_ratio))

  for(col in c("total_probs", "total_probs_bounded", "total_probs_bounded_cnsv")){

    # which names have top 10 highest values of col
    top_indices <- order(df[,col], decreasing = TRUE)[1:10]
    top_names <- df[top_indices, c("First","Last",col)]
    log[[paste0("top10_names_", col)]] <- top_names

    # which names have lowest 10 values of col
    bottom_indices <- order(df[,col], decreasing = FALSE)[1:10]
    bottom_names <- df[bottom_indices, c("First","Last",col)]
    log[[paste0("bottom10_names_", col)]] <- bottom_names

  }
  EDA_log[[state]] = log
}

saveRDS(EDA_log, file.path(finalprob_dir, "vf_final_jprobs-EDA_log.RDS"))
print("saved EDA log")
if(length(problem_log) > 0){
  saveRDS(problem_log, file.path(finalprob_dir, "vf_finalfile-problem_log.RDS"))
  print("saved problem log")
}

print("All Checks Run")
