#'
#' Author: Kyla Chasalow
#'
#'
#' FILE MANAGEMENT
#'
#' This script assumes that sampler has already been run and uses an absolute
#' filepath on cluster to access the files there
#
#' The purpose of this script is to impose bounds calculated using `get_bounds.R`
#' script to impose bounds on the sampler estimates, setting them to upper bound if they
#' are greater than upper and lower bound if they are less than lower
#'
#' The files are then re-saved with extra columns for the bound-corrected versions
#'
#'
library(tidyverse)
library(glue)

# Get directory where p(S|G,R=1) estimates are stored
state_surname_directory = fs::path("/n/imai_lab/Lab/lroyden/files/state_obitdata_surname_probs")
bound_directory = fs::path("/n/imai_lab/Lab/lroyden/files/state_obitdata_surname_probs/bounds/surname_prob_bounds.csv")

#state_surname_directory = "temp"
#bound_directory = "temp/surname_prob_bounds.csv"

state_csv_files <- list.files(state_surname_directory,
                              pattern = "^[A-Za-z]{2}_.*\\.csv$",
                              full.names = TRUE)
nf = length(state_csv_files)
print(glue("State files to run: {nf}"))

df_bounds = read.csv(bound_directory)
print("Loaded bounds")


for(file in state_csv_files){

  df = read.csv(file)
  state <- stringr::str_extract(fs::path_file(file), "^[A-Za-z]{2}")
  print(state)
  lowname = glue("lower_{state}")
  upname = glue("upper_{state}")
  lowname_c = glue("lower_{state}_c")
  upname_c = glue("upper_{state}_c")


  new <- left_join(
    select(df, surname, p_surname_gvn_state_R1),
    select(df_bounds, surname, all_of(c(upname, lowname, upname_c, lowname_c))),
    by = "surname"
  )


  # enforce  non-conservative bounds
  new$p_bounded = ifelse(new$p_surname_gvn_state_R1 > new[[upname]],
                                 new[[upname]], #if true, replaces by upper. If false, goes to lower check
                          ifelse(new$p_surname_gvn_state_R1 < new[[lowname]],
                                 new[[lowname]],
                                 new$p_surname_gvn_state_R1))
  new$p_bounded[is.na(new$p_bounded)] = new$p_surname_gvn_state_R1[is.na(new$p_bounded)]


  # checks and stats
  stopifnot(all(new$p_bounded <= new[,upname], na.rm = T))
  stopifnot(all(new$p_bounded >= new[,lowname], na.rm = T))
  print(glue("Non-conservative upper bounds enforced: {sum(new$p_bounded == new[,upname], na.rm = T)} of {nrow(df)}"))
  print(glue("Non-conservative lower bounds enforced: {sum(new$p_bounded == new[,lowname], na.rm = T)} of {nrow(df)}"))
  df$numerator_bounded = new$p_bounded * df$p_R1_gvn_state


  # enforce conservative bounds
  new$p_bounded_cnsv = ifelse(new$p_surname_gvn_state_R1 > new[[upname_c]],
                          new[[upname_c]], #if true, replaces by upper. If false, goes to lower check
                          ifelse(new$p_surname_gvn_state_R1 < new[[lowname_c]],
                                 new[[lowname_c]],
                                 new$p_surname_gvn_state_R1))
  new$p_bounded_cnsv[is.na(new$p_bounded_cnsv)] = new$p_surname_gvn_state_R1[is.na(new$p_bounded_cnsv)]
  stopifnot(all(new$p_bounded_cnsv <= new[,upname_c], na.rm = T))
  stopifnot(all(new$p_bounded_cnsv >= new[,lowname_c], na.rm = T))
  print(glue("Conservative upper bounds enforced: {sum(new$p_bounded_cnsv == new[,upname_c], na.rm = T)} of {nrow(df)}"))
  print(glue("Conservative lower bounds enforced: {sum(new$p_bounded_cnsv == new[,lowname_c], na.rm = T)} of {nrow(df)}"))
  df$numerator_bounded_cnsv = new$p_bounded_cnsv * df$p_R1_gvn_state

  # re-save file
  write.csv(df, file = fs::path(file), row.names = F)
  print(glue("Resaved file for {state}"))
}

