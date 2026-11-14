#'
#' Kyla Chasalow
#'
#' FILE MANAGEMENT
#'
#' This script gets P(S|G) from the voter file outputs
#' Note that for this it uses an absolute path to these files on the cluster
#'
#' This script uses P(R=1|G) and P(R=0|G) upper and lower bounds from the AJPP estimates.
#' Note that for this it uses a relative path to where these are supposed to be stored
#' in the Data folder of this repo. These can be re-generated using
#' clean_county_data.R in /src/1_data_processing/ but after running that script,
#' are stored in Data/DataObjects
#
#' The purpose of this script is to combine these two objects
#' to calculate bounds on each $P(S|G,R=1)$
#'
#' These are also stored on the cluster
#'
library(tidyverse)
library(glue)

# Load Brandeis files
p_R0_gvn_state = readRDS(fs::path("../../../Data/DataObjects/p_R0_gvn_state.RDS"))
p_R1_gvn_state = readRDS(fs::path("../../../Data/DataObjects/p_R1_gvn_state.RDS"))

p_R0_gvn_state_lower = readRDS(fs::path("../../../Data/DataObjects/p_R0_gvn_state_lower.RDS"))
p_R0_gvn_state_upper = readRDS(fs::path("../../../Data/DataObjects/p_R0_gvn_state_upper.RDS"))

p_R1_gvn_state_lower = readRDS(fs::path("../../../Data/DataObjects/p_R1_gvn_state_lower.RDS"))
p_R1_gvn_state_upper = readRDS(fs::path("../../../Data/DataObjects/p_R1_gvn_state_upper.RDS"))
print("Loaded Brandeis files")


state_surname_directory = fs::path("/n/imai_lab/Lab/lroyden/files/jewish_probs/state_vf_with_surname_only_jprobs")
#state_surname_directory = "temp"
state_csv_files <- list.files(state_surname_directory,
                              pattern = "^vf_surname_jprobs_[A-Za-z]{2}\\.csv$",
                              full.names = TRUE)
nf = length(state_csv_files)
print("Loaded state surname directories")
print(glue("Files to process:{nf}"))


# Place to save bounds
bound_directory = fs::path("/n/imai_lab/Lab/lroyden/files/state_obitdata_surname_probs/bounds/")
#bound_directory = "temp"
fs::dir_create(bound_directory)


# create a holder for only the surnames that are in the obituary data will fill it in state by state
obit_path = fs::path("/n/imai_lab/Lab/lroyden/files/original_obituary_data/obit_list_092424-fulldehyph-check.csv")
#obit_path = fs::path("../../../Data/ObitData/obit_list_092424-fulldehyph-check.csv")
df_obit = read.csv(obit_path)
lnames = sort(unique(df_obit$last))
holder = data.frame("surname" = lnames)
print(glue("Obtained obituary names:{length(lnames)}"))

for(file in state_csv_files){

  df = read.csv(file)
  state <- stringr::str_extract(file, "[A-Za-z]{2}(?=\\.csv$)")
  print(state)
  lowname = glue("lower_{state}")
  upname = glue("upper_{state}")
  lowname_c = glue("lower_{state}_c")
  upname_c = glue("upper_{state}_c")

  # Fixing possible encoding issues using iconv
  # df$surname <- iconv(df$surname, from = "latin1", to = "UTF-8", sub = "byte")

  # capitalize for consistency with other files
  df$surname = toupper(df$surname)

  # Non-CONSERVATIVE VERSION of bounds that will definitely fix > 1 problem
  df[,lowname] = pmax((df$denom - p_R0_gvn_state[state])/p_R1_gvn_state[state],0)
  df[,upname]  = pmin(df$denom / p_R1_gvn_state[state],1)
  stopifnot(all(df[,upname] >= df[,lowname]))
  stopifnot(all(df[,upname] >= df$denom))

    # CONSERVATIVE VERSION of bounds that may not fix > 1 problem
  df[,lowname_c] = pmax((df$denom - p_R0_gvn_state_upper[state])/p_R1_gvn_state_upper[state],0)
  df[,upname_c]  = pmin(df$denom / p_R1_gvn_state_lower[state],1)
  stopifnot(all(df[,upname_c] >= df[,lowname_c]))
  stopifnot(all(df[,upname_c] >= df$denom))

  dfsub <- df %>%
    filter(surname %in% lnames) %>%
    select(surname, .data[[lowname]], .data[[upname]], .data[[lowname_c]], .data[[upname_c]])

  holder = left_join(holder, dfsub, by = "surname")

  # If a surname does not occur in that state in the voter file, get a NA here. This is ok.

}

write.csv(holder,
          file.path(bound_directory, "surname_prob_bounds.csv"),
          row.names = F)
print("saved file")


# Check if any value in columns starting with "lower_" is greater than 0
any_above_zero <- holder %>%
  select(starts_with("lower_")) %>%
  summarise(across(everything(), ~ any(. > 0))) %>%
  unlist() %>%
  any()

# Count the number of values equal to 1 in columns starting with "upper_"
count_upper_ones <- holder %>%
  select(starts_with("upper_")) %>%
  summarise(across(everything(), ~ sum(. == 1, na.rm = TRUE))) %>%
  unlist() %>%
  sum()

# Count the number of values greater than or equal to 0.5 in columns starting with "upper_"
count_high <- holder %>%
  select(starts_with("upper_")) %>%
  summarise(across(everything(), ~ sum(. >= 0.5, na.rm = TRUE))) %>%
  unlist() %>%
  sum()

count_mediumhigh <- holder %>%
  select(starts_with("upper_")) %>%
  summarise(across(everything(), ~ sum(. >= 0.25, na.rm = TRUE))) %>%
  unlist() %>%
  sum()


print(glue("Lower bounds > 0? {any_above_zero}"))
print(glue("Upper bounds >= 0.5? {count_high}"))
print(glue("Upper bounds >= 0.25? {count_mediumhigh}"))
print(glue("Upper bounds = 1? {count_upper_ones}"))
