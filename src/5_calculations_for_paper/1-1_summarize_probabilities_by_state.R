#'
#'
#' NOTE: A subset of this code appears in long full code.R
#'
#' This script assumes that has already been run so that the sampling probs actually appear
#' in the vf_fullnames_with_jprobs and calculates some summaries for the paper
#'

library(tidyverse)
library(glue)

# For encoding bin counts below
encode_vec <- function(x) {
  paste(paste0(names(x), "=", x), collapse = "|")
}
decode_vec <- function(s) {
  parts <- strsplit(s, "\\|")[[1]]
  kv <- strsplit(parts, "=")
  vals <- as.numeric(sapply(kv, `[`, 2))
  names(vals) <- sapply(kv, `[`, 1)
  vals
}

#####################################################################
# Calculating Normalizing Denom and State probability totals
#####################################################################

finalprob_dir <- fs::path("/n/imai_lab/Lab/lroyden/files/vf_fullnames_with_jprobs")

# Get all CSV files with exactly two-letter names
state_csv_files <- list.files(finalprob_dir,
                              pattern = "^[A-Za-z]{2}\\.csv$",
                              full.names = TRUE)

nf = length(state_csv_files)
prob_matrix <- matrix(NA, nrow = nf, ncol = 13)

colnames(prob_matrix) <- c(
  "state",
  "prob_totals",
  "prob_squared_totals",
  "prob_4thpower_totals",
  "prob_min",
  "prob_max",
  "state_counts",
  "state_nonzeroprob_counts",
  "num_with_bounds_applied",
  "num_with_trunc_at_1",
  "total_probs_bounded_bins",
  "sampling_prob_bins"
)


for (i in 1:nf){

  file = state_csv_files[i]
  prob_matrix[i, "state"] = stringr::str_extract(file, "[A-Za-z]{2}(?=\\.csv$)")

  print(prob_matrix[i, "state"])
  df = read.csv(file)

  # Get rows with a non-NA obit data based probability and
  # make sure none missing total or a count
  sub = filter(df, !is.na(df$total_probs_bounded))
  stopifnot(sum(is.na(sub$total_probs_bounded)) == 0)
  stopifnot(sum(is.na(sub$Count)) == 0)
  if(nrow(sub)<nrow(df)){print(glue("{states[i]}: {nrow(df)-nrow(sub)} people with NA prob dropped"))}

  # Get Total of Probabilities truncating any > 1 (from first name ratio mult) as do in real calc - get \pi_g
  prob_matrix[i,"prob_totals"]  = sum(pmin(sub$total_probs_bounded,1) * sub$Count)

  # Get sum of squared probs
  prob_matrix[i,"prob_squared_totals"] = sum(pmin(sub$total_probs_bounded,1)^2 * sub$Count)

  # Get sum of 4th power of probs
  prob_matrix[i,"prob_4thpower_totals"]  = sum(pmin(sub$total_probs_bounded,1)^4 * sub$Count)

  # Get min prob >0
  prob_matrix[i,"prob_min"]  = min(sub$total_probs_bounded[sub$total_probs_bounded > 0])

  # Get max prob
  prob_matrix[i,"prob_max"]  = max(sub$total_probs_bounded)

    # Get total number of people
  prob_matrix[i,"state_counts"] = sum(df$Count) #includes NA people

  # Get total number of people with non-zero prob of being sampled
  prob_matrix[i,"state_nonzeroprob_counts"] = sum(sub$Count[sub$total_probs_bounded > 0])

  # Get number in which total_probs_bounded > 1 (b/c of first name ratio > 1). Can get frac using frame_count or total_count
  prob_matrix[i,"num_with_trunc_at_1"] = sum(sub$total_probs_bounded > 1)

  # Get binned break-down
  breaks <- seq(0, 1, by = 1/100)
  prob_matrix[i, "total_probs_bounded_bins"] <- encode_vec(table(
    cut(
      pmin(sub$total_probs_bounded, 1),
      breaks = breaks,
      include.lowest = TRUE,
      right = FALSE
    )
  ))
  prob_matrix[i,"sampling_prob_bins"] <- encode_vec(table(
    cut(
      sub$sampling_prob,
      breaks = breaks,
      include.lowest = TRUE,
      right = FALSE
    )
  ))

}


df_tot = data.frame(prob_matrix)


write.csv(df_tot,  "probability_summaries_by_state.csv",
          #file.path(finalprob_dir, "probability_summaries_by_state.csv"),
          row.names = F)
