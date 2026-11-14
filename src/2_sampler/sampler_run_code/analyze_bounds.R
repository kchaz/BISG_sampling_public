#'
#' Author: Kyla Chasalow
#'
#' Small script for analyzing bounds I downloaded off cluster and put in
#' obitdata_outputs together with the sampler output I downloaded off cluster
#' Analyzing how often the bounds were actually needed
#'
# Libraries
library(readr)
library(glue)
library(tidyverse)

# Scripts
source("../../dataprep/data_processing_helpers.R")
source("../HBayesSampler/sampler_analysis_helpers.R")

out = readRDS("obitdata_outputs/20241127011108_realdf_45000.RDS")
Theta_means = out$theta_sums/out$sums_denom

#NAs just mean that surname-state combo wasn't observed in voter file
bounds = read.csv("obitdata_outputs/surname_prob_bounds.csv")

states = colnames(Theta_means)

num_bounded = rep(NA, length(states))
names(num_bounded) = states

upper_cols = rep(NA, length(states))
names(upper_cols) = states


# Just looking at strict, not conservative right now
for (state in states){

  low_name = glue("lower_{state}")
  high_name = glue("upper_{state}")

  if(any(bounds[,low_name]!=0, na.rm = T)){
    print(glue("{state} has lower bounds > 0"))
  }

  # get all info together
  bounds_sub = bounds[,c("surname",high_name)]
  df <- data.frame(surname = names(Theta_means[,state]), est_prob = Theta_means[,state], row.names = NULL)
  combo = left_join(bounds_sub, df, by = "surname")

  # which ones have val > bound
  num_bounded[state] = sum(combo$est_prob > combo[,high_name], na.rm = T)

  upper_cols[state] = high_name


}

print(num_bounded)
print(sum(num_bounded))
sum(num_bounded)/prod(dim(Theta_means))
sum(is.na(bounds[,upper_cols]))/prod(dim(Theta_means))
