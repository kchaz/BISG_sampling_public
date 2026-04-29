library(dplyr)


get_per_state_targets_strat <- function(Target, sampling_frame, p_R1_gvn_state){

  state_counts = table(sampling_frame$state)
  num = sqrt(p_R1_gvn_state[names(state_counts)]) * state_counts
  sample_frac = num/sum(num)
  Tg = round(sample_frac * Target)

  Tg = as.vector(Tg)
  names(Tg) = names(sample_frac)

  stopifnot(abs(sum(Tg) - Target) <= 10)

  return(Tg)
}


get_per_state_targets_pois <- function(Target, sampling_frame, p_state_gvn_R1){

  # Get prob counts and prob squared counts
  df = sampling_frame %>%
    group_by(state) %>%
    summarise(
      sum_prob = sum(p_R1_gvn_surname_state, na.rm = TRUE),
      sum_prob_sq = sum(p_R1_gvn_surname_state^2, na.rm = TRUE),
      .groups = "drop"
    )
  # Merge in p_state_gvn_R1
  p_state_gvn_R1 = data.frame("state" = names(p_state_gvn_R1),
                              "p_state_gvn_R1" = p_state_gvn_R1)
  df = left_join(df, p_state_gvn_R1, by = "state")

  # Allocation formula
  numerators = df$p_state_gvn_R1 / sqrt(df$sum_prob_sq / df$sum_prob)
  sample_fracs = numerators/sum(numerators)
  Tg = round(Target*sample_fracs)
  stopifnot(abs(sum(Tg) - Target) <= 10)

  names(Tg) = df$state

  return(Tg)
}





get_multipliers <- function(sampling_frame, group_column, probability_column, targets) {

  multiplier_df = sampling_frame %>%
    group_by(.data[[group_column]]) %>%
    summarise(probsum = sum(.data[[probability_column]], na.rm = TRUE), .groups = "drop")

  multiplier_df$multiplier = targets[multiplier_df[[group_column]]] / multiplier_df$probsum
  multiplier_df$probsum = NULL

  # if any are Inf, set to 0 - means probsum was already 0...so just keep them 0
  # Means will just not sample from that group
  multiplier_df$multiplier[multiplier_df$multiplier == Inf] = 0

  return(multiplier_df)
}



get_rescaled_sampling_probs <- function(sampling_frame, group_column, group_targets,
                                        probability_column, new_prob_column_name,
                                        remove_multiplier_column = T,
                                        fix_greater_than_1_probs = F){

  stopifnot(group_column %in% colnames(sampling_frame))

  # Get multipliers for rescaling
  multipliers = get_multipliers(sampling_frame, group_column, probability_column, group_targets)

  # Join multipliers temporarily to data frame
  sampling_frame = left_join(sampling_frame, multipliers, by = group_column)

  # Apply multipliers
  sampling_frame[[new_prob_column_name]] = sampling_frame[[probability_column]] * sampling_frame[["multiplier"]]

  # Remove unneeded column
  if (remove_multiplier_column){
    sampling_frame[["multiplier"]] = NULL
  }
  sampling_frame = as.data.frame(sampling_frame)

  #Including this option but warning: this could mask issues like being further off of target T
  if(fix_greater_than_1_probs){

    # set to max value less than 1
    max_val_lq_1 = max(sampling_frame[[new_prob_column_name]][sampling_frame[[new_prob_column_name]] < 1])
    sampling_frame[[new_prob_column_name]][sampling_frame[[new_prob_column_name]] >= 1] <- max_val_lq_1
  }


  return(sampling_frame)
}



# q can be a vector - creates a column for each
get_filtered_sampling_probs <- function(sampling_frame,
                                        group_column,
                                        group_targets,
                                        probability_column,
                                        qvals){

  stopifnot(all(qvals >= 0 & qvals <= 1))

  cutoffs = quantile(sampling_frame[[probability_column]], qvals)

  for (i in seq_along(qvals)){

    c = cutoffs[i]

    colname = glue("filtered_prob_{qvals[i]}")
    sampling_frame[[colname]] = sampling_frame[[probability_column]]
    sampling_frame[sampling_frame[[probability_column]] <= c, colname] = 0

    sampling_frame = get_rescaled_sampling_probs(sampling_frame,
                                                 group_targets = group_targets,
                                                 group_column = group_column,
                                                 probability_column = colname,
                                                 new_prob_column_name = colname
    )
  }

  return(sampling_frame)

}


#--------------------
make_precise_prob <- function(sampling_frame, u, Target, group_column, group_targets) {

  l = 1-u

  colname <- glue("precise_prob_{u}")
  sampling_frame[[colname]] <- NA_real_

  # For rare_status = 1
  idx1 <- which(sampling_frame$rare_status == 1)
  sampling_frame[[colname]][idx1] <- ifelse(runif(length(idx1)) < u, u, l)

  # For rare_status = 0
  idx0 <- which(sampling_frame$rare_status == 0)
  sampling_frame[[colname]][idx0] <- ifelse(runif(length(idx0)) < u, l, u)

  sampling_frame = get_rescaled_sampling_probs(sampling_frame,
                                               group_targets = group_targets,
                                               group_column = group_column,
                                               probability_column = colname,
                                               new_prob_column_name = colname
  )

  return(sampling_frame)
}




run_checks <- function(sampling_frame, Target, prob_columns, tol = 5){

  for (col in prob_columns){
    if (!col %in% colnames(sampling_frame)){
      stop(glue("Column {col} not present in sampling frame"))
    }
    if (any(is.na(sampling_frame[col]))){
      stop(glue("{col} contains NAs"))
    }
    if (abs(sum(sampling_frame[col]) - Target) > tol){
      stop(glue("Sum of {col} is >{tol} away from target"))
    }
    if(any(sampling_frame[col] < 0)){
      stop(glue("Column {col} has negative values"))
    }
    if(any(sampling_frame[col] > 1)){
      stop(glue("Column {col} has values > 1"))
    }
  }

}
