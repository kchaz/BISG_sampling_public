library(readxl)
library(readr) # for loading data
library(gtools) # for rdirichlet
library(tidyverse)

##########################################################################
#' Get list of States and their abbreviations and get pre-gen list of names
##########################################################################
# df_states <- read_csv("../../Data/states.csv")
# df_states = df_states[order(df_states$Abbreviation),] #make sure alphabetical order by abbreviation, not state name!
# #write_csv(df_states, "../../Data/states.csv") #save so that this is the default elsewhere
#
# n_states <- dim(df_states)[1]
# n_names <- 100
# surnames <- readRDS(paste0("SimulationObjects/surnames10000.RDS"))
# surnames = surnames[1:n_names]


##########################################################################
#' COUNTY-JEWISH DATA
##########################################################################
# filename <- "../../DATA/JewishCountyCounts/brandeis_fips_pop_counts-proc-V2.csv"
# county_counts2 <- read_csv(filename)
# fips2state = readRDS("../../Data/DataObjects/fips2state.RDS")
#
# fips_codes <- as.character(county_counts2$fips)
#
# # see generate_fips_state_proportions.R for code to generate these
# pR1 = readRDS("../../Data/DataObjects/pR1.RDS")
# p_R1_gvn_fips = readRDS("../../Data/DataObjects/p_R1_gvn_fips.RDS")
# p_R0_gvn_fips = readRDS("../../Data/DataObjects/p_R0_gvn_fips.RDS")
# p_fips_gvn_R1 = readRDS("../../Data/DataObjects/p_fips_gvn_R1.RDS")
# p_fips_gvn_R0 = readRDS("../../Data/DataObjects/p_fips_gvn_R0.RDS")
#
# p_state_gvn_R1 = readRDS("../../Data/DataObjects/p_state_gvn_R1.RDS")
# p_state_gvn_R0 = readRDS("../../Data/DataObjects/p_state_gvn_R0.RDS")
# p_R1_gvn_state = readRDS("../../Data/DataObjects/p_R1_gvn_state.RDS")
# p_R0_gvn_state = readRDS("../../Data/DataObjects/p_R0_gvn_state.RDS")
#

##########################################################################
#' OPTIONAL: SPECIFY STATE SUBSET
#' Optionally, reduce number of states used (for simulation purposes,
#' so can work with lower dimensional datasets)
##########################################################################
#TODO: decide if want to write a function to add this back
# use_state_subset = F
#
# if (use_state_subset){
#   #get a subset of some high-prev, some low-prev, some medium-prev states
#   #and include FL and CA which have high state | R=1
#   ordered_states = names(sort(p_R1_gvn_state))
#   state_subset = unique(c(ordered_states[3:6],
#                    ordered_states[23:27],
#                    ordered_states[48:51],
#                    "CA","FL"))
#
#   df_states = filter(df_states, Abbreviation %in% state_subset)
#   n_states = length(state_subset)
#
#   county_counts2 = filter(county_counts2, state %in% state_subset)
#   fips_codes = as.character(county_counts2$fips)
#
#   p_R1_gvn_fips = p_R1_gvn_fips[as.character(fips_codes)]
#   p_R0_gvn_fips = 1 - p_R1_gvn_fips
#   p_fips_gvn_R1 = p_fips_gvn_R1[as.character(fips_codes)]
#   p_fips_gvn_R0 = p_fips_gvn_R0[as.character(fips_codes)]
#   p_fips_gvn_R1 = p_fips_gvn_R1/sum(p_fips_gvn_R1)
#   p_fips_gvn_R0 = p_fips_gvn_R0/sum(p_fips_gvn_R0)
#
#   p_R1_gvn_state = p_R1_gvn_state[state_subset]
#   p_R0_gvn_state = 1 - p_R1_gvn_state
#   p_state_gvn_R1 = p_state_gvn_R1[state_subset]
#   p_state_gvn_R0 = p_state_gvn_R0[state_subset]
#   p_state_gvn_R1 = p_state_gvn_R1/sum(p_state_gvn_R1)
#   p_state_gvn_R0 = p_state_gvn_R0/sum(p_state_gvn_R0)
#
#   saveRDS(state_subset,"StateSubset.RDS")
# }

##########################################################################
#' SURNAME-STATE DATA
#' Note: generating a surname distribution for each state, not for each FIPS
#' but prevalence will be based on FIPS in simulated data
##########################################################################

# Generate P(S|state,R=1) varying for each state (assume constant by FIPS).
get_pmat <- function(state_names, surnames, rate = 1){

  n_names = length(surnames)
  n_states = length(state_names)

  gamma <- rexp(n_names, rate = rate)
  pmat <- t(rdirichlet(n_states, gamma)) # each row is a name, each col a state
  colnames(pmat) <- state_names
  rownames(pmat) <- surnames
  return(pmat)
}


# Get P(S|R=1) using the real P(G|R=1)
get_p_surname_gvn_R1 <- function(pmat, p_state_gvn_R1){
  p <- apply(t(pmat) * p_state_gvn_R1, 2, sum) #not matrix mult, each row (state) gets multiplied by its P(G|R=1)
  p <- t(as.matrix(p)) # same format as others
  return(p)
}




#' ######################################
#' FUNCTION TO GENERATE DATA FRAME
#' ######################################
#' @ps - vector of probs for each surname OR matrix of probs for each surname-state combo
#' @n - total sample size to generate in the end
#' @probs - either vector of P(state|R=1) or P(fips|R=1)

# functions to get a count based on props
get_counts_by_state_probs <- function(n, props) {
  counts <- rmultinom(1, n, props)
  rownames(counts) <- names(props)
  return(counts)
}
get_counts_by_fips_props <- function(n, props) {
  stopifnot(length(props) == length(fips_codes))
  fips_counts <- rmultinom(1, n, props)

  counts <- data.frame(
    fips_count = fips_counts,
    state = county_counts2$state
  )
  # state_counts = data.frame(count = fips_counts, state = county_counts2$state) %>%
  #   group_by(state) %>%
  #   summarize(state_total = sum(count))
  # state_vec = state_counts$state_total
  # names(state_vec) = state_counts$state
  return(counts)
}
# test = get_counts_by_fips_props(n = 1000, props = p_fips_gvn_R1)

# MAIN FUNCTION
create_df <- function(pmat, n, props, fips_mode = T) {

  state_names = colnames(pmat)
  surnames = rownames(pmat)

  # 1. Draw num obs for each fips (else each state) using props
  if (fips_mode) {
    counts <- get_counts_by_fips_props(n, props)
  } else {
    counts <- get_counts_by_state_probs(n = n, props = props)
  }

  # 2. Draw S count according to P(S|G,R=1) for each count > 0
  if (dim(pmat)[1] == 1) { # single P(S|R=1) dist
    name_counts <- sapply(counts, rmultinom, n = 1, prob = ps)
  } else {
    if (fips_mode) {
      name_counts <- sapply(
        fips_codes,
        function(j) {
          return(rmultinom(1,
            size = counts[j, "fips_count"],
            prob = pmat[, counts[j, "state"]]
          )) # state-specific prob
        }
      )
      colnames(name_counts) <- rownames(counts)
    }
    if (!fips_mode) {
      name_counts <- sapply(
        state_names,
        function(j) {
          rmultinom(1, size = counts[j, ], prob = pmat[, j])
        }
      )
      colnames(name_counts) <- state_names
    }
  }

  # create these rows in the data using rep function applied to each col
  datlist <- apply(name_counts, 2, function(c) {
    if (sum(c) > 0) {
      rep(surnames, times = c)
    } else {
      NULL
    }
  })

  # create dataset from list
  df <- data.frame()
  if (fips_mode) {
    for (f in fips_codes) {
      if (!is.null(datlist[[f]])) {
        rows <- data.frame(
          fips = rep(f, sum(name_counts[, f])),
          surname = datlist[[f]],
          state = counts[f, "state"]
        )
        df <- rbind(df, rows)
      }
    }
  } else {
    for (s in state_names) {
      if (!is.null(datlist[[s]])) {
        rows <- data.frame(
          state = rep(s, sum(name_counts[, s])),
          surname = datlist[[s]]
        )
        df <- rbind(df, rows)
      }
    }
  }
  df
}


# Generate P(S|state,R=0) that diverges from P(S|state, R=1)
flip_probabilities <- function(pmat, sd, beta = 10){
  pmat_flip <- apply(pmat, 2, function(p) {
    e <- rnorm(length(p), sd = sd)
    #pnew <- exp(-beta*sqrt(p) + e) / sum(exp(-beta*sqrt(p) + e))
    pnew <- exp(-beta*p + e) / sum(exp(-beta*p + e))
  })
}

# Calculate P(R=1|S,fips) and P(R=1|S,state) using Bayes rule and the R=0 and R=1 values
get_p_R1_gvn_state_surname <- function(pmat, pmat_flip, p_R1_gvn_state, p_R0_gvn_state){

  # assumed orientation of pmat
  n_names = nrow(pmat)
  n_states = ncol(pmat)
  state_names = colnames(pmat)
  surnames = rownames(pmat)

  # checks
  stopifnot(n_states == length(p_R1_gvn_state))
  stopifnot(n_states == length(p_R0_gvn_state))

  # P(R=1|S,state)
  output = matrix(NA, nrow = n_names, ncol = n_states)
  colnames(output) = state_names
  rownames(output) = surnames
  for (st in state_names){
    num = pmat[,st]*p_R1_gvn_state[st]
    den = pmat[,st]*p_R1_gvn_state[st] + pmat_flip[,st]*p_R0_gvn_state[st]
    output[,st] <- num/den
  }
  return(output)
}


#Not incorporating fips right now
generate_sampling_frame <- function(pmat, pmat_flip, n1, pR1, p_state_gvn_R1,
                                    p_R1_gvn_state, p_R0_gvn_state){


  dfrare <- create_df(pmat, n = n1, props = p_state_gvn_R1, fips_mode = F)
  dfrare$rare_status = rep(1, n1)

  n0 = round(n1/pR1)
  dfnotrare <- create_df(pmat_flip, n = n0, props = p_state_gvn_R0, fips_mode = F)
  dfnotrare$rare_status = rep(0, n0)

  # create single dataframe
  sampling_frame = rbind(dfrare, dfnotrare)

  # join in the true probabilities of R1
  truep = get_p_R1_gvn_state_surname(pmat, pmat_flip, p_R1_gvn_state, p_R0_gvn_state)
  temp = data.frame(truep, surname = rownames(truep)) %>%
    gather(key = "state", value = "p_R1_gvn_surname_state", colnames(truep))
  sampling_frame = left_join(sampling_frame, temp)

  return(sampling_frame)

}













#TODO: add FIPs version to get_p_R1_gvn_state_surname
# uses assumption S \perp Fips | state, R  (which is built to hold here)
#------------------------------------------------------
# # Calculate overall P(R=1|S,fips)
# final <- matrix(NA, nrow = n_names, ncol = length(fips_codes))
# colnames(final) <- fips_codes
# rownames(final) <- surnames
# for (f in fips_codes) {
#   st <- fips2state[[f]]
#   num <- ps3mat[, st] * p_R1_gvn_fips[[f]]
#   den <- ps3mat[, st] * p_R1_gvn_fips[[f]] + pmat_flip[, st] * (1 - p_R1_gvn_fips[[f]])
#   final[, f] <- num / den
# }
#
# # sanity check - what are highest probability fips?
# print(final[which.max(final)])
# i <- which(apply(final >= final[which.max(final)]-.1, 1, sum) != 0)
# j <- which(apply(final >= final[which.max(final)]-.1, 2, sum) != 0)
# states <- unique(fips2state[names(j)])
# for (st in states) {
#   print(st)
#   print(ps3mat[i, st])
#   print(pmat_flip[i, st])
# }
# filter(county_counts2, fips %in% !!names(j))

#saveRDS(final, paste0(folder_path, "/p_R1_gvn_S_fips.RDS"))


# Could turn these matrices into dataframes to save but file sizes bigger!
#-------------------------------------------------------------------------
# df_final = final %>% data.frame(surname = rownames(final))
# colnames(df_final) = colnames(final)
# df_final = gather(df_final, key = "fips", value = "prob", colnames(final))
# write_csv(df_final, paste0(folder_path, "/p_R1_gvn_S_fips.csv"))
# df_final2 = final2 %>% data.frame(surname = rownames(final2)) %>%
#   gather(key = "state", value = "prob", colnames(final2))
# write_csv(df_final2,paste0(folder_path, "/p_R1_gvn_S_state.csv"))

