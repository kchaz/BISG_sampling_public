# Script to generate and save all relevant proportions
# of form Pr(R=r|State) and P(state|R=1) were R=1 is Jewish
# and R=0 is not Jewish

filename <- "../../data/AJPP_data/brandeis_fips_pop_counts-proc.csv"
county_counts1 <- read_csv(filename)

state_counts <- group_by(county_counts1, state) %>%
  summarize(
    jewish_adults_lower = sum(jewish_adults_lowci),
    jewish_adults = sum(jewish_adults),
    jewish_adults_upper = sum(jewish_adults_highci),
    all_adults = sum(all_adults),
    non_jewish_adults = sum(all_adults) - sum(jewish_adults)
  )

# sort by alphabetical order
state_counts = state_counts[order(state_counts$state),]

#' Value ~ the .024 mentioned here
#' https://www.pewresearch.org/religion/2021/05/11/the-size-of-the-u-s-jewish-population/
pR1 <- sum(county_counts1$jewish_adults) / sum(county_counts1$all_adults)

#====================================================================
# calculate proportions for JEWISH | STATE and NOT JEWISH | STATE
#====================================================================
# R=1 case
p_R1_gvn_state <- state_counts$jewish_adults / state_counts$all_adults
names(p_R1_gvn_state) = state_counts$state

# R=0 case
p_R0_gvn_state <- 1 - p_R1_gvn_state
names(p_R0_gvn_state) = state_counts$state

#checks
stopifnot(all(p_R1_gvn_state + p_R0_gvn_state == 1))

# calculate proportions for STATE | JEWISH and STATE | NOT JEWISH
#====================================================================
#Note: not calculating upper/lower for these right now because not using them

p_state_gvn_R1 <- state_counts$jewish_adults / sum(state_counts$jewish_adults)
names(p_state_gvn_R1) = state_counts$state

p_state_gvn_R0 <- state_counts$non_jewish_adults / sum(state_counts$non_jewish_adults)
names(p_state_gvn_R0) = state_counts$state


# saving
saveRDS(pR1, "../../data/AJPP_objects/pR1.RDS")
saveRDS(p_state_gvn_R1, "../../data/AJPP_objects/p_state_gvn_R1.RDS")
saveRDS(p_state_gvn_R0, "../../data/AJPP_objects/p_state_gvn_R0.RDS")
saveRDS(p_R1_gvn_state, "../../data/AJPP_objects/p_R1_gvn_state.RDS")
saveRDS(p_R0_gvn_state, "../../data/AJPP_objects/p_R0_gvn_state.RDS")
