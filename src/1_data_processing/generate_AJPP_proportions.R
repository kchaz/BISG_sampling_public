#
# Author: Kyla Chasalow
#
# Script to generate and save all relevant proportions
# of form Pr(R=1|state), Pr(state|R=1), Pr(fips|R=1),Pr(R=1|fips)
# and the R=0 versions

filename <- "../../data/JewishCountyCounts/brandeis_fips_pop_counts-proc.csv"
AJPP_data1 <- read_csv(filename)
filename <- "../../data/JewishCountyCounts/brandeis_fips_pop_counts-proc-V2.csv"
AJPP_data2 <- read_csv(filename)

fips_codes <- as.character(AJPP_data2$fips)

# Value approx equals the .024 mentioned here
# https://www.pewresearch.org/religion/2021/05/11/the-size-of-the-u-s-jewish-population/
pR1 <- sum(AJPP_data1$jewish_adults) / sum(AJPP_data1$all_adults)

# Pr(R=1|G = fips) and Pr(R=0|G=fips)
p_R1_gvn_fips <- AJPP_data2$jewish_adult_given_county
names(p_R1_gvn_fips) <- fips_codes
p_R0_gvn_fips = AJPP_data2$non_jewish_adult_given_county
names(p_R0_gvn_fips) = fips_codes

# Pr(G = fips |R=1) and Pr(R=0|G=fips)
p_fips_gvn_R1 <- AJPP_data2$county_given_jewish_adult
names(p_fips_gvn_R1) <- fips_codes
p_fips_gvn_R0 = AJPP_data2$county_given_non_jewish_adult
names(p_fips_gvn_R0) = fips_codes

# Pr(state|R=1) and Pr(state|R=0)
state_counts <- group_by(AJPP_data1, state) %>%
  summarize(
    jewish_adults_lower = sum(jewish_adults_lowci),
    jewish_adults = sum(jewish_adults),
    jewish_adults_upper = sum(jewish_adults_highci),
    #jewish_all = sum(jewish_all),                    #Note: not using b/c data do not give total pop...only "all_adults"?
    #jewish_all_upper = sum(jewish_allhighci),
    #jewish_all_lower = sum(jewish_alllowci),
    all_adults = sum(all_adults),
    non_jewish_adults = sum(all_adults) - sum(jewish_adults)
  )

# sort by alphabetical order
state_counts = state_counts[order(state_counts$state),]

#====================================================================
# calculate proportions for JEWISH | STATE and NOT JEWISH | STATE
#====================================================================
# R=1 case
p_R1_gvn_state <- state_counts$jewish_adults / state_counts$all_adults
names(p_R1_gvn_state) = state_counts$state

p_R1_gvn_state_upper <- state_counts$jewish_adults_upper / state_counts$all_adults
names(p_R1_gvn_state_upper) = state_counts$state

p_R1_gvn_state_lower <- state_counts$jewish_adults_lower / state_counts$all_adults
names(p_R1_gvn_state_lower) = state_counts$state

# R=0 case
p_R0_gvn_state <- 1 - p_R1_gvn_state
names(p_R0_gvn_state) = state_counts$state

p_R0_gvn_state_upper <- 1 - p_R1_gvn_state_lower
names(p_R0_gvn_state_upper) = state_counts$state

p_R0_gvn_state_lower <- 1 - p_R1_gvn_state_upper
names(p_R0_gvn_state) = state_counts$state

stopifnot(all(p_R1_gvn_state_upper >= p_R1_gvn_state))
stopifnot(all(p_R1_gvn_state_lower <= p_R1_gvn_state))
stopifnot(all(p_R0_gvn_state_upper >= p_R0_gvn_state))
stopifnot(all(p_R0_gvn_state_lower <= p_R0_gvn_state))
stopifnot(all(p_R1_gvn_state + p_R0_gvn_state == 1))

# calculate proportions for STATE | JEWISH and STATE | NOT JEWISH
#====================================================================
#Note: not calculating upper/lower for these right now because not using them

p_state_gvn_R1 <- state_counts$jewish_adults / sum(state_counts$jewish_adults)
names(p_state_gvn_R1) = state_counts$state

p_state_gvn_R0 <- state_counts$non_jewish_adults / sum(state_counts$non_jewish_adults)
names(p_state_gvn_R0) = state_counts$state

# lookup table for FIPS to State
state_lookup <- AJPP_data2$state
names(state_lookup) <- AJPP_data2$fips

# saving
saveRDS(pR1, "../../data/DataObjects/pR1.RDS")
saveRDS(p_R1_gvn_fips, "../../data/DataObjects/p_R1_gvn_fips.RDS")
saveRDS(p_R0_gvn_fips, "../../data/DataObjects/p_R0_gvn_fips.RDS") #this is 1 - (previous)
saveRDS(p_fips_gvn_R1, "../../data/DataObjects/p_fips_gvn_R1.RDS")
saveRDS(p_fips_gvn_R0, "../../data/DataObjects/p_fips_gvn_R0.RDS")
saveRDS(p_state_gvn_R1, "../../data/DataObjects/p_state_gvn_R1.RDS")
saveRDS(p_state_gvn_R0, "../../data/DataObjects/p_state_gvn_R0.RDS")
saveRDS(p_R1_gvn_state, "../../data/DataObjects/p_R1_gvn_state.RDS")
saveRDS(p_R0_gvn_state, "../../data/DataObjects/p_R0_gvn_state.RDS")
saveRDS(state_lookup, "../../data/DataObjects/fips2state.RDS")

# saving the upper and lower bounds for P(R|G)
saveRDS(p_R1_gvn_state_upper, "../../data/DataObjects/p_R1_gvn_state_upper.RDS")
saveRDS(p_R1_gvn_state_lower, "../../data/DataObjects/p_R1_gvn_state_lower.RDS")
saveRDS(p_R0_gvn_state_upper, "../../data/DataObjects/p_R0_gvn_state_upper.RDS")
saveRDS(p_R0_gvn_state_lower, "../../data/DataObjects/p_R0_gvn_state_lower.RDS")
