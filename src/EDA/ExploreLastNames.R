#
# Author: Kyla Chasalow
#

library(tidyverse)
library(readr)
df = read_csv("../../data/ObitData/obit_list_092424-fulldehyph-check.csv")

# MISSINGNESS
apply(df, 2, function(x){sum(is.na(x))})

# STATES
state.counts <- table(df$state)
length(state.counts)
hist(state.counts)
sum(state.counts > 10)
sum(state.counts > 500)
sum(state.counts > 5000)
sum(state.counts > 10000)


# LAST NAMES
last.counts <- table(df$last)
length(last.counts)

just.counts <- sort(unname(last.counts), decreasing = T)

#how many names have >= x count?
#-----------------------------
sum(just.counts > 1000)
sum(just.counts > 100 & just.counts <= 1000)
sum(just.counts > 10 & just.counts <= 100)
sum(just.counts >= 2 & just.counts <= 10)
sum(just.counts == 1)

hist(just.counts,
     main = "Distribution of Counts for Distinct Last Names",
     xlab = "Last Name Count",
     ylab = "Frequency",
     breaks = 50)


which(last.counts > 1000)
which(last.counts > 100 & last.counts <= 1000)
which(last.counts > 10 & last.counts <= 100)
which(last.counts == 1)

#Uncommon names with common versions
#--------------------------------------
df %>% filter(last == "BERLOWITZ")
berwitz_df <- df %>% filter(grepl("BER", last) & grepl("WITZ", last))
table(berwitz_df$last)

stein_df <- df %>% filter(grepl("STEIN", last))
dim(stein_df)
sort(table(stein_df$last), decreasing = T)


# for names that occur more than once, how many states do they occur in
#----------------------------------------
last.multis <- last.counts[which(last.counts > 1)]
last.multis.names <- as.vector(names(last.multis))
multis_df <- df[df$last %in% last.multis.names,]

n_states_df <- multis_df %>%
  group_by(last) %>%
  summarize(num_states = n_distinct(state))
  #mutate(num_states = n_distinct(state)) %>%
  #ungroup()
n_states_df$count = as.vector(last.multis)

#uses the mutate() version to look at which states the
#one-state-only names come from
#table(n_states_df[n_states_df$num_states == 1,"state"])

#intuitively, the more occurrences, the more states occurred in
#but there are names that occur quite a few times and only in a few states
plot(n_states_df$count, n_states_df$num_states,
     xlab = "Num. Occurences of Last Name",
     ylab = "Num. States Occurred In")


#names that occur in only one state
one_state <- n_states_df %>% filter(num_states == 1)
one_state <- one_state[order(one_state$count, decreasing = T),]
#df %>% filter(last == "WEINGARDEN")
#df %>% filter(last == "ABELOW")
head(one_state,20)


#names that occur in only two states
two_state <- n_states_df %>% filter(num_states == 2)
two_state <- two_state[order(two_state$count, decreasing = T),]
df %>% filter(last == "GARCIA") %>% select(state)
df %>% filter(last == "RABINOVITZ") %>% select(state)
head(two_state,20)

