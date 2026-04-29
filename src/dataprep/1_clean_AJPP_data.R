library(tidyverse)

# Load original AJPP data
original_filename <- "../../data/AJPP_data/brandeis_fips_pop_counts.csv"
df <- read_csv(original_filename)

# Fix missing state entries
df[10, "state"] = "VA" #manully checked using df[10,"countyname"]
df[128, "state"] = "SD" #manually checked using df[128, "countyname"]

# Check # of states is as expected
states = unique(df$state)
paste0("Num unique states: ",length(states))
stopifnot(length(states) == 51)

# Missing checks
miss_list = list()
for (col in colnames(df)){
  miss_list[[col]] = sum(is.na(df[,col]))
}
stopifnot(miss_list$state == 0)
stopifnot(miss_list$countyname == 0)
stopifnot(miss_list$jewish_adults == 0)
stopifnot(miss_list$fips1 == 0) #all at least one


#------------
# ADDITIONS
#------------

# Proportion Jewish Adults
df$jewish_adult_given_county = df$jewish_adults/df$all_adults
df$county_given_jewish_adult = df$jewish_adults/sum(df$jewish_adults)

# Proportion non-Jewish adults
df$non_jewish_adult_given_county = 1 - df$jewish_adult_given_county
df$county_given_non_jewish_adult = (df$all_adults - df$jewish_adults)/sum((df$all_adults - df$jewish_adults))

# Proportion of adult population
df$adult_prop = df$all_adults/sum(df$all_adults)

#------------------------------------------------------------
# SAVE with corrections and additions
#------------------------------------------------------------
new_filename <- sub("\\.csv$", "-proc.csv", original_filename)
write_csv(df, file = new_filename)
