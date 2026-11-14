#
# Author: Kyla Chasalow
#
# Load AJPP file and do some cleaning and checks
#

library(tidyverse)


# Load original data
original_filename <- "../../data/JewishCountyCounts/brandeis_fips_pop_counts.csv"
df <- read_csv(original_filename)

# Fix missing state entries
#------------------------------------------------------------
df[10, "state"] = "VA" #df[10,"countyname"] #manually checked
df[128, "state"] = "SD" #df[128, "countyname"] #manually checked

# Check # of states is as expected
states = unique(df$state)
paste0("Num unique states:",length(states))
stopifnot(length(states) == 51)


# Check that state named in county name is same as state in state column
#------------------------------------------------------------
test_consistency <- function(df){
  issues = 0
  for (s in states){
    present = grepl(s,df$countyname)
    if (!all(df$state[present] == s)){
      print(paste0("inconsistency for",s))
      issues = issues + 1
    }
  }
  if (issues == 0){print("no inconsistencies")}
}
test_consistency(df)


# Fix case where state name does not match - looked it up and it is in WI
# Redundant because remove these later but still want consistency first
present = grepl("MI",df$countyname)
which(df$state[present] == "WI")
sub = df[present,]
sub[15,"countyname"] #incorrectly labeled MI instead of WI
i = which(df$countyname == "Wausau & Nearby Counties, MI")
df[i,"countyname"] = "Wausau & Nearby Counties, WI"

# Recheck consistency
test_consistency(df)


# Something could do if wanted to isolate county name
# # Remove state name from county name since this is in other column
# for (s in states){
#   df$countyname = gsub(paste0(", ",s),"", df$countyname)
# }

# Are there cases with same exact name and state combination?
#------------------------------------------------------------
n = sum(table(df$state, df$countyname) > 1)
paste0("Number of non-unique state-county name pairs: ",n)
stopifnot(n==0)


# Missingness checks
#-----------------
miss_list = list()
for (col in colnames(df)){
  miss_list[[col]] = sum(is.na(df[,col]))
}
print(paste0("Rows missing state:",miss_list$state))
print(paste0("Rows missing countyname:",miss_list$countyname))
print(paste0("Rows missing jewish_adults:",miss_list$jewish_adults))
stopifnot(miss_list$state == 0)
stopifnot(miss_list$countyname == 0)
stopifnot(miss_list$jewish_adults == 0)
stopifnot(miss_list$fips1 == 0) #all at least one


# ADDITIONS
#------------

# Proportion Jewish Adults Columns
df$jewish_adult_given_county = df$jewish_adults/df$all_adults
df$county_given_jewish_adult = df$jewish_adults/sum(df$jewish_adults)

# same For Non-Jewish
df$non_jewish_adult_given_county = 1 - df$jewish_adult_given_county
df$county_given_non_jewish_adult = (df$all_adults - df$jewish_adults)/sum((df$all_adults - df$jewish_adults))

# Proportion of adult population represented by each FIPs
df$adult_prop = df$all_adults/sum(df$all_adults)

# How many FIPS codes associated with each row
fips_col = colnames(df)[grepl("fips",colnames(df))]
df$num_fips = rowSums(!is.na(df[,fips_col]))


# SAVE original file but with corrections and additions
#------------
new_filename <- sub("\\.csv$", "-proc.csv", original_filename)
write_csv(df, file = new_filename)


# Reshape data so there is just a FIPS columns
#----------------
df2 = gather(df, "fips_num_index","fips", fips_col)

#remove the NA's this creates
df2 = df2[!is.na(df2$fips),]

#check
stopifnot(sum(!is.na(df[,fips_col])) == dim(df2)[1])
stopifnot(length(unique(df2$fips)) == dim(df2)[1])


#drop fips_num_index - assume meaningless which is fips1 vs fips2 etc.
df2 = select(df2, -fips_num_index)


# save new file
new_filename <- sub("\\.csv$", "-proc-V2.csv", original_filename)
write_csv(df2, file = new_filename)




