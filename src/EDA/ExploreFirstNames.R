#
# Author: Kyla Chasalow
#
library(tidyverse)
library(readr)
df = read_csv("../../data/ObitData/obit_list_092424-fulldehyph-check.csv")


##############################
# First names in obit data
##############################
unique_names <- unique(df$first)
length(unique_names)
first_names <- df$first

#there are some names that appear to still have an initial like S Barbara
space_names = first_names[grepl("\\s", first_names)]
length(space_names)

# Check for the pattern "[single letter][space][rest of letters]" and look at which
# do not meet this
pattern <- "^[a-zA-Z]\\s.*$"
matches <- grepl(pattern, space_names)
space_names[!matches]
space_names <- space_names[matches]
filtered_space_names <- gsub("^[a-zA-Z]\\s", "", space_names)

# For any with [single letter][space][rest of letters], replace with just the part
# after the space
df$first = sapply(df$first, function(name){
  present = grepl("^[a-zA-Z]\\s", name)
  if(!present){return(name)}

  name = gsub("^[a-zA-Z]\\s","", name)
  return(name)
})

# any names that are only of length 1 or 2?
nchars = nchar(df$first)
unique(df$first[nchars == 2])
unique(df$first[nchars <= 1]) #TODO: should remove these for name first name analysis


# Look at post-processed first names
unique_names2 <- unique(df$first)
length(unique_names2)
first_names2 <- df$first

#create table of first name counts
first.tab = table(first_names2)
first.tab = sort(first.tab, decreasing = T)
hist(first.tab)
sum(first.tab == 1)
sum(first.tab >= 2 & first.tab <= 10)
sum(first.tab > 10 & first.tab <= 100)
sum(first.tab > 100 & first.tab <= 1000)
sum(first.tab > 1000)


#####################################
# ChatGPT obit extraction dataset
#####################################
df_scr = read_csv("../../data/ObitData/clean_first_names.csv")


# drop rows with no first name
df_scr = df_scr[!is.na(df_scr$name),]

# drop rows with first name "-"
df_scr = df_scr[df_scr$name != "–",]
df_scr = df_scr[df_scr$name != "T+",]

# do any of the names have space in them? No (maybe extracted as separate)
space_names = df_scr$name[grepl("\\s", df_scr$name)]
length(space_names)

# mistake names
for (string in c("MOTHER","FATHER","STEPSON","DAUGHTER",
                 "GRANDFATHER","GRANDPA","GRANDMA",
                 "GREAT-", "STEP-",
                 "UNCLE","COUSIN$")){
  mis = unique(df_scr$name[grepl(string, df_scr$name)])
  print(mis)
  if (length(mis)>0){print(sum(df_scr$name %in% mis))}
  df_scr = df_scr[!df_scr$name %in% mis,]
}

# do any of the names have length <= 2
nchars = nchar(df_scr$name)
hist(nchars)
unique(df_scr$name[nchars <= 2])
sum(nchars <= 2)

# are any unusually long and if so, do they look ok?
unique(df_scr$name[nchars >= 14])

#what about hyphens
hyphenated = df_scr$name[grep("-",df_scr$name)]
counts = table(hyphenated)
counts[counts > 5]

# number of unique first names after cleaning
length(unique(df_scr$name))

# number of total obs
nrow(df_scr)

# distribution of name counts
counts = table(df_scr$name)
hist(counts, breaks = 100)

length(counts[counts>100])
length(counts[counts>1000])


##### other variables
# relationship variable is very messy
length(unique(df_scr$relationship))



#####################################
# Compare the two
#####################################
obit_names = df$first
scrape_names = df_scr$name
obit_names.u = unique(obit_names)
scrape_names.u = unique(scrape_names)

#how much overlap?
sum(obit_names.u %in% scrape_names.u)/length(obit_names.u)

#which do not?
sub = obit_names.u[!obit_names.u %in% scrape_names.u]
#count = table(obit_names[obit_names %in% sub])


#################
# Explore some particular names
sum(obit_names == "MILDRED")/length(obit_names)
sum(scrape_names == "MILDRED")/length(scrape_names)

sum(obit_names == "CALEB")/length(obit_names)
sum(scrape_names == "CALEB")/length(scrape_names)



