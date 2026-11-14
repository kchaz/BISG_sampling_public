#
# Author: Kyla Chasalow
#
# Take initial first name file and do some cleaning
# Create an overall pooled file of first names
#

library(tidyverse)


count_tracker = list()

get_num_unique <- function(df){
  return(length(unique(df$first)))
}

# Script to extract cleaned first name lists from both original obit data
# and the scraped data

df = read_csv("../../data/ObitData/obit_list_092424-fulldehyph-check.csv")
df_scr = read_csv("../../data/ObitData/clean_first_names.csv") #scraped
df_scr = rename(df_scr, first = name)


count_tracker[["step 0 - obit data"]] = c(nrow(df),get_num_unique(df))
count_tracker[["step 0 - gpt data"]] = c(nrow(df_scr),get_num_unique(df_scr))


# Replace any of pattern[single letter][space][rest of letters] with just the part after the space
space_replacer <- function(firsts){
  s = gsub("^[a-zA-Z]\\s","", firsts)
  s = gsub("\\s[a-zA-Z]$","", s)
  return(s)}
df$first = space_replacer(df$first)
df_scr$first = space_replacer(df_scr$first)
stopifnot(sum(grepl("^[a-zA-Z]\\s", df$first))==0)
stopifnot(sum(grepl("^[a-zA-Z]\\s", df_scr$first))==0)

count_tracker[["step 1 - obit data"]] = c(nrow(df),get_num_unique(df))
count_tracker[["step 1 - gpt data"]] = c(nrow(df_scr),get_num_unique(df_scr))



# drop rows with no first name
df = df[!is.na(df$first),]
df_scr = df_scr[!is.na(df_scr$first),]

n_pre = nrow(df_scr)

# remove things that are not likely to actually be a first name
for (string in c("MOTHER","FATHER","STEPSON",
                 "DAUGHTER","^PARENT$",
                 "STEPFATHER", "STEPMOTHER", "STEP-",
                 "GRANDFATHER","GRANDPA","GRANDMA",
                 "^GREAT-", "^GREAT$",
                 "^UNCLE$","^AUNT$","^COUSIN$", "^COUSINS$", "^CANTOR$",
                 "^NEE[^a-zA-Z]$", "FIANCE","FIANCé","SPOUSE",
                 "^FRIEND", "^DOCTOR", "^PROFESSOR$",
                 "^MR$", "^MS$", "^MISS$","^MRS$","^JR$",
                 "^DR$", "^PROF$", "^SR$", "^JR$", "^MD$", #abbrev below suggestedby GPT
                 "^PHD$", "^REV$", "^HON$", "^ADM$",
                 "^GEN$", "^LT$", "^COL$", "^SGT$",
                 "^CPT$", "^ESQ$", "^ST$", "^FR$", "^PR$",
                 "^CEO$", "^CFO$", "^COO$", "^HUSBAND$",
                 "^DEVOTED$", "^CHERISHED$", "^LOVING$", "^DEAREST$",
                 "^DEARLY$", "^GRIEF$", "^MOURN",
                 "^ONE$", "^TWO$", "^THREE$", "^FOUR$", #Decision: numbers unlikely to be names
                 "^FIVE$", "^SIX$", "^SEVEN$", "^EIGHT$", "^NINE$", "^TEN$",
                 "^BRIEF$", "^SONS$"
                 )){
  mis = unique(df_scr$first[grepl(string, df_scr$first)])
  print("-------")
  print(string)
  print(mis)
  if (length(mis)>0){print(sum(df_scr$first %in% mis))}
  df_scr = df_scr[!df_scr$first %in% mis,]
}

count_tracker[["step 2 - obit data"]] = c(nrow(df),get_num_unique(df))
count_tracker[["step 2 - gpt data"]] = c(nrow(df_scr),get_num_unique(df_scr))


# Get rid of special characters
clean_characters <- function(df){

  #remove paranetheses
  df$first = gsub("\\(|\\)", "", df$first)

  #fix a / case that seems likely to be actually a full name
  df$first = gsub("PE/TER", "PETER", df$first)

  #for names split by /, take the first one (more likely to be official name? e.g., SORAH/SANDY)
  df$first[grepl("/", df$first)] = sapply(strsplit(df$first[grepl("/", df$first)],"/"), function(x){x[[1]]})

  #remove any & appearing purely at beginning or end
  df$first = gsub("^&", "", df$first)
  df$first = gsub("&$", "", df$first)

  #remove names that have any non-A-Z characters in it
  #TODO: might be able to improve this with converting things with apostrophes to not
  #but there are so few cases of that here (<10) that not bothering right now
  df = df[!grepl("[^a-zA-Z -]", df$first),]

  return(df)
}

df = clean_characters(df)
df_scr = clean_characters(df_scr)

count_tracker[["step 3 - obit data"]] = c(nrow(df),get_num_unique(df))
count_tracker[["step 3 - gpt data"]] = c(nrow(df_scr),get_num_unique(df_scr))



# drop rows with single character first name
single_character_dropper <- function(df){
  names = df$first
  nchars = sapply(names, nchar)
  to_remove = names[nchars <= 1]

  df = df[nchars > 1,]
  return(list("df" = df,
              "removed" = to_remove))
}
out = single_character_dropper(df)
df = out$df
out = single_character_dropper(df_scr)
df_scr = out$df


count_tracker[["step 4 - obit data"]] = c(nrow(df),get_num_unique(df))
count_tracker[["step 4 - gpt data"]] = c(nrow(df_scr),get_num_unique(df_scr))



#Decision: leaving names with hyphens. Could consider doing something with them


write.csv(df[,c("first","state")], "../../Data/ObitData/obit_list_092424_first_names_only.csv")
write.csv(df_scr, "../../Data/ObitData/clean_first_names-v2.csv")

#temp = read.csv("../../Data/ObitData/obit_list_092424_first_names_only.csv")
