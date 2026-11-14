#
# Author: Kyla Chasalow
#
# Take initial obit file and do some cleaning related to hypens
#
###############################
# Load Data
###############################
library(tidyverse)
library(readr)
df = read_csv("../../data/ObitData/obit_list_092424.csv")


##############################
# Hyphenated Names
##############################
hyphenated = grep("-",df$last)
df.hyphenated = df[hyphenated,]
num.hyphens = sapply(df.hyphenated$last, function(name){length(gregexpr("-", name)[[1]])})
counts = table(df.hyphenated$last)

flattened = unlist(strsplit(unique(df.hyphenated$last),"-"))
counts.flattened = table(flattened)
unique_flattened = unique(flattened)

unique_last = unique(df$last)

#find overlap - those in flattened that are in last
overlaps = unique_flattened[unique_flattened %in% unique_last]


summary = list()
summary[["num.unique.hyphenated.names"]] = length(unique(df.hyphenated$last))
summary[["num.hyphenated.names"]] = sum(counts)
summary[["number.of.hyphens.table"]] = table(num.hyphens)
summary[["max.hyphenated.name.count"]] = max(counts)
summary[["num.flattened.hyphen.componants"]] = length(flattened)
summary[["num.unique.flattened.hyphen.componants"]] = length(unique(flattened))
summary[["num.componants.in.main"]] = sum(unique_flattened %in% unique_last)
summary[["num.componants.not.in.main"]] = sum(!unique_flattened %in% unique_last)


#print(unique_flattened[unique_flattened %in% unique_last])
#print(unique_flattened[!unique_flattened %in% unique_last])

# examine counts of those with hyphenated names
last.counts = table(df$last)
ind = which(unique_flattened %in% unique_last)
count.componants = last.counts[unique_flattened[ind]]
hist(count.componants)
sum(count.componants <= 5)


count.componants["BERCOVITZ"]
count.componants["GOTFRIED"]
count.componants["ELBIN"]
count.componants["BARCHAM"]
count.componants["ROSENS"]
filter(df, last == "BERCOVITZ")
df.hyphenated[grep("BERCOVITZ", df.hyphenated$last),]



# Code which, for AAA-BBB, creates AAA, and BBB,
# checks if they appear elsewhere in data and adds only those that do.
expand_row <- function(row, last_id, nc, check_overlap_status){

  #get last and split
  last = as.vector(row[[last_id]])
  last.vec = strsplit(last, "-")[[1]]

  # If true, check if at least one component
  # appears elsewhere and if not, keep original hyphenated version.
  # If false, all hyphens get split and all componants added
  if (check_overlap_status){
      status = last.vec %in% overlaps
      if(sum(status) == 0){return(row)}
      last.vec = last.vec[status]
  }
  nh = length(last.vec)
  # repeat # of hypens times
  d = matrix(rep(row[-last_id], nh),
             ncol = nc-1,
             nrow = nh, byrow = T)
  return(cbind(last = last.vec, d))
}

# options
check_overlap_status = T
remove_hyphenated_rows = T

c_names = colnames(df)
last_id = which(colnames(df) == "last")
nc = ncol(df)
df.hyphen.flattened = apply(df.hyphenated, 1, expand_row, last_id = last_id, nc = nc,
                            check_overlap_status = check_overlap_status)
df.hyphen.flattened = as.data.frame(do.call(rbind, df.hyphen.flattened))
colnames(df.hyphen.flattened) = colnames(df)

#how many hyphens left
leftovers = df.hyphen.flattened[grepl("-",df.hyphen.flattened$last),]
dim(leftovers)


# Remove hyphenated rows from original dataset
# (some may be kept as part of df.hyphen.flattened, depending on choices above)
if (remove_hyphenated_rows){
  df.dehyph = df[!grepl("-",df$last),]
} else{
  df.dehyph = df
}



# new dataset merged
df.new = rbind(df.dehyph, df.hyphen.flattened)
if (check_overlap_status & remove_hyphenated_rows){
  name = "fulldehyph-check.csv"
} else if (check_overlap_status & !remove_hyphenated_rows){
  name = "partdehyph-check.csv"
} else if (!check_overlap_status & remove_hyphenated_rows) {
  name = "fulldehyph-nocheck.csv"
} else if (!check_overlap_status & !remove_hyphenated_rows) {
  name = "partdehyph-nocheck.csv"
}
write.csv(df.new, file.path("../../data/ObitData/", paste0("obit_list_092424-", name)))


