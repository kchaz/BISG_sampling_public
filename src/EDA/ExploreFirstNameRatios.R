#
# Author: Kyla Chasalow
#

# This script is for exploring first name ratios
# As created on the cluster using ratio of proportion in
# obituary data to proportion in voter file. Larger proportion
# means more common in obituary data than voter file.
#
# Note: code for creating firstnames.csv is in voterfile/pi.R
library(tidyverse)

dfratio = read_csv("../../data/cluster_output/firstnames_jratio.csv")
n = nrow(dfratio)

# Old code when ratio was calculated for everyone
#-------------------------------------------------------------------
# sum(is.na(dfratio$ratio))/n
# sum(is.na(dfratio$ratio)) == sum(is.na(dfratio$prop_j)) + sum(is.na(dfratio$prop_vf))
#
# #majority of names are NA because do not appear in obit data. That is ok.
# sum(is.na(dfratio$prop_j))/n
# head(dfratio[is.na(dfratio$prop_j),"first"], 20)
#
# #some names do not appear in voter file despite appearing in obit data but this less common.
# sum(is.na(dfratio$prop_vf))
# head(dfratio[is.na(dfratio$prop_vf),"first"], 20) #see ACRHIE, misspelling of Archie probably
# tail(dfratio[is.na(dfratio$prop_vf),"first"], 20)


#how many names in dfratio$first have spaces in them? about 14%.
# possible these should have been processed when calculating ratios
# but not entirely clear that is right choice. Obit data names were cleaned up a bit
space_names = dfratio$first[grepl(" ", dfratio$first)]
length(space_names)/n

#how many have single letter, space, name pattern? these should be cleaned
# initial_names = dfratio$first[grepl("^[a-zA-Z]\\s", dfratio$first)]
# initial_names2 = dfratio$first[grepl("\\s[a-zA-Z]$", dfratio$first)]
# (length(initial_names) + length(initial_names2))/n

hyph_names = dfratio$first[grepl("-", dfratio$first)]
length(hyph_names)/n


# Prevalence
#------------------------------------------
mean(dfratio$ratio >= 1, na.rm = T)
mean(dfratio$ratio < 1, na.rm = T)

# Plot histogram of ratios
hist(log(dfratio$ratio[!is.na(dfratio$ratio)]),
     breaks=100, main="Histogram of First Name Ratios",
     xlab="Ratio of Proportions: Obit Data / Voter File")

sum(dfratio$ratio >= 500, na.rm = T)
sum(dfratio$ratio >= 100 & dfratio$ratio < 500, na.rm = T)
sum(dfratio$ratio >= 10 & dfratio$ratio < 100, na.rm = T)
sum(dfratio$ratio >= 1 & dfratio$ratio < 10, na.rm = T)
sum(dfratio$ratio <= 1 & dfratio$ratio > .1, na.rm = T)
sum(dfratio$ratio <= .1 & dfratio$ratio > .01, na.rm = T)
sum(dfratio$ratio <= .01 & dfratio$ratio > .001, na.rm = T)
sum(dfratio$ratio < .001, na.rm = T)

#-------------------------------------------
#  observing some mess like "DEVOTED", "HUSBAND", "CHERISHED" but those shouldn't cause
#  problems if do not actually occur in voter file much
dfratio$first[which(dfratio$ratio>500)]
dfratio$first[which(dfratio$ratio>250 & dfratio$ratio<=500)]
dfratio$first[which(dfratio$ratio>150 & dfratio$ratio<=250)]
dfratio$first[which(dfratio$ratio>125 & dfratio$ratio<=150)]
dfratio$first[which(dfratio$ratio>100 & dfratio$ratio<=125)]
dfratio$first[which(dfratio$ratio<1)]

# These seem to make sense mostly?
dfratio$first[which(dfratio$ratio<.01)]


# Make sense to be <1
dfratio %>% filter(first == "KYLA") %>% select(first, count_j, count_vf, ratio)
dfratio %>% filter(first == "JADE") %>% select(first, count_j, count_vf, ratio)
dfratio %>% filter(first == "CHRISTIAN") %>% select(first, count_j, count_vf, ratio)
dfratio %>% filter(first == "JOSIAH") %>% select(first, count_j, count_vf, ratio)
dfratio %>% filter(first == "JAMES") %>% select(first, ratio)
dfratio %>% filter(first == "JOHAN") %>% select(first, ratio)
dfratio %>% filter(first == "PABLO") %>% select(first, count_j, count_vf, ratio)
dfratio %>% filter(first == "ENRIQUE") %>% select(first, ratio)
dfratio %>% filter(first == "MARILYN") %>% select(first, ratio)
dfratio %>% filter(first == "SUN") %>% select(first, ratio)
dfratio %>% filter(first == "ANTONIO") %>% select(first, ratio)
dfratio %>% filter(first == "ANTONIA") %>% select(first, ratio)
dfratio %>% filter(first == "NELSON") %>% select(first, ratio)
dfratio %>% filter(first == "JESUS") %>% select(first, count_j, count_vf, ratio)
dfratio %>% filter(first == "HANS") %>% select(first, ratio)
dfratio %>% filter(first == "MOHAMMAD") %>% select(first, count_j, count_vf, ratio)
dfratio %>% filter(first == "ROBERT") %>% select(first, count_j, count_vf, ratio)
dfratio %>% filter(first == "PABLO") %>% select(first, count_j, count_vf, ratio)
dfratio %>% filter(first == "HUNTER") %>% select(first, count_j, count_vf, ratio)
dfratio %>% filter(first == "MARIA") %>% select(first, count_j, count_vf, ratio)
dfratio %>% filter(first == "JOHN") %>% select(first, count_j, count_vf, ratio)
dfratio %>% filter(first == "CHRISTINA") %>% select(first, count_j, count_vf, ratio)
dfratio %>% filter(first == "KATELYN") %>% select(first, count_j, count_vf, ratio)
dfratio %>% filter(first == "CRYSTAL") %>% select(first, count_j, count_vf, ratio)


# make sense to be > 1
dfratio %>% filter(first == "EITAN") %>% select(first, count_j, count_vf, ratio)
dfratio %>% filter(first == "ELIANA") %>% select(first, count_j, count_vf, ratio)
dfratio %>% filter(first == "JOSEPH") %>% select(first, count_j, count_vf, ratio)
dfratio %>% filter(first == "JACOB") %>% select(first, count_j, count_vf, ratio)
dfratio %>% filter(first == "DAVID") %>% select(first, ratio)
dfratio %>% filter(first == "ASHER") %>% select(first, ratio)
dfratio %>% filter(first == "EZRA") %>% select(first, ratio)
dfratio %>% filter(first == "DANIEL") %>% select(first, ratio)
dfratio %>% filter(first == "ITAI") %>% select(first, ratio)
dfratio %>% filter(first == "MIRIAM") %>% select(first, ratio)
dfratio %>% filter(first == "LEAH") %>% select(first, ratio)
dfratio %>% filter(first == "ESTHER") %>% select(first, ratio)
dfratio %>% filter(first == "AARON") %>% select(first, ratio)
dfratio %>% filter(first == "ABRAHAM") %>% select(first, ratio)
dfratio %>% filter(first == "ISAAC") %>% select(first, ratio)
dfratio %>% filter(first == "HYMEN") %>% select(first, ratio)
dfratio %>% filter(first == "HYMAN") %>% select(first, ratio)
dfratio %>% filter(first == "CHAIM") %>% select(first, ratio)
dfratio %>% filter(first == "ARI") %>% select(first, ratio)
dfratio %>% filter(first == "BENJAMIN") %>% select(first, ratio)
dfratio %>% filter(first == "IRVING") %>% select(first, ratio)
dfratio %>% filter(first == "FIVEL") %>% select(first, ratio)
dfratio %>% filter(first == "AVROMI") %>% select(first, count_j, count_vf, ratio)
dfratio %>% filter(first == "ISRAEL") %>% select(first, count_j, count_vf, ratio)

# formerly weird, now fixed
dfratio %>% filter(first == "SARAH") %>% select(first, ratio)
dfratio %>% filter(first == "SARA") %>% select(first, ratio)
dfratio %>% filter(first == "NOA") %>% select(first, ratio)
dfratio %>% filter(first == "SHIRA") %>% select(first, ratio)
dfratio %>% filter(first == "TALIA") %>% select(first, ratio)
dfratio %>% filter(first == "MOSHE") %>% select(first, ratio)
dfratio %>% filter(first == "RACHEL") %>% select(first, ratio)
dfratio %>% filter(first == "YITZCHAK") %>% select(first, ratio)
dfratio %>% filter(first == "YISROEL") %>% select(first, ratio)


# weird that >1?
dfratio %>% filter(first == "KOSUKE") %>% select(first, ratio)
dfratio %>% filter(first == "JAN") %>% select(first, ratio)



# unsure
dfratio %>% filter(first == "LAURA") %>% select(first, ratio)
dfratio %>% filter(first == "ZACHARY") %>% select(first, ratio)
dfratio %>% filter(first == "REBECCA") %>% select(first, ratio)
dfratio %>% filter(first == "OMER") %>% select(first, ratio)
dfratio %>% filter(first == "URI") %>% select(first, ratio)
dfratio %>% filter(first == "HARRY") %>% select(first, ratio)
dfratio %>% filter(first == "GEORGE") %>% select(first, ratio)
dfratio %>% filter(first == "GRACE") %>% select(first, ratio)
dfratio %>% filter(first == "ALEX") %>% select(first, count_j, count_vf, ratio)



#TODO: what kinds of counts do the extreme ratios have?
#TODO: what if subset to only those with counts > c in and else treat as NA?

table(dfratio$count_j)[1:30]
table(dfratio$count_vf)[1:30]

(12525)/sum(table(dfratio$count_j))
(1974650)/sum(table(dfratio$count_vf))


# Old code for exploring cut-offs
#------------------------------------------
# dfratiosub = dfratio %>% filter(!is.na(ratio)) %>% filter(count_j > 5 | count_vf > 10)
# max(dfratiosub$ratio)
# dfratiosub$first[which(dfratiosub$ratio>500)]
# dfratiosub$first[which(dfratiosub$ratio>250 & dfratiosub$ratio<=500)]
# dfratiosub$first[which(dfratiosub$ratio>150 & dfratiosub$ratio<=250)]
# dfratiosub$first[which(dfratiosub$ratio>125 & dfratiosub$ratio<=150)]
# dfratiosub$first[which(dfratiosub$ratio>2 & dfratiosub$ratio<=4)]
# dfratiosub$first[which(dfratiosub$ratio<1)]



################################################################################
#TODO:
# * process and combine names with a space?
# * subset just to those with count > c in at least one dataset
# * re-calculate ratios just for those using that pool of counts

dfratiosub = dfratio %>% filter(!is.na(ratio)) %>% filter(count_j > 10 | count_vf > 10)





# Analysis of size
bothnames = dfratio %>% filter(!is.na(count_j) & !is.na(count_vf))

n1_orig = bothnames %>% summarize(sum(count_j))
n2_orig = bothnames %>% summarize(sum(count_vf))


cvals = c(1, 5, 1:10*10, 1:5*100)
holder1 = rep(NA, length(cvals))
holder2 = rep(NA, length(cvals))
holder3 = rep(NA, length(cvals))
holder4 = rep(NA, length(cvals))

for (i in 1:length(cvals)){
  c = cvals[i]
  sub = bothnames %>% filter(count_j > c | count_vf > c)
  n1 = sub %>% summarize(sum(count_j))
  n2 = sub %>% summarize(sum(count_vf))
  holder1[i] = n2/n1
  holder2[i] = nrow(sub)
  holder3[i] = n1/n1_orig
  holder4[i] = n2/n2_orig
}
holder1
holder2

plot(cvals, holder1,"b",xlab = "cutoff")
par(mfrow=c(1,3))
plot(cvals, holder2,"b", xlab = "cutoff", ylab = "number of names")
plot(cvals, holder3,"b", xlab = "cutoff", ylab = "proportion of names count mass")
plot(cvals, holder4,"b", xlab = "cutoff", ylab = "proportion of names count mass")



#Restriction
par(mfrow = c(1,2))
c = 5
for (c in c(5,10)){
  subdfratio = dfratio %>% filter(!is.na(count_j) & !is.na(count_vf) & count_j > c & count_vf > c)
  subdfratio$ratio2 = (subdfratio$prop_j/sum(subdfratio$prop_j))/(subdfratio$prop_vf/sum(subdfratio$prop_vf))
  plot(subdfratio$ratio, subdfratio$ratio2)
  abline(a =0 , b=1, lty = 2)
  abline(h=100)
}
