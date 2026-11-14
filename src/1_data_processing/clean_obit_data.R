#
# Author: Kyla Chasalow
#
# Take initial obit file and do some correcting of errors via
# via manual lookup matching against known information
#

library(tidyverse)
library(readr)

####################
# Load Raw Data
####################
df = read_csv("../../data/ObitData/obit_list_091824.csv")

##########################
# Data Cleaning and Checks
##########################
#----------------------------------------------------------------------
# Deal with NA rows (two)
j = which(is.na(df$state)); print(j)
df[j[1],] #everything missing row
df[j[2],] #data error row, last is "NEW JERSEY JEWISH MEMORIAL CHAPEL BELOW"
df = df[-j,] #remove all-missing row at start (?)
#----------------------------------------------------------------------
# Check that there are no 0 or 1 letter last names
stopifnot(all(nchar(df$last) >=2))
#----------------------------------------------------------------------
# Examine two-letter last names
# Confirmed Issue: names with last name "ah" have actually just taken
# the 'a"h" short for "alav hashalom" as last name
# 25 instances
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/dr-howard-lerner-ah-11377173
i = which(df$last == "AH" & df$first == "HOWARD" & df$dod == "07/19/2023"); stopifnot(i > 0)
df[i,"last"] = "LERNER"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/yisroel-semmel-ah-10610990
i = which(df$last == "AH" & df$first == "YISROEL" & df$dod == "03/01/2022"); stopifnot(i > 0)
df[i,"last"] = "SEMMEL"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/miriam-kalker-ah-10714040
i = which(df$last == "AH" & df$first == "MIRIAM" & df$dod == "04/17/2022"); stopifnot(i > 0)
df[i,"last"] = "KALKER"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/mollie-freireich-ah-11086755
i = which(df$last == "AH" & df$first == "MOLLIE" & df$dod == "12/30/2022"); stopifnot(i > 0)
df[i,"last"] = "FREIREICH"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/mindy-zimmer-gunsburg-ah-11180955
i = which(df$last == "AH" & df$first == "MINDY" & df$dod == "03/03/2023"); stopifnot(i > 0)
df[i,"last"] = "ZIMMER-GUNSBURG"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/rifka-borenstein-ah-10700639
i = which(df$last == "AH" & df$first == "RIFKA" & df$dod == "04/07/2022"); stopifnot(i > 0)
df[i,"last"] = "BORENSTEIN"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/gila-rollhaus-ah-11135351
i = which(df$last == "AH" & df$first == "GILA" & df$dod == "02/03/2023"); stopifnot(i > 0)
df[i,"last"] = "ROLLHAUS"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/zena-barsky-ah-11261248
i = which(df$last == "AH" & df$first == "ZENA" & df$dod == "04/23/2023"); stopifnot(i > 0)
df[i,"last"] = "BARSKY"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/dr-reuven-rosenkranz-ah-10806035
i = which(df$last == "AH" & df$first == "REUVEN" & df$dod == "06/24/2022"); stopifnot(i > 0)
df[i,"last"] = "ROSENKRANZ"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/julius-leon-yazersky-ah-10477689
i = which(df$last == "AH" & df$first == "JULIUS" & df$dod == "12/06/2021"); stopifnot(i > 0)
df[i,"last"] = "YAZERSKY"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/ruth-weichselbaum-ah-10920319
i = which(df$last == "AH" & df$first == "RUTH" & df$dod == "09/08/2022"); stopifnot(i > 0)
df[i,"last"] = "WEICHSELBAUM"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/amelia-gazal-ah-11317317
i = which(df$last == "AH" & df$first == "AMELIA" & df$dod == "06/05/2023"); stopifnot(i > 0)
df[i,"last"] = "GAZAL"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/ann-bennett-ah-10666020
i = which(df$last == "AH" & df$first == "ANN" & df$dod == "03/30/2022"); stopifnot(i > 0)
df[i,"last"] = "BENNETT"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/dr-mandell-ganchrow-ah-10588124
i = which(df$last == "AH" & df$first == "MANDELL" & df$dod == "2/15/2022"); stopifnot(i > 0)
df[i,"last"] = "GANCHROW"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/aharon-gitchel-ah-11226544
i = which(df$last == "AH" & df$first == "AHARON" & df$dod == "03/29/2023"); stopifnot(i > 0)
df[i,"last"] = "GITCHEL"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/bella-grossman-ah-10584857
i = which(df$last == "AH" & df$first == "BELLA" & df$dod == "02/12/2022"); stopifnot(i > 0)
df[i,"last"] = "GROSSMAN"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/bryna-greenberg-ah-10548053
i = which(df$last == "AH" & df$first == "BRYNA" & df$dod == "01/25/2022"); stopifnot(i > 0)
df[i,"last"] = "GREENBERG"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/david-moskowitz-ah-11490064
i = which(df$last == "AH" & df$first == "DAVID" & df$dod == "10/09/2023"); stopifnot(i > 0)
df[i,"last"] = "MOSKOWITZ"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/dolores-stiefel-dolnane-ah-11263195
i = which(df$last == "AH" & df$first == "DOLORES" & df$dod == "4/24/2023"); stopifnot(i > 0)
df[i,"last"] = "STIEFEL-DOLNANE"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/joseph-adlerstein-ah-11080962
i = which(df$last == "AH" & df$first == "JOSEPH" & df$dod == "12/28/2022"); stopifnot(i > 0)
df[i,"last"] = "ALDERSTEIN"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/marilyn-frid-ah-10582666
i = which(df$last == "AH" & df$first == "MARILYN" & df$dod == "02/10/2022"); stopifnot(i > 0)
df[i,"last"] = "FRID"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/marshall-joffe-ah-10387939
i = which(df$last == "AH" & df$first == "MARSHALL" & df$dod == "10/05/2021"); stopifnot(i > 0)
df[i,"last"] = "JOFFE"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/rebecca-goldberg-ah-10614007
i = which(df$last == "AH" & df$first == "REBECCA" & df$dod == "03/03/2022"); stopifnot(i > 0)
df[i,"last"] = "GOLDBERG"
#https://www.dignitymemorial.com/obituaries/forest-hills-ny/sarah-barbanel-ah-11310745
i = which(df$last == "AH" & df$first == "SARAH" & df$dod == "05/30/2023"); stopifnot(i > 0)
df[i,"last"] = "BARBANEL"
#df[which(df$last == "AH"),]
#df = df[-which(df$last == "AH"),]


# Other 2-letter last names seem more realistic. Some paired with Asian first names
# e.g. Helen Ax is real https://www.legacy.com/us/obituaries/nytimes/name/hellen-ax-obituary?id=28464193
print(df[nchar(df$last) ==2,], n = 40)

# The really long names also look like real last names
print(df[nchar(df$last) > 15,], n = 20)
#----------------------------------------------------------------------
# Remove the one row from Canada
df = filter(df, state != "CN")

# Save final file
write_csv(df, "../../DATA/ObitData/obit_list_092424.csv")
