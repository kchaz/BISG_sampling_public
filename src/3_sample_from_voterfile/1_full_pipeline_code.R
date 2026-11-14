#'
#'
#' Author: Laura Royden, with a few parts written by Kyla Chasalow
#'
#'
#' Notes:
#' - script uses paths on cluster
#' - script assumes that code in 2_sampler/sampler_run_code has already been run
#'

library(haven)
library(readr)
library(tidyverse)
library(data.table)
library(glue)


# start of process that was formerly in vf_names_with_jewish_probs.R------------
#
#-------------------------------------------------------------------------------
# Generate tables from the voter files into state_name_counts tables
# containing: first | last | state | count
#-------------------------------------------------------------------------------

setwd("/files/Vm2Uniform/tabs")

# pick which files you want to do
all.files <- list.files(pattern="tab")

# loop over each file
for (file in all.files) {

  # read in the columns we need and store as data
  data <- fread(file, select = c("Voters_FirstName", "Voters_LastName"))

  # get counts
  data <- data %>%
    group_by(Voters_FirstName, Voters_LastName) %>%
    tally() %>%
    ungroup()

  # add the two digit state code from the file as a column
  data$State <- str_extract(file, "(?<=VM2Uniform--).{2}")

  # do some renameing
  data <- data %>%
    rename(First = Voters_FirstName,
           Last = Voters_LastName,
           Count = n) %>%
    select(First, Last, State, Count) %>%
    mutate(Last = str_replace_all(Last, "[^[:graph:]]", " "))

  # save these files
  write_csv(data,
            paste0("/files/fullname_state_counts/fullnames_",
                   str_extract(file, "(?<=VM2Uniform--).{2}"),
                   ".csv"))
  print(file)
}

rm(list = ls())


#-------------------------------------------------------------------------------
# Get P(S=s|G) in each state G and then use these and the P(R=1|G) to calculate
# bounds on what P(S=s|G,R=1) can be. There is a conservative and a non-conservative
# version that we carry through but we ultimately use the non-conservative version
# because it guarantees final P(R=1|S,G) probs are in [0,1]
#-------------------------------------------------------------------------------

# Load AJPP files
brandeis_directory = fs::path("/files/Brandeis_data/")
p_R0_gvn_state = readRDS(file.path(brandeis_directory, "p_R0_gvn_state.RDS"))
p_R1_gvn_state = readRDS(file.path(brandeis_directory, "p_R1_gvn_state.RDS"))

p_R0_gvn_state_lower = readRDS(file.path(brandeis_directory, "p_R0_gvn_state_lower.RDS"))
p_R0_gvn_state_upper = readRDS(file.path(brandeis_directory, "p_R0_gvn_state_upper.RDS"))

p_R1_gvn_state_lower = readRDS(file.path(brandeis_directory, "p_R1_gvn_state_lower.RDS"))
p_R1_gvn_state_upper = readRDS(file.path(brandeis_directory, "p_R1_gvn_state_upper.RDS"))

# Load obituary data and create holder with names in obit data to fill state by state
obit_directory = fs::path("/files/original_obituary_data/obit_list_092424-fulldehyph-check.csv")

df_obit = read.csv(obit_directory)
lnames = sort(unique(df_obit$last))
holder = data.frame("surname" = lnames)

# Prep to load state surname files
state_surname_directory = fs::path("/files/fullname_state_counts/")
bound_directory = fs::path("/files/fullname_state_counts/bounds/")
fs::dir_create(bound_directory)

state_csv_files <- list.files(state_surname_directory, pattern = "^fullnames_[A-Za-z]{2}\\.csv$", full.names = TRUE)
nf = length(state_csv_files)

# process each file
for(file in state_csv_files){

  df = read.csv(file)

  state <- stringr::str_extract(file, "[A-Za-z]{2}(?=\\.csv$)"); print(state)
  lowname = glue("lower_{state}")
  upname = glue("upper_{state}")
  lowname_c = glue("lower_{state}_c")
  upname_c = glue("upper_{state}_c")

  # capitalize for consistency with other files
  df$Last = toupper(df$Last)

  # Get proportions for each last name
  df <- df %>%
    group_by(Last) %>%
    summarise(Total = sum(Count), .groups = "drop") %>%
    mutate(prop = Total / sum(Total))

  # Non-CONSERVATIVE VERSION of bounds that will definitely fix > 1 problem
  df[,lowname] = pmax((df$prop - p_R0_gvn_state[state])/p_R1_gvn_state[state],0)
  df[,upname]  = pmin(df$prop / p_R1_gvn_state[state],1)
  stopifnot(all(df[,upname] >= df[,lowname]))
  stopifnot(all(df[,upname] >= df$prop))

  # CONSERVATIVE VERSION of bounds that may not fix > 1 problem
  df[,lowname_c] = pmax((df$prop - p_R0_gvn_state_upper[state])/p_R1_gvn_state_upper[state],0)
  df[,upname_c]  = pmin(df$prop / p_R1_gvn_state_lower[state],1)
  stopifnot(all(df[,upname_c] >= df[,lowname_c]))
  stopifnot(all(df[,upname_c] >= df$prop))

  # Reduce to only Lasts in obituary data
  dfsub <- df %>%
    filter(Last %in% lnames) %>%
    select(Last, .data[[lowname]], .data[[upname]], .data[[lowname_c]], .data[[upname_c]])

  # Update Holder
  # If a Last does not occur in that state in the voter file, get a NA here. This is ok.
  dfsub <- dfsub %>% rename(surname = Last)
  holder = left_join(holder, dfsub, by = "surname")
}
write.csv(holder, file.path(bound_directory, "surname_prob_bounds.csv"),row.names = F)


#-------------------------------------------------------------------------------
# Apply bounds on P(S=s|G,R=1)
# to obituary data estimation output from 2_sampler/sampler_run_code
#-------------------------------------------------------------------------------
# Get directory where p(S|G,R=1) estimates are stored
obit_outputs_directory = fs::path("/files/state_obitdata_surname_probs_new")

# Load files to process
state_csv_files <- list.files(obit_outputs_directory, pattern = "^[A-Za-z]{2}_.*\\.csv$", full.names = TRUE)
nf = length(state_csv_files)

# Load bounds
df_bounds = read.csv(file.path(bound_directory, "surname_prob_bounds.csv"))

# Process each file
for(file in state_csv_files){

  df = read.csv(file)

  state <- stringr::str_extract(fs::path_file(file), "^[A-Za-z]{2}"); print(state)
  lowname = glue("lower_{state}")
  upname = glue("upper_{state}")
  lowname_c = glue("lower_{state}_c")
  upname_c = glue("upper_{state}_c")

  new <- left_join(
    select(df, surname, p_surname_gvn_state_R1),
    select(df_bounds, surname, all_of(c(upname, lowname, upname_c, lowname_c))),
    by = "surname"
  )

  # enforce  non-conservative bounds
  new$p_bounded = ifelse(new$p_surname_gvn_state_R1 > new[[upname]],
                         new[[upname]], #if true, replaces by upper. If false, goes to lower check
                         ifelse(new$p_surname_gvn_state_R1 < new[[lowname]],
                                new[[lowname]],
                                new$p_surname_gvn_state_R1))
  new$p_bounded[is.na(new$p_bounded)] = new$p_surname_gvn_state_R1[is.na(new$p_bounded)]


  # checks and stats
  stopifnot(all(new$p_bounded <= new[,upname], na.rm = T))
  stopifnot(all(new$p_bounded >= new[,lowname], na.rm = T))
  print(glue("Non-conservative upper bounds enforced: {sum(new$p_bounded == new[,upname], na.rm = T)} of {nrow(df)}"))
  print(glue("Non-conservative lower bounds enforced: {sum(new$p_bounded == new[,lowname], na.rm = T)} of {nrow(df)}"))
  df$numerator_bounded = new$p_bounded * df$p_R1_gvn_state

  # enforce conservative bounds
  new$p_bounded_cnsv = ifelse(new$p_surname_gvn_state_R1 > new[[upname_c]],
                              new[[upname_c]], #if true, replaces by upper. If false, goes to lower check
                              ifelse(new$p_surname_gvn_state_R1 < new[[lowname_c]],
                                     new[[lowname_c]],
                                     new$p_surname_gvn_state_R1))
  new$p_bounded_cnsv[is.na(new$p_bounded_cnsv)] = new$p_surname_gvn_state_R1[is.na(new$p_bounded_cnsv)]
  stopifnot(all(new$p_bounded_cnsv <= new[,upname_c], na.rm = T))
  stopifnot(all(new$p_bounded_cnsv >= new[,lowname_c], na.rm = T))
  print(glue("Conservative upper bounds enforced: {sum(new$p_bounded_cnsv == new[,upname_c], na.rm = T)} of {nrow(df)}"))
  print(glue("Conservative lower bounds enforced: {sum(new$p_bounded_cnsv == new[,lowname_c], na.rm = T)} of {nrow(df)}"))
  df$numerator_bounded_cnsv = new$p_bounded_cnsv * df$p_R1_gvn_state

  # re-save file
  write.csv(df, file = fs::path(file), row.names = F)
  print(glue("Resaved file for {state}"))
}


# continue process that was formerly in vf_names_with_jewish_probs.R -----------

# -------------------------------------------------------------------------------
# get just first name counts from vf tables
# ------------------------------------------------------------------------------
setwd("/n/holylabs/LABS/imai_lab/Lab/lroyden")
all.files <- list.files("files/fullname_state_counts/", pattern = "\\.csv$")

firstnames <- data.frame()

for (file in all.files) {

  # read in each file with counts
  state <- read_csv(paste0("files/fullname_state_counts/", file))

  # grab just the first name and counts
  state <- state %>%
    select(First, Count) %>%
    group_by(First) %>%
    summarize(Count = sum(Count))

  firstnames <- rbind(firstnames, state)

  print(file)
}

firstnames <- firstnames %>%
  group_by(First) %>%
  summarize(Count = sum(Count))

write_csv(firstnames, "files/vf_firstnames.csv")



#-------------------------------------------------------------------------------
# calculate jewish probabilities for surnames in vf
#-------------------------------------------------------------------------------
# this joins together p(surname|state, jewish) * p(jewish|state) which is
# in state_obitdata_surname_probs.csv
# with p(surname|state) which comes from fullname_state_counts/
# the output is a list of surnames from each state's voterfile and the probability
# of Jewish given surname and geography (but not yet first name)

setwd("/")

# load in files
denom.files <- list.files("files/fullname_state_counts", pattern=".csv")
num.files <- list.files("files/state_obitdata_surname_probs_new", pattern=".csv")


counts <- data.frame(state=character(),
                     total_vf=numeric(),
                     in_obits=numeric(),
                     truncated=numeric(),
                     truncated_cnsv=numeric())

# for each state....
for (i in 1:length(denom.files)) {
  # read in the denominator data
  denom <- read.csv(paste0("files/fullname_state_counts/",denom.files[i]))

  # get in format of last | count | prop
  denom <- denom %>% group_by(Last) %>% summarize(Count = sum(Count)) %>%
    mutate(prop = Count/sum(Count)) %>%
    rename(LastName = Last)

  # clean and rename to match for joining
  denom$LastName <- str_replace_all(denom$LastName,"[^[:graph:]]", " ")
  denom$surname <- toupper(denom$LastName)

  # drop unnecessary column
  denom <- denom %>% select(-LastName)

  ## read in the numerator probability data
  num <- read.csv(paste0("files/state_obitdata_surname_probs_new/", num.files[i]), encoding = "UTF-8")

  ## join the data -- right join so only names in the state's vf remain
  join <- right_join(num, denom, by="surname") %>%
    rename(denom = prop, num = numerator, num_bounded = numerator_bounded,
           num_bounded_cnsv = numerator_bounded_cnsv) %>% select(-X)

  # put zero in for any missing in the numerator
  join$num <- ifelse(is.na(join$num), 0, join$num)
  join$num_bounded <- ifelse(is.na(join$num_bounded), 0, join$num_bounded)
  join$num_bounded_cnsv <- ifelse(is.na(join$num_bounded_cnsv), 0, join$num_bounded_cnsv)


  # put in state for any missing state
  join$state <- ifelse(is.na(join$state), na.omit(join$state), join$state)

  # calculate the bayesian probability
  join$total <- join$num/join$denom
  join$total_bounded <- join$num_bounded/join$denom
  join$total_bounded_cnsv <- join$num_bounded_cnsv/join$denom

  truncated <- sum(join$total > 1)
  truncated_cnsv <- sum(join$total_bounded_cnsv > 1)

  join$total <- if_else(join$total > 1, 1, join$total)
  join$total_bounded_cnsv <- if_else(join$total_bounded_cnsv > 1, 1, join$total_bounded_cnsv)

  # pick columns we want in the order we want
  join <- join %>% select(surname, state, Count, num, denom, total, total_bounded, total_bounded_cnsv)

  write_csv(join, paste0("files/jewish_probs/state_vf_with_surname_only_jprobs/vf_surname_jprobs_",
                         str_extract(denom.files[i], "(?<=fullnames_)[^\\.]+"), ".csv"))

  ## get counts of names used and total in voter file
  state_count <- data.frame(
    state = str_extract(denom.files[i], "(?<=fullnames_)[^\\.]+"),
    total_vf = nrow(denom),
    in_obits = sum(join$num > 0),
    truncated = truncated,
    truncated_cnsv = truncated_cnsv)

  counts <- rbind(counts, state_count)

  print(denom.files[i])
}

counts$prop <- counts$in_obits/counts$total_vf

write_csv(counts, "files/jewish_probs/surname_used_by_state.csv")

#-------------------------------------------------------------------------------
# calculate jewish ratios for first names
#-------------------------------------------------------------------------------
# this code generates the ratio of p(first|Jewish) to p(first) overall for each name
# only for names that occur at least once in each file and with
# and at least 10 instances of each name in one file

# takes in files/vf_firstnames.csv (count of first names in voter file),
# names from obits (obit_list_092424_first_names_only.csv)
# and names scraped from obits (clean_first_names-v2.csv)

setwd("/")

# read in Jewish names from scraped data and obit
jewish_names <- read_csv("files/clean_first_names-v2.csv")
obit_names <- read_csv("files/obit_list_092424_first_names_only.csv")

# combine into one thing of Jewish names
all_jewish_names <- c(jewish_names$first, obit_names$first)

rm(jewish_names, obit_names)

# make table of Jewish first name counts
j_firstnames <- as.data.frame(all_jewish_names) %>%
  group_by(all_jewish_names) %>%
  count() %>%
  ungroup() %>%
  rename(first = all_jewish_names,
         count_j = n)

rm(all_jewish_names)

vf_firstnames <- read_csv("files/vf_firstnames.csv")

vf_firstnames$First <- str_replace_all(vf_firstnames$First,"[^[:graph:]]", " ")

# make table of voter file first name counts
vf_firstnames <- vf_firstnames %>%
  mutate(first = toupper(First)) %>%
  select(first, Count) %>%
  rename(count_vf = Count)

jnames_leftover <- j_firstnames$first[is.na(match(j_firstnames$first, vf_firstnames$first))]

combined_fn <- full_join(j_firstnames, vf_firstnames, by="first")

rm(j_firstnames, vf_firstnames)

# fix apostrophe spacing problems, and fix extra letter
combined_fn <- combined_fn %>%
  mutate(first = str_replace_all(first, "\\s*'\\s*", "'")) %>%
  mutate(first = str_replace_all(first, "^[a-zA-Z]\\s|\\s[a-zA-Z]$", ""))

combined_fn[is.na(combined_fn)] <- 0

# turn nicknames into real names

# Create a named vector for nicknames
nickname_lookup <- c(
  "ABE" = "ABRAHAM",
  "ABBY" = "ABIGAIL",
  "AL" = "ALBERT",
  "ALLIE" = "ALLISON",
  "ALLY" = "ALLISON",
  "ANDI" = "ANDREA",
  "ANDIE" = "ANDREA",
  "ANDY" = "ANDREW",
  "ART" = "ARTHUR",
  "ARTIE" = "ARTHUR",
  "BARB" = "BARBARA",
  "BECCA" = "REBECCA",
  "BECKY" = "REBECCA",
  "BEN" = "BENJAMIN",
  "BENJI" = "BENJAMIN",
  "BENNY" = "BENJAMIN",
  "BILL" = "WILLIAM",
  "BILLY" = "WILLIAM",
  "BOB" = "ROBERT",
  "BOBBY" = "ROBERT",
  "BURT" = "BURTON",
  "CHARLEY" = "CHARLES",
  "CHARLIE" = "CHARLES",
  "DAN" = "DANIEL",
  "DANNY" = "DANIEL",
  "DAVE" = "DAVID",
  "DEBBIE" = "DEBORAH",
  "DEBBY" = "DEBORAH",
  "DEBI" = "DEBORAH",
  "DICK" = "RICHARD",
  "DON" = "DONALD",
  "DONNIE" = "DONALD",
  "DONNY" = "DONALD",
  "DOUG" = "DOUGLAS",
  "DOUGIE" = "DOUGLAS",
  "EDDIE" = "EDWARD",
  "FRANCIE" = "FRANCINE",
  "GEOFF" = "GEOFFREY",
  "HANK" = "HENRY",
  "HOWIE" = "HOWARD",
  "IRV" = "IRVING",
  "JAKE" = "JACOB",
  "JACKI" = "JACQUELINE",
  "JACKIE" = "JACQUELINE",
  "JACKY" = "JACQUELINE",
  "JEFF" = "JEFFREY",
  "JEN" = "JENNIFER",
  "JENNIE" = "JENNIFER",
  "JENNY" = "JENNIFER",
  "JESS" = "JESSICA",
  "JESSI" = "JESSICA",
  "JESSIE" = "JESSICA",
  "JESSY" = "JESSICA",
  "JIM" = "JAMES",
  "JIMMY" = "JAMES",
  "JOE" = "JOSEPH",
  "JOEY" = "JOSEPH",
  "JOHNNY" = "JOHN",
  "JON" = "JONATHAN",
  "JONNY" = "JONATHAN",
  "JOSH" = "JOSHUA",
  "KEN" = "KENNETH",
  "KENNY" = "KENNETH",
  "LAINIE" = "ELAINE",
  "LEN" = "LEONARD",
  "LENNY" = "LEONARD",
  "LIZ" = "ELIZABETH",
  "LIZZI" = "ELIZABETH",
  "LIZZIE" = "ELIZABETH",
  "LIZZY" = "ELIZABETH",
  "LOU" = "LOUIS",
  "LOUIE" = "LOUIS",
  "MADDIE" = "MADELINE",
  "MADDY" = "MADELINE",
  "MAURY" = "MAURICE",
  "MARTY" = "MARTIN",
  "MATT" = "MATTHEW",
  "MEG" = "MARGARET",
  "MIKE" = "MICHAEL",
  "MIKEY" = "MICHAEL",
  "MITCH" = "MITCHELL",
  "MORT" = "MORTON",
  "NICK" = "NICHOLAS",
  "NICKY" = "NICHOLAS",
  "PEG" = "MARGARET",
  "PEGGY" = "MARGARET",
  "PHIL" = "PHILIP",
  "RICH" = "RICHARD",
  "RICHIE" = "RICHARD",
  "RICK" = "RICHARD",
  "ROB" = "ROBERT",
  "ROBBIE" = "ROBERT",
  "ROBBY" = "ROBERT",
  "ROD" = "RODNEY",
  "RODDY" = "RODNEY",
  "RON" = "RONALD",
  "RUSS" = "RUSSELL",
  "RUTHIE" = "RUTH",
  "SANDI" = "SANDRA",
  "SANDEE" = "SANDRA",
  "SANDY" = "SANDRA",
  "SID" = "SIDNEY",
  "SOL" = "SOLOMON",
  "STEVE" = "STEVEN",
  "STAN" = "STANLEY",
  "SUE" = "SUSAN",
  "SUSIE" = "SUZANNE",
  "SUZIE" = "SUZANNE",
  "SUZY" = "SUZANNE",
  "TED" = "THEODORE",
  "TIM" = "TIMOTHY",
  "TIMMY" = "TIMOTHY",
  "TOM" = "THOMAS",
  "TOMMIE" = "THOMAS",
  "TOMMY" = "THOMAS",
  "TRISH" = "PATRICIA",
  "VAL" = "VALERIE",
  "VICKIE" = "VICTORIA",
  "VICKY" = "VICTORIA",
  "VIKKY" = "VICTORIA",
  "WES" = "WESLEY",
  "WILL" = "WILLIAM",
  "YOSSI" = "YOSEF",
  "ZAC" = "ZACHARY",
  "ZACH" = "ZACHARY",
  "ZACK" = "ZACHARY",
  "ZAK" = "ZACHARY"
)


combined_fn <- combined_fn %>%
  mutate(first = recode(first, !!!nickname_lookup))


# regroup by name and summarize counts
combined_fn <- combined_fn %>% group_by(first) %>%
  summarise(count_j = sum(count_j), count_vf = sum(count_vf))


# Select only rows for which (1) both names occur in the files and
#(2) at least one file has a count of >= 10 (does not have to be both)

combined_fn <- combined_fn %>% filter(count_j > 0, count_vf > 0) %>%
  filter(count_j > 9 | count_vf > 9)


# Calculate, using only those >= 10 rows
#(obit count / obit total) / (voter file count / voter file total)
combined_fn <- combined_fn %>%
  mutate(prop_j = count_j / sum(count_j),
         prop_vf = count_vf / sum(count_vf),
         ratio = prop_j/prop_vf)

# save first name ratio results
write_csv(combined_fn %>% filter(ratio > 10), "files/jewish_probs/firstnames/first-names-that-were-truncated.csv")

combined_fn <- combined_fn %>% mutate(ratio = if_else(ratio > 10, 10, ratio))

write_csv(combined_fn, "files/jewish_probs/firstnames/firstnames_jratio.csv")
write_csv(as.data.frame(jnames_leftover), "files/jewish_probs/firstnames/unused-jewish-first-names.csv")


#-------------------------------------------------------------------------------
# join voter file names with first name Jewish ratios and last name jewish probs
# to get final probs
#-------------------------------------------------------------------------------
# this code combines the first | last | state | count table
# with the ratio of the first name being Jewish and P(Jewish | Last, Geo)
# takes in firstnames_jratio.csv, fullname_state_counts/ (first | last | state | count), and
# files/jewish_probs/state_vf_with_surname_only_jprobs/surname_j_probs_

rm(list = ls())

setwd("/")


#load in first name ratios
firstnames <- read_csv("files/jewish_probs/firstnames/firstnames_jratio.csv")

#pick the only columns we need
firstnames <- firstnames %>%
  select(first, ratio) %>%
  rename(First = first)


# choose files to do and extract the state names from them
all.files <- str_extract(list.files("files/fullname_state_counts/", pattern=".csv"), "[A-Z]{2}(?=\\.csv)")

# create a lookup table for nicknames
nicknames <- c(
  "ABE" = "ABRAHAM", "ABBY" = "ABIGAIL", "AL" = "ALBERT",
  "ALLIE" = "ALLISON", "ALLY" = "ALLISON", "ANDI" = "ANDREA",
  "ANDIE" = "ANDREA", "ANDY" = "ANDREW", "ART" = "ARTHUR",
  "ARTIE" = "ARTHUR", "BARB" = "BARBARA", "BECCA" = "REBECCA",
  "BECKY" = "REBECCA", "BEN" = "BENJAMIN", "BENJI" = "BENJAMIN",
  "BENNY" = "BENJAMIN", "BILL" = "WILLIAM", "BILLY" = "WILLIAM",
  "BOB" = "ROBERT", "BOBBY" = "ROBERT", "BURT" = "BURTON",
  "CHARLEY" = "CHARLES", "CHARLIE" = "CHARLES", "DAN" = "DANIEL",
  "DANNY" = "DANIEL", "DAVE" = "DAVID", "DEBBIE" = "DEBORAH",
  "DEBBY" = "DEBORAH", "DEBI" = "DEBORAH", "DICK" = "RICHARD",
  "DON" = "DONALD", "DONNIE" = "DONALD", "DONNY" = "DONALD",
  "DOUG" = "DOUGLAS", "DOUGIE" = "DOUGLAS", "EDDIE" = "EDWARD",
  "FRANCIE" = "FRANCINE", "GEOFF" = "GEOFFREY", "HANK" = "HENRY",
  "HOWIE" = "HOWARD", "IRV" = "IRVING", "JAKE" = "JACOB",
  "JACKI" = "JACQUELINE", "JACKIE" = "JACQUELINE", "JACKY" = "JACQUELINE",
  "JEFF" = "JEFFREY", "JEN" = "JENNIFER", "JENNIE" = "JENNIFER",
  "JENNY" = "JENNIFER", "JESS" = "JESSICA", "JESSI" = "JESSICA",
  "JESSIE" = "JESSICA", "JESSY" = "JESSICA", "JIM" = "JAMES",
  "JIMMY" = "JAMES", "JOE" = "JOSEPH", "JOEY" = "JOSEPH",
  "JOHNNY" = "JOHN", "JON" = "JONATHAN", "JONNY" = "JONATHAN",
  "JOSH" = "JOSHUA", "KEN" = "KENNETH", "KENNY" = "KENNETH",
  "LAINIE" = "ELAINE", "LEN" = "LEONARD", "LENNY" = "LEONARD",
  "LIZ" = "ELIZABETH", "LIZZI" = "ELIZABETH", "LIZZIE" = "ELIZABETH",
  "LIZZY" = "ELIZABETH", "LOU" = "LOUIS", "LOUIE" = "LOUIS",
  "MADDIE" = "MADELINE", "MADDY" = "MADELINE", "MAURY" = "MAURICE",
  "MARTY" = "MARTIN", "MATT" = "MATTHEW", "MEG" = "MARGARET",
  "MIKE" = "MICHAEL", "MIKEY" = "MICHAEL", "MITCH" = "MITCHELL",
  "MORT" = "MORTON", "NICK" = "NICHOLAS", "NICKY" = "NICHOLAS",
  "PEG" = "MARGARET", "PEGGY" = "MARGARET", "PHIL" = "PHILIP",
  "RICH" = "RICHARD", "RICHIE" = "RICHARD", "RICK" = "RICHARD",
  "ROB" = "ROBERT", "ROBBIE" = "ROBERT", "ROBBY" = "ROBERT",
  "ROD" = "RODNEY", "RODDY" = "RODNEY", "RON" = "RONALD",
  "RUSS" = "RUSSELL", "RUTHIE" = "RUTH", "SANDI" = "SANDRA",
  "SANDEE" = "SANDRA", "SANDY" = "SANDRA", "SID" = "SIDNEY",
  "SOL" = "SOLOMON", "STEVE" = "STEVEN", "STAN" = "STANLEY",
  "SUE" = "SUSAN", "SUSIE" = "SUZANNE", "SUZIE" = "SUZANNE",
  "SUZY" = "SUZANNE", "TED" = "THEODORE", "TIM" = "TIMOTHY",
  "TIMMY" = "TIMOTHY", "TOM" = "THOMAS", "TOMMIE" = "THOMAS",
  "TOMMY" = "THOMAS", "TRISH" = "PATRICIA", "VAL" = "VALERIE",
  "VICKIE" = "VICTORIA", "VICKY" = "VICTORIA", "VIKKY" = "VICTORIA",
  "WES" = "WESLEY", "WILL" = "WILLIAM", "YOSSI" = "YOSEF",
  "ZAC" = "ZACHARY", "ZACH" = "ZACHARY", "ZACK" = "ZACHARY",
  "ZAK" = "ZACHARY"
)


### loop will start here

for (file in all.files) {


  vf_counts <- read_csv(paste0("files/fullname_state_counts/fullnames_", file, ".csv"))

  ## get all data ready for joins

  ### clean voter file aggregates
  vf_counts$Last <- str_replace_all(vf_counts$Last,"[^[:graph:]]", " ")
  vf_counts$Last <- toupper(vf_counts$Last)

  vf_counts$First <- str_replace_all(vf_counts$First,"[^[:graph:]]", " ")
  vf_counts$First <- toupper(vf_counts$First)

  # fix apostrophe spacing problems, and fix extra letter
  vf_counts <- vf_counts %>%
    mutate(First = str_replace_all(First, "\\s*'\\s*", "'")) %>%
    mutate(First = str_replace_all(First, "^[a-zA-Z]\\s|\\s[a-zA-Z]$", ""))

  # replace nicknames
  vf_counts <- vf_counts %>%
    mutate(First = recode(First, !!!nicknames))

  #recombine after cleaning names
  vf_counts <- vf_counts %>% group_by(First, Last, State) %>% summarize(Count = sum(Count))


  ### read in surname probs for the state
  surname_j_probs <- read_csv(paste0("files/jewish_probs/state_vf_with_surname_only_jprobs/vf_surname_jprobs_", file, ".csv"),
                              col_select = c(surname, total, total_bounded, total_bounded_cnsv))

  ### change column names to match
  surname_j_probs <- surname_j_probs %>%
    rename(Last = surname)

  ## join

  # join <- left_join(vf_counts, surname_j_probs, by="Last") %>%
  #   rename(last_probs = total) %>%
  #   left_join(firstnames, by="First") %>%
  #   rename(first_ratio = ratio) %>%
  #   mutate(total_probs = if_else(!is.na(first_ratio), last_probs*first_ratio, last_probs)) %>%
  #   select(First, Last, Count, State, first_ratio, last_probs, total_probs)

  # Convert to data.table
  vf_counts <- as.data.table(vf_counts)
  surname_j_probs <- as.data.table(surname_j_probs)
  firstnames <- as.data.table(firstnames)

  # Perform left joins
  join <- merge(vf_counts, surname_j_probs, by = "Last", all.x = TRUE)
  setnames(join, "total", "last_probs")  # Rename column
  setnames(join, "total_bounded", "last_probs_bounded")
  setnames(join, "total_bounded_cnsv", "last_probs_bounded_cnsv")


  join <- merge(join, firstnames, by = "First", all.x = TRUE)
  setnames(join, "ratio", "first_ratio")  # Rename column

  # Calculate total_probs efficiently
  join[, total_probs := fifelse(!is.na(first_ratio), last_probs * first_ratio, last_probs)]
  join[, total_probs_bounded := fifelse(!is.na(first_ratio), last_probs_bounded * first_ratio, last_probs_bounded)]
  join[, total_probs_bounded_cnsv := fifelse(!is.na(first_ratio), last_probs_bounded_cnsv * first_ratio, last_probs_bounded_cnsv)]


  # Select the required columns
  join <- join[, .(First, Last, Count, State, first_ratio, last_probs, total_probs, total_probs_bounded, total_probs_bounded_cnsv)]



  write_csv(join,
            paste0("/files/vf_fullnames_with_jprobs/", file, ".csv"))

  rm(join, vf_counts, surname_j_probs)
  print(file)
}


#-------------------------------------------------------------------------------
# Content originally in get_pi_total_across_states.R
#-------------------------------------------------------------------------------

#' Goal: calculate total sampling probability normalizing total pi(g) for each state
#' and then calculate sampling probabilities
#'
#' Total Calculations
#' ------------------
#' Decisions:
#'  1. use total_bounded to do this. These are not so different from
#'     total_bounded_cnsv (see output of run_checks.R) but have guarantees that at
#'     least the surname component in final files is <1
#'
#'  2. because of multiplication by first name ratios, some total_bounded_cnsv
#'     still end up >1 because we truncated first name ratios at 10 and if surname
#'     based probability was >1/10, multiplying by 10 could result in a number >1.
#'     From some exploration, this is not an extremely common problem however.
#'     Decision:  Truncate these to 1 before sum(prob * count) for the state.
#'     Reason: else final probabilities sometimes end up >1; also worry that
#'     the distribution of probabilities will end up too skewed.
#'
#' Optimal Sampling Calculations
#' ------------------------------
#' This script also calculates, for each state:
#'     * sum(count) = total number of people
#'     * sum(count[prob > 0]) = total number of people with positive probability in sampling frame
#'
#' It uses these and the known P(R=1|G) to do an optimal sample allocation to
#' calculate the T_g
#'
#' Normalizing Ratio Calculations
#' ------------------------------
#' After calculating these T_g, the script also calculates the ratio T_g/Total Prob needed for
#' adjusting all the state probabilities so that they sum to T_g and hence so that the
#' expected sample size from that state is Poisson sampling is T_g.
#'
#'
#' See get_final_sampling_probabilities.R for the use of these ratios to actually calculate
#' the final sampling probabilities



# Calculating Normalizing Denom and State Totals ------------------------

finalprob_dir <- fs::path("/files/vf_fullnames_with_jprobs")

# Get all CSV files with exactly two-letter names
state_csv_files <- list.files(finalprob_dir,
                              pattern = "^[A-Za-z]{2}\\.csv$",
                              full.names = TRUE)

nf = length(state_csv_files)

prob_totals = rep(NA, nf)
state_counts = rep(NA, nf)
state_frame_counts = rep(NA, nf)
states = rep(NA, nf)


for (i in 1:nf){

  file = state_csv_files[i]
  states[i] = stringr::str_extract(file, "[A-Za-z]{2}(?=\\.csv$)")

  print(states[i])
  df = read.csv(file)

  # Get rows with a non-NA obit data based probability and
  # make sure none missing total or a count
  sub = filter(df, !is.na(df$total_probs_bounded))
  stopifnot(sum(is.na(sub$total_probs_bounded)) == 0)
  stopifnot(sum(is.na(sub$Count)) == 0)
  if(nrow(sub)<nrow(df)){print(glue("{states[i]}: {nrow(df)-nrow(sub)} people with NA prob dropped"))}

  # Get Total of Probabilities
  prob_totals[i]  = sum(pmin(sub$total_probs_bounded,1) * sub$Count)

  # Get total number of people
  state_counts[i] = sum(df$Count) #includes NA people

  # Get total number of people with non-zero prob of being sampled
  state_frame_counts[i] = sum(sub$Count[sub$total_probs_bounded > 0])

}

df_tot = data.frame(state = states,
                    prob_total = prob_totals,
                    total_count = state_counts,
                    frame_count = state_frame_counts)


# Target Sample Size Allocation -------------------------------------------

# Set overall target
Target = 50000

# State-Specific Targets

p_R1_gvn_state = readRDS(fs::path("/files/Brandeis_data/p_R1_gvn_state.RDS"))


# August 11, 2025 Comment: the commented out line below is a bug - this is not how
# I intended to calculate the T_g targets - not sure how it happened!
# I am fixing this but leaving the old one uncommented as documentation of what was
# run originally when doing the sampling. I tried it out and the T_g calculated by
# Each are thankfully very correlated. This bug only affected the T_g's calculated
# And not any of the downstream calculations using T_g so thankfully no major impact and
# no impact on validity of resulting sample.
#numerators = sqrt(p_R1_gvn_state[df_tot$state] * df_tot$prob_total)

numerators = sqrt(p_R1_gvn_state[df_tot$state]) * df_tot$frame_count
sample_fracs = numerators/sum(numerators)
state_targets = round(Target*sample_fracs)
print(glue("Target Total:{sum(state_targets)}"))

# Add targets to data frame
df_tot$target = state_targets[df_tot$state]

# Add ratios to dataframe
df_tot$ratio = df_tot$target/df_tot$prob_total

write.csv(df_tot, file.path(finalprob_dir, "normalizing_totals_by_state.csv"),
          row.names = F)

#-------------------------------------------------------------------------------
# Content originally in get_final_sampling_probabilities.R
#-------------------------------------------------------------------------------
#'
#'  For calculating final probabilities using totals
#'  from get_pi_total_across_states.R for normalization. This is necessary
#'  to control the final sample size of the poisson sampling
#'
#'  V1 - calculate total across probs and normalize by target / total
#'  without doing anything to the probs that end up >1 because of high first name ratios
#'
#'  V2 - first take any probs > 1 and set them to 1. Then calculate total across probs
#'  and normalize by target / total. #THIS PROBABLY BETTER - else likely to get > 1
#'  final sampling probs.
#'
#'  Note: could in theory use overall target T/(overall total) to normalize
#'  but this results in > 1 probabilities. There is unfortunately no guarantee of this
#'  renormalizing resulting in >1 probabilities. Here, for small states, this >1 does
#'  happen. A solution is to specify state-specific targets. This does not guarantee a fix to the issue
#'  but it tends to help... long as states with lower probabilities / lower total number of
#'  possible people to sample also have lower (enough) targets.
#'

finalprob_dir <- fs::path("/files/vf_fullnames_with_jprobs")
#finalprob_dir = "test/final"

# Get all CSV files with exactly two-letter names
state_csv_files <- list.files(finalprob_dir,
                              pattern = "^[A-Za-z]{2}\\.csv$",
                              full.names = TRUE)

nf = length(state_csv_files)

# Load ratios to multiply by
df_total = read.csv(file.path(finalprob_dir, "normalizing_totals_by_state.csv"))

info = list()


for (i in 1:nf){

  file = state_csv_files[i]
  state <- stringr::str_extract(file, "[A-Za-z]{2}(?=\\.csv$)")
  info[[state]] = list()
  print(state)

  df = read.csv(file)
  df$sampling_prob = pmin(df$total_probs_bounded,1) * df_total[df_total$state == state, "ratio"]

  # Remove NA rows
  num_NA_removed = sum(is.na(df$sampling_prob))
  df = df[!is.na(df$sampling_prob),]

  info[[state]][["prob_total"]] = sum(df$sampling_prob * df$Count)
  info[[state]][["target"]] = df_total$target[df_total$state == state]

  # re-save df file, overwriting old
  write.csv(df, file, row.names = F)

  # Summary Info below here

  info[[state]][["max_probs"]] = max(df$sampling_prob, na.rm = TRUE)
  info[[state]][["min_probs"]] = min(df$sampling_prob, na.rm = TRUE)
  info[[state]][["n"]] = nrow(df)
  info[[state]][["n_NA_removed"]] = num_NA_removed

  thresholds <- c(0, 0.1, 0.01, 0.001)
  for (t in thresholds) {
    name <- glue(">{t}_probs")
    mask <- df$sampling_prob > t
    info[[state]][[name]] <- list(sum(mask, na.rm = TRUE), mean(mask, na.rm = TRUE))
  }
  info[[state]][["=max_probs"]] = list(sum(df$sampling_prob == max(df$sampling_prob, na.rm = TRUE), na.rm = TRUE),
                                       mean(df$sampling_prob == max(df$sampling_prob, na.rm = TRUE), na.rm = TRUE))
  # Glimpse of top and bottom names
  c = 6
  top <- df %>%
    arrange(desc(sampling_prob)) %>%
    slice_head(n = c) %>%
    select(First, Last, Count, sampling_prob, sampling_prob)
  info[[state]][[glue("top{c}")]] = top

  top_lq_max <- df %>%
    filter(sampling_prob < (max(sampling_prob, na.rm = T) - 1e-4))  %>%
    arrange(desc(sampling_prob)) %>%
    slice_head(n = c) %>%
    select(First, Last, Count, sampling_prob, sampling_prob)
  info[[state]][[glue("top_lq_max{c}")]] = top_lq_max

  bottom_gt_zero <- df %>%
    filter(sampling_prob > 0) %>%
    arrange(sampling_prob) %>%
    slice_head(n = c) %>%
    select(First, Last, Count, sampling_prob, sampling_prob)
  info[[state]][[glue("bottom_gt_zero{c}")]] = bottom_gt_zero

  saveRDS(info, file.path(finalprob_dir, "sampling_prob_info.RDS"))
}

#info = readRDS("sampling_prob_info.RDS")
#info


#-------------------------------------------------------------------------------
# Code formerly in sampling.R -------------------------------------------------
### this script takes voter files from /Vm2Uniform/tabs, cleans them,
### and combines them with sampling probabilities from /vf_fullnames_with_jprobs/
### then it samples according to those probabilities and saves the names and LALs
### then it pulls out the columns needed from the voter files from sampled LALs


rm(list = ls())

# sample LALs -------------------------------------------------------------



# read in LALVOTERID, Voters_FirstName, Voters_LastName from the voter file

setwd("/files/Vm2Uniform/tabs")


# pick which files you want to do
all.files <- list.files(pattern="tab")


file <- all.files

set.seed(02138)

# for all files....
for (file in all.files) {

  # store which state you're doing
  state <- str_extract(file, "(?<=VM2Uniform--).{2}")

  # read in voter ID and name
  data <- fread(file,
                select = c("Voters_FirstName", "Voters_LastName", "LALVOTERID"))

  # clean names according to prior procedure

  ## create lookup for nicknames
  nicknames <- c(
    "ABE" = "ABRAHAM", "ABBY" = "ABIGAIL", "AL" = "ALBERT",
    "ALLIE" = "ALLISON", "ALLY" = "ALLISON", "ANDI" = "ANDREA",
    "ANDIE" = "ANDREA", "ANDY" = "ANDREW", "ART" = "ARTHUR",
    "ARTIE" = "ARTHUR", "BARB" = "BARBARA", "BECCA" = "REBECCA",
    "BECKY" = "REBECCA", "BEN" = "BENJAMIN", "BENJI" = "BENJAMIN",
    "BENNY" = "BENJAMIN", "BILL" = "WILLIAM", "BILLY" = "WILLIAM",
    "BOB" = "ROBERT", "BOBBY" = "ROBERT", "BURT" = "BURTON",
    "CHARLEY" = "CHARLES", "CHARLIE" = "CHARLES", "DAN" = "DANIEL",
    "DANNY" = "DANIEL", "DAVE" = "DAVID", "DEBBIE" = "DEBORAH",
    "DEBBY" = "DEBORAH", "DEBI" = "DEBORAH", "DICK" = "RICHARD",
    "DON" = "DONALD", "DONNIE" = "DONALD", "DONNY" = "DONALD",
    "DOUG" = "DOUGLAS", "DOUGIE" = "DOUGLAS", "EDDIE" = "EDWARD",
    "FRANCIE" = "FRANCINE", "GEOFF" = "GEOFFREY", "HANK" = "HENRY",
    "HOWIE" = "HOWARD", "IRV" = "IRVING", "JAKE" = "JACOB",
    "JACKI" = "JACQUELINE", "JACKIE" = "JACQUELINE", "JACKY" = "JACQUELINE",
    "JEFF" = "JEFFREY", "JEN" = "JENNIFER", "JENNIE" = "JENNIFER",
    "JENNY" = "JENNIFER", "JESS" = "JESSICA", "JESSI" = "JESSICA",
    "JESSIE" = "JESSICA", "JESSY" = "JESSICA", "JIM" = "JAMES",
    "JIMMY" = "JAMES", "JOE" = "JOSEPH", "JOEY" = "JOSEPH",
    "JOHNNY" = "JOHN", "JON" = "JONATHAN", "JONNY" = "JONATHAN",
    "JOSH" = "JOSHUA", "KEN" = "KENNETH", "KENNY" = "KENNETH",
    "LAINIE" = "ELAINE", "LEN" = "LEONARD", "LENNY" = "LEONARD",
    "LIZ" = "ELIZABETH", "LIZZI" = "ELIZABETH", "LIZZIE" = "ELIZABETH",
    "LIZZY" = "ELIZABETH", "LOU" = "LOUIS", "LOUIE" = "LOUIS",
    "MADDIE" = "MADELINE", "MADDY" = "MADELINE", "MAURY" = "MAURICE",
    "MARTY" = "MARTIN", "MATT" = "MATTHEW", "MEG" = "MARGARET",
    "MIKE" = "MICHAEL", "MIKEY" = "MICHAEL", "MITCH" = "MITCHELL",
    "MORT" = "MORTON", "NICK" = "NICHOLAS", "NICKY" = "NICHOLAS",
    "PEG" = "MARGARET", "PEGGY" = "MARGARET", "PHIL" = "PHILIP",
    "RICH" = "RICHARD", "RICHIE" = "RICHARD", "RICK" = "RICHARD",
    "ROB" = "ROBERT", "ROBBIE" = "ROBERT", "ROBBY" = "ROBERT",
    "ROD" = "RODNEY", "RODDY" = "RODNEY", "RON" = "RONALD",
    "RUSS" = "RUSSELL", "RUTHIE" = "RUTH", "SANDI" = "SANDRA",
    "SANDEE" = "SANDRA", "SANDY" = "SANDRA", "SID" = "SIDNEY",
    "SOL" = "SOLOMON", "STEVE" = "STEVEN", "STAN" = "STANLEY",
    "SUE" = "SUSAN", "SUSIE" = "SUZANNE", "SUZIE" = "SUZANNE",
    "SUZY" = "SUZANNE", "TED" = "THEODORE", "TIM" = "TIMOTHY",
    "TIMMY" = "TIMOTHY", "TOM" = "THOMAS", "TOMMIE" = "THOMAS",
    "TOMMY" = "THOMAS", "TRISH" = "PATRICIA", "VAL" = "VALERIE",
    "VICKIE" = "VICTORIA", "VICKY" = "VICTORIA", "VIKKY" = "VICTORIA",
    "WES" = "WESLEY", "WILL" = "WILLIAM", "YOSSI" = "YOSEF",
    "ZAC" = "ZACHARY", "ZACH" = "ZACHARY", "ZACK" = "ZACHARY",
    "ZAK" = "ZACHARY"
  )

  ## strip weird characters and make uppercase
  data$Voters_LastName <- str_replace_all(data$Voters_LastName,"[^[:graph:]]", " ")
  data$Voters_LastName <- toupper(data$Voters_LastName)

  data$Voters_FirstName <- str_replace_all(data$Voters_FirstName,"[^[:graph:]]", " ")
  data$Voters_FirstName <- toupper(data$Voters_FirstName)

  ## fix apostrophe spacing problems, and fix extra letter
  data <- data %>%
    mutate(Voters_FirstName = str_replace_all(Voters_FirstName, "\\s*'\\s*", "'")) %>%
    mutate(Voters_FirstName = str_replace_all(Voters_FirstName, "^[a-zA-Z]\\s|\\s[a-zA-Z]$", ""))



  ## replace nicknames
  data <- data %>%
    mutate(Voters_FirstName = recode(Voters_FirstName, !!!nicknames))

  # read in vf_probs
  vf_probs <- read_csv(paste0("/files/vf_fullnames_with_jprobs/", state, ".csv"))

  # rename so columns match, drop 0 sampling prob
  vf_probs <- vf_probs %>% select(First, Last, sampling_prob) %>%
    rename(Voters_FirstName = First, Voters_LastName = Last) %>%
    filter(sampling_prob > 0)


  # merge into vf

  data <- inner_join(data, vf_probs, by = c("Voters_FirstName", "Voters_LastName"))


  # sample

  data$sampled <- rbinom(nrow(data), 1, data$sampling_prob)

  sample <- data %>% filter(sampled == 1)

  rm(data, vf_probs)


  write_csv(sample,
            paste0("/files/samples/",
                   str_extract(file, "(?<=VM2Uniform--).{2}"),
                   ".csv"))

  rm(sample)
}



# combine sampled LALs into one file --------------------------------------

# combine into one sample file
setwd("/files/samples")

all.files <- list.files(pattern=".csv")

full.sample <- data.frame()

for (file in all.files){

  sample <- read_csv(file)

  sample$state <- str_extract(file, "[A-Z]{2}(?=\\.csv$)")

  full.sample <- rbind(full.sample, sample)

}



# match VF columns from LAL sample --------------------------------------------


setwd("/files/Vm2Uniform/tabs")


# pick which files you want to do
all.files <- list.files(pattern="tab")

# for all files....

for (file in all.files) {

  # store which state you're doing
  state.sample <- str_extract(file, "(?<=VM2Uniform--).{2}")

  LALs <- full.sample %>% filter(state == state.sample) %>%
    select(LALVOTERID)

  # read in voter ID and name
  data <- fread(file,
                select = c("Voters_FirstName", "Voters_LastName", "LALVOTERID",
                           "Voters_Gender", "Voters_BirthDate",
                           "Ethnic_Description",
                           "EthnicGroups_EthnicGroup1Desc",
                           "Parties_Description",
                           "VoterParties_Change_Changed_Party",
                           "Residence_Addresses_AddressLine",
                           "Residence_Addresses_ExtraAddressLine",
                           "Residence_Addresses_City",
                           "Residence_Addresses_State",
                           "Residence_Addresses_Zip",
                           "Mailing_Addresses_AddressLine" ,
                           "Mailing_Addresses_ExtraAddressLine" ,
                           "Mailing_Addresses_City" ,
                           "Mailing_Addresses_State",
                           "Mailing_Addresses_Zip",
                           "General_2024",
                           "Primary_2024",
                           "PresidentialPrimary_2024",
                           "General_2022",
                           "Primary_2022",
                           "General_2020",
                           "Primary_2020",
                           "PresidentialPrimary_2020"
                ))

  sampled.people <- data %>% filter(LALVOTERID %in% LALs$LALVOTERID)

  write.csv(sampled.people,
            paste0("/files/sampled_people/",
                   state.sample, ".csv"))

  rm(sampled.people, LALs, state.stample, data)
}



# combine sampled.people into one file

setwd("/files/sampled_people")


all.files <- list.files(pattern=".csv")
vf.sample <- data.frame()

for (file in all.files){
  sample <- read_csv(file)
  sample$state <- str_extract(file, "[A-Z]{2}(?=\\.csv$)")
  vf.sample <- rbind(vf.sample, sample)
  rm(sample)
}

rm(file, all.files)

vf.sample <- vf.sample %>% select(-`...1`)

write.csv(vf.sample, "files/sample_apr23.csv")
