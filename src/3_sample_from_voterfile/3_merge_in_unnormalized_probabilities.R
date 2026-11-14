library(tidyverse)
library(glue)

#'
#' Post-hoc script to merge unnormalized probabilities into final sample
#' for purpose of post-survey analysis
#'

path_to_sample = fs::path("/n/holylabs/LABS/imai_lab/Lab/lroyden/files/sample_May02.csv")
finalprob_dir <- fs::path("/n/imai_lab/Lab/lroyden/files/vf_fullnames_with_jprobs")


df_sample = read.csv(path_to_sample)
nr = nrow(df_sample)
nc = ncol(df_sample)

# Create copies of columns that are synced up to name and format in VF files
df_sample <- df_sample %>%
  mutate(
    First = Voters_FirstName,
    Last  = Voters_LastName,
    State = state
  )

# Repeat the processing steps (particularly nick names) needed to
# get matches on the vf_fullnames_with_jprobs files.  Copied from long full code.R
#-------------------------

## create look-up for nicknames
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
df_sample$Last <- str_replace_all(df_sample$Last,"[^[:graph:]]", " ")
df_sample$Last<- toupper(df_sample$Last)

df_sample$First <- str_replace_all(df_sample$First,"[^[:graph:]]", " ")
df_sample$First <- toupper(df_sample$First)

## fix apostrophe spacing problems, and fix extra letter
df_sample <- df_sample %>%
  mutate(First = str_replace_all(First, "\\s*'\\s*", "'")) %>%
  mutate(First = str_replace_all(First, "^[a-zA-Z]\\s|\\s[a-zA-Z]$", ""))

## replace nicknames
df_sample <- df_sample %>% mutate(First = recode(First, !!!nicknames))


#-------------------------

# Get all CSV files with exactly two-letter names
state_csv_files <- list.files(finalprob_dir,
                              pattern = "^[A-Za-z]{2}\\.csv$",
                              full.names = TRUE)

nf = length(state_csv_files)

# Column where probabilities will go
df_sample$Prob <- NA_real_

# Load each state file in turn
for (i in 1:nf){

  file = state_csv_files[i]
  state_name = stringr::str_extract(file, "[A-Za-z]{2}(?=\\.csv$)")
  print(state_name)

  df = read.csv(file) %>%
    select(First, Last, State, total_probs_bounded) %>%
    rename(Prob_new = total_probs_bounded)

  # Join in probabilities, merging together with probs from previous states
  df_sample <- left_join(df_sample, df, by = c("First", "Last", "State")) %>%
    mutate(Prob = coalesce(Prob, Prob_new)) %>%
    select(-Prob_new) #remove temporary column

}

# Remove temp columns
df_sample <- df_sample %>%
  select(-First, -Last, -State)

# Restore prob col name
df_sample <- df_sample %>%
  rename(total_probs_bounded = Prob)

# Checks
#---------
if (nrow(df_sample) != nr) {
  stop(glue("Row count changed: expected {nr}, got {nrow(df_sample)}"))
}
if (ncol(df_sample) != (nc + 1)) {
  stop(glue("Column count changed: expected {nc + 1}, got {ncol(df_sample)}"))
}


# Info check
num_na <- sum(is.na(df_sample$total_probs_bounded))
print(glue("Number of missing total_probs_bounded values: {num_na}"))


write.csv(df_sample,
          "/n/holylabs/LABS/imai_lab/Lab/lroyden/files/sample_May02-with_total_probs_bounded.csv",
          row.names = FALSE)
