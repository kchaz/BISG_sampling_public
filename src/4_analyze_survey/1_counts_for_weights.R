
library(haven)
library(readr)
library(tidyverse)
library(data.table)
library(glue)


setwd("/files/Vm2Uniform/tabs")


# pick which files you want to do
all.files <- list.files(pattern="tab")


# for loop for each file

for (file in all.files) {
  
  # read in the columns we need and store as data  
  data <- fread(file, 
                select = c("Voters_FirstName", "Voters_LastName", "Voters_Gender", "Voters_Age", "Parties_Description"))
  
  # add the two digit state code from the file as a column
  data$State <- str_extract(file, "(?<=VM2Uniform--).{2}")
  
  data <- data %>%
    rename(First = Voters_FirstName, 
           Last = Voters_LastName, 
           Gender = Voters_Gender,
           Age = Voters_Age,
           Party = Parties_Description) %>%
    select(First, Last, State, Gender, Age, Party) %>%
    mutate(Last = str_replace_all(Last, "[^[:graph:]]", " "))
  
  
  ### clean voter file aggregates
  data$Last <- str_replace_all(data$Last,"[^[:graph:]]", " ") 
  data$Last <- toupper(data$Last)
  
  data$First <- str_replace_all(data$First,"[^[:graph:]]", " ") 
  data$First <- toupper(data$First)
  
  # fix apostrophe spacing problems, and fix extra letter
  data <- data %>% 
    mutate(First = str_replace_all(First, "\\s*'\\s*", "'")) %>%
    mutate(First = str_replace_all(First, "^[a-zA-Z]\\s|\\s[a-zA-Z]$", ""))
  
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
  
  
  
  # replace nicknames
  data <- data %>%
    mutate(First = recode(First, !!!nicknames))
  
  state_code <- str_extract(file, "(?<=VM2Uniform--).{2}")
  
  # Construct the file path and read it
  state_prob <- read_csv(
    file.path("/files/vf_fullnames_with_jprobs", paste0(state_code, ".csv"))
  )
  
  
  state_prob <- state_prob %>% select(First, Last, total_probs_bounded)
  
  left_join(data, state_prob, by=c("First", "Last")) -> data
  
  write_csv(data, 
            paste0("files/counts_for_weights/", 
                   str_extract(file, "(?<=VM2Uniform--).{2}"),
                   ".csv"))
  
  print(file)
  
}

rm(list = ls())





