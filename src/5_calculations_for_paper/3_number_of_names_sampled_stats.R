#'
#' Get some stats on the sample file's names and probabilities
#'
#' Author: Kyla Chasalow
#'
#' AI Use: Script written originally by me, checked by GPT-6.1 Sol Medium,
#' which caught no major issues
#'

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Get info on number of names sampled
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
df_sample = read.csv("/n/holylabs/LABS/imai_lab/Lab/lroyden/files/sample_May02-with_total_probs_bounded.csv",)
df_obit = read.csv("/n/holylabs/LABS/imai_lab/Lab/lroyden/files/original_obituary_data/obit_list_092424-fulldehyph-check.csv")

#-------------------------------------------------------------------------------
# Merge info on how common each name was in obit data
#-------------------------------------------------------------------------------
#df_obit = read.csv("../../data/ObitData/obit_list_092424-fulldehyph-check.csv")
obit_counts <- table(df_obit$last)
obit_counts <- data.frame(
  last = names(obit_counts),
  obit_count = as.integer(obit_counts)
)

obit_counts$obit_rank <- dplyr::dense_rank(dplyr::desc(obit_counts$obit_count))
obit_counts <- obit_counts[order(obit_counts$obit_rank),]

df_sample$Voters_LastName <- toupper(df_sample$Voters_LastName)

df_sample <- merge(
  df_sample,
  obit_counts,
  by.x = "Voters_LastName",
  by.y = "last",
  all.x = TRUE,
  sort = FALSE
)

stopifnot(!anyNA(df_sample$obit_rank)) #check everyone matched and got a rank

#-------------------------------------------------------------------------------
# Main summary stats
#-------------------------------------------------------------------------------
n <- nrow(df_sample)

unique_lastnames <- unique(df_sample$Voters_LastName)
unique_firstnames <- unique(df_sample$Voters_FirstName)

last_counts <- sort(table(df_sample$Voters_LastName), decreasing = TRUE)
first_counts <- sort(table(df_sample$Voters_FirstName), decreasing = TRUE)

num_last_names <- length(unique_lastnames)
num_first_names <- length(unique_firstnames)

summary_stats <- data.frame(
  Statistic = c(
    "Number of people sampled",
    "Number of unique last names",
    "Number of unique first names",
    "Most common last name",
    "Largest count of a last name",
    "2nd most common last name",
    "2nd largest count of a last name",
    "Most common first name",
    "Largest count of a first name",
    "2nd most common first name",
    "2nd largest count of a first name",
    "Mean count per last name",
    "Mean count per first name",
    "Median count per last name",
    "Median count per first name",
    "Number of sampled people with names in top 10",
    "Number of obit names in top 10",
    "Number of sampled people with names in top 25",
    "Number of obit names in top 25",
    "Number of sampled people with names in top 100",
    "Number of obit names in top 100",
    "Number of sampled people with names in top 150",
    "Number of obit names in top 150",
    "Number of sampled people with names in top 200",
    "Number of obit names in top 200",
    "Number of sampled people with names in top 250",
    "Number of obit names in top 250",
    "Number of sampled people with names below top 250",
    "Number of obit names below top 250",
    "Mean_obituary_name_rank",
    "Median_obituary_name_rank",
    "Max_obituary_name_rank",
    "Min_obituary_name_rank"
  ),
  Value = c(
    n,
    num_last_names,
    num_first_names,
    names(last_counts)[1],
    last_counts[1],
    names(last_counts)[2],
    last_counts[2],
    names(first_counts)[1],
    first_counts[1],
    names(first_counts)[2],
    first_counts[2],
    mean(last_counts),
    mean(first_counts),
    median(last_counts),
    median(first_counts),
    sum(df_sample$obit_rank <= 10),
    sum(obit_counts$obit_rank <= 10), #to capture ties
    sum(df_sample$obit_rank <= 25),
    sum(obit_counts$obit_rank <= 25),
    sum(df_sample$obit_rank <= 100),
    sum(obit_counts$obit_rank <= 100),
    sum(df_sample$obit_rank <= 150),
    sum(obit_counts$obit_rank <= 150),
    sum(df_sample$obit_rank <= 200),
    sum(obit_counts$obit_rank <= 200),
    sum(df_sample$obit_rank <= 250),
    sum(obit_counts$obit_rank <= 250),
    sum(df_sample$obit_rank >  250),
    sum(obit_counts$obit_rank > 250),
    mean(df_sample$obit_rank),
    median(df_sample$obit_rank),
    max(df_sample$obit_rank),
    min(df_sample$obit_rank)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_stats,
  fs::path("/n/holylabs/LABS/imai_lab/Lab/lroyden/files/paper_stats/sampled_data_name_summary_statistics.csv"),
  row.names = FALSE
)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Summarize the BISG and sampling probabilities of the sampled set
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# First, cap at 1 (do this elsehwere too)
df_sample$total_probs_bounded <- pmin(df_sample$total_probs_bounded, 1)

summarize_values <- function(x) {
  data.frame(
    n = length(x),
    n_missing = sum(is.na(x)),
    min = min(x, na.rm = TRUE),
    Q1 = unname(quantile(x, 0.25, na.rm = TRUE)),
    median = median(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    Q3 = unname(quantile(x, 0.75, na.rm = TRUE)),
    max = max(x, na.rm = TRUE)
  )
}

summary_stats <- rbind(
  data.frame(
    Variable = "total_probs_bounded",
    summarize_values(df_sample$total_probs_bounded)
  ),
  data.frame(
    Variable = "sampling_prob",
    summarize_values(df_sample$sampling_prob)
  )
)

write.csv(
  summary_stats,
  fs::path("/n/holylabs/LABS/imai_lab/Lab/lroyden/files/paper_stats/sampled_data_prob_summary_statistics.csv"),
  row.names = FALSE
)


