#'
#' Processing for public release
#'

library(readxl)
df <- read_excel("../../data/SurveyPrivate/survey_for_analysis_100125_with_total_prob.xlsx")

# Subset to those that responded - have no info for these columns for the rest anyway
stopifnot(sum(!is.na(df$recordeddate))== 1765)
df = df[!is.na(df$recordeddate),]

# Cut-off at 1, which also did in analysis
df$total_probs_bounded = pmin(df$total_probs_bounded, 1)

# Rename to interpretable name
df$BISG = df$total_probs_bounded
df$total_probs_bounded = NULL

# Store originals temporarily
originals = df$BISG

# Round to two significant figures
df$BISG = signif(df$BISG, digits = 2)
df$sampling_prob = signif(df$sampling_prob, digits = 2)

# Checks
stopifnot(all(df$BISG > 0, na.rm = TRUE))
stopifnot(all(df$BISG <= 1, na.rm = TRUE))
stopifnot(sum(is.na(df$BISG)) == 6) #the 6 people who could not be linked back to sampling frame

stopifnot(all(df$sampling_prob > 0, na.rm = TRUE))
stopifnot(all(df$sampling_prob < 1, na.rm = TRUE))
stopifnot(sum(is.na(df$sampling_prob)) == 6) #the 6 people who could not be linked back to sampling frame

write.csv(df, "../../data/SurveyPublic/survey_for_analysis_100125_for_release.csv")
