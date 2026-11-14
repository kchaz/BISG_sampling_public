#'
#' Estimating raking weights for age, sex, party using BISG probs
#'
#' Author: Kyla Chasalow
#'
#' AI Use: Wrote the set-up part by hand and then detailed psuedocode for each
#' step below. Then prompted ChatGPT 6.1 Sol Medium to write the code. Then
#' copied it in in parts, checked it line by line and made some tweaks.
#' Some of the original psuedocode is preserved as comments
#'
#' Written: October 23, 2026
#'
options(scipen = 999)

library(dplyr)
# NOTE: I had to run this to get birdie to install correctly
# install.packages(
#   c("RcppParallel", "birdie"),
#   repos = "https://cloud.r-project.org",
#   type = "binary"
# )
library(birdie)


# #--------------------------------------------------------------
# # CREATE A DUMMY DATASET TO ILLUSTRATE THE INTENDED FORMAT
# #--------------------------------------------------------------
# set.seed(222)
# n <- 10000000
#
# # Some fake surnames, each a three-letter string
# surname_pool <- sample(
#   apply(expand.grid(LETTERS, LETTERS, LETTERS), 1, paste0, collapse = ""),
#   size = 25
# )
#
# # 25 hypothetical surname groups with different Jewish probabilities
# p_by_surname <- seq(0.05, 0.95, length.out = 25)
# surname_id <- sample(seq_along(p_by_surname), n, replace = TRUE)
# total_prob_bounded <- p_by_surname[surname_id]
#
# # Helper: randomly replace approximately 10% of values with NA
# add_missing <- function(x, missing_rate = 0.10) {
#   x[runif(length(x)) < missing_rate] <- NA
#   x
# }
# df <- data.frame(
#   total_prob_bounded = total_prob_bounded,
#   X = rbinom(n, 25, .5),
#   state = "MA",
#   Age = add_missing(sample(18:90, n, replace = TRUE)),
#   Party = factor(add_missing(sample(
#     c("Democrat", "Republican", "Independent"),
#     n, replace = TRUE, prob = c(0.45, 0.25, 0.30)
#   ))),
#   Sex = factor(add_missing(sample(
#     c("Female", "Male"),
#     n, replace = TRUE, prob = c(0.52, 0.48)
#   )))
# )
# df = data.table::setDT(df)
# state_name = "MA"
# #--------------------------------------------------------------


#-----------
# Settings
#-----------
vf_file_directory = "/n/holylabs/imai_lab/Lab/lroyden/files/counts_for_weights/"
save_directory = "/n/holylabs/imai_lab/Lab/lroyden/files/counts_for_weights/estimated_weights/"
fs::dir_create(save_directory)

X_columns_for_weighted_estimator = c("Age", "Gender")
X_columns_for_birdie_estimator = c("Party")
X_columns = c(
  X_columns_for_birdie_estimator,
  X_columns_for_weighted_estimator
)
R1_prob_column_name = "total_probs_bounded"


#----------------------------
# Helper
#----------------------------
# Return NA when there are no observed values
safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

#----------------------------

#Load each file in turn of format <state abbreviation>.csv
state_csv_files <- list.files(
  vf_file_directory,
  pattern = "^[A-Za-z]{2}\\.csv$",
  full.names = TRUE)

missingness_stats = list()
estimates_list = list()

# Load each state file in turn
for (i in seq_along(state_csv_files)){

  file = state_csv_files[i]
  state_name = stringr::str_extract(file, "[A-Za-z]{2}(?=\\.csv$)")
  print(state_name)

  # Load only columns needed
  df <- data.table::fread(
    file,
    select = c(X_columns, R1_prob_column_name),
    data.table = TRUE
  )
  stopifnot(R1_prob_column_name %in% colnames(df))

  # Cap probabilities at 1 in df so both estimators use corrected values
  # Note: also did this in the actual sampling prob calculations so this is
  # consistent with there. Some can end up > 1 (but not >10) from multiplication
  # by first name r_f ratio
  df[, (R1_prob_column_name) := pmin(get(R1_prob_column_name), 1)]

  # Get fraction missing for each X variable and the average
  # R1 probability among the missing and among the non-missing for each X variable
  # Put them in a named vector and add to a list so that can consolidate into a df below
  state_stats <- unlist(lapply(X_columns, function(x) {
    missing <- is.na(df[[x]]) | df[[x]] == ""
    stats <- c(
      fraction_missing = mean(missing),
      mean_R1prob_missing = safe_mean(df[[R1_prob_column_name]][missing]),
      mean_R1prob_nonmissing = safe_mean(df[[R1_prob_column_name]][!missing])
    )
    setNames(stats, paste(x, names(stats), sep = "_"))
  }), use.names = TRUE)
  missingness_stats[[state_name]] <- as.list(state_stats)


  # Turn all missing X values to "Missing" so gets encoded as a level
  # Convert NA to an explicit "Missing" category in every X column.
  # Also make sure all of the X columns are encoded as factors.
  df[, (X_columns) := lapply(.SD, function(x) {
    values <- as.character(x)
    values[is.na(values) | values == ""] <- "Missing"
    factor(values)
  }), .SDcols = X_columns]


  # Get prob tables with some checks that are valid. Table will only be used
  # in Birdie explicitly but checks matter for all
  r1_prob <- df[[R1_prob_column_name]]
  if (any(!is.na(r1_prob) &
          (!is.finite(r1_prob) | r1_prob < 0 | r1_prob > 1))) {
    stop("Invalid R1 probabilities in state: ", state_name)
  }
  r_probs <- data.frame(
    pr_R1 = r1_prob,    #Note to self: BIRDIE wants columns starting with pr_
    pr_R0 = 1 - r1_prob
  )

  #----------------------
  # Weighted estimator
  #----------------------
  # Get denom of weighted est as just sum of the R1 probs
  denominator = sum(df[[R1_prob_column_name]], na.rm = TRUE)

  # For each X column in X_columns_for_weighted_estimator, calculate the sum
  # of the R1 probabilities for each level of X and then divide by denominator to get
  # weighted estimator for that level of X variable.
  # Want to add to estimates_list list a row that is
  # (X_variable, Variable_Level, State, Estimate)
  state_estimates <- lapply(X_columns_for_weighted_estimator, function(x) {
    result <- df[
      ,
      .(
        Estimate = if (denominator > 0) {
          sum(get(R1_prob_column_name), na.rm = TRUE) / denominator
        } else {
          NA_real_
        }
      ),
      by = .(Variable_Level = get(x)) #this ensures group by X levels and also creates X var level column
    ]
    # Also save the variable, variable level, and state and estimate type with each estimate
    result[, X_variable := x]
    result[, Variable_Level := as.character(Variable_Level)]
    result[, State := state_name]
    result[, Est_type := "weighted"]
    data.table::setcolorder(
      result,
      c("X_variable", "Variable_Level", "State", "Estimate", "Est_type")
    )
    result
  })

  # Append estimates to the existing estimates table for this state
  estimates_list[[state_name]] <- data.table::rbindlist(state_estimates)

  #----------
  # BIRDIE
  #----------

  # Run Birdie fit simple X ~ 1 model for categorical outcome, no geography
  # since already running this separately in each state file
  # NA values of X should be encoded as a level of X called missing
  state_estimates <- lapply(X_columns_for_birdie_estimator, function(x) {

    # Keep missing X values; exclude only missing probabilities
    keep <- !is.na(r1_prob)


    fit <- birdie::birdie(
      r_probs = r_probs[keep, , drop = FALSE],
      formula = reformulate("1", response = x),
      data = df[keep, ],
      family = birdie::cat_dir(),
      algorithm = "em"
    )

    # Extract Pr(X = level | R = R1)
    estimates <- coef(fit)

    if (any(!is.finite(estimates[, "R1"]))) {
      print(summary(fit))
      stop("Non-finite BIRDiE estimates for ", x, " in ", state_name)
    }

    result <- data.table::data.table(
      Variable_Level = rownames(estimates),
      Estimate = as.numeric(estimates[, "R1"])
    )
    result[, X_variable := x]
    result[, State := state_name]
    result[, Est_type := "birdie"]

    data.table::setcolorder(
      result,
      c("X_variable", "Variable_Level", "State", "Estimate", "Est_type")
    )

    result
  })

  # Append estimates to the existing estimates table for this state
  estimates_list[[state_name]] <- data.table::rbindlist(
    c(list(estimates_list[[state_name]]), state_estimates),
    use.names = TRUE
  )

}#i

#-------------------------------------------
# Consolidate across states and save tables
#-------------------------------------------
missingness_df <- data.table::rbindlist(
  missingness_stats,
  idcol = "State"
)
estimates_df <- data.table::rbindlist(
  estimates_list,
  use.names = TRUE
)

data.table::fwrite(
  missingness_df,
  file.path(save_directory, "missingness_by_state.csv")
)
data.table::fwrite(
  estimates_df,
  file.path(save_directory, "estimated_weights_by_state.csv")
)
