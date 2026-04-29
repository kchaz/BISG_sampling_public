library(tidyverse)


#' @setting_list - a named list of vectors
#'
#' Generates csv.file containing all combinations of settings in
#' @setting_list with columns named according to names of @setting_list
#'
#' If @path_to_existing_df is NULL:
#'  * creates table of settings
#'  * adds a columnn num2run containing @num2run
#'  * returns list with elements "settings" containing table and
#'    "return_code" set to 0
#'
#' If @path_to_existing_df is not NULL,
#'   * tries to load existing file via read.csv - if cannot, raises an error
#'   * creates table of settings
#'   * checks how many times each setting has already been run (see @run_grouping below)
#'   * adds num2run column where num2run = max(@num2run - times run, 0)
#'   * returns list with elements "settings" containing table and
#'     "return_code" set to 1
#'
#' In both cases, also also adds a column run_ID. If there was no
#' existing file to adjust for, all the run_ID's are 1. If there were,
#' the run_ID is <number of runs already run> + 1. When the simulation
#' is next run on that setting, this run_ID can then be saved as a way of
#' marking the distinct run.
#'
#'
#' @run_grouping
#'   * This is relevant if there is an existing file to process. Else it is ignored
#'
#'   * Default is 1. If this is set to NULL, will set it back to 1 (strong default)
#'
#'   * If this is specified and is an INTEGER,
#'     then the function assumes that across all settings, a single
#'     run generates that number of rows in the existing file of simulation output
#'
#'   * If this is specified and is a STRING,
#'     then the function assumes that that string represents a run identifier
#'     column that it can use to identify, for a given setting, distinct runs
#'
#'
#'
#'
#' WARNING
#' ---------
#' This works well for sequential processing but not for parallel processing
#' where there cannot be order-dependent IDs

setting_updater <- function(setting_list,
                            num2run,
                            existing_df = NULL, #this gets priority if given
                            path_to_existing_df = NULL,
                            run_grouping = 1,
                            return_num_runs = F
){

  if(is.null(run_grouping)){run_grouping = 1}
  setting_names = names(setting_list)

  # Generate table of settings just based on setting_list
  settings = expand.grid(setting_list)
  settings[,sort(colnames(settings))] # consistent order

  if(is.null(existing_df)){
      # Try to load or else let existing_df == NULL
      if(!is.null(path_to_existing_df)){
        existing_df <- tryCatch(
          {
            read.csv(path_to_existing_df)
          },
          error = function(e) {
            warning(paste(#"Failed to load path_to_existing_df as csv file:", path_to_existing_df,
              #"\nError:", e$message,
              "\nINSTEAD: loading full settings table from beginning"))
            NULL
          }
        )
      }
  }

  # If no existing file given or if failed to load, simply generate table and return
  if(is.null(path_to_existing_df) | is.null(existing_df)){
    settings$num2run = num2run
    settings$run_ID = rep(1, nrow(settings))

    if(return_num_runs){
      settings$num_runs = rep(0, nrow(settings))
    }

    return(list("settings" = settings,
                "return_code" = 0))
  }

  # If here, having an existing file. Check that it contains setting names
  if(!all(setting_names %in% colnames(existing_df))){
    stop(glue("Some settings missing from existing file: {setdiff(setting_names, colnames(existing_df))}"))
  }


  # Calculate how often each setting already run has been run
  # Option 1: count occurrences of each setting, divide by # outputs per setting run
  if(is.numeric(run_grouping)){
    num_runs_df <- existing_df %>%
      select(all_of(setting_names)) %>%
      group_by(across(setting_names)) %>%
      summarise(num_runs = n() / run_grouping, .groups = "drop")
    # Option 2: use specified column to group on within setting
    # Does not require that each setting has same number of rows
  }else{
    if(!run_grouping %in% colnames(existing_df)){
      stop(glue("run_grouping = {run_grouping} column not found"))
    }
    num_runs_df <- existing_df %>%
      group_by(across(all_of(c(setting_names, run_grouping)))) %>%
      summarise(count = n(), .groups = "drop") %>%  # up to this point gets # of rows per unique setting + run_ID combo
      group_by(across(all_of(setting_names))) %>%   # below this part then gets # of runs for that setting
      summarize(num_runs = n(),
                all_same_count = n_distinct(count), .groups = "drop")

    # Warning if for same setting, dif runs have dif # of rows
    if (!all(num_runs_df$all_same_count == 1)) {
      warning(glue::glue("There exist {sum(num_runs_df$all_same_count > 1)} settings in existing df with different numbers of rows for the same exact setting"))
    }
    num_runs_df$all_same_count = NULL
  }

  # join to settings data and fill in NAs with 0
  settings = settings %>%  left_join(num_runs_df, by = setting_names)
  settings$num_runs[is.na(settings$num_runs)] = 0
  settings$num2run = pmax(num2run - settings$num_runs, 0)  # Compute remaining runs to do

  # add a default run ID for potential use when storing output of next run
  settings$run_ID = settings[,"num_runs"] + 1

  # Remove internal column
  if(!return_num_runs){
    settings$num_runs = NULL
  }

  return(list("settings" = settings,
              "return_code" = 1))
}

# Dev Examples and for testing
#---------------------------
# # Ex1: one output per run
# setting_list = list("a" = c(1,2,3),
#                     "b" = c(1,2,3))
# existing_df = data.frame("a" = c(1,2,3,1,2,3),
#                            "b" = c(1,1,1,1,1,1),
#                            "c" = c(1,2,3,4,5,6))
# write.csv(existing_df, "temp_existing_df.csv", row.names = F)
# out = setting_updater(setting_list,
#                 path_to_existing_df = NULL,
#                 num2run = 4)
# stopifnot(out$return_code == 0)
# stopifnot(all(out$settings$run_ID == 1))
# out = setting_updater(setting_list,
#                 path_to_existing_df = "temp_existing_df.csv",
#                 num2run = 4,
#                 run_grouping = 1)
# stopifnot(out$return_code == 1)
# stopifnot(sum(out$settings$run_ID == 3) == 3)
# stopifnot(sum(out$settings$run_ID == 1) == 6)
# print(out$settings)
# # Should raise a warning
# out = setting_updater(setting_list,
#                 path_to_existing_df = "blahblahblah",
#                 num2run = 4,
#                 run_grouping = 1)
# # Ex2: multiple outputs per run
# setting_list = list("a" = c(1,2),
#                     "b" = c(1,2))
# existing_df = data.frame("a" = c(1,1,1,1,2,2, 2,2),
#                            "b" = c(1,1,2,2,1,1, 2,2),
#                            "run_ID" = c(1,1,1,1,1,1,1,1))
# write.csv(existing_df, "temp_existing_df.csv", row.names = F)
# out = setting_updater(setting_list,
#                 path_to_existing_df = NULL,
#                 num2run = 2)
# stopifnot(out$return_code == 0)
# stopifnot(all(out$settings$num2run == 2))
# stopifnot(all(out$settings$run_ID == 1))
# out = setting_updater(setting_list,
#                 path_to_existing_df = "temp_existing_df.csv",
#                 num2run = 2,
#                 run_grouping = 2)
# stopifnot(out$return_code == 1)
# stopifnot(all(out$settings$num2run == 1))
# stopifnot(all(out$settings$run_ID == 2))
# out = setting_updater(setting_list,
#                       path_to_existing_df = "temp_existing_df.csv",
#                       num2run = 2,
#                       run_grouping = "run_ID")
# stopifnot(out$return_code == 1)
# stopifnot(all(out$settings$num2run == 1))
# stopifnot(all(out$settings$run_ID == 2))
# # Ex 3: multiple outputs per run, but number differing by setting
# setting_list = list("a" = c(1,2),
#                     "b" = c(1,2))
# existing_df = data.frame("a" = c(1,1,1,1,1,2,2, 2,2),
#                          "b" = c(1,1,1,2,2,1,1, 2,2),
#                          "run_ID" = c(1,1,1,1,1,1,1,1,1))
# write.csv(existing_df, "temp_existing_df.csv", row.names = F)
# out = setting_updater(setting_list,
#                       path_to_existing_df = "temp_existing_df.csv",
#                       num2run = 2,
#                       run_grouping = "run_ID")
# stopifnot(out$return_code == 1)
# stopifnot(all(out$settings$num2run == 1))
# stopifnot(all(out$settings$run_ID == 2))
# file.remove("temp_existing_df.csv")


# Function for filtering out incomplete groups if have a set number of
# repetitions that expect for each setting and something has gone wrong
filter_complete_groups <- function(df, setting_names,
                                   output_per_group = 1,
                                   return_filtered = TRUE) {
  df_with_flag <- df %>%
    group_by(across(all_of(setting_names))) %>%
    mutate(is_complete_group = n() == output_per_group) %>%
    ungroup()

  message("Number of incomplete-group rows: ", sum(!df_with_flag$is_complete_group))
  message("Original dim: ", paste(dim(df), collapse = " x "))

  if (return_filtered) {
    df_filtered <- df[df_with_flag$is_complete_group, ]
    message("Filtered dim: ", paste(dim(df_filtered), collapse = " x "))
    return(df_filtered)
  } else {
    return(df_with_flag)
  }
}

# # Create toy dataframe
# toy_df <- data.frame(
#   setting1 = c("A", "A", "A", "B", "B", "C", "C", "C"),
#   setting2 = c(1, 1, 1, 2, 2, 3, 3, 3),
#   value = 1:8
# )
# setting_names <- c("setting1", "setting2")
# result_df <- filter_complete_groups(toy_df, setting_names, output_per_group = 3)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# System 2
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#' This system is meant to work even if using parallel processing where
#' can't be sure in what order files will be run and hence cannot assume
#' e.g., that if you see 2 replicates of a setting, these will be the
#' run_ID = 1 and run_ID = 2 runs of this setting. This sytem requires
#' that in the file of existing settings, a filename column has been saved
#' that records the exact location of a file containing the output from
#' that setting run. It will generate all the filenames of all the settings
#' to be run based on num2run and then it will see which of those filenames
#' (if any) are present in the already run df

# helper: convert a single value to a clean string (no scientific notation)
val_to_str <- function(v) {
  stopifnot(length(v) == 1)
  v <- v[[1]] #unwrap, just in case

  # Warning just means it is not numeric
  suppressWarnings({
    num <- as.numeric(v)
  })
  if (!is.na(num) && is.finite(num)) {
    return(format(num, scientific = FALSE, trim = TRUE))
  }
  v <- as.character(v)
  return(v)
}

get_filename_from_settings <- function(row, setting_names, path = NULL, prefix = NULL) {

  # ensure consistent order
  setting_names <- sort(setting_names)

  # 2) build "<setting value>_<setting value>_... .csv"
  values <- sapply(setting_names, function(v){val_to_str(row[[v]])})
  fname  <- paste0(paste(values, collapse = "_"), ".csv")

  # optional prefix at start
  if (!is.null(prefix) && nzchar(prefix)) { #if non-empty string
    fname <- paste0(prefix, fname)
  }

  # optional directory path
  if (!is.null(path) && nzchar(path)) {
    return(fs::path(path, fname))
  }
  fname
}


# MAIN FUNCTION
filename_based_setting_updater <- function(setting_list,
                                           num2run,
                                           df_already_run,
                                           prefix #used in filenames
                                           ){

  # Variables that form grid
  setting_names = names(setting_list)

  # Generate table of settings just based on setting_list
  df_settings = expand.grid(setting_list)
  df_settings$num2run = num2run

  # Expand to repeat each num2run times
  df_expanded <- df_settings %>%
    uncount(num2run)

  # Create run ID that is unique within each group of settings
  df_expanded <- df_expanded %>%
    group_by(across(all_of(setting_names))) %>%
    mutate(run_ID = row_number())

  # Generate tempfile that includes run ID
  get_filename <- function(row, path){
    num_sur_str <- format(as.numeric(row[['num_sur']]),  scientific = FALSE, trim = TRUE)
    n_str <- format(as.numeric(row[['n']]), scientific = FALSE, trim = TRUE)
    filepath <- fs::path(glue("{path}/tf_{num_sur_str}_{n_str}_{row[['alpha_init_option']]}_{row[['run_ID']]}.csv"))
    return(filepath)
  }

  # Generate file paths for each setting, including run ID
  df_expanded$filename = apply(df_expanded, 1, get_filename_from_settings,
                               setting_names = c(setting_names,"run_ID"),
                               path = tempfile_path, prefix = prefix)

  # Run only on those not yet run
  if(!is.null(df_already_run)){
    to_run = setdiff(df_expanded$filename, df_already_run$filename)
    df_expanded = df_expanded[df_expanded$filename %in% to_run, ]
  }

  return(df_expanded)
}


# FOR GETTING df of ALREADY RUN
# Consolidate all temp files that might create separately
# at first for for ease of parallel stuff
# merge_tf_csvs <- function(dir, pattern, output_file = NULL) {
#   files <- list.files(
#     dir,
#     pattern = pattern,
#     full.names = TRUE
#   )
#   if (length(files) == 0) {
#     return(NULL)
#   }
#   merged <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
#
#   if (!is.null(output_file)){
#     write.csv(merged, file = output_file, row.names = FALSE)
#   }
#   return(merged)
# }

merge_tf_csvs <- function(dir, pattern, output_file = NULL) {
  files <- list.files(
    dir,
    pattern = pattern,
    full.names = TRUE
  )

  if (length(files) == 0) {
    return(NULL)
  }

  dfs <- lapply(files, read.csv, stringsAsFactors = FALSE)

  merged <- dplyr::bind_rows(dfs)

  if (!is.null(output_file)) {
    write.csv(merged, file = output_file, row.names = FALSE)
  }

  merged
}


