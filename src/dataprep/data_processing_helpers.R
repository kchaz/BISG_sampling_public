
#' ----------------------------------------------------------------
#' Function to get a D count matrix from a dataset df containing a column
#' called 'state' and a column called 'surnames' and each entry a
#' count n_{s,g}
#' ----------------------------------------------------------------
get_count_mat <- function(df, loc_names, surnames, surname_colname = "surname") {
  # set-up frame
  n_loc <- length(loc_names)
  n_sur <- length(surnames)
  D <- matrix(0, nrow = n_sur, ncol = n_loc)
  colnames(D) <- loc_names
  rownames(D) <- surnames
  
  # tabulate - result possibly not full dimension in rows/col
  D2 <- table(df[[surname_colname]], df$state)
  
  
  # find any rows in D that are not in D2 and if any,
  # pad D2 with them with 0 counts and make sure order same as D
  added_sur <- surnames[!surnames %in% rownames(D2)]
  additionD2 <- matrix(0,
                       ncol = length(colnames(D2)),
                       nrow = length(added_sur)
  )
  colnames(additionD2) <- colnames(D2)
  rownames(additionD2) <- added_sur
  D2 <- rbind(D2, additionD2) # merge them
  D2 <- D2[rownames(D), ] # same order as D
  
  
  
  if (length(setdiff(colnames(D2), colnames(D))) != 0){
    stop("df contains rows with state values that are not one of the 50 state acronyms or DC")
  }
  
  
  # fill col of D2 into D so that maintain 0 count col
  # for any col of D that do not appear in D2
  D[, colnames(D2)] <- D2
  
  #make sure D order is same order as location names
  D = D[,loc_names]
  
  return(D)
}

