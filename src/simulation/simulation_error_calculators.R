#'
#' Functions for calculating mean squared error
#'


#' Looking for scale invariance in following sense:
#' Let a and b be probability vectors that each sum to 1 and are of dimension k
#' For c>1, let a' = rep(a/c, c) and b' = rep(b/c, c) so that a' and b' each still sum to 1
#' (because a componant that usms to 1/c is repeated c times) but their elements are of
#' a systematically smaller magnitude
#'
#' Then want dist(a,b) = dist(a',b')
#'
#' This is satisfied by KL, TV, and rlogMSE but not by rMSE
#'
#' Downside of KL and rlogMSE is that if an entry in a or b (or just b for KL) is 0,
#' it explodes to infinity
#'
#' TV does not have this issue - choosing this
#'


rMSE <- function(a, b) {
  sqrt(mean((a - b)^2))
}
rlogMSE <- function(a, b) {
  sqrt(mean((log(a) - log(b))^2))
}
TV <- function(a, b) {
  0.5 * sum(abs(a - b))
}
KL <- function(a, b) {
  sum(a * (log(a) - log(b)))
}


# a1 <- c(.15, .05, .2, .3, .2, .1)
# stopifnot(sum(a1) == 1)
# b1 <- a1 + rnorm(length(a1), 0, .01)
# b1 <- b1 / sum(b1)
#
# r <- 10
# a2 <- rep(a1 / r, r) # same as above only one scale down, repeat 10 times
# stopifnot(sum(a2) == 1)
# b2 <- rep(b1 / r, r)
#
# rMSE(a1, b1)
# rMSE(a2, b2) # smaller b/c smaller scale
#
# rlogMSE(a1, b1)
# rlogMSE(a2, b2) # now they are similar, 2 case larger more often?
#
# TV(a1, b1)
# TV(a2, b2) # also similar scale?
#
# KL(a1, b1)
# KL(a2, b2) # same scale






#----------------------------------------
# Plotting helper
#----------------------------------------
plot_error_comparison <- function(postmean_errors,
                                  mle_errors,
                                  parameter_name) {
  l <- min(postmean_errors, mle_errors)
  u <- max(postmean_errors, mle_errors)
  absmax <- max(abs(l), abs(u))
  lim <- c(-absmax, absmax) # make it symmetric

  par(mfrow = c(1, 3))
  breaks <- 50

  hist(postmean_errors,
    main = paste0(glue("Distribution of posterior mean estimate errors for\n all {parameter_name} parameters")),
    breaks = breaks,
    xlim = lim,
    xlab = "Posterior Mean Estimation Error"
  )
  abline(v = 0, lwd = 2)

  hist(mle_errors,
    main = paste0(glue("Distribution of MLE estimate errors for\n all {parameter_name} parameters")),
    breaks = breaks,
    xlim = lim,
    xlab = "MLE Estimation Error"
  )
  abline(v = 0, lwd = 2)

  plot(postmean_errors, mle_errors,
    ylab = "MLE Estimation Error",
    xlab = "Posterior Mean Estimation Error",
    ylim = lim,
    xlim = lim,
    main = "Comparison of Errors"
  )
  abline(v = 0, lwd = 1, col = "black", lty = 2)
  abline(h = 0, lwd = 1, col = "black", lty = 2)
}



#----------------------------------------
# Function to calculate error for alphas
#----------------------------------------
get_alpha_error_analysis <- function(D, postmean, true_prob_marg) {
  ns <- apply(D, 1, sum)
  mle <- matrix(ns / sum(ns), ncol = 1)

  # make sure all a matrix even single vector
  if (is.null(dim(postmean))) {
    postmean <- matrix(postmean, ncol = 1)
  }
  true_prob_marg <- matrix(true_prob_marg, ncol = 1)
  stopifnot(dim(mle)[1] == length(true_prob_marg))
  stopifnot(dim(postmean[1]) == length(true_prob_marg))

  # Subtract off of each column
  # postmean_errors <- sweep(postmean, 1, true_prob_marg, "-")
  # mle_errors <- sweep(mle, 1, true_prob_marg, "-")


  outnames <- colnames(postmean)
  if (is.null(outnames)) {
    outnames <- as.character(seq_along(ncol(postmean)))
  }
  if (length(outnames) == 1) {
    outnames <- ""
  }

  outlist <- list()
  outlist["alpha_mle_rmse"] <- rMSE(mle, true_prob_marg)
  outlist["alpha_mle_TV"] <- TV(mle, true_prob_marg)
  for (i in seq_along(outnames)) {
    outlist[glue("alpha_postmean{outnames[i]}_rmse")] <- rMSE(postmean[, i], true_prob_marg)
    outlist[glue("alpha_postmean{outnames[i]}_TV")] <- TV(postmean[, i], true_prob_marg)
  }
  return(outlist)
}



#----------------------------------------
# Function to calculate MSEs for thetas
#----------------------------------------
get_theta_error_analysis <- function(postmean, D, true_prob_mat, burnin = 1,
                                     mle_smoothing = 0.000001) {
  if (any(dim(postmean) != dim(true_prob_mat))) {
    stop("Error: true_prob_mat dimensions and postmean matrix dimensions do not match")
  }
  if (any(dim(postmean) != dim(D))) {
    stop("Error: D dimensions and postmean matrix dimensions do not match")
  }

  # Calculate MLE with tiny bit of smoothing to avoid 0's
  mle <- apply(D + mle_smoothing, 2, function(v) {
    return(v / sum(v))
  })

  # plotting
  # if (plot){
  #   postmean_errors = postmean - true_prob_mat
  #   mle_errors = mle - true_prob_mat
  #   plot_error_comparison(postmean_errors, mle_errors, parameter_name = "theta")
  # }

  return(list(
    "theta_postmean_rmse" = rMSE(postmean, true_prob_mat),
    "theta_mle_rmse" = rMSE(mle, true_prob_mat),
    "theta_postmean_TV" = TV(postmean, true_prob_mat) / ncol(D), # make it avg
    "theta_mle_TV" = TV(mle, true_prob_mat) / ncol(D)
  ))
}


# All in one wrapper for when multiple cutoffs
analyze_theta_for_multiple_cutoffs <- function(D,
                                               ng,
                                               out,
                                               true_prob_mat,
                                               burnin,
                                               upper_cutoffs
                                               ) {
  theta_list <- list()
  if (is.null(upper_cutoffs)){
      niter = length(out$eta_chain)
      upper_cutoffs = niter
  }

  for (i in seq_along(upper_cutoffs)) {
    u <- upper_cutoffs[i]
    theta_means <- R_get_thetameans1( #TODO: non-R version?
      D,
      ng,
      out,
      burnin = burnin,
      upper = u)
    theta_output <- get_theta_error_analysis(
      theta_means,
      D,
      true_prob_mat,
      burnin = burnin
    )
    if (i == 1) {
      theta_list["theta_mle_rmse"] <- theta_output$theta_mle_rmse
      theta_list["theta_mle_TV"] <- theta_output$theta_mle_TV
    }
    theta_list[[glue("theta_postmean{u}_rmse")]] <- theta_output$theta_postmean_rmse
    theta_list[[glue("theta_postmean{u}_TV")]] <- theta_output$theta_postmean_TV
  }
  return(theta_list)
}
