library(glue)

#-------------------------------------------------------
# Acceptance rate stats
#-------------------------------------------------------

get_sampler1_accept_stats <- function(out, burnin = 1){

  niter = length(out$eta_chain) - 1

  eta_rate <- out$num_accept_eta / niter

    a_info <- out$num_accept_per_a_cycle_iter
  n_groups <- out$num_groups_per_a_cycle
  a_mean_num_accept <- mean(a_info[burnin:niter])
  a_mean_prop_accept <- mean(a_info[burnin:niter] / n_groups)
  a_max_accept_prop <- max(a_info[burnin:niter] / n_groups)
  a_min_accept_prop <- min(a_info[burnin:niter] / n_groups)

  return(list("eta_rate" = eta_rate,
              "a_mean_num_accept" = a_mean_num_accept,
              "a_mean_prop_accept" = a_mean_prop_accept,
              "a_max_accept_prop" = a_max_accept_prop,
              "a_min_accept_prop" = a_min_accept_prop))
}



#-------------------------------------------
# Chain Plotters
#-------------------------------------------
# plot_eta_chain <- function(out, burnin = 1){
#   niter = length(out$eta_chain)-1
#   eta_chain = out$eta_chain[burnin:(niter+1)]
#   postmean = mean(eta_chain)
#   par(mfrow = c(1, 1))
#   plot(burnin:(niter + 1), eta_chain,
#        main = "eta",
#        "l",
#        xlab = "iterations",
#        ylab = "eta")
#   abline(h = postmean,col = "red")
#   return(postmean)
# }

plot_eta_chain <- function(out, burnin = 1, xlab = "Index") {
  niter <- length(out$eta_chain) - 1
  eta_chain <- out$eta_chain[burnin:(niter + 1)]
  postmean <- mean(eta_chain)
  par(
    mar = c(4.5, 4.5, 2.5, 1),
    las = 1,
    bty = "l"
  )
  plot(
    burnin:(niter + 1), eta_chain,
    type = "l",
    lwd = 1.5,
    col = "gray30",
    xlab = xlab,
    ylab = expression(eta),
    main = expression("Trace plot of " * eta),
    cex.main = 2,
    cex.lab = 2
  )
  abline(
    h = postmean,
    col = "firebrick",
    lwd = 1,
    lty = 1
  )
  invisible(postmean)
}





plot_a_chain = function(out, indices, burnin = 1, labels = NULL, xlab = "Index"){
  par(mfrow = c(3,3))
  niter = length(out$eta_chain) - 1
  a_chain = out$a_chain[,burnin:(niter+1)]

  if (is.null(labels)){
    titles = rownames(out$a_chain)[indices]
  } else {
    titles = labels
    stopifnot(length(labels) == length(indices))
  }
  lower <- min(a_chain)
  upper <- max(a_chain)

  for (j in 1:length(indices)){
    ind = indices[j]
    plot(burnin:(niter + 1), a_chain[ind, ], "l",
         main = titles[j],
         xlab = xlab,
         ylab = expression(alpha[s]),
         cex.lab = 1.5,
         cex.main = 1.5)
  }
}

#plots a single entry of the matrix
plot_theta_chain <- function(out, i, j, burnin = 1){

  niter = length(out$eta_chain) - 1

  #get rid of burn-in elements
  theta_chain = out$theta_chain[burnin:(niter+1)]

  # get values
  v <- sapply(theta_chain, function(mat) {mat[i, j]})
  lower <- min(v)
  upper <- max(v)
  plot(burnin:(niter + 1), v, "l",
       main = paste0("theta_",i,j),
       ylim = c(lower, upper))
}


#-----------------------------------------------
# Func to get posterior means for a and c and look at a MSE
# relative to true pop probs
#-----------------------------------------------
get_alpha_and_eta_postmeans <- function(out, burnin = 1, upper_cutoffs = NULL){

  if(is.null(upper_cutoffs)){

    l = length(out$eta_chain)

    eta_chain = out$eta_chain[burnin:l]
    a_chain = out$a_chain[,burnin:l]

    a_means = apply(a_chain, 1, mean)
    eta_mean = mean(eta_chain)
    return(list("a_means" = a_means,
                "eta_means" = eta_mean))

  }
  else{

    upper_cutoffs = sort(upper_cutoffs)

    num_sur = dim(out$a_chain)[1]

    nc = length(upper_cutoffs)
    a_mean_holder = matrix(NA, ncol = nc, nrow = num_sur)
    eta_mean_holder = matrix(NA, ncol = nc, nrow = 1)

    for(i in seq_along(upper_cutoffs)){

      u = upper_cutoffs[i]

      if (u < burnin){
        stop("upper cutoff is smaller than burnin in get_alpha_and_eta_postmeans")
      }

      eta_chain = out$eta_chain[burnin:u]
      a_chain = out$a_chain[,burnin:u, drop = F]

      a_mean_holder[,i] = apply(a_chain, 1, mean)
      eta_mean_holder[,i] = mean(eta_chain)

    }
    colnames(a_mean_holder) = upper_cutoffs
    rownames(a_mean_holder) = rownames(out$a_chain)
    colnames(eta_mean_holder) = upper_cutoffs
    return(list("a_means" = a_mean_holder,
                "eta_means" = eta_mean_holder))

  }

}


#---------------------------------
# Sampler 1 get postmean for theta
# These two have Rcpp versions too
# but there seems to still be a bug that sometimes makes them crash...
#---------------------------------

# Old less efficient version
# R_get_postmean_mat1_old <- function(D, c, a, ng){
#   num = apply(D, 2, function(col){col + c * a})
#   means = t(t(num) / (ng + c)) #divide each column
#   return(means)
# }


R_get_postmean_mat1 <- function(D, eta, alpha, ng) {
  # D: |S| x |G|
  # alpha: length |S| (rows)
  # ng: length |G| (cols)
  sweep(D + eta * alpha, 2, ng + eta, "/") # for each column of (arg1) divide by corresponding value of (arg2)
}


R_get_thetameans1 <- function(D, ng, out, burnin = 1, upper = NULL){

  holder = matrix(0, ncol = ncol(D), nrow = nrow(D))

  if(is.null(upper)){
    upper = length(out$eta_chain)
  } else{
    if(burnin > upper){
      stop("Burnin is greater than upper limit")
    }
  }

  for (t in burnin : upper){
    a = out$a_chain[,t]
    eta = out$eta_chain[t]
    holder = holder + get_postmean_mat1(D, eta, a, ng)
  }
  return(holder/(upper - burnin + 1))
}



#------------------------------------
# analyze eta via rho
#------------------------------------

#' calculate posterior distribution of rho based on eta chain
analyze_rho <- function(out, ng, burnin = 1, upper = NULL){

  if (is.null(upper)) {
    upper <- length(out$eta_chain)
  } else {
    if (burnin > upper) stop("Burnin is greater than upper limit")
    if (upper > length(out$eta_chain)) stop("Upper limit exceeds length of eta_chain")
  }
  if (burnin < 1) stop("burnin must be >= 1")

  # Subset chain indices used for posterior summaries
  idx <- burnin:upper

  # Ensure ng is named (for column names)
  if (is.null(names(ng)) || any(names(ng) == "")) {
    stop("ng must be a named vector (names used as state/column names).")
  }


  # Extract eta chain (assumed numeric vector)
  eta_chain <- out$eta_chain
  if (!is.numeric(eta_chain)) stop("out$eta_chain must be numeric.")
  if (length(eta_chain) == 0) stop("out$eta_chain is empty.")

  # create a matrix rho_mat with column for each state (with names from names of ng vector)
  # and row for each element of chain and value is (eta/ng) / (1+ eta/ng)
  #
  # This simplifies to: eta / (eta + ng)
  rho_mat <- outer(
    eta_chain[idx],
    ng,
    FUN = function(eta, ng_i) eta / (eta + ng_i)
  )
  colnames(rho_mat) <- names(ng)

  # obtain posterior mean of rho for each state by averaging over chain, subsetting to burnin:upper
  rho_postmean <- colMeans(rho_mat)

  # obtain the mean rho across each state for each iteration
  overall_rho_mean_chain = rowMeans(rho_mat)
  overall_rho_sd_chain = apply(rho_mat, 1, sd)
  overall_rho_min_chain = apply(rho_mat, 1, min)
  overall_rho_max_chain = apply(rho_mat, 1, max)

  postmean_rho_mean = mean(overall_rho_mean_chain)
  postmean_rho_sd = mean(overall_rho_sd_chain)
  postmean_rho_min = mean(overall_rho_min_chain)
  postmean_rho_max = mean(overall_rho_max_chain)

  # Weighted mean
  rho_postmean_over_individuals = sum(rho_postmean * ng) / sum(ng)

  # return rho mat and the posterior means
  list(
    rho_mat = rho_mat,
    rho_postmean = rho_postmean,
    postmean_rho_mean = postmean_rho_mean,
    postmean_rho_sd  =  postmean_rho_sd,
    postmean_rho_min = postmean_rho_min,
    postmean_rho_max = postmean_rho_max,
    rho_postmean_over_individuals = rho_postmean_over_individuals
  )
}




plot_rho_chain <- function(rho_mat, state, xlab = "Index") {
  niter <- length(out$eta_chain) - 1
  par(
    mar = c(4.5, 4.5, 2.5, 1),
    las = 1,
    bty = "l"
  )
  plot(
    1:nrow(rho_mat), rho_mat[,state],
    type = "l",
    lwd = 1.5,
    col = "gray30",
    xlab = xlab,
    ylab = expression(rho),
    main = bquote("Trace plot of " ~ rho ~ " for " ~ .(state)),
    cex.main = 2,
    cex.lab = 2
  )
  abline(
    h = mean(rho_mat[,state]),
    col = "firebrick",
    lwd = 1,
    lty = 1
  )
}




