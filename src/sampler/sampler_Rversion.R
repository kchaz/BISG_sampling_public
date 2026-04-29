#' This script is for general helpers used by multiple samplers
#'
#' R_ before the function indicates I have also implemented an Rcpp version
#'
#' By: Kyla Chasalow
#------------------------------------------------------------

#' L is vector of indices (possibly in shuffled order)
#' r is position
#' Returns a pair of indices to run on
#'
R_get_ij <- function(L, r) {
  if (r + 1 <= length(L)) {
    i <- L[r]
    j <- L[r + 1]
  } else {
    i <- L[r]
    j <- L[1] # if odd # surnames, group last with first elem of L
  }
  return(c(i, j))
}

#' ----------------------------------------------------------------
#' alpha proposal generator and calculates likelihood ratio
#' ----------------------------------------------------------------
R_calculate_distance <- function(x1, y1, x2, y2) {
  return(sqrt((x2 - x1)^2 + (y2 - y1)^2))
}

R_sample_bounded_line <- function(ai, aj,
                                li = 0,
                                lj = 0,
                                ui = 1,
                                uj = 1){

  tot = ai+aj
  #slope = (tot - li - lj)/(li - tot + lj)
  #intercept = (tot - li) - slope * li
  slope = (-ai + li)/(ai-li)
  intercept = aj - ai*slope

  #Note: really slope is always -1 but below I just have this as a variable
  # So e.g. (uj - intercept)/slope could also just be intercept - uj

  # find intersection point where y = uj
  xcross_uj = (uj - intercept)/slope
  ycross_uj = uj#slope * xcross_uj + intercept

  # find the intersection point y = lj
  xcross_lj = (lj - intercept)/slope
  ycross_lj = lj#slope * xcross_lj + intercept

  # find intersection point x = li
  xcross_li = li
  ycross_li = slope * xcross_li + intercept

  # find intersection point x = ui
  xcross_ui = ui
  ycross_ui = slope * xcross_ui + intercept

  #Two candidates for upper are where hits li or where hits uj
  xcross_upper = max(xcross_li, xcross_uj)
  ycross_upper = min(ycross_li, ycross_uj)

  #Two candidates for lwoer are where hits lj or where hits ui
  xcross_lower = min(xcross_lj, xcross_ui)
  ycross_lower = max(ycross_lj, ycross_ui)

  #distances that can sample from
  distance_neg <- R_calculate_distance(xcross_upper, ycross_upper, ai, aj)
  distance_pos <- R_calculate_distance(xcross_lower, ycross_lower, ai, aj)

  # seems to work reasonably in practice
  sigma = mean(c(ai,aj))

  e = EnvStats::rnormTrunc(n = 1, mean = 0, sd = sigma, min = -distance_neg, max = distance_pos)
  theta = atan(abs(slope))
  x = e * cos(theta)
  y = e * sin(theta)
  ainew = ai + x #if e>0, should be -|e|cos(theta)  so do +
  ajnew = aj - y #if e<0, should be +|e|cos(theta)  so do -


  p_old2new = EnvStats::dnormTrunc(x = e, mean = 0, sd = sigma,
                                   min = -distance_neg, max = distance_pos)


  distance_neg <- R_calculate_distance(xcross_upper, ycross_upper, ainew, ajnew)
  distance_pos <- R_calculate_distance(xcross_lower, ycross_lower, ainew, ajnew)
  p_new2old = EnvStats::dnormTrunc(-e,0, sigma,
                                   min = -distance_neg, max = distance_pos)

  lratio = log(p_new2old) - log(p_old2new)

  return(list("new" = c(ainew, ajnew),
              "lratio" = lratio))
}


#' -------------------------------------------------------------------
#' c scaling up/down proposal generator, evaluater + eta prior function
#' -------------------------------------------------------------------
R_dcprior <- function(eta, lambda1 = 1, lambda2 = 1 / 10, log = T) {
  return(dgamma(eta, shape = lambda1, rate = lambda2, log = log))
}

R_eta_proposal_generator <- function(eta) {
  etanew <- eta * exp(rnorm(1, mean = 0, sd = 1))
  return(etanew)
}

R_eta_proposal_evaluator <- function(eta, etanew, log = T) {
  dnorm(log(etanew / eta), mean = 0, sd = 1, log = log)
}

#' ----------------------------------------------------------------
#' Accept ratios for marginal posterior (mp) (alpha c)
#' ----------------------------------------------------------------
# D can be Dij or D
# a can be aij or a
R_get_prod1 <- function(D, eta, a) {
  sum <- 0

  svals <- 1:length(a) # for alpha ratio only s in {i,j}
  gvals <- 1:dim(D)[2]

  for (s in svals) {
    for (g in gvals) {
      nsg <- D[s, g]
      if (nsg > 0) {
        v <- sum(log(0:(nsg - 1) + eta * a[s]))
      } else {
        v <- 0
      }
      sum <- sum + v
    }
  }
  return(sum)
}
R_get_prod2 <- function(D, eta, a) {
  gsums <- 0
  # get the sum of logs for each s for all the g's
  for (k in 1:length(a)) {
    gsums <- gsums + lgamma(D[k, ] + eta * a[k]) - lgamma(eta * a[k])
  }
  # sum over the g's at end
  return(sum(gsums))
}


R_mp_eta_accept_ratio <- function(D, ng, etanew, eta, a, lambda1, lambda2,
                                log = T, option = 2) {
  # preliminaries
  stopifnot(eta > 0 & etanew > 0)

  # prior part
  priorpart <- R_dcprior(etanew,
    lambda1 = lambda1,
    lambda2 = lambda2,
    log = T
  ) - R_dcprior(eta, lambda1 = lambda1, lambda2 = lambda2, log = T)

  # gamma part
  gammapart <- lgamma(etanew) - lgamma(eta)

  # probability ratio likelihood part 1
  if (option == 1) {
    p1 <- R_get_prod1(D, etanew, a) - R_get_prod1(D, eta, a)
  }
  if (option != 1) {
    p1 <- R_get_prod2(D, etanew, a) - R_get_prod2(D, eta, a)
  }

  # probability ratio likelihood part 2
  lognum <- sum(-lgamma(ng + etanew))
  logden <- sum(-lgamma(ng + eta))
  p2 <- lognum - logden

  lr1 <- priorpart + gammapart + p1 + p2

  # proposal ratio part
  lognum <- R_eta_proposal_evaluator(eta = etanew, etanew = eta, log = T)
  logden <- R_eta_proposal_evaluator(eta = eta, etanew = etanew, log = T)
  lr2 <- lognum - logden

  s <- lr1 + lr2
  if (log) {
    return(s)
  }
  if (!log) {
    return(exp(s))
  }
}



#'-------------------------------------------------------------------
#' SAMPLER 1 Main Helpers
#'-------------------------------------------------------------------
#' Group indices are indices of first entry in each pair (or more generally,
#' could imagine implementing something with groups)

R_run_alpha_cycle <- function(D, a, gam, eta, u, L, group_indices,
                              testmode, option = 2){

  num_accept <- 0

  for (r in group_indices) {

    # set i and j for this group
    ij <- R_get_ij(L, r)
    i <- ij[1]
    j <- ij[2]

    # Get current alphas
    ai <- a[i]
    aj <- a[j]
    gi <- gam[i]
    gj <- gam[j]
    Dij <- D[c(i, j), ]
    aij <- c(ai, aj)

    # sample new alphas. Comes with accept ratio in output
    proposal = R_sample_bounded_line(
                                ai = ai,
                                aj = aj,
                                li = 0,
                                lj = 0,
                                ui = 1,
                                uj = 1
                                )
    aijnew <- proposal$new
    ainew = aijnew[1]
    ajnew = aijnew[2]
    lratio_alpha_part = proposal[["lratio"]]

    if (testmode) {
      if (sum(a == 0) > 0) {
        warning("a has a 0 entry in call to mp_alpha_accept_ratio()")
      }
      stopifnot(i > 0 & i <= length(a))
      stopifnot(j > 0 & j <= length(a))
      stopifnot(round(sum(aijnew), 8) == round(sum(aij), 8))
    }

    # calculate ratio components from probability distribution
    #---------------------------------------------------------
    p1 <- (gi - 1) * (log(ainew) - log(ai)) + (gj - 1) * (log(ajnew) - log(aj))
    # Option 1 - slower but less risk of numerical overflow
    if (option == 1) {
      p2 <- R_get_prod1(Dij, eta, aijnew) - R_get_prod1(Dij, eta, aij)
    }
    # Option 2 - faster but more risk of numerical overflow
    if (option != 1) {
      p2 <- R_get_prod2(Dij, eta, aijnew) - R_get_prod2(Dij, eta, aij)
    }

    # calculate log(transition prob)
    lratio <- lratio_alpha_part + p1 + p2
    if (testmode) {
      if (is.nan(lratio)) {
        stop("a ratio is NaN!")
      }
    }

    # accept/reject with this prob and if accept, update a
    # if reject, do nothing because not saving within-cycle
    v <- log(runif(1))
    if (v < min(0, lratio)) {
      a[i] <- ainew
      a[j] <- ajnew
      num_accept <- num_accept + 1
    }
  }
  return(list("num_accept" = num_accept,
              "a" = a))
}


R_run_eta_update <- function(D, ng, a, eta, lambda1, lambda2,
                         testmode, option){

  accept_eta = 0

  # propose new value
  etanew <- R_eta_proposal_generator(eta)

  # calculate accept ratio
  lratio <- R_mp_eta_accept_ratio(D = D, ng = ng, etanew = etanew, eta = eta, a = a,
                              lambda1 = lambda1, lambda2 = lambda2, log = T,
                              option = option)
  if (testmode) {
    if (is.nan(lratio)) {
      stop("c ratio is NaN!")
    }
  }
  if (is.nan(lratio)) {
    lratio <- -1
  } # if  0/0 ratio, just reject?   #TODO - consider this

  #accept or reject with appropriate probability
  v <- log(runif(1))
  if (v < min(0, lratio)) {
    eta <- etanew
    accept_eta = 1
  }

  return(list(
              "accept_eta" = accept_eta,
              "eta" = eta))

}

#'----------------------------------------------------------------------
#' Main Function
#'----------------------------------------------------------------------
#' Note that initialization is saved as part of chain - first element
R_runMCMC1 <- function(D, ng, niter, init_a, init_eta,
                     lambda1, lambda2, gam,
                     testmode = F, verbose = F,
                     option = 2,
                     increment = 1, #only store every <increment> element of chains
                     means_burnin = 1 #at what point to start accumulating the posterior means - init alway discarded and not part of this count
                     ) {

  stopifnot(round(increment) == increment)

  num_sur <- nrow(D)
  group_indices <- seq(1, num_sur, 2) # entries mark first in every group from L below
  ngroups <- length(group_indices)
  L <- 1:num_sur

  # Calculate how big holders should be based on increment
  hsize = length(seq(1,niter, by = increment))

  # Holders for output, with initialization saved as part of chain
  eta_chain <- c(init_eta, rep(NA, hsize))
  a_chain <- matrix(NA, ncol = hsize + 1, nrow = num_sur)
  a_chain[, 1] <- init_a

  # Holders for accumulating theta posterior means
  theta_sums = matrix(0, nrow = num_sur, ncol = length(ng))
  sums_denom = 0

  # Trackers
  num_accept_a_per_iter <- c(ngroups, rep(0, niter))
  num_accept_eta = 0

  # Temporary holders
  alpha_old = init_a
  eta_old = init_eta
  ic = 1 #increment counter, start at 1 to account for time 0 storage

  # New element of the chain
  for (t in 2:(niter + 1)) {
    if (verbose) {
      if ((t - 1) %% 50 == 0) {
        print(paste("Iter:", t - 1))
      }
    }

    # alpha cycle update
    #--------------------
    out_a = R_run_alpha_cycle(D = D,
                            a = alpha_old, gam = gam,
                            eta = eta_old,
                            L = L, group_indices = group_indices,
                            testmode = testmode,
                            option = option)

    alpha_old = out_a[["a"]]
    num_accept_a_per_iter[t] <- out_a[["num_accept"]]
    if (testmode) {
      stopifnot(round(sum(alpha_old), 10) == 1)
    }

    # c update
    #---------
    out_eta = R_run_eta_update(D = D, ng = ng,
                         a = alpha_old ,
                         eta = eta_old,
                         lambda1 = lambda1,
                         lambda2 = lambda2,
                         testmode = testmode,
                         option = option)

    eta_old = out_eta[["eta"]]
    num_accept_eta = num_accept_eta + out_eta[["accept_eta"]]
    if (testmode) {
      stopifnot(out_eta[["eta"]] > 0)
    }

    # if hit an increment where going to save, save
    if ((t-2) %% increment == 0){  #b/c start at 2
      ic = ic + 1
      a_chain[,ic] <- alpha_old
      eta_chain[ic] <- eta_old
    }

    # Theta posterior means addition update - will never include init value
    #--------------------------------------
    if (t-2 >= means_burnin){
      numerator = apply(D, 2, function(col){col + eta_old * alpha_old})
      theta_sums = theta_sums +  t(t(numerator) / (ng + eta_old))
      sums_denom = sums_denom + 1
    }

    # shuffle order of L
    L <- sample(L)

  }
  return(list(
    "a_chain" = a_chain,
    "eta_chain" = eta_chain,
    "num_groups_per_a_cycle" = ngroups,
    "num_accept_per_a_cycle_iter" = num_accept_a_per_iter,
    "num_accept_eta" = num_accept_eta,
    "theta_sums" = theta_sums,
    "sums_denom" = sums_denom
  ))
}
