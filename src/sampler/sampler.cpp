#include <Rcpp.h>
using namespace Rcpp;



//---------------
// Minor helpers
//---------------
// [[Rcpp::export]]
NumericVector get_ij(NumericVector L, int r) {
  int n = L.size();
  int i, j;

  if (r + 1 <= n) {
    i = L[r - 1];  // Adjust for 0-based index in C++
    j = L[r];      // Adjust for 0-based index in C++
  } else {
    i = L[r - 1];  // Adjust for 0-based index in C++
    j = L[0];      // Get the first element if r + 1 exceeds length of L
  }

  return NumericVector::create(i, j);
}

// if niter/increment has remainder, need + 1 for that
// [[Rcpp::export]]
int calculate_numsave(int niter, int increment) {
  int result = niter / increment;
  //Rcpp::Rcout << result << std::endl;
  if (niter % increment != 0){
    result = result + 1;
  }
  return result;
}

//--------------------------------------------------------------
// Functions to initialize various kinds of vectors and matrices
//--------------------------------------------------------------

// function to create a vector of length n with elements e  - Note: this is kind of redundant
// [[Rcpp::export]]
NumericVector create_vec(double e, int n) {
  NumericVector vec(n, e);
  return vec;
}


// function to create a vector of length n with elements 1:n
// [[Rcpp::export]]
NumericVector create_vec_1ton(int n) {
  NumericVector vec(n);
  for (int i = 0; i < n; ++i) {
    vec[i] = i + 1;
  }
  return vec;
}

// function to create a vector of length n with elements 1:n by 2's
// [[Rcpp::export]]
NumericVector create_vec_1ton_odds(int n) {

  int l;
  if (n % 2 == 0){
    l = n/2;
  } else{
    l = (n+1)/2;
  }
  NumericVector vec(l);
  for (int i = 0; i < l; ++i) {
       vec[i] = 2*i+1;
  }
  return vec;
}


// function to create a matrix of NAs
// [[Rcpp::export]]
NumericMatrix create_empty_mat(int nrow,int ncol) {

  NumericMatrix mat(nrow, ncol); //create
  std::fill(mat.begin(), mat.end(), NA_REAL); //fill
  return mat;
}

// function to create a matrix of 0s
// [[Rcpp::export]]
NumericMatrix create_zero_mat(int nrow,int ncol) {

  NumericMatrix mat(nrow, ncol); //create
  std::fill(mat.begin(), mat.end(), 0); //fill
  return mat;
}


// function to randomly shuffle a vector
// [[Rcpp::export]]
NumericVector shuffle_vector(NumericVector L) {
  NumericVector shuffled = Rcpp::sample(L, L.size(), false, R_NilValue);
  return shuffled;
}


// [[Rcpp::export]]
NumericMatrix create_zero_matrix(int n, int m) {
  NumericMatrix mat(n, m); // Initialize a NumericMatrix of size n by m

  // Fill the matrix with 0's
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      mat(i, j) = 0.0;
    }
  }

  return mat;
}



// alpha proposal generator and evaluator
// ---------------------------------------------------

// [[Rcpp::export]]
double calculate_distance(double x1, double y1, double x2, double y2) {
  return std::sqrt(std::pow(x2 - x1, 2.0) + std::pow(y2 - y1, 2.0));
}


// Create rnorm (generator) and dnorm (evaluator) functions for truncated normal

//TODO: UNDERSTAND EVERYTHING BELOW

// Quick wrappers for standard normal pnorm, qunorm, dnorm
static inline double pnorm_std(double x) {
  return R::pnorm(x, 0.0, 1.0, 1, 0); // value, mean, sd, lower_tail = 1 means P(Z<= z), log_p =0 means return prob, not log
}

static inline double qnorm_std(double p) {
  return R::qnorm(p, 0.0, 1.0, 1, 0);
}

static inline double dnorm_std(double x) {
  return R::dnorm(x, 0.0, 1.0, 0);
}

// Random generator from truncated normal
static inline double rnorm_trunc(double mu,
                                 double sigma,
                                 double lower,
                                 double upper) {
  // Standardize the bounds to N(0,1) space
  double a = (lower - mu) / sigma;
  double b = (upper - mu) / sigma;

  // Calculate P(Z<=a) and P(Z<=b)
  double pa = pnorm_std(a);
  double pb = pnorm_std(b);

  // Guard against edge cases like if sigma is 0 so that a and b end up infinite
  if (!R_finite(pa) || !R_finite(pb) || pb <= pa) {
    return mu;
  }

  // Use PIT (U~Unif -> F^{-1}(U) ~ F) - see my truncnorm notes for why this works
  double u = R::runif(pa, pb);
  double z = qnorm_std(u);
  return mu + sigma * z; // Go back to original scale
}

// Calculate density for truncnorm
static inline double dnorm_trunc(double x,
                                 double mu,
                                 double sigma,
                                 double lower,
                                 double upper) {
  double a = (lower - mu) / sigma;
  double b = (upper - mu) / sigma;

  double pa = pnorm_std(a);
  double pb = pnorm_std(b);

  double C = pb - pa; // normalizing constant to divide by
  if (!R_finite(C) || C <= 0.0) return NA_REAL;

  double z = (x - mu) / sigma;
  return (dnorm_std(z) / sigma) / C; // comes from change of var formula
}

// [[Rcpp::export]]
Rcpp::List sample_bounded_line(double ai, double aj,
                                     double li = 0.0, double lj = 0.0,
                                     double ui = 1.0, double uj = 1.0) {

  double tot = ai + aj;
  //double slope = (tot - li - lj) / (li - tot + lj);
  //double intercept = (tot - li) - slope * li;
  double slope = (-ai + li)/(ai-li);
  double intercept = aj - ai*slope;

  double xcross_uj = (uj - intercept) / slope;
  double ycross_uj = uj; //slope * xcross_uj + intercept;

  double xcross_lj = (lj - intercept) / slope;
  double ycross_lj = lj;// * xcross_lj + intercept;

  double xcross_li = li;
  double ycross_li = slope * xcross_li + intercept;

  double xcross_ui = ui;
  double ycross_ui = slope * xcross_ui + intercept;

  double xcross_upper = std::max(xcross_li, xcross_uj);
  double ycross_upper = std::min(ycross_li, ycross_uj);

  double xcross_lower = std::min(xcross_lj, xcross_ui);
  double ycross_lower = std::max(ycross_lj, ycross_ui);

  double distance_neg = calculate_distance(xcross_upper, ycross_upper, ai, aj);
  double distance_pos = calculate_distance(xcross_lower, ycross_lower, ai, aj);

  double sigma = 0.5 * (ai + aj);

  double e = rnorm_trunc(0.0, sigma, -distance_neg, distance_pos);
  double theta = std::atan(std::abs(slope));

  double x = e * std::cos(theta);
  double y = e * std::sin(theta);

  double ainew = ai + x;
  double ajnew = aj - y;

  double p_old2new = dnorm_trunc(e, 0.0, sigma, -distance_neg, distance_pos);

  distance_neg = calculate_distance(xcross_upper, ycross_upper, ainew, ajnew);
  distance_pos = calculate_distance(xcross_lower, ycross_lower, ainew, ajnew);

  double p_new2old = dnorm_trunc(-e, 0.0, sigma, -distance_neg, distance_pos);

  double lratio = std::log(p_new2old) - std::log(p_old2new);

  return Rcpp::List::create(
    Rcpp::Named("new") = Rcpp::NumericVector::create(ainew, ajnew),
    Rcpp::Named("lratio") = lratio
  );
}



// eta proposal generator and evaluator and prior
// ---------------------------------------------------

// gamma prior, default exp(1/10)
// [[Rcpp::export]]
double dcprior(double eta, double lambda1 = 1, double lambda2 = 1/10, bool log = true){
  // int l;
  // if (log){
  //   l = 1;
  // } else {
  //   l = 0;
  // }
  int log_flag = log ? 1 : 0;
  double p = R::dgamma(eta, lambda1, 1/lambda2, log_flag); //seems to use different parameterization
  return(p);
}

// making sigma larger decreases amount of movement in chain but makes it converge faster
// < 1 seemed to converge too slowly in some tests
// [[Rcpp::export]]
double eta_proposal_generator(double eta) {
  double etanew = eta * exp(R::rnorm(0, 1)); // Generating random number with mean 0 and sd 1
  return etanew;
}

// [[Rcpp::export]]
double eta_proposal_evaluator(double eta, double etanew, bool log) {
  int log_flag = log ? 1 : 0;
  double log_density = R::dnorm(std::log(etanew / eta), 0, 1, log_flag);
  return log_density;
}



// Helpers for evaluating accept ratios
//-------------------------------
// D can be Dij or D
// a can be aij or a
// Version 1 uses simplifications of the ratios of gamma functions
// Version 2 just does the (log) gamme functions
//
// This is all done entirely on log scale
// [[Rcpp::export]]
double get_prod1(NumericMatrix D, double eta, NumericVector a) {
  int n = D.nrow();
  int m = D.ncol();

  double sum = 0.0;

  // for each row (surname)
  for (int s = 0; s < n; ++s) {

    // for each column (state)
    for (int g = 0; g < m; ++g) {

      // want to go for k = 1,...,n_{sg} count for that pair
      // but actually since have n_{sg}-k, translates also to doing
      // k = 0,...,n_{sg}-1 which happens naturally below because
      // of 0-indexing of C++
      int nsg = D(s, g);
      double v;

      if (nsg > 0) {
        v = 0.0;
        for (int k = 0; k < nsg; ++k) {
          v += std::log(k + eta * a[s]);
        }
      } else {
        v = 0.0;
      }
      sum += v;
    }
  }
  return sum;
}


// [[Rcpp::export]]
double get_prod2(NumericMatrix D,
                 double eta,
                 NumericVector a) {

  int n = D.nrow();
  int m = D.ncol();

  double gsums = 0.0;

  // for each fixed s (row entry)
  for (int k = 0; k < n; ++k) {
    double sum_log_gamma = 0.0;

    // for each g (state entry)
    for (int j = 0; j < m; ++j) {
      sum_log_gamma += std::lgamma(D(k, j) + eta * a[k]) - std::lgamma(eta * a[k]);
      //Rcpp::Rcout << "sum_log_gamma: " << sum_log_gamma << std::endl;

    }

    gsums += sum_log_gamma;
  }

  return gsums;
}



// [[Rcpp::export]]
double mp_eta_accept_ratio(NumericMatrix D, NumericVector ng,
                         double etanew,
                         double eta,
                         NumericVector a,
                         double lambda1,
                         double lambda2,
                         bool log =  true,
                         int option = 2){

  double priorpart = dcprior(etanew, lambda1, lambda2, log) - dcprior(eta, lambda1, lambda2, log);
  double gammapart = std::lgamma(etanew) - std::lgamma(eta);

  // probability like ratio part 1
  double p1;
  if (option == 1) { // slightly slower, but less risk of overflow
    p1 = get_prod1(D, etanew, a) - get_prod1(D, eta, a);
  } else { //slightly faster, more risk of overflow
    p1 = get_prod2(D, etanew, a) - get_prod2(D, eta, a);
  }

  // probability like ratio part 2
  double lognum = -sum(Rcpp::lgamma(ng + etanew));                           //NTS: so you can use vectorized lgamma...
  double logden = -sum(Rcpp::lgamma(ng + eta));
  double p2 = lognum - logden;

  double lr1 = priorpart + gammapart + p1 + p2;
  //Rcpp::Rcout << "lr1: " << lr1 << std::endl;

  //NTS: could cut this - I have symmetric proposal where this always 0
  // proposal ratio part
  double lognum2 = eta_proposal_evaluator(etanew, eta, true);
  double logden2 = eta_proposal_evaluator(eta, etanew, true);
  double lr2 = lognum2 - logden2;
  //Rcpp::Rcout << "lr2: " << lr2 << std::endl;

  // putting them together
  double result = lr1 + lr2;
  if (log) {
    return(result);
  } else {
    return(std::exp(result));
  }
}

// Functions for running components of Sampler 1
//------------------------------------------------
// [[Rcpp::export]]
List run_alpha_cycle(NumericMatrix D,
                     NumericVector a,
                     NumericVector gam,
                     double eta,
                     NumericVector L,
                     NumericVector group_indices,
                     bool testmode = false,
                     int option = 2){

  if (testmode) {
    if (sum(a == 0) > 0) {
      Rcpp::warning("a has a 0 entry in call to run_alpha_cycle");
    }
  }

  int num_accept = 0;

  for (int r : group_indices) {

    // set and i and j for this group
    NumericVector ij = get_ij(L,r);
    int i = ij[0];
    int j = ij[1];

    // adjust for 0-based indexing
    i--;
    j--;

    // Get current alphas
    double ai = a[i];
    double aj = a[j];
    double gi = gam[i];
    double gj = gam[j];
    NumericVector aij = NumericVector::create(ai, aj);

    // bit more effort than in R version to get i and j rows of D
    int p = D.ncol();
    NumericMatrix Dij(2, p);
    for (int col = 0; col < p; ++col) {
      Dij(0, col) = D(i, col);
      Dij(1, col) = D(j, col);
    }

    // sample new alphas. Comes with accept ratio in output
    // For now just use 0 and 1 bounds on alpha - have not yet implemented
    // custom bounds
    List proposal = sample_bounded_line(ai, aj, 0.0, 0.0, 1.0, 1.0);

    NumericVector aijnew = proposal["new"];
    double ainew = aijnew[0];
    double ajnew = aijnew[1];
    double lratio_alpha_part = proposal["lratio"]; //log accept ratio componant

    if (testmode) {
      if (std::round( (ainew + ajnew) * 1e8) != std::round( (ai+aj)*1e8)) {
        Rcpp::warning("Round equality check failed within run_alpha_cycle()");
      }
    }

    // calculate ratio components from probability distribution
    double p1 =
      (gi - 1.0) * (std::log(ainew) - std::log(ai)) +
      (gj - 1.0) * (std::log(ajnew) - std::log(aj));

    double p2;
    if (option == 1) {
      p2 = get_prod1(Dij, eta, aijnew) - get_prod1(Dij, eta, aij);
    } else {
      p2 = get_prod2(Dij, eta, aijnew) - get_prod2(Dij, eta, aij);
    }

    double lratio = lratio_alpha_part + p1 + p2;

    if (testmode) {
      if (R_IsNaN(lratio)) {
        Rcpp::stop("a ratio is NaN!");
      }
    }

    // Accept/Reject with prob = min(0, lratio) (because on log scale)
    // if reject, don't have to save because not saving within-cycle updates
    double v = std::log(R::runif(0,1));

    if (v < std::min(0.0, lratio)) {
       a[i] = ainew; // WARNING - this modifies the underlying a object, not just within func
       a[j] = ajnew;
       num_accept += 1;
    }
  }

  return List::create(
         _["num_accept"] = num_accept,
         _["a"] = a
       );
}


// [[Rcpp::export]]
List run_eta_update(NumericMatrix D,
                  NumericVector ng,
                  NumericVector a,
                  double eta,
                  double lambda1, double lambda2,
                  bool testmode = false,
                  int option = 2){


  // propose a new value
  double etanew = eta_proposal_generator(eta);
  //Rcpp::Rcout << "etanew:" << etanew << std::endl;

  // calculate ratio
  double lratio = mp_eta_accept_ratio(D,ng,etanew,eta,a,lambda1,lambda2,true,option);

  // deal with NaN case if 0/0 ratio for some reason
  if (NumericVector::is_na(lratio)){
    Rf_warning("eta ratio is NaN!");
    return List::create(
      _["accept_eta"] = 0,
      _["eta"] = eta
    );
  }
  // accept or reject based on lratio
  double v = std::log(R::runif(0,1));
  if (v < std::min(0.0, lratio)) {
    return List::create(
      _["accept_eta"] = 1,
      _["eta"] = etanew
    );
  }
  // else return the old one
  return List::create(
    _["accept_eta"] = 0,
    _["eta"] = eta
  );

}

//-------------------------------------------------------------------------
// Main Function: Sampler 1
//-------------------------------------------------------------------------
// [[Rcpp::export]]
List runMCMC1(NumericMatrix D,
              NumericVector ng,
              int niter,
              NumericVector init_a,
              double init_eta,
              double lambda1, double lambda2,
              NumericVector gam,
              bool testmode = false,
              bool verbose = false,
              int option = 2,
              int increment = 1,
              int means_burnin = 1){

  // initialize meta parameters
  int num_sur = D.nrow();
  int num_state = D.ncol();

  // create group_indices as seq(1, num_sur, 2) and L as 1:num_sur
  NumericVector group_indices = create_vec_1ton_odds(num_sur);
  int ngroups = group_indices.size();
  NumericVector L = create_vec_1ton(num_sur);

  // figure out length of chain based on niter and increment
  int hsize = calculate_numsave(niter, increment);

  // set up c chain
  NumericVector eta_chain(1 + hsize); // +1 for the initial
  eta_chain[0] = init_eta;
  std::fill(eta_chain.begin() + 1, eta_chain.end(), NA_REAL); //fill eta_chain with NA for rest

  // set up a chains
  NumericMatrix a_chain = create_empty_mat(num_sur, hsize + 1);
  for (int i = 0; i < num_sur; ++i) { //fill in initial
      a_chain(i, 0) = init_a[i];
  }

  // Trackers
  NumericVector num_accept_a_per_iter = create_vec(0,niter);
  int num_accept_eta = 0;

  // Holders for Accumulating Sums and denominator to calculate mean
  int sums_denom = 0;
  NumericMatrix theta_sums = create_zero_mat(num_sur, num_state);

  // Temporary Holders
  double eta_old = init_eta;
  NumericVector alpha_old = init_a;
  int ic = 0;

  // MAIN
  for (int t = 1; t <= niter; ++t) { //starts at 1 because 0 indexing and then <= niter is enough to go niter
    //Rcpp::Rcout << "Lstart" << L << std::endl;

    // counter
    if (verbose) {
      if (t % 50 == 0) {
        std::cout << "Iter: " << t << std::endl;
      }
    }

    // alpha cycle run on most recent a value
    List out_a = run_alpha_cycle(D, alpha_old, gam, eta_old, L, group_indices, testmode, option);

    // update 'old' alpha
    alpha_old = out_a["a"];
    num_accept_a_per_iter[t-1] = out_a["num_accept"];

    // c update
    List out_eta = run_eta_update(D, ng, alpha_old, eta_old, lambda1, lambda2, testmode, option);
    eta_old = out_eta["eta"];
    int accept_eta = out_eta["accept_eta"];
    num_accept_eta = num_accept_eta + accept_eta;  //for some reason cannot smoosh with prev

    // if increment is right, store values of a and eta
    if ((t-1) % increment == 0){
      ic = ic + 1;
      // store a value
      NumericVector atemp = out_a["a"];
      for (int i = 0; i < num_sur; ++i) { //fill in initial
        a_chain(i, ic) = atemp[i];
      }
      // store c value
      eta_chain[ic] = out_eta["eta"];
    }

    // accumulate mean if past burnin period
    if (t-1 >= means_burnin){
       for(int r = 0; r < num_sur; ++r){
         for(int f = 0; f < num_state; ++f){
            theta_sums(r,f) = theta_sums(r,f) + (D(r,f) + eta_old * alpha_old[r])/(ng[f] + eta_old);
         }
       }
       sums_denom = sums_denom + 1;
    }

    // shuffle L
    L = shuffle_vector(L);
  }

  return List::create(
    _["a_chain"] = a_chain,
    _["eta_chain"] = eta_chain,
    _["num_groups_per_a_cycle"] = ngroups,
    _["num_accept_per_a_cycle_iter"] = num_accept_a_per_iter,
    _["num_accept_eta"] = num_accept_eta,
    _["theta_sums"] = theta_sums,
    _["sums_denom"] = sums_denom
  );
}







//-----------------------------------------------------------------
// Output analysis functions
//-----------------------------------------------------------------

// Function to compute postmean matrix
// [[Rcpp::export]]
NumericMatrix get_postmean_mat1(const NumericMatrix& D,
                                double eta,
                                const NumericVector& a,
                                const NumericVector& ng) {
  int nrow = D.nrow();
  int ncol = D.ncol();

  NumericMatrix num(nrow, ncol);

  // D is num_sur x num_state
  // Compute num = D + eta * a (broadcast a to each column of D)
  for (int j = 0; j < ncol; ++j) {
    for (int i = 0; i < nrow; ++i) {
      num(i, j) = D(i, j) + eta * a[i];
    }
  }

  // Compute means by dividing each column by (ng + c)
  NumericMatrix means(nrow, ncol);
  for (int j = 0; j < ncol; ++j) {
    double divisor = ng[j] + eta;
    for (int i = 0; i < nrow; ++i) {
      means(i, j) = num(i, j) / divisor;
    }
  }
  return means;
}

// [[Rcpp::export]]
NumericMatrix get_thetameans1(const NumericMatrix& D,
                              const NumericVector& ng,
                              const List& out,
                              int burnin = 1) {
  int nrow_D = D.nrow();
  int ncol_D = D.ncol();

  // get objects from output
  NumericMatrix a_chain = out["a_chain"];
  NumericVector eta_chain = out["eta_chain"];
  int niter = eta_chain.size();

  // Initialize holder to zero
  NumericMatrix holder(nrow_D, ncol_D);
  holder.fill(0);
  holder.attr("rownames") = D.attr("rownames");
  holder.attr("colnames") = D.attr("colnames"); // This isn't working for some reason

  // Loop over the iterations
  for (int t = burnin - 1; t < niter; ++t) {
    NumericVector a = a_chain(_, t);
    double eta = eta_chain[t];
    NumericMatrix post_mean_mat = get_postmean_mat1(D, eta, a, ng);

    for (int i = 0; i < nrow_D; ++i){
      for (int j = 0; j < ncol_D; ++j){
        holder(i,j) = holder(i,j) + post_mean_mat(i,j);
      }
    }
  }

  // Normalize the holder matrix to get overall mean
  int num_samples = niter - burnin + 1;

  for (int i = 0; i < nrow_D; ++i){
    for (int j = 0; j < ncol_D; ++j){
      holder(i,j) = holder(i,j)/num_samples;
    }
  }

  return holder;
}
