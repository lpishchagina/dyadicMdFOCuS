#include <Rcpp.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullVertexSet.h>
#include <libqhullcpp/PointCoordinates.h>
using namespace Rcpp;
using namespace orgQhull;

/*Function returns integer number:
 1  - gaussian cost, FOCuS0;
 11 - gaussian cost, FOCuS;
 2  - poisson cost, FOCuS0;
 22 - poisson cost, FOCuS;
 */
unsigned int get_costType(std::string method, 
                          std::string cost) {
  unsigned int costType = 0;
  if (cost == "gauss") costType = 1;
  if (cost == "poisson") costType = 2;
  if (method == "FOCuS") costType = costType * 11;
  return costType;
}

/*Function calculates matrix cumsumMatrix ( n+1 rows x (p+1) columns )
 * first line by default, (0,...,0)
 * row i :(i, (x_1^k+...x_i^k)), k=1,..,p i=1,..., n
 * */
void get_cumsumMatrix(unsigned int p, 
                      unsigned int length, 
                      double** &cumsumMatrix,  
                      Rcpp::NumericMatrix data) {
  for (size_t k = 0; k < p; k++) cumsumMatrix[0][k + 1] =  data(0, k); 
  for (size_t i = 0; i < length; i++) cumsumMatrix[i][0] = i + 1; // column 0 : (1,..,n)
  //cumsum
  for (size_t i = 1; i < length; i++) 
    for (size_t k = 1; k <= p; k++) 
      cumsumMatrix[i][k] = cumsumMatrix[i-1][k] + data(i, k - 1);
  
 /*Print------------------------------------
  for (size_t i = 0; i < length; i++) {
   Rcpp::Rcout<< endl;
   for (size_t k = 0; k <= p; k++) 
     Rcpp::Rcout<<cumsumMatrix[i][k]<< " ";
  }
  ----------------------------------------*/
}



/* Function returns the value of cost function for FOCuS0 or FOCuS */
double get_cost(unsigned int typeCost, 
                unsigned int p, 
                unsigned int change, 
                unsigned int t, 
                double** &cumsumMatrix) {
  //Focus0 gauss
  double cost =  0;
  //memory
  double* right_cumsum = new double[p+1];
  //[t-change, S(t)-S(change)]
  for (size_t k = 0; k <= p; k++)
    right_cumsum[k] = cumsumMatrix[t][k] - cumsumMatrix[change][k];
  //Gaussian model
  if (typeCost == 1 || typeCost == 11) {
    //FOCuS0 +
    for (size_t k = 1; k <= p; k++)
      cost = cost - right_cumsum[k] * right_cumsum[k] / right_cumsum[0];
    //FOCuS
    if (typeCost == 11)
      for (size_t k = 1; k <= p; k++)
        cost = cost - cumsumMatrix[change][k] * cumsumMatrix[change][k] / (cumsumMatrix[change][0]);
  }
  //Poisson model :check!
  if (typeCost == 2 || typeCost == 22) {
    //FOCuS0
    for (size_t k = 1; k <= p; k++)
      if(right_cumsum[k] != 0)
        cost = cost + right_cumsum[k] * ( 1 - log(right_cumsum[k] / right_cumsum[0]));
      //FOCuS
      if (typeCost == 22)
        for (size_t k = 1; k <= p; k++)
          if(cumsumMatrix[change][k] != 0)
            cost = cost + cumsumMatrix[change][k] * (1 - log(cumsumMatrix[change][k] / (cumsumMatrix[change][0])));
        cost = 2*cost; 
  }
  //clean memory
  delete []right_cumsum;
  return cost;
}

/* Function obtains candidate's indices, 
 * builds the convex hull 
 * and removes candidate's indices that are not vertices of convex hull by dyadic procedure*/
void doPruning_dyadic(unsigned int qmin,
                      unsigned int n, 
                      unsigned int p, 
                       std::list<int> &candidates,              //i in Tau 
                       bool* & isPruning,                       //if true -remove i
                       double** & cumsumMatrix) { 
std::vector<coordT> points;                           // points P(i) in one line
  unsigned int q = qmin;
  unsigned int left_bound_chunk = n;
  unsigned int nb_points = 0;
  int index;
  int chunk_size = 0;
  //loop
  
  while ((n % (1 << q) == 0)  && (q <= int(log2(n))) ) {  // + 1 because counting starts from 1 and not from 0
    nb_points = 0;
    points.clear();
    left_bound_chunk  = n- (std::pow(2, q));
    //pre-processing for Qhull
   // flag: do pruning for all candidates at time t 
    auto rit = candidates.rbegin();
    while ((*rit >= left_bound_chunk) && (rit != candidates.rend())) {
      isPruning[*rit] = true;  // flag: do pruning for all candidates at time t in chunk
      for (unsigned int k = 0; k <= p; ++k) 
        points.push_back(cumsumMatrix[*rit][k]); //get points P(i) i in Tau'
      ++nb_points;
      ++rit;
    }
    //we need to check that the number of points in two chunks is more then in the first chunk after qhull
    if (chunk_size < nb_points) {     //do pruning
      //do Qhull
      orgQhull::Qhull qhull;
      qhull.runQhull("", p + 1, nb_points, points.data(), "s");
      //get all vertices of Qhull
      const orgQhull::QhullVertexList& vertices = qhull.vertexList();
      //Attention :  vertex.point().coordinates()[0] is always candidate's "index+1" 
      for (const orgQhull::QhullVertex& vertex : vertices) {
        const double* coords = vertex.point().coordinates();
        index = int(vertex.point().coordinates()[0]) - 1;
        if (isPruning[index] == true) isPruning[index] = false;// change pruning's flag : this candidate exists at time t
      }
      //remove interior points
      auto pruning_it = candidates.rbegin();
      while ((*pruning_it >= left_bound_chunk) && (pruning_it != candidates.rend())) {
        if (isPruning[*pruning_it]) {
          nb_points = nb_points - 1;
          pruning_it = decltype(pruning_it)(candidates.erase(std::next(pruning_it).base()));
        } else {
          ++pruning_it;
        }
      }
    }
    chunk_size = nb_points; 
    q = q + 1;
  }
}




//' @title dyadMdFOCuS
//'
//' @description Single Changepoint Detection using methods: MdFOCuS0 and MdFOCuS.
//' @param data is a matrix of data (n rows x p-columns, where p - dimension, n - length).
//' @param method is the algorithm: 'FOCuS0'or 'FOCuS'.
//' @param cost is the type of cost: 'gauss' or 'poisson'.
//' @param common_difference_step is the beta parameter for 'next'.
//' @param common_ratio_step is the alpha parameter for 'next'.
//' @param first_step_qhull is the parameter for the first next (next = first_step_qhull + p+1).
//' @param cand_nb is the boolean parameter ('cand_nb=true' => get candidate number)
//' @param opt_changes is the boolean parameter ('opt_changes=true' => get optimal changes)
//' @param opt_costs is the boolean parameter ('opt_costs=true' => get optimal costs)
//' @param cands is the boolean parameter ('cands=true' => get candidates)
//'
//' @return a list of  elements  = (changes, means, UnpenalizedCost, NumberOfCandidates).
//'
//' \describe{
//' \item{\code{change}}{is the changepoint.}
//' \item{\code{cost}}{is a number equal to the segmentation cost.}
//' \item{\code{pre_change_mean}}{is the vector of pre-change means for the p-variate time series.}
//' \item{\code{post_change_mean}}{is the vector of post-change means for the p-variate time series.}
//' \item{\code{nb_candidates}}{is a number of candidates at time n.}
//' \item{\code{candidates}}{is a list of candidates at each at time n.}
//' \item{\code{nb_candidates_at_time}}{is a list of candidatenumber at each iteration.}
//' \item{\code{opt_changes_at_time}}{is a list of optimal changes at each iteration.}
//' \item{\code{opt_costs_at_time}}{is a list of coptimal costs (without constant B') at each iteration.}
//' }
//' @examples
//' 
//' set.seed(21)
//' N <- 1000
//' P <- 2
//' Change <- N %/% 2
//' theta0 <-  rep(0, P)
//' #Data generation
//' ##the time series with one change
//' ts_gauss <- generate_ts(type = "gauss", p = P, n = N, changes = Change, means = matrix(c(theta0, theta0 + 5), nrow = P))
//' ts_poisson <- generate_ts(type = "poisson", p = P, n = N, changes = Change,  means = matrix(c(theta0 + 1, theta0 + 5), nrow = P))
//' ##the time series with one change
//' ts_gauss0 <- generate_ts(type = "gauss", p = P, n = N, changes = NULL, means = matrix(0, ncol = 1, nrow = P))
//' ts_poisson0 <- generate_ts(type = "poisson", p = P, n = N, changes = NULL, means = matrix(1, ncol = 1, nrow = P))
//' dyadMdFOCuS(ts_gauss, method ='FOCuS0', qmin = 2)

// [[Rcpp::export]]
List dyadMdFOCuS(Rcpp::NumericMatrix data, 
                     std::string method = "FOCuS0", 
                     std::string cost = "gauss",
                     bool cand_nb = false, 
                     bool opt_changes = false, 
                     bool opt_costs = false,
                     bool cands = false,
                     unsigned int qmin = 3
                       ) {
  //data parameters
  unsigned int p = (unsigned int)data.ncol();
  unsigned int length = (unsigned int) data.nrow();
  //stop
  if (method != "FOCuS0" && method !="FOCuS") {throw std::range_error("Parameter method should be FOCuS0 or FOCuS !");}
  if (cost != "gauss" && cost != "poisson") {throw std::range_error("Parameter cost should be gauss or poisson !");}
  if (qmin < log2(p + 2)) {throw std::range_error("The integer parameter qmin must be greater than or equal to log2(p+2)!");}
  
  // cost and method
  unsigned int typeCost = get_costType(method, cost);
  
  //memory help tables and vectors
  bool* isPruning = new bool[length - 1];   //for pruning's flag
  double** cumsumMatrix = new double*[length];
  for (unsigned int i = 0; i < length; i++) cumsumMatrix[i] = new double[p + 1];
  
  //pre-processing
  for (int i  = 0;  i < (length - 1); i++) isPruning[i] = false;
  get_cumsumMatrix(p, length, cumsumMatrix, data);

  //initialization
  std::list<int> candidates;     // change point candidates at time t
  std::vector<unsigned int> nb_at_time;   // vector of candidate number at time t
  std::vector<double>  opt_cost;          // optimal cost at time t
  std::vector<int>  opt_change;            // optimal change at time t
  std::vector<unsigned int> candidate_set; // candidates at time t=length-1
  for (unsigned int i = 0; i < length; i++) {
    nb_at_time.push_back(0);
    opt_cost.push_back(INFINITY);
    opt_change.push_back(0);
  }
  opt_cost[0] = 0;
  double candidate_cost;
  //loop
  for (unsigned int i = 1; i < length; i++) {
    candidates.push_back(i-1);
    nb_at_time[i] = candidates.size();
    // optimal cost and change
    for(auto it = candidates.begin(); it != candidates.end(); it++){
      candidate_cost = get_cost(typeCost, p, (*it), i, cumsumMatrix);
      if (candidate_cost < opt_cost[i]) {                      //find min
        opt_cost[i] = candidate_cost;
        opt_change[i] = (*it) + 1; //add 1 for R
      }
    }
    //pruning
    doPruning_dyadic(qmin, i, p, candidates, isPruning, cumsumMatrix);  // i not (i-1) because counting starts from 1 and not from 0
  }
  //get res
  unsigned int change  = opt_change[length-1];
  std::vector<double> pre_change_mean;
  std::vector<double> post_change_mean;
  for(unsigned int k = 1; k <= p; k++) {
    pre_change_mean.push_back(cumsumMatrix[change-1][k]/change); 
    post_change_mean.push_back((cumsumMatrix[length-1][k] - cumsumMatrix[change-1][k])/(length - change)); 
  }
  //add the last nb_cand
  nb_at_time.push_back(candidates.size());
  //memory
  for (unsigned int i = 0; i < length; i++) delete(cumsumMatrix[i]);
  delete [] isPruning;
  delete [] cumsumMatrix ;
  isPruning = NULL;
  cumsumMatrix = NULL;
  //constant B'(x) = sum_{i=1}^length x²
  double constant = 0;
  for (size_t i = 0; i < length; i++)
    for (size_t k = 0; k < p; k++) 
      constant = constant + data(i, k)*data(i, k);
  List res;
  res["change"] = change;
  res["cost"] = opt_cost[length-1] + constant;
  res["pre_change_mean"] = pre_change_mean;
  res["post_change_mean"] = post_change_mean;
  
  res["nb_candidates"] = nb_at_time[length];
  
  if (cands){
    std::vector<int> cands_n;
    auto pruning_it = candidates.begin();
    while (pruning_it != candidates.end()) {
      cands_n.push_back((*pruning_it)+1);//for R
      ++pruning_it;
    }
    res["candidates"] = cands_n;
  } 
  if (cand_nb) res["nb_candidates_at_time"] = nb_at_time;
  if (opt_changes) res["opt_changes_at_time"] = opt_change;
  if (opt_costs) res["opt_costs_at_time"] = opt_cost;
  return res;
}
