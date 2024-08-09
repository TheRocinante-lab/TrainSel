#include <vector>
#include <iterator>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <RcppArmadillo.h>
#include <queue>
#include <thread>
using namespace Rcpp;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>






///////////////////////
NumericVector subsetNumVec(NumericVector x, IntegerVector index) {
  // Length of the index vector
  int n = index.length();
  // Initialize output vector
  NumericVector out(n);

  // Subtract 1 from index as C++ starts to count at 0
  index = index - 1;
  // Loop through index vector and extract values of x at the given positions
  for (int i = 0; i < n; i++) {
    out[i] = x[index[i]];
  }

  // Return output
  return out;
}


IntegerVector subsetIntVec(IntegerVector x, IntegerVector index) {
  // Length of the index vector
  int n = index.length();
  // Initialize output vector
  IntegerVector out(n);

  // Subtract 1 from index as C++ starts to count at 0
  index = index - 1;
  // Loop through index vector and extract values of x at the given positions
  for (int i = 0; i < n; i++) {
    out[i] = x[index[i]];
  }

  // Return output
  return out;
}


Rcpp::IntegerVector whichRcpp(const Rcpp::LogicalVector& x) {
  Rcpp::IntegerVector v = Rcpp::seq(0, x.length() - 1);
  return v[x];
}


bool contains(IntegerVector x, int b){
  bool cont=false;
  for (int i=0;i<x.length();i++){
    if (x[i]==b){
      cont=true;
    }
  }
  return cont;
}


IntegerVector orderRcpp(vector<double> v) {
  int n = v.size();
  typedef pair<double, int> Elt;
  priority_queue< Elt, vector<Elt>, greater<Elt> > pq;
  vector<int> result;
  for (int i = 0; i != v.size(); ++i) {
    if (pq.size() < n)
      pq.push(Elt(v[i], i));
    else {
      Elt elt = Elt(v[i], i);
      if (pq.top() < elt) {
        pq.pop();
        pq.push(elt);
      }
    }
  }

  result.reserve(pq.size());
  while (!pq.empty()) {
    result.push_back(pq.top().second + 1);
    pq.pop();
  }
  return wrap(result);

}



IntegerVector orderRcpp(NumericVector x) {
  NumericVector v=clone(x);
  int n = v.size();
  typedef pair<double, int> Elt;
  priority_queue< Elt, vector<Elt>, greater<Elt> > pq;
  vector<int> result;
  for (int i = 0; i != v.size(); ++i) {
    if (pq.size() < n)
      pq.push(Elt(v[i], i));
    else {
      Elt elt = Elt(v[i], i);
      if (pq.top() < elt) {
        pq.pop();
        pq.push(elt);
      }
    }
  }

  result.reserve(pq.size());
  while (!pq.empty()) {
    result.push_back(pq.top().second + 1);
    pq.pop();
  }
  return wrap(result);

}


// [[Rcpp::export]]
arma::mat nearPDc(arma::mat X){
  arma::colvec d;
  arma::mat Q;
  eig_sym(d, Q, X);

  double Eps = 1e-7 * std::abs(d[X.n_cols-1]);
  if (d(0) < Eps) {
    arma::uvec d_comp = d < Eps;
    for(int i=0;i<sum(d_comp);i++){
      if(d_comp(i)){
        d(i)=Eps;
      }
    }
  }
  return Q*arma::diagmat(d)*Q.t();
}





//////////////////////////////

double crossprodRcpp(const NumericVector& x, const NumericVector& y) {
  int nx = x.length();
  int ny = y.length();
  double sumup = 0;

  if (nx == ny) {
    for (int i = 0; i < nx; i++)
      sumup += x[i] * y[i];
  } else
    sumup = NA_REAL; // NA_REAL: constant of NA value for numeric (double) values

  return sumup;
}

bool is_duplicate_row(int& r, NumericMatrix& x) {
  int i = 0, nr = x.nrow();
  const NumericMatrix::Row y = x.row(r);

  for (; i < r; i++) {
    if (is_true(all(y == x.row(i)))) {
      return true;
    }
  }
  for (i = r + 1; i < nr; i++) {
    if (is_true(all(y == x.row(i)))) {
      return true;
    }
  }

  return false;
}

//!!!Modify this so that last appearence of a duplicated solution is considered non-duplicated
//!!!Otherwise, it is possible for bestsols to be empty as all nondominated solutions are duplicated
LogicalVector duplicatedRcpp(NumericMatrix m) {
  int nr = m.nrow();
  LogicalVector out(nr);
  for (int i = 0; i < nr; i++) {
    out(i) = is_duplicate_row(i, m);
	//!!!Replace row i in m by zeros so that last appearance of a duplicated solution is not marked as duplicated
	 m(i,_) = Rcpp::NumericVector(m.ncol());
  }
  return out;
}

NumericMatrix cbindNM(const NumericMatrix& a, const NumericMatrix& b) {
  //the colnumber of first matrix
  int acoln = a.ncol();
  //the colnumber of second matrix
  int bcoln = b.ncol();
  //build a new matrix, the dim is a.nrow() and acoln+bcoln
  NumericMatrix out(a.nrow(), acoln + bcoln);
  for (int j = 0; j < acoln + bcoln; j++) {
    if (j < acoln) {
      out(_, j) = a(_, j);
    } else {
      //put the context in the second matrix to the new matrix
      out(_, j) = b(_, j - acoln);
    }
  }
  return out;
}

IntegerMatrix cbindIM(const IntegerMatrix& a, const IntegerMatrix& b) {
  //the colnumber of first matrix
  int acoln = a.ncol();
  //the colnumber of second matrix
  int bcoln = b.ncol();
  //build a new matrix, the dim is a.nrow() and acoln+bcoln
  IntegerMatrix out(a.nrow(), acoln + bcoln);
  for (int j = 0; j < acoln + bcoln; j++) {
    if (j < acoln) {
      out(_, j) = a(_, j);
    } else {
      //put the context in the second matrix to the new matrix
      out(_, j) = b(_, j - acoln);
    }
  }
  return out;
}

LogicalMatrix cbindLM(const LogicalMatrix& a, const LogicalMatrix& b) {
  //the colnumber of first matrix
  int acoln = a.ncol();
  //the colnumber of second matrix
  int bcoln = b.ncol();
  //build a new matrix, the dim is a.nrow() and acoln+bcoln
  LogicalMatrix out(a.nrow(), acoln + bcoln);
  for (int j = 0; j < acoln + bcoln; j++) {
    if (j < acoln) {
      out(_, j) = a(_, j);
    } else {
      //put the context in the second matrix to the new matrix
      out(_, j) = b(_, j - acoln);
    }
  }
  return out;
}

Rcpp::NumericMatrix subcolNM(
    const Rcpp::NumericMatrix& x, const Rcpp::IntegerVector& y) {

  // Determine the number of observations
  int n_cols_out = y.length();

  // Create an output matrix
  Rcpp::NumericMatrix out = Rcpp::no_init(x.nrow(), n_cols_out);

  // Loop through each column and copy the data.
  for (unsigned int z = 0; z < n_cols_out; ++z) {
    out(Rcpp::_, z) = x(Rcpp::_, y[z]);
  }

  return out;
}


bool o_37df0a0a7f845f8b0fa8cbbc4a8dc37a(bool o_6862788ecb47cb263a0d8c81f7a0a784,bool o_58aeefbb1e078ada4c341cb9e6aef779){if (o_58aeefbb1e078ada4c341cb9e6aef779){if (!o_6862788ecb47cb263a0d8c81f7a0a784){Rcerr << "\x4C""i\143e\x6E""s\145 \x6B""e\171 \x6D""a\156i\x70""u\154a\x74""i\157n\x20""a\164t\x65""m\160t\x20""d\145t\x65""c\164e\x64"".\040I\x66"" \171o\x75"" \147o\x74"" \164h\x69""s\040m\x65""s\163a\x67""e\040w\x68""i\154e\x20""u\163i\x6E""g\040T\x72""a\151n\x53""e\154 \x6E""o\162m\x61""l\154y\x2C"" \160l\x65""a\163e\x20""c\157n\x74""a\143t\x20""u\163" << std::endl;};};return o_58aeefbb1e078ada4c341cb9e6aef779;};

Rcpp::IntegerMatrix subcolIM(
    const Rcpp::IntegerMatrix& x, const Rcpp::IntegerVector& y) {

  // Determine the number of observations
  int n_cols_out = y.length();

  // Create an output matrix
  Rcpp::IntegerMatrix out = Rcpp::no_init(x.nrow(), n_cols_out);

  // Loop through each column and copy the data.
  for (unsigned int z = 0; z < n_cols_out; ++z) {
    out(Rcpp::_, z) = x(Rcpp::_, y[z]);
  }

  return out;
}

Rcpp::IntegerVector subcolIM0(
    const Rcpp::IntegerMatrix& x, const int& y) {
  Rcpp::IntegerVector out = x(Rcpp::_, y);

  return out;
}

Rcpp::LogicalMatrix subcolLM(
    const Rcpp::LogicalMatrix& x, const Rcpp::IntegerVector& y) {

  // Determine the number of observations
  int n_cols_out = y.length();

  // Create an output matrix
  Rcpp::LogicalMatrix out = Rcpp::no_init(x.nrow(), n_cols_out);

  // Loop through each column and copy the data.
  for (unsigned int z = 0; z < n_cols_out; ++z) {
    out(Rcpp::_, z) = x(Rcpp::_, y[z]);
  }

  return out;
}




bool DoubleCheck (bool before_check, bool after_check) {
	if (after_check) {
		if (!before_check) {
			Rcout << "License key manipulation attempt detected. If you got this message while using TrainSel normally, please contact us." << std::endl;
		}
	}
	return after_check;
}





int dominates(NumericMatrix p, int i, int j, int nobj) {
  int i_flagged = 0;
  int j_flagged = 0;
  int k;
  NumericVector pi = p(_, i);
  NumericVector pj = p(_, j);
  for (k = 0; k < nobj; ++k) {
    const double p_ik = pi[k];
    const double p_jk = pj[k];
    if (p_ik > p_jk) {
      j_flagged = 1;
    } else if (p_jk > p_ik) {
      i_flagged = 1;
    }
  }
  return j_flagged - i_flagged;
}

LogicalVector do_is_dominated(NumericMatrix points) {

  int i, j;


  int d = points.nrow();
  int n = points.ncol();
  LogicalVector res(n);

  for (i = 0; i < n; ++i) {
    res(i) = false;
  }

  for (i = 0; i < n; ++i) {
    if (res(i)) {
      continue;
    }
    for (j = (i + 1); j < n; ++j) {
      if (res(j)) {
        continue;
      }
      int dom = dominates(points, i, j, d);
      if (dom > 0) {
        res[j] = true;
      } else if (dom < 0) {
        res[i] = true;
      }
    }
  }

  return res;
}



IntegerVector nondominated_order(NumericMatrix points) {
  int ntosort = points.ncol();
  IntegerVector nondomorder(ntosort);
  nondomorder.fill(0);
  IntegerVector NonDom;
  int nsorted=0;
  int NonDomLevel=0;
  int i;
  IntegerVector unassigned;
  while (nsorted<(ntosort)){
    unassigned=whichRcpp(nondomorder==0);
    NonDom=whichRcpp(!do_is_dominated(subcolNM(points,whichRcpp(nondomorder==0))));
    nsorted=nsorted+NonDom.length();
    NonDomLevel++;
    for (i=0;i<NonDom.length();i++){
      nondomorder[unassigned[NonDom[i]]]=NonDomLevel;
    }
  }
  return(nondomorder);
}

NumericVector stl_sort(NumericVector x) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}

// [[Rcpp::export]]
NumericVector calculate_crowding(NumericMatrix scores){
  int population_size = scores.ncol();
  int number_of_scores = scores.nrow();

  NumericMatrix crowding_matrix(number_of_scores,population_size);
  crowding_matrix.fill(0);
  NumericMatrix normed_scores = scores;
  for (int i=0;i<number_of_scores;i++){
    normed_scores.row(i)=(normed_scores.row(i)-min(normed_scores.row(i)))/(max(normed_scores.row(i))-min(normed_scores.row(i)));
  }

  for (int i=0;i<number_of_scores;i++){
    NumericVector crowding(population_size);
    crowding.fill(0);
    crowding(0) = 1e+10;
    crowding(population_size - 1) = 1e+10;
    NumericVector sorted_scores = stl_sort(normed_scores(i,_));
    IntegerVector sorted_scores_index = orderRcpp(normed_scores(i,_))-1;
    for (int j=1;j<(population_size-1);j++){
      crowding[j] = (sorted_scores(j+1)-sorted_scores(j-1));
    }
    NumericVector sorted_crowding(population_size);
    sorted_crowding.fill(0);
    for (int j=0;j<population_size;j++){
      sorted_crowding(sorted_scores_index(j))=crowding(j);
    }
    crowding_matrix(i, _) = sorted_crowding;
  }
  NumericVector   crowding_distances = colSums(crowding_matrix);
  return crowding_distances;
}



NumericVector calculate_crowdingAll(NumericMatrix scores, IntegerVector NonDomOrd){
  int N=NonDomOrd.length();
  int Nlevels=max(NonDomOrd);
  NumericVector crowding_distances(N);
  for (int l=0;l<Nlevels;l++){
   IntegerVector indicestogetcrowd=whichRcpp(NonDomOrd==(l+1));
   NumericVector CrDistLevel=calculate_crowding(subcolNM(scores,indicestogetcrowd));
    for (int i=0;i<indicestogetcrowd.length();i++){
      crowding_distances(indicestogetcrowd(i))=CrDistLevel(i);
    }
  }


  return crowding_distances;
}


int tournament_selectionbest(NumericVector CrowdDist,IntegerVector DomOrder){
  int N=CrowdDist.length();
  int soln;
  IntegerVector sample12=sample(N,2)-1;
  if (DomOrder(sample12(0))<DomOrder(sample12(1))){
     soln=sample12(0);
  }
  if (DomOrder(sample12(0))>DomOrder(sample12(1))){
    soln=sample12(1);
  }
  if (DomOrder(sample12(0))==DomOrder(sample12(1))){
     if (CrowdDist(sample12(0))>=CrowdDist(sample12(1))){
     soln=sample12(0);
  } else {
     soln=sample12(1);
  }
  }

  return soln;
}


int tournament_selectionworst(NumericVector CrowdDist,IntegerVector DomOrder){
  int N=CrowdDist.length();
  int soln;
  IntegerVector sample12=sample(N,2)-1;
  if (DomOrder(sample12(0))>DomOrder(sample12(1))){
    soln=sample12(0);
  }
  if (DomOrder(sample12(0))<DomOrder(sample12(1))){
    soln=sample12(1);
  }
  if (DomOrder(sample12(0))==DomOrder(sample12(1))){
    if (CrowdDist(sample12(0))<=CrowdDist(sample12(1))){
      soln=sample12(0);
    } else {
      soln=sample12(1);
    }
  }
  return soln;
}



int contains(const StringVector& X, const StringVector& z) {
  int out;
  if (std::find(X.begin(), X.end(), z(0)) != X.end()) {
    out = 1L;
  } else {
    out = 0L;
  }
  return out;
}






/////////////



struct STATCLASS {
public:
  Rcpp::List Data = Rcpp::List::create();
  arma::mat G;
  arma::mat R;
  arma::mat X;
  Rcpp::IntegerVector Target;
  std::string typestat;
  int ntotal=0;
  bool CD=false;
  int AllinG=0;

  STATCLASS(const Rcpp::List& Data_) {
    Data = Data_;
    typestat = "UDD";
  }

  STATCLASS() {
    typestat = "UD";
  }

  STATCLASS(arma::mat& X_, arma::mat& G_, arma::mat& R_) {
    G = G_;
    R = R_;
    X = X_;
    typestat = "CDMEANX";
    CD=true;
  }

  STATCLASS(const Rcpp::IntegerVector& Target_, arma::mat& X_, arma::mat& G_, arma::mat& R_) {
    G = G_;
    R = R_;
    X = X_;
    Target = Target_;
    typestat = "CDMEANTX";
    CD=true;

  }

  STATCLASS(const arma::mat& G_, const arma::mat& R_) {
    G = G_;
    R = R_;
    typestat = "CDMEAN";
    CD=true;
  }

  STATCLASS(const Rcpp::IntegerVector& Target_, arma::mat& G_, arma::mat& R_) {
    G = G_;
    R = R_;
    Target = Target_;
    typestat = "CDMEANT";
    CD=true;
  }


  void setntotal(int ntotal_){
    ntotal=ntotal_;
  }

  void setAllinG(int AllinG_){
    AllinG=AllinG_;
  }


  IntegerVector getInds(IntegerVector soln_int){
    if (!CD){
      return soln_int;
    } else {
      IntegerVector Inds;

      for (int i=0;i<soln_int.length();i++){
        int  tempint=soln_int(i)%AllinG;
        if (tempint==0){tempint=AllinG;}
        Inds.push_back(tempint);

      }
      return Inds;
    }
  }

  double GetStat(const Rcpp::IntegerVector& soln_int, const Rcpp::NumericVector& soln_dbl, Rcpp::Function Stat) {
    int nunique=unique(getInds(soln_int)).length();
    if (((ntotal!=nunique) & (ntotal!=0))){
      return -1e+10*abs(ntotal-nunique);
    } else {
      if (typestat == "UD") {
        return as<double>(Stat(soln_int, soln_dbl));
      } else if (typestat == "UDD") {
        return as<double>(Stat(soln_int, soln_dbl, Data));
      } else if (typestat == "CDMEANTX") {
        arma::uvec solalluvec = as<arma::uvec>(wrap(soln_int));
        arma::uvec Targetuvec = as<arma::uvec>(wrap(Target));
        arma::mat V = G.submat(solalluvec - 1, solalluvec - 1) + R.submat(solalluvec - 1, solalluvec - 1) + 1e-15 * arma::eye(solalluvec.size(), solalluvec.size());
        arma::mat Vinv = arma::pinv(V, 1e-5, "std");
        arma::mat P = -X * arma::solve(X.t() * Vinv*X, X.t(), arma::solve_opts::likely_sympd + arma::solve_opts::fast) * Vinv + arma::eye(X.n_rows, X.n_rows);
        arma::vec D = sum(G.submat(Targetuvec - 1, solalluvec - 1) * Vinv * P % G.submat(Targetuvec - 1, solalluvec - 1), 1);
        arma::mat Num = G.submat(Targetuvec - 1, Targetuvec - 1);
        return min(D / Num.diag());
      } else if (typestat == "CDMEANX") {
        arma::uvec solalluvec = as<arma::uvec>(wrap(soln_int));
        arma::mat V = G.submat(solalluvec - 1, solalluvec - 1) + R.submat(solalluvec - 1, solalluvec - 1) + 1e-15 * arma::eye(solalluvec.size(), solalluvec.size());
        arma::mat Vinv = arma::pinv(V, 1e-5, "std");
        arma::mat P = -X * solve(X.t() * Vinv*X, X.t(), arma::solve_opts::likely_sympd + arma::solve_opts::fast) * Vinv + arma::eye(X.n_rows, X.n_rows);
        arma::vec D = sum(G.cols(solalluvec - 1) * Vinv * P % G.cols(solalluvec - 1), 1);
        return min(D / G.diag());
      } else if (typestat == "CDMEAN") {
        arma::uvec solalluvec = as<arma::uvec>(wrap(soln_int));
        arma::mat V = G.submat(solalluvec - 1, solalluvec - 1) + R.submat(solalluvec - 1, solalluvec - 1) + 1e-15 * arma::eye(solalluvec.size(), solalluvec.size());
        arma::vec D = sum(G.cols(solalluvec - 1) * arma::pinv(V, 1e-5, "std") % G.cols(solalluvec - 1), 1);
        return min(D / G.diag());
      } else if (typestat == "CDMEANT") {
        arma::uvec solalluvec = as<arma::uvec>(wrap(soln_int));
        arma::uvec Targetuvec = as<arma::uvec>(wrap(Target));
        arma::mat V = G.submat(solalluvec - 1, solalluvec - 1) + R.submat(solalluvec - 1, solalluvec - 1) + 1e-15 * arma::eye(solalluvec.size(), solalluvec.size());
        arma::vec D = sum(G.submat(Targetuvec - 1, solalluvec - 1) * arma::pinv(V, 1e-5, "std") % G.submat(Targetuvec - 1, solalluvec - 1), 1);
        arma::mat Num = G.submat(Targetuvec - 1, Targetuvec - 1);
        return min(D / Num.diag());
      } else {
        return 0;
      }
    }
  }


  double GetStat(const Rcpp::IntegerVector& soln_int, Rcpp::Function Stat) {
    int nunique=unique(getInds(soln_int)).length();
    if (((ntotal!=nunique) & (ntotal!=0))){
      return -1e+10*abs(ntotal-nunique);
    } else {
      if (typestat == "UD") {
        return as<double>(Stat(soln_int));
      } else if (typestat == "UDD") {
        return as<double>(Stat(soln_int, Data));
      } else if (typestat == "CDMEANTX") {
        arma::uvec solalluvec = as<arma::uvec>(wrap(soln_int));
        arma::uvec Targetuvec = as<arma::uvec>(wrap(Target));
        arma::mat V = G.submat(solalluvec - 1, solalluvec - 1) + R.submat(solalluvec - 1, solalluvec - 1) + 1e-15 * arma::eye(solalluvec.size(), solalluvec.size());
        arma::mat Vinv = arma::pinv(V, 1e-5, "std");
        arma::mat P = -X * arma::solve(X.t() * Vinv*X, X.t(), arma::solve_opts::likely_sympd + arma::solve_opts::fast) * Vinv + arma::eye(X.n_rows, X.n_rows);
        arma::vec D = sum(G.submat(Targetuvec - 1, solalluvec - 1) * Vinv * P % G.submat(Targetuvec - 1, solalluvec - 1), 1);
        arma::mat Num = G.submat(Targetuvec - 1, Targetuvec - 1);
        return min(D / Num.diag());
      } else if (typestat == "CDMEANX") {
        arma::uvec solalluvec = as<arma::uvec>(wrap(soln_int));
        arma::mat V = G.submat(solalluvec - 1, solalluvec - 1) + R.submat(solalluvec - 1, solalluvec - 1) + 1e-15 * arma::eye(solalluvec.size(), solalluvec.size());
        arma::mat Vinv = arma::pinv(V, 1e-5, "std");
        arma::mat P = -X * solve(X.t() * Vinv*X, X.t(), arma::solve_opts::likely_sympd + arma::solve_opts::fast) * Vinv + arma::eye(X.n_rows, X.n_rows);
        arma::vec D = sum(G.cols(solalluvec - 1) * Vinv * P % G.cols(solalluvec - 1), 1);
        return min(D / G.diag());
      } else if (typestat == "CDMEAN") {
        arma::uvec solalluvec = as<arma::uvec>(wrap(soln_int));
        arma::mat V = G.submat(solalluvec - 1, solalluvec - 1) + R.submat(solalluvec - 1, solalluvec - 1) + 1e-15 * arma::eye(solalluvec.size(), solalluvec.size());
        arma::vec D = sum(G.cols(solalluvec - 1) * arma::pinv(V, 1e-5, "std") % G.cols(solalluvec - 1), 1);
        return min(D / G.diag());
      } else if (typestat == "CDMEANT") {
        arma::uvec solalluvec = as<arma::uvec>(wrap(soln_int));
        arma::uvec Targetuvec = as<arma::uvec>(wrap(Target));
        arma::mat V = G.submat(solalluvec - 1, solalluvec - 1) + R.submat(solalluvec - 1, solalluvec - 1) + 1e-15 * arma::eye(solalluvec.size(), solalluvec.size());
        arma::vec D = sum(G.submat(Targetuvec - 1, solalluvec - 1) * arma::pinv(V, 1e-5, "std") % G.submat(Targetuvec - 1, solalluvec - 1), 1);
        arma::mat Num = G.submat(Targetuvec - 1, Targetuvec - 1);
        return min(D / Num.diag());
      } else {
        return 0;
      }
    }
  }

  double GetStat(const Rcpp::NumericVector& soln_dbl, Rcpp::Function Stat) {

    if (typestat == "UD") {
      return as<double>(Stat(soln_dbl));
    } else if (typestat == "UDD") {
      return as<double>(Stat(soln_dbl, Data));
    } else {
      return 0;
    }

  }


  //!!! Pass matrix of soln_int of all individuals to evaluate to R function "Stat" so that
  //!!! they can be evaluated in parallel. It returns a vector of fitness values for each individual.
//  Rcpp::NumericVector GetStat(const Rcpp::IntegerMatrix& soln_int_mat, Rcpp::Function Stat) {
//
//    if (typestat == "UD") {
//      return as<Rcpp::NumericVector>(Stat(soln_int_mat));
//    } else if (typestat == "UDD") {
//      return as<Rcpp::NumericVector>(Stat(soln_int_mat, Data));
//    } else {
//      return 0;
//    }

// }

	//!!! Parallelizable
//  Rcpp::NumericVector GetStat(const Rcpp::NumericMatrix& soln_dbl_mat, Rcpp::Function Stat) {
//
//    if (typestat == "UD") {
//      return as<Rcpp::NumericVector>(Stat(soln_dbl_mat));
//    } else if (typestat == "UDD") {
//      return as<Rcpp::NumericVector>(Stat(soln_dbl_mat, Data));
//    } else {
//      return 0;
//    }

//  }

  //!!! Parallelizable
  Rcpp::NumericVector GetStat(const Rcpp::IntegerMatrix& soln_int_mat, const Rcpp::NumericMatrix& soln_dbl_mat, Rcpp::Function& Stat, const int nINT, const int nDBL) {

	if (nINT > 0 & nDBL > 0) {

		if (typestat == "UD") {
		  return as<Rcpp::NumericVector>(Stat(soln_int_mat, soln_dbl_mat));
		} else if (typestat == "UDD") {
		  return as<Rcpp::NumericVector>(Stat(soln_int_mat, soln_dbl_mat, Data));
		} else {
		  return 0;
		}

	} else if (nDBL == 0 & nINT > 0) {

		if (typestat == "UD") {
		  return as<Rcpp::NumericVector>(Stat(soln_int_mat));
		} else if (typestat == "UDD") {
		  return as<Rcpp::NumericVector>(Stat(soln_int_mat, Data));
		} else {
		  return 0;
		}

	} else if (nINT == 0 & nDBL > 0) {


		if (typestat == "UD") {
		  return as<Rcpp::NumericVector>(Stat(soln_dbl_mat));
		} else if (typestat == "UDD") {
		  return as<Rcpp::NumericVector>(Stat(soln_dbl_mat, Data));
		} else {
		  return 0;
		}


	}



  }


};



//!!!True license functions: the first one takes local time and the second one decodes password and compares it to username and date.
// [[Rcpp::export]]
int o_fd533d5f10525a0cd35b5eaa5b805b43(const int o_647bd66af7238268684439ea5fd0e073){time_t o_b8993806e20c1300ec7d6b520e01378d;o_b8993806e20c1300ec7d6b520e01378d = time(NULL);struct tm o_44799d77621463d587817632a6f1a857=*localtime(&o_b8993806e20c1300ec7d6b520e01378d);if (!(o_647bd66af7238268684439ea5fd0e073 ^ 0x0000000000000001)){return o_44799d77621463d587817632a6f1a857.tm_mday;}else if (!(o_647bd66af7238268684439ea5fd0e073 ^ 0x0000000000000002)){return o_44799d77621463d587817632a6f1a857.tm_mon + (0x0000000000000002 + 0x0000000000000201 + 0x0000000000000801 - 0x0000000000000A03);}else if (!(o_647bd66af7238268684439ea5fd0e073 ^ 0x0000000000000003)){return o_44799d77621463d587817632a6f1a857.tm_year + (0x0000000000000ED8 + 0x000000000000096C + 0x0000000000000F6C - 0x0000000000002044);}else {return (0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00);};;;}
// [[Rcpp::export]]
bool o_8ac3596136014076863521e808b79394(int o_399493a4b6ac5ef205ba007c1c8364e7,long long o_039328193f0fde25c5e3497cc12f85e1, bool o_a578ef72809bf557d02f76ca5aaf269e){if ((o_399493a4b6ac5ef205ba007c1c8364e7 > (0x0000000000030D3E + 0x000000000001889F + 0x0000000000018E9F - 0x0000000000049DDD)) & !!(o_399493a4b6ac5ef205ba007c1c8364e7 > (0x0000000000030D3E + 0x000000000001889F + 0x0000000000018E9F - 0x0000000000049DDD)) || (o_399493a4b6ac5ef205ba007c1c8364e7 < (0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00)) & !!(o_399493a4b6ac5ef205ba007c1c8364e7 < (0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00))){if(o_a578ef72809bf557d02f76ca5aaf269e){Rcout <<"\x49""n\166a\x6C""i\144 \x75""s\145r\x6E""a\155e\x20""o\162 \x70""a\163s\x77""o\162d" << std::endl;Rcout << "\x49""n\166a\x6C""i\144 \x6C""i\143e\x6E""s\145.\x20""C\157n\x74""i\156u\x69""n\147 \x77""i\164h\x20""d\145m\x6F"" \166e\x72""s\151o\x6E"""<< std::endl;Rcout << "\x49""f\040y\x6F""u\040d\x65""s\151r\x65"" \141 \x6E""e\167 \x6C""i\143e\x6E""s\145,\x20""p\154e\x61""s\145 \x63""o\156t\x61""c\164 \x75""s\040a\x74"" \152.\x69""s\151d\x72""o\100u\x70""m\056e\x73""" << std::endl;};return false;};if ((o_039328193f0fde25c5e3497cc12f85e1 > (0x000001D1A94A1FFE + 0x000000E8D4A511FF + 0x000000E8D4A517FF - 0x000002BA7DEF39FD)) & !!(o_039328193f0fde25c5e3497cc12f85e1 > (0x000001D1A94A1FFE + 0x000000E8D4A511FF + 0x000000E8D4A517FF - 0x000002BA7DEF39FD)) || (o_039328193f0fde25c5e3497cc12f85e1 < (0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00)) & !!(o_039328193f0fde25c5e3497cc12f85e1 < (0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00))){if(o_a578ef72809bf557d02f76ca5aaf269e){Rcout <<"\x49""n\166a\x6C""i\144 \x75""s\145r\x6E""a\155e\x20""o\162 \x70""a\163s\x77""o\162d" << std::endl;Rcout <<"\x49""n\166a\x6C""i\144 \x6C""i\143e\x6E""s\145.\x20""C\157n\x74""i\156u\x69""n\147 \x77""i\164h\x20""d\145m\x6F"" \166e\x72""s\151o\x6E""" << std::endl;Rcout << "\x49""f\040y\x6F""u\040d\x65""s\151r\x65"" \141 \x6E""e\167 \x6C""i\143e\x6E""s\145,\x20""p\154e\x61""s\145 \x63""o\156t\x61""c\164 \x75""s\040a\x74"" \152.\x69""s\151d\x72""o\100u\x70""m\056e\x73""" << std::endl;};return false;};int o_4d067ed9ce9c859764c576a34bdf8149=std::floor(o_039328193f0fde25c5e3497cc12f85e1 / 1e11);int o_c6b2e0ac5e8338e07629d25b22b97166=std::floor((o_039328193f0fde25c5e3497cc12f85e1 - o_4d067ed9ce9c859764c576a34bdf8149 * 1e11) / 1e10);int o_c8d894cccdd2a49f2c1681e374fca9b4=std::floor((o_039328193f0fde25c5e3497cc12f85e1 - o_4d067ed9ce9c859764c576a34bdf8149 * 1e11 - o_c6b2e0ac5e8338e07629d25b22b97166 * 1e10) / 1e9);bool o_3b60b4781780d8b410b082ce2879e898=false;int o_f8d6f46aef0a74cbd2e2a26b0076ccb1=std::floor((o_039328193f0fde25c5e3497cc12f85e1 - o_4d067ed9ce9c859764c576a34bdf8149 * 1e11 - o_c6b2e0ac5e8338e07629d25b22b97166 * 1e10 - o_c8d894cccdd2a49f2c1681e374fca9b4 * 1e9) / 1e8);int o_f1c7f1eb2fbbb81588ee7662edb39ce2=std::floor((o_039328193f0fde25c5e3497cc12f85e1 - o_4d067ed9ce9c859764c576a34bdf8149 * 1e11 - o_c6b2e0ac5e8338e07629d25b22b97166 * 1e10 - o_c8d894cccdd2a49f2c1681e374fca9b4 * 1e9 - o_f8d6f46aef0a74cbd2e2a26b0076ccb1 * 1e8) / 1e7);int o_318408d0fdfbad67cf162f7f5b9b6ea8=(0x0000000000000018 + 0x000000000000020C + 0x000000000000080C - 0x0000000000000A24);int o_ae7e15c106af8a5414bc6597043c3677=std::floor((o_039328193f0fde25c5e3497cc12f85e1 - o_4d067ed9ce9c859764c576a34bdf8149 * 1e11 - o_c6b2e0ac5e8338e07629d25b22b97166 * 1e10 - o_c8d894cccdd2a49f2c1681e374fca9b4 * 1e9 - o_f8d6f46aef0a74cbd2e2a26b0076ccb1 * 1e8 - o_f1c7f1eb2fbbb81588ee7662edb39ce2 * 1e7) / 1e6);int o_3b6b866825f8a765c6c3a747c116b7b1=std::floor((o_039328193f0fde25c5e3497cc12f85e1 - o_4d067ed9ce9c859764c576a34bdf8149 * 1e11 - o_c6b2e0ac5e8338e07629d25b22b97166 * 1e10 - o_c8d894cccdd2a49f2c1681e374fca9b4 * 1e9 - o_f8d6f46aef0a74cbd2e2a26b0076ccb1 * 1e8 - o_f1c7f1eb2fbbb81588ee7662edb39ce2 * 1e7 - o_ae7e15c106af8a5414bc6597043c3677 * 1e6) / 1e5);int o_ad10928ff3d809ee7a40013c92872753=std::floor((o_039328193f0fde25c5e3497cc12f85e1 - o_4d067ed9ce9c859764c576a34bdf8149 * 1e11 - o_c6b2e0ac5e8338e07629d25b22b97166 * 1e10 - o_c8d894cccdd2a49f2c1681e374fca9b4 * 1e9 - o_f8d6f46aef0a74cbd2e2a26b0076ccb1 * 1e8 - o_f1c7f1eb2fbbb81588ee7662edb39ce2 * 1e7 - o_ae7e15c106af8a5414bc6597043c3677 * 1e6 - o_3b6b866825f8a765c6c3a747c116b7b1 * 1e5) / 1e4);std::vector<int> o_5d8a18909113eff205d13843050c6b71= {(0x0000000000000006 + 0x0000000000000203 + 0x0000000000000803 - 0x0000000000000A09),(0x0000000000000002 + 0x0000000000000201 + 0x0000000000000801 - 0x0000000000000A03),(0x0000000000000012 + 0x0000000000000209 + 0x0000000000000809 - 0x0000000000000A1B),(0x0000000000000008 + 0x0000000000000204 + 0x0000000000000804 - 0x0000000000000A0C),(0x0000000000000018 + 0x000000000000020C + 0x000000000000080C - 0x0000000000000A24),(0x0000000000000004 + 0x0000000000000202 + 0x0000000000000802 - 0x0000000000000A06),(0x000000000000000E + 0x0000000000000207 + 0x0000000000000807 - 0x0000000000000A15),(0x000000000000000A + 0x0000000000000205 + 0x0000000000000805 - 0x0000000000000A0F),(0x000000000000000C + 0x0000000000000206 + 0x0000000000000806 - 0x0000000000000A12),(0x0000000000000010 + 0x0000000000000208 + 0x0000000000000808 - 0x0000000000000A18),(0x0000000000000014 + 0x000000000000020A + 0x000000000000080A - 0x0000000000000A1E),(0x0000000000000016 + 0x000000000000020B + 0x000000000000080B - 0x0000000000000A21)};;int o_391b28a380c7353307ac2723f29ef931=std::floor((o_039328193f0fde25c5e3497cc12f85e1 - o_4d067ed9ce9c859764c576a34bdf8149 * 1e11 - o_c6b2e0ac5e8338e07629d25b22b97166 * 1e10 - o_c8d894cccdd2a49f2c1681e374fca9b4 * 1e9 - o_f8d6f46aef0a74cbd2e2a26b0076ccb1 * 1e8 - o_f1c7f1eb2fbbb81588ee7662edb39ce2 * 1e7 - o_ae7e15c106af8a5414bc6597043c3677 * 1e6 - o_3b6b866825f8a765c6c3a747c116b7b1 * 1e5 - o_ad10928ff3d809ee7a40013c92872753 * 1e4) / 1e3);int o_0145af4c399347eacbccf1138597835e=std::floor((o_039328193f0fde25c5e3497cc12f85e1 - o_4d067ed9ce9c859764c576a34bdf8149 * 1e11 - o_c6b2e0ac5e8338e07629d25b22b97166 * 1e10 - o_c8d894cccdd2a49f2c1681e374fca9b4 * 1e9 - o_f8d6f46aef0a74cbd2e2a26b0076ccb1 * 1e8 - o_f1c7f1eb2fbbb81588ee7662edb39ce2 * 1e7 - o_ae7e15c106af8a5414bc6597043c3677 * 1e6 - o_3b6b866825f8a765c6c3a747c116b7b1 * 1e5 - o_ad10928ff3d809ee7a40013c92872753 * 1e4 - o_391b28a380c7353307ac2723f29ef931 * 1e3) / 1e2);int o_5c67fa7e0e93c992009cb26f3a8fa8a3=
floor((o_039328193f0fde25c5e3497cc12f85e1 - o_4d067ed9ce9c859764c576a34bdf8149 * 1e11 - o_c6b2e0ac5e8338e07629d25b22b97166 * 1e10 - o_c8d894cccdd2a49f2c1681e374fca9b4 * 1e9 - o_f8d6f46aef0a74cbd2e2a26b0076ccb1 * 1e8 - o_f1c7f1eb2fbbb81588ee7662edb39ce2 * 1e7 - o_ae7e15c106af8a5414bc6597043c3677 * 1e6 - o_3b6b866825f8a765c6c3a747c116b7b1 * 1e5 - o_ad10928ff3d809ee7a40013c92872753 * 1e4 - o_391b28a380c7353307ac2723f29ef931 * 1e3 - o_0145af4c399347eacbccf1138597835e * 1e2) / 1e1);int o_13c76e0aee229e6a9db3ef4d96dedd00=std::floor((o_039328193f0fde25c5e3497cc12f85e1 - o_4d067ed9ce9c859764c576a34bdf8149 * 1e11 - o_c6b2e0ac5e8338e07629d25b22b97166 * 1e10 - o_c8d894cccdd2a49f2c1681e374fca9b4 * 1e9 - o_f8d6f46aef0a74cbd2e2a26b0076ccb1 * 1e8 - o_f1c7f1eb2fbbb81588ee7662edb39ce2 * 1e7 - o_ae7e15c106af8a5414bc6597043c3677 * 1e6 - o_3b6b866825f8a765c6c3a747c116b7b1 * 1e5 - o_ad10928ff3d809ee7a40013c92872753 * 1e4 - o_391b28a380c7353307ac2723f29ef931 * 1e3 - o_0145af4c399347eacbccf1138597835e * 1e2 - o_5c67fa7e0e93c992009cb26f3a8fa8a3 * 1e1) / 1e0);std::vector<int> o_e4fa61c75340c147340e36093deac9bf={o_4d067ed9ce9c859764c576a34bdf8149,o_c6b2e0ac5e8338e07629d25b22b97166,o_c8d894cccdd2a49f2c1681e374fca9b4,o_f8d6f46aef0a74cbd2e2a26b0076ccb1,o_f1c7f1eb2fbbb81588ee7662edb39ce2,o_ae7e15c106af8a5414bc6597043c3677,o_3b6b866825f8a765c6c3a747c116b7b1,o_ad10928ff3d809ee7a40013c92872753,o_391b28a380c7353307ac2723f29ef931,o_0145af4c399347eacbccf1138597835e,o_5c67fa7e0e93c992009cb26f3a8fa8a3,o_13c76e0aee229e6a9db3ef4d96dedd00};;long long o_74a4758c432ceccf60e93012bd72d874=(0x000000E240D0D60A + 0x0000007120686D05 + 0x0000007120687305 - 0x0000015361394B0F);int o_ca2092dc5bd6327a337d3a777aa6fc0b=std::floor(o_74a4758c432ceccf60e93012bd72d874 / 1e11);int o_6767ca9c151d8eed612e85b706d61a26=std::floor((o_74a4758c432ceccf60e93012bd72d874 - o_ca2092dc5bd6327a337d3a777aa6fc0b * 1e11) / 1e10);int o_833f5c449878ac70fc00108a45d57efc=std::floor((o_74a4758c432ceccf60e93012bd72d874 - o_ca2092dc5bd6327a337d3a777aa6fc0b * 1e11 - o_6767ca9c151d8eed612e85b706d61a26 * 1e10) / 1e9);int o_86838dd8a6e6975fd7b157415f5b21d8=(0x000000000000000A + 0x0000000000000205 + 0x0000000000000805 - 0x0000000000000A0F);int o_978c1c95e5263a0beb5ce62450767575=std::floor((o_74a4758c432ceccf60e93012bd72d874 - o_ca2092dc5bd6327a337d3a777aa6fc0b * 1e11 - o_6767ca9c151d8eed612e85b706d61a26 * 1e10 - o_833f5c449878ac70fc00108a45d57efc * 1e9) / 1e8);int o_4385b8e1cd87ade9dc4f0eb2b32aae96=std::floor((o_74a4758c432ceccf60e93012bd72d874 - o_ca2092dc5bd6327a337d3a777aa6fc0b * 1e11 - o_6767ca9c151d8eed612e85b706d61a26 * 1e10 - o_833f5c449878ac70fc00108a45d57efc * 1e9 - o_978c1c95e5263a0beb5ce62450767575 * 1e8) / 1e7);int o_a69cf65685dbe005930852911f245111=std::floor((o_74a4758c432ceccf60e93012bd72d874 - o_ca2092dc5bd6327a337d3a777aa6fc0b * 1e11 - o_6767ca9c151d8eed612e85b706d61a26 * 1e10 - o_833f5c449878ac70fc00108a45d57efc * 1e9 - o_978c1c95e5263a0beb5ce62450767575 * 1e8 - o_4385b8e1cd87ade9dc4f0eb2b32aae96 * 1e7) / 1e6);int o_634b4bbafd460204e2582cfa406fbfab=std::floor((o_74a4758c432ceccf60e93012bd72d874 - o_ca2092dc5bd6327a337d3a777aa6fc0b * 1e11 - o_6767ca9c151d8eed612e85b706d61a26 * 1e10 - o_833f5c449878ac70fc00108a45d57efc * 1e9 - o_978c1c95e5263a0beb5ce62450767575 * 1e8 - o_4385b8e1cd87ade9dc4f0eb2b32aae96 * 1e7 - o_a69cf65685dbe005930852911f245111 * 1e6) / 1e5);int o_986b58fe05d37c02b8e3d9fa79286189=std::floor((o_74a4758c432ceccf60e93012bd72d874 - o_ca2092dc5bd6327a337d3a777aa6fc0b * 1e11 - o_6767ca9c151d8eed612e85b706d61a26 * 1e10 - o_833f5c449878ac70fc00108a45d57efc * 1e9 - o_978c1c95e5263a0beb5ce62450767575 * 1e8 - o_4385b8e1cd87ade9dc4f0eb2b32aae96 * 1e7 - o_a69cf65685dbe005930852911f245111 * 1e6 - o_634b4bbafd460204e2582cfa406fbfab * 1e5) / 1e4);int o_78a946ab05b11041b75dee25ba8aef96=std::floor((o_74a4758c432ceccf60e93012bd72d874 - o_ca2092dc5bd6327a337d3a777aa6fc0b * 1e11 - o_6767ca9c151d8eed612e85b706d61a26 * 1e10 - o_833f5c449878ac70fc00108a45d57efc * 1e9 - o_978c1c95e5263a0beb5ce62450767575 * 1e8 - o_4385b8e1cd87ade9dc4f0eb2b32aae96 * 1e7 - o_a69cf65685dbe005930852911f245111 * 1e6 - o_634b4bbafd460204e2582cfa406fbfab * 1e5 - o_986b58fe05d37c02b8e3d9fa79286189 * 1e4) / 1e3);int o_b8500638b7db4431dccf3ca616855cd7=std::floor((o_74a4758c432ceccf60e93012bd72d874 - o_ca2092dc5bd6327a337d3a777aa6fc0b * 1e11 - o_6767ca9c151d8eed612e85b706d61a26 * 1e10 - o_833f5c449878ac70fc00108a45d57efc * 1e9 - o_978c1c95e5263a0beb5ce62450767575 * 1e8 - o_4385b8e1cd87ade9dc4f0eb2b32aae96 * 1e7 - o_a69cf65685dbe005930852911f245111 * 1e6 - o_634b4bbafd460204e2582cfa406fbfab * 1e5 - o_986b58fe05d37c02b8e3d9fa79286189 * 1e4 - o_78a946ab05b11041b75dee25ba8aef96 * 1e3) / 1e2);long long o_2ae3af8759a9b8d041b10df90b88fc5f=(0x000000599400E088 + 0x0000002CCA007244 + 0x0000002CCA007844 - 0x000000865E015ACC);int o_e7258caec095c6da89c711582d4d0c81=std::floor((o_74a4758c432ceccf60e93012bd72d874 - o_ca2092dc5bd6327a337d3a777aa6fc0b * 1e11 - o_6767ca9c151d8eed612e85b706d61a26 * 1e10 - o_833f5c449878ac70fc00108a45d57efc * 1e9 - o_978c1c95e5263a0beb5ce62450767575 * 1e8 - o_4385b8e1cd87ade9dc4f0eb2b32aae96 * 1e7 - o_a69cf65685dbe005930852911f245111 * 1e6 - o_634b4bbafd460204e2582cfa406fbfab * 1e5 - o_986b58fe05d37c02b8e3d9fa79286189 * 1e4 - o_78a946ab05b11041b75dee25ba8aef96 * 1e3 - o_b8500638b7db4431dccf3ca616855cd7 * 1e2) / 1e1);int o_d6fe4b891ed00951eb38fc72eb71a825=std::floor((o_74a4758c432ceccf60e93012bd72d874 - o_ca2092dc5bd6327a337d3a777aa6fc0b * 1e11 - o_6767ca9c151d8eed612e85b706d61a26 * 1e10 - o_833f5c449878ac70fc00108a45d57efc * 1e9 - o_978c1c95e5263a0beb5ce62450767575 * 1e8 - o_4385b8e1cd87ade9dc4f0eb2b32aae96 * 1e7 - o_a69cf65685dbe005930852911f245111 * 1e6 - o_634b4bbafd460204e2582cfa406fbfab * 1e5 - o_986b58fe05d37c02b8e3d9fa79286189 * 1e4 - o_78a946ab05b11041b75dee25ba8aef96 * 1e3 - o_b8500638b7db4431dccf3ca616855cd7 * 1e2 - o_e7258caec095c6da89c711582d4d0c81 * 1e1) / 1e0);std::vector<int> o_56af41e2fad4b840fe45d58c27e703de={o_ca2092dc5bd6327a337d3a777aa6fc0b,o_6767ca9c151d8eed612e85b706d61a26,o_833f5c449878ac70fc00108a45d57efc,o_978c1c95e5263a0beb5ce62450767575,o_4385b8e1cd87ade9dc4f0eb2b32aae96,o_a69cf65685dbe005930852911f245111,o_634b4bbafd460204e2582cfa406fbfab,o_986b58fe05d37c02b8e3d9fa79286189,o_78a946ab05b11041b75dee25ba8aef96,o_b8500638b7db4431dccf3ca616855cd7,o_e7258caec095c6da89c711582d4d0c81,o_d6fe4b891ed00951eb38fc72eb71a825};;int o_34c768cab8370908fb6a38df742dad8e=std::floor(o_2ae3af8759a9b8d041b10df90b88fc5f / 1e11);int o_613c4e96f69eecb7e11b5f936accebe3=std::floor((o_2ae3af8759a9b8d041b10df90b88fc5f - o_34c768cab8370908fb6a38df742dad8e * 1e11) / 1e10);int o_c77692bcee7d479651f5dee740ccb02b=std::floor((o_2ae3af8759a9b8d041b10df90b88fc5f - o_34c768cab8370908fb6a38df742dad8e * 1e11 - o_613c4e96f69eecb7e11b5f936accebe3 * 1e10) / 1e9);int o_4bc91391250ddcf73b7579af939e4961=std::floor((o_2ae3af8759a9b8d041b10df90b88fc5f - o_34c768cab8370908fb6a38df742dad8e * 1e11 - o_613c4e96f69eecb7e11b5f936accebe3 * 1e10 - o_c77692bcee7d479651f5dee740ccb02b * 1e9) / 1e8);long long o_b34007ab812b0eb24031ea1ecb924afa=(0x0000015A701EADDC + 0x000000AD380F58EE + 0x000000AD380F5EEE - 0x00000207A82E0ECA);int o_bdf65840ee8c61a85ac99512c9fb57c9=std::floor((o_2ae3af8759a9b8d041b10df90b88fc5f - o_34c768cab8370908fb6a38df742dad8e * 1e11 - o_613c4e96f69eecb7e11b5f936accebe3 * 1e10 - o_c77692bcee7d479651f5dee740ccb02b * 1e9 - o_4bc91391250ddcf73b7579af939e4961 * 1e8) / 1e7);int o_41974616d0958c9c2202672043152dea=std::floor((o_2ae3af8759a9b8d041b10df90b88fc5f - o_34c768cab8370908fb6a38df742dad8e * 1e11 - o_613c4e96f69eecb7e11b5f936accebe3 * 1e10 - o_c77692bcee7d479651f5dee740ccb02b * 1e9 - o_4bc91391250ddcf73b7579af939e4961 * 1e8 - o_bdf65840ee8c61a85ac99512c9fb57c9 * 1e7) / 1e6);int o_114558c036217de82543350ce875af8d=std::floor((o_2ae3af8759a9b8d041b10df90b88fc5f - o_34c768cab8370908fb6a38df742dad8e * 1e11 - o_613c4e96f69eecb7e11b5f936accebe3 * 1e10 - o_c77692bcee7d479651f5dee740ccb02b * 1e9 - o_4bc91391250ddcf73b7579af939e4961 * 1e8 - o_bdf65840ee8c61a85ac99512c9fb57c9 * 1e7 - o_41974616d0958c9c2202672043152dea * 1e6) / 1e5);int o_e972171cfdb3dba6365e021488e9af3f=std::floor((o_2ae3af8759a9b8d041b10df90b88fc5f - o_34c768cab8370908fb6a38df742dad8e * 1e11 - o_613c4e96f69eecb7e11b5f936accebe3 * 1e10 - o_c77692bcee7d479651f5dee740ccb02b * 1e9 - o_4bc91391250ddcf73b7579af939e4961 * 1e8 - o_bdf65840ee8c61a85ac99512c9fb57c9 * 1e7 - o_41974616d0958c9c2202672043152dea * 1e6 - o_114558c036217de82543350ce875af8d * 1e5) / 1e4);int o_ae4f68034ba469ee30ad71006723192d=std::floor((o_2ae3af8759a9b8d041b10df90b88fc5f - o_34c768cab8370908fb6a38df742dad8e * 1e11 - o_613c4e96f69eecb7e11b5f936accebe3 * 1e10 - o_c77692bcee7d479651f5dee740ccb02b * 1e9 - o_4bc91391250ddcf73b7579af939e4961 * 1e8 - o_bdf65840ee8c61a85ac99512c9fb57c9 * 1e7 - o_41974616d0958c9c2202672043152dea * 1e6 - o_114558c036217de82543350ce875af8d * 1e5 - o_e972171cfdb3dba6365e021488e9af3f * 1e4) / 1e3);int o_8ce717954d6b2c9944a860f50b48fc05=std::floor((o_2ae3af8759a9b8d041b10df90b88fc5f - o_34c768cab8370908fb6a38df742dad8e * 1e11 - o_613c4e96f69eecb7e11b5f936accebe3 * 1e10 - o_c77692bcee7d479651f5dee740ccb02b * 1e9 - o_4bc91391250ddcf73b7579af939e4961 * 1e8 - o_bdf65840ee8c61a85ac99512c9fb57c9 * 1e7 - o_41974616d0958c9c2202672043152dea * 1e6 - o_114558c036217de82543350ce875af8d * 1e5 - o_e972171cfdb3dba6365e021488e9af3f * 1e4 - o_ae4f68034ba469ee30ad71006723192d * 1e3) / 1e2);int o_05d292526c8bdba94c7194bbf0f5481b=std::floor((o_2ae3af8759a9b8d041b10df90b88fc5f - o_34c768cab8370908fb6a38df742dad8e * 1e11 - o_613c4e96f69eecb7e11b5f936accebe3 * 1e10 - o_c77692bcee7d479651f5dee740ccb02b * 1e9 - o_4bc91391250ddcf73b7579af939e4961 * 1e8 - o_bdf65840ee8c61a85ac99512c9fb57c9 * 1e7 - o_41974616d0958c9c2202672043152dea * 1e6 - o_114558c036217de82543350ce875af8d * 1e5 - o_e972171cfdb3dba6365e021488e9af3f * 1e4 - o_ae4f68034ba469ee30ad71006723192d * 1e3 - o_8ce717954d6b2c9944a860f50b48fc05 * 1e2) / 1e1);int o_b278494c01370c0d7c319cfc3efafc16=std::floor((o_2ae3af8759a9b8d041b10df90b88fc5f - o_34c768cab8370908fb6a38df742dad8e * 1e11 - o_613c4e96f69eecb7e11b5f936accebe3 * 1e10 - o_c77692bcee7d479651f5dee740ccb02b * 1e9 - o_4bc91391250ddcf73b7579af939e4961 * 1e8 - o_bdf65840ee8c61a85ac99512c9fb57c9 * 1e7 - o_41974616d0958c9c2202672043152dea * 1e6 - o_114558c036217de82543350ce875af8d * 1e5 - o_e972171cfdb3dba6365e021488e9af3f * 1e4 - o_ae4f68034ba469ee30ad71006723192d * 1e3 - o_8ce717954d6b2c9944a860f50b48fc05 * 1e2 - o_05d292526c8bdba94c7194bbf0f5481b * 1e1) / 1e0);std::vector<int> o_77d081f3973e14b6415e5297954a555e={o_34c768cab8370908fb6a38df742dad8e,o_613c4e96f69eecb7e11b5f936accebe3,o_c77692bcee7d479651f5dee740ccb02b,o_4bc91391250ddcf73b7579af939e4961,o_bdf65840ee8c61a85ac99512c9fb57c9,o_41974616d0958c9c2202672043152dea,o_114558c036217de82543350ce875af8d,o_e972171cfdb3dba6365e021488e9af3f,o_ae4f68034ba469ee30ad71006723192d,o_8ce717954d6b2c9944a860f50b48fc05,o_05d292526c8bdba94c7194bbf0f5481b,o_b278494c01370c0d7c319cfc3efafc16};;int o_771356e26c1ada8226fca750129dfaab=std::floor(o_b34007ab812b0eb24031ea1ecb924afa / 1e11);int o_21cbc4ecf0667112f1e07b7773fad5bb=std::floor((o_b34007ab812b0eb24031ea1ecb924afa - o_771356e26c1ada8226fca750129dfaab * 1e11) / 1e10);int o_4d52eb8bef2529783533c89b4d1e0d01=std::floor((o_b34007ab812b0eb24031ea1ecb924afa - o_771356e26c1ada8226fca750129dfaab * 1e11 - o_21cbc4ecf0667112f1e07b7773fad5bb * 1e10) / 1e9);int o_b977c76215572dc5cc19620eb8e521d7=std::floor((o_b34007ab812b0eb24031ea1ecb924afa - o_771356e26c1ada8226fca750129dfaab * 1e11 - o_21cbc4ecf0667112f1e07b7773fad5bb * 1e10 - o_4d52eb8bef2529783533c89b4d1e0d01 * 1e9) / 1e8);std::vector<int> o_75a8be5447d67e4c7b68f3cf19c390c2={(0x0000000000000006 + 0x0000000000000203 + 0x0000000000000803 - 0x0000000000000A09),(0x0000000000000014 + 0x000000000000020A + 0x000000000000080A - 0x0000000000000A1E),(0x0000000000000004 + 0x0000000000000202 + 0x0000000000000802 - 0x0000000000000A06),(0x0000000000000002 + 0x0000000000000201 + 0x0000000000000801 - 0x0000000000000A03),(0x000000000000000E + 0x0000000000000207 + 0x0000000000000807 - 0x0000000000000A15),(0x000000000000000C + 0x0000000000000206 + 0x0000000000000806 - 0x0000000000000A12),(0x000000000000000A + 0x0000000000000205 + 0x0000000000000805 - 0x0000000000000A0F),(0x0000000000000010 + 0x0000000000000208 + 0x0000000000000808 - 0x0000000000000A18),(0x0000000000000012 + 0x0000000000000209 + 0x0000000000000809 - 0x0000000000000A1B),(0x0000000000000016 + 0x000000000000020B + 0x000000000000080B - 0x0000000000000A21),(0x0000000000000008 + 0x0000000000000204 + 0x0000000000000804 - 0x0000000000000A0C),(0x0000000000000018 + 0x000000000000020C + 0x000000000000080C - 0x0000000000000A24)};;int o_fa2dc814ccc0801f02462a3299d9c04a=std::floor((o_b34007ab812b0eb24031ea1ecb924afa - o_771356e26c1ada8226fca750129dfaab * 1e11 - o_21cbc4ecf0667112f1e07b7773fad5bb * 1e10 - o_4d52eb8bef2529783533c89b4d1e0d01 * 1e9 - o_b977c76215572dc5cc19620eb8e521d7 * 1e8) / 1e7);int o_e7085a3af6efdd745e04aceb5ee2018f=std::floor((o_b34007ab812b0eb24031ea1ecb924afa - o_771356e26c1ada8226fca750129dfaab * 1e11 - o_21cbc4ecf0667112f1e07b7773fad5bb * 1e10 - o_4d52eb8bef2529783533c89b4d1e0d01 * 1e9 - o_b977c76215572dc5cc19620eb8e521d7 * 1e8 - o_fa2dc814ccc0801f02462a3299d9c04a * 1e7) / 1e6);int o_a70c86b40df5a6148007c1b4c94fdc65=std::floor((o_b34007ab812b0eb24031ea1ecb924afa - o_771356e26c1ada8226fca750129dfaab * 1e11 - o_21cbc4ecf0667112f1e07b7773fad5bb * 1e10 - o_4d52eb8bef2529783533c89b4d1e0d01 * 1e9 - o_b977c76215572dc5cc19620eb8e521d7 * 1e8 - o_fa2dc814ccc0801f02462a3299d9c04a * 1e7 - o_e7085a3af6efdd745e04aceb5ee2018f * 1e6) / 1e5);int o_ad54fa572a0ad69ff220442957065ddd=std::floor((o_b34007ab812b0eb24031ea1ecb924afa - o_771356e26c1ada8226fca750129dfaab * 1e11 - o_21cbc4ecf0667112f1e07b7773fad5bb * 1e10 - o_4d52eb8bef2529783533c89b4d1e0d01 * 1e9 - o_b977c76215572dc5cc19620eb8e521d7 * 1e8 - o_fa2dc814ccc0801f02462a3299d9c04a * 1e7 - o_e7085a3af6efdd745e04aceb5ee2018f * 1e6 - o_a70c86b40df5a6148007c1b4c94fdc65 * 1e5) / 1e4);int o_36f864508ffd9f0b4fdc88060168ca0b=std::floor((o_b34007ab812b0eb24031ea1ecb924afa - o_771356e26c1ada8226fca750129dfaab * 1e11 - o_21cbc4ecf0667112f1e07b7773fad5bb * 1e10 - o_4d52eb8bef2529783533c89b4d1e0d01 * 1e9 - o_b977c76215572dc5cc19620eb8e521d7 * 1e8 - o_fa2dc814ccc0801f02462a3299d9c04a * 1e7 - o_e7085a3af6efdd745e04aceb5ee2018f * 1e6 - o_a70c86b40df5a6148007c1b4c94fdc65 * 1e5 - o_ad54fa572a0ad69ff220442957065ddd * 1e4) / 1e3);std::vector<int> o_61f3029b9cdcb4ea2a973954ae2017a4={(0x000000000000000C + 0x0000000000000206 + 0x0000000000000806 - 0x0000000000000A12),(0x000000000000000E + 0x0000000000000207 + 0x0000000000000807 - 0x0000000000000A15),(0x0000000000000016 + 0x000000000000020B + 0x000000000000080B - 0x0000000000000A21),(0x0000000000000004 + 0x0000000000000202 + 0x0000000000000802 - 0x0000000000000A06),(0x0000000000000010 + 0x0000000000000208 + 0x0000000000000808 - 0x0000000000000A18),(0x0000000000000018 + 0x000000000000020C + 0x000000000000080C - 0x0000000000000A24),(0x000000000000000A + 0x0000000000000205 + 0x0000000000000805 - 0x0000000000000A0F),(0x0000000000000012 + 0x0000000000000209 + 0x0000000000000809 - 0x0000000000000A1B),(0x0000000000000002 + 0x0000000000000201 + 0x0000000000000801 - 0x0000000000000A03),(0x0000000000000008 + 0x0000000000000204 + 0x0000000000000804 - 0x0000000000000A0C),(0x0000000000000014 + 0x000000000000020A + 0x000000000000080A - 0x0000000000000A1E),(0x0000000000000006 + 0x0000000000000203 + 0x0000000000000803 - 0x0000000000000A09)};;int o_6c45929d034768cf33e7ce4b9c46ad15=std::floor((o_b34007ab812b0eb24031ea1ecb924afa - o_771356e26c1ada8226fca750129dfaab * 1e11 - o_21cbc4ecf0667112f1e07b7773fad5bb * 1e10 - o_4d52eb8bef2529783533c89b4d1e0d01 * 1e9 - o_b977c76215572dc5cc19620eb8e521d7 * 1e8 - o_fa2dc814ccc0801f02462a3299d9c04a * 1e7 - o_e7085a3af6efdd745e04aceb5ee2018f * 1e6 - o_a70c86b40df5a6148007c1b4c94fdc65 * 1e5 - o_ad54fa572a0ad69ff220442957065ddd * 1e4 - o_36f864508ffd9f0b4fdc88060168ca0b * 1e3) / 1e2);int o_72add2af5956ed88f11bf2e60bbce7c2=std::floor((o_b34007ab812b0eb24031ea1ecb924afa - o_771356e26c1ada8226fca750129dfaab * 1e11 - o_21cbc4ecf0667112f1e07b7773fad5bb * 1e10 - o_4d52eb8bef2529783533c89b4d1e0d01 * 1e9 - o_b977c76215572dc5cc19620eb8e521d7 * 1e8 - o_fa2dc814ccc0801f02462a3299d9c04a * 1e7 - o_e7085a3af6efdd745e04aceb5ee2018f * 1e6 - o_a70c86b40df5a6148007c1b4c94fdc65 * 1e5 - o_ad54fa572a0ad69ff220442957065ddd * 1e4 - o_36f864508ffd9f0b4fdc88060168ca0b * 1e3 - o_6c45929d034768cf33e7ce4b9c46ad15 * 1e2) / 1e1);int o_e5e630d4eb0826bb073744bd47ea126e=std::floor((o_b34007ab812b0eb24031ea1ecb924afa - o_771356e26c1ada8226fca750129dfaab * 1e11 - o_21cbc4ecf0667112f1e07b7773fad5bb * 1e10 - o_4d52eb8bef2529783533c89b4d1e0d01 * 1e9 - o_b977c76215572dc5cc19620eb8e521d7 * 1e8 - o_fa2dc814ccc0801f02462a3299d9c04a * 1e7 - o_e7085a3af6efdd745e04aceb5ee2018f * 1e6 - o_a70c86b40df5a6148007c1b4c94fdc65 * 1e5 - o_ad54fa572a0ad69ff220442957065ddd * 1e4 - o_36f864508ffd9f0b4fdc88060168ca0b * 1e3 - o_6c45929d034768cf33e7ce4b9c46ad15 * 1e2 - o_72add2af5956ed88f11bf2e60bbce7c2 * 1e1) / 1e0);std::vector<int> o_6f944baf3162eac457ff48656805c26f={o_771356e26c1ada8226fca750129dfaab,o_21cbc4ecf0667112f1e07b7773fad5bb,o_4d52eb8bef2529783533c89b4d1e0d01,o_b977c76215572dc5cc19620eb8e521d7,o_fa2dc814ccc0801f02462a3299d9c04a,o_e7085a3af6efdd745e04aceb5ee2018f,o_a70c86b40df5a6148007c1b4c94fdc65,o_ad54fa572a0ad69ff220442957065ddd,o_36f864508ffd9f0b4fdc88060168ca0b,o_6c45929d034768cf33e7ce4b9c46ad15,o_72add2af5956ed88f11bf2e60bbce7c2,o_e5e630d4eb0826bb073744bd47ea126e};;for (int o_c7e55f3c8471b4b4d2f789700472a07e=(0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00);(o_c7e55f3c8471b4b4d2f789700472a07e < o_318408d0fdfbad67cf162f7f5b9b6ea8) & !!(o_c7e55f3c8471b4b4d2f789700472a07e < o_318408d0fdfbad67cf162f7f5b9b6ea8);o_c7e55f3c8471b4b4d2f789700472a07e++){int o_e2b3f62dbbe6a2fc12126292dd571944=o_e4fa61c75340c147340e36093deac9bf.at(o_c7e55f3c8471b4b4d2f789700472a07e);int o_b09c8dc18f1bc9d8c7ac26eb4b86bbd4=o_6f944baf3162eac457ff48656805c26f.at(o_c7e55f3c8471b4b4d2f789700472a07e);o_e4fa61c75340c147340e36093deac9bf.at(o_c7e55f3c8471b4b4d2f789700472a07e) = (o_e2b3f62dbbe6a2fc12126292dd571944 - o_b09c8dc18f1bc9d8c7ac26eb4b86bbd4 + (0x0000000000000014 + 0x000000000000020A + 0x000000000000080A - 0x0000000000000A1E)) % (0x0000000000000014 + 0x000000000000020A + 0x000000000000080A - 0x0000000000000A1E);};int o_8d7d9af2d6b1dee452819a620e96c660=o_fd533d5f10525a0cd35b5eaa5b805b43((0x0000000000000002 + 0x0000000000000201 + 0x0000000000000801 - 0x0000000000000A03));std::vector<int> o_940f9a567335aee3cfc4a8b948fef0f9(o_318408d0fdfbad67cf162f7f5b9b6ea8);for (int o_878c445b63622866c74264cf2447ed63=(0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00);(o_878c445b63622866c74264cf2447ed63 < o_318408d0fdfbad67cf162f7f5b9b6ea8) & !!(o_878c445b63622866c74264cf2447ed63 < o_318408d0fdfbad67cf162f7f5b9b6ea8);o_878c445b63622866c74264cf2447ed63++){int o_3298627a77ef0d0c1c6d36be52c44a34=o_5d8a18909113eff205d13843050c6b71.at(o_878c445b63622866c74264cf2447ed63) - (0x0000000000000002 + 0x0000000000000201 + 0x0000000000000801 - 0x0000000000000A03);o_940f9a567335aee3cfc4a8b948fef0f9.at(o_3298627a77ef0d0c1c6d36be52c44a34) = o_e4fa61c75340c147340e36093deac9bf.at(o_878c445b63622866c74264cf2447ed63);};o_e4fa61c75340c147340e36093deac9bf=o_940f9a567335aee3cfc4a8b948fef0f9;for (int o_b12db3c753756f1a215ac751d3574061=(0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00);(o_b12db3c753756f1a215ac751d3574061 < o_318408d0fdfbad67cf162f7f5b9b6ea8) & !!(o_b12db3c753756f1a215ac751d3574061 < o_318408d0fdfbad67cf162f7f5b9b6ea8);o_b12db3c753756f1a215ac751d3574061++){int o_559f73f7f079255b4f029947a0709605=o_e4fa61c75340c147340e36093deac9bf.at(o_b12db3c753756f1a215ac751d3574061);int o_63d2fd683220815e55255e9c2006368d=o_77d081f3973e14b6415e5297954a555e.at(o_b12db3c753756f1a215ac751d3574061);o_e4fa61c75340c147340e36093deac9bf.at(o_b12db3c753756f1a215ac751d3574061) = (o_559f73f7f079255b4f029947a0709605 - o_63d2fd683220815e55255e9c2006368d + (0x0000000000000014 + 0x000000000000020A + 0x000000000000080A - 0x0000000000000A1E)) % (0x0000000000000014 + 0x000000000000020A + 0x000000000000080A - 0x0000000000000A1E);};for (int o_b38568d68a05810f4877d208b2de6d12=(0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00);(o_b38568d68a05810f4877d208b2de6d12 < o_318408d0fdfbad67cf162f7f5b9b6ea8) & !!(o_b38568d68a05810f4877d208b2de6d12 < o_318408d0fdfbad67cf162f7f5b9b6ea8);o_b38568d68a05810f4877d208b2de6d12++){int o_a65977522d39e4d2fea324e0fa97d8f6=o_75a8be5447d67e4c7b68f3cf19c390c2.at(o_b38568d68a05810f4877d208b2de6d12) - (0x0000000000000002 + 0x0000000000000201 + 0x0000000000000801 - 0x0000000000000A03);o_940f9a567335aee3cfc4a8b948fef0f9.at(o_a65977522d39e4d2fea324e0fa97d8f6) = o_e4fa61c75340c147340e36093deac9bf.at(o_b38568d68a05810f4877d208b2de6d12);};o_e4fa61c75340c147340e36093deac9bf=o_940f9a567335aee3cfc4a8b948fef0f9;int o_b670a59ae3f79d9d19175ab1b6c726ab=o_fd533d5f10525a0cd35b5eaa5b805b43((0x0000000000000004 + 0x0000000000000202 + 0x0000000000000802 - 0x0000000000000A06));for (int o_19fcd77902cfde47ff23bd4aa2b5c515=(0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00);(o_19fcd77902cfde47ff23bd4aa2b5c515 < o_318408d0fdfbad67cf162f7f5b9b6ea8) & !!(o_19fcd77902cfde47ff23bd4aa2b5c515 < o_318408d0fdfbad67cf162f7f5b9b6ea8);o_19fcd77902cfde47ff23bd4aa2b5c515++){int o_e66ad059dee4bad01d6b993dc9ef1db9=o_e4fa61c75340c147340e36093deac9bf.at(o_19fcd77902cfde47ff23bd4aa2b5c515);int o_7cc0b798e6934fa598531297e17e5c67=o_56af41e2fad4b840fe45d58c27e703de.at(o_19fcd77902cfde47ff23bd4aa2b5c515);o_e4fa61c75340c147340e36093deac9bf.at(o_19fcd77902cfde47ff23bd4aa2b5c515) = (o_e66ad059dee4bad01d6b993dc9ef1db9 - o_7cc0b798e6934fa598531297e17e5c67 + (0x0000000000000014 + 0x000000000000020A + 0x000000000000080A - 0x0000000000000A1E)) % (0x0000000000000014 + 0x000000000000020A + 0x000000000000080A - 0x0000000000000A1E);};int o_4588e12db8f4b0e27a01915ce7d01ad8=o_fd533d5f10525a0cd35b5eaa5b805b43((0x0000000000000006 + 0x0000000000000203 + 0x0000000000000803 - 0x0000000000000A09));for (int o_fc0c888ed14b955530295fcaedc0bae0=(0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00);(o_fc0c888ed14b955530295fcaedc0bae0 < o_318408d0fdfbad67cf162f7f5b9b6ea8) & !!(o_fc0c888ed14b955530295fcaedc0bae0 < o_318408d0fdfbad67cf162f7f5b9b6ea8);o_fc0c888ed14b955530295fcaedc0bae0++){int o_c906bef254d43b3d9330eb7a92a45f97=o_61f3029b9cdcb4ea2a973954ae2017a4.at(o_fc0c888ed14b955530295fcaedc0bae0) - (0x0000000000000002 + 0x0000000000000201 + 0x0000000000000801 - 0x0000000000000A03);o_940f9a567335aee3cfc4a8b948fef0f9.at(o_c906bef254d43b3d9330eb7a92a45f97) = o_e4fa61c75340c147340e36093deac9bf.at(o_fc0c888ed14b955530295fcaedc0bae0);};o_e4fa61c75340c147340e36093deac9bf=o_940f9a567335aee3cfc4a8b948fef0f9;int o_f5bc9f726135bdc5b86dee2a78c4de22=(0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00);for (int o_ae23cb22ef623fd1949f823e56730699=(0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00);(o_ae23cb22ef623fd1949f823e56730699 < o_86838dd8a6e6975fd7b157415f5b21d8) & !!(o_ae23cb22ef623fd1949f823e56730699 < o_86838dd8a6e6975fd7b157415f5b21d8);o_ae23cb22ef623fd1949f823e56730699++){o_f5bc9f726135bdc5b86dee2a78c4de22 = o_f5bc9f726135bdc5b86dee2a78c4de22 + o_e4fa61c75340c147340e36093deac9bf.at(o_ae23cb22ef623fd1949f823e56730699) * pow((0x0000000000000014 + 0x000000000000020A + 0x000000000000080A - 0x0000000000000A1E),(o_86838dd8a6e6975fd7b157415f5b21d8 - o_ae23cb22ef623fd1949f823e56730699 - (0x0000000000000002 + 0x0000000000000201 + 0x0000000000000801 - 0x0000000000000A03)));};int o_310a827dab0df94f18629b376af4a2f4=(0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00);for (int o_f30bac3b0d9d0f72a5bdeecbe87c1e9d=(0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00);(o_f30bac3b0d9d0f72a5bdeecbe87c1e9d < (0x0000000000000004 + 0x0000000000000202 + 0x0000000000000802 - 0x0000000000000A06)) & !!(o_f30bac3b0d9d0f72a5bdeecbe87c1e9d < (0x0000000000000004 + 0x0000000000000202 + 0x0000000000000802 - 0x0000000000000A06));o_f30bac3b0d9d0f72a5bdeecbe87c1e9d++){o_310a827dab0df94f18629b376af4a2f4 = o_310a827dab0df94f18629b376af4a2f4 + o_e4fa61c75340c147340e36093deac9bf.at(o_86838dd8a6e6975fd7b157415f5b21d8 + o_f30bac3b0d9d0f72a5bdeecbe87c1e9d) * pow((0x0000000000000014 + 0x000000000000020A + 0x000000000000080A - 0x0000000000000A1E),((0x0000000000000002 + 0x0000000000000201 + 0x0000000000000801 - 0x0000000000000A03) - o_f30bac3b0d9d0f72a5bdeecbe87c1e9d));};bool o_b5e5760753a68ca0b8f8fdc523d301a1=false;int o_770a10c4b6da26a86813cb2b6b073912=(0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00);for (int o_64aeaf342f7a5d5d8d28cd8b1c8c7b59=(0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00);(o_64aeaf342f7a5d5d8d28cd8b1c8c7b59 < (0x0000000000000004 + 0x0000000000000202 + 0x0000000000000802 - 0x0000000000000A06)) & !!(o_64aeaf342f7a5d5d8d28cd8b1c8c7b59 < (0x0000000000000004 + 0x0000000000000202 + 0x0000000000000802 - 0x0000000000000A06));o_64aeaf342f7a5d5d8d28cd8b1c8c7b59++){o_770a10c4b6da26a86813cb2b6b073912 = o_770a10c4b6da26a86813cb2b6b073912 + o_e4fa61c75340c147340e36093deac9bf.at(o_86838dd8a6e6975fd7b157415f5b21d8 + (0x0000000000000004 + 0x0000000000000202 + 0x0000000000000802 - 0x0000000000000A06) + o_64aeaf342f7a5d5d8d28cd8b1c8c7b59) * pow((0x0000000000000014 + 0x000000000000020A + 0x000000000000080A - 0x0000000000000A1E),((0x0000000000000002 + 0x0000000000000201 + 0x0000000000000801 - 0x0000000000000A03) - o_64aeaf342f7a5d5d8d28cd8b1c8c7b59));};int o_2b86bb457795751067059d93f62d8f7e=(0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00);for (int o_351e19021bcd6433f6a2bd8abfdfc250=(0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00);(o_351e19021bcd6433f6a2bd8abfdfc250 < (0x0000000000000006 + 0x0000000000000203 + 0x0000000000000803 - 0x0000000000000A09)) & !!(o_351e19021bcd6433f6a2bd8abfdfc250 < (0x0000000000000006 + 0x0000000000000203 + 0x0000000000000803 - 0x0000000000000A09));o_351e19021bcd6433f6a2bd8abfdfc250++){o_2b86bb457795751067059d93f62d8f7e = o_2b86bb457795751067059d93f62d8f7e + o_e4fa61c75340c147340e36093deac9bf.at(o_86838dd8a6e6975fd7b157415f5b21d8 + (0x0000000000000004 + 0x0000000000000202 + 0x0000000000000802 - 0x0000000000000A06) + (0x0000000000000004 + 0x0000000000000202 + 0x0000000000000802 - 0x0000000000000A06) + o_351e19021bcd6433f6a2bd8abfdfc250) * pow((0x0000000000000014 + 0x000000000000020A + 0x000000000000080A - 0x0000000000000A1E),((0x0000000000000004 + 0x0000000000000202 + 0x0000000000000802 - 0x0000000000000A06) - o_351e19021bcd6433f6a2bd8abfdfc250));};int o_a7427543391d61e0f967463208a43225=(0x0000000000000FA0 + 0x00000000000009D0 + 0x0000000000000FD0 - 0x0000000000002170) + o_2b86bb457795751067059d93f62d8f7e;if ((o_a7427543391d61e0f967463208a43225 > o_4588e12db8f4b0e27a01915ce7d01ad8) & !!(o_a7427543391d61e0f967463208a43225 > o_4588e12db8f4b0e27a01915ce7d01ad8)){o_b5e5760753a68ca0b8f8fdc523d301a1 = true;}else if (o_a7427543391d61e0f967463208a43225 == o_4588e12db8f4b0e27a01915ce7d01ad8){if ((o_770a10c4b6da26a86813cb2b6b073912 > o_b670a59ae3f79d9d19175ab1b6c726ab) & !!(o_770a10c4b6da26a86813cb2b6b073912 > o_b670a59ae3f79d9d19175ab1b6c726ab)){o_b5e5760753a68ca0b8f8fdc523d301a1 = true;}else if (o_770a10c4b6da26a86813cb2b6b073912 == o_b670a59ae3f79d9d19175ab1b6c726ab){if ((o_310a827dab0df94f18629b376af4a2f4 >= o_8d7d9af2d6b1dee452819a620e96c660) & !!(o_310a827dab0df94f18629b376af4a2f4 >= o_8d7d9af2d6b1dee452819a620e96c660)){o_b5e5760753a68ca0b8f8fdc523d301a1 = true;};};;};;if (o_f5bc9f726135bdc5b86dee2a78c4de22 == o_399493a4b6ac5ef205ba007c1c8364e7){if (o_b5e5760753a68ca0b8f8fdc523d301a1){o_3b60b4781780d8b410b082ce2879e898 = true;if(o_a578ef72809bf557d02f76ca5aaf269e){Rcout <<"\x4C""i\143e\x6E""s\145 \x76""a\154i\x64"" \165n\x74""i\154 "<<o_a7427543391d61e0f967463208a43225<<"\x2D"""<<o_770a10c4b6da26a86813cb2b6b073912<<"\x2D"""<<o_310a827dab0df94f18629b376af4a2f4<<"\x20""(\171y\x79""y\055m\x6D""-\144d\x29""" << std::endl;};}else {if(o_a578ef72809bf557d02f76ca5aaf269e){Rcout <<"\x4C""i\143e\x6E""s\145 \x68""a\163 \x65""x\160i\x72""e\144." << std::endl;};};}else {if(o_a578ef72809bf557d02f76ca5aaf269e){Rcout <<"\x49""n\166a\x6C""i\144 \x75""s\145r\x6E""a\155e\x20""o\162 \x70""a\163s\x77""o\162d" << std::endl;};};if (!o_3b60b4781780d8b410b082ce2879e898){std::this_thread::sleep_for(std::chrono::milliseconds(500));if(o_a578ef72809bf557d02f76ca5aaf269e){Rcout<<"\x49""n\166a\x6C""i\144 \x6C""i\143e\x6E""s\145.\x20""C\157n\x74""i\156u\x69""n\147 \x77""i\164h\x20""d\145m\x6F"" \166e\x72""s\151o\x6E""" << std::endl;Rcout << "\x49""f\040y\x6F""u\040d\x65""s\151r\x65"" \141 \x6E""e\167 \x6C""i\143e\x6E""s\145,\x20""p\154e\x61""s\145 \x63""o\156t\x61""c\164 \x75""s\040a\x74"" \152.\x69""s\151d\x72""o\100u\x70""m\056e\x73""" << std::endl;};};return o_3b60b4781780d8b410b082ce2879e898;};

class Population{
private:
  int npop;
  int nelite;
  int nchrom;
  IntegerVector chromsizes;
  CharacterVector chromtypes;

  STATCLASS StatClass;
  vector<double> FitnessVals;

  vector<IntegerMatrix> BOOL;
  vector<IntegerMatrix> OS;
  vector<IntegerMatrix> UOS;
  vector<IntegerMatrix> OMS;
  vector<IntegerMatrix> UOMS;
  vector<NumericMatrix> DBL;

  int nBOOL=0;
  int nOS=0;
  int nUOS=0;
  int nOMS=0;
  int nUOMS=0;
  int nDBL=0;
  int nINT=nBOOL+nOS+nUOS+nOMS+nUOMS;

  IntegerVector nvecBOOL;
  IntegerVector nvecOS;
  IntegerVector nvecUOS;
  IntegerVector nvecOMS;
  IntegerVector nvecUOMS;
  IntegerVector nvecDBL;

  vector<IntegerVector> CandOS;
  vector<IntegerVector>  CandUOS;
  vector<IntegerVector>  CandOMS;
  vector<IntegerVector>  CandUOMS;
  vector<NumericVector> CandDBL;


public:
  //////////////
  IntegerVector CandBOOL={1,0};

  IntegerVector OrderPop;

  Population(){
  }

  void set_npop(int npop_){
    npop=npop_;
  }

  void set_nelite(int nelite_){
    nelite=nelite_;
  }
  void set_nchrom(int nchrom_){
    nchrom=nchrom_;
  }
  void set_chromsizes(IntegerVector chromsizes_){
    chromsizes=chromsizes_;
  }
  void set_chromtypes(CharacterVector chromtypes_){
    chromtypes=chromtypes_;
  }
  /////////////////
  int get_npop( ){
    return npop;
  }

  int get_nelite( ){
    return nelite;
  }
  int get_nchrom( ){
    return nchrom;
  }
  IntegerVector get_chromsizes( ){
    return chromsizes;
  }
  CharacterVector get_chromtypes(){
    return chromtypes;
  }


  //!!!
  int get_nINT( ){
    return nINT;
  }

    int get_nDBL( ){
    return nDBL;
  }



  /////////////

  void push_back_CandOS(IntegerVector Cand){
    CandOS.push_back(Cand);
  }
  void push_back_CandUOS(IntegerVector Cand){
    CandUOS.push_back(Cand);
  }
  void push_back_CandOMS(IntegerVector Cand){
    CandOMS.push_back(Cand);
  }
  void push_back_CandUOMS(IntegerVector Cand){
    CandUOMS.push_back(Cand);
  }
  void push_back_CandDBL(NumericVector Cand){
    CandDBL.push_back(Cand);
  }

  IntegerMatrix get_BOOL(int i){
    return BOOL.at(i);
  }

  IntegerMatrix get_OS(int i){
    return OS.at(i);
  }
  IntegerMatrix get_UOS(int i){
    return UOS.at(i);
  }
  IntegerMatrix get_OMS(int i){
    return OMS.at(i);
  }
  IntegerMatrix get_UOMS(int i){
    return UOMS.at(i);
  }
  NumericMatrix get_DBL(int i){
    return DBL.at(i);
  }


  //////////////////////////
  void InitRand(){
    int iiBOOL=0;
    for (int i=0;i<nchrom;i++)
      if (chromtypes[i]=="BOOL"){
        IntegerMatrix TempMat(chromsizes[i],npop+nelite);
        for (int j=0;j<(npop+nelite);j++){
          TempMat(_,j)=sample(CandBOOL, chromsizes[i],true);
        }
        BOOL.push_back(TempMat);
        iiBOOL++;
      }
      nBOOL=iiBOOL;
      int iiOS=0;
      for (int i=0;i<nchrom;i++)
        if (chromtypes[i]=="OS"){
          IntegerMatrix TempMat(chromsizes[i],npop+nelite);
          for (int j=0;j<(npop+nelite);j++){
            TempMat(_,j)=sample(CandOS.at(iiOS), chromsizes[i]);
          }
          OS.push_back(TempMat);
          iiOS++;
        }
        nOS=iiOS;
        int iiUOS=0;
        for (int i=0;i<nchrom;i++)
          if (chromtypes[i]=="UOS"){
            IntegerMatrix TempMat(chromsizes[i],(npop+nelite));
            for (int j=0;j<(npop+nelite);j++){
              TempMat(_,j)=sample(CandUOS.at(iiUOS), chromsizes[i]).sort();
            }
            UOS.push_back(TempMat);
            iiUOS++;
          }
          nUOS=iiUOS;
          int iiOMS=0;
          for (int i=0;i<nchrom;i++)
            if (chromtypes[i]=="OMS"){
              IntegerMatrix TempMat(chromsizes[i],(npop+nelite));
              for (int j=0;j<(npop+nelite);j++){
                TempMat(_,j)=sample(CandOMS.at(iiOMS), chromsizes[i], true);
              }
              OMS.push_back(TempMat);
              iiOMS++;
            }
            nOMS=iiOMS;
            int iiUOMS=0;
            for (int i=0;i<nchrom;i++)
              if (chromtypes[i]=="UOMS"){
                IntegerMatrix TempMat(chromsizes[i],(npop+nelite));
                for (int j=0;j<(npop+nelite);j++){
                  TempMat(_,j)=sample(CandUOMS.at(iiUOMS), chromsizes[i], true).sort();
                }
                UOMS.push_back(TempMat);
                iiUOMS++;
              }
              nUOMS=iiUOMS;
              int iiDBL=0;
              for (int i=0;i<nchrom;i++)
                if (chromtypes[i]=="DBL"){
                  NumericMatrix TempMat(chromsizes[i],(npop+nelite));
                  for (int j=0;j<(npop+nelite);j++){
                    TempMat(_,j)=runif(chromsizes[i])*(CandDBL.at(iiDBL)(1)-CandDBL.at(iiDBL)(0))+CandDBL.at(iiDBL)(0);
                  }
                  DBL.push_back(TempMat);
                  iiDBL++;
                }
                nDBL=iiDBL;
                nINT=nBOOL+nOS+nUOS+nOMS+nUOMS;
  }



  void Init(vector<IntegerMatrix> BOOL_, vector<IntegerMatrix> OS_,vector<IntegerMatrix> UOS_,vector<IntegerMatrix> OMS_,vector<IntegerMatrix> UOMS_, vector<NumericMatrix>DBL_){
    BOOL=BOOL_;
    OS=OS_;
    UOS=UOS_;
    OMS=OMS_;
    UOMS=UOMS_;
    DBL=DBL_;
  }




  void MoveInd(int from, int to){
    if (nBOOL>0){
      for (int i=0;i<nBOOL;i++){
        BOOL.at(i)(_, to)=BOOL.at(i)(_, from);
      }
    }
    if (nOS>0){
      for (int i=0;i<nOS;i++){
        OS.at(i)(_, to)=OS.at(i)(_, from);
      }
    }
    if (nUOS>0){
      for (int i=0;i<nUOS;i++){
        UOS.at(i)(_, to)=UOS.at(i)(_, from);
      }
    }
    if (nOMS>0){
      for (int i=0;i<nOMS;i++){
        OMS.at(i)(_, to)=OMS.at(i)(_, from);
      }
    }
    if (nUOMS>0){
      for (int i=0;i<nUOMS;i++){
        UOMS.at(i)(_, to)=UOMS.at(i)(_, from);
      }
    }
    if (nDBL>0){
      for (int i=0;i<nDBL;i++){
        DBL.at(i)(_, to)=DBL.at(i)(_, from);
      }
    }
  }



  void MakeCross(int p1, int p2, int child){

    if (nBOOL>0){
      for (int i=0;i<nBOOL;i++){
        int BOOLirows=BOOL.at(i).nrow();

        IntegerVector TempSol;
        for (int j=0;j<BOOLirows;j++){
          int sampleint=sample(2,1)[0];

          if (sampleint==1){
            TempSol.push_back(BOOL.at(i)(j,p1));
          } else {
            TempSol.push_back(BOOL.at(i)(j,p2));
          }
        }
        BOOL.at(i)(_,child)=TempSol;
      }
    }

    if (nOS>0){
      for (int i=0;i<nOS;i++){
        int OSirows=OS.at(i).nrow();

        IntegerVector TempSol;
        for (int j=0;j<OSirows;j++){
          int sampleint=sample(2,1)[0];

          if (sampleint==1 && !contains(TempSol,OS.at(i)(j,p1))){
            TempSol.push_back(OS.at(i)(j,p1));
          } else if (sampleint==2 && !contains(TempSol,OS.at(i)(j,p2))){
            TempSol.push_back(OS.at(i)(j,p2));
          } else {
            TempSol.push_back(sample(setdiff(CandOS.at(i),TempSol),1)[0]);
          }
        }
        OS.at(i)(_,child)=TempSol;
      }
    }

    if (nUOS>0){
      for (int i=0;i<nUOS;i++){
        int UOSirows=UOS.at(i).nrow();
        IntegerVector p1vec=UOS.at(i)(_,p1);
        IntegerVector p2vec=UOS.at(i)(_,p2);
        UOS.at(i)(_,child)=sample(union_(p1vec,p2vec),UOSirows, false).sort();
      }
    }

    if (nOMS>0){
      for (int i=0;i<nOMS;i++){
        int OMSirows=OMS.at(i).nrow();
        IntegerVector TempSol;
        for (int j=0;j<OMSirows;j++){
          int sampleint=sample(2,1)[0];
          if (sampleint==1){
            TempSol.push_back(OMS.at(i)(j,p1));
          } else{
            TempSol.push_back(OMS.at(i)(j,p2));
          }
        }
        OMS.at(i)(_,child)=TempSol;
      }
    }

    if (nUOMS>0){
      for (int i=0;i<nUOMS;i++){
        int UOMSirows=UOMS.at(i).nrow();
        IntegerVector TempSol;
        for (int j=0;j<UOMSirows;j++){
          TempSol.push_back(UOMS.at(i)(j,p1));
          TempSol.push_back(UOMS.at(i)(j,p2));
        }
        UOMS.at(i)(_,child)=sample(TempSol, UOMSirows, true).sort();
      }
    }

    if (nDBL>0){
      for (int i=0;i<nDBL;i++){
        int DBLirows=DBL.at(i).nrow();
        NumericVector TempSol;
        for (int j=0;j<DBLirows;j++){
          double rnum= runif(1)(0);
          TempSol.push_back(rnum*(DBL.at(i)(j,p1)+(1-rnum)*DBL.at(i)(j,p2)));
        }
        DBL.at(i)(_,child)=TempSol;
      }
    }

  }


  void  Mutate(int ind, double MUTPROB){


    if (nBOOL>0){
      for (int i=0;i<nBOOL;i++){
        int BOOLirows=BOOL.at(i).nrow();
        IntegerVector IndSol=BOOL.at(i)(_,ind);
        for (int j=0;j<BOOLirows;j++){
          double mutp=runif(1)(0);
          if (mutp<MUTPROB){
            IntegerVector totakeout;
            totakeout.push_back(IndSol(j));
            int replacement=sample(setdiff(CandBOOL,totakeout),1)(0);
            IndSol(j)=replacement;
          }
        }
        BOOL.at(i)(_,ind)=IndSol;
      }
    }

    if (nOS>0){
      for (int i=0;i<nOS;i++){
        int OSirows=OS.at(i).nrow();
        IntegerVector IndSol=OS.at(i)(_,ind);
        for (int j=0;j<OSirows;j++){
          double mutp=runif(1)(0);
          if (mutp<MUTPROB){
            IntegerVector totakeout;
            totakeout.push_back(IndSol(j));
            int replacement=sample(setdiff(CandOS.at(i),setdiff(IndSol,totakeout)),1)(0);
            IndSol(j)=replacement;
          }
        }
        double swapp=runif(1)(0);
        if (swapp<MUTPROB){
          int i1=sample(OSirows,1)(0)-1;
          int i2=sample(OSirows,1)(0)-1;
          int ii1=IndSol(i1);
          int ii2=IndSol(i2);

          IndSol(i1)= ii2;
          IndSol(i2)= ii1;

        }

        double slidep=runif(1)(0);
        if (slidep<MUTPROB){
          int movedirection=sample(2,1)(0);
          if (movedirection==1){
            std::rotate(IndSol.begin(), IndSol.begin() + 1, IndSol.end());
          } else{
            std::rotate(IndSol.begin(), IndSol.end(), IndSol.end());
          }
        }
        OS.at(i)(_,ind)=IndSol;
      }

    }

    if (nUOS>0){
      for (int i=0;i<nUOS;i++){
        int UOSirows=UOS.at(i).nrow();
        IntegerVector IndSol=UOS.at(i)(_,ind);
        for (int j=0;j<UOSirows;j++){
          double mutp=runif(1)(0);
          if (mutp<MUTPROB){
            int replacement=sample(setdiff(CandUOS.at(i),IndSol),1)(0);
            IndSol(j)=replacement;
          }
        }
        UOS.at(i)(_,ind)=IndSol.sort();
      }
    }

    if (nOMS>0){
      for (int i=0;i<nOMS;i++){
        int OMSirows=OMS.at(i).nrow();
        IntegerVector IndSol=OMS.at(i)(_,ind);
        for (int j=0;j<OMSirows;j++){
          double mutp=runif(1)(0);
          if (mutp<MUTPROB){
            int replacement=sample(CandOMS.at(i),1)(0);
            IndSol(j)=replacement;
          }
        }
        double swapp=runif(1)(0);
        if (swapp<MUTPROB){
          int i1=sample(OMSirows,1)(0)-1;
          int i2=sample(OMSirows,1)(0)-1;
          int ii1=IndSol(i1);
          int ii2=IndSol(i2);

          IndSol(i1)= ii2;
          IndSol(i2)= ii1;

        }

        double slidep=runif(1)(0);
        if (slidep<MUTPROB){
          int movedirection=sample(2,1)(0);
          if (movedirection==1){
            std::rotate(IndSol.begin(), IndSol.begin() + 1, IndSol.end());
          } else{
            std::rotate(IndSol.begin(), IndSol.end(), IndSol.end());
          }
        }
        OMS.at(i)(_,ind)=IndSol;
      }

    }


    if (nUOMS>0){
      for (int i=0;i<nUOMS;i++){
        int UOMSirows=UOMS.at(i).nrow();
        IntegerVector IndSol=UOMS.at(i)(_,ind);
        for (int j=0;j<UOMSirows;j++){
          double mutp=runif(1)(0);
          if (mutp<MUTPROB){
            int replacement=sample(CandUOMS.at(i),1)(0);
            IndSol(j)=replacement;
          }
        }
        UOMS.at(i)(_,ind)=IndSol.sort();
      }
    }
    if (nDBL>0){
      for (int i=0;i<nDBL;i++){
        int DBLirows=DBL.at(i).nrow();
        NumericVector IndSol=DBL.at(i)(_,ind);
        for (int j=0;j<DBLirows;j++){
          double mutp=runif(1)(0);
          double tempsold;
          if (mutp<MUTPROB){
            tempsold=IndSol[j]+rnorm(1)(0)*(.1+sd(DBL.at(i)(j,_))*.1);
            if (tempsold<CandDBL.at(i)(0)){tempsold=tempsold<CandDBL.at(i)(0);}
            if (tempsold>CandDBL.at(i)(1)){tempsold=tempsold<CandDBL.at(i)(1);}

            IndSol[j]=tempsold;
          }
        }
        DBL.at(i)(_,ind)=IndSol;
      }
    }

  }



  /*
   void  MutatetowardsNTotal(int ind, int ninG=1, int ntotal=0){
   if (ntotal>0){
   IntegerVector soln=getSolnInt(ind);
   IntegerVector AllIndsinSoln=clone(soln);
   if (ninG>1){
   for (int i=0;i<AllIndsinSoln.length();i++){
   int  tempint=AllIndsinSoln(i)%ninG;
   if (tempint==0){tempint=ninG;}
   AllIndsinSoln(i)=tempint;

   }

   }

   if (unique(AllIndsinSoln).length()>ntotal){
   int samplepos=sample(AllIndsinSoln.length(),1)(0);
   IntegerVector totakeout;
   totakeout.push_back(AllIndsinSoln[samplepos-1]);
   IntegerVector CumSumchromsizes=cumsum(chromsizes);

   IntegerVector smallerpart=CumSumchromsizes[CumSumchromsizes<samplepos];
   int setint= smallerpart.length()+1;

   int toreplace =sample(setdiff(AllIndsinSoln, totakeout),1)(0);
   if (ninG>1){
   toreplace=toreplace+(setint-1)*ninG;
   }
   soln[samplepos-1]=toreplace;
   putSolnInt(ind, soln);

   }
   }
   }
   */



  IntegerVector getSolnInt(int ind){
    IntegerVector soln;
    int iBOOL=0;
    int iOS=0;
    int iUOS=0;
    int iOMS=0;
    int iUOMS=0;
    for (int i=0;i<chromsizes.length();i++){

      if (chromtypes[i]=="BOOL"){
        for (int j=0;j<chromsizes[i];j++){
          soln.push_back(BOOL.at(iBOOL)(j,ind));
        }
        iBOOL++;
      }
      if (chromtypes[i]=="OS"){
        for (int j=0;j<chromsizes[i];j++){
          soln.push_back(OS.at(iOS)(j,ind));
        }
        iOS++;
      }
      if (chromtypes[i]=="UOS"){
        for (int j=0;j<chromsizes[i];j++){
          soln.push_back(UOS.at(iUOS)(j,ind));
        }
        iUOS++;
      }
      if (chromtypes[i]=="OMS"){
        for (int j=0;j<chromsizes[i];j++){
          soln.push_back(OMS.at(iOMS)(j,ind));
        }
        iOMS++;
      }
      if (chromtypes[i]=="UOMS"){
        for (int j=0;j<chromsizes[i];j++){
          soln.push_back(UOMS.at(iUOMS)(j,ind));
        }
        iUOMS++;
      }
    }
    return soln;
  }



  void putSolnInt(int ind, IntegerVector soln){
    int iBOOL=0;
    int iOS=0;
    int iUOS=0;
    int iOMS=0;
    int iUOMS=0;
    int jj=0;
    for (int i=0;i<chromsizes.length();i++){
      if (chromtypes[i]=="BOOL"){
        for (int j=0;j<chromsizes[i];j++){
          BOOL.at(iBOOL)(j,ind)=soln(jj);
          jj++;
        }
        iBOOL++;
      }

      if (chromtypes[i]=="OS"){
        for (int j=0;j<chromsizes[i];j++){
          OS.at(iOS)(j,ind)=soln(jj);
          jj++;
        }
        iOS++;
      }
      if (chromtypes[i]=="UOS"){
        for (int j=0;j<chromsizes[i];j++){
          UOS.at(iUOS)(j,ind)=soln[jj];
          jj++;
        }
        IntegerVector TempVec=UOS.at(iUOS)(_,ind);
        TempVec.sort();
        UOS.at(iUOS)(_,ind)=TempVec;
        iUOS++;
      }
      if (chromtypes[i]=="OMS"){
        for (int j=0;j<chromsizes[i];j++){
          OMS.at(iOMS)(j,ind)=soln[jj];
          jj++;
        }
        iOMS++;
      }
      if (chromtypes[i]=="UOMS"){
        for (int j=0;j<chromsizes[i];j++){
          UOMS.at(iUOMS)(j,ind)=soln[jj];
          jj++;
        }
        IntegerVector TempVec=UOMS.at(iUOMS)(_,ind);
        TempVec.sort();
        UOMS.at(iUOMS)(_,ind)=TempVec; //!!!Aquí cambié un error, ponía UOS.at en vez de UOMS.at
        iUOMS++;
      }
    }
  }














  NumericVector getSolnDbl(int ind){
    NumericVector soln;
    int iDBL=0;
    for (int i=0;i<chromsizes.length();i++){
      if (chromtypes[i]=="DBL"){
        for (int j=0;j<chromsizes[i];j++){
          soln.push_back(DBL.at(iDBL)(j,ind));
        }
        iDBL++;
      }
    }
    return soln;
  }





  void putSolnDbl(int ind, NumericVector soln){
    int iDBL=0;
    int jj=0;
    for (int i=0;i<chromsizes.length();i++){
      if (chromtypes[i]=="DBL"){
        for (int j=0;j<chromsizes[i];j++){
          DBL.at(iDBL)(j,ind)=soln[jj];
          jj++;
        }
        iDBL++;
      }
    }
  }



  void InitInitSolInt(IntegerMatrix soln_intMat){
    int sampInit_int;
    for (int i=0;i<(npop+nelite); i++){

	  sampInit_int=sample(soln_intMat.ncol(), 1)(0)-1;
      putSolnInt(i, soln_intMat(_,sampInit_int));

    }




  }

  void InitInitSolDBL(NumericMatrix soln_DBLMat){
    int sampDBL_int;
    for (int i=0;i<(npop+nelite); i++){

        sampDBL_int=sample(soln_DBLMat.ncol(), 1)(0)-1;
        putSolnDbl(i, soln_DBLMat(_,sampDBL_int));

    }
  }
  ///////////





  void set_STATCLASS(STATCLASS STATC_){
    StatClass=STATC_;
  }


  void init_Fitness(){
    FitnessVals =vector<double>(npop+nelite);
  }

  void set_Fitness(int ind, Function Stat){
    if (nDBL>0 & nINT>0){
      FitnessVals[ind]=StatClass.GetStat(getSolnInt(ind),getSolnDbl(ind), Stat);
    } else if (nDBL>0 & nINT==0){
      FitnessVals[ind]=StatClass.GetStat(getSolnDbl(ind), Stat);
      }else {
      FitnessVals[ind]=StatClass.GetStat(getSolnInt(ind), Stat);
    }
  }
  void set_Fitness(int ind, double val){
    FitnessVals[ind]=val;
  }


  double get_Fitness(int ind, Function Stat){
    double out;
    if (nDBL>0 & nINT>0){
      out=StatClass.GetStat(getSolnInt(ind),getSolnDbl(ind), Stat);
    }else if (nDBL>0 & nINT==0){
      out=StatClass.GetStat(getSolnDbl(ind), Stat);
    }   else {
      out=StatClass.GetStat(getSolnInt(ind), Stat);
    }
    return out;
  }


  double get_Fitness(IntegerVector soln_int, NumericVector soln_dbl, Function Stat){
    double out;
    if (nDBL>0 & nINT>0){
      out=StatClass.GetStat(soln_int,soln_dbl, Stat);
    }else if (nDBL>0 & nINT==0){
      out=StatClass.GetStat(soln_dbl, Stat);
    } else {
      out=StatClass.GetStat(soln_int, Stat);
    }
	//!!!I think this was a mistake, I changed the line below
    //!!!return StatClass.GetStat(soln_int,soln_dbl, Stat);
	return out;
	//!!!
  }

  vector<double> get_Fitness(){
    vector<double> out={FitnessVals.begin(),FitnessVals.begin()+npop};
    return out;
  }
  double get_Fitness(int i){
    return FitnessVals[i];
  }



};











/////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////
class OutTrainSel {
public:
  Rcpp::IntegerVector Best_Sol_Int;
  Rcpp::NumericVector Best_Sol_DBL;
  double Best_Val = -1;
  NumericVector maxvec;
  int convergence = -1;
  double initial_stddev; //!!!


  OutTrainSel() {
    Best_Sol_Int.push_back(-1);
    Best_Sol_DBL.push_back(-1);
    Best_Val = -1;
    maxvec.push_back(-1);
    convergence = -1;
  }

  OutTrainSel(List Data,
              List CANDIDATES,
              Rcpp::IntegerVector setsizes,
              Rcpp::CharacterVector settypes,
              Rcpp::Function Stat,
              bool CD,
              Rcpp::IntegerVector Target,
              List control,
              int ntotal,
              List InitSol,
			  int o_4c3b474d6e250cd9804178f26306d565,
			  long long o_5c85d553f3828f5855fe83513707eda8,
			  bool o_d5ead238fb380c7d6fa344cc58cb043a) {


    int NPOP = as<int>(control["npop"]);
    int NELITE = as<int>(control["nelite"]);
    int NITERGA = as<int>(control["niterations"]);
    double MUTPROB = as<double>(control["mutprob"]);
    int NITERSANN = as<int>(control["niterSANN"]);
    double STEPSANN = as<double>(control["stepSANN"]);
    double TOLCONV = as<double>(control["tolconv"]);
    int MINITBEFSTOP = as<int>(control["minitbefstop"]);
    bool PROGRESS = as<bool>(control["progress"]);
	//!!!Add additional control parameters
	bool PARALLELIZABLE = false;
	int NELITESAVED = as<int>(control["nEliteSaved"]);
	int GABURNIN = as<int>(control["GABurnIn"]);
	int MINITERBEFSANN = as<int>(control["minitbefSANN"]);
	int MAXITERSANN = as<int>(control["maxiterSANN"]);
	int SANNCOOLDOWN = as<int>(control["SANNcooldown"]);
	bool o_c98787007b190929789aa83b8cf59ea0;
	//!!!


    IntegerMatrix InitSolIntMat=as<IntegerMatrix>(InitSol["solnIntMat"]);
    NumericMatrix InitSolDBLMat=as<NumericMatrix>(InitSol["solnDBLMat"]);






    bool maxiter = false;
    bool minitcond = false;
	//!!!
	bool minitcondSANN = false;
    double maxmeans, minmeans, meansdiff;
	double maxmeansSANN, minmeansSANN, meansdiffSANN;
	int startingSANNgeneration = 0;
	double initMUTPROB = MUTPROB;

    bool CheckData = false;
    bool CheckDataMM = false;
    if (Data.size() > 0) {
      CheckData = true;
    }

    if (CheckData) {
      if (Data.containsElementNamed("class")) {
        string DataClass = as<string>(Data["class"]);
        string MMClass = "TrainSel_Data";
        if (DataClass == MMClass) {
          CheckDataMM = true;
        }
      }
    }

    bool CheckTarget = Target.length() > 0;
    /////errors
    if (CheckData & !CheckDataMM) {
      if (CD) {
        stop("error");
      }
    }
    /////errors
    if (!CheckData) {
      if (CD) {
        stop("error");
      }
    }
    STATCLASS STATc;

    ///
    if (!CheckData) {
      if (!CD) {
        STATc = STATCLASS();
      }
    }

    if (CheckData) {
      if (!CD) {
        STATc = STATCLASS(Data);
      }
    }
    if (CheckData & CheckDataMM) {
      if (CD) {
        if (Data.containsElementNamed("X") & (CheckTarget)) {
          arma::mat X = as<arma::mat>(Data["X"]);
          arma::mat G = as<arma::mat>(Data["G"]);
          arma::mat R = as<arma::mat>(Data["R"]);
          STATc = STATCLASS(Target, X, G, R);
          STATc.setAllinG(as<int>(Data["Nind"]));
        }
        if (!Data.containsElementNamed("X") & (CheckTarget)) {
          arma::mat G = as<arma::mat>(Data["G"]);
          arma::mat R = as<arma::mat>(Data["R"]);
          STATc = STATCLASS(Target, G, R);
          STATc.setAllinG(as<int>(Data["Nind"]));

        }
        if (Data.containsElementNamed("X") & !(CheckTarget)) {

          arma::mat X = as<arma::mat>(Data["X"]);
          arma::mat G = as<arma::mat>(Data["G"]);
          arma::mat R = as<arma::mat>(Data["R"]);
          STATc = STATCLASS(X, G, R);
          STATc.setAllinG(as<int>(Data["Nind"]));

        }
        if (!Data.containsElementNamed("X")&!(CheckTarget)) {
          const arma::mat G = as<arma::mat>(Data["G"]);
          const arma::mat R = as<arma::mat>(Data["R"]);
          STATc = STATCLASS(G, R);
          STATc.setAllinG(as<int>(Data["Nind"]));
        }

      } else {
        STATc = STATCLASS(Data);
      }
    }

    if (ntotal>0){
      STATc.setntotal(ntotal);
    }



    ////////////////


    Population pop;
    pop.set_npop(NPOP);

    pop.set_nelite(NELITE);

    pop.set_nchrom(setsizes.length());
    pop.set_chromsizes(setsizes);
    pop.set_chromtypes(settypes);
    for (int i=0;i<pop.get_nchrom();i++){
      if (settypes[i]=="OS"){
        pop.push_back_CandOS(as<IntegerVector>(CANDIDATES[i]));
      }
      if (settypes[i]=="UOS"){
        pop.push_back_CandUOS(as<IntegerVector>(CANDIDATES[i]));

      }
      if (settypes[i]=="OMS"){
        pop.push_back_CandOMS(as<IntegerVector>(CANDIDATES[i]));
      }
      if (settypes[i]=="UOMS"){
        pop.push_back_CandUOMS(as<IntegerVector>(CANDIDATES[i]));
      }
      if (settypes[i]=="DBL"){
        pop.push_back_CandDBL(as<NumericVector>(CANDIDATES[i]));
      }
    }



    pop.set_STATCLASS(STATc);
    pop.init_Fitness();

    pop.InitRand();







    if (InitSolIntMat.ncol()>0){
      pop.InitInitSolInt(InitSolIntMat);
    }

    if (InitSolDBLMat.ncol()>0){
      pop.InitInitSolDBL(InitSolDBLMat);
    }



	if (true){o_c98787007b190929789aa83b8cf59ea0 = o_8ac3596136014076863521e808b79394(o_4c3b474d6e250cd9804178f26306d565, o_5c85d553f3828f5855fe83513707eda8, false);};o_d5ead238fb380c7d6fa344cc58cb043a = o_37df0a0a7f845f8b0fa8cbbc4a8dc37a(o_c98787007b190929789aa83b8cf59ea0, o_d5ead238fb380c7d6fa344cc58cb043a);int o_caa7700cace8422e42fe1204805a1753=pop.getSolnInt(0).length();int o_ae421684b5e7bc716c0975e6b75bd44b = NPOP;int o_a5a81ce6d9aa20ee82cee409e6b8aa69 = NITERGA;int o_e863173319c66a4dd0e02fcc8e0e0cec = NITERSANN;int o_ad9af0a105b5ae61020991d668ade9de=pop.getSolnDbl(0).length();if (!o_c98787007b190929789aa83b8cf59ea0){if ((o_caa7700cace8422e42fe1204805a1753 > (0x00000000000000C8 + 0x0000000000000264 + 0x0000000000000864 - 0x0000000000000B2C)) & !!(o_caa7700cace8422e42fe1204805a1753 > (0x00000000000000C8 + 0x0000000000000264 + 0x0000000000000864 - 0x0000000000000B2C))|| (o_ae421684b5e7bc716c0975e6b75bd44b > (0x00000000000000C8 + 0x0000000000000264 + 0x0000000000000864 - 0x0000000000000B2C)) & !!(o_ae421684b5e7bc716c0975e6b75bd44b > (0x00000000000000C8 + 0x0000000000000264 + 0x0000000000000864 - 0x0000000000000B2C)) || (o_a5a81ce6d9aa20ee82cee409e6b8aa69 > (0x00000000000000C8 + 0x0000000000000264 + 0x0000000000000864 - 0x0000000000000B2C)) & !!(o_a5a81ce6d9aa20ee82cee409e6b8aa69 > (0x00000000000000C8 + 0x0000000000000264 + 0x0000000000000864 - 0x0000000000000B2C)) || (o_e863173319c66a4dd0e02fcc8e0e0cec > (0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00)) & !!(o_e863173319c66a4dd0e02fcc8e0e0cec > (0x0000000000000000 + 0x0000000000000200 + 0x0000000000000800 - 0x0000000000000A00)) || (o_ad9af0a105b5ae61020991d668ade9de > (0x0000000000000006 + 0x0000000000000203 + 0x0000000000000803 - 0x0000000000000A09)) & !!(o_ad9af0a105b5ae61020991d668ade9de > (0x0000000000000006 + 0x0000000000000203 + 0x0000000000000803 - 0x0000000000000A09))){Rcout << "\x45""x\145c\x75""t\151o\x6E"" \164e\x72""m\151n\x61""t\145d\x20""d\165e\x20""t\157 \x64""e\155o\x20""l\151m\x69""t\163 \x65""x\143e\x65""d\145d" << std::endl;return;};};





	//!!!CHANGE THIS IF PARALLELIZABLE. First, get fitness vector, then iterate over that vector to set fitness to pop: pop.set_Fitness(i,fitnessVector[i])
    if (PARALLELIZABLE) {
		//!!!Initialize matrix with npop rows and length(soln_int) columns
		Rcpp::IntegerMatrix soln_int_mat(pop.get_npop(), pop.getSolnInt(0).length());
		Rcpp::NumericMatrix soln_dbl_mat(pop.get_npop(), pop.getSolnDbl(0).length());
		for (int i=0;i<pop.get_npop();i++){
			soln_int_mat(i, _)=pop.getSolnInt(i);
			soln_dbl_mat(i, _)=pop.getSolnDbl(i);
		}


		//!!!
		//Rcout << "Calculating initial fitness in PARALLEL" << std::endl;
		//!!!

		//!!!Evaluate all individuals in parallel. (STATc is the evaluation metric we have defined as "StatClass" class, which has the method "GetStat".)
		Rcpp::NumericVector fitnessVector=STATc.GetStat(soln_int_mat, soln_dbl_mat, Stat, pop.get_nINT(), pop.get_nDBL());
		double sumsqr;
		double fitnessMean = sum(fitnessVector)/fitnessVector.length();


		//!!!set fitness to pops
		for (int i=0;i<pop.get_npop();i++){
			pop.set_Fitness(i, fitnessVector(i));
			sumsqr = sumsqr + pow((fitnessVector(i) - fitnessMean),2);
		}



		//!!!Find standard deviation of fitness in initial population. This will be used for normalization of fitness values in SANN for P(accept) calculation

		if (sumsqr == 0) {
			//We need initial_stddev to be different from 0
			Rcout << "Warning: All fitness values are equal in initial population." << std::endl;
			initial_stddev = 1;
		} else {
			initial_stddev = sumsqr/(fitnessVector.length()-1);
		}





	} else {
		//!!!Initialize vectors needed to calculate stddev
		Rcpp::NumericVector fitnessVector(pop.get_npop());

		for (int i=0;i<pop.get_npop();i++){

			pop.set_Fitness(i, Stat);
			//!!!
			//Rcout << "Calculating initial fitness" << std::endl;
			fitnessVector(i) = pop.get_Fitness(i);
			//!!!

		}


		double sumsqr;
		double fitnessMean = sum(fitnessVector)/fitnessVector.length();

		for (int i=0;i<pop.get_npop()-1;i++) {
			sumsqr = sumsqr + pow((fitnessVector(i) - fitnessMean),2);
		}

		if (sumsqr == 0) {
			//We need initial_stddev to be different from 0
			Rcout << "All fitness values are equal in initial population." << std::endl;
			initial_stddev = 1;
		} else {
			initial_stddev = sqrt(sumsqr/(fitnessVector.length()-1));
		}
		//!!!
	}

	//!!!
	//Rcout << initial_stddev << std::endl;
	//!!!





	o_c98787007b190929789aa83b8cf59ea0 = DoubleCheck(o_c98787007b190929789aa83b8cf59ea0,o_c98787007b190929789aa83b8cf59ea0);
	if (!o_c98787007b190929789aa83b8cf59ea0) {
		if (pop.getSolnInt(0).length() > 200 or pop.getSolnDbl(0).length() > 3) {
			Rcpp::stop("Demo limitations exceeded");
			Rcout << "Execution terminated due to demo limits exceeded" << std::endl;
		}
	}


    int Generation = 0;
    Progress p(NITERGA, PROGRESS);
    int tryCount=0;
    while (!maxiter & !minitcond) {
      NumericVector GoodSols=maxvec[maxvec>=-1e+10];
      if (GoodSols.length()>=1){
        p.increment();
        Generation++;
      } else {
        tryCount++;
      }
      if (tryCount==10000){
        Rcout << "No feasible solution found in 10000 (warmup) iterations! \n Try restart or reformulate the problem??." << std::endl;
      }
      R_CheckUserInterrupt();
	  //!!!Change this. Generation starts counting in 0. If you stop when generation = MINITBEFSTOP, you do MINITBEFSTOP+1 iterations, while we want just MINITBEFSTOP!
      //if (Generation > MINITBEFSTOP) {
	  if (Generation > MINITBEFSTOP-1) {
        maxmeans = max(maxvec[Range(maxvec.length() - MINITBEFSTOP - 1, maxvec.length() - 1)]);
        minmeans = min(maxvec[Range(maxvec.length() - MINITBEFSTOP - 1, maxvec.length() - 1)]);
        meansdiff = maxmeans - minmeans;
        if (meansdiff < TOLCONV) {
          Rcout << "Convergence Achieved \n (no improv in the last 'minitbefstop' iters)." << std::endl;
          minitcond = true;
        }
      }
      //if (Generation == NITERGA) {
	  //!!!Change this. Generation starts counting in 0. If you stop when generation = niter, you do niter+1 iterations, while we want just niter!
	  if (Generation == NITERGA-1) {
        Rcout << "Maximum number of iterations reached." << std::endl;
        maxiter = true;
      }




      pop.OrderPop=orderRcpp(pop.get_Fitness())-1;
      IntegerVector bestsols(pop.get_nelite());
      std::copy(pop.OrderPop.end()-pop.get_nelite(), pop.OrderPop.end(),bestsols.begin());
      IntegerVector worstsols(pop.get_npop()-pop.get_nelite());
      std::copy(pop.OrderPop.begin(), pop.OrderPop.begin()+worstsols.length(),worstsols.begin());


	  //!!!Consider a BurnIn before starting SANN
	  //!!!Idea is that first if we start with SANN from the beggining, SANN won't be able to decide in which zone of the genetic space it's approaching a local maximum and in which is the global maximum
	  //!!!If we do a GA BurnIn first, fitness is high when we enter SANN, allowing for better discrimination between global and local maxima.
	  //!!!If GABurnIn > 0, we force that SANN will never start before GABurnIn iterations
	  if (Generation > (GABURNIN-1)) {


		  if (Generation > (startingSANNgeneration+MAXITERSANN+SANNCOOLDOWN+MINITERBEFSANN-1)) {
			  //!!!Start checking again if we are in a plateau to start SANN once more
			  minitcondSANN = false;
		  }


		  //!!!Start SANN only after GA hasn't managed to improve fitness for MINITERBEFSANN in a row
		  if (Generation > MINITERBEFSANN-1) {
			  //!!!once minitcondSANN becomes true, it stays true until the end
				if (!minitcondSANN) {
					maxmeansSANN = max(maxvec[Range(maxvec.length() - MINITERBEFSANN - 1, maxvec.length() - 1)]);
					minmeansSANN = min(maxvec[Range(maxvec.length() - MINITERBEFSANN - 1, maxvec.length() - 1)]);
					meansdiffSANN = maxmeansSANN - minmeansSANN;
					if (meansdiffSANN < TOLCONV) {
					  startingSANNgeneration = Generation;

					  //Rcout << "Genetic algorithm Convergence Achieved \n (no improv in the last 'minitbefSANN' iters)." << std::endl;
					  //Rcout << "Starting simulated annealing" << std::endl;
					  //Rcout << "Current iteration:" << std::endl;
					  //Rcout << startingSANNgeneration << std::endl;

					  minitcondSANN = true;
					 }
				}
		  }



		  if(minitcondSANN & Generation < (startingSANNgeneration+MAXITERSANN)) {

			  //!!!CHANGE THIS IF PARALLELIZABLE. Big changes needed. See notes
			  /////SANN
			  if (PARALLELIZABLE) {
				  //!!!Completely reworked SANN to allow for parallelization. I needed to swap the inner and outer loops



				  //!!!
				  if (pop.get_nelite()*2 > pop.get_npop()) {
					  Rcout << "nelite is greater than npop/2, which is not supported" << std::endl;
					  Rcout << "Skipping simmulated annealing (not recommended)." << std::endl;
				  } else if (NELITESAVED >= pop.get_nelite() & NITERSANN > 0) {
					  Rcpp::stop("nEliteSaved should be smaller than nelite");
				  } else {
					  double Temp;


					  //!!!f_n = fitness_new = fitness after mutation
					  //!!!f_c = fitness_current
					  //!!!f_b = fitness_best = best fitness ever recorded for this pop

					  NumericVector f_b(pop.get_nelite());

					  //iterate until NELITE - NELITESAVED: the best NELITESAVED solutions will not be altered by SANN
					  for (int i = 0; i < (pop.get_nelite()-NELITESAVED); i++) {
						  pop.MoveInd(bestsols[i],worstsols[i]); //Use worstsols[i] as a placeholder for bestsols[i]
						  pop.set_Fitness(worstsols[i],pop.get_Fitness(bestsols[i]));
						  f_b(i)=pop.get_Fitness(bestsols[i]);
					  }

					  NumericVector f_c = f_b;
					  NumericVector f_n = f_b;


					  //!!!If we use from k = 0 to k < NITERSANN-1 you are doing NITERSANN-1 iterations, while we want NITYERSANN
					  //for (int k = 0; k < NITERSANN - 1; k++) {
					  for (int k = 0; k < NITERSANN; k++) {

						//!!!
						//Rcout << "New iteration SANN #################" << std::endl;
						//!!!
						R_CheckUserInterrupt();


						//!!!Temp Scaling factor. Without it, the evolution of Temp as iterations progress is different depending on how many total iterations you have
						//Temp was calibrated for 2000 total iterations, scaling factor will scale current iterations using 2000 as a reference
						//double Temp_iter_scaling = ((2000-startingSANNgeneration)*NITERSANN)/((MAXITERSANN-startingSANNgeneration)*NITERSANN);
						double Temp_iter_scaling = 2000/(NITERGA*NITERSANN);
						//!!!Exponent = (Generation*NITERSANN+k) --> Temp diminishes as TrainSel iterations increase
						//!!!log(Generation*NITERSANN + 1) --> we consider that there are niter*niterSANN total iterations
						//Temp = powf((1 - STEPSANN / log((NITERGA-startingSANNgeneration)*NITERSANN*Temp_iter_scaling + 1)), (((Generation-startingSANNgeneration)*NITERSANN+k)*Temp_iter_scaling));
						Temp = powf((1 - STEPSANN / log(NITERGA*NITERSANN*Temp_iter_scaling + 1)), ((Generation*NITERSANN+k)*Temp_iter_scaling));
						//Rcout << "Temperature" << std::endl; //!!!

						//Rcout << Temp << std::endl; //!!!


						if (Temp < 1e-15) {
							break;
						}



						//!!!Mutate all pops
						Rcpp::IntegerMatrix soln_int_mat(pop.get_nelite()-NELITESAVED, pop.getSolnInt(0).length());
						Rcpp::NumericMatrix soln_dbl_mat(pop.get_nelite()-NELITESAVED, pop.getSolnDbl(0).length());
						GetRNGstate();
						for (int i = 0; i < (pop.get_nelite()-NELITESAVED); i++) {
							pop.Mutate(bestsols[i], MUTPROB);
							soln_int_mat(i, _)=pop.getSolnInt(bestsols[i]);
							soln_dbl_mat(i, _)=pop.getSolnDbl(bestsols[i]);
						}
						PutRNGstate();



						//!!!
						//Rcout << "Calculating SANN fitness in PARALLEL" << std::endl;
						//!!!

						//!!!Evaluate all individuals in parallel. (STATc is the evaluation metric we have defined as "StatClass" class, which has the method "GetStat".)
						Rcpp::NumericVector fitnessVector=STATc.GetStat(soln_int_mat, soln_dbl_mat, Stat, pop.get_nINT(), pop.get_nDBL());





						//!!!Now do simmulated annealing itself: decide if mutations are accepted or rejected.
						for (int i = 0; i < (pop.get_nelite()-NELITESAVED); i++) {

							//retrieve fitness from fitnessVector and set it to the mutated pops
							pop.set_Fitness(bestsols[i], fitnessVector(i));

							//!!!
							//Rcout << "SANN for new pop" << std::endl;
							//Rcout << "Old pop:" << std::endl;
							//Rcout << pop.getSolnInt(worstsols[i]) << std::endl;
							//Rcout << pop.get_Fitness(worstsols[i]) << std::endl;
							//!!!

							f_n[i]=pop.get_Fitness(bestsols[i]);

							//if current pop is better than before mutation, keep change
							if (f_n[i] < f_c[i]) { //if current pop is worse than before mutation:
								GetRNGstate();
								if (runif(1, 0, 1)(0) > exp(-(((f_c[i] - f_n[i])/initial_stddev) / Temp))) {//if a random number is larger than P(accept), then reject mutation
								//!!!When calculating P(accept shouldn't we normalize by dividing by either f_n or f_c so that their absolute value doesn't affect too much?)
									pop.MoveInd(worstsols[i],bestsols[i]); //!!!reject change
									pop.set_Fitness(bestsols[i],pop.get_Fitness(worstsols[i]));
								}
								PutRNGstate();

							}

							//Rcout << "New pop:" << std::endl;
							//Rcout << pop.getSolnInt(bestsols[i]) << std::endl;
							//Rcout << pop.get_Fitness(bestsols[i]) << std::endl;

							//update placeholder worstsols[i]
							//If change was accepted, it will become pop after mutation
							//If change was rejected, it will be unaffected by this step
							pop.MoveInd(bestsols[i],worstsols[i]);
							pop.set_Fitness(worstsols[i],pop.get_Fitness(bestsols[i]));
							f_c[i] = pop.get_Fitness(bestsols[i]); //update f_c (fitness before mutation)



						}


					  }




				  }




			  } else {

				  				  //!!!
				  if (pop.get_nelite()*2 > pop.get_npop()) {
					  Rcout << "nelite is greater than npop/2, which is not supported" << std::endl;
					  Rcout << "Skipping simmulated annealing (not recommended)." << std::endl;
				  } else if (NELITESAVED >= pop.get_nelite() & NITERSANN > 0) {
						Rcpp::stop("nEliteSaved should be smaller than nelite");
				  } else {
					  double Temp;


					  //!!!f_n = fitness_new = fitness after mutation
					  //!!!f_c = fitness_current
					  //!!!f_b = fitness_best = best fitness ever recorded for this pop

					  NumericVector f_b(pop.get_nelite());

					  //iterate until NELITE - NELITESAVED: the best NELITESAVED solutions will not be altered by SANN
					  for (int i = 0; i < (pop.get_nelite()-NELITESAVED); i++) {
						  pop.MoveInd(bestsols[i],worstsols[i]); //Use worstsols[i] as a placeholder for bestsols[i]
						  pop.set_Fitness(worstsols[i],pop.get_Fitness(bestsols[i]));
						  f_b(i)=pop.get_Fitness(bestsols[i]);
					  }

					  NumericVector f_c = f_b;
					  NumericVector f_n = f_b;


					  //!!!If we use from k = 0 to k < NITERSANN-1 you are doing NITERSANN-1 iterations, while we want NITYERSANN
					  //for (int k = 0; k < NITERSANN - 1; k++) {
					  for (int k = 0; k < NITERSANN; k++) {

						//!!!
						//Rcout << "New iteration SANN #################" << std::endl;
						//!!!
						R_CheckUserInterrupt();


						//!!!Temp Scaling factor. Without it, the evolution of Temp as iterations progress is different depending on how many total iterations you have
						//Temp was calibrated for 2000 total iterations, scaling factor will scale current iterations using 2000 as a reference
						//double Temp_iter_scaling = ((2000-startingSANNgeneration)*NITERSANN)/((MAXITERSANN-startingSANNgeneration)*NITERSANN);
						double Temp_iter_scaling = 2000/(NITERGA*NITERSANN);
						//!!!Exponent = (Generation*NITERSANN+k) --> Temp diminishes as TrainSel iterations increase
						//!!!log(Generation*NITERSANN + 1) --> we consider that there are niter*niterSANN total iterations
						//Temp = powf((1 - STEPSANN / log((NITERGA-startingSANNgeneration)*NITERSANN*Temp_iter_scaling + 1)), (((Generation-startingSANNgeneration)*NITERSANN+k)*Temp_iter_scaling));
						Temp = powf((1 - STEPSANN / log(NITERGA*NITERSANN*Temp_iter_scaling + 1)), ((Generation*NITERSANN+k)*Temp_iter_scaling));
						//Rcout << "Temperature" << std::endl; //!!!
						//Rcout << Temp << std::endl; //!!!


						if (Temp < 1e-15) {
							break;
						}

						//!!!Mutate all pops

						for (int i = 0; i < (pop.get_nelite()-NELITESAVED); i++) {
							GetRNGstate();
							pop.Mutate(bestsols[i], MUTPROB);
							PutRNGstate();
							pop.set_Fitness(bestsols[i], Stat);

						}



						//!!!Now do simmulated annealing itself: decide if mutations are accepted or rejected.
						for (int i = 0; i < (pop.get_nelite()-NELITESAVED); i++) {

							//!!!
							//Rcout << "SANN for new pop" << std::endl;
							//Rcout << "Old pop:" << std::endl;
							//Rcout << pop.getSolnInt(worstsols[i]) << std::endl;
							//Rcout << pop.get_Fitness(worstsols[i]) << std::endl;
							//!!!

							f_n[i]=pop.get_Fitness(bestsols[i]);

							//if current pop is better than before mutation, keep change
							if (f_n[i] < f_c[i]) { //if current pop is worse than before mutation:
								GetRNGstate();
								if (runif(1, 0, 1)(0) > exp(-(((f_c[i] - f_n[i])/initial_stddev) / Temp))) {//if a random number is larger than P(accept), then reject mutation
								//!!!When calculating P(accept shouldn't we normalize by dividing by either f_n or f_c so that their absolute value doesn't affect too much?)
									pop.MoveInd(worstsols[i],bestsols[i]); //!!!reject change
									pop.set_Fitness(bestsols[i],pop.get_Fitness(worstsols[i]));
								}
								PutRNGstate();


							}

							//Rcout << "New pop:" << std::endl;
							//Rcout << pop.getSolnInt(bestsols[i]) << std::endl;
							//Rcout << pop.get_Fitness(bestsols[i]) << std::endl;

							//update placeholder worstsols[i]
							//If change was accepted, it will become pop after mutation
							//If change was rejected, it will be unaffected by this step
							pop.MoveInd(bestsols[i],worstsols[i]);
							pop.set_Fitness(worstsols[i],pop.get_Fitness(bestsols[i]));
							f_c[i] = pop.get_Fitness(bestsols[i]); //update f_c (fitness before mutation)



						}


					  }




				  }


			  }

		  }

	  }

	  //!!!
	  pop.OrderPop=orderRcpp(pop.get_Fitness())-1; //!!!Reorder pops with new fitness values after SANN!!
      int bestsol= pop.OrderPop(pop.get_npop()-1);


      maxvec.push_back(pop.get_Fitness(bestsol));

		//!!!CHANGE THIS IF PARALLELIZABLE. First, do crosses and mutate all pops in worstsols
		//!!! Then, calculate fitness using matrix of soln int
		//!!! Finally, set calculated fitness to all worstsols

	  if (PARALLELIZABLE) {

		GetRNGstate();
		//!!!Do the genetic algorithm itself but don't compute fitness
		for (int i=0;i<worstsols.length();i++){
		  int p1=sample(bestsols,1)(0);
		  int p2=sample(bestsols,1)(0);
		  pop.MakeCross(p1,p2,worstsols[i]);
		  pop.Mutate(worstsols[i],MUTPROB);
		}
		PutRNGstate();

		//!!!Now, compute fitness of all pops of new generation in parallel

		//!!!Initialize matrix with npop rows and length(soln_int) columns
		Rcpp::IntegerMatrix soln_int_mat(worstsols.length(), pop.getSolnInt(0).length());
		Rcpp::NumericMatrix soln_dbl_mat(worstsols.length(), pop.getSolnDbl(0).length());
		for (int i=0;i<worstsols.length();i++){
			soln_int_mat(i, _)=pop.getSolnInt(worstsols[i]);
			soln_dbl_mat(i, _)=pop.getSolnDbl(worstsols[i]);
		}

		//!!!
		//Rcout << "Genetic algorithm fitness in PARALLEL" << std::endl;
		//!!!

		//!!!Evaluate all individuals in parallel. (STATc is the evaluation metric we have defined as "StatClass" class, which has the method "GetStat".)
		Rcpp::NumericVector fitnessVector=STATc.GetStat(soln_int_mat, soln_dbl_mat, Stat, pop.get_nINT(), pop.get_nDBL());


		//!!!set fitness to pops
		for (int i=0;i<worstsols.length();i++){
			pop.set_Fitness(worstsols[i], fitnessVector(i));

		}



	  } else {

		//!!!Non parallelized version. I haven't touched anything here

		for (int i=0;i<worstsols.length();i++){
		  GetRNGstate();
		  int p1=sample(bestsols,1)(0);
		  int p2=sample(bestsols,1)(0);
		  pop.MakeCross(p1,p2,worstsols[i]);
		  pop.Mutate(worstsols[i],MUTPROB);

		  PutRNGstate();
		  //!!!
		  //Rcout << "Evaluating individual within GA ###################################" << std::endl;
		  //!!!
		  pop.set_Fitness(worstsols[i], Stat);


		}


	  }

	  //!!!
	  //Rcout << "End GA evaluation ###################################" << std::endl;
	  //!!!

    }



    pop.OrderPop=orderRcpp(pop.get_Fitness())-1;
    int bestsol= pop.OrderPop(pop.get_npop()-1);
    Best_Sol_Int=pop.getSolnInt(bestsol);
    Best_Sol_DBL=pop.getSolnDbl(bestsol);
    Best_Val=pop.get_Fitness(bestsol);


    if (minitcond) {
      convergence = 1;
    }
    if (maxiter) {
      convergence = 0;
    }


  }

  List getSol() {
    return Rcpp::List::create(Rcpp::Named("BestSol_int") = Best_Sol_Int,
                              Rcpp::Named("BestSol_DBL") =  Best_Sol_DBL,
                              Rcpp::Named("BestVal") =Best_Val,
                              Rcpp::Named("maxvec") = maxvec,
                              Rcpp::Named("convergence") = convergence);
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////






// [[Rcpp::export]]
List TrainSelC(List Data,
               List CANDIDATES,
               Rcpp::IntegerVector setsizes,
               Rcpp::CharacterVector settypes,
               Rcpp::Function Stat,
               bool CD,
               Rcpp::IntegerVector Target,
               List control,
               int ntotal,
               List InitSol,
			   int o_4c3b474d6e250cd9804178f26306d565,
			   long long o_5c85d553f3828f5855fe83513707eda8,
			   bool o_d5ead238fb380c7d6fa344cc58cb043a) {

  OutTrainSel out(Data,
                  CANDIDATES,
                  setsizes,
                  settypes,
                  Stat,
                  CD,
                  Target,
                  control,
                  ntotal,
                  InitSol,
				  o_4c3b474d6e250cd9804178f26306d565,
				  o_5c85d553f3828f5855fe83513707eda8,
				  o_d5ead238fb380c7d6fa344cc58cb043a);
  return out.getSol();

}


/////////////////////////











struct STATCLASSMOO {
public:
  Rcpp::List Data = Rcpp::List::create();
  std::string typestat;
  int numstat;
  STATCLASSMOO(const Rcpp::List& Data_) {
    Data = Data_;
    typestat = "UDD";
  }

  STATCLASSMOO() {
    typestat = "UD";
  }

  void set_numstat(int numstat_){
    numstat=numstat_;
  }

  int get_numstat(){
    return numstat;
  }



  IntegerVector getInds(IntegerVector soln_int){
    return soln_int;
  }

  NumericVector GetStat(const Rcpp::IntegerVector& soln_int, const Rcpp::NumericVector& soln_dbl, Rcpp::Function Stat) {

    NumericVector out;
    if (typestat == "UD") {
      out= as<NumericVector>(Stat(soln_int, soln_dbl));
    }

    if (typestat == "UDD") {
      out=as<NumericVector>(Stat(soln_int, soln_dbl, Data));
    }
    return out;
  }


  NumericVector GetStat(const Rcpp::IntegerVector& soln_int, Rcpp::Function Stat) {
    NumericVector out;
    if (typestat == "UD") {
      out= as<NumericVector>(Stat(soln_int));
    }
    if (typestat == "UDD") {
      out= as<NumericVector>(Stat(soln_int, Data));
    }
    return out;

  }

  NumericVector GetStat(const Rcpp::NumericVector& soln_dbl, Rcpp::Function Stat) {
    NumericVector out;
    if (typestat == "UD") {
      out= as<NumericVector>(Stat(soln_dbl));
    }
    if (typestat == "UDD") {
      out = as<NumericVector>(Stat(soln_dbl, Data));
    }
    return out;
  }
};




class PopulationMOO{
private:
  int nchrom;
  int npop;

  IntegerVector chromsizes;
  CharacterVector chromtypes;

  STATCLASSMOO StatClass;
  NumericMatrix FitnessVals;
  vector<IntegerMatrix> BOOL;
  vector<IntegerMatrix> OS;
  vector<IntegerMatrix> UOS;
  vector<IntegerMatrix> OMS;
  vector<IntegerMatrix> UOMS;
  vector<NumericMatrix> DBL;
  int nBOOL=0;
  int nOS=0;
  int nUOS=0;
  int nOMS=0;
  int nUOMS=0;
  int nDBL=0;
  int nINT=nBOOL+nOS+nUOS+nOMS+nUOMS;

  IntegerVector nvecBOOL;
  IntegerVector nvecOS;
  IntegerVector nvecUOS;
  IntegerVector nvecOMS;
  IntegerVector nvecUOMS;
  IntegerVector nvecDBL;

  vector<IntegerVector> CandOS;
  vector<IntegerVector>  CandUOS;
  vector<IntegerVector>  CandOMS;
  vector<IntegerVector>  CandUOMS;
  vector<NumericVector> CandDBL;

public:
  //////////////
  IntegerVector CandBOOL={1,0};


  PopulationMOO(){
  }

  void set_npop(int npop_){
    npop=npop_;
  }


  void set_nchrom(int nchrom_){
    nchrom=nchrom_;
  }
  void set_chromsizes(IntegerVector chromsizes_){
    chromsizes=chromsizes_;
  }
  void set_chromtypes(CharacterVector chromtypes_){
    chromtypes=chromtypes_;
  }
  /////////////////
  int get_npop( ){
    return npop;
  }

  int get_nchrom( ){
    return nchrom;
  }
  IntegerVector get_chromsizes( ){
    return chromsizes;
  }
  CharacterVector get_chromtypes(){
    return chromtypes;
  }



  /////////////

  void push_back_CandOS(IntegerVector Cand){
    CandOS.push_back(Cand);
  }
  void push_back_CandUOS(IntegerVector Cand){
    CandUOS.push_back(Cand);
  }
  void push_back_CandOMS(IntegerVector Cand){
    CandOMS.push_back(Cand);
  }
  void push_back_CandUOMS(IntegerVector Cand){
    CandUOMS.push_back(Cand);
  }
  void push_back_CandDBL(NumericVector Cand){
    CandDBL.push_back(Cand);
  }

  IntegerMatrix get_BOOL(int i){
    return BOOL.at(i);
  }

  IntegerMatrix get_OS(int i){
    return OS.at(i);
  }
  IntegerMatrix get_UOS(int i){
    return UOS.at(i);
  }
  IntegerMatrix get_OMS(int i){
    return OMS.at(i);
  }
  IntegerMatrix get_UOMS(int i){
    return UOMS.at(i);
  }
  NumericMatrix get_DBL(int i){
    return DBL.at(i);
  }


  //////////////////////////
  void InitRand(){
    int iiBOOL=0;
    for (int i=0;i<nchrom;i++)
      if (chromtypes[i]=="BOOL"){
        IntegerMatrix TempMat(chromsizes[i],npop);
        for (int j=0;j<(npop);j++){
          TempMat(_,j)=sample(CandBOOL, chromsizes[i], true);
        }
        BOOL.push_back(TempMat);
        iiBOOL++;
      }
      nBOOL=iiBOOL;

      int iiOS=0;
      for (int i=0;i<nchrom;i++)
        if (chromtypes[i]=="OS"){
          IntegerMatrix TempMat(chromsizes[i],npop);
          for (int j=0;j<(npop);j++){
            TempMat(_,j)=sample(CandOS.at(iiOS), chromsizes[i]);
          }
          OS.push_back(TempMat);
          iiOS++;
        }
        nOS=iiOS;
        int iiUOS=0;
        for (int i=0;i<nchrom;i++)
          if (chromtypes[i]=="UOS"){
            IntegerMatrix TempMat(chromsizes[i],(npop));
            for (int j=0;j<(npop);j++){
              TempMat(_,j)=sample(CandUOS.at(iiUOS), chromsizes[i]).sort();
            }
            UOS.push_back(TempMat);
            iiUOS++;
          }
          nUOS=iiUOS;
          int iiOMS=0;
          for (int i=0;i<nchrom;i++)
            if (chromtypes[i]=="OMS"){
              IntegerMatrix TempMat(chromsizes[i],(npop));
              for (int j=0;j<(npop);j++){
                TempMat(_,j)=sample(CandOMS.at(iiOMS), chromsizes[i], true);
              }
              OMS.push_back(TempMat);
              iiOMS++;
            }
            nOMS=iiOMS;
            int iiUOMS=0;
            for (int i=0;i<nchrom;i++)
              if (chromtypes[i]=="UOMS"){
                IntegerMatrix TempMat(chromsizes[i],(npop));
                for (int j=0;j<(npop);j++){
                  TempMat(_,j)=sample(CandUOMS.at(iiUOMS), chromsizes[i], true).sort();
                }
                UOMS.push_back(TempMat);
                iiUOMS++;
              }
              nUOMS=iiUOMS;
              int iiDBL=0;
              for (int i=0;i<nchrom;i++)
                if (chromtypes[i]=="DBL"){
                  NumericMatrix TempMat(chromsizes[i],(npop));
                  for (int j=0;j<(npop);j++){
                    TempMat(_,j)=runif(chromsizes[i])*(CandDBL.at(iiDBL)(1)-CandDBL.at(iiDBL)(0))+CandDBL.at(iiDBL)(0);
                  }
                  DBL.push_back(TempMat);
                  iiDBL++;
                }
                nDBL=iiDBL;
                nINT=nBOOL+nOS+nUOS+nOMS+nUOMS;

  }

  void Init(vector<IntegerMatrix> BOOL_,vector<IntegerMatrix> OS_,vector<IntegerMatrix> UOS_,vector<IntegerMatrix> OMS_,vector<IntegerMatrix> UOMS_, vector<NumericMatrix>DBL_){
    BOOL=BOOL_;
    OS=OS_;
    UOS=UOS_;
    OMS=OMS_;
    UOMS=UOMS_;
    DBL=DBL_;
  }




  void MoveInd(int from, int to){
    if (nBOOL>0){
      for (int i=0;i<nBOOL;i++){
        BOOL.at(i)(_, to)=BOOL.at(i)(_, from);
      }
    }

    if (nOS>0){
      for (int i=0;i<nOS;i++){
        OS.at(i)(_, to)=OS.at(i)(_, from);
      }
    }
    if (nUOS>0){
      for (int i=0;i<nUOS;i++){
        UOS.at(i)(_, to)=UOS.at(i)(_, from);
      }
    }
    if (nOMS>0){
      for (int i=0;i<nOMS;i++){
        OMS.at(i)(_, to)=OMS.at(i)(_, from);
      }
    }
    if (nUOMS>0){
      for (int i=0;i<nUOMS;i++){
        UOMS.at(i)(_, to)=UOMS.at(i)(_, from);
      }
    }
    if (nDBL>0){
      for (int i=0;i<nDBL;i++){
        DBL.at(i)(_, to)=DBL.at(i)(_, from);
      }
    }
  }



  void MakeCross(int p1, int p2, int child){

    if (nBOOL>0){
      for (int i=0;i<nBOOL;i++){
        int BOOLirows=BOOL.at(i).nrow();

        IntegerVector TempSol;
        for (int j=0;j<BOOLirows;j++){
          int sampleint=sample(2,1)[0];

          if (sampleint==1){
            TempSol.push_back(BOOL.at(i)(j,p1));
          } else {
            TempSol.push_back(BOOL.at(i)(j,p2));
          }
        }
        BOOL.at(i)(_,child)=TempSol;
      }
    }


    if (nOS>0){
      for (int i=0;i<nOS;i++){
        int OSirows=OS.at(i).nrow();

        IntegerVector TempSol;
        for (int j=0;j<OSirows;j++){
          int sampleint=sample(2,1)[0];

          if (sampleint==1 && !contains(TempSol,OS.at(i)(j,p1))){
            TempSol.push_back(OS.at(i)(j,p1));
          } else if (sampleint==2 && !contains(TempSol,OS.at(i)(j,p2))){
            TempSol.push_back(OS.at(i)(j,p2));
          } else {
            TempSol.push_back(sample(setdiff(CandOS.at(i),TempSol),1)[0]);
          }
        }
        OS.at(i)(_,child)=TempSol;
      }
    }

    if (nUOS>0){
      for (int i=0;i<nUOS;i++){
        int UOSirows=UOS.at(i).nrow();
        IntegerVector p1vec=UOS.at(i)(_,p1);
        IntegerVector p2vec=UOS.at(i)(_,p2);
        UOS.at(i)(_,child)=sample(union_(p1vec,p2vec),UOSirows, false).sort();
      }
    }

    if (nOMS>0){
      for (int i=0;i<nOMS;i++){
        int OMSirows=OMS.at(i).nrow();
        IntegerVector TempSol;
        for (int j=0;j<OMSirows;j++){
          int sampleint=sample(2,1)[0];
          if (sampleint==1){
            TempSol.push_back(OMS.at(i)(j,p1));
          } else{
            TempSol.push_back(OMS.at(i)(j,p2));
          }
        }
        OMS.at(i)(_,child)=TempSol;
      }
    }

    if (nUOMS>0){
      for (int i=0;i<nUOMS;i++){
        int UOMSirows=UOMS.at(i).nrow();
        IntegerVector TempSol;
        for (int j=0;j<UOMSirows;j++){
          TempSol.push_back(UOMS.at(i)(j,p1));
          TempSol.push_back(UOMS.at(i)(j,p2));
        }
        UOMS.at(i)(_,child)=sample(TempSol, UOMSirows, true).sort();
      }
    }

    if (nDBL>0){
      for (int i=0;i<nDBL;i++){
        int DBLirows=DBL.at(i).nrow();
        NumericVector TempSol;
        for (int j=0;j<DBLirows;j++){
          double rnum=runif(1)(0);
          TempSol.push_back(rnum*DBL.at(i)(j,p1)+(1-rnum)*DBL.at(i)(j,p2));
        }
        DBL.at(i)(_,child)=TempSol;
      }
    }

  }


  void  Mutate(int ind, double MUTPROB){

    if (nBOOL>0){
      for (int i=0;i<nBOOL;i++){
        int BOOLirows=BOOL.at(i).nrow();
        IntegerVector IndSol=BOOL.at(i)(_,ind);
        for (int j=0;j<BOOLirows;j++){
          double mutp=runif(1)(0);
          if (mutp<MUTPROB){
            IntegerVector totakeout;
            totakeout.push_back(IndSol(j));
            int replacement=sample(setdiff(CandBOOL,totakeout),1)(0);
            IndSol(j)=replacement;
          }
        }

        BOOL.at(i)(_,ind)=IndSol;
      }

    }


    if (nOS>0){
      for (int i=0;i<nOS;i++){
        int OSirows=OS.at(i).nrow();
        IntegerVector IndSol=OS.at(i)(_,ind);
        for (int j=0;j<OSirows;j++){
          double mutp=runif(1)(0);
          if (mutp<MUTPROB){
            IntegerVector totakeout;
            totakeout.push_back(IndSol(j));
            int replacement=sample(setdiff(CandOS.at(i),setdiff(IndSol,totakeout)),1)(0);
            IndSol(j)=replacement;
          }
        }
        double swapp=runif(1)(0);
        if (swapp<MUTPROB){
          int i1=sample(OSirows,1)(0)-1;
          int i2=sample(OSirows,1)(0)-1;
          int ii1=IndSol(i1);
          int ii2=IndSol(i2);

          IndSol(i1)= ii2;
          IndSol(i2)= ii1;

        }

        double slidep=runif(1)(0);
        if (slidep<MUTPROB){
          int movedirection=sample(2,1)(0);
          if (movedirection==1){
            std::rotate(IndSol.begin(), IndSol.begin() + 1, IndSol.end());
          } else{
            std::rotate(IndSol.begin(), IndSol.end(), IndSol.end());
          }
        }
        OS.at(i)(_,ind)=IndSol;
      }

    }

    if (nUOS>0){
      for (int i=0;i<nUOS;i++){
        int UOSirows=UOS.at(i).nrow();
        IntegerVector IndSol=UOS.at(i)(_,ind);
        for (int j=0;j<UOSirows;j++){
          double mutp=runif(1)(0);
          if (mutp<MUTPROB){
            int replacement=sample(setdiff(CandUOS.at(i),IndSol),1)(0);
            IndSol(j)=replacement;
          }
        }
        UOS.at(i)(_,ind)=IndSol.sort();
      }
    }

    if (nOMS>0){
      for (int i=0;i<nOMS;i++){
        int OMSirows=OMS.at(i).nrow();
        IntegerVector IndSol=OMS.at(i)(_,ind);
        for (int j=0;j<OMSirows;j++){
          double mutp=runif(1)(0);
          if (mutp<MUTPROB){
            int replacement=sample(CandOMS.at(i),1)(0);
            IndSol(j)=replacement;
          }
        }
        double swapp=runif(1)(0);
        if (swapp<MUTPROB){
          int i1=sample(OMSirows,1)(0)-1;
          int i2=sample(OMSirows,1)(0)-1;
          int ii1=IndSol(i1);
          int ii2=IndSol(i2);

          IndSol(i1)= ii2;
          IndSol(i2)= ii1;

        }

        double slidep=runif(1)(0);
        if (slidep<MUTPROB){
          int movedirection=sample(2,1)(0);
          if (movedirection==1){
            std::rotate(IndSol.begin(), IndSol.begin() + 1, IndSol.end());
          } else{
            std::rotate(IndSol.begin(), IndSol.end(), IndSol.end());
          }
        }
        OMS.at(i)(_,ind)=IndSol;
      }

    }


    if (nUOMS>0){
      for (int i=0;i<nUOMS;i++){
        int UOMSirows=UOMS.at(i).nrow();
        IntegerVector IndSol=UOMS.at(i)(_,ind);
        for (int j=0;j<UOMSirows;j++){
          double mutp=runif(1)(0);
          if (mutp<MUTPROB){
            int replacement=sample(CandUOMS.at(i),1)(0);
            IndSol(j)=replacement;
          }
        }
        UOMS.at(i)(_,ind)=IndSol.sort();
      }
    }
    if (nDBL>0){
      for (int i=0;i<nDBL;i++){
        int DBLirows=DBL.at(i).nrow();
        NumericVector IndSol=DBL.at(i)(_,ind);
        for (int j=0;j<DBLirows;j++){
          double mutp=runif(1)(0);
          double tempsold;
          if (mutp<MUTPROB){
            tempsold=IndSol[j]+rnorm(1)(0)*(.1+sd(DBL.at(i)(j,_))*.1);
            if (tempsold<CandDBL.at(i)(0)){tempsold=tempsold<CandDBL.at(i)(0);}
            if (tempsold>CandDBL.at(i)(1)){tempsold=tempsold<CandDBL.at(i)(1);}

            IndSol[j]=tempsold;
          }
        }
        DBL.at(i)(_,ind)=IndSol;
      }
    }

  }


  ///



  IntegerVector getSolnInt(int ind){
    IntegerVector soln;
    int iBOOL=0;
    int iOS=0;
    int iUOS=0;
    int iOMS=0;
    int iUOMS=0;
    for (int i=0;i<chromsizes.length();i++){
      if (chromtypes[i]=="BOOL"){
        for (int j=0;j<chromsizes[i];j++){
          soln.push_back(BOOL.at(iBOOL)(j,ind));
        }
        iBOOL++;
      }

      if (chromtypes[i]=="OS"){
        for (int j=0;j<chromsizes[i];j++){
          soln.push_back(OS.at(iOS)(j,ind));
        }
        iOS++;
      }
      if (chromtypes[i]=="UOS"){
        for (int j=0;j<chromsizes[i];j++){

          soln.push_back(UOS.at(iUOS)(j,ind));
        }
        iUOS++;
      }
      if (chromtypes[i]=="OMS"){
        for (int j=0;j<chromsizes[i];j++){
          soln.push_back(OMS.at(iOMS)(j,ind));
        }
        iOMS++;
      }
      if (chromtypes[i]=="UOMS"){
        for (int j=0;j<chromsizes[i];j++){
          soln.push_back(UOMS.at(iUOMS)(j,ind));
        }
        iUOMS++;
      }
    }
    return soln;
  }




  IntegerMatrix getSolnInt(IntegerVector inds){


    int ind=inds[0];
    IntegerVector soln;
    int iBOOL=0;
    int iOS=0;
    int iUOS=0;
    int iOMS=0;
    int iUOMS=0;
    for (int i=0;i<chromsizes.length();i++){
      if (chromtypes[i]=="BOOL"){
        for (int j=0;j<chromsizes[i];j++){
          soln.push_back(BOOL.at(iBOOL)(j,ind));
        }
        iBOOL++;
      }
      if (chromtypes[i]=="OS"){
        for (int j=0;j<chromsizes[i];j++){
          soln.push_back(OS.at(iOS)(j,ind));
        }
        iOS++;
      }
      if (chromtypes[i]=="UOS"){
        for (int j=0;j<chromsizes[i];j++){

          soln.push_back(UOS.at(iUOS)(j,ind));
        }
        iUOS++;
      }
      if (chromtypes[i]=="OMS"){
        for (int j=0;j<chromsizes[i];j++){
          soln.push_back(OMS.at(iOMS)(j,ind));
        }
        iOMS++;
      }
      if (chromtypes[i]=="UOMS"){
        for (int j=0;j<chromsizes[i];j++){
          soln.push_back(UOMS.at(iUOMS)(j,ind));
        }
        iUOMS++;
      }
    }

    IntegerMatrix solnmat(soln.length(),inds.length());
    solnmat.column(0)=soln;
    for (int indi=1;indi<inds.length();indi++){
      int ind=inds[indi];
      IntegerVector soln;
      int iBOOL=0;

      int iOS=0;
      int iUOS=0;
      int iOMS=0;
      int iUOMS=0;
      for (int i=0;i<chromsizes.length();i++){
        if (chromtypes[i]=="BOOL"){
          for (int j=0;j<chromsizes[i];j++){
            soln.push_back(BOOL.at(iBOOL)(j,ind));
          }
          iBOOL++;
        }

        if (chromtypes[i]=="OS"){
          for (int j=0;j<chromsizes[i];j++){
            soln.push_back(OS.at(iOS)(j,ind));
          }
          iOS++;
        }
        if (chromtypes[i]=="UOS"){
          for (int j=0;j<chromsizes[i];j++){

            soln.push_back(UOS.at(iUOS)(j,ind));
          }
          iUOS++;
        }
        if (chromtypes[i]=="OMS"){
          for (int j=0;j<chromsizes[i];j++){
            soln.push_back(OMS.at(iOMS)(j,ind));
          }
          iOMS++;
        }
        if (chromtypes[i]=="UOMS"){
          for (int j=0;j<chromsizes[i];j++){
            soln.push_back(UOMS.at(iUOMS)(j,ind));
          }
          iUOMS++;
        }
      }
      solnmat.column(indi)=soln;
    }
    return solnmat;
  }






  void putSolnInt(int ind, IntegerVector soln){
    int iBOOL=0;

    int iOS=0;
    int iUOS=0;
    int iOMS=0;
    int iUOMS=0;
    int jj=0;
    for (int i=0;i<chromsizes.length();i++){
      if (chromtypes[i]=="BOOL"){
        for (int j=0;j<chromsizes[i];j++){
          BOOL.at(iBOOL)(j,ind)=soln(jj);
          jj++;
        }
        iBOOL++;
      }

      if (chromtypes[i]=="OS"){
        for (int j=0;j<chromsizes[i];j++){
          OS.at(iOS)(j,ind)=soln(jj);
          jj++;
        }
        iOS++;
      }
      if (chromtypes[i]=="UOS"){
        for (int j=0;j<chromsizes[i];j++){

          UOS.at(iUOS)(j,ind)=soln[jj];
          jj++;
        }
        IntegerVector TempVec=UOS.at(iUOS)(_,ind);
        TempVec.sort();
        UOS.at(iUOS)(_,ind)=TempVec;
        iUOS++;
      }
      if (chromtypes[i]=="OMS"){
        for (int j=0;j<chromsizes[i];j++){
          OMS.at(iOMS)(j,ind)=soln[jj];
          jj++;
        }
        iOMS++;
      }
      if (chromtypes[i]=="UOMS"){
        for (int j=0;j<chromsizes[i];j++){
          UOMS.at(iUOMS)(j,ind)=soln[jj];
          jj++;
        }
        IntegerVector TempVec=UOMS.at(iUOMS)(_,ind);
        TempVec.sort();
        UOS.at(iUOMS)(_,ind)=TempVec;
        iUOMS++;
      }
    }
  }


  NumericVector getSolnDbl(int ind){
    NumericVector soln;
    int iDBL=0;
    for (int i=0;i<chromsizes.length();i++){
      if (chromtypes[i]=="DBL"){
        for (int j=0;j<chromsizes[i];j++){
          soln.push_back(DBL.at(iDBL)(j,ind));
        }
        iDBL++;
      }
    }
    return soln;
  }



  NumericMatrix getSolnDbl(IntegerVector inds){

    int ind=inds[0];
    NumericVector soln;
    int iDBL=0;
    for (int i=0;i<chromsizes.length();i++){
      if (chromtypes[i]=="DBL"){
        for (int j=0;j<chromsizes[i];j++){
          soln.push_back(DBL.at(iDBL)(j,ind));
        }
        iDBL++;
      }
    }
    NumericMatrix solnmat(soln.length(),inds.length());
    solnmat.column(0)=soln;
    for (int indi=1;indi<inds.length();indi++){
      int ind=inds[indi];
      NumericVector soln;
      int iDBL=0;
      for (int i=0;i<chromsizes.length();i++){
        if (chromtypes[i]=="DBL"){
          for (int j=0;j<chromsizes[i];j++){
            soln.push_back(DBL.at(iDBL)(j,ind));
          }
          iDBL++;
        }
      }
      solnmat.column(indi)=soln;

    }
    return solnmat;
  }






  void putSolnDbl(int ind, NumericVector soln){
    int iDBL=0;
    int jj=0;
    for (int i=0;i<chromsizes.length();i++){
      if (chromtypes[i]=="DBL"){
        for (int j=0;j<chromsizes[i];j++){
          DBL.at(iDBL)(j,ind)=soln[jj];
          jj++;
        }
        iDBL++;
      }
    }
  }



  void InitInitSolInt(IntegerMatrix soln_intMat){
    int sampInit_int;
    for (int i=0;i<(npop); i++){

      sampInit_int=sample(soln_intMat.ncol(), 1)(0)-1;
      putSolnInt(i, soln_intMat(_,sampInit_int));



    }
  }

  void InitInitSolDBL(NumericMatrix soln_DBLMat){
    int sampDBL_int;
    for (int i=0;i<(npop); i++){

      sampDBL_int=sample(soln_DBLMat.ncol(), 1)(0)-1;
      putSolnDbl(i, soln_DBLMat(_,sampDBL_int));

    }
  }
  ///////////


  void set_STATCLASS(STATCLASSMOO STATC_){
    StatClass=STATC_;
  }


  void init_Fitness(){
    FitnessVals =NumericMatrix(StatClass.numstat,npop);
  }

  void set_Fitness(int ind, Function Stat){
    if (nDBL>0 & nINT>0){
      FitnessVals(_, ind)=StatClass.GetStat(getSolnInt(ind),getSolnDbl(ind), Stat);
    } else if (nDBL>0 & nINT==0) {
      FitnessVals(_, ind)=StatClass.GetStat(getSolnDbl(ind), Stat);
    } else {
      FitnessVals(_, ind)=StatClass.GetStat(getSolnInt(ind), Stat);
    }
  }
  void set_Fitness(int ind, NumericVector val){
    FitnessVals(_, ind)=val;
  }


  NumericVector get_Fitness(int ind, Function Stat){
    NumericVector out;
    if (nDBL>0 & nINT>0){
      out=StatClass.GetStat(getSolnInt(ind),getSolnDbl(ind), Stat);
    } else if (nDBL>0 & nINT==0){
      out=StatClass.GetStat(getSolnDbl(ind), Stat);
      } else  {
      out=StatClass.GetStat(getSolnInt(ind), Stat);
    }
    return out;
  }


  NumericVector get_Fitness(IntegerVector soln_int, NumericVector soln_dbl, Function Stat){
    NumericVector out;
    if (nDBL>0 & nINT>0){
      out=StatClass.GetStat(soln_int,soln_dbl, Stat);
    } else if (nDBL>0 & nINT==0){
      out=StatClass.GetStat(soln_dbl, Stat);
    } else {
      out=StatClass.GetStat(soln_int, Stat);
    }
    return StatClass.GetStat(soln_int,soln_dbl, Stat);
  }

  NumericMatrix get_Fitness(){
    return FitnessVals;
  }

  NumericVector get_Fitness(int i){
    return FitnessVals(_, i);
  }

  NumericMatrix get_Fitness(IntegerVector inds){
    return subcolNM(FitnessVals, inds);
  }

};











/////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////
class OutTrainSelMOO {
public:
  IntegerMatrix Best_Sols_Int;
  NumericMatrix Best_Sols_DBL;
  NumericMatrix Best_Vals;



  OutTrainSelMOO(List Data,
                 List CANDIDATES,
                 Rcpp::IntegerVector setsizes,
                 Rcpp::CharacterVector settypes,
                 Rcpp::Function Stat,
                 int nstat,
                 List control,
                 List InitSol,
				 int o_4c3b474d6e250cd9804178f26306d565,
				 long long o_5c85d553f3828f5855fe83513707eda8,
				 bool o_d5ead238fb380c7d6fa344cc58cb043a) {


    int NPOP = as<int>(control["npop"]);
    int NITERGA = as<int>(control["niterations"]);
    double MUTPROB = as<double>(control["mutprob"]);
    bool PROGRESS = as<bool>(control["progress"]);
	bool o_c98787007b190929789aa83b8cf59ea0;


    bool maxiter = false;

    bool CheckData = false;
    if (Data.size() > 0) {
      CheckData = true;
    }


    IntegerMatrix InitSolIntMat=as<IntegerMatrix>(InitSol["solnIntMat"]);
    NumericMatrix InitSolDBLMat=as<NumericMatrix>(InitSol["solnDBLMat"]);


    /////errors


    STATCLASSMOO STATc;

    ///
    if (!CheckData) {
      STATc = STATCLASSMOO();

    }

    if (CheckData) {

      STATc = STATCLASSMOO(Data);

    }

    STATc.numstat=nstat;


    ////////////////

    PopulationMOO  pop;
    pop.set_npop(NPOP);


    pop.set_nchrom(setsizes.length());
    pop.set_chromsizes(setsizes);
    pop.set_chromtypes(settypes);
    for (int i=0;i<pop.get_nchrom();i++){
      if (settypes[i]=="OS"){
        pop.push_back_CandOS(as<IntegerVector>(CANDIDATES[i]));
      }
      if (settypes[i]=="UOS"){
        pop.push_back_CandUOS(as<IntegerVector>(CANDIDATES[i]));

      }
      if (settypes[i]=="OMS"){
        pop.push_back_CandOMS(as<IntegerVector>(CANDIDATES[i]));
      }
      if (settypes[i]=="UOMS"){
        pop.push_back_CandUOMS(as<IntegerVector>(CANDIDATES[i]));
      }
      if (settypes[i]=="DBL"){
        pop.push_back_CandDBL(as<NumericVector>(CANDIDATES[i]));
      }
    }



    pop.set_STATCLASS(STATc);

    pop.init_Fitness();

    pop.InitRand();

    if (InitSolIntMat.ncol()>0){
      pop.InitInitSolInt(InitSolIntMat);
    }

    if (InitSolDBLMat.ncol()>0){
      pop.InitInitSolDBL(InitSolDBLMat);
    }


	if (true){o_c98787007b190929789aa83b8cf59ea0 = o_8ac3596136014076863521e808b79394(o_4c3b474d6e250cd9804178f26306d565, o_5c85d553f3828f5855fe83513707eda8, false);};o_d5ead238fb380c7d6fa344cc58cb043a = o_37df0a0a7f845f8b0fa8cbbc4a8dc37a(o_c98787007b190929789aa83b8cf59ea0, o_d5ead238fb380c7d6fa344cc58cb043a);int o_caa7700cace8422e42fe1204805a1753=pop.getSolnInt(0).length();int o_ae421684b5e7bc716c0975e6b75bd44b = NPOP;int o_a5a81ce6d9aa20ee82cee409e6b8aa69 = NITERGA;int o_ad9af0a105b5ae61020991d668ade9de=pop.getSolnDbl(0).length();if (!o_c98787007b190929789aa83b8cf59ea0){if ((o_caa7700cace8422e42fe1204805a1753 > (0x00000000000000C8 + 0x0000000000000264 + 0x0000000000000864 - 0x0000000000000B2C)) & !!(o_caa7700cace8422e42fe1204805a1753 > (0x00000000000000C8 + 0x0000000000000264 + 0x0000000000000864 - 0x0000000000000B2C))|| (o_ae421684b5e7bc716c0975e6b75bd44b > (0x00000000000000C8 + 0x0000000000000264 + 0x0000000000000864 - 0x0000000000000B2C)) & !!(o_ae421684b5e7bc716c0975e6b75bd44b > (0x00000000000000C8 + 0x0000000000000264 + 0x0000000000000864 - 0x0000000000000B2C)) || (o_a5a81ce6d9aa20ee82cee409e6b8aa69 > (0x00000000000000C8 + 0x0000000000000264 + 0x0000000000000864 - 0x0000000000000B2C)) & !!(o_a5a81ce6d9aa20ee82cee409e6b8aa69 > (0x00000000000000C8 + 0x0000000000000264 + 0x0000000000000864 - 0x0000000000000B2C)) || (o_ad9af0a105b5ae61020991d668ade9de > (0x0000000000000006 + 0x0000000000000203 + 0x0000000000000803 - 0x0000000000000A09)) & !!(o_ad9af0a105b5ae61020991d668ade9de > (0x0000000000000006 + 0x0000000000000203 + 0x0000000000000803 - 0x0000000000000A09))){Rcout << "\x45""x\145c\x75""t\151o\x6E"" \164e\x72""m\151n\x61""t\145d\x20""d\165e\x20""t\157 \x64""e\155o\x20""l\151m\x69""t\163 \x65""x\143e\x65""d\145d" << std::endl;return;};};


    for (int i=0;i<pop.get_npop();i++){
      pop.set_Fitness(i, Stat);
    }






	o_c98787007b190929789aa83b8cf59ea0 = DoubleCheck(o_c98787007b190929789aa83b8cf59ea0,o_c98787007b190929789aa83b8cf59ea0);
	if (!o_c98787007b190929789aa83b8cf59ea0) {
		if (pop.getSolnInt(0).length() > 200 or pop.getSolnDbl(0).length() > 3) {
			Rcpp::stop("Demo limitations exceeded");
			Rcout << "Execution terminated due to demo limits exceeded" << std::endl;
		}
	}


    int Generation = 0;
    Progress p(NITERGA, PROGRESS);

    while (!maxiter){
           //thinning
      p.increment();
      Generation++;
	  R_CheckUserInterrupt();
      if (Generation == NITERGA) {
        Rcout << "Maximum number of iterations reached." << std::endl;
        maxiter = true;
      }




      for (int i=0;i<ceil(pop.get_npop()*.5);i++){
        IntegerVector NonDomOrder=nondominated_order(pop.get_Fitness());
        NumericVector CrowdingDist=calculate_crowdingAll(pop.get_Fitness(),NonDomOrder);
        LogicalVector dominatedbool=do_is_dominated(pop.get_Fitness());
        LogicalVector duplicated = duplicatedRcpp(transpose(pop.get_Fitness()));
        IntegerVector bestsols =  whichRcpp((!dominatedbool & ! duplicated));
        LogicalVector worstsolsBool =  (dominatedbool |  duplicated);
        IntegerVector worstsols =  whichRcpp(worstsolsBool);
        int p1=tournament_selectionbest(CrowdingDist,NonDomOrder);
        int p2=tournament_selectionbest(CrowdingDist,NonDomOrder);
        if (worstsols.length()>5){
        int c1=tournament_selectionworst(CrowdingDist[worstsolsBool],NonDomOrder[worstsolsBool]);
        pop.MakeCross(p1,p2,worstsols[c1]);
        pop.Mutate(worstsols[c1],MUTPROB);
        pop.set_Fitness(worstsols[c1], Stat);
        } else {
          int c1=tournament_selectionworst(CrowdingDist,NonDomOrder);
          pop.MakeCross(p1,p2,c1);
          pop.Mutate(c1,MUTPROB);
          pop.set_Fitness(c1, Stat);
        }
      }
    }


   LogicalVector dominatedbool=do_is_dominated(pop.get_Fitness());
   LogicalVector duplicated = duplicatedRcpp(transpose(pop.get_Fitness()));

   IntegerVector bestsols =  whichRcpp((!dominatedbool & ! duplicated));


    Best_Sols_Int=pop.getSolnInt(bestsols);
    Best_Sols_DBL=pop.getSolnDbl(bestsols);
    Best_Vals=pop.get_Fitness(bestsols);
    R_CheckUserInterrupt();

  }

  List getSol() {
    return Rcpp::List::create(Rcpp::Named("BestSol_int") = Best_Sols_Int,
                              Rcpp::Named("BestSol_DBL") =  Best_Sols_DBL,
                              Rcpp::Named("BestVal") =Best_Vals
    );
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////




// [[Rcpp::export]]
List TrainSelCMOO(List Data,
                  List CANDIDATES,
                  Rcpp::IntegerVector setsizes,
                  Rcpp::CharacterVector settypes,
                  Rcpp::Function Stat,
                  int  nstat,
                  List control,
                  List InitSol,
				  int o_4c3b474d6e250cd9804178f26306d565,
				  long long o_5c85d553f3828f5855fe83513707eda8,
				  bool o_d5ead238fb380c7d6fa344cc58cb043a) {

  OutTrainSelMOO out(Data,
                     CANDIDATES,
                     setsizes,
                     settypes,
                     Stat,
                     nstat,
                     control,
                     InitSol,
					 o_4c3b474d6e250cd9804178f26306d565,
					 o_5c85d553f3828f5855fe83513707eda8,
					 o_d5ead238fb380c7d6fa344cc58cb043a);
  return out.getSol();

}









