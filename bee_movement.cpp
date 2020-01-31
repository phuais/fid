#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <ctime>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
List bee_movement(NumericMatrix bees, arma::mat mosaic, std::vector<int> ibees, std::vector<int> target, int n_sampled, float first_approach, float t) {

  srand(time(NULL));
  for(int i = 0; i < bees.nrow(); i++){
    
    int b_row = bees(i, 0);
    int b_col = bees(i, 1);
    // Subset mosaic to the surrounding of the individual bee position
    arma::mat m = mosaic.submat(b_row-1, b_col-1, b_row+1, b_col+1);
    int s = arma::accu(m);
    
    NumericVector m_pos = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };  
    NumericVector probs;
    
    for(int ii = 0; ii < 3; ii++){
      for(int jj = 0; jj < 3; jj++)
        probs.push_back(m(jj, ii)/s);
    }  
    
   // Randomly or not randomly (i.e. with certain probabilities) choose a new cell to move
   // movements ids:
   // 1  4  7
   // 2  5  8
   // 3  6  9
    NumericVector tomove = RcppArmadillo::sample(m_pos, 1, 0, probs);
    float move = tomove(0);

    // new individual bee row position
    bees(i, 0) = bees(i, 0) - (2 - (move - 3*(ceil(move/3)-1)));
    // new individual bee col position
    bees(i, 1) = bees(i, 1) - (2 - ceil(move/3));
    
    std::vector<int> bee_position; 
    for(int w = 0; w < 2; w++)
      bee_position.push_back(bees(i, w));

    // When a bee hits the target, n_sampled increase a unit
    if((std::find(ibees.begin(), ibees.end(), i) == ibees.end()) && (bee_position == target)){
      n_sampled = n_sampled + 1;
      ibees.push_back(i);
    }
    
    // Saves the time (i.e. t) of the first interaction
    if(std::isnan(first_approach)){
      if(bee_position == target){
        first_approach = t;
      }
    }
}

  return List::create(bees, n_sampled, first_approach, ibees);
}