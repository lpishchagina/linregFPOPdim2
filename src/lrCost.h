#include <Rcpp.h>
using namespace Rcpp;

#ifndef LRCOST_H
#define LRCOST_H

#include <vector>
/*
 Class lrCost
 -------------------------------------------------------------------------------
 Description: 
 y = kx+a +e
 The cost function for the interval (i,t) in 2-dimension for time series (x,y)
 
 Parameters:
 "cnst" = (t - i + 1);
 "mi1p" - sum of the value of the optimal cost at moment (i-1) and penalty;
  a,b  - regression coefficients.
 -------------------------------------------------------------------------------
 */

class lrCost{
private:
  unsigned int cnst;
  double mi1p; 
  double k;
  double a;
  double value_min;
  
public:
  lrCost(){};
  lrCost(unsigned int i, unsigned int t, double* Si1, double* St, double mi_1pen);//S = vector of sum xy,x,y,x2,y2
  
  unsigned int get_cnst() const;
  double get_mi1p() const;
  double get_k() const;
  double get_a() const;
  double get_min() const;
};

#endif // LRCOST_H