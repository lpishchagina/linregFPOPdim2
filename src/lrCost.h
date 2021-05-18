#include <Rcpp.h>
using namespace Rcpp;

#ifndef LRCOST_H
#define LRCOST_H

#include <vector>
/*
 Class lrCost
 -------------------------------------------------------------------------------
 Description: 
 Cost function of the time series (x_{i:t},y_{i:t}) when y_j = k*x_j+a +ej, j=i,t.  
 a,k  - regression coefficients.
 -------------------------------------------------------------------------------
 */
class lrCost{
private:
  double mi1p;  //sum of the value of the optimal cost at moment (i-1) and penalty;
  unsigned int cnst;  // (t - i + 1)
  // Fun(a,k) = 0, Fun(k,a) = A*k^2+2*B*k*a+C*a^2+2*D*k+2*E*a+F = 0 
  //Parameters
  double A;
  double B;
  double C;
  double D;
  double E;
  double F;
  //Estimations
  double k; 
  double a;
  //value_min = q_i^t(a,k)= mi1p + sum_{j = i}^t (yj -(k*x_j + a))^2;
  double value_min;
public:
  lrCost(){};
  lrCost(unsigned int i, unsigned int t, double* Si1, double* St, double mi_1pen);//S = vector of sum xy,x,y,x2,y2
  //accessory
  double g_A() const;
  double g_B() const;
  double g_C() const;
  double g_D() const;
  double g_E() const;
  double g_F() const;
  
  unsigned int g_cnst() const;
  double g_mi1p() const;
  
  double g_k() const;
  double g_a() const;
  
  double g_min() const;
};

#endif // LRCOST_H