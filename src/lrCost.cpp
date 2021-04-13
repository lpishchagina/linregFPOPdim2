#include "lrCost.h"

#include <iostream>
#include "math.h"

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
//The method of least squares
//k =  [(t-i+1)sum_{j = i}^t (x_j*y_j) - sum_{j = i}^t (x_j)*sum_{j = i}^t(y_j)]/[(t-i+1)*sum_{j = i}^t (x_j^2) - (sum_{j = i}^t x_j)^2]
//a = [sum_{j = i}^t (y_j^2)  - b* sum_{j = i}^t (x_j)]/(t-i+1)
//value_min - cost function minimum

// qmin = q(a,k)

//S[0]=xy,S[1]=x,S[2]=y,S[3]=x^2,S[4]=y^2

//constructor*******************************************************************
lrCost::lrCost(unsigned int i, unsigned int t,double* Si1, double* St, double mi_1pen){
  cnst = t - i + 1;
  mi1p = mi_1pen;   //S[0]=xy,S[1]=x,S[2]=y,S[3]=x^2,S[4]=y^2
  k = (cnst* (St[0] - Si1[0]) - (St[2] - Si1[2])*(St[1] - Si1[1]))/(cnst*(St[3] - Si1[3]) - (St[1] - Si1[1])*(St[1] - Si1[1]));
  a = ((St[2] - Si1[2]) - k*(St[1] - Si1[1]))/cnst; 
  
  value_min = mi1p + cnst*a*a - 2*k*(St[0] - Si1[0]) +2*a*k*(St[1] - Si1[1]) - 2*a*(St[2] - Si1[2]) + k*k*(St[3] - Si1[3]) + (St[4] - Si1[4]);
}
//accessory*********************************************************************
unsigned int lrCost::get_cnst()const{return cnst;}

double lrCost::get_mi1p()const{return mi1p;}

double lrCost::get_k()const{return k;}

double lrCost::get_a()const{return a;}

double lrCost::get_min()const{return value_min;}

//******************************************************************************