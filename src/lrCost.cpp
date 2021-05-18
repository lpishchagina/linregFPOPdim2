#include "lrCost.h"

#include <iostream>
#include "math.h"

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//constructor*******************************************************************
lrCost::lrCost(unsigned int i, unsigned int t,double* Si1, double* St, double mi_1pen){
  mi1p = mi_1pen;   
  cnst = t - i + 1;
  
  /*comment---------------------------------------------------------------------
  Notations: Sj[0]=sum_j(x_jy_j),Sj[1]=sum_j(x_j),Sj[2]=sum_j(y_j),Sj[3]=sum_j(x_j^2),S_j[4]=sum(y_j^2);
  
  Estimations:
    k =  [(t-i+1)sum_{j = i}^t (x_j*y_j) - sum_{j = i}^t (x_j)*sum_{j = i}^t(y_j)]/[(t-i+1)*sum_{j = i}^t (x_j^2) - (sum_{j = i}^t x_j)^2];
    a = [sum_{j = i}^t (y_j^2)  - k* sum_{j = i}^t (x_j)]/(t-i+1);
  ----------------------------------------------------------------------------*/
  k = (cnst* (St[0] - Si1[0]) - (St[2] - Si1[2])*(St[1] - Si1[1]))/(cnst*(St[3] - Si1[3]) - (St[1] - Si1[1])*(St[1] - Si1[1]));
  a = ((St[2] - Si1[2]) - k*(St[1] - Si1[1]))/cnst; 
  
  /*comment---------------------------------------------------------------------
    q_i^t(k, a) =  k*k*sum_{j=i}^t(x_j^2) + 2*a*k*sum_{j=i}^t(x_j) + (t-i+1)*a*a - 2*k*sum_{j=i}^t(x_j*y_j) - 2*a*sum_{j=i}^t(y_j) + sum_{j=i}^t(y_j^2) + mi1p;
  
    value_min = q_i^t(k_est, a_est)=  k*k*(St[3]-Si1[3]) + 2*a*k*(St[1]-Si1[1]) + cnst*a*a + 2*k*(Si1[0]-St[0]) + 2*a*(Si1[2]-St[2]) + (St[4]-Si1[4]) + mi1p;
  
    value_min = A*k^2 + 2*B*k*a + C*a^2 + 2*D*k + 2*E*a + F + mi1p;
  ----------------------------------------------------------------------------*/
  //Parameters
  A = St[3] - Si1[3]; 
  B = St[1] - Si1[1];
  C = cnst;
  D = Si1[0] - St[0]; 
  E = Si1[2] - St[2];
  F = St[4] - Si1[4];
  //value_min
  if (i == t){ value_min = INFINITY;} else { value_min = A*k*k + 2*a*k*B + C*a*a + 2*k*D + 2*a*E + F + mi1p;}
}
//accessory*********************************************************************
unsigned int lrCost::g_cnst()const{return cnst;}
double lrCost::g_mi1p()const{return mi1p;}

double lrCost::g_k()const{return k;}
double lrCost::g_a()const{return a;}

double lrCost::g_A()const{return A;}
double lrCost::g_B()const{return B;}
double lrCost::g_C()const{return C;}
double lrCost::g_D()const{return D;}
double lrCost::g_E()const{return E;}
double lrCost::g_F()const{return F;}

double lrCost::g_min()const{return value_min;}
//******************************************************************************