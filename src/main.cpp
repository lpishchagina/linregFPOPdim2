#include "OP.h"
#include "Rect.h"
#include "lrCost.h"
#include "Ellipse.h"
#include "Geom0.h"
#include "Geom1.h"

#include "math.h"

#include <Rcpp.h>
//using namespace Rcpp;
//using namespace std;

//' @title FPOP2D 
//'                                                                                                        
//' @description Detecting changepoints using the functional pruning optimal partitioning method (fpop) in linear regression y = kCoef*x + aCoef.                         
//'                                                                                                       
//' @param x.                                
//' @param y.                                
//' @param penalty is a value of penalty (a non-negative real number).                                        
//' @param type is a value defining the  type of geometry for FPOP-pruning: type=1: ("intersection" of sets), approximation - rectangle.       
//'                                                                                                          
//' @return a list of 4 elements  = (chpts, kCoef, aCoef, globalCost).                    
//'  
//' \describe{
//' \item{\code{chpts}}{is the vector of changepoints.}
//' \item{\code{kCoef}}{is the vector of k.}
//' \item{\code{aCoef}}{is the vector of a.}
//' \item{\code{globalCost}}{is a number equal to the global cost.}
//' }                                                                                                                                                                             #     
//'             
//' @examples 
// [[Rcpp::export]]
List lrFPOP2D(std::vector<double> x, std::vector<double> y, double penalty, int type)
{
  //----------stop--------------------------------------------------------------
  if(x.size() != y.size()){throw std::range_error("x and y have different length");}
  if(penalty < 0) {throw std::range_error("penalty should be a non-negative number");}
  if(type < 0 || type > 3)
  {throw std::range_error("type must be one of: 0,1");}
  //----------------------------------------------------------------------------
  List res;
  bool test;
  test = false;
  if (type == 0)
  {
    //test = true;//
    OP<Geom0> N = OP<Geom0>(x, y, penalty);
    N.algoFPOP(x, y, type, test);  
    res["chpts"] = N.get_chpts();
    res["kCoef"] = N.get_kCoef();
    res["aCoef"] = N.get_aCoef();
    res["globalCost"] = N.get_globalCost();
  }
  if (type == 1)
  {
    //test = true;//
    OP<Geom1> N = OP<Geom1>(x, y, penalty);
    N.algoFPOP(x, y, type, test);  
    res["chpts"] = N.get_chpts();
    res["kCoef"] = N.get_kCoef();
    res["aCoef"] = N.get_aCoef();
    res["globalCost"] = N.get_globalCost();
  }
  return res;
}
