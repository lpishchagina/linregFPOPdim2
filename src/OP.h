#ifndef OP_H
#define OP_H
#include "Geom0.h"
#include "Geom3.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <iterator>
#include <string> 

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

/*
 Class OP
 -------------------------------------------------------------------------------
 Description: 
 Template for the realization of FPOP-Algorithm in 2-dimension. 
 
 Parameters:
 "n" - data length;
 "penalty" - value of penalty;
 "Sum" - matrix(n+1x2*p) of sum:x1:xp, x1^2:xp^2;
 "chpts" - vector of changepoints;
 "kCoef" - k of changepoints; 
 "aCoef" - a of changepoints;  
 "globalCost" - value of global cost.
 --------------------------------------------------------------------------------
 */

template <class GeomX>
class OP{
private:
  double penalty; //value of penalty 
  unsigned int n; //data length
  double** Sum;  // vector sum xy,x,y, x^2, y^2
  std::vector< unsigned int> chpts;    //changepoints vector 
  std::vector<double> kCoef;         //kCoef
  std::vector<double> aCoef;         //aCoef        
  double globalCost;                  //value of global cost
public:
  OP<GeomX>(){};
  //----------------------------------------------------------------------------
  OP<GeomX> (std::vector<double>& x, std::vector<double>&y, double beta){
    n = x.size();
    penalty = beta;
    Sum = new double*[n+1]; 
    for(unsigned int i = 0; i < n + 1; i++) {Sum[i] = new double[5];}
  }
  //----------------------------------------------------------------------------
  ~OP<GeomX>(){
    for(unsigned int i = 0; i < n+1; i++) {delete(Sum[i]);}
    delete [] Sum;
    Sum = NULL;
  }
  //----------------------------------------------------------------------------
  std::vector< unsigned int > get_chpts() const {return chpts;}
  std::vector< double > get_kCoef() const{return kCoef;}
  std::vector< double > get_aCoef() const{return aCoef;}
  double get_globalCost() const {return globalCost;}
  unsigned int get_n() const {return n;}
  double get_penalty() const {return penalty;}
  //----------------------------------------------------------------------------
  double** vect_Sum(std::vector<double>& x, std::vector<double>& y) {
    //Sum[0]=xy,Sum[1]=x,Sum[2]=y,Sum[3]=x^2,Sum[4]=y^2
    for (unsigned i = 0; i < 5; i++){ Sum[0][i] = 0;}
    for (unsigned int j = 1; j < (n + 1); j++){
      Sum[j][0] = Sum[j - 1][0] + x[j - 1]*y[j - 1];
      Sum[j][1] = Sum[j - 1][1] + x[j - 1];
      Sum[j][2] = Sum[j - 1][2] + y[j - 1];
      Sum[j][3] = Sum[j - 1][3] + x[j - 1] * x[j - 1];
      Sum[j][4] = Sum[j - 1][4] + y[j - 1] * y[j - 1];
    }
    return(Sum);
  }
  //----------------------------------------------------------------------------
  void algoFPOP(std::vector<double>& x, std::vector<double>& y, int type, bool test_mode){
    //preprocessing-------------------------------------------------------------
    Sum = vect_Sum(x, y); 
    double* m = new double[n + 1];        // "globalCost" = m[n] - chpts.size()*penalty
    m[0] = 0;  
    double** Chpt_k_a = new double*[3];// vectors of best last changepoints, k and a
    for(unsigned int i = 0; i < 3; i++) {Chpt_k_a[i] = new double[n];}
    
    std::ofstream test_file;
    if (test_mode == true){test_file.open("test.txt");}
    
    //Algorithm-----------------------------------------------------------------
    double min_val;                                 //if qit(a,k)> mt+penalty =>PELT
    double kTemp;
    double aTemp;
    double mdif;
    unsigned int lbl; 
    unsigned int u; 
    lrCost cost;
    GeomX geom = GeomX(0);
    std::list<GeomX> list_geom;    //list of geometry
    std::list<Ellps> list_Ellps;//list of active ellipses(t-1)
    for (unsigned int t = 0; t < n ; t++){
      cost = lrCost(t, t, Sum[t], Sum[t+1], m[t]);
      min_val = cost.get_min();              
      kTemp =  cost.get_k();   
      aTemp = cost.get_a(); 
      lbl = t;
      list_Ellps.clear();

      //First run: searching min------------------------------------------------
      typename std::list<GeomX>::reverse_iterator rit_geom = list_geom.rbegin();
      while(rit_geom!= list_geom.rend()){
        u = rit_geom -> get_label_t(); 
        // Searching: min
        cost = lrCost(u, t, Sum[u], Sum[t + 1], m[u]);
        if( min_val >= cost.get_min()){
          lbl = u;
          min_val = cost.get_min();
          kTemp = cost.get_k();
          aTemp = cost.get_a();  
        }
        //list of active Ellpss(t-1)
        cost = lrCost(u, t-1, Sum[u], Sum[t], m[u]);
        list_Ellps.push_back(Ellps(cost));
        ++rit_geom;
      }
      //best last changepoints and means
      Chpt_k_a[0][t] = lbl;       //vector of best last chpt
      Chpt_k_a[1][t] = kTemp;     //vector of k
      Chpt_k_a[2][t] = aTemp;     //vector of a
      //new min 
      m[t + 1] = min_val + penalty;
      
      //Initialisation of geometry----------------------------------------------
      geom.InitialGeometry(t, list_Ellps);
      list_geom.push_back(geom);
      
      //Second run: Update list of geometry-------------------------------------
      typename std::list<GeomX>::iterator it_geom = list_geom.begin();
      while (it_geom != list_geom.end()){
        lbl = it_geom -> get_label_t();
        cost = lrCost(lbl, t, Sum[lbl], Sum[t + 1], m[lbl]);
        mdif = m[t + 1] - m[lbl] - cost.get_min(); //if qit(a,k)> mt+penalty =>PELT
        //PELT
        if (mdif <= 0){it_geom = list_geom.erase(it_geom); --it_geom;}
        //FPOP
        if (mdif > 0){
          it_geom -> UpdateGeometry(Ellps(cost));
          if (it_geom -> EmptyGeometry()){it_geom = list_geom.erase(it_geom);--it_geom;}
          else {if (test_mode == true && (type == 2 || type == 3)){ test_file << it_geom ->get_label_t() << " "<< it_geom ->get_ellps().size() << " ";}}
        }//else
        ++it_geom;
      }
      if (test_mode == true){test_file << "\n";} 
    }
    if (test_mode == true){test_file.close();}
    //Result vectors------------------------------------------------------------
    unsigned int chp = n;
    while (chp > 0){
      chpts.push_back(chp);
      kCoef.push_back(Chpt_k_a[1][chp-1]);
      aCoef.push_back(Chpt_k_a[2][chp-1]);
      chp = Chpt_k_a[0][chp-1];
    }
    reverse(chpts.begin(), chpts.end());
    chpts.pop_back();                
    reverse(kCoef.begin(), kCoef.end());
    reverse(aCoef.begin(), aCoef.end());
    globalCost = m[n] - penalty * chpts.size();  
    //memory--------------------------------------------------------------------
    for(unsigned int i = 0; i < 3; i++) {delete(Chpt_k_a[i]);}
    delete [] Chpt_k_a;
    Chpt_k_a = NULL;
    delete [] m;
    m = NULL;
  }
  //----------------------------------------------------------------------------
};

#endif //OP_H      
    