#ifndef GEOM3_H
#define GEOM3_H

#include <iostream>
#include <vector>
#include <list>
#include <iterator>

#include "Ellps.h"
#include "lrCost.h"
#include "Mat2x2.h"

//Class Geom3
//------------------------------------------------------------------------------
//Description of geometry "Geom3": 
//Geometry for FPOP-Algorithm in 2-dimension. 
//Parameters of geometry:"label_t" - moment of time, "ellps" - list of ellipses(t-1)
//The updated geometry is a ellipse that approximates (ellipse at the moment t) minus (list of ellps) .
//Check for emptiness - Points extrem inside Ellipse_t. 
//------------------------------------------------------------------------------
class Geom3
{
  private:
    unsigned int label_t;       //time moment 
    std::list<Ellps> ellps;  //list of ellipses(t-1)
    bool fl_empty;
    
  public:
    Geom3();
    Geom3(unsigned int t);
    
    unsigned int get_label_t() const;
    std::list<Ellps> get_ellps() const;
    
    unsigned int InterKalman(double c1, double c2, const double* &A, double d1, double d2, const double* &B);
    unsigned int ExclKalman(double c1, double c2, const double* &A, double d1, double d2, const double* &B);
    unsigned int seqShturm(double l, double C0, double C1, double C2, double C3);
    
//    double* eigenValues(const Mat2X2 &A, const Mat2X2 &B);
   
   
    void InitialGeometry(unsigned int i, const std::list<Ellps> &ellpses);
    void UpdateGeometry(const Ellps &Et);
    bool EmptyGeometry();
};
#endif //GEOM3_H