#ifndef GEOM3_H
#define GEOM3_H

#include <iostream>
#include <vector>
#include <list>
#include <iterator>

#include "Ellps.h"
#include "lrCost.h"

//Class Geom3
//------------------------------------------------------------------------------
//Description of geometry "Geom3": 
//Geometry for FPOP-Algorithm in 2-dimension. 
//Parameters of geometry:"label_t" - moment of time, "ellpst1" - list of ellipses(t-1)
//The updated geometry is a ellipse that approximates (ellipse at the moment t) minus (list of ellpst1) .
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
    
    double Dist(double a1, double a2, double b1, double b2);
    
    void InitialGeometry(unsigned int i, const std::list<Ellps> &ellpses);
    void UpdateGeometry(Ellps &elt);
    bool EmptyGeometry();
};
#endif //GEOM3_H