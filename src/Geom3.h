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
    
    double dsX_to_x(double dsx, double dsy, double dx, double dy, double angle);
    double dsY_to_y(double dsx, double dsy, double dx, double dy, double angle);
    
    double x_to_dsX(double x, double y, double dx, double dy, double angle);
    double y_to_dsY(double x, double y, double dx, double dy, double angle);
    
    bool pnt_insd_E(double x, double y, const Ellps E);
    
    double Dist(double a1, double a2, double b1, double b2);
    
    int filtKalman(double c1, double c2, const Mat2X2 &A, double d1, double d2, const Mat2X2 &B);
   
    void InitialGeometry(unsigned int i, const std::list<Ellps> &ellpses);
    void UpdateGeometry(const Ellps &Et);
    bool EmptyGeometry();
};
#endif //GEOM3_H