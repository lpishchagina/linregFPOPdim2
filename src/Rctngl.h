#ifndef RCTNGL_H
#define RCTNGL_H

#include "math.h"
#include "Ellps.h"

#include <vector>


/*
 Class Rctngl
 -------------------------------------------------------------------------------
 Description: 
 Rectangle in 2-dimension. 
 
 Parameters:
 -------------------------------------------------------------------------------
 */
class Rctngl{
private:
  double x1, y1;
  double x2, y2;                     
  double x3, y3;
  double x4, y4;
public:
  
  Rctngl():x1(-INFINITY), y1(-INFINITY), x2(-INFINITY), y2(-INFINITY), x3(-INFINITY), y3(-INFINITY),x4(-INFINITY), y4(-INFINITY){}
  Rctngl(double xx1, double yy1, double xx2, double yy2, double xx3, double yy3, double xx4, double yy4):x1(xx1), y1(yy1), x2(xx2), y2(yy2), x3(xx3), y3(yy3),x4(xx4), y4(yy4){}
  
  double get_x1() const;
  double get_y1() const;
  double get_x2() const;
  double get_y2() const;
  double get_x3() const;
  double get_y3() const;
  double get_x4() const;
  double get_y4() const;
  
  double min_ab(double a, double b);
  double max_ab(double a, double b);
  
  double dsX_to_x(double dsx, double dsy, double dx, double dy, double angle);
  double dsY_to_y(double dsx, double dsy, double dx, double dy, double angle);
  
  double x_to_dsX(double x, double y, double dx, double dy, double angle);
  double y_to_dsY(double x, double y, double dx, double dy, double angle);
  
  double dst_pnts(double a1, double b1, double a2, double b2);
  
  bool insd_pnt(double x, double y);
  bool pnt_insd_trngl(double x, double y, double a1, double b1, double a2, double b2, double c1, double c2);
  double get_r(double x, double y, const Ellps &E);
  
  bool EmptyIntersection(const Ellps &E);
  
  bool IsEmpty_Rctngl();
  void Exclusion_disk(const Ellps &E){};
  void Intersection_disk(const Ellps &E){};
  
  
};

#endif //RCTNGL_H