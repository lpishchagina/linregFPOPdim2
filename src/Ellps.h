#ifndef ELLPS_H
#define ELLPS_H

#include "lrCost.h"
#include <vector>
/*
 Class Ellps
 -------------------------------------------------------------------------------
 Description: 
 Ellps  that corresponds of Fun(k,a) = A*k^2+2*B*k*a+C*a^2+2*D*k+2*E*a+F = 0
 (X,Y) =>(x,y)
 (x/a)^2 +(y/b)^2 = c^2
 (A*k^2+2*B*k*a+C*a^2+2*D*k+2*E*a+F=0) => (displacement, slope) => (x/a)^2 +(y/b)^2 = c^2
 --------------------------------------------------------------------------------
 */
class Ellps
{
private:
  //displacement (x0,y0)=>(0,0)  
  double d1, d2;
  //roots lmbd1>=lmbd2
  double lmbd1,lmbd2;
  //slope
  double k1, k2;
  //angle
  double angl;
  //coefficients
  double a, b, c2;
  //focus F1 (-F, 0), F2(F, 0)
  double F;
public:
  Ellps(){};
  Ellps(lrCost Q);  
  //accessory
  double get_d1() const;
  double get_d2() const;
 
  double get_a() const;
  double get_b() const;
  double get_c2() const;
  
  double get_F() const;
  
  double get_lmbd1() const;
  double get_lmbd2() const;
  
  double get_k1() const;
  double get_k2() const;
  double get_angl() const;
  
  bool insd_pnt(double x, double y);
  bool insd_ellps(Ellps &El);
  
  double dst_pnts(double x1, double y1, double x2, double y2);
  double dst_ellps_pnt(double x1, double y1);
 
  double x_XY(double X, double Y);
  double y_XY(double X, double Y);
    
  double X_xy(double x, double y);
  double Y_xy(double x, double y);
  
  double get_r(double x, double y);
  
}; 

#endif //ELLPS_H
//------------------------------------------------------------------------------