#ifndef ELLPS_H
#define ELLPS_H

#include "lrCost.h"
#include "Mat2x2.h"

#include <vector>
/*
 Class Ellps
 -------------------------------------------------------------------------------
 Description: 
 Ellps  that corresponds of Fun(k,a) = A*k^2+2*B*k*a+C*a^2+2*D*k+2*E*a+F = 0
 (X,Y) =>(x,y)
 (x/a)^2 +(y/b)^2 = c^2
 (A*k^2+2*B*k*a+C*a^2+2*D*k+2*E*a+F=0) => (displacement, slope) => (x/a)^2 +(y/b)^2 = c^2
 -------------------------------------------------------------------------------
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
  double a, b;
  //focus Fs1 (-Fs, 0), Fs2(Fs, 0)
  double Fs;
  //Matrix (A,B//B,C)
  Mat2X2 Sigma;
  
public:
  Ellps(){};
  Ellps(lrCost Q, double r2);
  //accessory
  double get_d1() const;
  double get_d2() const;
 
  double get_a() const;
  double get_b() const;
  
  double get_Fs() const;
  
  double get_lmbd1() const;
  double get_lmbd2() const;
  
  double get_k1() const;
  double get_k2() const;
  double get_angl() const;
  
  Mat2X2 get_Sigma() const;
  
  bool insd_pnt(double x, double y);
  double dst_pnts(double x1, double y1, double x2, double y2);
 
  double x_XY(double X, double Y);
  double y_XY(double X, double Y);
    
  double X_xy(double x, double y);
  double Y_xy(double x, double y);
  
  double get_r(double x, double y);
  
  int testInter(const Ellps &E);
}; 

#endif //ELLPS_H
//------------------------------------------------------------------------------