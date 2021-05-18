#ifndef ELLPS_H
#define ELLPS_H

#include "lrCost.h"
#include "Mat2x2.h"

#include <vector>
/*
 Class Ellps
 -------------------------------------------------------------------------------
 Description: 
 Ellps  that corresponds of Fun(k,a) = A*k^2+2*B*k*a+C*a^2+2*D*k+2*E*a+fr_mb = 0
(x0,y0) - center
 M = (m11, m12,m22) of matrix (m11,m12//m12,m22)
 -------------------------------------------------------------------------------
 */
class Ellps
{
private:
  //center
  double x0, y0;
  //elements of matrix
  double* M; //3 values A, B and C of Matrix (A,B//B,C)
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
  //Mat2X2 Sigma;
 
public:
  Ellps(){};
  Ellps(const Ellps & E);
  ~Ellps();
  Ellps(lrCost Q, double mdif);
  //accessory
  double* get_M() const;
  double get_x0() const;
  double get_y0() const;
 
  double get_a() const;
  double get_b() const;
  double get_Fs() const;
  double get_lmbd1() const;
  double get_lmbd2() const;
  double get_k1() const;
  double get_k2() const;
  double get_angl() const;
  
  //Mat2X2 get_Sigma() const;
  
  bool insd_pnt(double x, double y);
  double dst_pnts(double x1, double y1, double x2, double y2);
  //transfer of coordinates
  double x_XY(double X, double Y);
  double y_XY(double X, double Y);
    
  double X_xy(double x, double y);
  double Y_xy(double x, double y);
  
  double get_r(double x, double y);
}; 

#endif //ELLPS_H
//------------------------------------------------------------------------------