#include "Ellps.h"
#include "math.h"

#include <Rcpp.h>
using namespace Rcpp;
//constructor*******************************************************************
Ellps::Ellps(lrCost Q, double r2)
{
  //r2 = (gCost(t) + beta - gCost(i-1) - beta) => free member fr_mb = F - r2          
  double fr_mb = Q.get_F() - r2;//free member
  
  //invariant condition: (theta1^2/C1 + theta2^2/A1 = -F) <=> (A1+C1=inv1) && (A1C1=inv2) && (A1C1F1=inv3) => A1 =lmbd1; C1=inv1-lmbd1; F1=inv3/inv2 ; 
  //a =sqrt(-F1/C1); b = sqrt(-F1/A1); c = 1; angle = atan((A1-A)/B) 
  double inv1 = Q.get_A() + Q.get_C();
  double inv2 = Q.get_A()*Q.get_C() - Q.get_B()*Q.get_B();
  double inv3 = Q.get_A()*Q.get_C()*fr_mb + 2*Q.get_B()*Q.get_D()*Q.get_E() - Q.get_C()*Q.get_D()*Q.get_D() - Q.get_A()*Q.get_E()*Q.get_E() - Q.get_B()*Q.get_B()*fr_mb;
    
  //Matrix transformation: tr_cnst = (CD^2- 2BDE + AE^2)/(AC-B^2) - fr_mb => Sigma = (1/tr_cnst)*(A B//B C)  
  double tr_cnst = ( Q.get_C()*Q.get_D()*Q.get_D() - 2*Q.get_B()*Q.get_D()*Q.get_E() +  Q.get_A()*Q.get_E()*Q.get_E())/(inv2) - fr_mb;
  Sigma = Mat2X2(Q.get_A()/tr_cnst, Q.get_B()/tr_cnst, Q.get_B()/tr_cnst, Q.get_C()/tr_cnst);
  
  //Center transformation: x0 = (BE-DC)/(AC-B^2); y0 = (BD-AE)/(AC-B^2); 
  x0 = (Q.get_B()*Q.get_E() - Q.get_D()*Q.get_C())/inv2;
  y0 = (Q.get_B()*Q.get_D() - Q.get_A()*Q.get_E())/inv2;
    
  //Eigenvalues: lmbd1 >=lmbd2 
  lmbd1 = (inv1 - sqrt(inv1*inv1 - 4*inv2))/2;
  lmbd2 = (inv1 + sqrt(inv1*inv1 - 4*inv2))/2;
  if(lmbd2 > lmbd1){double rt = lmbd1; lmbd1 = lmbd2; lmbd2 = rt;}
    
  //Slope
  k1 = (lmbd1 - Q.get_A())/Q.get_B();
  k2 = (lmbd2 - Q.get_A())/Q.get_B();
  angl = atan((lmbd1 - Q.get_A())/Q.get_B());
    
  //a,b
  a = sqrt(-(inv3/inv2)/lmbd1);     
  b = sqrt(-(inv3/inv2)/(inv1 - lmbd1));
  //Focus
  Fs = sqrt(a*a - b*b);
}
//accessory*********************************************************************
double Ellps::get_x0() const {return x0;}
double Ellps::get_y0() const {return y0;}

double Ellps::get_k1() const {return k1;}
double Ellps::get_k2() const {return k2;}
double Ellps::get_angl() const {return angl;}

double Ellps::get_lmbd1() const {return lmbd1;}
double Ellps::get_lmbd2() const {return lmbd2;}

double Ellps::get_a() const {return a;}
double Ellps::get_b() const {return b;}
double Ellps::get_Fs() const {return Fs;}
Mat2X2 Ellps::get_Sigma() const {return Sigma;}

//dist_pnts********************************************************************
double Ellps::dst_pnts(double x1, double y1, double x2, double y2){return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));}

//insd_pnt********************************************************************
bool  Ellps::insd_pnt(double x, double y){ if ((x*x/a*a + y*y/b*b) < 1) {return true;} else {return false;}}

//get_r*************************************************************************
double Ellps::get_r(double x, double y){      //r = |2rx-x0|*|2ry-y0|/sqrt((2rx-x0)^2sin^2(phi) + (2ry-y0)^2cos^2(phi)), x0,y0 = (0,0)
  double csns = sqrt(1/((y/x)*(y/x)+1));
  double sns = sqrt(1-csns*csns);
  return (4*abs(a*b)/sqrt(2*(b*b*csns*csns + a*a*sns*sns)));
}

//******************************************************************************
//(X,Y)- Descart coord
double  Ellps::X_xy(double x, double y){return x*cos(angl) - y*sin(angl) +x0;}
double  Ellps::Y_xy(double x, double y){return x*sin(angl) + y*cos(angl) +y0;}
//(x,y) - modification
double  Ellps::x_XY(double X, double Y){return (X - x0)*cos(angl) + (Y - y0)*sin(angl);}
double  Ellps::y_XY(double X, double Y){return -(X - x0)*sin(angl) + (Y - y0)*cos(angl);}


//******************************************************************************
