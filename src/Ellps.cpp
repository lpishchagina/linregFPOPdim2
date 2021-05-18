#include "Ellps.h"
#include "math.h"

#include <Rcpp.h>
using namespace Rcpp;
//constructor*******************************************************************
Ellps::Ellps(lrCost Q, double mdif)
{
  /*comment-------------------------------------------------
   * forme: Ak^2+2Bak+Ca^2+2Dk+2Ea+fr_mb
   * when free member is fr_mb = F - mdif 
   */         
  double fr_mb = Q.g_F() - mdif;
  
  /*comment------------------------------------------------
   *  invariants of forme
   */
  double inv1 = Q.g_A() + Q.g_C();
  double inv2 = Q.g_A()*Q.g_C() - Q.g_B()*Q.g_B();
  double inv3 = Q.g_A()*Q.g_C()*fr_mb + 2*Q.g_B()*Q.g_D()*Q.g_E() - Q.g_C()*Q.g_D()*Q.g_D() - Q.g_A()*Q.g_E()*Q.g_E() - Q.g_B()*Q.g_B()*fr_mb;
  
  /*comment------------------------------------------------
   *  Matrix transformation:
   *  before: (A B//B C) => after: (A/tr_cnst B/tr_cnst // C/tr_cnst B/tr_cnst)
   *  when tr_cnst = (CD^2 - 2BDE+AE^2)/(AC-B^2) - fr_mb
   *  => M = (A/tr_cnst, B/tr_cnst, C/tr_cnst)
   */
 
  double tr_cnst = (Q.g_C()*Q.g_D()*Q.g_D() - 2*Q.g_B()*Q.g_D()*Q.g_E() +  Q.g_A()*Q.g_E()*Q.g_E())/inv2 - fr_mb;
  M = new double[3];
  M[0] = Q.g_A()/tr_cnst;
  M[1] = Q.g_B()/tr_cnst;
  M[2] = Q.g_C()/tr_cnst;
    
  // - Sigma = Mat2X2(Q.g_A()/tr_cnst, Q.g_B()/tr_cnst, Q.g_B()/tr_cnst, Q.g_C()/tr_cnst);
  
  /*comment------------------------------------------------
   *  Center coordinates of ellipse
   *  x0 = (BE-DC)/(AC-B^2);
   *  y0 = (BD-AE)/(AC-B^2); 
   */
  x0 = (Q.g_B()*Q.g_E() - Q.g_D()*Q.g_C())/inv2;
  y0 = (Q.g_B()*Q.g_D() - Q.g_A()*Q.g_E())/inv2;
  
  /*comment------------------------------------------------
   *  Eigenvalues
   *  lmbd1 >=lmbd2 
   */
  lmbd1 = (inv1 - sqrt(inv1*inv1 - 4*inv2))/2;
  lmbd2 = (inv1 + sqrt(inv1*inv1 - 4*inv2))/2;
  if(lmbd2 > lmbd1){double rt = lmbd1; lmbd1 = lmbd2; lmbd2 = rt;}
  
  /*comment------------------------------------------------
   *  invariant condition:
   *  (theta1^2/C1 + theta2^2/A1 = -F1) <=> (A1+C1=inv1) && (A1C1=inv2) && (A1C1F1=inv3)
   *   => A1 =lmbd1; C1=inv1-lmbd1; F1=inv3/inv2 ; 
   *   a =sqrt(-F1/C1); b = sqrt(-F1/A1); c = 1; angl = atan((A1-A)/B), Slope: k1,k2, Focus:Fs
   */
  k1 = (lmbd1 - Q.g_A())/Q.g_B();
  k2 = (lmbd2 - Q.g_A())/Q.g_B();
  angl = atan((lmbd1 - Q.g_A())/Q.g_B());
  a = sqrt(-(inv3/inv2)/lmbd1);     
  b = sqrt(-(inv3/inv2)/(inv1 - lmbd1));
  Fs = sqrt(a*a - b*b);
}

//constructor copy**************************************************************
Ellps::Ellps(const Ellps & E){
  M = new double[3];
  for(unsigned int i = 0; i < 3; i++){M[i] = E.M[i];}
  x0 = E.x0;
  y0 = E.y0;
  //
  lmbd1 = E.lmbd1;
  lmbd2 = E.lmbd2;
  k1 = E.k1;
  k2 = E.k2;
  angl = E.angl;
  a = E.a;     
  b = E.b;
  Fs = E.Fs;
}

//destructor********************************************************************
Ellps::~Ellps(){
  delete[]M;
  M = NULL;
}

//accessory*********************************************************************
double Ellps::g_x0() const {return x0;}
double Ellps::g_y0() const {return y0;}
double* Ellps::g_M() const {return M;}

double Ellps::g_k1() const {return k1;}
double Ellps::g_k2() const {return k2;}
double Ellps::g_angl() const {return angl;}

double Ellps::g_lmbd1() const {return lmbd1;}
double Ellps::g_lmbd2() const {return lmbd2;}

double Ellps::g_a() const {return a;}
double Ellps::g_b() const {return b;}
double Ellps::g_Fs() const {return Fs;}
// - Mat2X2 Ellps::g_Sigma() const {return Sigma;}

//dist_pnts********************************************************************
double Ellps::dst_pnts(double x1, double y1, double x2, double y2){
  return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

//insd_pnt********************************************************************
bool  Ellps::insd_pnt(double x, double y){ 
  if ((x*x/a*a + y*y/b*b) < 1) {return true;} 
  else {return false;}
}

//g_r*************************************************************************
double Ellps::g_r(double x, double y){      //r = |2rx-x0|*|2ry-y0|/sqrt((2rx-x0)^2sin^2(phi) + (2ry-y0)^2cos^2(phi)), x0,y0 = (0,0)
  double csns = sqrt(1/((y/x)*(y/x)+1));
  double sns = sqrt(1-csns*csns);
  return (4*abs(a*b)/sqrt(2*(b*b*csns*csns + a*a*sns*sns)));
}

//transfer**********************************************************************
//(X,Y)- Descart coord
double  Ellps::X_xy(double x, double y){return x*cos(angl) - y*sin(angl) +x0;}
double  Ellps::Y_xy(double x, double y){return x*sin(angl) + y*cos(angl) +y0;}
//(x,y) - modification
double  Ellps::x_XY(double X, double Y){return (X - x0)*cos(angl) + (Y - y0)*sin(angl);}
double  Ellps::y_XY(double X, double Y){return -(X - x0)*sin(angl) + (Y - y0)*cos(angl);}


//******************************************************************************
