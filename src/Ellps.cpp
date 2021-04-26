#include "Ellps.h"
#include "math.h"

#include <Rcpp.h>
using namespace Rcpp;
//constructor*******************************************************************
Ellps::Ellps(lrCost Q, double r2)
{
  Sigma = Mat2X2(Q.get_A(), Q.get_B(), Q.get_B(), Q.get_C());
  
  double inv1 = Q.get_A()+Q.get_C();
  double inv2 = Q.get_A()*Q.get_C() - Q.get_B()*Q.get_B();
  double inv3 = Q.get_A()*Q.get_C()*(Q.get_F() + r2) + 2*Q.get_B()*Q.get_D()*Q.get_E() - Q.get_C()*Q.get_D()*Q.get_D() - Q.get_A()*Q.get_E()*Q.get_E() - Q.get_B()*Q.get_B()*(Q.get_F() + r2);
    
  //invariant condition: (theta1^2/C1 + theta2^2/A1 = -F) <=> (A1+C1=inv1) && (A1C1=inv2) && (A1C1F1=inv3) 
  // ((inv2 > 0) && (inv3 != 0) && (inv1*inv3 < 0)) ellipse condition
    //A1 =lmbd1; C1=inv1-lmbd1; F1=inv3/inv2 ; 
    //a =sqrt(-F1/C1); b = sqrt(-F1/A1); c = 1; angle = atan((A1-A)/B) 
  
    //center
    d1 = (Q.get_B()*Q.get_E() - Q.get_D()*Q.get_C())/inv2;
    d2 = (Q.get_D()*Q.get_B() - Q.get_A()*Q.get_E())/inv2;
    
    //lmbd1 >=lmbd2 
    lmbd1 = (inv1 - sqrt(inv1*inv1 - 4*inv2))/2;
    lmbd2 = (inv1 + sqrt(inv1*inv1 - 4*inv2))/2;
    if(lmbd2 > lmbd1){double rt = lmbd1; lmbd1 = lmbd2; lmbd2 = rt;}
    
    //slope
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
double Ellps::get_d1() const {return d1;}
double Ellps::get_d2() const {return d2;}

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
double  Ellps::X_xy(double x, double y){return x*cos(angl) - y*sin(angl) +d1;}
double  Ellps::Y_xy(double x, double y){return x*sin(angl) + y*cos(angl) +d2;}
//(x,y) - modification
double  Ellps::x_XY(double X, double Y){return (X - d1)*cos(angl) + (Y - d2)*sin(angl);}
double  Ellps::y_XY(double X, double Y){return -(X - d1)*sin(angl) + (Y - d2)*cos(angl);}


//******************************************************************************

int Ellps::testInter(const Ellps &E2)
{
  //Ellipse E1
  Mat2X2 Sig1 = Sigma;
  double* c1 = new double[2];
  c1[0] = d1;
  c1[1] = d2;
  //Ellipse E2
  Mat2X2 Sig2 = E2.get_Sigma();
  double* c2 = new double[2];
  c2[0] = E2.get_d1();
  c2[1] = E2.get_d2();
  
  int res; //  if isn't intersection res = 2, else res = 0, 1 
  
  double* m1 = Sig1.MatMultVec(c1);
  double* m2 = Sig2.MatMultVec(c2);
  m1[0] = m1[0] - m2[0];
  m1[1] = m1[1] - m2[1];
  
  Mat2X2 N = (Sig1.DifMat(Sig2)).Obrat();
  
  Mat2X2 M = (Sig2.MultMat(N)).MultNb(-1);
  
  Mat2X2 Told = Mat2X2();
  Mat2X2 I = Mat2X2(1,0,0,1);
  Mat2X2 T = Mat2X2();
  
  double* P = new double [3];
  for (unsigned int k = 0; k<2; k++){P[k] = 0;}
  
  double* ai = new double [3];
  ai[2] = 1
  unsigned int i = 2;
  
  while (i > 0){
    T = (M.MultMat(Told)).SumMat(I.MultNb(ai[i]));
    
    P[i-1] = P[i-1] + T.VecMultMat(N.VecMultMat(m2))[0]*m2[0] + T.VecMultMat(N.VecMultMat(m2))[1]*m2[1];
    P[i] = P[i] + 2*(T.VecMultMat(N.VecMultMat(m1))[0]*m2[0] + T.VecMultMat(N.VecMultMat(m1))[1]*m2[1]);
    P[i+1] = P[i+1] + + T.VecMultMat(N.VecMultMat(m1))[0]*m1[0] + T.VecMultMat(N.VecMultMat(m1))[1]*m1[1];
    
    Told = T;
    i = i-1;
    ai[i] = (M.MultMat(T)).Trace()/(i - 2);
  }
  
  //obtenir p(lambda) => Shturm
}

