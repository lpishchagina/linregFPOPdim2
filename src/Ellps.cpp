#include "Ellps.h"
#include "math.h"


#include <Rcpp.h>
using namespace Rcpp;
//constructor*******************************************************************
Ellps::Ellps(lrCost Q)
{
  double inv1 = Q.get_A()+Q.get_C();
  double inv2 = Q.get_A()*Q.get_C() - Q.get_B()*Q.get_B();
  double inv3 = Q.get_A()*Q.get_C()*Q.get_F() + 2*Q.get_B()*Q.get_D()*Q.get_E() - Q.get_C()*Q.get_D()*Q.get_D() - Q.get_A()*Q.get_E()*Q.get_E() - Q.get_B()*Q.get_B()*Q.get_F();
    
  //invariant condition: (theta1^2/C1 + theta2^2/A1 = -F) <=> (A1+C1=inv1) && (A1C1=inv2) && (A1C1F1=inv3) 
  if ((inv2 > 0) && (inv3 != 0) && (inv1*inv3 < 0)){//ellipse condition
    //A1 =lmbd1; C1=inv1-lmbd1; F1=inv3/inv2 ; 
    //a =sqrt(-F1/C1); b = sqrt(-F1/A1); c = 1; angle = atan((A1-A)/B) 
  
    //displacement 
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
    F = sqrt(a*a - b*b);
  } else{ d1 = 0; d1 = 0; a = 0; b = 0; angl = 0; k1 =0; k2 = 0; lmbd1 = 0; lmbd2 = 0; F = 0;}
  Rcpp::Rcout<<" A="<< Q.get_A()<<" B="<< Q.get_B()<<" C="<< Q.get_C()<<" D="<< Q.get_D()<<" E="<< Q.get_E()<<" F="<< Q.get_F()<<" a="<< a<<" b="<< b<<" d1="<< d1<<" d2="<<d2<<" k1="<< k1<<" k2="<< k2<<std::endl;
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
double Ellps::get_F() const {return F;}

//dist_pnts********************************************************************
double Ellps::dst_pnts(double x1, double y1, double x2, double y2){return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));}

//dst_ellps_pnt*****************************************************************
//|F1P|+|F2P| = 2(a+d) => d = (|F1P|+|F2P|)/2 -a
double Ellps::dst_ellps_pnt(double x1, double y1){ return  (abs((dst_pnts(x1,y1,(-F),0) + dst_pnts(x1,y1,(F),0))/2 - a));}


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
bool  Ellps::insd_ellps(Ellps & El){
  //(0,0)
  if (!insd_pnt(x_XY(El.get_d1(), El.get_d2()), y_XY(El.get_d1(), El.get_d2()))){return false;}
  //(-a,0)
  double X = El.X_xy((-El.get_a()), 0);
  double Y = El.Y_xy((-El.get_a()), 0);
  if (!insd_pnt(x_XY(X,Y),y_XY(X,Y))){return false;}
  //(a,0) 
  X = El.X_xy(El.get_a(), 0);
  Y = El.Y_xy(El.get_a(), 0);
  if (!insd_pnt(x_XY(X,Y),y_XY(X,Y))){return false;}
  //(0,-b) 
  X = El.X_xy(0, (-El.get_b()));
  Y = El.Y_xy(0, (-El.get_b()));
  if (!insd_pnt(x_XY(X,Y), y_XY(X,Y))){return false;}
  //(0,-b) 
  X = El.X_xy(0, El.get_b());
  Y = El.Y_xy(0, El.get_b());
  if (!insd_pnt(x_XY(X,Y),y_XY(X,Y))){return false;}
  return true;
}
//****************************************************************************** 
