#include "Rctngl.h"

#include <Rcpp.h>
using namespace Rcpp;

//accessory*********************************************************************
double Rctngl::get_x1() const {return x1;}
double Rctngl::get_y1() const {return y1;}
double Rctngl::get_x2() const {return x2;}
double Rctngl::get_y2() const {return y2;}
double Rctngl::get_x3() const {return x3;}
double Rctngl::get_y3() const {return y3;}
double Rctngl::get_x4() const {return x4;}
double Rctngl::get_y4() const {return y4;}

double Rctngl::min_ab(double a, double b){if (a < b) {return a;} else {return b;}}
double Rctngl::max_ab(double a, double b){if (a > b) {return a;} else {return b;}}

//(dsX,dsY)- Descart coord, (x,y) - modification********************************
//from decart(x,y) to modif (X,Y)
double Rctngl::dsX_to_x(double dsx, double dsy, double dx, double dy, double angle){return (dsx - dx)*cos(angle) + (dsy - dy)*sin(angle);}
double Rctngl::dsY_to_y(double dsx, double dsy, double dx, double dy, double angle){return -(dsx - dx)*sin(angle) + (dsy - dy)*cos(angle);}
//from modif (X,Y) to decart(x,y) 
double Rctngl::x_to_dsX(double x, double y, double dx, double dy, double angle){return x*cos(angle) - y*sin(angle) + dx;}
double Rctngl::y_to_dsY(double x, double y, double dx, double dy, double angle){return x*sin(angle) + y*cos(angle) + dy;}

//IsEmpty_rect******************************************************************
bool Rctngl::IsEmpty_Rctngl(){ if ((x1 >= x4)&&(x2 >= x3) || ((y1 >= y2)&&(y4 >= y3))) {return true;} else { return false;}}

//pnt_insd_trngl******************************************************************
bool Rctngl::pnt_insd_trngl(double x, double y, double a1, double a2, double b1, double b2, double c1, double c2){
  double a = (a1 - x) * (b2 - y) - (b1 - a1) * (a2 - y);
  double b = (b1 - x) * (c2 - y) - (c1 - b1) * (b2 - y);
  double c = (c1 - x) * (a2 - y) - (a1 - c1) * (c2 - y);
  if ((a >= 0 && b >= 0 && c >= 0) || (a <= 0 && b <= 0 && c <= 0)){return true;} else {return false;}
}

//insd_pnt********************************************************************
bool  Rctngl::insd_pnt(double x, double y){ if ((x*x/a*a + y*y/b*b) < 1) {return true;} else {return false;}}

//dist_pnts********************************************************************
double Rctngl::dst_pnts(double a1, double b1, double a2, double b2){return sqrt((a1-a2)*(a1-a2) + (b1-b2)*(b1-b2));}

//get_r*************************************************************************
double Rctngl::get_r(double x, double y, const Ellps &E){      //r = |2rx-x0|*|2ry-y0|/sqrt((2rx-x0)^2sin^2(phi) + (2ry-y0)^2cos^2(phi)), x0,y0 = (0,0)
  double csns = sqrt(1/((dsY_to_y(x,y,E.get_d1(),E.get_d2(),E.get_angl())/dsX_to_x(x,y,E.get_d1(),E.get_d2(),E.get_angl()))*(dsY_to_y(x,y,E.get_d1(),E.get_d2(),E.get_angl())/dsX_to_x(x,y,E.get_d1(),E.get_d2(),E.get_angl()))+1));
  double sns = sqrt(1-csns*csns);
  return (4*abs(E.get_a()*E.get_b())/sqrt(2*(E.get_b()*E.get_b()*csns*csns + E.get_a()*E.get_a()*sns*sns)));
}

//EmptyIntersection*************************************************************
bool Rctngl::EmptyIntersection(const Ellps &E){
  //Ellipse center inside rectangle => intersection
  //we divide the rectangle (1,2,3,4)into 2 triangles (1,2,3) and (1,3,4) 
  if (pnt_insd_trngl( E.get_d1(), E.get_d2(), x1, y1, x2, y2, x3, y3)){return false;}
  if (pnt_insd_trngl( E.get_d1(), E.get_d2(), x1, y1, x3, y3, x4, y4)){return false;}
  
  //Ellipse center outside rectangle
  //ERROR!!!!переделать.. надо искать кратчайшие расстояния до граней
  //Ellipse center outside rectangle => no intersection if all distance (vertex, center) > Ellipse radius
  if ( get_r(x1, y1, E)>= dst_pnts(x1, y1, E.get_d1(),E.get_d2())){return false;}
  if ( get_r(x2, y2, E)>= dst_pnts(x2, y2, E.get_d1(),E.get_d2())){return false;}
  if ( get_r(x3, y3, E)>= dst_pnts(x3, y3, E.get_d1(),E.get_d2())){return false;}
  if ( get_r(x4, y4, E)>= dst_pnts(x4, y4, E.get_d1(),E.get_d2())){return false;}
  return true;
}