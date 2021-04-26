#include "Geom3.h"

#include <iostream>

#include <list>
#include <iterator>
#include <math.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//constructor*******************************************************************
Geom3::Geom3(){
  label_t = 0;
  fl_empty = false;
}

Geom3::Geom3(unsigned  int t){
  label_t = t; 
  fl_empty = false;
}

//accessory*********************************************************************
unsigned int Geom3::get_label_t() const {return label_t;}
std::list<Ellps> Geom3::get_ellps()const {return ellps;}

//(dsX,dsY)- Descart coord, (x,y) - modification********************************
//from decart(x,y) to modif (X,Y)
double  Geom3::dsX_to_x(double dsx, double dsy, double dx, double dy, double angle){return (dsx - dx)*cos(angle) + (dsy - dy)*sin(angle);}
double  Geom3::dsY_to_y(double dsx, double dsy, double dx, double dy, double angle){return -(dsx - dx)*sin(angle) + (dsy - dy)*cos(angle);}
//from modif (X,Y) to decart(x,y) 
double  Geom3::x_to_dsX(double x, double y, double dx, double dy, double angle){return x*cos(angle) - y*sin(angle) + dx;}
double  Geom3::y_to_dsY(double x, double y, double dx, double dy, double angle){return x*sin(angle) + y*cos(angle) + dy;}

//Dist pnt_pnt******************************************************************
double Geom3::Dist(double a1, double a2, double b1, double b2){ return sqrt((a1 - b1)*(a1 - b1) +(a2 - b2)*(a2 - b2));}

//pnt_insd_E********************************************************************
bool Geom3::pnt_insd_E(double x, double y, const Ellps E){ 
  if ((x*x/E.get_a()*E.get_a() + y*y/E.get_b()*E.get_b()) < 1) { return true;} 
  else {return false;}
}

//InitialGeometry***************************************************************
void Geom3::InitialGeometry(unsigned int i, const std::list<Ellps> &ellpses){
  label_t = i;
  ellps.clear();
  ellps = ellpses;
  fl_empty = false;
}

//UpdateGeometry****************************************************************
void Geom3::UpdateGeometry (const Ellps &Et){
  double d_Ei_cEt, d_Et_cEi;// distance from first ellipse to center  of second ellipse
  std::list<Ellps>::iterator Ei = ellps.begin();
  /*
  while( Ei != ellps.end()){
    //if Et center inside (*Ei)
    if (pnt_insd_E(dsX_to_x(Et.get_d1(), Et.get_d2(), (*Ei).get_d1(), (*Ei).get_d2(),(*Ei).get_angl()), dsY_to_y(Et.get_d1(),Et.get_d2(),(*Ei).get_d1(), (*Ei).get_d2(), (*Ei).get_angl()), (*Ei))){
      if (E1_insd_E2(Et, (*Ei))){
        fl_empty = true;
        return;
      }
      else{++Ei;}     // => exist only intersection next *Ei 
    }
    else{
      //if (*Ei) center inside Et
      if (pnt_insd_E(dsX_to_x((*Ei).get_d1(), (*Ei).get_d2(), Et.get_d1(), Et.get_d2(), Et.get_angl()), dsY_to_y((*Ei).get_d1(), (*Ei).get_d2(), Et.get_d1(), Et.get_d2(), Et.get_angl()), Et))
        { ++Ei;}// => exist intersection next *Ei   
      //(*Ei)center outside Et and Et center outside (*Ei)
      else{  
     
        if (condition) {Ei = ellps.erase(Ei);} //=>empty_intersection
        else{++Ei;} // => exist intersection next *Ei 
     }
    }
  }*/
} 

//EmptyGeometry*****************************************************************
bool Geom3::EmptyGeometry() {return fl_empty;}

//******************************************************************************

int Geom3::testIntersection(const Ellps &E1,const Ellps &E2)
{
  int res = E1.testInter(E2);
  /*
  //Ellipse E1
  Mat2X2 Sig1 = E1.get_Sigma();
  double* c1 = new double[2];
  c1[0] = E1.get_d1();
  c1[1] = E1.get_d2();
  //Ellipse E2
  Mat2X2 Sig2 = E2.get_Sigma();
  double* c2 = new double[2];
  c2[0] = E2.get_d1();
  c2[1] = E2.get_d2();
  
  int res; //  if isn't intersection res = 2, else res = 0, 1 
  
   */
}
