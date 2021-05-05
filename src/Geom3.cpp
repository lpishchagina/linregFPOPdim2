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
   if (filtKalman(Et.get_x0(), Et.get_y0(), Et.get_Sigma(), Ei.get_x0(), Ei.get_y0(), Ei.get_Sigma()) ==2) {Ei = ellps.erase(Ei);} 
   else{
     if (condition excl) {
       fl_empty = true;
       return;
     } 
    else{++Ei;} // => exist intersection next *Ei 
   }
  }
  */
} 

//EmptyGeometry*****************************************************************
bool Geom3::EmptyGeometry() {return fl_empty;}

//******************************************************************************
unsigned int Geom3::filtKalman(double c1, double c2, const Mat2X2 &A, double d1, double d2, const Mat2X2 &B)
{ 
  // E(lambda) = lambda*A +(1-lambda)*B
  //polynome: det(E(lambda)) = k2*lambda^2 + k1*lambda +k0
  double k_cnst = A.get_el(1,1)*B.get_el(0,0) + A.get_el(0,0)*B.get_el(1,1) - 2*A.get_el(0,1)*B.get_el(0,1);
  
  double k0 = B.get_el(0,0)*B.get_el(1,1) - B.get_el(0,1)*B.get_el(0,1);
  double k1 = k_cnst = k_cnst - 2* k0;
  double k2 = A.get_el(0,0)*A.get_el(1,1) - A.get_el(0,1)*A.get_el(0,1) + k0 - k_cnst;
  
  //det(E(lambda))K(lambda) = det(E(lambda)) + (lambda^2-lambda)((d-c)^T )B(E(lambda)^(-1))A(d-c) =>
  //det(E(lambda)) *K(lambda) = z1*lambda^3 + (k2+z2-z1)*lambda^2 + (k1-z2)*lambda + k0
    
  //(d-c)^T*B = (l1,l2) => l1 = b11d1-b11c1+b12d2-b12c2; l2 = b12d1-b12c1+b22d2-b22c2
  double l1 = B.get_el(0,0)*(d1 - c1) + B.get_el(0,1)*(d2 - c2);
  double l2 = B.get_el(0,1)*(d1 - c1) + B.get_el(1,1)*(d2 - c2);
  
  //A*(d-c) = (n1,n2)^T => n1 = a11d1-a11c1+a12d2-a12c2; n2 = a12d1-a12c1+a22d2-a22c2
  double n1 = A.get_el(0,0)*(d1 - c1) + A.get_el(0,1)*(d2 - c2);
  double n2 = A.get_el(0,1)*(d1 - c1) + A.get_el(1,1)*(d2 - c2);
  
  //z1 = n1l1a22 - n1l1b22 + n1l2a12 +n2l1b12 -n2l1a12 +n2l2a11 - n2l2b11
  //z2 = n1l1b22 - n1l2b12 + n2l2b11 - n2l1b12
  
  double z1 = n1*l1*A.get_el(1,1) - n1*l1*B.get_el(1,1) + n1*l2*A.get_el(0,1) +n2*l1*B.get_el(0,1) -n2*l1*A.get_el(0,1) +n2*l2*A.get_el(0,0) - n2*l2*B.get_el(0,0);
  double z2 = n1*l1*B.get_el(1,1) - n1*l2*B.get_el(0,1) + n2*l2*B.get_el(0,0) - n2*l1*B.get_el(0,1);
  
  /*Shturm sequence: f_1 = f' , fi = -(f_(i-2) mod f_(i-1)), i = 2,3 
  det(E(lambda)) *K(lambda) = C3*lambda^3 + C2*lambda^2 + C1*lambda + C0
   */
  double C0 = k0; 
  double C1 = k1 - z1; 
  double C2 = k2 + z2 - z1; 
  double C3 = z1; 
  double CC0 = C1*C2/(9*C3) - C0; 
  double CC1 =  2*C2*C2/(9*C3)  - 2*C1/3;
  
  unsigned int count = seqShturm(0, C0, C1, C2, C3,CC0, CC1) - seqShturm(1, C0, C1, C2, C3,CC0, CC1);
  return count;
}
//------------------------------------------------------------------------------
unsigned int Geom3::seqShturm(double lmb, double C0, double C1, double C2, double C3,double CC0,double CC1)
{
  double* Sh_f = new double [4];
  Sh_f[0] = C3*lmb*lmb*lmb + C2*lmb*lmb + C1* lmb +C0;
  Sh_f[1] = 3*C3*lmb*lmb +2*C2*lmb + C1;
  Sh_f[2] = CC1*lmb + CC0;
  Sh_f[3] = CC0/CC1* (2*C2 - 3*C3*CC0/C1) - C1;
  unsigned int count = 0;
  //count sign
  for (unsigned int i = 0; i < 3 ; i++){
    if (Sh_f[i] < 0 && Sh_f[i+1] > 0 || Sh_f[i] > 0 && Sh_f[i+1] < 0) {count++;}
  }
  //memory
  delete []Sh_f;
  Sh_f = NULL;
  return count;
}
//#Stock########################################################################

/*
 functions: lambda = 0, lambda = 1,
 f_0 = C0;
 f_1 = C3+C2+C1+C0;
 f1_0 = C1;
 f1_1 = 3*C3 +2*C2 +C1;
 f2_0 = CC0;
 f2_1 = CC1 +CC0;
 f3_01 = CC0*(2*C2-(3*C3*CC0)/C1)/CC1 - C1; 
 */
/*
 double* Sh_f0 = new double [4];
 double* Sh_f1 = new double [4];
 Sh_f0[0] = k0;
 Sh_f1[0] = k2 + z2 + k1 - z1 + k0;
 Sh_f0[1] = k1 - z1;
 Sh_f1[1] = 2*(k2 + z2) + k1;
 Sh_f0[2] = (k1 - z1)*(k2 + z2 - z1)/(9*z1) - k0;
 Sh_f1[2] = 2*(k2 + z2 - z1)*(k2 + z2 - z1)/(9*z1)  - 2*(k1 - z1)/3 + (k1 - z1)*(k2 + z2 - z1)/(9*z1) - k0;
 Sh_f0[3] = ((k1 - z1)*(k2 + z2 - z1)/(9*z1) - k0) * (2*(k2 + z2 - z1)-(3*z1*((k1 - z1)*(k2 + z2 - z1)/(9*z1) - k0))/(k1 - z1))/(2*(k2 + z2 - z1)*(k2 + z2 - z1)/(9*z1)  - 2*(k1 - z1)/3) - (k1 - z1); 
 Sh_f1[3] =  Sh_f0[3];
 unsigned int count0 = 0;
 unsigned int count1 = 0;
 
 //count sign
 for (unsigned int i = 0; i < 3 ; i++){
 if (Sh_f0[i] < 0 && Sh_f0[i+1] > 0 || Sh_f0[i] > 0 && Sh_f0[i+1] < 0) {count0++;}
 if (Sh_f1[i] < 0 && Sh_f1[i+1] > 0 || Sh_f1[i] > 0 && Sh_f1[i+1] < 0) {count1++;}
 }
 //memory  
 delete []Sh_f0;
 delete []Sh_f1;
 Sh_f0 = NULL;
 Sh_f1 = NULL;
 return (count0-count1); 
 */
//##############################################################################