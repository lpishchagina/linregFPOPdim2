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

//EmptyGeometry*****************************************************************
bool Geom3::EmptyGeometry() {return fl_empty;}

//seqShturm*********************************************************************
/* INPUT : l -  specified value
 * C0,C1,C2,C3 - coefficient of equation f(x) = C3*x^3+C2*x^2+C1*x+C0  
 * 
 * OUTPUT : The function returns the number of changes of signs in the Shturm sequence at the specified value "l"
 */
unsigned int Geom3::seqShturm(double l, double C0, double C1, double C2, double C3){
  /*comment------------------------------------------------
   * Shturm sequence
   * f(l) = C3*l^3+C2*l^2+C1*l+C0
   * f1(l) = f'(l) =3*C3*l^2+2*C2*l+C1
   * f2(l) = C5*l + C4
   * f3(l) = (C4/C5)*(2*C2 - 3*C3*C4/C1) - C1
   * when 
   * C4 = ((1/9)*(C1*C2)/C3 - C0)
   * C5 = ((2/9)*(C2*C2/C3) - (2/3)*C1)
   */
  double C4 = (1/9)*(C1*C2)/C3 - C0;
  double C5 = (2/9)*(C2*C2)/C3 - (2/3)*C1;
  
  /*comment------------------------------------------------
   * Values of Shturm sequence
   * Sh_f = (f(l),f1(l),f2(l), f3(l))
   */
  double* Sh_f = new double [4];
  
  Sh_f[0] = C3*l*l*l + C2*l*l + C1*l + C0;
  Sh_f[1] = 3*C3*l*l + 2*C2*l + C1;
  Sh_f[2] = C5*l + C4;
  Sh_f[3] = C4/C5* (2*C2 - 3*C3*C4/C1) - C1;
  
  /*comment------------------------------------------------
   * count - number of changes of signs in the Shturm sequence 'Sh_f'
   */
  unsigned int count = 0;
  for (unsigned int i = 0; i < 3 ; i++){
    if (Sh_f[i] < 0 && Sh_f[i+1] > 0 || Sh_f[i] > 0 && Sh_f[i+1] < 0) {count++;}
  }
  /*memory cleaning-------------------------------------*/
  delete []Sh_f;
  Sh_f = NULL;
  
  return count;
}

//InterKalman*********************************************************************
/* INPUT : parameters of 2 ellipses:
 * 1 : (c1, c2) -center (x0,y0), A -  M-values of matrix (a11,a12, a12,a22)
 * 2 : (d1, d2) -center (x0,y0), B -  M-values of matrix (b11,b12, b12,b22)
 * OUTPUT : The function returns the number of roots of the Kalman filter
 * NOTE : if roots = 2 => exist intersection
 */
unsigned int Geom3::InterKalman(double c1, double c2, const double* &A, double d1, double d2, const double* &B){ 
  /*comment------------------------------------------------
   * E - linear combination of matrix  2 ellipses
   * E(l) = l*A +(1-l)*B
   * detE(l) =  k2*l^2 + k1*l + k0
   * when
   * k2 = (b11*b22-b12^2)+(b11*b22-b12^2) - (a22b11+a11b22-2a12b12) = detA + detB - detcnst
   * k1 = (a22b11+a11b22-2a12b12) - 2(b11*b22-b12^2) = detcnst - 2*detB
   * k0 = b11*b22-b12^2 = detB
   */
  double detA = A[0]*A[2] - A[1]*A[1];
  double detB = B[0]*B[2] - B[1]*B[1];
  double detcnst = A[2]*B[0] + A[0]*B[2] - 2* A[1]*B[1];

  /*comment------------------------------------------------
   * Kalman filter k(l) = 1- l(1-l)(d-c)^T*B*E^(-1)*A(d-c) =>
   * detE(l)*K(l) = detE(l) + detE(l)*(l^2-l)(d-c)^T*B*E(l)^(-1)*A(d-c)
   *======>*/

  /*comment------------------------------------------------
   * (d-c)^T*B = (l1,l2) => l1 = b11(d1-c1) + b12(d2-c2); l2 = b12(d1-c1) + b22(d2-c2)
   */
  double l1 = B[0]*(d1-c1) + B[2]*(d2-c2);
  double l2 = B[2]*(d1-c1) + B[3]*(d2-c2);
  
  /*comment------------------------------------------------
   * A*(d-c) = (n1,n2)^T =>  n1 = a11(d1-c1) + a12(d2-c2); n2 = a12(d1-c1) + a22(d2-c2)
   */
  double n1 = A[0]*(d1-c1) + A[2]*(d2-c2);
  double n2 = A[1]*(d1-c1) + A[2]*(d2-c2);
  
  /*comment------------------------------------------------
   * detE(l)*(l1,l2)*E(l)^(-1)*(n1, n2)^T = z1*l + z2
   * z1 = n1*l1(a22-b22) + n2*l2(a11-b11)   +  (a12-b12)(n1*l2 - n2*l1) 
   * z2 = n1(l1*b22-l2*b12) + n2(l2*b11-l1*b12)
   */
  double z1 = n1*l1*(A[2]-B[2]) + n2*l2*(A[0]-B[0]) + (A[1]-B[1])*(n1*l2 - n2*l1);
  double z2 = n1*(l1*B[2]-l2*B[1]) + n2*(l2*B[0]-l1*B[1]);
  
  /*comment------------------------------------------------
   * detE(l)*K(l) = detE(l) + (l^2-l)(z1*l +z2) = (k2*l^2 + k1*l + k0) + (z1*l^3 + (z2 -z1)*l^2 - z2*l) =>
   * detE(l)*K(l) = C3*l^3 + C2*l^2 + C1*l + C0
   * when
   * C3 = z1
   * C2 = k2 + z2 - z1
   * C1 = k1 - z2
   * C0 = k0
   */
  double C0 = detB; 
  double C1 = detcnst - 2*detB - z2; 
  double C2 = detA + detB - detcnst + z2 - z1; 
  double C3 = z1; 
  unsigned int count = seqShturm(0, C0, C1, C2, C3) - seqShturm(1, C0, C1, C2, C3);
  return count;
}

//ExclKalman*********************************************************************
/* INPUT : parameters:
 * 1 : (c1, c2) -center (x0,y0), A -  M-values of matrix (a11,a12, a12,a22)
 * 2 : (d1, d2) -center (x0,y0), B -  M-values of matrix (b11,b12, b12,b22)
 * OUTPUT : The function returns the number of roots of the Kalman filter for lmb1 and lmb2
 * NOTE : if roots = 4 => exist exclusion
 */
unsigned int Geom3::ExclKalman(double c1, double c2, const double* &A, double d1, double d2, const double* &B){ 
  /*comment------------------------------------------------
   * E - linear combination: E(l) = l*A +(1-l)*(-B) = l*A +(l-1)*B => detE(l) =  k2*l^2 + k1*l + k0
   * when k2 = detA + detB - detcnst;  k1 = detcnst - 2*detB; k0 = detB
   */
  double detA = A[0]*A[2] - A[1]*A[1];
  double detB = B[0]*B[2] - B[1]*B[1];
  double detcnst = 2*A[1]*B[1] - A[2]*B[0] - A[0]*B[2];
  
  /*comment------------------------------------------------
   * Kalman filter k(l) = 1- l(1-l)(d-c)^T*(-B)*E^(-1)*A(d-c) =>
   * detE(l)*K(l) = detE(l) + detE(l)*(l^2-l)(d-c)^T*(-B)*E(l)^(-1)*A(d-c)
   *======>*/
  
  /*comment------------------------------------------------
   * (d-c)^T*(-B) = (l1,l2) => 
   * l1 = b11(d1-c1) + b12(d2-c2); l2 = b12(d1-c1) + b22(d2-c2)
   */
  double l1 = B[0]*(c1-d1) + B[2]*(c2-d2);
  double l2 = B[2]*(c1-d1) + B[3]*(c2-d2);
  
  /*comment------------------------------------------------
   * A*(d-c) = (n1,n2)^T =>  n1 = a11(d1-c1) + a12(d2-c2); n2 = a12(d1-c1) + a22(d2-c2)
   */
  double n1 = A[0]*(d1-c1) + A[2]*(d2-c2);
  double n2 = A[1]*(d1-c1) + A[2]*(d2-c2);
  
  /*comment------------------------------------------------
   * detE(l)*(l1,l2)*E(l)^(-1)*(n1, n2)^T = z1*l + z2
   * z1 = n1*l1(a22+b22) + n2*l2(a11+b11) - (a12+b12)(n1*l2+n2*l1) 
   * z2 = n1(l2*b12-l1*b22) + n2(l1*b12-l2*b11)
   */
  double z1 = n1*l1*(A[2]+B[2]) + n2*l2*(A[0]+B[0]) - (A[1]-B[1])*(n1*l2 + n2*l1);
  double z2 = n1*(l2*B[1]-l1*B[2]) + n2*(l1*B[1]-l2*B[0]);
  
  /*comment------------------------------------------------
   * detE(l)*K(l) = detE(l) + (l^2-l)(z1*l +z2) = (k2*l^2 + k1*l + k0) + (z1*l^3 + (z2 -z1)*l^2 - z2*l) =>
   * detE(l)*K(l) = C3*l^3 + C2*l^2 + C1*l + C0
   * when C3 = z1; C2 = k2 + z2 - z1; C1 = k1 - z2; C0 = k0
   */
  double C0 = detB; 
  double C1 = detcnst - 2*detB - z2; 
  double C2 = detA + detB - detcnst + z2 - z1; 
  double C3 = z1; 
  
  /*comment------------------------------------------------
   * Eigenvalues
   */
  double detApB = (B[0]+A[0])*(A[2]+B[2]) - (A[1]+B[1])*(A[1]+B[1]);
  
  double b = 2*B[1]*(A[1]+B[1]) - B[2]*(A[0]+B[0]) - B[0]*(A[2]+B[2]); 
  
  double D = b*b - 4*detB*detApB;
  
  double lb1 = ((-1)*b - sqrt(D))/(2*detApB);
  double lb2 = ((-1)*b + sqrt(D))/(2*detApB);
  
  if (lb1 > lb2) {D = lb1; lb1 = lb2; lb2 = lb1;} //(lb1<lb2)&&((lb1*lb2)>0)
  
  unsigned int count = (seqShturm(0,C0,C1,C2,C3)-seqShturm(lb1,C0,C1,C2,C3))+(seqShturm(lb2,C0,C1,C2,C3)-seqShturm(1,C0,C1,C2,C3));
  return count;
}
/*
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
  
} 

 */
//##############################################################################