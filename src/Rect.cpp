#include "Rect.h"
#include "Ellipse.h"
#include "math.h"
#include<iostream>
#include <vector>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//accessory*********************************************************************
double Rect::get_rectx0() const {return rectx0;}

double Rect::get_recty0() const {return recty0;}

double Rect::get_rectx1() const {return rectx1;}

double Rect::get_recty1() const {return recty1;}

double Rect::min_ab(double a, double b){if (a < b) {return a;} else {return b;}}

double Rect::max_ab(double a, double b){if (a > b) {return a;} else {return b;}}

//IsEmpty_rect******************************************************************
bool Rect::IsEmpty_rect(){ if (rectx0 >= rectx1 || recty0 >= recty1) {return true;} else { return false;}}

//EmptyIntersection*************************************************************
bool Rect::EmptyIntersection(const Ellipse &ellipse){
  double c1 = ellipse.get_center1(); 
  double c2 = ellipse.get_center2();
  //point_min-------------------------------------------------------------------
  double pnt_min1 = c1;
  double pnt_min2 = c2;
  if (c1 <= rectx0){ pnt_min1 = rectx0;}
  if (c1 >= rectx1){ pnt_min1 = rectx1;}
  if (c2 <= recty0){ pnt_min2 = recty0;}
  if (c2 >= recty1){ pnt_min2 = recty1;}
  //distance--------------------------------------------------------------------
  if (sqrt((pnt_min1 - c1)*(pnt_min1 - c1) + (pnt_min2 - c2)*(pnt_min2 - c2)) >= ellipse.get_radius()) {return true;}
  else {return false;}
}

//Intersection_Ellipse*************************************************************
void Rect::Intersection_Ellipse(const Ellipse &ellipse){
  double r = ellipse.get_radius();        
  double c1 = ellipse.get_center1(); 
  double c2 = ellipse.get_center2();
  //point_min-------------------------------------------------------------------
  double pnt_min1 = c1;
  double pnt_min2 = c2;
 
  if (c1 <= rectx0){ pnt_min1 = rectx0;}
  if (c1 >= rectx1){ pnt_min1 = rectx1;}
  if (c2 <= recty0){ pnt_min2 = recty0;}
  if (c2 >= recty1){ pnt_min2 = recty1;}
  
  //discriminant----------------------------------------------------------------
  double dxy2_1 =  r * r - (pnt_min2 - c2) * (pnt_min2 - c2);
  double dxy2_2 =  r * r - (pnt_min1 - c1) * (pnt_min1 - c1);
  //----------------------------------------------------------------------------
  if ((dxy2_1 <= 0)||(dxy2_2 <= 0)){rectx1 = rectx0;}
  else{
    if (dxy2_1 > 0){
      rectx0 = max_ab(rectx0, c1 - sqrt(dxy2_1)); 
      rectx1 = min_ab(rectx1, c1 + sqrt(dxy2_1));
    }
    if (dxy2_2 > 0){
      recty0 = max_ab(recty0, c2 - sqrt(dxy2_2)); 
      recty1 = min_ab(recty1, c2 + sqrt(dxy2_2));
    }
  }
}

//Exclusion_Ellipse****************************************************************
void Rect::Exclusion_Ellipse(const Ellipse &ellipse){
  double r = ellipse.get_radius();        
  double c1 = ellipse.get_center1(); 
  double c2 = ellipse.get_center2();
  //-point_max------------------------------------------------------------------
  double pnt_max1;
  double pnt_max2;
  if (abs(c1 - rectx1) >= abs(c1 - rectx0)) {pnt_max1 = rectx1;} else{pnt_max1 = rectx0;}
  if (abs(c2 - recty1) >= abs(c2 - recty0)) {pnt_max2 = recty1;} else{pnt_max2 = recty0;}
  //discriminant----------------------------------------------------------------
  double dxy2_1 =  r * r - (pnt_max2 - c2) * (pnt_max2 - c2);
  double dxy2_2 =  r * r - (pnt_max1 - c1) * (pnt_max1 - c1);
  //----------------------------------------------------------------------------
  if (dxy2_1 > 0){
    if ((pnt_max1 == rectx0) && (rectx1 <= c1 + sqrt(dxy2_1))) {rectx1 = min_ab(rectx1, c1 - sqrt(dxy2_1));}
    if ((pnt_max1 == rectx1) && (rectx0 >= c1 - sqrt(dxy2_1))) {rectx0 = max_ab(rectx0, c1 + sqrt(dxy2_1));}
  }
  if (dxy2_2 > 0){
    if ((pnt_max2 == recty0) && (recty1 <= c2 + sqrt(dxy2_2))) {recty1 = min_ab(recty1, c2 - sqrt(dxy2_2));}
    if ((pnt_max2 == recty1) && (rectx0 >= c2 - sqrt(dxy2_2))) {recty0 = max_ab(recty0, c2 + sqrt(dxy2_2));}
  }
}
//******************************************************************************



