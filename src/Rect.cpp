#include "Rect.h"
#include "Ellipse.h"
#include "math.h"
#include<iostream>
#include <vector>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//accessory*********************************************************************
double Rect::get_x0() const {return x0;}

double Rect::get_y0() const {return y0;}

double Rect::get_dx() const {return dx;}

double Rect::get_dy() const {return dy;}

double Rect::min_ab(double a, double b){if (a < b) {return a;} else {return b;}}

double Rect::max_ab(double a, double b){if (a > b) {return a;} else {return b;}}


//correct

//IsEmpty_rect******************************************************************
bool Rect::IsEmpty_rect(){  {return false;} //correct

//EmptyIntersection*************************************************************
bool Rect::EmptyIntersection(const Ellipse &ellipse)
{return false;}

//Intersection_Ellipse*************************************************************
void Rect::Intersection_Ellipse(const Ellipse &ellipse)
{}

//Exclusion_Ellipse****************************************************************
void Rect::Exclusion_Ellipse(const Ellipse &ellipse)
{}
//******************************************************************************



