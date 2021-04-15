#include "Geom3.h"

#include <iostream>
#include <iterator>
#include <list>
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

//Dist**************************************************************************
double Geom3::Dist(double a1, double a2, double b1, double b2){ return sqrt((a1 - b1)*(a1 - b1) +(a2 - b2)*(a2 - b2));}

//InitialGeometry***************************************************************
void Geom3::InitialGeometry(unsigned int i, const std::list<Ellps> &ellpses){
  label_t = i;
  ellps.clear();
  ellps = ellpses;
  fl_empty = false;
}

//UpdateGeometry****************************************************************
void Geom3::UpdateGeometry(const Ellps &elt){
  std::list<Ellps>::iterator iter = ellps.begin();
  while( iter != ellps.end()){
    if ((*iter).insd_ellps(elt)) {
      fl_empty = true;
      return;
    }
    //if (elt.empty_intersection(*iter) ) {iter = disks_t_1.erase(iter);}
    //else {++iter;}
  }
} 

//EmptyGeometry*****************************************************************
bool Geom3::EmptyGeometry() {return fl_empty;}

//******************************************************************************
