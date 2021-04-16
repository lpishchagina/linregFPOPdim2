#include "Geom0.h"

#include <iostream>
#include <iterator>
#include <list>
#include <math.h>

using namespace std;

//constructor*******************************************************************
Geom0::Geom0(){
  label_t = 0;
}
Geom0::Geom0(unsigned  int t){
  label_t = t; 
}
//accessory*********************************************************************
unsigned int Geom0::get_label_t() const {return label_t;}
std::list<Ellps> Geom0::get_ellps() const {return ellps;}

//InitialGeometry***************************************************************
void Geom0::InitialGeometry(unsigned int i, const std::list<Ellps> &ellpses){label_t = i;}

//UpdateGeometry****************************************************************
void Geom0::UpdateGeometry(Ellps &elt){} 

//EmptyGeometry*****************************************************************
bool Geom0::EmptyGeometry() {return false;}

//******************************************************************************
