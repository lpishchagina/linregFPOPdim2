#include "Geom0.h"

#include <iostream>
#include <iterator>
#include <list>
#include <math.h>

using namespace std;

//constructor*******************************************************************
Geom0::Geom0()
{
  label_t = 0;
}
Geom0::Geom0(unsigned  int t)
{
  label_t = t; 
}



//accessory*********************************************************************
unsigned int Geom0::get_label_t() const {return label_t;}
std::list<Disk> Geom0::get_disks_t_1() const {return disks_t_1;}

//InitialGeometry***************************************************************
void Geom0::InitialGeometry(unsigned int i, const std::list<Disk> &disks){label_t = i;}

//UpdateGeometry****************************************************************
void Geom0::UpdateGeometry(const Disk &disk_t){} 

//EmptyGeometry*****************************************************************
bool Geom0::EmptyGeometry() {return false;}

//******************************************************************************
