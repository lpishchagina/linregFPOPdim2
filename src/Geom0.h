#ifndef GEOM0_H
#define GEOM0_H

#include <iostream>
#include <vector>
#include <list>
#include <iterator>

#include "Ellps.h"

//Class Geom0
//------------------------------------------------------------------------------
//PELT
//------------------------------------------------------------------------------
class Geom0
{
  private:
    unsigned int label_t;       //time moment 
    std::list<Ellps> ellps;  //list of ellipses(t-1)
    
  public:
    Geom0();
    Geom0(unsigned int t);
    
    unsigned int get_label_t() const;
    std::list<Ellps> get_ellps() const;
    
    void InitialGeometry(unsigned int i, const std::list<Ellps> &ellpses);
    void UpdateGeometry(const Ellps & elt);
    bool EmptyGeometry();
};
#endif //GEOM0_H