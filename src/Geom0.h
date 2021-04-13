#ifndef GEOM0_H
#define GEOM0_H

#include <iostream>
#include <vector>
#include <list>
#include <iterator>

#include "Ellipse.h"

//Class Geom0
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
class Geom0
{
  private:
    unsigned int label_t;       //time moment 
    std::list<Ellipse> ellipse_t_1;  //list of disks(t-1)
    
  public:
    Geom0();
    Geom0(unsigned int t);
    
    unsigned int get_label_t() const;
    std::list<Ellipse> get_ellipse_t_1() const;
    
    void InitialGeometry(unsigned int i, const std::list<Ellipse> &ellipses);
    void UpdateGeometry(const Ellipse & ellipse_t);
    bool EmptyGeometry();
};
#endif //GEOM0_H