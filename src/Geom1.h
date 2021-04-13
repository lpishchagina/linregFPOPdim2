#ifndef GEOM1_H
#define GEOM1_H

#include <iostream>
#include <vector>
#include <list>

#include "Ellipse.h"
#include "Rect.h"
#include "lrCost.h"

//Class Geom1
//-----------------------------------------------------------------------------------------------
//Description of geometry "Geom1": 
//Geometry for FPOP-Algorithm in 2-dimension. 
//Parameters of geometry: rectangle "rect t" - approximated set, "label_t" - moment of time
//The updated geometry is a rectangle that approximates the intersection of the rectangle and disK.
//Check for emptiness - correct rectangle coordinates. 
//-----------------------------------------------------------------------------------------------
class Geom1{
private:
  unsigned int label_t; //time moment 
  Rect rect_t;          //approx rectangle
public:
  Geom1():label_t (0), rect_t(Rect()){}
  Geom1(unsigned int t):label_t(t), rect_t(Rect()){}
 
  unsigned int get_label_t() const;
  Rect get_rect_t() const;
  std::list<Ellipse> get_Ellipse_t_1() const;
  
  void InitialGeometry(unsigned int i, const std::list<Ellipse> &ellipses);
  void UpdateGeometry(const Ellipse &ellipse_t);
  bool EmptyGeometry();
};
#endif //GEOM1_H
