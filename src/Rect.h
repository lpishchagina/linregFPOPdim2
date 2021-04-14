#ifndef RECT_H
#define RECT_H

#include "math.h"
#include "Ellipse.h"

#include <vector>


/*
 Class Rect
 -------------------------------------------------------------------------------
 Description: 
 Rectangle in 2-dimension. 
 
 Parameters:
 "(rectx0,recty0)" - the coordinates of the bottom left vertex;
 "(dx,dy)".
 -------------------------------------------------------------------------------
 */
class Rect{
private:
  double x0;                     
  double y0;
  double dx;                      
  double dy;
public:
  
  Rect():x0(-INFINITY), y0(-INFINITY), dx(INFINITY), dy(INFINITY){}
  Rect(double rx0, double ry0, double rdx, double rdy):x0(rx0), y0(ry0), dx(rdx), dy(rdy){}

  double get_x0() const;
  double get_y0() const;
  double get_dx() const;
  double get_dy() const;

  double min_ab(double a, double b);
  double max_ab(double a, double b);
  
  bool EmptyIntersection(const Ellipse &ellipse);

  bool IsEmpty_rect();
  void Exclusion_Ellipse(const Ellipse &ellipse);
  void Intersection_Ellipse(const Ellipse &ellipse);
  

};

#endif //RECT_H