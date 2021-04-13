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
 "(rectx1,recty1)" - the coordinates of the top right vertex.
 -------------------------------------------------------------------------------
 */
class Rect{
private:
  double rectx0;                     
  double recty0;
  double rectx1;                      
  double recty1;
public:
  
  Rect():rectx0(-INFINITY), recty0(-INFINITY), rectx1(INFINITY), recty1(INFINITY){}
  Rect(double x0, double y0, double x1, double y1):rectx0(x0), recty0(y0), rectx1(x1), recty1(y1){}

  double get_rectx0() const;
  double get_recty0() const;
  double get_rectx1() const;
  double get_recty1() const;

  double min_ab(double a, double b);
  double max_ab(double a, double b);
  
  bool EmptyIntersection(const Ellipse &ellipse);

  bool IsEmpty_rect();
  void Exclusion_Ellipse(const Ellipse &ellipse);
  void Intersection_Ellipse(const Ellipse &ellipse);
  

};

#endif //RECT_H