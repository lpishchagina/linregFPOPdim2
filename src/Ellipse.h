#ifndef ELLIPSE_H
#define ELLIPSE_H

#include <vector>
/*
 Class Ellipse
 -------------------------------------------------------------------------------
 Description: 
 Disk in 2-dimension. 
 
 Parameters:
 "(center1, center2)"  - the disk  center coordinates;
 "radius" - the value of the disk radius.
 -------------------------------------------------------------------------------
 */
class Ellipse
{
private:
  double center1;                           
  double center2;  
  double radius;                                     

public:
  Ellipse(){};
  Ellipse(double c1, double c2, double r):center1(c1), center2(c2), radius(r){}  

  double get_radius() const;
  double get_center1() const;
  double get_center2() const;
}; 

#endif //ELLIPSE_H
//------------------------------------------------------------------------------