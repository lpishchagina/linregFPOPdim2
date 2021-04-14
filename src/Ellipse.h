#ifndef ELLIPSE_H
#define ELLIPSE_H

#include <vector>
/*
 Class Ellipse
 -------------------------------------------------------------------------------
 Description: 
 Disk in 2-dimension. 
 
 Parameters:
------------------------------------
 */
class Ellipse
{
private:
  double a;                           
  double k;  
  double dif;                                     

public:
  Ellipse(){};
  Ellipse(double aa, double kk, double diff):a(aa), k(kk), dif(diff){}  

  double get_a() const;
  double get_k() const;
  double get_dif() const;
}; 

#endif //ELLIPSE_H
//------------------------------------------------------------------------------