#include <Rcpp.h>
using namespace Rcpp;

#ifndef MAT2X2_H
#define MAT2X2_H

#include <vector>
/*
 Class Mat2X2
 -------------------------------------------------------------------------------
 Description: 
 
 -------------------------------------------------------------------------------
 */
class Mat2X2{
  
private:
  double **mat;
public:
  
  Mat2X2(); 
  Mat2X2(double a11, double a12, double a21, double a22); 
  ~Mat2X2();
  Mat2X2(const Mat2X2 & M);
  Mat2X2& operator=(const Mat2X2& M);          
  
  Mat2X2 SumMat( const Mat2X2&M2);      
  Mat2X2 DifMat(const Mat2X2&M2);      
  Mat2X2 MultNb(double nb);           
  Mat2X2 MultMat(const Mat2X2& M2);           
  double* MatMultVec(double* & vec);  
  double* VecMultMat(double* & vec);
  
  double Det();
  Mat2X2 Obrat(); 
  Mat2X2 Trans(); 
  double Trace(); 
  
  double** get_mat() const;
  double get_el(unsigned int i, unsigned int j) const;
};

#endif // MAT2X2_H
