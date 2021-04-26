#include "Mat2X2.h"
#include "math.h"

//constructor*******************************************************************
Mat2X2::Mat2X2(){
  mat = new double*[2]; 
  for (unsigned int i = 0; i < 2; i++) {
    mat[i] = new double[2]; 
    mat[i][0] = 0;
    mat[i][1] = 0;
  }
}

Mat2X2::Mat2X2(double a11, double a12, double a21, double a22){
  mat = new double*[2]; 
  for (unsigned int i = 0; i < 2; i++) {mat[i] = new double[2];}
  mat[0][0] = a11;
  mat[0][1] = a12;
  mat[1][0] = a21;
  mat[1][1] = a22;
}

//constructor copy**************************************************************
Mat2X2::Mat2X2(const Mat2X2 & M){
  mat = new double*[2]; 
  for(unsigned int i = 0; i < 2; i++) {
    mat[2] = new double[2];
    mat[i][0] = M.mat[i][0];
    mat[i][1] = M.mat[i][1];
  }
}

//destructor********************************************************************
Mat2X2::~Mat2X2(){
  for(unsigned int i = 0; i < 2; i++) {delete[]mat[i];}
  delete[]mat;
  mat = NULL;
}

//operator = *******************************************************************
Mat2X2 &Mat2X2::operator=(const Mat2X2 &M){
  mat = new double*[2]; 
  for(unsigned int i = 0; i < 2; i++) {
    mat[2] = new double[2];
    mat[i][0] = M.mat[i][0];
    mat[i][1] = M.mat[i][1];
  }
  return * this;
}
//SumMat************************************************************************
Mat2X2 Mat2X2::SumMat(const Mat2X2&M2){
  Mat2X2 res;
  for(unsigned int i = 0; i < 2; i++) {
    res.mat[i][0] = mat[i][0] + M2.mat[i][0];
    res.mat[i][1] = mat[i][1] + M2.mat[i][1];
  }
  return res;
}
//DifMat************************************************************************
Mat2X2 Mat2X2::DifMat(const Mat2X2&M2){
  Mat2X2 res;
  for(unsigned int i = 0; i < 2; i++) {
    res.mat[i][0] = mat[i][0] - M2.mat[i][0];
    res.mat[i][1] = mat[i][1] - M2.mat[i][1];
  }
  return res;
}
//MultNb************************************************************************
Mat2X2 Mat2X2::MultNb(double nb){
  Mat2X2 res;
  for(unsigned int i = 0; i < 2; i++) {
    res.mat[i][0] = nb*mat[i][0];
    res.mat[i][1] = nb*mat[i][1];
  }
  return res;
}

//MultMat***********************************************************************
Mat2X2 Mat2X2::MultMat(const Mat2X2& M2){
  Mat2X2 res;
  for(unsigned int i = 0; i < 2; i++) {
    for(unsigned int j = 0; j < 2; i++) {
      for(unsigned int k = 0; k < 2; i++) {
         res.mat[i][j] = res.mat[i][j] + mat[i][k]*M2.mat[k][j];
      }
    }
  }
  return res;
}
//******************************************************************************
// Mat x vec
double* Mat2X2::MatMultVec(double* & vec){
  double* res = new double[2];
  for (unsigned int i = 0; i < 2; i++){ res[i] = mat[i][0]*vec[0] + mat[i][1]*vec[1];}
  return res;
}
//******************************************************************************
//  vec x Mat 
double* Mat2X2::VecMultMat(double* & vec){
  double* res = new double[2];
  for (unsigned int i =0; i<2; i++){ res[i] = mat[0][i]*vec[0] + mat[1][i]*vec[1];}
  return res;
}

//det***************************************************************************
double Mat2X2::Det(){ 
  double res = mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
  return res;
}

//Obrat*************************************************************************
Mat2X2 Mat2X2::Obrat(){
  Mat2X2 res;
  double det = mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
  res.mat[0][0] = det*(mat[1][1]);
  res.mat[0][1] = -det*(mat[1][0]);
  res.mat[1][0] = -det*(mat[0][1]);
  res.mat[1][1] = det*(mat[0][0]);
  return res;
}

//******************************************************************************
Mat2X2  Mat2X2::Trans(){
  Mat2X2 res;
  for (unsigned int i = 0; i < 2; i++){
    for (unsigned int j = 0; j < 2; j++){res.mat[i][j] = mat[j][i];}
  }
  return res;
}
//******************************************************************************
double  Mat2X2::Trace(){return mat[0][0] +m at[1][1];}
//******************************************************************************
double** Mat2X2::get_mat() const {return mat;} 
double Mat2X2::get_element(unsigned int i, unsigned int j) const {return mat[i][j];}

