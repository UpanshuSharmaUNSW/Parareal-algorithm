#ifndef MATRICE_HPP
#define MATRICE_HPP
#include <iostream>
#include <fstream>
#include "vecteur_FL.hpp"

using namespace std;

class mat
{
public:

  // champs
  int row, col;
  double * GScoord;

  // constructeurs, destructeurs
  mat();
  ~mat();
  mat(int,int);
  mat(const mat &);
  double & operator() (int,int);
  void set_size(int,int);
  void zeros();
  
  // operations surchargees
  mat & operator= (const mat &);
  mat operator+ (const mat &);
  mat operator- (const mat &);
  mat operator* (const mat &);
  mat operator* (const double); 
  vec mult(const vec &);
  double trace();
  
};

#endif

