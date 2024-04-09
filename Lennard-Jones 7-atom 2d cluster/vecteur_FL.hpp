#ifndef VECTEUR_HPP
#define VECTEUR_HPP
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
using namespace std;

///// This is the vecteur class of matrice_FL
///// The independent class belongs to Gabriel

///////////////////////////////////////////////
//// CLASSE vecteur

class vec
{
public:

  int FLdim;            // Dimensions of the vector
  double* FLcoord;      // FLcoord is a pointer to doubles

  vec();
    // default constructor

  vec(int);
  //build vector of given size

  vec(const vec &);
  // copying operator

  ~vec();
  // destructor

  vec & operator= (const vec &);
  // overload assignemnt for =

  vec operator+ (const vec &);
  // overload assignemnt for +

  vec operator/ (const double a);
  //division by a scalar

  double & operator() (int);
  // return a component of the vector

  double & operator[] (int);
  // same as ()
  
  void set_size(int);
  // set current vector to the desired length

  void zeros();
  // make a vector 0
  
  //  vec mult(double);
  // multiplication of the current vector by a scalar and return thr result
  
  double scal(const vec &);
  // scalar product of two vectors

  double norm2();
  // L2 norm/ Cartesian norm of a vector
  
};

vec operator*(double t, const vec &v);
  // scalar multiplication of a vector 'v' by salar 't'

#endif

