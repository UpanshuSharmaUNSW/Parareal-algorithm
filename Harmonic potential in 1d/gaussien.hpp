#ifndef GAUSSIEN_
#define GAUSSIEN_

#include<iostream>
#include "MersenneTwister.hpp"

class gaussien
{
public :
  int random_flag;
  double r1;
  double r2;
  MTRand mtrand1;
  
  gaussien();
  double random_number();
  double rand_number(double,double);
  double rand_number_tronque(double, double, double);
};
#endif




