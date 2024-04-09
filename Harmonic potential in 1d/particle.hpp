#ifndef PARTICLE_HPP
#define PARTICLE_HPP
#include"input.hpp"
#include "vecteur_FL.hpp"
class particle
{
public:
  double q;  // position
  double p; // impulsion
  
  particle(){};
  particle(const input & I);
  particle & operator= (const particle &);
  // affectation surchargée
  
};

#endif

