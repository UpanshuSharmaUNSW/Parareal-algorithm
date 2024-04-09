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

///// On reprend ici la classe vecteur de matrice_FL
///// La classe matrice deveient independante, c'est celle de Gabriel

///////////////////////////////////////////////
//// CLASSE vecteur

class vec
{
public:

  int FLdim;
  double* FLcoord;

  vec();
  // constructeur par defaut

  vec(int);
  //construit un vec de taille donnée

  vec(const vec &);
  // operateur de recopie

  ~vec();
  // destruction

  vec & operator= (const vec &);
  // affectation surchargée

  vec operator+ (const vec &);
  // addition surchargée

  vec operator/ (const double a);
  //division par un scalaire

  double & operator() (int);
  // renvoie la composante d'un vec

  double & operator[] (int);
  // idem que ()
  
  void set_size(int);
  // met le vec courant à la longueur souhaitée

  void zeros();
  // met un vecteur a 0
  
  //  vec mult(double);
  // multiplication du vec courant par un scalaire
  
  double scal(const vec &);
  // produit scalaire

  double norm2();
  // norme au carre du vecteur
  
};

vec operator*(double t, const vec &v);
  // multiplication scalaire par vecteur

#endif

