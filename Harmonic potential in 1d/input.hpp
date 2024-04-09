#ifndef INPUT_HPP
#define INPUT_HPP
#include<iostream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string.h>
#include "vecteur_FL.hpp"
#include "gaussien.hpp"
using namespace std;

///////////////////////////////////////////////////////////////////////////////
/////////////////////////////// INPUT PARAMETERS //////////////////////////////
////////////////// OBTAINED FROM THE INPUT FILE "input_file" //////////////////
///////////////////////////////////////////////////////////////////////////////


class input
{
public:

  int Nb_para; // nb d'iterations parareel
  int N_Coarse_steps; // nb d'intervalles macro
  int N_Fine_steps; // nb de pas de temps dans un intervalle macro
  double t_step;  // pas de temps
  int type_algo; // type de l'algo

  double beta; // inverse de la temperature
  double gamma; // friction dans langevin

  template <class T>
  void read_item(istream& ,const char *, T *); /// reading function for 1 data (int or double)
  // void read_item_string (istream& ,char *, char *); /// reading function for 1 string data
  void read(int & Nb_para_, int & N_Coarse_step_, int & N_Fine_step_, double& t_step_, int& type_algo_, double& beta_, double& gamma_);
  void load(void);	/// copying datas from input_file to the parameters class

};
#endif
