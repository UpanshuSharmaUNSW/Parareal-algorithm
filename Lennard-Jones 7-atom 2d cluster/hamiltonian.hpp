#ifndef HAM_HPP
#define HAM_HPP

#include"particle.hpp"
#include"vecteur_FL.hpp"
#include <math.h>
#include"matrice_GS.hpp"


class hamiltonian
{
public:
    double * force_x;
    double * force_y;

    virtual ~hamiltonian(){}
    virtual void compute_force(particle *, int ){}
};




class hamiltonian_full : public hamiltonian
{
public:
    hamiltonian_full();  // default constructor
    hamiltonian_full(int ); // makes force_x,force_y into arrays of size int

    void compute_force(particle * , int );    // Calculates the force acting on the full system
    double pot_ener(particle * , int ); // returns the potential energy of system
    double kin_ener(particle * , int ); // returns the kinetic energy of system
    double ener(particle * , int ); // returns the total energy of system
};



class hamiltonian_cg : public hamiltonian
{
public:
    hamiltonian_cg();  // constructor
    hamiltonian_cg(int ); // the force_x,force_y are made arrays

    void compute_force(particle * , int ); // computes force using harmonic approximation of the 1d LJ potential
    void compute_force(particle * , int ,particle *); // computes force using harmonic approximation of the initial well (in natoms*2 dimensions)
};

#endif
