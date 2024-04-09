#ifndef HAM_HPP
#define HAM_HPP

#include"particle.hpp"
#include"vecteur_FL.hpp"
#include <math.h>

class hamiltonian
{
public:
    double force;

    virtual ~hamiltonian(){};
    virtual void compute_force(double & ){};
};


class hamiltonian_full : public hamiltonian
{
public:
    hamiltonian_full(){ ; }
    void compute_force(double &);
    void compute_force_adaptive(double &, double );
    double ener(double & x, double & v);
    double kin_ener(double & x, double & v);
    double pot_ener(double & x, double & v);
};


class hamiltonian_cg : public hamiltonian
{
public:
    hamiltonian_cg(){ ; }
    void compute_force(double &);
};

#endif
