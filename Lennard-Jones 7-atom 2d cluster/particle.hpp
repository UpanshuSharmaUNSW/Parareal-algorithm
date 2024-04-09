#ifndef PARTICLE_HPP
#define PARTICLE_HPP
#include"input.hpp"
#include "vecteur_FL.hpp"

class particle
{
public:
    // position
    double q_x;
    double q_y;
    // momentum
    double p_x;
    double p_y;

    // default constructor - assigns pre-set values to particle
    particle();

    // assigns p=my_p and q=my_q
    particle(double my_q_x, double my_q_y, double my_p_x, double my_p_y);

    // overload assignment
    particle & operator= (const particle &);    
};

#endif
