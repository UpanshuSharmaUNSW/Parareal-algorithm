#include"hamiltonian.hpp"
#include"gaussien.hpp"

// fine force
void hamiltonian_full::compute_force(double & x)
{

    // ** double well: energy= (x^2-1)^2
    // force = -4.*(x*x-1.)*x;

    // ** single well: energy= x^2/2, force = - energy' = -x
    force = - x;

    // ** Lennard-Jones: energy x^{-12}-2x^{-6}, force = - energy' = 12(x^{-13}-x^{-7})
    // force = 12. * (pow(x,-13) - pow(x,-7));
}

// fine energy
double hamiltonian_full::ener(double & x, double & v)
{
    // ** double well
    // double e = pow(x*x-1.,2) + 0.5*v*v;

    // ** single well
    double e = pow(x,2)/2 + 0.5*v*v;

    // ** Lennard-Jones
    // double e = (pow(x,-12)-2.*pow(x,-6)) + 0.5*v*v;

    return e;
}

// fine kinetic energy
double hamiltonian_full::kin_ener(double & x, double & v)
{
    // ** does not depend on choice of potential
    double e = 0.5*v*v;
    return e;
}

// fine potential energy
double hamiltonian_full::pot_ener(double & x, double & v)
{
    // ** double well
    // double e = pow(x*x-1.,2);

    // ** single well
    double e = pow(x,2)/2;

    // ** Lennard-Jones
    // double e = pow(x,-12)-2.*pow(x,-6);

    return e;
}


// coarse force
void hamiltonian_cg::compute_force(double & x)
{
    // ** double well
    // ** energy: if x>0, 4(x-1)^2; if x<0, 4(x+1)^2;
    // if (x>0)
    // {
    //     force = -8*(x-1.);
    // }
    // else
    // {
    //     force = -8*(x+1.);
    // }

    // ** single well - use coarse potential to be a single well 1/2.x^2
    double w = 0.1;
    force = - w * x;

    // ** harmonic approximation of Lennard-Jones around minimum x=1
    // ** energy: -1+36(x-1)^2, force = -energy'
    // force = 72.*(1. - x);
}
