#ifndef ALGORITHM_NATOMS_HPP
#define ALGORITHM_NATOMS_HPP

#include "gaussien.hpp"
#include "particle.hpp"
#include "hamiltonian.hpp"
#include "matrice_GS.hpp"
// #include <Eigen/Dense>

// using namespace Eigen;

//#include <gsl/gsl_math.h>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_odeiv2.h>
//#include <gsl/gsl_sf.h>

class Algorithm_langevin
{
public:

    hamiltonian_full * H_full;
    hamiltonian_cg * H_cg;

    mat * noise_c_x; // Coarse noise for variable p_x
    mat * noise_c_y; // Coarse noise for variable p_y
    mat * noise_f_x; // Fine noise for vaiable p_x
    mat * noise_f_y; // Fine noise for variable p_y

    // Default constructor
    Algorithm_langevin(input &, int );

    void construct_noise(input & , gaussien & , int );

    void construct_noise_invariant(input & , gaussien & , int );

    // Fine propogator
    particle * prop_fine(input & , vec * , vec * , particle * , int );
    particle * prop_fine(input & , vec * , vec * , particle * , int , double & , double &, double &);
    particle * prop_fine_inv(input & , vec * , vec * , particle * , int ); // fixes atom 0 to origin and atom 1 to the x-axis

    // Coarse propoagtor
    particle * prop_coarse(input & , vec * , vec *,  particle * , int, particle * well);// coarse-propoagtor that uses harmonic approximation of the initial well (last input)
    particle * prop_coarse_inv(input & , vec * , vec *,  particle * , int, particle * well);// coarse-propoagtor that uses harmonic approximation of the initial well (last input) and fixes atom 0 to origin and atom 1 to x-axis

    // Main parareal routine
    void compute_traj(input & , gaussien & , particle * , int);

    void compute_traj_AdSlab(input & , gaussien & , particle * , int);

    // Well-id routines
    double  well_id(particle * , int , int ); // second int keeps track of the Verlet iteration calling gradient descent algorithm
    double  well_id(particle * , int , int, int ); // last int keeps track of the parareal iteration calling gradient descent algorithm
    particle * well_id(particle * atom_inp, int natoms); // Same as above - returns position of the particles at the well minima
    particle * LineSearch(particle * current_solution, double * current_direction_x, double * current_direction_y, double slope, double mu, int natoms);// Line search algorithm

    // // Compare hessain calculated analytically and using finite difference
    // void hessian_comp(particle *, int ); // Computes the hessian for LJ7 using analytical expression
};

#endif
