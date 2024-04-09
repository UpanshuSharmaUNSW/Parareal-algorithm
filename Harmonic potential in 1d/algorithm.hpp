#ifndef ALGORITH_HPP
#define ALGORITH_HPP

#include"gaussien.hpp"
#include "particle.hpp"
#include "hamiltonian.hpp"
#include "matrice_GS.hpp"

class Algorithm_langevin
{
public:

    hamiltonian_full * H_full;
    hamiltonian_cg * H_cg;

    mat noise_c;
    mat noise_f;

    // Default constructor
    Algorithm_langevin(input &);

    void construct_noise(input& I, gaussien & g);

    // Fine propogator
    vec prop_fine(input& I, vec & noise, double & q_old, double & p_old);
    vec prop_fine(input& I, vec & noise, double & q_old, double & p_old, double & moy_p2);

    // Coarse propogator
    vec prop_coarse(input& I, vec & noise, double & q_old, double & p_old);

    // Main parareal routine
    void compute_traj(gaussien & g, input &, particle & X0);

    // Main parareal routine with split slab and comparision to fine for entire trajectory
    void compute_traj_AdSlab(gaussien & g, input &, particle & X0);

    // Main parareal routine with split slab and comparision to fine initialised for each slab
    // void compute_traj_slab_fineInit(gaussien & g, input &, particle & X0);
};

#endif
