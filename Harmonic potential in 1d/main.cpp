#include"algorithm.hpp"

int main(void)
{
    //  srand(time(NULL));
    //srand48(1);
    //  cerr<<RAND_MAX;

    gaussien g;

    input I;
    I.load();

    // initial conditions for the particle are drawn from N(0,1/beta) with beta=3 (see burn-in code with seed 120)
    particle X0(I);
    X0.q = -0.386593;
    X0.p = 0.48453;

    if (I.type_algo > 1)
    {
        cerr<<"I.type_algo > 1: on va comme si I.type_algo == 0 dans le full"<<endl;
    }

    Algorithm_langevin* A;
    A = new Algorithm_langevin(I);
    A->H_full=new hamiltonian_full();
    A->H_cg=new hamiltonian_cg();
    // A->compute_traj(g,I,X0);
    // A->compute_traj_slab(g,I,X0);
    A->compute_traj_AdSlab(g,I,X0);
    // A->compute_traj_slab_fineInit(g,I,X0);
}
