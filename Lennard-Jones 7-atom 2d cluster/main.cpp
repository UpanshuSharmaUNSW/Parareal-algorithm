#include"algorithm_natoms.hpp"

/*
Initial condition is the stable well (taken from BinderLelievreSimpson15):  one particle at origin (0,0), rest placed on the unit circle \pi/3 apart - (1,0), (-1,0), (1/2,\sqrt{3}/2), (-1/2,\sqrt{3}/2), (1/2,-\sqrt{3}/2), (-1/2,-\sqrt{3}/2)
*/

int main()
{
    gaussien g;

    input I;              // This creates an instance I of class input
    I.load();             // I gets its variables from the input_file -- see input.cpp

    int natoms = 7;
    // int natoms = 2; // For testing

    particle initial_condition [natoms];
    // Assign value to initial_condition
    initial_condition[0].q_x = 0.124867;
    initial_condition[0].q_y = -0.262327;
    initial_condition[0].p_x = 0.164891;
    initial_condition[0].p_y = -0.743497;

    initial_condition[1].q_x = 1.13063;
    initial_condition[1].q_y = -0.10536;
    initial_condition[1].p_x = -0.033897;
    initial_condition[1].p_y = -0.0320293;

    initial_condition[2].q_x = -0.852266;
    initial_condition[2].q_y = -0.428867;
    initial_condition[2].p_x = -0.528037;
    initial_condition[2].p_y = 0.227548;

    initial_condition[3].q_x = 0.396524;
    initial_condition[3].q_y = 0.726056;
    initial_condition[3].p_x = 0.00710503;
    initial_condition[3].p_y = -0.113494;

    initial_condition[4].q_x = -0.560437;
    initial_condition[4].q_y = 0.503632;
    initial_condition[4].p_x = -0.27061;
    initial_condition[4].p_y = 0.25531;

    initial_condition[5].q_x = 0.815286;
    initial_condition[5].q_y = -1.07102;
    initial_condition[5].p_x = 0.0889535;
    initial_condition[5].p_y = -0.0183866;

    initial_condition[6].q_x = -0.202437;
    initial_condition[6].q_y = -1.26057;
    initial_condition[6].p_x = -0.0415541;
    initial_condition[6].p_y = -0.357407;

    particle * initial_atoms = NULL;
    initial_atoms = new particle[natoms];

    for (int i =0 ; i < natoms ; i++)
    {
        initial_atoms[i].q_x=initial_condition[i].q_x;
        initial_atoms[i].q_y=initial_condition[i].q_y;
        initial_atoms[i].p_x=initial_condition[i].p_x;
        initial_atoms[i].p_y=initial_condition[i].p_y;
    }

    Algorithm_langevin * A;   // create pointer to object of class Algorithm_Langevin
    A = new Algorithm_langevin(I,natoms);
    A->H_full=new hamiltonian_full(natoms);
    A->H_cg=new hamiltonian_cg(natoms);

    // A->compute_traj(I, g, initial_atoms, natoms);
    A->compute_traj_AdSlab(I, g, initial_atoms, natoms);
    // A->hessian_comp(initial_atoms, natoms);
    delete A;
    return 0;
}
