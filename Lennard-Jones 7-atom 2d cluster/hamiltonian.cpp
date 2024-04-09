  #include"hamiltonian.hpp"
#include"matrice_GS.hpp"
#include"gaussien.hpp"
#include <cmath>

/*
macros that calculate powers of a number
*/

#define X2(a)   (a)*(a)
#define X3(a)   X2(a)*(a)
#define X4(a)   X2(a)*X2(a)
#define X6(a)   X4(a)*X2(a)
#define X12(a) X6(a)*X6(a)

////////////////////////////////
//////// Fine routines ////////
///////////////////////////////

// default constructor
hamiltonian_full::hamiltonian_full()
{
    force_x = NULL;
    force_y = NULL;
}

// initialise force_x,force_y arrays to size 'natoms' and value 0
hamiltonian_full::hamiltonian_full(int natoms)
{
    force_x = new double[natoms];
    force_y = new double[natoms];
    for (int i=0 ; i < natoms ; i++)
    {
        force_x[i]=0.0;
        force_y[i]=0.0;
    }
}


/*
Calculates the fine force (-nabla V) on the full system - stored as an array
- r_{ij}_2 = (atom[i].q_x - atom[j].q_x)^2 + (atom[i].q_y - atom[j].q_y)^2
- force_x on i-th particle = \sum_{j\neq i} 6(r_ij^{-12} - r_ij^{-6})*(atom[i].q_x - atom[j].q_x)/(r_ij_2)
- force_y on i-th particle = \sum_{j\neq i} 6(r_ij^{-12} - r_ij^{-6})*(atom[i].q_y - atom[j].q_y)/(r_ij_2)
- inputs: array of particles of size natoms; natoms - number of atoms
*/
void hamiltonian_full::compute_force(particle * atoms, int natoms)
{
    double dx_i; double dy_i;
    double dx_j; double dy_j;
    double r_ij_2; double r_ij;
    double de = 0.0;

    // initialise forces to zero
    for (int i=0 ; i < natoms ; i++)
    {
        force_x[i]=0.0;
        force_y[i]=0.0;
    }

    for (int i=0 ; i < natoms ; i++)
    {
        dx_i = atoms[i].q_x;
        dy_i = atoms[i].q_y;

        for (int j=(i+1) ; j < natoms ; j++)
        {
            dx_j = atoms[j].q_x;
            dy_j = atoms[j].q_y;
            r_ij_2 = X2(dx_i - dx_j) + X2(dy_i - dy_j);
            r_ij = sqrt(X2(dx_i - dx_j) + X2(dy_i - dy_j));

            // Lennard-Jones Potential
            de = 12. * ( 1.0/(X6(r_ij_2)) - 1.0/(X3(r_ij_2)) ) / (r_ij_2);
            force_x[i] += de * (dx_i - dx_j); // bound ij: update forces on atoms #i and #j
            force_y[i] += de * (dy_i - dy_j);
            force_x[j] -= de * (dx_i - dx_j);
            force_y[j] -= de * (dy_i - dy_j);
        }
    }
}

// Compute potential energy for the fine trajectory
double hamiltonian_full::pot_ener(particle * atoms, int natoms)
{
    double dx_i; double dy_i;
    double dx_j; double dy_j;
    double r_ij_2;
    double p_ene = 0.;

    for (int i=0 ; i < natoms; i++)
    {
        dx_i = atoms[i].q_x;
        dy_i = atoms[i].q_y;

        for (int j=(i+1) ; j < natoms ; j++)
        {
            dx_j = atoms[j].q_x;
            dy_j = atoms[j].q_y;
            r_ij_2 = X2(dx_i - dx_j) + X2(dy_i - dy_j);

            //Len-Jon pot = r^{-12} - 2r^{-6}
            p_ene += 1.0/(X6(r_ij_2)) - 2.0/(X3(r_ij_2));
        }

    }

    return p_ene;
}

// Compute kinetic energy for fine trajectory
double hamiltonian_full::kin_ener(particle * atoms, int natoms)
{
    double k_ene = 0.;
    for (int i=0 ; i < natoms; i++)
    {
        k_ene += (X2(atoms[i].p_x))/2.0 + (X2(atoms[i].p_y))/2.0;
    }
    return k_ene;
}


/*
Calculates the total energy of the system ener = k_ene + p_ene
- Inputs: array of particles of size natoms; natoms - number of atoms
- Output: total energy (as a double)
*/
double hamiltonian_full::ener(particle * atoms, int natoms)
{
    double dx_i; double dy_i;
    double dx_j; double dy_j;
    double r_ij_2;
    double p_ene = 0.;
    double k_ene = 0.;
    double ene = 0.;

    //Potential energy
    for (int i=0 ; i < natoms; i++)
    {
        dx_i = atoms[i].q_x;
        dy_i = atoms[i].q_y;

        for (int j=(i+1) ; j < natoms ; j++)
        {
            dx_j = atoms[j].q_x;
            dy_j = atoms[j].q_y;
            r_ij_2 = X2(dx_i - dx_j) + X2(dy_i - dy_j);

            //Len-Jon pot = r^{-12} - 2r^{-6}
            p_ene += 1.0/(X6(r_ij_2)) - 2.0/(X3(r_ij_2));
        }

    }

    //Kinetic energy
    for (int i=0 ; i < natoms; i++)
    {
        k_ene += (X2(atoms[i].p_x))/2.0 + (X2(atoms[i].p_y))/2.0;
    }

    // Full energy
    ene = p_ene + k_ene;

    return ene;
}


/////////////////////////////////
//////// Coarse routines ////////
/////////////////////////////////

// default constructor
hamiltonian_cg::hamiltonian_cg()
{
    force_x = NULL;
    force_y = NULL;
}

// initialise force_x,force_y arrays to size 'natoms' and value 0
hamiltonian_cg::hamiltonian_cg(int natoms)
{
    force_x = new double[natoms];
    force_y = new double[natoms];
    for (int i=0 ; i < natoms ; i++)
    {
        force_x[i]=0.0;
        force_y[i]=0.0;
    }
}


/*
Calculates CG force (-nabla V_c) using Harmonic approximation of 1d LJ potential
// r_{ij} = sqrt((atom[i].q_x - atom[j].q_x)^2 + (atom[i].q_y - atom[j].q_y)^2)
// Two choices of coarse potential: (1) Harmonic approximation and (2) Harmonic approximation with cutoff
// (1)  \phi(r) = -1+72(r-1)^2
// (2) \phi(r) = { -1+72(r-1)^2 when r < 1+1/sqrt(72} and { 0 else }
//  \phi'(r) = { 72*2(r-1)} and {0}
// force_x on i-th particle = { -\sum_{j\neq i} 72(1-1/r_ij)*(atom[i].q_x - atom[j].q_x)} and {0}
// force_y on i-th particle = { -\sum_{j\neq i} 72(1-1/r_ij)*(atom[i].q_y - atom[j].q_y)} and {0}
// inputs: array of particles of size natoms; natoms - number of atoms
*/
void hamiltonian_cg::compute_force(particle * atoms, int natoms)
{
    double dx_i; double dy_i;
    double dx_j; double dy_j;
    double r_ij; double r_ij_2;
    double de = 0.0;

    // initialise to 0
    for (int i=0 ; i < natoms ; i++)
    {
        force_x[i]=0.0;
        force_y[i]=0.0;
    }

    for (int i=0 ; i < natoms ; i++)
    {
        dx_i = atoms[i].q_x;
        dy_i = atoms[i].q_y;

        for (int j=0 ; j < natoms ; j++)
        {
            if (i==j) continue;

            dx_j = atoms[j].q_x;
            dy_j = atoms[j].q_y;
            r_ij = sqrt(X2(dx_i - dx_j) + X2(dy_i - dy_j));
            r_ij_2 = X2(dx_i - dx_j) + X2(dy_i - dy_j);

            //// harmonic approximation
            de = 72. * (1.0/r_ij - 1.0);

            //// harmonic approximation + r^{-12} repulsive
            // de = 72. * (1.0/r_ij - 1.0) + 6 * (1./( X12(r_ij) * X2(r_ij) ) );

            // //// cubic approximation for r<=1 and harmonic app. for r>1
            // if (r_ij <= 1.)
            // {
            //     de = 72. * (1.0/r_ij - 1.0) + 756 * 3 * X2(r_ij-1) * (1.0/r_ij);
            // }
            // else
            // {
            //     de = 72. * (1.0/r_ij - 1.0);
            // }

            //// harmonic approximation with cutoff
            // if (r_ij < 1.+ 1./(sqrt(72.)))
            // {
            //     de = 72. * (1.0/r_ij - 1.0);
            // }
            // else
            // {
            //     de = 0.;
            // }

            // sanity check - use fine potential
            // de = 12. * ( 1.0/(X6(r_ij_2)) - 1.0/(X3(r_ij_2)) ) / (r_ij_2);

            force_x[i] += de * (dx_i - dx_j);
            force_y[i] += de * (dy_i - dy_j) ;
        }
    }
}

/*
Computes CG force using harmonic approximation of the initial well (in natoms*2 dimensions)
*/
void hamiltonian_cg::compute_force(particle * atoms, int natoms, particle * well)
{
    double dx_i = 0.; double dy_i = 0.;
    double dx_ref_i = 0.; double dy_ref_i = 0.;
    double dx_j = 0.; double dy_j = 0.;
    double dx_l = 0.; double dy_l = 0.;
    double r_ij_2 = 0.;
    double r_il_2 = 0.;
    double theta_ij = 0.; double theta_r_ij = 0.;
    double theta_il = 0.; double theta_r_il = 0.;

    particle atom_min [natoms];
    for ( int i = 0 ; i < natoms ; i++)
    {
        atom_min[i] = well[i];
    }

    // initialise the forces to zero
    for (int i=0 ; i < natoms ; i++)
    {
        force_x[i]=0.0;
        force_y[i]=0.0;
    }

    /////////////////////
    // compute hessian //
    /////////////////////

    mat hessian_xx,hessian_yy,hessian_xy;
    hessian_xx.set_size(natoms,natoms);
    hessian_yy.set_size(natoms,natoms);
    hessian_xy.set_size(natoms,natoms);
    for (int i=0 ; i < natoms ; i++)
    {
        for (int l = 0 ; l < natoms ; l++)
        {
            hessian_xx(i,l) = 0.;
            hessian_yy(i,l) = 0.;
            hessian_xy(i,l) = 0.;
        }
    }


    for (int i=0 ; i < natoms ; i++)
    {
        dx_i = atom_min[i].q_x;
        dy_i = atom_min[i].q_y;

        for (int l = 0 ; l < natoms ; l++)
        {
            if (l==i)
            {
                for (int j = 0 ; j < natoms ; j++)
                {
                    if (j!=i)
                    {
                        dx_j = atom_min[j].q_x;
                        dy_j = atom_min[j].q_y;

                        r_ij_2 = X2(dx_i - dx_j)+X2(dy_i - dy_j);
                        theta_ij = - (12./r_ij_2) * (1./(X6(r_ij_2)) - 1./(X3(r_ij_2))); // theta(r_ij)
                        theta_r_ij = (1./(X2(r_ij_2))) * ( (12.*14.)/(X6(r_ij_2)) - (12.*8.)/(X3(r_ij_2)) );   // theta'(r_ij)/r_ij

                        hessian_xx(i,i) += theta_r_ij*(dx_i - dx_j)*(dx_i - dx_j) + theta_ij;
                        hessian_yy(i,i) += theta_r_ij*(dy_i - dy_j)*(dy_i - dy_j) + theta_ij;
                        hessian_xy(i,i) += theta_r_ij*(dx_i - dx_j)*(dy_i - dy_j);
                    }
                }
            }
            else    // l \neq i
            {
                dx_l = atom_min[l].q_x;
                dy_l = atom_min[l].q_y;

                r_il_2 = X2(dx_i - dx_l)+X2(dy_i - dy_l);

                theta_il = - (12./r_il_2) * (1./(X6(r_il_2)) - 1./(X3(r_il_2))); // theta(r_il)
                theta_r_il = (1./(X2(r_il_2))) * ( (12.*14.)/(X6(r_il_2)) - (12.*8.)/(X3(r_il_2)) );   // theta'(r_il)/r_il

                hessian_xx(i,l) += theta_r_il*(dx_l - dx_i)*(dx_i - dx_l) - theta_il;
                hessian_yy(i,l) += theta_r_il*(dy_l - dy_i)*(dy_i - dy_l) - theta_il;
                hessian_xy(i,l) += theta_r_il*(dx_l - dx_i)*(dy_i - dy_l);
            }
        }
    }

    /////////////////////////////////
    // compute force using hessian //
    /////////////////////////////////

    for (int l = 0 ; l < natoms ; l++)
    {
        for (int i=0 ; i < natoms ; i++)
        {
            dx_i = atoms[i].q_x;
            dy_i = atoms[i].q_y;
            dx_ref_i = atom_min[i].q_x;
            dy_ref_i = atom_min[i].q_y;
            force_x[l] -= (dx_i - dx_ref_i)*hessian_xx(i,l) + (dy_i - dy_ref_i)*hessian_xy(i,l);
            force_y[l] -= (dx_i - dx_ref_i)*hessian_xy(i,l) + (dy_i - dy_ref_i)*hessian_yy(i,l);
        }
    }
}
