#include "algorithm_natoms.hpp"
#include "gaussien.hpp"
// #include <Eigen/Dense>
// #include <eigen3/Eigen/Dense>

/*
macros for taking powers
*/
#define X2(a)   (a)*(a)
#define X3(a)   X2(a)*(a)
#define X4(a)   X2(a)*X2(a)
#define X6(a)   X4(a)*X2(a)
#define X12(a) X6(a)*X6(a)

// The numerical schemes used here are taken from Joubaud-Stoltz-12(MMS, Sec 4.2, page 11)

/*
Initialise object from class Algorithm_langevin
Input - (1) input file with preset parameters, (2) # of atoms
Output - Set the size of the noise matrices
*/
Algorithm_langevin::Algorithm_langevin(input& I, int natoms)
{
    H_full = NULL;
    H_cg   = NULL;

    noise_c_x = new mat[natoms];
    noise_c_y = new mat[natoms];
    noise_f_x = new mat[natoms];
    noise_f_y = new mat[natoms];

    for (int i=0 ; i < natoms ; i++)
    {
        noise_c_x[i].set_size(I.N_Coarse_steps,I.N_Fine_steps);
        noise_c_y[i].set_size(I.N_Coarse_steps,I.N_Fine_steps);
        noise_f_x[i].set_size(I.N_Coarse_steps,I.N_Fine_steps);
        noise_f_y[i].set_size(I.N_Coarse_steps,I.N_Fine_steps);
    }
}


/*
Create noise matrices
Input - (1) input file with preset parameters, (2)object of class gaussien, (3) # of atoms
Output - (1) fill the coarse noise matrix using 'g', (2) fine noise = coarse noise
*/
void Algorithm_langevin::construct_noise(input & I, gaussien & g, int natoms)
{
    for(int j=0 ; j < I.N_Coarse_steps; j++)
    {
        for(int k=0 ; k < I.N_Fine_steps ; k++)
        {
            for (int i=0 ; i < natoms ; i++)
            {
                noise_c_x[i](j,k) = g.rand_number(0.,1.); // draws a gaussian with mean=0 and variance=1
                noise_c_y[i](j,k) = g.rand_number(0.,1.);

                // use same noise for coarse and fine dynamics
                noise_f_x[i](j,k) = noise_c_x[i](j,k);
                noise_f_y[i](j,k) = noise_c_y[i](j,k);
            }
        }
    }
}

/*
Create noise matrices with zero mean in x and y direction. This implies that the systems is translation invariant.
Input - (1) input file with preset parameters, (2)object of class gaussien, (3) # of atoms
Output - (1) fill the coarse noise matrix using 'g', (2) fine noise = coarse noise
*/
void Algorithm_langevin::construct_noise_invariant(input & I, gaussien & g, int natoms)
{
    // generate noise
    for(int j=0 ; j < I.N_Coarse_steps; j++)
    {
        for(int k=0 ; k < I.N_Fine_steps ; k++)
        {
            for (int i=0 ; i < natoms ; i++)
            {
                noise_c_x[i](j,k) = g.rand_number(0.,1.); // draws a gaussian with mean = 0 and variance = 1
                noise_c_y[i](j,k) = g.rand_number(0.,1.);
            }// end i<natoms
        }//end k < I.N_Fine_steps
    }// end j < I.N_Coarse_steps

    // make noise zero mean in x,y direction
    for(int j=0 ; j < I.N_Coarse_steps; j++)
    {
        for(int k=0 ; k < I.N_Fine_steps ; k++)
        {
            double x_mean = 0.;
            double y_mean = 0.;
            // calculate x,y mean
            for(int i=0 ; i < natoms ; i++)
            {
                x_mean += noise_c_x[i](j,k);
                y_mean += noise_c_y[i](j,k);
            }// end i<natoms

            for (int i=0 ; i < natoms ; i++)
            {
                noise_c_x[i](j,k) = noise_c_x[i](j,k) - (1./7.) * x_mean;
                noise_c_y[i](j,k) = noise_c_y[i](j,k) - (1./7.) * y_mean;

                // use same noise for coarse and fine dynamics
                noise_f_x[i](j,k) = noise_c_x[i](j,k);
                noise_f_y[i](j,k) = noise_c_y[i](j,k);
            }// end i<natoms
        }//end k < I.N_Fine_steps
    }// end j < I.N_Coarse_steps
}


/*
Fine Langevin propagator
Input - (1) input file ;  (2) 'noise_x' is a pointer to an array of class 'vec' of size 'natoms' - i-th element in the array and k-th row in the vec contains the noise for the i-th particle and k-th row, i.e. noise_x[i](j) = noise_f_x [i] (j,l) for 0<= l < N_Fine_steps - in the current input file 'N_Fine_steps=1' ; (3) same as (2) ; (4) 'atom_inp' is a pointer to an array of class particle of size 'natoms' and contains the inout configuration, (5) natoms - # of atoms
Output - array of size 'natoms' of class 'particle'
Remark - To avoid mistakes with pointers, I have created new variables in this function. This makes this function memory inefficient.
*/
particle * Algorithm_langevin::prop_fine(input& I, vec * noise_x, vec * noise_y, particle * atom_inp, int natoms)
{
    double dt = I.t_step;
    particle * atom_out = NULL;
    atom_out = new particle[natoms];

    double alpha = exp(-I.gamma*dt);

    particle atom_old [natoms];
    for(int i=0 ; i < natoms ; i++)
    {
        atom_old[i] = atom_inp[i];
    }

    particle atom_new [natoms];

    if (I.type_algo == 0)
    {
        for(int step = 0; step < I.N_Fine_steps ; step++)
        {
            // verlet:
            H_full->compute_force(atom_old,natoms);

            for(int i=0; i < natoms; i++)
            {
                //p_x^{n+1/2}, p_y^{n+1/2}
                atom_new[i].p_x = atom_old[i].p_x + (dt/2.) * H_full->force_x[i];
                atom_new[i].p_y = atom_old[i].p_y + (dt/2.) * H_full->force_y[i];

                //q_x^{n+1}, q_y^{n+1}
                atom_new[i].q_x = atom_old[i].q_x + dt * atom_new[i].p_x;
                atom_new[i].q_y = atom_old[i].q_y + dt * atom_new[i].p_y;
            }

            // recalculate the force
            H_full->compute_force(atom_new,natoms);

            for(int i=0; i < natoms; i++)
            {
                //p_x^{n+1}, p_y^{n+1}
                atom_new[i].p_x += (dt/2.) * H_full->force_x[i];
                atom_new[i].p_y += (dt/2.) * H_full->force_y[i];

                // OU part:
                atom_new[i].p_x = atom_new[i].p_x * alpha + sqrt((1.-alpha*alpha)/I.beta) * noise_x[i](step);
                atom_new[i].p_y = atom_new[i].p_y * alpha + sqrt((1.-alpha*alpha)/I.beta) * noise_y[i](step);

                // we are done -- copy into atom_old
                atom_old[i].q_x = atom_new[i].q_x;
                atom_old[i].q_y = atom_new[i].q_y;
                atom_old[i].p_x = atom_new[i].p_x;
                atom_old[i].p_y = atom_new[i].p_y;
            }
        }
    }

    if (I.type_algo != 0)
    {
        cerr<<"This scheme is not programmed"<<endl;
        exit(-1);
    }

    for(int i = 0; i < natoms ; i++)
    {
        atom_out[i]=atom_new[i];
    }

    return atom_out;
}


/*
Same routine as above -- we fix the atom 0 to the origin and atom 1 to the x-axis
// Here inv in routine title stands for rotation-translation invariant
*/
particle * Algorithm_langevin::prop_fine_inv(input& I, vec * noise_x, vec * noise_y, particle * atom_inp, int natoms)
{
    double dt = I.t_step;
    particle * atom_out = NULL;
    atom_out = new particle[natoms];

    double alpha = exp(-I.gamma*dt);

    particle atom_old [natoms];
    for(int i=0 ; i < natoms ; i++)
    {
        atom_old[i] = atom_inp[i];
    }

    particle atom_new [natoms];

    if (I.type_algo == 0)
    {
        for(int step = 0; step < I.N_Fine_steps ; step++)
        {
            // verlet:
            H_full->compute_force(atom_old,natoms);

            for(int i=0; i < natoms; i++)
            {
                if(i==0) // atom 0 is fixed at the origin
                {
                    atom_new[i].q_x = 0.0; atom_new[i].p_x = 0.0;
                    atom_new[i].q_y = 0.0; atom_new[i].p_y = 0.0;
                }
                else if(i==1) // atom 1 only moves in the x-direction and its y-components are 0
                {
                    atom_new[i].p_y = 0.0; atom_new[i].q_y = 0.0;
                    atom_new[i].p_x = atom_old[i].p_x + (dt/2.) * H_full->force_x[i];
                    atom_new[i].q_x = atom_old[i].q_x + dt * atom_new[i].p_x;
                }
                else
                {
                    //p_x^{n+1/2}, p_y^{n+1/2}
                    atom_new[i].p_x = atom_old[i].p_x + (dt/2.) * H_full->force_x[i];
                    atom_new[i].p_y = atom_old[i].p_y + (dt/2.) * H_full->force_y[i];

                    //q_x^{n+1}, q_y^{n+1}
                    atom_new[i].q_x = atom_old[i].q_x + dt * atom_new[i].p_x;
                    atom_new[i].q_y = atom_old[i].q_y + dt * atom_new[i].p_y;
                }
            }

            // recalculate the force
            H_full->compute_force(atom_new,natoms);

            for(int i=0; i < natoms; i++)
            {
                //p_x^{n+1}, p_y^{n+1}
                if(i==0)
                {
                    atom_new[i].p_x=0.0; atom_new[i].p_y=0.0;
                }
                else if (i==1)
                {
                    atom_new[i].p_y=0.0;
                    atom_new[i].p_x += (dt/2.) * H_full->force_x[i];
                    atom_new[i].p_x = atom_new[i].p_x * alpha + sqrt((1.-alpha*alpha)/I.beta) * noise_x[i](step);
                }
                else
                {
                    atom_new[i].p_x += (dt/2.) * H_full->force_x[i];
                    atom_new[i].p_y += (dt/2.) * H_full->force_y[i];

                    // OU part:
                    atom_new[i].p_x = atom_new[i].p_x * alpha + sqrt((1.-alpha*alpha)/I.beta) * noise_x[i](step);
                    atom_new[i].p_y = atom_new[i].p_y * alpha + sqrt((1.-alpha*alpha)/I.beta) * noise_y[i](step);
                }

                // we are done -- copy into atom_old
                atom_old[i].q_x = atom_new[i].q_x;
                atom_old[i].q_y = atom_new[i].q_y;
                atom_old[i].p_x = atom_new[i].p_x;
                atom_old[i].p_y = atom_new[i].p_y;
            }
        }
    }

    if (I.type_algo != 0)
    {
        cerr<<"This scheme is not programmed"<<endl;
        exit(-1);
    }

    for(int i = 0; i < natoms ; i++)
    {
        atom_out[i]=atom_new[i];
    }

    return atom_out;
}


/*
Fine Langevin propagator (overloaded - same as above) -- to calculate additional objects
Input - (1) input file;  (2) 'noise_x' is a pointer to an array of class 'vec' of size 'natoms' - i-th element in the array and k-th row in the vec contains the noise for the i-th particle and k-th row, i.e. noise_x[i](j) = noise_f_x [i] (j,l) for 0<= l < N_Fine_steps - in the current input file 'N_Fine_steps=1' ; (3) same as (2) ; (4) 'atom_inp' is a pointer to an array of class particle of size 'natoms' and contains the inout configuration, (5) natoms - # of atoms
Output - array of size 'natoms' of class 'particle'
*/
particle * Algorithm_langevin::prop_fine(input& I, vec * noise_x, vec * noise_y, particle * atom_inp, int natoms, double & aver_p2, double & aver_p4, double & aver_config_temp)
{
    double dt = I.t_step;
    particle * atom_out = NULL;
    atom_out = new particle[natoms];

    double alpha = exp(-I.gamma*dt);

    particle atom_old [natoms];
    for(int i=0 ; i < natoms ; i++)
    {
        atom_old[i] = atom_inp[i];
    }

    particle atom_new [natoms];

    if (I.type_algo == 0)
    {
        for(int step = 0; step < I.N_Fine_steps ; step++)
        {
            // verlet:
            H_full->compute_force(atom_old,natoms);

            for(int i=0; i < natoms; i++)
            {
                //p_x^{n+1/2}, p_y^{n+1/2}
                atom_new[i].p_x = atom_old[i].p_x + (dt/2.) * H_full->force_x[i];
                atom_new[i].p_y = atom_old[i].p_y + (dt/2.) * H_full->force_y[i];

                //q_x^{n+1}, q_y^{n+1}
                atom_new[i].q_x = atom_old[i].q_x + dt * atom_new[i].p_x;
                atom_new[i].q_y = atom_old[i].q_y + dt * atom_new[i].p_y;
            }

            // recalculate the force
            H_full->compute_force(atom_new,natoms);

            for(int i=0; i < natoms; i++)
            {
                //p_x^{n+1}, p_y^{n+1}
                atom_new[i].p_x += (dt/2.) * H_full->force_x[i];
                atom_new[i].p_y += (dt/2.) * H_full->force_y[i];

                // OU part:
                atom_new[i].p_x = atom_new[i].p_x * alpha + sqrt((1.-alpha*alpha)/I.beta) * noise_x[i](step);
                atom_new[i].p_y = atom_new[i].p_y * alpha + sqrt((1.-alpha*alpha)/I.beta) * noise_y[i](step);

                // we are done -- copy in atom_old
                atom_old[i] = atom_new[i];

                // Calculate additional objects
                // aver_p2 += atom_old[i].p_x * atom_old[i].p_x + atom_old[i].p_y * atom_old[i].p_y;
                aver_p2 += X2(atom_old[i].p_x) + X2(atom_old[i].p_y);
                aver_p4 += X4(atom_old[i].p_x) + X4(atom_old[i].p_y);
                aver_config_temp -= atom_old[i].q_x * H_full->force_x[i];
                aver_config_temp -= atom_old[i].q_y * H_full->force_y[i];
            }
        }
    }

    if (I.type_algo != 0)
    {
        cerr<<"This scheme is not programmed"<<endl;
        exit(-1);
    }

    for(int i = 0; i < natoms ; i++)
    {
        atom_out[i]=atom_new[i];
    }

    return atom_out;
}




/*
Coarse Langevin propagator which uses harmonic approximation of initial well
Inputs - (1) same as above, (2) additional input of initial well
*/
particle * Algorithm_langevin::prop_coarse(input& I, vec * noise_x, vec * noise_y, particle * atom_inp, int natoms, particle * well)
{
    double dt = I.t_step;
    particle * atom_out = NULL;
    atom_out = new particle[natoms];

    double alpha = exp(-I.gamma*dt);

    particle atom_old [natoms];
    for(int i=0 ; i < natoms ; i++)
    {
        atom_old[i] = atom_inp[i];
    }

    particle atom_new [natoms];

    if (I.type_algo == 0)
    {
        for(int step = 0; step < I.N_Fine_steps ; step++)
        {
            // verlet:
            H_cg->compute_force(atom_old,natoms,well);  // harmonic approximation around initial well

            for(int i=0; i < natoms; i++)
            {
                //p_x^{n+1/2}, p_y^{n+1/2}
                atom_new[i].p_x = atom_old[i].p_x + (dt/2.) * H_cg->force_x[i];
                atom_new[i].p_y = atom_old[i].p_y + (dt/2.) * H_cg->force_y[i];

                //q_x^{n+1}, q_y^{n+1}
                atom_new[i].q_x = atom_old[i].q_x + dt * atom_new[i].p_x;
                atom_new[i].q_y = atom_old[i].q_y + dt * atom_new[i].p_y;
            }

            // recalculate the force
            H_cg->compute_force(atom_new,natoms,well); // harmonic approximation around initial well

            for(int i=0; i < natoms; i++)
            {
                //p_x^{n+1}, p_y^{n+1}
                atom_new[i].p_x += (dt/2.) * H_cg->force_x[i];
                atom_new[i].p_y += (dt/2.) * H_cg->force_y[i];

                // OU part:
                atom_new[i].p_x = atom_new[i].p_x * alpha + sqrt((1.-alpha*alpha)/I.beta) * noise_x[i](step);
                atom_new[i].p_y = atom_new[i].p_y * alpha + sqrt((1.-alpha*alpha)/I.beta) * noise_y[i](step);

                atom_old[i] = atom_new[i];
            }
        }
    }

    if (I.type_algo != 0)
    {
        cerr<<"This scheme is not programmed"<<endl;
        exit(-1);
    }

    for(int i = 0; i < natoms ; i++)
    {
        atom_out[i]=atom_new[i];
    }

    return atom_out;
}

/*
Coarse Langevin propagator which uses harmonic approximation of initial well (same as above), and fixes atom 0 to origin and atom 1 to the x-axis
*/
particle * Algorithm_langevin::prop_coarse_inv(input& I, vec * noise_x, vec * noise_y, particle * atom_inp, int natoms, particle * well)
{
    double dt = I.t_step;
    particle * atom_out = NULL;
    atom_out = new particle[natoms];

    double alpha = exp(-I.gamma*dt);

    particle atom_old [natoms];
    for(int i=0 ; i < natoms ; i++)
    {
        atom_old[i] = atom_inp[i];
    }

    particle atom_new [natoms];

    if (I.type_algo == 0)
    {
        for(int step = 0; step < I.N_Fine_steps ; step++)
        {
            // verlet:
            H_cg->compute_force(atom_old,natoms,well);  // harmonic approximation around initial well

            for(int i=0; i < natoms; i++)
            {
                if(i==0) // atom 0 is fixed at the origin
                {
                    atom_new[i].p_x = 0.; atom_new[i].p_y = 0.;
                    atom_new[i].q_x = 0.; atom_new[i].q_y = 0.;
                }
                else if(i==1) // atom 1 only moves in the x-direction and its y-components are 0
                {
                    atom_new[i].p_y = 0.; atom_new[i].q_y = 0.;
                    atom_new[i].p_x = atom_old[i].p_x + (dt/2.) * H_cg->force_x[i];
                    atom_new[i].q_x = atom_old[i].q_x + dt * atom_new[i].p_x;
                }
                else
                {
                    //p_x^{n+1/2}, p_y^{n+1/2}
                    atom_new[i].p_x = atom_old[i].p_x + (dt/2.) * H_cg->force_x[i];
                    atom_new[i].p_y = atom_old[i].p_y + (dt/2.) * H_cg->force_y[i];

                    //q_x^{n+1}, q_y^{n+1}
                    atom_new[i].q_x = atom_old[i].q_x + dt * atom_new[i].p_x;
                    atom_new[i].q_y = atom_old[i].q_y + dt * atom_new[i].p_y;
                }
            }

            // recalculate the force
            H_cg->compute_force(atom_new,natoms,well);  // harmonic approximation around initial well

            for(int i=0; i < natoms; i++)
            {
                //p_x^{n+1}, p_y^{n+1}
                if(i==0)
                {
                    atom_new[i].p_x=0.; atom_new[i].p_y=0.;
                }
                else if (i==1)
                {
                    atom_new[i].p_y=0.;
                    atom_new[i].p_x += (dt/2.) * H_cg->force_x[i];
                    atom_new[i].p_x = atom_new[i].p_x * alpha + sqrt((1.-alpha*alpha)/I.beta) * noise_x[i](step);
                }
                else
                {
                    atom_new[i].p_x += (dt/2.) * H_cg->force_x[i];
                    atom_new[i].p_y += (dt/2.) * H_cg->force_y[i];

                    // OU part:
                    atom_new[i].p_x = atom_new[i].p_x * alpha + sqrt((1.-alpha*alpha)/I.beta) * noise_x[i](step);
                    atom_new[i].p_y = atom_new[i].p_y * alpha + sqrt((1.-alpha*alpha)/I.beta) * noise_y[i](step);
                }
                atom_old[i] = atom_new[i];
            }
        }
    }

    if (I.type_algo != 0)
    {
        cerr<<"This scheme is not programmed"<<endl;
        exit(-1);
    }

    for(int i = 0; i < natoms ; i++)
    {
        atom_out[i]=atom_new[i];
    }

    return atom_out;
}

/*
//// *** MAIN ROUTINE *** ////
Compute parareal and reference Langevin trajectory (and well-id)
Input - (1) input file;  (2) object of class 'gaussien' ; (3) 'atom_inp' is a pointer to an array of class particle of size 'natoms' and contains the input configuration, (4) natoms - # of atoms
*/
void Algorithm_langevin::compute_traj(input & I, gaussien & g, particle * atom_inp, int natoms)
{
    particle atom_cur [natoms];
    for(int i=0 ; i < natoms ; i++)
    {
        atom_cur[i] = atom_inp[i];
    }

    particle * well_init;
    well_init = new particle [natoms];
    // well_init = well_id (atom_inp , natoms); well corresponing to initial condition
    for (int i = 0 ; i < natoms ; i++)
    {
        well_init[i] = atom_inp[i];
    }


    // (x,y) position and momentum matrices for the i-th particle, with k-th row corresponding to k-th parareal iteration
    mat all_q_x [natoms];
    mat all_q_y [natoms];
    mat all_p_x [natoms];
    mat all_p_y [natoms];
    //Initialisation
    for(int i=0 ; i < natoms ; i++)
    {
        all_q_x[i].set_size(I.Nb_para,1+I.N_Coarse_steps);
        all_q_y[i].set_size(I.Nb_para,1+I.N_Coarse_steps);
        all_p_x[i].set_size(I.Nb_para,1+I.N_Coarse_steps);
        all_p_y[i].set_size(I.Nb_para,1+I.N_Coarse_steps);
    }

    mat err(I.Nb_para,1+I.N_Coarse_steps);
    vec max_err(I.Nb_para);
    mat energy(I.Nb_para,1+I.N_Coarse_steps);    // stores energy of the system at each time-step

    // initial data (t=0) for every parareal iteration is given
    for (int k=0 ; k < I.Nb_para ; k++)
    {
        for (int i=0 ; i < natoms ; i++)
        {
            all_q_x[i](k,0) = atom_inp[i].q_x;
            all_q_y[i](k,0) = atom_inp[i].q_y;
            all_p_x[i](k,0) = atom_inp[i].p_x;
            all_p_y[i](k,0) = atom_inp[i].p_y;
        }
        energy(k,0) = H_full->ener(atom_inp,natoms); // energy of initial configuration
    }

    // construct_noise(I, g, natoms);    // create the noise matrices noise_c_x and noise_c_y
    construct_noise_invariant(I, g, natoms);    // create the noise matrices noise_c_x and noise_c_y with zero mean in x,y

    vec * noise_cx = NULL;
    vec * noise_cy = NULL;
    noise_cx = new vec[natoms];
    noise_cy = new vec[natoms];

    vec * noise_fx = NULL;
    vec * noise_fy = NULL;
    noise_fx = new vec[natoms];
    noise_fy = new vec[natoms];

    // initialise of noise vectors
    for (int i=0 ; i < natoms ; i++)
    {
        noise_cx[i].FLdim = I.N_Fine_steps;
        noise_cx[i].FLcoord = new double[I.N_Fine_steps];
        noise_cy[i].FLdim = I.N_Fine_steps;
        noise_cy[i].FLcoord = new double[I.N_Fine_steps];
        noise_fx[i].FLdim = I.N_Fine_steps;
        noise_fx[i].FLcoord = new double[I.N_Fine_steps];
        noise_fy[i].FLdim = I.N_Fine_steps;
        noise_fy[i].FLcoord = new double[I.N_Fine_steps];
    }

    ////////////////////////////////
    // Initial parareal iteration //
    ////////////////////////////////
    particle * atom_new = NULL;
    atom_new = new particle[natoms];

    for(int step_m = 0; step_m < I.N_Coarse_steps ; step_m++)
    {
        // this ensures that we get the correct noise
        for(int step = 0; step < I.N_Fine_steps ; step++)
        {
            for(int i=0 ; i < natoms ; i++)
            {
                noise_cx[i](step) = noise_c_x[i](step_m,step);
                noise_cy[i](step) = noise_c_y[i](step_m,step);
            }
        }

        //// Choose coarse propogator ////
        // atom_new = prop_coarse(I,noise_cx,noise_cy,atom_cur,natoms,well_init); // harmonic app.
        atom_new = prop_coarse_inv(I,noise_cx,noise_cy,atom_cur,natoms,well_init); // harmonic app. with atom 0,1 fixed at origin and x-axis

        energy(0,1+step_m) = H_full->ener(atom_new,natoms);

        for(int i=0 ; i < natoms ; i++)
        {
            atom_cur[i] = atom_new[i];
            all_q_x[i](0,1+step_m) = atom_cur[i].q_x;
            all_q_y[i](0,1+step_m) = atom_cur[i].q_y;
            all_p_x[i](0,1+step_m) = atom_cur[i].p_x;
            all_p_y[i](0,1+step_m) = atom_cur[i].p_y;
        }

        //// Uncomment to calculate well-id
        // well_id(atom_new,natoms,step_m,0);
    }
    cerr << "Finished parareal iteration 0" << endl;

    ///////////////////////////////////
    // Remaining parareal iterations//
    /////////////////////////////////

    particle * atom_prev_c = NULL;
    atom_prev_c = new particle[natoms];
    particle * atom_prev_f = NULL;
    atom_prev_f = new particle[natoms];

    for (int k = 1 ; k < I.Nb_para ; k++)
    {
        // initial data is given by 'atom_inp' for every parareal iteration
        for(int i = 0 ; i < natoms ; i++)
        {
            atom_cur[i] = atom_inp[i];
        }

        // This for-loop calculates trajectory at the 'k-th' parareal iteration
        for(int step_m = 0; step_m < I.N_Coarse_steps ; step_m++)
        {
            // this ensures that we use the correct noise
            for(int step = 0; step < I.N_Fine_steps ; step++)
            {
                for(int i = 0 ; i < natoms ; i++)
                {
                    noise_cx[i](step) = noise_c_x[i](step_m,step);
                    noise_cy[i](step) = noise_c_y[i](step_m,step);
                    noise_fx[i](step) = noise_f_x[i](step_m,step);
                    noise_fy[i](step) = noise_f_y[i](step_m,step);
                }
            }

            // These variables store the state from 'k-1' parareal iteration
            particle atom_prev [natoms];
            for(int i = 0 ; i < natoms ; i++)
            {
                atom_prev[i].q_x = all_q_x[i](k-1,step_m);
                atom_prev[i].q_y = all_q_y[i](k-1,step_m);
                atom_prev[i].p_x = all_p_x[i](k-1,step_m);
                atom_prev[i].p_y = all_p_y[i](k-1,step_m);
            }

            // ** Parareal iteration -- choose coarse/fine propogators ** //

            // atom_prev_c = prop_coarse(I,noise_cx,noise_cy,atom_prev,natoms,well_init); // harmonic app.
            atom_prev_c = prop_coarse_inv(I,noise_cx,noise_cy,atom_prev,natoms,well_init); //harmonic app. with atom 0,1 fixed at origin and x-axis

            // atom_prev_f = prop_fine(I,noise_fx,noise_fy,atom_prev,natoms);
            atom_prev_f = prop_fine_inv(I,noise_fx,noise_fy,atom_prev,natoms); // fine prop. with atom 0,1 fixed at origin and x-axis

            // atom_new = prop_coarse(I,noise_cx,noise_cy,atom_cur,natoms,well_init);// harmonic app.
            atom_new = prop_coarse_inv(I,noise_cx,noise_cy,atom_cur,natoms,well_init); //harmonic app. with atom 0,1 fixed at origin and x-axis

            // **  Parareal update for the k-th iteration ** //
            for(int i = 0 ; i < natoms ; i++)
            {
                atom_cur[i].q_x = atom_new[i].q_x + atom_prev_f[i].q_x - atom_prev_c[i].q_x;
                atom_cur[i].q_y = atom_new[i].q_y + atom_prev_f[i].q_y - atom_prev_c[i].q_y;
                atom_cur[i].p_x = atom_new[i].p_x + atom_prev_f[i].p_x - atom_prev_c[i].p_x;
                atom_cur[i].p_y = atom_new[i].p_y + atom_prev_f[i].p_y - atom_prev_c[i].p_y;
            }
            energy(k,1+step_m) = H_full->ener(atom_cur,natoms);

            // update the matrices storing states
            for(int i = 0 ; i < natoms ; i++)
            {
                all_q_x[i](k,1+step_m) = atom_cur[i].q_x;
                all_q_y[i](k,1+step_m) = atom_cur[i].q_y;
                all_p_x[i](k,1+step_m) = atom_cur[i].p_x;
                all_p_y[i](k,1+step_m) = atom_cur[i].p_y;
            }

            //// Uncomment to calculate well-id for the last parareal iteration
            // if (k == I.Nb_para-1)
            // well_id(atom_new,natoms,step_m,k);
        }
        cerr<<"Finished parareal iteration " << k <<endl;
    }

    //////////////////////////////
    //// Reference trajectory ////
    //////////////////////////////

    //// Used to calculate the first exit event -- parareal should work well before this exit event
    ofstream well_id_ref;
    well_id_ref.open("well_id_ref");
    well_id_ref << setprecision(10);

    // ofstream coord_pot_ener; // potential energy of reference trajectory
    // coord_pot_ener.open("coord_pot_ener");
    // coord_pot_ener << setprecision(10);

    // ofstream coord_force; // force = -nabla V of reference trajectory
    // coord_force.open("coord_force");
    // coord_force << setprecision(10);


    for(int i=0 ; i < natoms ; i++)
    {
        atom_cur[i] = atom_inp[i];
    }

    // to store reference trajectory
    vec ref_q_x [natoms];
    vec ref_q_y [natoms];
    vec ref_p_x [natoms];
    vec ref_p_y [natoms];

    // initialization
    for(int i = 0 ; i < natoms ; i++)
    {
        ref_q_x[i].set_size(1+I.N_Coarse_steps);
        ref_q_y[i].set_size(1+I.N_Coarse_steps);
        ref_p_x[i].set_size(1+I.N_Coarse_steps);
        ref_p_y[i].set_size(1+I.N_Coarse_steps);
    }

    // initial data is given
    for (int i=0 ; i < natoms ; i++)
    {
        ref_q_x[i](0) = atom_inp[i].q_x;
        ref_q_y[i](0) = atom_inp[i].q_y;
        ref_p_x[i](0) = atom_inp[i].p_x;
        ref_p_y[i](0) = atom_inp[i].p_y;
    }


    for(int step_m = 0; step_m < I.N_Coarse_steps ; step_m++)
    {
        // ensures that we use the correct noise
        for(int step = 0; step < I.N_Fine_steps ; step++)
        {
            for(int i=0 ; i < natoms ; i++)
            {
                noise_fx[i](step) = noise_f_x[i](step_m,step);
                noise_fy[i](step) = noise_f_y[i](step_m,step);
            }
        }

        //// force for initial configuration
        // if (step_m == 0)
        // {
        //     double f_init = 0.;
        //     H_full->compute_force(atom_inp,natoms);
        //     for (int i =0 ; i < natoms ; i++)
        //     {
        //         f_init += X2(H_full->force_x[i])+X2(H_full->force_y[i]);
        //     }
        //     coord_force << "0 "<< f_init  << endl;
        // }

        particle * atom_new_ref;
        atom_new_ref = new particle[natoms];

        //// Choose fine propogator ////
        // atom_new_ref = prop_fine(I,noise_fx,noise_fy,atom_cur,natoms);
        atom_new_ref = prop_fine_inv(I,noise_fx,noise_fy,atom_cur,natoms); // atom 0,1 fixed at origin and x-axis

        for(int i=0 ; i < natoms ; i++)
        {
            atom_cur[i] = atom_new_ref[i];

            ref_q_x[i](1+step_m) = atom_cur[i].q_x;
            ref_q_y[i](1+step_m) = atom_cur[i].q_y;
            ref_p_x[i](1+step_m) = atom_cur[i].p_x;
            ref_p_y[i](1+step_m) = atom_cur[i].p_y;
        }

        /*
        Compute objects for reference trjectory
        */
        //// well-id
        well_id(atom_new_ref,natoms, step_m);
        //// potential energy
        // double cur_pot_ener = 0.;
        // cur_pot_ener = H_full->pot_ener(atom_new_ref,natoms);
        // coord_pot_ener << 1+step_m << " " << cur_pot_ener << endl;
        // double f2 = 0.;
        // H_full->compute_force(atom_new_ref,natoms);
        // for (int i = 0 ; i < natoms ; i++)
        // {
        //     f2 += X2(H_full->force_x[i])+X2(H_full->force_y[i]);
        // }
        // coord_force << 1+step_m << " " << f2 << endl;
    }


    ////////////////////////////////
    ////////// * Errors * //////////
    ////////////////////////////////

    // Errors as function of k
    // ofstream para_err; // sum of error for each time-step
    // para_err.open("para_err");
    // para_err << setprecision(10);
    // ofstream para_err_ter; // error at the terminal time-step (I.N_Coarse_steps)
    // para_err_ter.open("para_err_ter");
    // para_err_ter << setprecision(10);
    ofstream para_rel_err_ter; // relative error at the terminal time-step (I.N_Coarse_steps)
    para_rel_err_ter.open("para_rel_err_ter");
    para_rel_err_ter << setprecision(10);
    ofstream para_rel_err_traj; // relative error at the entire trajectory
    para_rel_err_traj.open("para_rel_err_traj");
    para_rel_err_traj << setprecision(10);

    //// Error for each time-step
    // for (int k=0; k<I.Nb_para; k++)
    // {
    //     double err=0.;
    //     for(int step_m = 1; step_m < I.N_Coarse_steps ; step_m++)
    //     {
    //         for(int i=0 ; i < natoms ; i++)
    //         {
    //             err += fabs(all_q_x[i](k,step_m) - ref_q_x[i](step_m));
    //             err += fabs(all_q_y[i](k,step_m) - ref_q_y[i](step_m));
    //         }
    //     }
    //     para_err << k << " " << err << endl;
    // }

    //// Relative error at terminal point
    for (int k=0; k<I.Nb_para; k++)
    {
        double err = 0.;
        double ref = 0.;

        for(int i = 0 ; i < natoms ; i++)
        {
            err += fabs(all_q_x[i](k,I.N_Coarse_steps-1) - ref_q_x[i](I.N_Coarse_steps-1));
            err += fabs(all_q_y[i](k,I.N_Coarse_steps-1) - ref_q_y[i](I.N_Coarse_steps-1));
            ref += fabs(ref_q_x[i](I.N_Coarse_steps-1)) + fabs(ref_q_y[i](I.N_Coarse_steps-1));
        }
        // para_err_ter << k << " " << err << endl;
        para_rel_err_ter << k << " " << err/ref << endl;
    }

    //// Relative error for entire trajectory
    for (int k=0; k<I.Nb_para; k++)
    {
        double err_traj = 0.;
        double ref_traj = 0.;
        for(int step_m = 1; step_m < I.N_Coarse_steps ; step_m++)
        {
            for(int i = 0 ; i < natoms ; i++)
            {
                err_traj += fabs(all_q_x[i](k,step_m) - ref_q_x[i](step_m));
                err_traj += fabs(all_q_y[i](k,step_m) - ref_q_y[i](step_m));
                ref_traj += fabs(ref_q_x[i](step_m)) + fabs(ref_q_y[i](step_m));
            }
        }
        para_rel_err_traj << k << " " << err_traj/ref_traj << endl;
    }

    //// Return the first para-iteration when relative traj-error is greater than threshold
    for (int k=0; k<I.Nb_para; k++)
    {
        double err_traj = 0.;
        double ref_traj = 0.;
        for(int step_m = 1; step_m < I.N_Coarse_steps ; step_m++)
        {
            for(int i = 0 ; i < natoms ; i++)
            {
                err_traj += fabs(all_q_x[i](k,step_m) - ref_q_x[i](step_m));
                err_traj += fabs(all_q_y[i](k,step_m) - ref_q_y[i](step_m));
                ref_traj += fabs(ref_q_x[i](step_m)) + fabs(ref_q_y[i](step_m));
            }
        }
        if (err_traj/ref_traj > 1e+3)
        {
            cout << "Error larger than 1e+3 at para-iteration #" << k << endl;
            break;
        }
    }


    // Relative terminal error as function of N for fixed k's //
    ofstream para0_rel_err_ter; // k=0
    para0_rel_err_ter.open("para0_rel_err_ter");
    para0_rel_err_ter << setprecision(10);
    ofstream para5_rel_err_ter; // k=5
    para5_rel_err_ter.open("para5_rel_err_ter");
    para5_rel_err_ter << setprecision(10);
    ofstream para14_rel_err_ter; // k=14
    para14_rel_err_ter.open("para14_rel_err_ter");
    para14_rel_err_ter << setprecision(10);
    ofstream para35_rel_err_ter; // k=35
    para35_rel_err_ter.open("para35_rel_err_ter");
    para35_rel_err_ter << setprecision(10);

    double err_k0 = 0.;
    double err_k5 = 0.;
    double err_k14 = 0.;
    double err_k35 = 0.;
    double ref = 0.;

    if (I.Nb_para > 35)
    {
        for(int step_m = 0 ; step_m < I.N_Coarse_steps ; step_m++)
        {
            for(int i = 0 ; i < natoms ; i++)
            {
                err_k0 += fabs(all_q_x[i](0,1+step_m) - ref_q_x[i](1+step_m));
                err_k0 += fabs(all_q_y[i](0,1+step_m) - ref_q_y[i](1+step_m));

                err_k5 += fabs(all_q_x[i](5,1+step_m) - ref_q_x[i](1+step_m));
                err_k5 += fabs(all_q_y[i](5,1+step_m) - ref_q_y[i](1+step_m));

                err_k14 += fabs(all_q_x[i](14,1+step_m) - ref_q_x[i](1+step_m));
                err_k14 += fabs(all_q_y[i](14,1+step_m) - ref_q_y[i](1+step_m));

                err_k35 += fabs(all_q_x[i](35,1+step_m) - ref_q_x[i](1+step_m));
                err_k35 += fabs(all_q_y[i](35,1+step_m) - ref_q_y[i](1+step_m));

                ref += fabs(ref_q_x[i](1+step_m)) + fabs(ref_q_y[i](1+step_m));
            }
            para0_rel_err_ter <<  step_m+1 << " " << err_k0/ref << endl;
            para5_rel_err_ter <<  step_m+1 << " " << err_k5/ref << endl;
            para14_rel_err_ter <<  step_m+1 << " " << err_k14/ref << endl;
            para35_rel_err_ter <<  step_m+1 << " " << err_k35/ref << endl;
        }
        err_k0 = 0.;
        err_k5 = 0.;
        err_k14 = 0.;
        err_k35 = 0.;
        ref = 0.;
    }

    ///////////////////////////////////////////////
    ////////// * Efficiency/Gain Plots * //////////
    ///////////////////////////////////////////////

    // Calculate efficiency plots for full trajectory error
    ofstream plot_traj_eff;
    plot_traj_eff.open("plot_traj_eff");
    plot_traj_eff << setprecision(10);

    for (double time_iter = 20., incre = 20. ; time_iter < 1 + I.N_Coarse_steps ; time_iter += incre)
    {
        double err;
        double ref;

        double max_err = 1e-6;
        for (int k = 1 ; k < I.Nb_para + 1 ; k++)
        {
            err = 0.;
            ref = 0.;

            for(int step_m = 1; step_m < time_iter ; step_m++) // these two loops calculates the error at iteration k
            {
                for(int i=0 ; i < natoms ; i++)
                {
                    err += X2(all_q_x[i](k,step_m) - ref_q_x[i](step_m));
                    err += X2(all_q_y[i](k,step_m) - ref_q_y[i](step_m));
                    ref += X2(ref_q_x[i](step_m)) + X2(ref_q_y[i](step_m));
                }
            }
            if (err/ref < max_err)
            {
                plot_traj_eff << time_iter << " " << k/time_iter << endl;
                break;
            }
        }
    }

    // Calculate efficiency plots for terminal time-iteration
    ofstream plot_ter_eff; // Plots k^*/n where k^* is the first para-iteration when error is less than prescribed error
    plot_ter_eff.open("plot_ter_eff");
    plot_ter_eff << setprecision(10);

    for (double time_iter = 20., incre = 20. ; time_iter < 1 + I.N_Coarse_steps ; time_iter += incre)
    {
        double err_ter;
        double ref_ter;

        double max_err = 1e-6;
        for (int k = 1 ; k < I.Nb_para + 1 ; k++)
        {
            err_ter = 0.;
            ref_ter = 0.;

            for(int i=0 ; i < natoms ; i++)
            {
                err_ter += X2(all_q_x[i](k,time_iter) - ref_q_x[i](time_iter));
                err_ter += X2(all_q_y[i](k,time_iter) - ref_q_y[i](time_iter));
                ref_ter += X2(ref_q_x[i](time_iter)) + X2(ref_q_y[i](time_iter));
            }
            if (err_ter/ref_ter < max_err)
            {
                plot_ter_eff << time_iter << " " << k/time_iter << endl;
                break;
            }
        }
    }
}




////** Adaptive slab **////
void Algorithm_langevin::compute_traj_AdSlab(input & I, gaussien & g, particle * atom_inp, int natoms)
{
    particle * well_init;
    well_init = new particle [natoms];
    for (int i = 0 ; i < natoms ; i++)
    {
        well_init[i] = atom_inp[i];
    }// end for (i=0;i<natoms)

    particle * atom_new = NULL;
    atom_new = new particle[natoms];
    particle * atom_prev_c = NULL;
    atom_prev_c = new particle[natoms];
    particle * atom_prev_f = NULL;
    atom_prev_f = new particle[natoms];

    // to store reference trajectory
    vec ref_q_x [natoms];
    vec ref_q_y [natoms];
    vec ref_p_x [natoms];
    vec ref_p_y [natoms];
    // initialization
    for(int i = 0 ; i < natoms ; i++)
    {
        ref_q_x[i].set_size(1+I.N_Coarse_steps);
        ref_q_y[i].set_size(1+I.N_Coarse_steps);
        ref_p_x[i].set_size(1+I.N_Coarse_steps);
        ref_p_y[i].set_size(1+I.N_Coarse_steps);
    }// end for (i=0;i<natoms)

    // construct_noise(I, g, natoms);
    construct_noise_invariant(I, g, natoms);    // create the noise matrices noise_c_x and noise_c_y with zero mean in x,y

    particle atom_cur [natoms];

    int cost = 0; // Wall-clock time neglecting coarse integrations

    //////////////////////////
    // Reference trajectory //
    //////////////////////////
    for(int i=0 ; i < natoms ; i++)
    {
        atom_cur[i] = atom_inp[i];
    }// end for (i=0;i<natoms)

    for (int i=0 ; i < natoms ; i++)
    {
        ref_q_x[i](0) = atom_inp[i].q_x;
        ref_q_y[i](0) = atom_inp[i].q_y;
        ref_p_x[i](0) = atom_inp[i].p_x;
        ref_p_y[i](0) = atom_inp[i].p_y;
    }// end for (i=0;i<natoms)

    vec * noise_fx = NULL;
    vec * noise_fy = NULL;
    noise_fx = new vec[natoms];
    noise_fy = new vec[natoms];
    for (int i=0 ; i < natoms ; i++) // initialise noise vectors
    {
        noise_fx[i].FLdim = I.N_Fine_steps;
        noise_fx[i].FLcoord = new double[I.N_Fine_steps];
        noise_fy[i].FLdim = I.N_Fine_steps;
        noise_fy[i].FLcoord = new double[I.N_Fine_steps];
    }// end for (i=0;i<natoms)

    // Reference trajectory
    for(int step_m = 0; step_m < I.N_Coarse_steps ; step_m++)
    {
        // ensure correct noise
        for(int step = 0; step < I.N_Fine_steps ; step++)
        {
            for(int i=0 ; i < natoms ; i++)
            {
                noise_fx[i](step) = noise_f_x[i](step_m,step);
                noise_fy[i](step) = noise_f_y[i](step_m,step);
            }// end for (i=0;i<natoms)
        }// end for(step=0;step < I.N_Fine_steps)

        //// Choose fine propogator ////
        atom_new = prop_fine(I,noise_fx,noise_fy,atom_cur,natoms);
        // atom_new = prop_fine_inv(I,noise_fx,noise_fy,atom_cur,natoms); // atom 0,1 fixed at origin and x-axis

        for(int i=0 ; i < natoms ; i++)
        {
            atom_cur[i] = atom_new[i];

            ref_q_x[i](1+step_m) = atom_cur[i].q_x;
            ref_q_y[i](1+step_m) = atom_cur[i].q_y;
            ref_p_x[i](1+step_m) = atom_cur[i].p_x;
            ref_p_y[i](1+step_m) = atom_cur[i].p_y;
        }// end for (i=0;i<natoms)
        well_id(atom_new,natoms, step_m); // calculate well-id for reference trajectory
    }// end for(step_m = 0; step_m < I.N_Coarse_steps)

    ////////////////////////////////////////
    // Adaptive Time-slab para trajectory //
    ////////////////////////////////////////
    int max_m = I.N_Coarse_steps; // # total time-iterations
    int init_m = 0; // # counter to keep track of initial
    int N_slab = max_m;

    // Error bounds
    double Rel_Err = 0.001; // to start with relative error is in [Rel_Err_min,Rel_Err_max]
    double delta = 1.;
    double Rel_Err_max = delta; // assume explosion when this threshold is crossed
    double Rel_Err_max2 = delta; // this to dertermine N_err -- time iteration where we stop
    double Rel_Err_min = 1.e-5; // criterion of cnvergence
    // we require Rel_Err_min <= Rel_Err_max2 = Rel_Err_max
    int N_err = 0; // time-step at which things go wrong
    int slab_counter = 1; // keeps track of # of slabs
    int k_counter = 1;

    // store gain as a function of Rel_Err_max = Rel_Err_max2
    ofstream gain_Adap;
    gain_Adap.open("gain_Adap", ofstream::app);
    gain_Adap << setprecision(10);

    ofstream gain_Adap_delta;
    gain_Adap_delta.open("gain_Adap_delta", ofstream::app);
    gain_Adap_delta << setprecision(10);

    ofstream slab_size; // keeps track of slab-size of each slab
    slab_size.open("slab_size");
    slab_size << setprecision(10);


    // ofstream Adap_fine;
    // Adap_fine.open("Adap_fine", ofstream::app);
    // Adap_fine << setprecision(10);


    // Initialise (q,p) for [initial_m,max_m] using coarse integrator
    // here q_traj (k= #para-iteration, step_m = #time-iteration)
    vec para_q_x [natoms];
    vec para_q_y [natoms];
    vec para_p_x [natoms];
    vec para_p_y [natoms];

    vec para_old_q_x [natoms];
    vec para_old_q_y [natoms];
    vec para_old_p_x [natoms];
    vec para_old_p_y [natoms];

    // initialization
    for(int i = 0 ; i < natoms ; i++)
    {
        para_q_x[i].set_size(1 + max_m);
        para_q_y[i].set_size(1 + max_m);
        para_p_x[i].set_size(1 + max_m);
        para_p_y[i].set_size(1 + max_m);

        para_old_q_x[i].set_size(1 + max_m);
        para_old_q_y[i].set_size(1 + max_m);
        para_old_p_x[i].set_size(1 + max_m);
        para_old_p_y[i].set_size(1 + max_m);
    }

    for (int i=0 ; i < natoms ; i++)
    {
        para_q_x[i](0) = atom_inp[i].q_x;
        para_q_y[i](0) = atom_inp[i].q_y;
        para_p_x[i](0) = atom_inp[i].p_x;
        para_p_y[i](0) = atom_inp[i].p_y;
    }

    for(int i=0 ; i < natoms ; i++)
    {
        atom_cur[i] = atom_inp[i];
    }

    vec * noise_cx = NULL;
    vec * noise_cy = NULL;
    noise_cx = new vec[natoms];
    noise_cy = new vec[natoms];
    for (int i=0 ; i < natoms ; i++)
    {
        noise_cx[i].FLdim = I.N_Fine_steps;
        noise_cx[i].FLcoord = new double[I.N_Fine_steps];
        noise_cy[i].FLdim = I.N_Fine_steps;
        noise_cy[i].FLcoord = new double[I.N_Fine_steps];
    }// end for(i=0li<natoms)

    // here init_m = 0, N_slab = max_m
    for(int step_m = init_m; step_m < init_m + N_slab; step_m++)
    {
        // ensure correct noise
        for(int step = 0; step < I.N_Fine_steps ; step++)
        {
            for(int i=0 ; i < natoms ; i++)
            {
                noise_cx[i](step) = noise_c_x[i](step_m,step);
                noise_cy[i](step) = noise_c_y[i](step_m,step);
            }// end for(i=0;i<natoms)
        }// end for(step=0;step<I.N_Fine_steps)

        //// Choose coarse propogator ////
        atom_new = prop_coarse(I,noise_cx,noise_cy,atom_cur,natoms,well_init); // harmonic app.
        // atom_new = prop_coarse_inv(I,noise_cx,noise_cy,atom_cur,natoms,well_init); // harmonic app. with atom 0,1 fixed at origin and x-axis

        for(int i=0 ; i < natoms ; i++)
        {
            atom_cur[i] = atom_new[i];
            para_q_x[i](1+step_m) = atom_cur[i].q_x;
            para_q_y[i](1+step_m) = atom_cur[i].q_y;
            para_p_x[i](1+step_m) = atom_cur[i].p_x;
            para_p_y[i](1+step_m) = atom_cur[i].p_y;
        }// end for(i=0;i<natoms)
    }//end for(step_m)

    cerr << endl;

    while (init_m < max_m)
    {
        cerr<<"Slab " << slab_counter << " : (init_m,N_slab) = ("<<init_m <<","<<N_slab << ")" << endl;

        while ( (Rel_Err < Rel_Err_max && Rel_Err > Rel_Err_min) || isnan(Rel_Err) )
        // while (Rel_Err < Rel_Err_max && Rel_Err > Rel_Err_min)
        {
            cost++; // we add to cost everytime parareal procedure is called

            // store previous parareal iteration
            for (int step_m = init_m; step_m < init_m + N_slab; step_m++)
            {
                for(int i=0 ; i < natoms ; i++)
                {
                    para_old_q_x[i](step_m) = para_q_x[i](step_m);
                    para_old_q_y[i](step_m) = para_q_y[i](step_m);
                    para_old_p_x[i](step_m) = para_p_x[i](step_m);
                    para_old_p_y[i](step_m) = para_p_y[i](step_m);
                }// end for(i=0;i<natoms)
            }//end for(init_m <= step_m <= init_m + N_slab)

            particle atom_prev [natoms];

            for(int i = 0 ; i < natoms ; i++)
            {// initial data is same for every para-iteration
                atom_cur[i].q_x = para_old_q_x[i](init_m);
                atom_cur[i].q_y = para_old_q_y[i](init_m);
                atom_cur[i].p_x = para_old_p_x[i](init_m);
                atom_cur[i].p_y = para_old_p_y[i](init_m);
            }// end for(i=0;i<natoms)

            //** Compute para-traj **//
            // compute next parareal-iteration on time interval [init_m,init_m+N_slab]
            // This should be done in parallel
            for(int step_m = init_m; step_m < init_m + N_slab; step_m++)
            {
                // ensure correct noise
                for(int step = 0; step < I.N_Fine_steps ; step++)
                {
                    for(int i = 0 ; i < natoms ; i++)
                    {
                        noise_cx[i](step) = noise_c_x[i](step_m,step);
                        noise_cy[i](step) = noise_c_y[i](step_m,step);
                        noise_fx[i](step) = noise_f_x[i](step_m,step);
                        noise_fy[i](step) = noise_f_y[i](step_m,step);
                    }// end for(i=0;i<natoms)
                }// end for for(step)

                for(int i = 0 ; i < natoms ; i++)
                {
                    atom_prev[i].q_x = para_old_q_x[i](step_m);
                    atom_prev[i].q_y = para_old_q_y[i](step_m);
                    atom_prev[i].p_x = para_old_p_x[i](step_m);
                    atom_prev[i].p_y = para_old_p_y[i](step_m);
                }

                atom_prev_c = prop_coarse(I,noise_cx,noise_cy,atom_prev,natoms,well_init); // harmonic app.
                // atom_prev_c = prop_coarse_inv(I,noise_cx,noise_cy,atom_prev,natoms,well_init); //harmonic app. with atom 0,1 fixed at origin and x-axis

                atom_prev_f = prop_fine(I,noise_fx,noise_fy,atom_prev,natoms);
                // atom_prev_f = prop_fine_inv(I,noise_fx,noise_fy,atom_prev,natoms); // fine prop. with atom 0,1 fixed at origin and x-axis

                atom_new = prop_coarse(I,noise_cx,noise_cy,atom_cur,natoms,well_init);// harmonic app.
                // atom_new = prop_coarse_inv(I,noise_cx,noise_cy,atom_cur,natoms,well_init); //harmonic app. with atom 0,1 fixed at origin and x-axis

                // **  Parareal update for the k-th iteration ** //
                for(int i = 0 ; i < natoms ; i++)
                {
                    atom_cur[i].q_x = atom_new[i].q_x + atom_prev_f[i].q_x - atom_prev_c[i].q_x;
                    atom_cur[i].q_y = atom_new[i].q_y + atom_prev_f[i].q_y - atom_prev_c[i].q_y;
                    atom_cur[i].p_x = atom_new[i].p_x + atom_prev_f[i].p_x - atom_prev_c[i].p_x;
                    atom_cur[i].p_y = atom_new[i].p_y + atom_prev_f[i].p_y - atom_prev_c[i].p_y;
                }// end for(i=0;i<natoms)

                for(int i = 0 ; i < natoms ; i++)
                {
                    para_q_x[i](1+step_m) = atom_cur[i].q_x;
                    para_q_y[i](1+step_m) = atom_cur[i].q_y;
                    para_p_x[i](1+step_m) = atom_cur[i].p_x;
                    para_p_y[i](1+step_m) = atom_cur[i].p_y;
                }// end for(i=0;i<natoms)

            }//end for(init_m <= step_m <= init_m + N_slab) which calculates para-traj

            //** Compute relative-error comparing previous and current para-iteration **//
            double err_num = 0.;
            double err_den = 0.;
            int m = init_m;
            Rel_Err = 0.;
            // variables for relative error w.r.t. fine trajectory
            double err_num_fine = 0.;
            double err_den_fine = 0.;
            double Rel_Err_fine = 0.;

            // compute relative error as long as Rel_Err < Rel_Err_max2
            while (m < init_m + N_slab && Rel_Err < Rel_Err_max2)
            {
                for(int i = 0 ; i < natoms ; i++)
                {
                    err_num += fabs(para_q_x[i](m) - para_old_q_x[i](m));
                    err_num += fabs(para_q_y[i](m) - para_old_q_y[i](m));
                    err_den += fabs(para_old_q_x[i](m)) + fabs(para_old_q_y[i](m));

                    // errors relative to fine traj
                    err_num_fine += fabs(para_q_x[i](m) - ref_q_x[i](m));
                    err_num_fine += fabs(para_q_y[i](m) - ref_q_y[i](m));
                    err_den_fine += fabs(ref_q_x[i](m))+fabs(ref_q_y[i](m));
                }// end for(i=0;i<natoms)

                Rel_Err = err_num/err_den;
                Rel_Err_fine = err_num_fine/err_den_fine;
                m++;
            }//end while(m)

            if(Rel_Err >= Rel_Err_max2)
            {
                N_err = m-1; // first step when Rel_Err >= Rel_Err_max2
                while (m < init_m + N_slab)  // this while loop completes the computation of relative error
                {
                    for(int i = 0 ; i < natoms ; i++)
                    {
                        err_num += fabs(para_q_x[i](m) - para_old_q_x[i](m));
                        err_num += fabs(para_q_y[i](m) - para_old_q_y[i](m));
                        err_den += fabs(para_old_q_x[i](m)) + fabs(para_old_q_y[i](m));

                        // errors relative to fine traj
                        err_num_fine += fabs(para_q_x[i](m) - ref_q_x[i](m));
                        err_num_fine += fabs(para_q_y[i](m) - ref_q_y[i](m));
                        err_den_fine += fabs(ref_q_x[i](m))+fabs(ref_q_y[i](m));
                    }// end for(i=0;i<natoms)
                    Rel_Err = err_num/err_den;
                    Rel_Err_fine = err_num_fine/err_den_fine;
                    m++;
                }//end while(m)

            }// end if(Rel_Err >= Rel_Err_max2)

            cerr << "Parareal iterations #" << k_counter << " : Rel_Err = "<< Rel_Err << ", Rel_Err_fine = " << Rel_Err_fine << endl;
            // Adap_fine << k_counter << " " << Rel_Err_fine << endl;
            k_counter++;
        }//end while(Rel_Err < Rel_Err_max && Rel_Err > Rel_Err_min)

        if(Rel_Err <= Rel_Err_min)
        {
            // while loop concluded with convergence
            cerr << "Convergence on Slab " << slab_counter << " : (max_m,init_m,N_slab) = (" << max_m << "," << init_m << ","<< N_slab<<")" << endl;
            slab_size << slab_counter << " " << N_slab << endl;
            cerr << endl;
            slab_counter++;
            k_counter = 1;

            init_m = init_m + N_slab;
            if(init_m <= max_m)
            {
                N_slab = max_m-init_m;

                // Initialise (q,p) on [init_m, init_m + N_slab]
                for(int i=0 ; i < natoms ; i++)
                {
                    atom_cur[i].q_x = para_q_x[i](init_m);
                    atom_cur[i].q_y = para_q_y[i](init_m);
                    atom_cur[i].p_x = para_p_x[i](init_m);
                    atom_cur[i].p_y = para_p_y[i](init_m);
                }// end for(i=0;i<natoms)

                for(int step_m = init_m; step_m < init_m + N_slab; step_m++)
                {
                    // ensure correct noise
                    for(int step = 0; step < I.N_Fine_steps ; step++)
                    {
                        for(int i=0 ; i < natoms ; i++)
                        {
                            noise_cx[i](step) = noise_c_x[i](step_m,step);
                            noise_cy[i](step) = noise_c_y[i](step_m,step);
                        }// end for(i=0;i<natoms)
                    }// end for(step)

                    //// Choose coarse propogator ////
                    atom_new = prop_coarse(I,noise_cx,noise_cy,atom_cur,natoms,well_init); // harmonic app.
                    // atom_new = prop_coarse_inv(I,noise_cx,noise_cy,atom_cur,natoms,well_init); // harmonic app. with atom 0,1 fixed at origin and x-axis

                    for(int i=0 ; i < natoms ; i++)
                    {
                        atom_cur[i] = atom_new[i];
                        para_q_x[i](1+step_m) = atom_cur[i].q_x;
                        para_q_y[i](1+step_m) = atom_cur[i].q_y;
                        para_p_x[i](1+step_m) = atom_cur[i].p_x;
                        para_p_y[i](1+step_m) = atom_cur[i].p_y;
                    }// end for(i=0;i<natoms)
                }//end for(step_m)
            }// end if(init_m <= max_m)
        }//end if(Rel_Err <= Rel_Err_min)
        else
        { // while loop concluded with convergence --> modfiy N_slab and continue parareal
            // cerr << "No Convergence on Slab " << slab_counter << " : (max_m,init_m,N_slab,N_err) = (" << max_m << "," << init_m << ","<< N_slab<<","<<N_err<<")" << endl;
            cerr << "No Convergence on slab : N_err = " << N_err<< endl;
            cerr << endl;
            N_slab = N_err-init_m;
            k_counter = 1;
        }//end else if(Rel_Err > Rel_Err_min)

        Rel_Err = 0.001; // reset Rel_Err to between the [Rel_Err_min,Rel_Err_max]

    }//end while(init_m <= max_m)

    // cerr<<"Cost = "<<cost<<", max_m = "<<max_m<<", parareal gain = "<< (double)max_m/cost<<endl;
    cerr << "max_m (# fine-iterations without parareal) = " << max_m << endl;
    cerr << "Cost (# parareal routine called in Adap_para) = " << cost << endl;
    cerr << "Adaptive parareal gain (assuming zero cost for coarse integration) = max_m/Cost = " <<  (double)max_m/(double)cost << endl;

    // store gain vs Rel_Err_max
    gain_Adap << max_m << " " << (double)max_m/(double)cost << endl;
    gain_Adap_delta << delta << " " << (double)max_m/(double)cost << endl;


    // relative error between final adaptive parareal and fine trajectory
    double fine_para_num = 0.;
    double fine_para_den = 0.;
    double fine_para_err = 0.;

    for(int step_m = 0; step_m < I.N_Coarse_steps ; step_m++)
    {
        for(int i = 0 ; i < natoms ; i++)
        {
            fine_para_num += fabs(para_q_x[i](1+step_m) - ref_q_x[i](1+step_m));
            fine_para_num += fabs(para_q_y[i](1+step_m) - ref_q_y[i](1+step_m));
            fine_para_den += fabs(ref_q_x[i](1+step_m))+fabs(ref_q_y[i](1+step_m));
        }// end for(i=0;i<natoms)
    }// end for
    fine_para_err = fine_para_num/fine_para_den;
    cerr << "Relative error between fine and Adap_para trajectory = " << fine_para_err << endl;
    cerr << endl;

    /////////////////////////////
    // Vanilla para trajectory //
    /////////////////////////////

    // store gain as function of theta for vanilla parareal
    ofstream gain_theta;
    gain_theta.open("gain_theta", ofstream::app);
    gain_theta << setprecision(10);
    // store gain as function of gamma for vanilla parareal
    ofstream gain_n;
    gain_n.open("gain_n", ofstream::app);
    gain_n << setprecision(10);

    // store relative error between consecutive vanilla para trajectories
    ofstream van_err;
    van_err.open("van_err");
    van_err << setprecision(10);
    // store relative error between vanilla para and fine traj
    ofstream van_fine;
    van_fine.open("van_fine");
    van_fine << setprecision(10);

    // // track center of mass of k=100
    // ofstream cent_mass_k90;
    // cent_mass_k90.open("cent_mass_k90");
    // cent_mass_k90 << setprecision(10);


    /////  **** For approach
    // (x,y) position and momentum matrices for the i-th particle, with k-th row corresponding to k-th parareal iteration
    mat all_q_x [natoms];
    mat all_q_y [natoms];
    mat all_p_x [natoms];
    mat all_p_y [natoms];
    //Initialisation
    for(int i=0 ; i < natoms ; i++)
    {
        all_q_x[i].set_size(I.Nb_para,1+I.N_Coarse_steps);
        all_q_y[i].set_size(I.Nb_para,1+I.N_Coarse_steps);
        all_p_x[i].set_size(I.Nb_para,1+I.N_Coarse_steps);
        all_p_y[i].set_size(I.Nb_para,1+I.N_Coarse_steps);
    }

    // initial data (t=0) for every parareal iteration is given
    for (int k=0 ; k < I.Nb_para ; k++)
    {
        for (int i=0 ; i < natoms ; i++)
        {
            all_q_x[i](k,0) = atom_inp[i].q_x;
            all_q_y[i](k,0) = atom_inp[i].q_y;
            all_p_x[i](k,0) = atom_inp[i].p_x;
            all_p_y[i](k,0) = atom_inp[i].p_y;
        }
    }

    ////////////////////////////////
    // Initial parareal iteration //
    ////////////////////////////////
    // particle * atom_new = NULL;
    // atom_new = new particle[natoms];
    for(int i=0 ; i < natoms ; i++)
    {
        atom_cur[i] = atom_inp[i];
    }

    for(int step_m = 0; step_m < I.N_Coarse_steps ; step_m++)
    {
        // this ensures that we get the correct noise
        for(int step = 0; step < I.N_Fine_steps ; step++)
        {
            for(int i=0 ; i < natoms ; i++)
            {
                noise_cx[i](step) = noise_c_x[i](step_m,step);
                noise_cy[i](step) = noise_c_y[i](step_m,step);
            }
        }

        //// Choose coarse propogator ////
        atom_new = prop_coarse(I,noise_cx,noise_cy,atom_cur,natoms,well_init); // harmonic app.
        // atom_new = prop_coarse_inv(I,noise_cx,noise_cy,atom_cur,natoms,well_init); // harmonic app. with atom 0,1 fixed at origin and x-axis

        for(int i=0 ; i < natoms ; i++)
        {
            atom_cur[i] = atom_new[i];
            all_q_x[i](0,1+step_m) = atom_cur[i].q_x;
            all_q_y[i](0,1+step_m) = atom_cur[i].q_y;
            all_p_x[i](0,1+step_m) = atom_cur[i].p_x;
            all_p_y[i](0,1+step_m) = atom_cur[i].p_y;
        }
    }

    ///////////////////////////////////
    // Remaining parareal iterations//
    /////////////////////////////////
    int k_break = 0;
    // int nan_check  = 0;

    for (int k = 1 ; k < I.Nb_para ; k++)
    {
        // initial data is given by 'atom_inp' for every parareal iteration
        for(int i = 0 ; i < natoms ; i++)
        {
            atom_cur[i] = atom_inp[i];
        }

        // This for-loop calculates trajectory at the 'k-th' parareal iteration
        for(int step_m = 0; step_m < I.N_Coarse_steps ; step_m++)
        {
            // this ensures that we use the correct noise
            for(int step = 0; step < I.N_Fine_steps ; step++)
            {
                for(int i = 0 ; i < natoms ; i++)
                {
                    noise_cx[i](step) = noise_c_x[i](step_m,step);
                    noise_cy[i](step) = noise_c_y[i](step_m,step);
                    noise_fx[i](step) = noise_f_x[i](step_m,step);
                    noise_fy[i](step) = noise_f_y[i](step_m,step);
                }
            }//end for(step)

            // These variables store the state from 'k-1' parareal iteration
            particle atom_prev [natoms];
            for(int i = 0 ; i < natoms ; i++)
            {
                atom_prev[i].q_x = all_q_x[i](k-1,step_m);
                atom_prev[i].q_y = all_q_y[i](k-1,step_m);
                atom_prev[i].p_x = all_p_x[i](k-1,step_m);
                atom_prev[i].p_y = all_p_y[i](k-1,step_m);
            }

            // ** Parareal iteration -- choose coarse/fine propogators ** //

            atom_prev_c = prop_coarse(I,noise_cx,noise_cy,atom_prev,natoms,well_init); // harmonic app.
            // atom_prev_c = prop_coarse_inv(I,noise_cx,noise_cy,atom_prev,natoms,well_init); //harmonic app. with atom 0,1 fixed at origin and x-axis

            atom_prev_f = prop_fine(I,noise_fx,noise_fy,atom_prev,natoms);
            // atom_prev_f = prop_fine_inv(I,noise_fx,noise_fy,atom_prev,natoms); // fine prop. with atom 0,1 fixed at origin and x-axis

            atom_new = prop_coarse(I,noise_cx,noise_cy,atom_cur,natoms,well_init);// harmonic app.
            // atom_new = prop_coarse_inv(I,noise_cx,noise_cy,atom_cur,natoms,well_init); //harmonic app. with atom 0,1 fixed at origin and x-axis

            // **  Parareal update for the k-th iteration ** //
            for(int i = 0 ; i < natoms ; i++)
            {
                atom_cur[i].q_x = atom_new[i].q_x + atom_prev_f[i].q_x - atom_prev_c[i].q_x;
                atom_cur[i].q_y = atom_new[i].q_y + atom_prev_f[i].q_y - atom_prev_c[i].q_y;
                atom_cur[i].p_x = atom_new[i].p_x + atom_prev_f[i].p_x - atom_prev_c[i].p_x;
                atom_cur[i].p_y = atom_new[i].p_y + atom_prev_f[i].p_y - atom_prev_c[i].p_y;
            }

            // update the matrices storing states
            for(int i = 0 ; i < natoms ; i++)
            {
                all_q_x[i](k,1+step_m) = atom_cur[i].q_x;
                all_q_y[i](k,1+step_m) = atom_cur[i].q_y;
                all_p_x[i](k,1+step_m) = atom_cur[i].p_x;
                all_p_y[i](k,1+step_m) = atom_cur[i].p_y;
            }
        }// end for(step_m = 0; step_m < I.N_Coarse_steps)

        // calculate relative errors
        // calculate relative error between consecutive para trajectories
        double num = 0.;
        double den = 0.;
        // double rel_err_van = 0.;

        double num_fine = 0.;
        double den_fine = 0.;
        // double rel_err_van_fine = 0.;

        for(int step_m = 0; step_m < I.N_Coarse_steps ; step_m++)
        {
            for(int i=0 ; i < natoms ; i++)
            {
                num += fabs(all_q_x[i](k,step_m) - all_q_x[i](k-1,step_m));
                num += fabs(all_q_y[i](k,step_m) - all_q_y[i](k-1,step_m));
                den += fabs(all_q_x[i](k-1,step_m)) + fabs(all_q_y[i](k-1,step_m));

                num_fine += fabs(all_q_x[i](k,step_m) - ref_q_x[i](step_m));
                num_fine += fabs(all_q_y[i](k,step_m) - ref_q_y[i](step_m));
                den_fine += fabs(ref_q_x[i](step_m))+fabs(ref_q_y[i](step_m));
            }
        }
        // if(isnan(num/den))
        // {
        //     nan_check  = 1;
        //     van_fine << k << " " << num_fine/den_fine << endl;
        // }
        // else
        // {
        //     van_fine << k << " " << num_fine/den_fine << endl;
        // }
        van_fine << k << " " << num_fine/den_fine << endl;
        van_err << k << " " << num/den << endl;
        // cerr<<"Finished parareal iteration " << k << " with relative errors w.r.t (k-1,fine) = " << num/den << ","<< num_fine/den_fine <<endl;

        // // store center of mass for iteration k=90
        // if(k == 90)
        // {
        //     double cent = 0.;
        //     for(int step_m = 0; step_m < I.N_Coarse_steps ; step_m++)
        //     {
        //         for(int i=0 ; i < natoms ; i++)
        //         {
        //             cent += all_q_x[i](k,step_m);
        //         }
        //         cent_mass_k90 << step_m << " " << fabs(cent)/7. << endl;
        //         cent = 0.;
        //     }
        // }

        //stopping criterion
        if(isnan(num/den)==false)
        {
            if(num/den <= 1.e-5)
            {
                    k_break = k;// breaks for (k) if relative error is <= 1.e-5
                    break;
            }
        }
    }// end for(k = 1 ; k < I.Nb_para)
    cerr << endl;
    cerr << "max_m (# fine-iterations without parareal) = " << max_m << endl;
    cerr << "Cost (# fine integrations called in Vanilla parareal) = " << k_break << endl;
    cerr << "Parareal gain (assuming zero cost for coarse integration) = max_m/Cost = " <<  (double)max_m/(double)k_break << endl;
    cerr << endl;
    gain_theta << I.beta << " " << (double)max_m/(double)k_break << endl;
    gain_n << max_m << " " << (double)max_m/(double)k_break << endl;

    // n << max_m << " " << (double)max_m/(double)k_break << endl;
    //n_G << I.gamma << " " << (double)max_m/(double)k_break << endl;

    // if(nan_check==1)
    // {
    //     n << max_m << " " << "nan" << endl;
    //     n_G << I.gamma << " " << "nan" << endl;
    // }
    // else
    // {
    //     n << max_m << " " << (double)max_m/(double)k_break << endl;
    //     n_G << I.gamma << " " << (double)max_m/(double)k_break << endl;
    // }

    // relative error between final parareal and fine trajectory
    double van_para_num = 0.;
    double van_para_den = 0.;

    for(int step_m = 0; step_m < I.N_Coarse_steps ; step_m++)
    {
        for(int i = 0 ; i < natoms ; i++)
        {
            van_para_num += fabs(all_q_x[i](k_break,1+step_m) - ref_q_x[i](1+step_m));
            van_para_num += fabs(all_q_y[i](k_break,1+step_m) - ref_q_y[i](1+step_m));
            van_para_den += fabs(ref_q_x[i](1+step_m))+fabs(ref_q_y[i](1+step_m));
        }// end for(i;i<natoms)
    }// end for(step_m)
    cerr << "Relative error between fine and vanilla para trajectory = " << van_para_num/van_para_den << endl;

}


/*
Compute well-id using Conjugate gradient descent
Inputs:
- atom_inp is the configuration of particles for which we want to find the well-id
- natoms is the number of atoms
- Ver_Iter keeps track of the Langevin iteration from which it was called
Output: Returns the well-id (potential energy of the well) of atom_inp
Strategy: first iteration follows simple gradient descent and the remaining iterations follow conjugate gradient descent. Step size is computed using line search
// Observation - Pollack-Ribiere method strongly outperforms Fletcher-Reeves at the saddle points
// Note - For LJ7 cluster in 2d there are 4 wells (Fig 14 BinderLelievreSimpson15)
*/
double Algorithm_langevin::well_id(particle * atom_inp, int natoms, int Ver_Iter)
{
    int max_iter = 5000000;              /* Maximum number of iterations allowed */
    double err_tol = 1e-10;           /* Magnitude of zero gradient tolerance */
    //      double tau = 1e-5;

    particle atom_old [natoms];
    particle atom_cur [natoms];
    particle * atom_new;
    atom_new = new particle[natoms];

    double h_old_x [natoms];
    double h_old_y [natoms];
    double h_cur_x [natoms];
    double h_cur_y [natoms];

    double force_old_2;
    double force_cur_2;


    double gamma_cur = 0.;

    //    double pot_ene_cur = 0.;
    double pot_ene_new = 0.;

    double well_id =0.;
    int iter_num = 0;

    ofstream coord_well;
    coord_well.open("coord_well", ofstream::app);
    coord_well << setprecision(8);
    ofstream coord_well_bis;
    coord_well_bis.open("coord_well_bis", ofstream::app);
    coord_well_bis << setprecision(8);

    ofstream coord_well_ter;
    coord_well_ter.open("coord_well_ter", ofstream::app);
    coord_well_ter << setprecision(10);

    double mu = 0.025; // Line search parameter
    // Original mu = 0.025; decreasing mu makes things slower; increasing to 0.25 makes it faster (but any faster does not change things)
    double residual = 0.;

    ////////// Iteration 0 - Simple gradient descent //////////

    for (int i = 0 ; i < natoms ; i++)
    {
        atom_cur[i].q_x = atom_inp[i].q_x;
        atom_cur[i].q_y = atom_inp[i].q_y;
        atom_cur[i].p_x = 0.;
        atom_cur[i].p_y = 0.;
    }

    H_full->compute_force(atom_cur,natoms);

    double slope = 0.;
    for (int i = 0 ; i < natoms ; i++)
    {
        h_cur_x[i] = H_full->force_x[i];
        h_cur_y[i] = H_full->force_y[i];
        slope += X2(h_cur_x[i]) + X2(h_cur_y[i]);
    }
    slope = - slope;      // To do: Is this slope or sqrt(slope)
    atom_new = LineSearch(atom_cur, h_cur_x, h_cur_y, slope, mu, natoms);

    ////////// Remaining Iterations - Conjugate gradient descent //////////
    for (int n = 0 ; n < max_iter ; n++)
    {
        for (int i = 0 ; i < natoms ; i++)
        {
            atom_old[i].q_x = atom_cur[i].q_x;
            atom_old[i].q_y = atom_cur[i].q_y;
            atom_old[i].p_x = 0.;
            atom_old[i].p_y = 0.;

            h_old_x[i] = h_cur_x[i];
            h_old_y[i] = h_cur_y[i];

            atom_cur[i].q_x = atom_new[i].q_x;
            atom_cur[i].q_y = atom_new[i].q_y;
            atom_old[i].p_x = 0.;
            atom_old[i].p_y = 0.;
        }

        H_full->compute_force(atom_old,natoms);
        force_old_2 = 0.;
        double force_old_x [natoms];
        double force_old_y [natoms];
        for (int i = 0 ; i < natoms ; i++)
        {
            force_old_x[i] = H_full->force_x[i];
            force_old_y[i] = H_full->force_y[i];
            // force_old_2 += X2(H_full->force_x[i]) + X2(H_full->force_y[i]);
            force_old_2 += X2(force_old_x[i]) + X2(force_old_y[i]);
        }

        H_full->compute_force(atom_cur,natoms);
        //        pot_ene_cur = H_full->pot_ener(atom_cur,natoms);
        force_cur_2 = 0.;
        double force_cur_x [natoms];
        double force_cur_y [natoms];
        for (int i = 0 ; i < natoms ; i++)
        {
            force_cur_x[i] = H_full->force_x[i];
            force_cur_y[i] = H_full->force_y[i];
            // force_cur_2 += X2(H_full->force_x[i]) + X2(H_full->force_y[i]);
            force_cur_2 += X2(force_cur_x[i]) + X2(force_cur_y[i]);
        }

        // Fletcher Reeves method
        // gamma_cur = force_cur_2/force_old_2;

        // Pollack-Ribiere method
        double gamma_aux = 0.;
        for (int i=0 ; i < natoms ; i++)
        {
            gamma_aux += force_old_x[i] * force_cur_x[i];
            gamma_aux += force_old_y[i] * force_cur_y[i];
        }
        gamma_cur = (force_cur_2 - gamma_aux)/force_old_2;

        slope = 0.;
        for (int i = 0 ; i < natoms ; i++)
        {
            h_cur_x[i] = H_full->force_x[i] + gamma_cur * h_old_x[i];
            h_cur_y[i] = H_full->force_y[i] + gamma_cur * h_old_y[i];
            slope += h_cur_x[i]*H_full->force_x[i] + h_cur_y[i]*H_full->force_y[i];
        }
        slope = -slope;
        atom_new = LineSearch(atom_cur, h_cur_x, h_cur_y, slope, mu, natoms);

        pot_ene_new = H_full->pot_ener(atom_new,natoms);
        H_full->compute_force(atom_new,natoms);
        residual = 0.;
        for (int i = 0 ; i < natoms ; i++)
        {
            residual += X2(H_full->force_x[i]);
            residual += X2(H_full->force_y[i]);
        }

        iter_num = iter_num + 1;

        // first algo: get out when we do not move anymore
        // if ( fabs(pot_ene_new - pot_ene_cur) < err_tol)
        //   {
        //     break;
        //   }
        // second algo: get out when small gradient
        if ( residual < err_tol ) //1.e-12
        {
            break;
        }

    }

    coord_well << natoms << endl;
    coord_well << "Reference Iteration # = " << 1+Ver_Iter << ". Found well_id in " << iter_num << " iterations; potential energy at min is " << pot_ene_new << " and residual = " << residual << endl;
    coord_well_bis << "Langevin iteration: "<< 1+Ver_Iter << ". Well id: " << pot_ene_new << ", found in " << iter_num << " iterations" << endl;
    coord_well_ter << 1+Ver_Iter << " " << pot_ene_new << endl;
    for (int i=0; i < natoms ; i++)
    {
        coord_well << i+1 << " " << atom_new[i].q_x << " " << atom_new[i].q_y << " " << 0. << endl;
    }

    well_id = pot_ene_new;

    return well_id;
}

/*
Compute well-id (same as above) and returns the configuration at the minima
Output: Returns positions of atoms at minima (instead of well-id above)
*/
particle * Algorithm_langevin::well_id(particle * atom_inp, int natoms)
{
    int max_iter = 5000000;              /* Maximum number of iterations allowed */
    double err_tol = 1e-10;           /* Magnitude of zero gradient tolerance */
    //      double tau = 1e-5;

    particle atom_old [natoms];
    particle atom_cur [natoms];
    particle * atom_new;
    atom_new = new particle[natoms];

    double h_old_x [natoms];
    double h_old_y [natoms];
    double h_cur_x [natoms];
    double h_cur_y [natoms];

    double force_old_2;
    double force_cur_2;


    double gamma_cur = 0.;

    //    double pot_ene_cur = 0.;
    //    double pot_ene_new = 0.;

    // double well_id =0.;
    int iter_num = 0;

    //    ofstream coord_well;
    // coord_well.open("coord_well", ofstream::app);
    //coord_well << setprecision(8);
    //    ofstream coord_well_bis;
    //    coord_well_bis.open("coord_well_bis", ofstream::app);
    //    coord_well_bis << setprecision(8);

    double mu = 0.025; // Line search parameter
    // Original mu = 0.025; decreasing mu makes things slower; increasing to 0.25 makes it faster (but any faster does not change things)
    double residual = 0.;

    ////////// Iteration 0 - Simple gradient descent //////////

    for (int i = 0 ; i < natoms ; i++)
    {
        atom_cur[i].q_x = atom_inp[i].q_x;
        atom_cur[i].q_y = atom_inp[i].q_y;
        atom_cur[i].p_x = 0.;
        atom_cur[i].p_y = 0.;
    }

    H_full->compute_force(atom_cur,natoms);

    double slope = 0.;
    for (int i = 0 ; i < natoms ; i++)
    {
        h_cur_x[i] = H_full->force_x[i];
        h_cur_y[i] = H_full->force_y[i];
        slope += X2(h_cur_x[i]) + X2(h_cur_y[i]);
    }
    slope = - slope;      // To do: Is this slope or sqrt(slope)
    atom_new = LineSearch(atom_cur, h_cur_x, h_cur_y, slope, mu, natoms);

    ////////// Remaining Iterations - Conjugate gradient descent //////////
    for (int n = 0 ; n < max_iter ; n++)
    {
        for (int i = 0 ; i < natoms ; i++)
        {
            atom_old[i].q_x = atom_cur[i].q_x;
            atom_old[i].q_y = atom_cur[i].q_y;
            atom_old[i].p_x = 0.;
            atom_old[i].p_y = 0.;

            h_old_x[i] = h_cur_x[i];
            h_old_y[i] = h_cur_y[i];

            atom_cur[i].q_x = atom_new[i].q_x;
            atom_cur[i].q_y = atom_new[i].q_y;
            atom_old[i].p_x = 0.;
            atom_old[i].p_y = 0.;
        }

        H_full->compute_force(atom_old,natoms);
        force_old_2 = 0.;
        double force_old_x [natoms];
        double force_old_y [natoms];
        for (int i = 0 ; i < natoms ; i++)
        {
            force_old_x[i] = H_full->force_x[i];
            force_old_y[i] = H_full->force_y[i];
            // force_old_2 += X2(H_full->force_x[i]) + X2(H_full->force_y[i]);
            force_old_2 += X2(force_old_x[i]) + X2(force_old_y[i]);
        }

        H_full->compute_force(atom_cur,natoms);
        // pot_ene_cur = H_full->pot_ener(atom_cur,natoms);
        force_cur_2 = 0.;
        double force_cur_x [natoms];
        double force_cur_y [natoms];
        for (int i = 0 ; i < natoms ; i++)
        {
            force_cur_x[i] = H_full->force_x[i];
            force_cur_y[i] = H_full->force_y[i];
            // force_cur_2 += X2(H_full->force_x[i]) + X2(H_full->force_y[i]);
            force_cur_2 += X2(force_cur_x[i]) + X2(force_cur_y[i]);
        }

        // Fletcher Reeves method
        // gamma_cur = force_cur_2/force_old_2;

        // Pollack-Ribiere method
        double gamma_aux = 0.;
        for (int i=0 ; i < natoms ; i++)
        {
            gamma_aux += force_old_x[i] * force_cur_x[i];
            gamma_aux += force_old_y[i] * force_cur_y[i];
        }
        gamma_cur = (force_cur_2 - gamma_aux)/force_old_2;

        slope = 0.;
        for (int i = 0 ; i < natoms ; i++)
        {
            h_cur_x[i] = H_full->force_x[i] + gamma_cur * h_old_x[i];
            h_cur_y[i] = H_full->force_y[i] + gamma_cur * h_old_y[i];
            slope += h_cur_x[i]*H_full->force_x[i] + h_cur_y[i]*H_full->force_y[i];
        }
        slope = -slope;
        atom_new = LineSearch(atom_cur, h_cur_x, h_cur_y, slope, mu, natoms);

        //        pot_ene_new = H_full->pot_ener(atom_new,natoms);
        H_full->compute_force(atom_new,natoms);
        residual = 0.;
        for (int i = 0 ; i < natoms ; i++)
        {
            residual += X2(H_full->force_x[i]);
            residual += X2(H_full->force_y[i]);
        }

        iter_num = iter_num + 1;

        // first algo: get out when we do not move anymore
        // if ( fabs(pot_ene_new - pot_ene_cur) < err_tol)
        //   {
        //     break;
        //   }
        // second algo: get out when small gradient
        if ( residual < err_tol ) //1.e-12
        {
            break;
        }
    }
    return atom_new;
}

/*
Compute well-id for parareal iterations (same as original) and keep track of parareal iteration
Inputs: First three are the same as above; Para_Iter keeps track of the parareal iteration from which it was called
*/
double Algorithm_langevin::well_id(particle * atom_inp, int natoms, int Ver_Iter, int Para_Iter)
{
    int max_iter = 100000;              /* Maximum number of iterations allowed */
    double err_tol = 1e-8;           /* Magnitude of zero gradient tolerance */
    //      double tau = 1e-5;

    particle atom_old [natoms];
    particle atom_cur [natoms];
    particle * atom_new;
    atom_new = new particle[natoms];

    double h_old_x [natoms];
    double h_old_y [natoms];
    double h_cur_x [natoms];
    double h_cur_y [natoms];

    double force_old_2;
    double force_cur_2;


    double gamma_cur = 0.;

    //    double pot_ene_cur = 0.;
    double pot_ene_new = 0.;

    double well_id =0.;
    int iter_num = 0;

    ofstream para_coord_well;
    para_coord_well.open("para_coord_well", ofstream::app);
    para_coord_well << setprecision(8);
    ofstream para_coord_well_bis;
    para_coord_well_bis.open("para_coord_well_bis", ofstream::app);
    para_coord_well_bis << setprecision(8);

    double mu = 0.025; // Line search parameter - different value for paprareal
    double residual = 0.;

    ////////// Iteration 0 - Simple gradient descent //////////

    for (int i = 0 ; i < natoms ; i++)
    {
        atom_cur[i].q_x = atom_inp[i].q_x;
        atom_cur[i].q_y = atom_inp[i].q_y;
        atom_cur[i].p_x = 0.;
        atom_cur[i].p_y = 0.;
    }

    H_full->compute_force(atom_cur,natoms);

    double slope = 0.;
    for (int i = 0 ; i < natoms ; i++)
    {
        h_cur_x[i] = H_full->force_x[i];
        h_cur_y[i] = H_full->force_y[i];
        slope += X2(h_cur_x[i]) + X2(h_cur_y[i]);
    }
    slope = - slope;      // To do: Is this slope or sqrt(slope)
    atom_new = LineSearch(atom_cur, h_cur_x, h_cur_y, slope, mu, natoms);

    ////////// Remaining Iterations - Conjugate gradient descent //////////
    for (int n = 0 ; n < max_iter ; n++)
    {
        for (int i = 0 ; i < natoms ; i++)
        {
            atom_old[i].q_x = atom_cur[i].q_x;
            atom_old[i].q_y = atom_cur[i].q_y;
            atom_old[i].p_x = 0.;
            atom_old[i].p_y = 0.;

            h_old_x[i] = h_cur_x[i];
            h_old_y[i] = h_cur_y[i];

            atom_cur[i].q_x = atom_new[i].q_x;
            atom_cur[i].q_y = atom_new[i].q_y;
            atom_old[i].p_x = 0.;
            atom_old[i].p_y = 0.;
        }

        H_full->compute_force(atom_old,natoms);
        force_old_2 = 0.;
        double force_old_x [natoms];
        double force_old_y [natoms];
        for (int i = 0 ; i < natoms ; i++)
        {
            force_old_x[i] = H_full->force_x[i];
            force_old_y[i] = H_full->force_y[i];
            // force_old_2 += X2(H_full->force_x[i]) + X2(H_full->force_y[i]);
            force_old_2 += X2(force_old_x[i]) + X2(force_old_y[i]);
        }

        H_full->compute_force(atom_cur,natoms);
        // pot_ene_cur = H_full->pot_ener(atom_cur,natoms);
        force_cur_2 = 0.;
        double force_cur_x [natoms];
        double force_cur_y [natoms];
        for (int i = 0 ; i < natoms ; i++)
        {
            force_cur_x[i] = H_full->force_x[i];
            force_cur_y[i] = H_full->force_y[i];
            // force_cur_2 += X2(H_full->force_x[i]) + X2(H_full->force_y[i]);
            force_cur_2 += X2(force_cur_x[i]) + X2(force_cur_y[i]);
        }

        // Fletcher Reeves method
        // gamma_cur = force_cur_2/force_old_2;

        // Pollack-Ribiere method
        double gamma_aux = 0.;
        for (int i=0 ; i < natoms ; i++)
        {
            gamma_aux += force_old_x[i] * force_cur_x[i];
            gamma_aux += force_old_y[i] * force_cur_y[i];
        }
        gamma_cur = (force_cur_2 - gamma_aux)/force_old_2;

        slope = 0.;
        for (int i = 0 ; i < natoms ; i++)
        {
            h_cur_x[i] = H_full->force_x[i] + gamma_cur * h_old_x[i];
            h_cur_y[i] = H_full->force_y[i] + gamma_cur * h_old_y[i];
            slope += h_cur_x[i]*H_full->force_x[i] + h_cur_y[i]*H_full->force_y[i];
        }
        slope = -slope;
        atom_new = LineSearch(atom_cur, h_cur_x, h_cur_y, slope, mu, natoms);

        pot_ene_new = H_full->pot_ener(atom_new,natoms);
        H_full->compute_force(atom_new,natoms);
        residual = 0.;
        for (int i = 0 ; i < natoms ; i++)
        {
            residual += X2(H_full->force_x[i]);
            residual += X2(H_full->force_y[i]);
        }

        iter_num = iter_num + 1;

        // first algo: get out when we do not move anymore
        // if ( fabs(pot_ene_new - pot_ene_cur) < err_tol)
        //   {
        //     break;
        //   }
        // second algo: get out when small gradient
        if ( residual < err_tol ) //1.e-12
        {
            break;
        }

    }

    para_coord_well << natoms << endl;
    para_coord_well << "Parareal Iteration # = " << Para_Iter << ". Langevin Iteration # = " << 1+Ver_Iter << ". Found well_id in " << iter_num << " iterations; potential energy at min is " << pot_ene_new << " and residual = " << residual << endl;
    para_coord_well_bis << "Parareal iteration: "<< Para_Iter<< ". Langevin iteration: "<< 1+Ver_Iter << ". Well id: " << pot_ene_new << ", found in " << iter_num << " iterations" << endl;
    for (int i=0; i < natoms ; i++)
    {
        para_coord_well << i+1 << " " << atom_new[i].q_x << " " << atom_new[i].q_y << " " << 0. << endl;
    }

    well_id = pot_ene_new;

    return well_id;
}


/*
LineSearch:
Inputs:
- current solution corresponds to the current configuration of the particles
- current_direction is the direction of descent (-nabla V in this case)
- slope is the slope in the direction of descent
- mu is a parameter (est un parametre dont j'ai oublié la signification)
*/
particle * Algorithm_langevin::LineSearch(particle * current_solution, double * current_direction_x, double * current_direction_y, double slope, double mu, int natoms)
{
    double FLT_MIN = 1.e-30; // flotant minimal

    int tst = 0;					// useful vars
    double alpha2=0, alpha_tmp=0, alpha_prev=0;
    double alpha_prev2=0;
    double f1=0, f2=0, fprev=0;
    double a=0, b=0;
    double abs_a;
    double c=0, cm11=0, cm12=0, cm21=0, cm22=0;
    double disc=0;
    double new_m=0, old_m=0;

    // We copy current_solution into new_solution, then we update new_solution which is returned at the end
    particle new_solution[natoms];
    for (int i = 0 ; i < natoms ; i++)
    {
        new_solution[i].q_x = current_solution[i].q_x;
        new_solution[i].q_y = current_solution[i].q_y;
        new_solution[i].p_x = 0.;
        new_solution[i].p_y = 0.;
    }

    old_m = H_full->pot_ener(current_solution,natoms);
    // Evaluation at the current solution

    int iterMax = 50;
    int iterNum = 0;                           // iteration counter
    iterNum++;
    double alpha = 1.e-4; // updating step: we start with something small on purpose
    // For Fletcher-Reeves, best alpha = 0.01. Performs very well except for the saddle points
    // For Polak-Rib , best alpha = 0.003 (around 1.e-3). Does very well on saddle points.

    for (int i = 0 ; i < natoms ; i++)
    {
        new_solution[i].q_x += alpha*current_direction_x[i];
        new_solution[i].q_y += alpha*current_direction_y[i];
        new_solution[i].p_x = 0.;
        new_solution[i].p_y = 0.;
    }
    new_m = H_full->pot_ener(new_solution,natoms);
    // Evaluation at the new solution (proposed by line search with step alpha = 1e-5)

    iterNum++;

    // Implementing Goldstein's test for alpha too small
    while (new_m < old_m + (1. - mu)*alpha*slope && iterNum< iterMax)
    {
        //	  cout<<"iteration LS "<<iterNum<<", alpha "<<alpha<<" is too small"<<endl;
        alpha *= 3;
        // on extrapole en alpha
        for (int i = 0 ; i < natoms ; i++)
        {
            new_solution[i].q_x = current_solution[i].q_x + alpha*current_direction_x[i];
            new_solution[i].q_y = current_solution[i].q_y + alpha*current_direction_y[i];
            new_solution[i].p_x = 0.;
            new_solution[i].p_y = 0.;
        }
        new_m = H_full->pot_ener(new_solution,natoms);
        iterNum++;
    }
    if (iterNum == iterMax)
    {
        cout << "Alpha overflowed! Iteration limit on line search = "<<iterMax<<endl;
    }


    // Armijo's test for alpha too large
    alpha_prev = alpha; // H.L. Deng, 6/13/95
    while (new_m > old_m + mu*alpha*slope && iterNum < iterMax)
    {
        //	  cout<<"iteration LS "<<iterNum<<", alpha "<<alpha<<" is too big"<<endl;
        alpha2 = alpha * alpha;
        f1 = new_m - old_m - slope * alpha;

        if (tst == 0)
        {
            alpha_tmp = -slope * alpha2 / (f1 * 2.);
            // tentative alpha
            tst = 1;
        }
        else
        {
            alpha_prev2 = alpha_prev * alpha_prev;
            f2 = fprev - old_m - alpha_prev * slope;

            c = 1. / (alpha - alpha_prev);
            cm11 = 1. / alpha2;
            cm12 = -1. / alpha_prev2;
            cm21 = -alpha_prev / alpha2;
            cm22 = alpha / alpha_prev2;

            a = c * (cm11 * f1 + cm12 * f2);
            b = c * (cm21 * f1 + cm22 * f2);
            disc = b * b - 3. * a * slope;

            abs_a = a;
            if (a < 0)
            {
                abs_a = -a;
            }
            if ((abs_a > FLT_MIN) && (disc > FLT_MIN))
            alpha_tmp = (-b + sqrt(disc)) / (3. * a);
            else
            alpha_tmp = slope * alpha2 / (2. * f1);

            if (alpha_tmp >= .5 * alpha)
            alpha_tmp = .5 * alpha;
        }
        alpha_prev = alpha;
        fprev = new_m;

        if (alpha_tmp < .1 * alpha)
        alpha *= .1;
        else
        alpha = alpha_tmp;
        // we have redefined alpha (the step from which we move forward) using cubic interpolation

        for (int i = 0 ; i < natoms ; i++)
        {
            new_solution[i].q_x = current_solution[i].q_x + alpha*current_direction_x[i];
            new_solution[i].q_y = current_solution[i].q_y + alpha*current_direction_y[i];
            new_solution[i].p_x = 0.;
            new_solution[i].p_y = 0.;
        }
        // calculation of x_n+1 according to this new step
        new_m = H_full->pot_ener(new_solution,natoms);
        // calculation of the new energy
        iterNum++;
    }
    if (iterNum == iterMax)
    {
        cout << "Alpha overflowed! Iteration limit on line search = "<<iterMax<<endl;
    }

    //cout<<"On sort de LS apres "<<iterNum<<" iterations: alpha "<<alpha<<" is acceptable: exit"<<endl;

    particle * new_sol;
    new_sol = new particle[natoms];
    for (int i = 0 ; i < natoms ; i++)
    {
        new_sol[i].q_x = new_solution[i].q_x;
        new_sol[i].q_y = new_solution[i].q_y;
        new_sol[i].p_x = 0.;
        new_sol[i].p_y = 0.;
    }
    return (new_sol);
}

// /*
// Compare hessain for LJ7 potential at initial-well using analytical formula and finite difference
// */
// void Algorithm_langevin::hessian_comp(particle * atoms, int natoms)
// {
//     //////////////////////////////////
//     //// ** Analytical formula ** ////
//     //////////////////////////////////
//     double dx_i = 0.; double dy_i = 0.;
//     double dx_j = 0.; double dy_j = 0.;
//     double dx_l = 0.; double dy_l = 0.;
//     double r_ij_2 = 0.;
//     double r_il_2 = 0.;
//     double theta_ij = 0.; double theta_r_ij = 0.;
//     double theta_il = 0.; double theta_r_il = 0.;

//     particle atom_ana [natoms];
//     for ( int i = 0 ; i < natoms ; i++)
//     {
//         atom_ana[i] = atoms[i];
//     }

//     mat hessian_xx,hessian_yy,hessian_xy;
//     hessian_xx.set_size(natoms,natoms);
//     hessian_yy.set_size(natoms,natoms);
//     hessian_xy.set_size(natoms,natoms);

//     // Initialise matrices
//     for (int i=0 ; i < natoms ; i++)
//     {
//         for (int l = 0 ; l < natoms ; l++)
//         {
//             hessian_xx(i,l) = 0.;
//             hessian_yy(i,l) = 0.;
//             hessian_xy(i,l) = 0.;
//         }
//     }

//     for (int i=0 ; i < natoms ; i++)
//     {
//         dx_i = atom_ana[i].q_x;
//         dy_i = atom_ana[i].q_y;

//         for (int l = 0 ; l < natoms ; l++)
//         {
//             if (l==i)
//             {
//                 for (int j = 0 ; j < natoms ; j++)
//                 {
//                     if (j!=i)
//                     {
//                         dx_j = atom_ana[j].q_x;
//                         dy_j = atom_ana[j].q_y;

//                         r_ij_2 = X2(dx_i - dx_j)+X2(dy_i - dy_j);
//                         theta_ij = - (12./r_ij_2) * (1./(X6(r_ij_2)) - 1./(X3(r_ij_2))); // theta(r_ij)
//                         theta_r_ij = (1./(X2(r_ij_2))) * ( (12.*14.)/(X6(r_ij_2)) - (12.*8.)/(X3(r_ij_2)) );   // theta'(r_ij)/r_ij

//                         hessian_xx(i,i) += theta_r_ij*(dx_i - dx_j)*(dx_i - dx_j) + theta_ij;
//                         hessian_yy(i,i) += theta_r_ij*(dy_i - dy_j)*(dy_i - dy_j) + theta_ij;
//                         hessian_xy(i,i) += theta_r_ij*(dx_i - dx_j)*(dy_i - dy_j);
//                     }
//                 }
//             }
//             else    // l \neq i
//             {
//                 dx_l = atom_ana[l].q_x;
//                 dy_l = atom_ana[l].q_y;

//                 r_il_2 = X2(dx_i - dx_l)+X2(dy_i - dy_l);

//                 theta_il = - (12./r_il_2) * (1./(X6(r_il_2)) - 1./(X3(r_il_2))); // theta(r_il)
//                 theta_r_il = (1./(X2(r_il_2))) * ( (12.*14.)/(X6(r_il_2)) - (12.*8.)/(X3(r_il_2)) );   // theta'(r_il)/r_il

//                 hessian_xx(i,l) += theta_r_il*(dx_l - dx_i)*(dx_i - dx_l) - theta_il;
//                 hessian_yy(i,l) += theta_r_il*(dy_l - dy_i)*(dy_i - dy_l) - theta_il;
//                 hessian_xy(i,l) += theta_r_il*(dx_l - dx_i)*(dy_i - dy_l);
//             }
//         }
//     }

//     //////////////////////////////////
//     //// ** Finite difference ** /////
//     //////////////////////////////////
//     mat FD_xx,FD_yy,FD_xy;
//     FD_xx.set_size(natoms,natoms);
//     FD_yy.set_size(natoms,natoms);
//     FD_xy.set_size(natoms,natoms);
//     double eps = 0.000001; // step size

//     // Initialise matrices
//     for (int i=0 ; i < natoms ; i++)
//     {
//         for (int l = 0 ; l < natoms ; l++)
//         {
//             FD_xx(i,l) = 0.;
//             FD_yy(i,l) = 0.;
//             FD_xy(i,l) = 0.;
//         }
//     }

//     particle atom_FD [natoms];
//     double f_x [natoms];
//     double f_y [natoms];
//     for (int j=0 ; j < natoms ; j++) // Initialise atom_FD
//     {
//         atom_FD[j] = atoms[j];
//         f_x[j] = 0.;
//         f_y[j] = 0.;
//     }
//     H_full->compute_force(atom_FD,natoms);
//     // f_x,f_y contain -force (=\nabla V) evaluated at atoms (given point in R^(7*7))
//     for (int j=0 ; j < natoms ; j++)
//     {
//         f_x[j] = - H_full->force_x[j];
//         f_y[j] = - H_full->force_y[j];
//     }


//     for (int i=0 ; i < natoms ; i++)
//     {
//         // Calculate FD_xx
//         for (int l = 0 ; l < natoms ; l++)
//         {
//             for (int j=0 ; j < natoms ; j++) // Re-initialise atom_FD
//             {
//                 atom_FD[j] = atoms[j];
//             }
//             atom_FD[l].q_x = atom_FD[l].q_x + eps;
//             H_full->compute_force(atom_FD,natoms);
//             FD_xx(i,l) = (-H_full->force_x[i] - f_x[i])/eps;
//         }

//         // Calculate FD_yy
//         for (int l = 0 ; l < natoms ; l++)
//         {
//             for (int j=0 ; j < natoms ; j++) // Re-initialise atom_FD
//             {
//                 atom_FD[j] = atoms[j];
//             }
//             atom_FD[l].q_y = atom_FD[l].q_y + eps;
//             H_full->compute_force(atom_FD,natoms);
//             FD_yy(i,l) = (-H_full->force_y[i] - f_y[i])/eps; // -force = \nabla V
//         }

//         // Calculate FD_xy
//         for (int l = 0 ; l < natoms ; l++)
//         {
//             for (int j=0 ; j < natoms ; j++) // Re-initialise atom_FD
//             {
//                 atom_FD[j] = atoms[j];
//             }
//             atom_FD[l].q_y = atom_FD[l].q_y + eps;
//             H_full->compute_force(atom_FD,natoms);
//             FD_xy(i,l) = (-H_full->force_x[i] - f_x[i])/eps;
//         }
//     }

    // ////////////////////////////////////
    // //// ** Compare the Hessians ** ////
    // ////////////////////////////////////
    // ofstream hess_comp; // stores the hessians from above two approaches
    // hess_comp.open("hess_comp");
    // hess_comp << setprecision(10);

    // for (int i = 0 ; i < natoms ; i++)
    // {
    //     for(int l = 0 ; l < natoms ; l++)
    //     {
    //         hess_comp << "Ana_xx (" << i << "," << l << ") = " << hessian_xx(i,l) << ", FD_xx (" << i << "," << l << ") = " <<FD_xx(i,l) << endl;
    //     }
    //     hess_comp << endl;
    //     for(int l = 0 ; l < natoms ; l++)
    //     {
    //         hess_comp << "Ana_yy (" << i << "," << l << ") = " << hessian_yy(i,l) << ", FD_yy (" << i << "," << l << ") = " <<FD_yy(i,l) << endl;
    //     }
    //     hess_comp << endl;
    //     for(int l = 0 ; l < natoms ; l++)
    //     {
    //         hess_comp << "Ana_xy (" << i << "," << l << ") = " << hessian_xy(i,l) << ", FD_xy (" << i << "," << l << ") = " <<FD_xy(i,l) << endl;
    //     }
    //     hess_comp << endl;
    // }

    // // ////////////////////////////////////////////////////////
    // //// ** Calculate spectrum (using Eigen Library) ** ////
    // ////////////////////////////////////////////////////////
    // ofstream hess_spectrum; // stores the hessians from above two approaches
    // hess_spectrum.open("hess_spectrum");
    // hess_spectrum << setprecision(10);

    // MatrixXd M(14,14); // construct matrix
    // for (int row = 0 ; row < natoms ; row++) // initialise matrix using analytical formula for hessian
    // {
    //     for (int col = 0 ; col < natoms ; col++)
    //     {
    //         M(row,col)= hessian_xx(row,col);
    //         M(row,col+natoms)= hessian_xy(row,col);
    //         M(row+natoms,col)= hessian_xy(row,col);
    //         M(row+natoms,col+natoms)= hessian_yy(row,col);
    //     }
    // }
    // EigenSolver<MatrixXd> es(M); // calculate the spectrum
    // hess_spectrum << "The eigenvalues of the 14x14 matrix are:" << endl << endl << es.eigenvalues() << endl << endl;
    // hess_spectrum << "The matrix of eigenvectors of the 14x14 matrix is:"<< endl << endl << es.eigenvectors() << endl << endl;
    // // print on screen
    // cout << "The eigenvalues of the 14x14 matrix are:" << endl << es.eigenvalues() << endl;
    // // cout << "The matrix of eigenvectors of the 14x14 matrix is:"<< endl << es.eigenvectors() << endl;
// }
