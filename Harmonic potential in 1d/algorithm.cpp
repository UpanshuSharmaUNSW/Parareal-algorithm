#include"algorithm.hpp"
#include "gaussien.hpp"

// The numerical schemes used here are taken from Joubaud-Stoltz-12(MMS, Sec 4.2, page 11)

Algorithm_langevin::Algorithm_langevin(input& I)
{
    noise_c.set_size(I.N_Coarse_steps,I.N_Fine_steps); // Coarse noise
    noise_f.set_size(I.N_Coarse_steps,I.N_Fine_steps); // Fine noise
}

/*
Create noise matrices
*/
void Algorithm_langevin::construct_noise(input& I, gaussien & g)
{
    for(int j=0 ; j < I.N_Coarse_steps; j++)
    {
        for(int k=0 ; k < I.N_Fine_steps ; k++)
        {
            noise_c(j,k) = g.rand_number(0.,1.);
            noise_f(j,k) = noise_c(j,k);
        }
    }
}


/*
Fine Langevin propoagtor
*/
vec Algorithm_langevin::prop_fine(input& I, vec & noise, double & q_, double & p_)
{

    double dt = I.t_step;

    vec new_particle;
    new_particle.set_size(2);

    double alpha = exp(-I.gamma*dt);
    double q_old = q_;
    double p_old = p_;
    double q_new;
    double p_new;

    if (I.type_algo == 0)
    {
        for(int step = 0; step < I.N_Fine_steps ; step++)
        {
            // verlet:
            H_full->compute_force(q_old);

            //p^{n+1/2}
            p_new = p_old + (dt/2.)*H_full->force;
            //q^{n+1}
            q_new = q_old + dt*p_new;
            // on recalcule les forces
            H_full->compute_force(q_new);
            //p^{n+1}
            p_new += (dt/2.)*H_full->force;

            // OU part:
            p_new = p_new*alpha + sqrt((1.-alpha*alpha)/I.beta)*noise(step);

            q_old = q_new;
            p_old = p_new;
        }
    }
    if (I.type_algo != 0)
    {
        cerr<<"this algorithm is not programmed"<<endl;
        exit(-1);
    }
    new_particle(0) = q_new;
    new_particle(1) = p_new;
    return new_particle;
}


/*
Fine Langevin propagator (overloaded - same as above) -- to calculate additional objects
*/
vec Algorithm_langevin::prop_fine(input& I, vec & noise, double & q_, double & p_, double & moy_p2_)
{

    double dt = I.t_step;

    vec new_particle;
    new_particle.set_size(2);

    double alpha = exp(-I.gamma*dt);

    double q_old = q_;
    double p_old = p_;

    double q_new;
    double p_new;

    if (I.type_algo == 0)
    {
        for(int step = 0; step < I.N_Fine_steps ; step++)
        {
            // verlet:
            H_full->compute_force(q_old);

            //p^{n+1/2}
            p_new = p_old + (dt/2.)*H_full->force;
            //q^{n+1}
            q_new = q_old + dt*p_new;
            // on recalcule les forces
            H_full->compute_force(q_new);
            //p^{n+1}
            p_new += (dt/2.)*H_full->force;

            // OU part:
            p_new = p_new*alpha + sqrt((1.-alpha*alpha)/I.beta)*noise(step);

            q_old = q_new;
            p_old = p_new;

            moy_p2_ += p_old*p_old;
        }
    }
    if (I.type_algo != 0)
    {
        cerr<<"this algorithm is not programmed"<<endl;
        exit(-1);
    }
    new_particle(0) = q_new;
    new_particle(1) = p_new;
    return new_particle;
}


/*
Coarse Langevin propoagtor
*/
vec Algorithm_langevin::prop_coarse(input& I, vec & noise, double & q_, double & p_)
{
    double dt = I.t_step;

    vec new_particle;
    new_particle.set_size(2);

    double alpha = exp(-I.gamma*dt);
    double q_old = q_;
    double p_old = p_;
    double q_new;
    double p_new;

    if (I.type_algo == 0)
    {
        for(int step = 0; step < I.N_Fine_steps ; step++)
        {
            // verlet:
            H_cg->compute_force(q_old);

            //p^{n+1/2}
            p_new = p_old + (dt/2.)*H_cg->force;
            //q^{n+1}
            q_new = q_old + dt*p_new;
            // on recalcule les forces
            H_cg->compute_force(q_new);
            //p^{n+1}
            p_new += (dt/2.)*H_cg->force;

            // OU part:
            p_new = p_new*alpha + sqrt((1.-alpha*alpha)/I.beta)*noise(step);

            q_old = q_new;
            p_old = p_new;
        }
    }
    if (I.type_algo != 0)
    {
        cerr<<"this algorithm is not programmed"<<endl;
        exit(-1);
    }
    new_particle(0) = q_new;
    new_particle(1) = p_new;
    return new_particle;
}


void Algorithm_langevin::compute_traj(gaussien & g, input & I, particle & X0)
{
    ofstream plot_eff_iter;
    plot_eff_iter.open("plot_eff_iter");
    plot_eff_iter << setprecision(10);

    double q0 = X0.q;
    double p0 = X0.p;

    double q_cur = q0;
    double p_cur = p0;

    double q_n_k;
    double p_n_k;

    vec new_part(2);
    vec new_part_c(2);
    vec new_part_f(2);
    vec new_part_f_Ham(2);

    // position and momentum matrices for the i-th particle, with k-th row corresponding to k-th parareal iteration
    mat all_q(1+I.Nb_para,1+I.N_Coarse_steps);
    mat all_p(1+I.Nb_para,1+I.N_Coarse_steps);

    mat err(1+I.Nb_para,1+I.N_Coarse_steps);
    mat energy(1+I.Nb_para,1+I.N_Coarse_steps);
    vec max_err(I.Nb_para);

    vec ref_q(1+I.N_Coarse_steps);
    vec ref_p(1+I.N_Coarse_steps);

    double moy_p2;

    for (int k=0;k<1+I.Nb_para;k++)
    {
        all_q(k,0) = q0;
        all_p(k,0) = p0;
        energy(k,0) = H_full->ener(q0,p0);
    }
    construct_noise(I,g);

    vec noise_c_vec(I.N_Fine_steps);
    vec noise_f_vec(I.N_Fine_steps);

    ////////////////////////////////
    // Initial parareal iteration //
    ////////////////////////////////
    for(int step_m = 0; step_m < I.N_Coarse_steps ; step_m++)
    {
        // ensure correct noise
        for(int step = 0; step < I.N_Fine_steps ; step++)
        {
            noise_c_vec(step) = noise_c(step_m,step);
        }

        new_part = prop_coarse(I,noise_c_vec,q_cur,p_cur); // Verlet

        q_cur = new_part(0);
        p_cur = new_part(1);
        all_q(0,1+step_m) = q_cur;
        all_p(0,1+step_m) = p_cur;
        energy(0,1+step_m) = H_full->ener(q_cur,p_cur);
    }

    cerr<<"Finished parareal iteration 0"<<endl;

    //////////////////////////////////
    // Remaining parareal iteration //
    //////////////////////////////////
    for (int k=1;k<1+I.Nb_para;k++)
    {
        q_cur = q0;
        p_cur = p0;

        for(int step_m = 0; step_m < I.N_Coarse_steps ; step_m++)
        {
            // ensure correct noise
            for(int step = 0; step < I.N_Fine_steps ; step++)
            {
                noise_c_vec(step) = noise_c(step_m,step);
                noise_f_vec(step) = noise_f(step_m,step);
            }

            q_n_k = all_q(k-1,step_m);
            p_n_k = all_p(k-1,step_m);

            new_part_c = prop_coarse(I,noise_c_vec,q_n_k,p_n_k);

            new_part_f = prop_fine(I,noise_f_vec,q_n_k,p_n_k);

            new_part = prop_coarse(I,noise_c_vec,q_cur,p_cur);

            q_cur = new_part(0) + new_part_f(0) - new_part_c(0);
            p_cur = new_part(1) + new_part_f(1) - new_part_c(1);

            all_q(k,1+step_m) = q_cur;
            all_p(k,1+step_m) = p_cur;
            energy(k,1+step_m) = H_full->ener(q_cur,p_cur);
        }
        cerr<<"Finished parareal iteration "<<k<<endl;
    }
    // //Remove
    // for (int k=0 ; k< 1+I.Nb_para ; k++)
    // {
    //     for(int step_m = 0; step_m < I.N_Coarse_steps+1 ; step_m++)
    //     cerr<< "para-iter= "<< k << ", time-step =" << step_m << ", (q,p)=" << "("<< all_q(k,step_m) << "," << all_p(k,step_m) << ")" << endl;
    // }

    //////////////////////////
    // Reference trajectory //
    //////////////////////////
    q_cur = q0;
    p_cur = p0;
    ref_q(0) = q_cur;
    ref_p(0) = p_cur;
    moy_p2 = 0.;

    for(int step_m = 0; step_m < I.N_Coarse_steps ; step_m++)
    {
        // ensure correct noise
        for(int step = 0; step < I.N_Fine_steps ; step++)
        {
            noise_f_vec(step) = noise_f(step_m,step);
        }

        new_part = prop_fine(I,noise_f_vec,q_cur,p_cur,moy_p2); // Verlet

        q_cur = new_part(0);
        p_cur = new_part(1);

        ref_q(1+step_m) = q_cur;
        ref_p(1+step_m) = p_cur;
    }
    // //Remove
    // cerr << endl;
    // for (int step_m = 0; step_m < I.N_Coarse_steps+1 ; step_m++)
    // cerr<< "Fine: time-step =" << step_m << ", (q,p)=" << "("<< ref_q(step_m) << "," << ref_p(step_m) << ")" << endl;


    moy_p2 = moy_p2/(I.N_Coarse_steps*I.N_Fine_steps);

    cerr<<"Average p^2 for fine trajectory = "<<moy_p2<<endl;

    ////// Efficiency plot calculation //////
    ofstream plot_eff_iter_rel;
    plot_eff_iter_rel.open("plot_eff_iter_rel");
    plot_eff_iter_rel << setprecision(10);


    for (double time_iter = 25., incre = 25. ; time_iter < 1 + I.N_Coarse_steps ; time_iter += incre)
    {
        double err;
        double max_err = 1e-6;
        for (int k = 1 ; k < I.Nb_para + 1 ; k++)
        {
            err = 0.;
            for(int step_m = 1; step_m < time_iter ; step_m++)
            {
                err += fabs(all_q(k,step_m) - ref_q(step_m));
            }

            if (err < max_err)
            {
                plot_eff_iter << time_iter << " " << k/time_iter << endl;
                break;
            }
        }
    }


    for (double time_iter = 25., incre = 25. ; time_iter < 1 + I.N_Coarse_steps ; time_iter += incre)
    {
        double err;
        double ref;
        double max_err = 1e-6;
        for (int k = 1 ; k < I.Nb_para + 1 ; k++)
        {
            err = 0.;
            ref = 0.;
            for(int step_m = 1; step_m < time_iter ; step_m++)
            {
                err += fabs(all_q(k,step_m) - ref_q(step_m));
                ref += fabs(ref_q(step_m));
            }

            if (err/ref < max_err)
            {
                plot_eff_iter_rel << time_iter << " " << k/time_iter << endl;
                break;
            }
        }
    }

    // Effciency plot for relative terminal error
    ofstream plot_eff_iter_rel_ter;
    plot_eff_iter_rel_ter.open("plot_eff_iter_rel_ter");
    plot_eff_iter_rel_ter << setprecision(10);

    for (double time_iter = 25., incre = 25. ; time_iter < 1 + I.N_Coarse_steps ; time_iter += incre)
    {
        double err_ter;
        double ref_ter;
        double max_err = 1e-6;

        for (int k = 1 ; k < I.Nb_para + 1 ; k++)
        {
            err_ter = 0.;
            ref_ter = 0.;

            err_ter = fabs(all_q(k,time_iter) - ref_q(time_iter));
            ref_ter = fabs(ref_q(time_iter));

            if (err_ter/ref_ter < max_err)
            {
                plot_eff_iter_rel_ter << time_iter << " " << k/time_iter << endl;
                break;
            }

        }
    }

    //////////////////////
    /////// Errors ///////
    //////////////////////

    // Relative error for full trajectory //
    ofstream para_rel_err_traj;
    para_rel_err_traj.open("para_rel_err_traj");
    para_rel_err_traj << setprecision(10);

    for (int k=0; k<I.Nb_para; k++)
    {
        double err=0.;
        double ref=0.;

        for(int step_m = 1; step_m < I.N_Coarse_steps ; step_m++)
        {
            err += fabs(all_q(k,step_m) - ref_q(step_m));
            ref += fabs(ref_q(step_m));
        }
        para_rel_err_traj << k << " " << err/ref << endl;
    }

    // Relative Terminal error = Terminal error/reference position at terminal point as a function of k //
    ofstream para_rel_err_ter;
    para_rel_err_ter.open("para_rel_err_ter");
    para_rel_err_ter << setprecision(10);

    for (int k=0; k<I.Nb_para; k++)
    {
        para_rel_err_ter << k << " " <<  fabs(ref_q(I.N_Coarse_steps-1)-all_q(k,I.N_Coarse_steps-1))/fabs(ref_q(I.N_Coarse_steps-1)) << endl;
    }

    ///// Sanity checks ////
    // fine trajectory is stable
    ofstream fine_pos;
    fine_pos.open("fine_pos");
    fine_pos << setprecision(10);

    ofstream fine_mom;
    fine_mom.open("fine_mom");
    fine_mom << setprecision(10);

    ofstream para_explode_k_pos;
    para_explode_k_pos.open("para_explode_k_pos");
    para_explode_k_pos << setprecision(10);

    for(int step_m = 1; step_m < I.N_Coarse_steps ; step_m++)
    {
        fine_pos << step_m << " " << ref_q(1+step_m) << endl;
        fine_mom << step_m << " " << ref_p(1+step_m) << endl;
    }
}



// Adaptive slab
void Algorithm_langevin::compute_traj_AdSlab(gaussien & g, input & I, particle & X0)
{
    double q0 = X0.q;
    double p0 = X0.p;

    double q_init = q0;
    double p_init = p0;

    vec new_part(2);
    vec new_part_c(2);
    vec new_part_f(2);
    vec new_part_f_Ham(2);

    // vector to store fine trajectory
    vec ref_q(1+I.N_Coarse_steps);
    vec ref_p(1+I.N_Coarse_steps);

    double moy_p2;

    construct_noise(I,g);

    double q_cur;
    double p_cur;

    int cost = 0; // Wall-clock time neglecting coarse integrations

    //////////////////////////
    // Reference trajectory //
    //////////////////////////
    q_cur = q0;
    p_cur = p0;
    ref_q(0) = q0;
    ref_p(0) = p0;
    moy_p2 = 0.;
    vec noise_f_vec(I.N_Fine_steps);

    // Reference trajectory
    for(int step_m = 0; step_m < I.N_Coarse_steps ; step_m++)
    {
        // ensure correct noise
        for(int step = 0; step < I.N_Fine_steps ; step++)
        {
            noise_f_vec(step) = noise_f(step_m,step);
        }

        new_part = prop_fine(I,noise_f_vec,q_cur,p_cur,moy_p2); // Verlet

        q_cur = new_part(0);
        p_cur = new_part(1);

        ref_q(1+step_m) = q_cur;
        ref_p(1+step_m) = p_cur;
    }

    // ofstream fine_q;
    // fine_q.open("fine_q");
    // fine_q << setprecision(10);
    // ofstream fine_p;
    // fine_p.open("fine_p");
    // fine_p << setprecision(10);

    // for (int step_m = 0 ; step_m < I.N_Coarse_steps ; step_m++)
    // {
    //     fine_q << step_m << " " << ref_q(step_m) << endl;
    //     fine_p << step_m << " " << ref_p(step_m) << endl;
    // }

    ////////////////////////////////////////
    // Adaptive Time-slab para trajectory //
    ////////////////////////////////////////
    int max_m = I.N_Coarse_steps; // # total time-iterations
    int init_m = 0; // # counter to keep track of initial
    int N_slab = max_m;

    // Error bounds
    double Rel_Err = 0.001; // to start with relative error is in [Rel_Err_min,Rel_Err_max]
    double delta = 9.;
    double Rel_Err_max = delta; // assume explosion when this threshold is crossed
    double Rel_Err_max2 = delta; // this to dertermine N_err -- time iteration where we stop
    double Rel_Err_min = 1.e-5; // criterion of cnvergence
    // we require Rel_Err_min <= Rel_Err_max2 < Rel_Err_max
    int N_err = 0; // time-step at which things go wrong
    int slab_counter = 1; // keeps track of # of slabs
    int k_counter = 1;

    // store gain as a function of N
    ofstream adap_gain;
    adap_gain.open("adap_gain", ofstream::app);
    adap_gain << setprecision(10);

    // ofstream slab_size; // keeps track of slab-size of each slab
    // slab_size.open("slab_size");
    // slab_size << setprecision(10);

    // Initialise (q,p) for [initial_m,max_m] using coarse integrator
    // here q_traj (k= #para-iteration, step_m = #time-iteration)
    vec para_q(1 + max_m);
    vec para_p(1 + max_m);
    vec para_q_old(1+max_m);
    vec para_p_old(1+max_m);
    para_q(0) = q_init;
    para_p(0) = p_init;
    q_cur = para_q(0);
    p_cur = para_p(0);
    vec noise_c_vec(I.N_Fine_steps);

    // here init_m = 0, N_slab = max_m
    for(int step_m = init_m; step_m < init_m + N_slab; step_m++)
    {
        // ensure correct noise
        for(int step = 0; step < I.N_Fine_steps ; step++)
        {
            noise_c_vec(step) = noise_c(step_m,step);
        }// end for(step)

        new_part = prop_coarse(I,noise_c_vec,q_cur,p_cur); // coarse-integration step
        q_cur = new_part(0);
        p_cur = new_part(1);

        // store (q,p) trajectory
        para_q(1+step_m) = q_cur;
        para_p(1+step_m) = p_cur;
    }//end for(step_m)

    // cerr << endl;

    while (init_m < max_m)
    {
        // cerr<<"Slab " << slab_counter << " : (init_m,N_slab) = ("<<init_m <<","<<N_slab << ")" << endl;

        while (Rel_Err < Rel_Err_max && Rel_Err > Rel_Err_min)
        {
            cost++; // we add to cost everytime parareal procedure is called

            //** Compute para-traj **//
            // store previous parareal iteration
            for (int step_m = init_m; step_m < init_m + N_slab; step_m++)
            {
                para_q_old(step_m) = para_q(step_m);
                para_p_old(step_m) = para_p(step_m);
            }//end for(init_m <= step_m <= init_m + N_slab)

            double q_old = 0.;
            double p_old = 0.;
            q_cur = para_q_old(init_m);
            p_cur = para_p_old(init_m);

            // compute next parareal-iteration on time interval [init_m,init_m+N_slab]
            // This should be done in parallel
            for(int step_m = init_m; step_m < init_m + N_slab; step_m++)
            {
                // ensure correct noise
                for(int step = 0; step < I.N_Fine_steps ; step++)
                {
                    noise_c_vec(step) = noise_c(step_m,step);
                    noise_f_vec(step) = noise_f(step_m,step);
                }// end for for(step)

                q_old = para_q_old(step_m);
                p_old = para_p_old(step_m);

                new_part_c = prop_coarse(I,noise_c_vec,q_old,p_old);

                new_part_f = prop_fine(I,noise_f_vec,q_old,p_old);

                new_part = prop_coarse(I,noise_c_vec,q_cur,p_cur);

                q_cur = new_part(0) + new_part_f(0) - new_part_c(0);
                p_cur = new_part(1) + new_part_f(1) - new_part_c(1);

                para_q(1+step_m) = q_cur;
                para_p(1+step_m) = p_cur;
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
                err_num += fabs(para_q(m) - para_q_old(m));
                err_den += fabs(para_q_old(m));
                Rel_Err = err_num/err_den;
                // calculate relative error w.r.t. fine trajectory
                err_num_fine += fabs(para_q(m) - ref_q(m));
                err_den_fine += fabs(ref_q(m));
                Rel_Err_fine = err_num_fine/err_den_fine;
                ///////////
                m++;
            }//end while(m)

            if(Rel_Err >= Rel_Err_max2)
            {
                N_err = m-1; // first step when Rel_Err >= Rel_Err_max2
                while (m < init_m + N_slab)  // this while loop completes the computation of relative error
                {
                    err_num += fabs(para_q(m) - para_q_old(m));
                    err_den += fabs(para_q_old(m));
                    Rel_Err = err_num/err_den;
                    // calculate relative error w.r.t. fine trajectory
                    err_num_fine += fabs(para_q(m) - ref_q(m));
                    err_den_fine += fabs(ref_q(m));
                    Rel_Err_fine = err_num_fine/err_den_fine;
                    //////////
                    m++;
                }//end while(m)

            }// end if(Rel_Err >= Rel_Err_max2)
            // cerr << "Parareal iterations : Slab " << slab_counter << ": Rel_Err = "<<Rel_Err << endl;
            // cerr << "Parareal iterations #" << k_counter << " : Rel_Err = "<< Rel_Err << ", Rel_Err_fine = " << Rel_Err_fine << endl;
            k_counter++;
        }//end while(Rel_Err < Rel_Err_max && Rel_Err > Rel_Err_min)

        if(Rel_Err <= Rel_Err_min)
        {
            // while loop concluded with convergence
            cerr << "Convergence on Slab " << slab_counter << " : (max_m,init_m,N_slab) = (" << max_m << "," << init_m << ","<< N_slab<<")" << endl;
            // slab_size << slab_counter << " " << N_slab << endl;
            cerr << endl;
            slab_counter++;
            k_counter = 1;

            init_m = init_m + N_slab;
            if(init_m <= max_m)
            {
                N_slab=max_m-init_m;

                // Initialise (q,p) on [init_m, init_m + N_slab]
                q_cur = para_q(init_m);
                p_cur = para_p(init_m);
                for(int step_m = init_m; step_m < init_m + N_slab; step_m++)
                {
                    // ensure correct noise
                    for(int step = 0; step < I.N_Fine_steps ; step++)
                    {
                        noise_c_vec(step) = noise_c(step_m,step);
                    }// end for(step)

                    new_part = prop_coarse(I,noise_c_vec,q_cur,p_cur); // coarse-integration step
                    q_cur = new_part(0);
                    p_cur = new_part(1);

                    // store (q,p) trajectory
                    para_q(1+step_m) = q_cur;
                    para_p(1+step_m) = p_cur;
                }//end for(step_m)
            }// end if(init_m <= max_m)
        }//end if(Rel_Err <= Rel_Err_min)
        else
        { // while loop concluded with convergence --> modfiy N_slab and continue parareal
            cerr << "No Convergence on Slab " << slab_counter << " : (max_m,init_m,N_slab,N_err) = (" << max_m << "," << init_m << ","<< N_slab<<","<<N_err<<")" << endl;
            cerr << "No Convergence on slab : N_err = " << N_err<< endl;
            cerr << endl;
            N_slab = N_err-init_m;
            k_counter = 1;
        }//end else if(Rel_Err > Rel_Err_min)

        Rel_Err = 0.001; // reset Rel_Err to between the [Rel_Err_min,Rel_Err_max]

    }//end while(init_m <= max_m)

    // cerr<<"Cost = "<<cost<<", max_m = "<<max_m<<", parareal gain = "<< (double)max_m/cost<<endl;
    cerr << "max_m (# fine-iterations without parareal) = " << max_m << endl;
    cerr << endl;
    cerr << "Cost (# parareal routine called in Adap_para) = " << cost << endl;

    cerr << "Parareal gain (assuming zero cost for coarse integration) = max_m/Cost = " <<  (double)max_m/(double)cost << endl;

    adap_gain << max_m << " " << (double)max_m/(double)cost << endl;

    // relative error between final parareal and fine trajectory
    double fine_para_num = 0.;
    double fine_para_den = 0.;
    double fine_para_err = 0.;

    for(int step_m = 0; step_m < I.N_Coarse_steps ; step_m++)
    {
        fine_para_num += fabs(para_q(1+step_m) - ref_q(1+step_m));
        fine_para_den += fabs(ref_q(1+step_m));
    }// end for
    fine_para_err = fine_para_num/fine_para_den;
    cerr << "Relative error between fine and Adap_para trajectory = " << fine_para_err << endl;
    cerr << endl;


    /////////////////////////////
    // Vanilla para trajectory //
    /////////////////////////////
    // store relative error
    ofstream van_err;
    van_err.open("van_err");
    van_err << setprecision(10);

    // store relative error w.r.t. fine
    ofstream van_fine;
    van_fine.open("van_fine");
    van_fine << setprecision(10);

    // store gain as a function of N,gamma
    ofstream van_gain;
    van_gain.open("van_gain", ofstream::app);
    van_gain << setprecision(10);

    // store parareal trajectory to compare with fine dynamics
    ofstream para_q_Int;
    para_q_Int.open("para_q_Int");
    para_q_Int << setprecision(10);
    ofstream para_p_Int;
    para_p_Int.open("para_p_Int");
    para_p_Int << setprecision(10);

    ofstream para_q_Int_80;
    para_q_Int_80.open("para_q_Int_80");
    para_q_Int_80 << setprecision(10);
    ofstream para_p_Int_80;
    para_p_Int_80.open("para_p_Int_80");
    para_p_Int_80 << setprecision(10);

    ofstream para_q_stab;
    para_q_stab.open("para_q_stab");
    para_q_stab << setprecision(10);
    ofstream para_p_stab;
    para_p_stab.open("para_p_stab");
    para_p_stab << setprecision(10);



    // variables
    Rel_Err = 0.001;
    int van_cost = 0;
    int k_van = 1;

    para_q(0) = q_init;
    para_p(0) = p_init;
    q_cur = para_q(0);
    p_cur = para_p(0);


    // initial parareal iterations
    for(int step_m = 0; step_m < max_m; step_m++)
    {
        // ensure correct noise
        for(int step = 0; step < I.N_Fine_steps ; step++)
        {
            noise_c_vec(step) = noise_c(step_m,step);
        }// end for(step)

        new_part = prop_coarse(I,noise_c_vec,q_cur,p_cur); // coarse-integration step
        q_cur = new_part(0);
        p_cur = new_part(1);

        // store (q,p) trajectory
        para_q(1+step_m) = q_cur;
        para_p(1+step_m) = p_cur;
    }//end for(step_m)

    // run para as long as error does not converge
    while (Rel_Err > Rel_Err_min) // Rel_Err_min = 1.e-5
    {
        van_cost++; // add cost each time para-routine is called
        //** Compute para-traj **//

        // store previous parareal iteration
        for (int step_m = 0; step_m < max_m; step_m++)
        {
            para_q_old(step_m) = para_q(step_m);
            para_p_old(step_m) = para_p(step_m);
        }//end for(step_m)

        // initialise k-th para iteration
        double q_old = 0.;
        double p_old = 0.;
        q_cur = para_q_old(0);
        p_cur = para_p_old(0);

        for(int step_m = 0; step_m < max_m; step_m++)
        {
            // ensure correct noise
            for(int step = 0; step < I.N_Fine_steps ; step++)
            {
                noise_c_vec(step) = noise_c(step_m,step);
                noise_f_vec(step) = noise_f(step_m,step);
            }// end for for(step)

            q_old = para_q_old(step_m);
            p_old = para_p_old(step_m);

            new_part_c = prop_coarse(I,noise_c_vec,q_old,p_old);

            new_part_f = prop_fine(I,noise_f_vec,q_old,p_old);

            new_part = prop_coarse(I,noise_c_vec,q_cur,p_cur);

            q_cur = new_part(0) + new_part_f(0) - new_part_c(0);
            p_cur = new_part(1) + new_part_f(1) - new_part_c(1);

            para_q(1+step_m) = q_cur;
            para_p(1+step_m) = p_cur;
            if(van_cost == 40)
            {
                para_q_Int << step_m << " " << q_cur << endl;
                para_p_Int << step_m << " " << p_cur << endl;
            }
            if(van_cost == 80)
            {
                para_q_Int_80 << step_m << " " << q_cur << endl;
                para_p_Int_80 << step_m << " " << p_cur << endl;
            }
            if (van_cost == 150)
            {
                para_q_stab << step_m << " " << q_cur << endl;
                para_p_stab << step_m << " " << p_cur << endl;
            }
        }//end for(init_m <= step_m <= init_m + N_slab) which calculates para-traj


        //** Compute relative-error comparing previous and current para-iteration **//
        double van_num = 0.;
        double van_den = 0.;
        int m = 0;
        Rel_Err = 0.;
        // variables for relative error w.r.t. fine trajectory
        double van_num_fine = 0.;
        double van_den_fine = 0.;
        double Rel_Err_fine = 0.;

        // compute relative error
        while (m < max_m)
        {
            van_num += fabs(para_q(m) - para_q_old(m));
            van_den += fabs(para_q_old(m));
            Rel_Err = van_num/van_den;

            // calculate relative error w.r.t. fine trajectory
            van_num_fine += fabs(para_q(m) - ref_q(m));
            van_den_fine += fabs(ref_q(m));
            Rel_Err_fine = van_num_fine/van_den_fine;

            m++;
        }//end while(m)

        van_fine << van_cost << " " << Rel_Err_fine << endl;
        van_err << van_cost << " " << Rel_Err << endl;
        // cerr << "Parareal iterations #" << k_van << "  : Rel_Err =  " << Rel_Err << ", Rel_Err_Fine = " << Rel_Err_fine << endl;
        k_van++;
    }//end while (Rel_Err > Rel_Err_min)

    cerr << "Cost (# fine integrations called in Vanilla parareal) = " << van_cost << endl;
    cerr << "Parareal gain (assuming zero cost for coarse integration) = max_m/Cost = " <<  (double)max_m/(double)van_cost << endl;
    cerr << endl;

    van_gain << max_m << " " << (double)max_m/(double)van_cost << endl;

    // relative error between final parareal and fine trajectory
    double van_para_num = 0.;
    double van_para_den = 0.;
    double van_para_err = 0.;

    for(int step_m = 0; step_m < I.N_Coarse_steps ; step_m++)
    {
        van_para_num += fabs(para_q(1+step_m) - ref_q(1+step_m));
        van_para_den += fabs(ref_q(1+step_m));
    }// end for
    van_para_err = van_para_num/van_para_den;
    cerr << "Relative error between fine and Adap_para trajectory = " << van_para_err << endl;

    // cerr << "Parareal gain_van = " <<  (double)max_m/(double)van_cost << endl;
    // cerr << "Parareal gain_adap = " <<  (double)max_m/(double)cost << endl;

}
