# Parareal-algorithm

## Project description 
This C++ repository deals with implements parareal algorithm for the (underdamped) Langevin equation where the coarse an fine integrators correspond to different force fields. This code was used to produce the results in [An Adaptive Parareal Algorithm: Application to the Simulation of Molecular Dynamics Trajectories](https://epubs.siam.org/doi/abs/10.1137/21M1412979) (see [hal-03189428](https://hal.science/hal-03189428) for a preprint). The repository contains two folders which discuss different examples: (1) harmonic potential in one dimension and (2) Lennard-Jones seven-atom cluster in two dimensions. 

## Instuctions
(1) Build a make file. For instance this can be achived by downloading the code and using the **makefile** command in terminal.

(2) Enter the parameters for the parareal algorithm in the **input_file** (see the following section for details).

(3) Run the executable file by entering **./CG++** in the terminal. 

(4) The outputs are stored in the following text files:
   
(a) File **van_err** stores the parareal iteration and corresponding trajectorial error between consecutive parareal trajectories. This error is called $E(k,N,\gamma,\beta)$ in the article (equation (3.1)). This file is overwritten during every run of the code. 

(b) File **van_fine** stores the parareal iteration and corresponding trajectorial error between parareal trajectory and the fine(/true) trajectory. This error is called $E_f(k,N,\gamma,\beta)$ in the article (equation (3.4)). Every run of the executable file overwrites this text file. 

(c) File **gain_n** stores the the total number of time steps and corresponidng gain for the classical parareal algorithm. The gain is defined in equation (3.3) in the article. Every run of the executable file appends to this text file.   

(d) File **gain_adap** stores the the total number of time steps and corresponidng gain for the adaptive parareal algorithm. The gain is defined in equation (3.3) in the article. Every run of the executable file appends to this text file.  


## Parameters for parareal algorithm

Each of the folder contains an identical file titled **input_file** which contains the various parameters used in the parareal algorithm, which we now describe. 

**WARNING: CHANGING THE FORMAT OF THE MAKE FILE WILL LEAD TO ERRONEOUS RESULTS. ONLY CHANGE THE NUMBERS IN THE input_file**   

(1) Number of para iterations: Number of parareal iterations to be used in the current run. This correspond to variable $K$ in the article. 

(2) Number of coarse intervals: Number of coarse-intervals to be used in the current run. Wile this number can be different from (1), in our experiments this is always chosen to be equal to input in (1). This correspond to variable $N$ in the article.  

(3) Number of fine steps: Number of fine steps within each coarse step. In our experiments this is always chosen to be 1, i.e. the fine and coarse time-step are exactly the same. For explanation see the article. 

(4) Time step: This is the time step (same for fine and coarse integration) that is used in the algorithm. This corresonds to $\Delta t$ in the article.

(5) Algorithm type: This number should always be **0**.  

(6) Beta: Value of inverse temperature in the Langevin dynamics (see $\beta$ in article).

(7) Gamma: Value of the friction coefficient in the Langevin dynamics (see $\gamma$ in the article).  


## Parameter for adaptive parareal algorith
The adaptive parareal algorithm uses an additional explosion parameter $\delta_{\mathrm{expl}}$ (see Algorithm 1 in the article). The value for this parameter is assigned to the variable **delta** in the **algorithm_natoms.cpp** file. 
