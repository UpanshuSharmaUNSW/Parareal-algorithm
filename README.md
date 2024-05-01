# Parareal-algorithm

## Project description 
This repository deals with implements parareal algorithm for the (underdamped) Langevin equation where the coarse an fine integrators correspond to different force fields. This code was used to produce the results in [An Adaptive Parareal Algorithm: Application to the Simulation of Molecular Dynamics Trajectories](https://epubs.siam.org/doi/abs/10.1137/21M1412979) (see [hal-03189428](https://hal.science/hal-03189428) for a preprint). The repository contains two folders which discuss different examples: (1) harmonic potential in one dimension and (2) Lennard-Jones seven-atom cluster in two dimensions. 

## Basic instuctions

### Parameters for parareal algorithm
Each of the folder contains an identical file titled **input_file** which contains the various parameters used in the parareal algorithm, which we now describe. 

**WARNING: CHANGING THE FORMAT OF THE MAKE FILE WILL LEAD TO ERRONEOUS RESULTS. ONLY CHANGE THE NUMBERS IN THE input_file**   

(1) Number of para iterations: Number of parareal iterations to be used in the current run. This correspond to variable $\mathbf{K}$ in the article above. 

(2) Number of coarse intervals: Number of coarse-intervals to be used in the current run. Wile this number can be different from (1), in our experiments this is always chosen to be equal to input in (1). For explanation see the article. 

(3) Number of fine steps: Number of fine steps within each coarse step. In our experiments this is always chosen to be 1, i.e. the fine and coarse time-step are exactly the same. For explanation see the article. 

(4) Time step: This is the time step (same for fine and coarse integration) that is used in the algorithm. This corresonds to $\Delta t$ in the article.

(5) Algorithm type: This number should always be **0**.  

(6) Beta: Value of inverse temperature in the Langevin dynamics (see $\beta$ in article).

(7) Gamma: Value of the friction coefficient in the Langevin dynamics (see $\gamma$ in the article).  
