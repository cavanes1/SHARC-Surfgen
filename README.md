# SHARC-Surfgen
Here are the very basics of using [SHARC](https://github.com/sharc-md/sharc) with COLUMBUS.

## Installation
1. `wget https://github.com/sharc-md/sharc/archive/refs/tags/v3.0.2.tar.gz` downloads SHARC
2. `cd source`
3. `vi Makefile` and change GNU to Intel and mkl for ARCH
4. `ml intel/2022.2-dpcpp` load Intel module
5. `make`
6. `make install`
7. Add `source ~/sharc/bin/sharcvars.sh` to .bash_profile

## Tutorial
This tutorial was provided to me by [Yuchen Wang](https://github.com/ywang312).

The files are located at .../scr16_dyarkon1/public/traj-sharc-tutorial

> This is a short tutorial to run sharc trajectory with different initial conditions

> There are three directories, trajdata for storing trajectory results,
> sharc-main is the main directory to run sharc. The QM(quantum mechanical) subdirectory inside it contains all PES and dipole/soc surfaces (if any)
> sharc-surfgen is the directory to compile the small program that reads PES and dipole/soc surfaces

> First, generate initial conditions. Strictly, the initial conditions are generated usually from Wigner distribution,
> the sharc website tutorial contain how to do this, I also have every input under sharc-main (the geometry and molden file),
> however, I want you to run a short test of pure random intial velocity (I want to see the difference it makes by comparing with my result)
> The way to do this is simply go inside sharc-main/inp, which contains the geninp.f90, change the value in line 10 to be "random x",
> where x is the value to vary, I recommend from 0.1 to 1 (this value is the energy on each atom, so the kinetic energy is 4x because NH3 has 4 atoms)
> and then run ./geninp.sh, this compiles geninp.f90 and make 100 initial conditions, remember to change 100 to at least 1000 in the bash script as 100 is too small for the ambition of PES.

> Second is to run the run.sh inside sharc-main. Also change the i value in run.sh to be the number of initial conditions in the first step.

> Third is to submit the job. Adjust the time accordingly.

> The command you needed (probably) are:

> `cd sharc-main/inp`

> #change the initial conditions here
> `vi geninp.f90`

> `./geninp.sh`

> #go back to sharc-main
> `cd ..`

> #change the number of trajectory here
> `vi run.sh`

> #change the run time of slurm script here
> `vi submit.sl`

> `sbatch submit.sl`

> When the trajectories are done, go into trajdata and analyze the result with `./collectgeominfo.sh`,
> this returns the state information when the bond distance are greater than a certain value.
> state1 and stata2 corresponds to the mch and diagonal representation (see sharc paper for this),
> since we use a adiabatic representation here, they should be equal.
> count the number of "State 1: 1" appears (I use grep). that's the ground state product ratio.
> and so for the excited state product. Note there is a chance (especially with lower energy) that it does not dissociate. That's why we check the bond distance.

> It is hoped that you can understand the scripts after you go through the entire process and make improvements!

> Note: You have to install sharc and have $SHARC point to the sharc program. You can modify and compile the f90 files in sharc-surfgen for other surfgen PES.
