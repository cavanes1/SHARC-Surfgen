# SHARC-Surfgen
This program was used for an upcoming research paper that will be published in the Journal of Chemical Physics.
An earlier version was used for [this research paper](https://doi.org/10.1021/acs.jpclett.4c02312) published in the Journal of Physical Chemistry Letters.
This code allows the [SHARC](https://sharc-md.org/) molecular dynamics program to use quantum chemistry data from potential energy surfaces constructed by Surfgen.
This repository includes the ammonia singlet and triplet surfaces generated using neural networks by the [Yarkony Group](https://github.com/yarkonygrp/).
Basic instructions are shown below for convenience, but see the [wiki](https://github.com/cavanes1/SHARC-Surfgen/wiki) for more details.

## Usage

1. Use wigner.py by running `$SHARC/wigner.py sharc-main/molden.freq -n 5000` (you may need more than 5000 if using energy bins)
2. Load the needed module `ml intel/2020.1` and compile in sharc-surfgen without the electric field with `make`. This step must be performed now if using energy bins. Otherwise, you can wait until step 6.
3. Make sure to `chmod +x sharc-surfgen/surfgen.x sharc-main/run.sh sharc-main/QM/runQM.sh sharc-main/QM/inQM.sh` if these do not already have executable permissions
5. Edit geninp.py as needed `vi geninp.py`
4. Run geninp.py using `sbatch geninp.sl`
6. If using an electric field, compile sharc-surfgen with the field enabled
7. Submit run.sh with `sbatch submit.sl` or `python parallel.py` in order to run SHARC
8. Make sure the desired number of trajectories actually completed using `grep finish slurm-* | wc -l`
9. Account for bad trajectories in the trajdata directory with `grep trajdata slurm-* | sort -u` and, if necessary, running `python errorcheck.py`
10. Perform analysis using the tools in the trajdata directory via `sbatch collectgeominfo.sh`
11. To clear up space after data is collected, delete the QM directories, e.g., `rm -r trajdata/traj*/QM`
