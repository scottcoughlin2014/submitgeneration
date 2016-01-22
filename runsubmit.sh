#!/bin/bash

module load python
module load mpi/openmpi-1.8.3-intel2015.0

source /projects/b1011/ligo_project/lsc/o1_lalinference_20151210/etc/lscsoftrc

./lalinference_mcmc_submit_new.py --inj /projects/b1011/spinning_runs/STT4injections.xml --event 11 --approx SpinTaylorT4 --lowM1 5.0 --lowM2 1.0 --dir /projects/b1011/scoughlin/test/ 
# --fix-distance --fix-rightascension --fix-distance --fix-costheta_jn
