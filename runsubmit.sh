#!/bin/bash

module load python
module load mpi/openmpi-1.8.3-intel2015.0

source /projects/b1011/ligo_project/lsc/o1_lalinference_20151210/etc/lscsoftrc

./lalinference_mcmc_submit_new.py --inj /projects/b1011/spinning_runs/STT2injections.xml --event 11 --approx SpinTaylorT4 --lowM1 1.0 --lowM2 5.0 --no-detector-frame
