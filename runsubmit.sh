#!/bin/bash

module load python
module load mpi/openmpi-1.8.3-intel2015.0

source /projects/b1011/ligo_project/lsc/o1_lalinference_20151210/etc/lscsoftrc

rm -rf ./condor

for event in 11 44;
		do ./lalinference_mcmc_submit_new.py --inj /projects/b1011/ligo_project/full_bns.xml --event ${event} --approx SpinTaylorT2 --lowM1 1.0 --lowM2 1.0 --dir /projects/b1011/spinning_runs/BNS/${event}/ --comp-min 1.0 --comp-max 5.0
done;
# --fix-distance --fix-rightascension --fix-declination --fix-costheta_jn
