#!/bin/bash

module load python
module load mpi/openmpi-1.8.3-intel2015.0

source /projects/b1011/ligo_project/lsc/freezingparams_20160402/etc/lscsoftrc

rm -rf ./condor

for event in 11 44;
		do ./lalinference_mcmc_submit_new.py --inj /projects/b1011/ligo_project/full_bns.xml --event ${event} --approx SpinTaylorT4 --lowM1 1.0 --lowM2 1.0 --dir /projects/b1011/spinning_runs/BNS/STT4/${event}/ --comp-min 1.0 --comp-max 5.0 --ppall --walltime 30:00:00:00
done;
# --fix-distance --fix-rightascension --fix-declination --fix-costheta_jn
