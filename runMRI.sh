#!/bin/bash

module load python
module load mpi/openmpi-1.8.3-intel2015.0

source /projects/b1011/ligo_project/lsc/kagra_o2_lalinference/etc/lscsoftrc

rm -rf ./condor

for event in 1 2 3;
        do ./lalinference_mcmc_submit_new.py --inj /projects/b1011/mri_runs/coinc.xml.gz --event ${event} --approx SpinTaylorT4threePointFivePN --lowM1 1.0 --lowM2 1.0 --walltime 21:00:00:00 --dir /projects/b1011/mri_runs/lalinference_mcmc_spin/scottyTry//${event}/ --ifo 'H1' 'L1'
done;
