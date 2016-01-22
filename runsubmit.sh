#!/bin/bash

module load python
module load mpi/openmpi-1.8.3-intel2015.0

source /projects/b1011/ligo_project/lsc/o1_lalinference_20151210/etc/lscsoftrc

rm -rf ./condor

for event in {0 11 44 101};
	do for combo in "none" "skyloc" "skyloc_dist" "skyloc_thetajn" "skyloc_thetajn_dist";
		if [ ${combo} = "skyloc" ]; then
                   addflags = "--fix-rightascension --fix-declination"
		else if [ ${combo} = "skyloc_dist" ]; then
		   addflags = "--fix-rightascension --fix-declination --fix-distance"
                else if [ ${combo} = "skyloc_thetajn_dist" ]; then
		   addflags = "--fix-rightascension --fix-declination --fix-costheta_jn"
                else if [ ${combo} = "skyloc__thetajn_dist" ]; then
                   addflags = "--fix-rightascension --fix-declination --fix-costheta_jn --fix-distance"
                else
		   addflags = ""
                fi
		do ./lalinference_mcmc_submit_new.py --inj /projects/b1011/spinning_runs/STT4injections.xml --event ${event} --approx SpinTaylorT4 --lowM1 5.0 --lowM2 1.0 --dir /projects/b1011/spinning_runs/freezingparams/${event}/${combo} ${addflags};
		done;
# --fix-distance --fix-rightascension --fix-declination --fix-costheta_jn
