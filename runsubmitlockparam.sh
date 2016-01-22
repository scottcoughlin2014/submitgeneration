#!/bin/bash

module load python
module load mpi/openmpi-1.8.3-intel2015.0

source /projects/b1011/ligo_project/lsc/o1_lalinference_20151210/etc/lscsoftrc

rm -rf ./condor

for event in 11 44;
	do for combo in "none" "skyloc" "skyloc_dist" "skyloc_thetajn" "skyloc_thetajn_dist";
                do if [ ${combo} = "skyloc" ]; then
                   addflags="--fix-rightascension --fix-declination"
                elif [ ${combo} = "skyloc_dist" ]; then
                   addflags="--fix-rightascension --fix-declination --fix-distance"
                elif [ ${combo} = "skyloc_thetajn" ]; then
                   addflags="--fix-rightascension --fix-declination --fix-costheta_jn"
                elif [ ${combo} = "skyloc_thetajn_dist" ]; then
                   addflags="--fix-rightascension --fix-declination --fix-costheta_jn --fix-distance"
                else
                   addflags=""
                fi;
		./lalinference_mcmc_submit_new.py --inj /projects/b1011/spinning_runs/STT4injections.xml --event ${event} --approx SpinTaylorT4 --lowM1 3.0 --lowM2 1.0 --dir /projects/b1011/spinning_runs/freezingparams/${event}/${combo} ${addflags};
	done;
done;
# --fix-distance --fix-rightascension --fix-declination --fix-costheta_jn
