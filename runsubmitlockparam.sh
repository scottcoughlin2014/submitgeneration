#!/bin/bash

module load python
module load mpi/openmpi-1.8.3-intel2015.0

source /projects/b1011/ligo_project/lsc/o1_lalinference_20160402/etc/lscsoftrc

rm -rf ./condor

for event in 11 23 44 53 63 67 85 98 102 105 112 118 124 129 132 138 171 194 218 236 240 243 244 262 273 281 292 320 339 346 359 370 376 378 386 391 392 404 416 423 450 469 484 486 498 529 536 553 566 572 576 610 636 640 647 650 654 655 666 677 680 698 716 729 733 770 775 796 802 804 820 821 824 843 848 851 876 881 896 901 906 908 916 921 935 946 953 973 979 981 985;
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
		./lalinference_mcmc_submit_new.py --inj /projects/b1011/spinning_runs/STT4injections.xml --event ${event} --approx SpinTaylorT4 --lowM1 9.0 --lowM2 1.0 --plot-2d --walltime 30:00:00:00 --dir /projects/b1011/spinning_runs/freezingparams_20160402/${event}/${combo} ${addflags} --compare;
	done;
done;
# --fix-distance --fix-rightascension --fix-declination --fix-costheta_jn
