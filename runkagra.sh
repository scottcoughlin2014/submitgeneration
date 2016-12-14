#!/bin/bash

module load python
module load mpi/openmpi-1.8.3-intel2015.0

source /projects/b1011/ligo_project/lsc/kagra_o2_lalinference/etc/lscsoftrc

rm -rf ./condor

for event in 11 23 44 53 63 67 85 98 102 105 112 118 124 129 132 138 171 194 218 236 240 243 244 262 273 281 292 320 339 346 359 370 376 378 386 391 392 404 416 423 450 469 484 486 498 529 536 553 566 572 576 610 636 640 647 650 654 655 666 677 680 698 716 729 733 770 775 796 802 804 820 821 824 843 848 851 876 881 896 901 906 908 916 921 935 946 953 973 979 981 985;
	./lalinference_mcmc_submit_new.py --inj /projects/b1011/spinning_runs/IMRfreezinginj.xml --event ${event} --approx IMRPhenomPv2pseudoFourPN --lowM1 1.0 --lowM2 1.0 --plot-2d --walltime 30:00:00:00 --dir /projects/b1011/kagra/kagra_o2_lalinference/${event}/ --ifo 'H1' 'L1' 'V1' 'K1'
# --fix-distance --fix-rightascension --fix-declination --fix-costheta_jn
