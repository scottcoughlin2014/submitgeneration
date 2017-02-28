#!/bin/bash

module load python
module load mpi/openmpi-1.10.5-intel2015.0

source /projects/b1011/ligo_project/lsc/o2_lalinference_GW170104_MPI_1_10_5/etc/lscsoftrc

rm -rf ./condor

for event in 11 23 44 53 63 67 85 98 102 112 118 124 129 132 171 194 218 236 244 262 273 281 292 339 359 376 378 386 392 416 423 450 484 486 498 536 553 566 572 610 636 640 647 654 655 666 677 680 716 733 770 796 802 804 820 821 824 848 851 881 896 901 906 908 916 921 935 953 973 979 981 985;
	do ./lalinference_mcmc_submit_new.py --inj /projects/b1011/spinning_runs/IMRfreezinginj.xml --event ${event} --approx IMRPhenomPv2pseudoFourPN --lowM1 2.0 --lowM2 2.0  --walltime 21:00:00:00 --dir /projects/b1011/kagra/HLVK/${event}/ --ifo 'H1' 'L1' 'V1' 'K1'  --plot-2d
done;
for event in 105 138 240 243 320 346 370 391 404 469 529 576 650 698 729 843 775 876 946;
	do ./lalinference_mcmc_submit_new.py --inj /projects/b1011/spinning_runs/IMRfreezinginj.xml --event ${event} --approx IMRPhenomPv2pseudoFourPN --lowM1 2.0 --lowM2 2.0  --walltime 21:00:00:00 --dir /projects/b1011/kagra/HLVK/${event}/ --ifo 'H1' 'L1' 'V1' 'K1' --distance-max 2000 --plot-2d
done;
# --fix-distance --fix-rightascension --fix-declination --fix-costheta_jn
