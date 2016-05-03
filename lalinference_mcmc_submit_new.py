#!/usr/bin/env python
import lalsimulation as lalsim

import os
import glob
import stat
import argparse
import getpass
import socket
import numpy as np

import math
import subprocess, shlex

from glue.ligolw import lsctables
from glue.ligolw import ligolw
from pylal import SimInspiralUtils

from lalsim_snr import *


def power_log(x):
    return 2**(math.ceil(math.log(x, 2)))

def exists(filename):
    """Check if filename exists"""
    try:
        with open(filename): return True
    except IOError:
        print "Warning: {} not found.".format(rc)
        return False

def check_for_arg_substring(arg_str, args):
    arg_check = ['-{}'.format(arg_str) in arg for arg in args]
    return True if True in arg_check else False

def num_of_arg_substring(arg_str, args):
    arg_check = ['-fix' in arg for arg in args]
    return sum(arg_check)

parser = \
  argparse.ArgumentParser(description='Generate a submit file for \
                          lalinference_mcmc on grail.')

msub = parser.add_argument_group('MSUB')
env = parser.add_argument_group('env')
li_mcmc = parser.add_argument_group('lalinference_mcmc')

msub.add_argument('--alloc', default='b1011',
        help='Allocation to charge SUs to (default=b1011).')
msub.add_argument('--queue', default='ligo',
        help='Queue for job (default=ligo).')
msub.add_argument('--jobName',
        help='Name of job, used for output file names and queue listing.')
msub.add_argument('--dir', default=os.getcwd(),
        help='Directory where submit script is written and \
              executed from (default=current directory).')
msub.add_argument('--dep', default=None,
        help='Wait for dependent job ID to complete (default=None).')
msub.add_argument('--name', default='submit',
        help='Name of submit file (default=submit).')
msub.add_argument('--walltime', default='2:00:00:00',
        help='Walltime for job (default=2:00:00:00).')
msub.add_argument('--cores-per-node', type=int, default=16,
        help='Number of cores per node (default=16).')
msub.add_argument('--multiple-nodes', default=False, action='store_true',
        help='If nChains > 16 then use more than one node.')
msub.add_argument('--nPar', default=None, type=int,
        help='Number of dimensions for MCMC.  Defaults for common templates \
              are set, assuming no PSD fitting.')
msub.add_argument('--homepath', \
        default = os.getcwd(),
        help='Where the compilation of msub commands should go')
msub.add_argument('--email', default = 'scottcoughlin2014@u.northwestern.edu',
        help='Email for when job starts, ends, or is aborted suddenly.')
msub.add_argument('--emailyes', default=False, action='store_true',
        help='Do you want email alert when job is finished?.')

env.add_argument('--branch', default='o1_lalinference_20160402',
        help='Branchname to use, assuming \
              /projects/p20251/USER/lsc/BRANCHNAME/etc/lscsoftrc \
              exists (default=master).')
env.add_argument('--rc', action='append',
        help='Specify direct path to rc files to be sourced (e.g. lscsoftrc). \
              /projects/b1011/non-lsc/lscsoft-user-env.sh added by default.')
env.add_argument('--sim-quest', default=False, action='store_true',
        help='Act as if on Quest.  Useful for setting up submit files on local\
              computer for uploading to Quest')
env.add_argument('--quest-iwd', help='Working directory for Quest submit script.')

li_mcmc.add_argument('--era', default='advanced',
        help='Era ("initial" or "advanced") of detector PSD for SNR \
              calculation. If no cache arguments given, this will add the \
              appropriate analytical PSD arguments to the submit file.')
li_mcmc.add_argument('--resume', default=False, action='store_true',
        help='What to resume?')
li_mcmc.add_argument('--psd', nargs='+',
        help='Pre-computed PSD(s), either as a single xml, or list of ascii\
              files, one for each IFO in the order the IFOs were specified. \
              XMLs are converted using lalinference_pipe_utils.')
li_mcmc.add_argument('--ifo', nargs='+', default=['H1', 'L1', 'V1'],
        help='IFOs for the analysis.')
li_mcmc.add_argument('--inj',
        help='Injection XML file.')
li_mcmc.add_argument('--event', type=int,
        help='Event number in XML to inject.')
li_mcmc.add_argument('--trigtime', type=str,
        help='Trigger time of event.  Automatically set when injecting.')
li_mcmc.add_argument('--approx', default='SpinTaylorT4',
        help='Specify a template approximant (default SpinTaylorT4).')
li_mcmc.add_argument('--amporder', default=None,
        help='Specify amplitude order of template.')
li_mcmc.add_argument('--flow', type=float,
        help='Lower frequency bound for all detectors (default=40).')
li_mcmc.add_argument('--lowM1', type=float,
        help='Lowest reasonable M1 of the binary your are searching for. Used to determine appropraite seglength.')
li_mcmc.add_argument('--lowM2', type=float,
        help='Lowest reasonable M2 of the binary your are searching for. Used to determine appropraite seglength.')
li_mcmc.add_argument('--fhigh', type=float,
        help='Upper frequency bound for all detectors, given as a fraction\
                of the injection\'s ISCO frequency.')
li_mcmc.add_argument('--srate', default=None, type=float,
        help='Sampling rate of the waveform.  If not provided and an injection\
              is peformed, it is set to be sufficient for the signal being \
              injected.  If no injection, it defaults to a sufficient value \
              for a 1.4-1.4 binary coalescence (expensive!).')
li_mcmc.add_argument('--seglen', default=None, type=float,
        help='Length of data segment used for likelihood compuatation. \
              Same default behavior as "--srate".')
li_mcmc.add_argument('--psdlength', default=None, type=float,
        help='Length of data segment to use for PSD estimation. \
              Defaults to 32*seglen.')
li_mcmc.add_argument('--psdstart', default=None, type=float,
        help='GPS time to start PSD calculation. \
              Defaults to trigtime - psdlength - seglen')
li_mcmc.add_argument('--temp-ladder-top-down', default=False, action='store_true',
        help='Build the temperature ladder from the bottom up, using an \
              analytic prescription for the spacing that should ensure \
              communication between chains.  Sets the number of cores so \
              that the hottest temperature should be sampling the prior.')
li_mcmc.add_argument('--noisy', default=False, action='store_true',
        help='Use a non-zero noise realization.')
li_mcmc.add_argument('--distance-max', type=float,
        help='Hard outer prior boundary on distance (default=1 Gpc iLIGO, 2 Gpc aLIGO).')
li_mcmc.add_argument('--no-malmquist', default=False, action='store_true',
        help='Do not use the Malmquist prior, which by default approximates the\
                selection effects imposed by the detection processes as a cut\
                in the SNR in the second loudest detector.')
li_mcmc.add_argument('--no-margtimephi', default=False, action='store_true',
        help='Do not use the time and phase marginalized likelihood function.')
li_mcmc.add_argument('--trigSNR', default=None, type=float,
        help='SNR of the trigger (calculated automatically if injecting).')
li_mcmc.add_argument('--tempMin', default=1.0, type=float,
        help='Temperature of coldest chain (default=1.0).')
li_mcmc.add_argument('--tempMax', type=float,
        help='Temperature of hotest chain.  Determined automatically if \
              injecting, or trigSNR is given.')
li_mcmc.add_argument('--Neff', default=1000, type=int,
        help='Requested number of independent samples.')
li_mcmc.add_argument('--Niter', default=1000000000, type=int,
        help='Maximum number of MCMC iterations to allow.')
li_mcmc.add_argument('--fix-rightascension', default=False, action='store_true',
	help='Fix RA')
li_mcmc.add_argument('--fix-declination', default=False, action='store_true',
        help='Fix DEC')
li_mcmc.add_argument('--fix-distance', default=False, action='store_true',
        help='Fix Distance')
li_mcmc.add_argument('--fix-costheta_jn', default=False, action='store_true',
        help='Fix costheta_jn')
li_mcmc.add_argument('--ppall', default=False, action='store_true',
        help='Post Process ALL chains')
li_mcmc.add_argument('--compare', default=False, action='store_true',
        help='Generate Script that will make comparison pages of all runs you just generated')
li_mcmc.add_argument('--plot-2d', default=False, action='store_true',
        help='Want 2D plots?')

##########################################################################
##########################################################################
################             MAIN                #########################
##########################################################################
##########################################################################

args,unknown = parser.parse_known_args()


# Assume all unknown arguments are meant for lalinference_mcmc
li_args = unknown

# If not on quser##, don't bother making a submit file
if 'quser' in socket.gethostname() or args.sim_quest:
    on_quest = True
else:
    on_quest = False

# Join name of submit file to PWD.

if args.sim_quest:
    if args.quest_iwd:
        out_dir = args.quest_iwd
        submitFilePath = os.path.join(args.dir, args.name)
    else:
        out_dir = os.path.abspath(args.dir)
        submitFilePath = os.path.join(out_dir, args.name)
else:
    out_dir = os.path.abspath(args.dir)
    submitFilePath = os.path.join(out_dir, args.name)

# Setup output directory for post processing
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

webdir = "%s/post"%out_dir
if not os.path.isdir(webdir):
   os.makedirs(webdir)

# Setup and check envirnment files to be sourced
rcs = args.rc if args.rc else []

# Add non-lsc standard location if one is not given
non_lsc_check = ['non-lsc' in rc_path for rc_path in rcs]
if True not in non_lsc_check and on_quest:
    rcs.append('/projects/b1011/non-lsc/lscsoft-user-env.sh')

# Assume a certain directory tree structure if branch name is given
if args.branch and on_quest:
    try:
        lscsoftrc = '/projects/b1011/ligo_project/lsc/{}/etc/lscsoftrc'.format(args.branch)
    except KeyError:
        lscsoftrc = '/projects/b1011/ligo_project/lsc/o1_lalinference_20160402/etc/lscsoftrc'

    rcs.append(lscsoftrc)

# Only include rc files that exist
if not args.sim_quest:
    rcs[:] = [rc for rc in rcs if exists(rc)]

# Necessary modules
modules = ['python','mpi/openmpi-1.8.3-intel2015.0']
unload_modules = []

# Determine which simulated PSD to use.
noise_psd_caches = {}
psd_files=None
print "Using {}-era PSDs.".format(args.era)
if args.era == 'initial':
    for _ifos, _cache in (
          (('H1', 'H2', 'L1', 'I1'), 'LALLIGO'),
          (('V1',), 'LALVirgo')):
            for _ifo in _ifos:
                noise_psd_caches[_ifo] = _cache

elif args.era == 'advanced':
    for _ifos, _cache in (
          (('H1', 'H2', 'L1', 'I1'), 'LALSimAdLIGO'),
          (('V1',), 'LALSimAdVirgo')):
            for _ifo in _ifos:
                noise_psd_caches[_ifo] = _cache


# Check if caches specified in extra arguments
caches_specified = check_for_arg_substring('cache', li_args)

# Check for params already marginalized by likelihood
if not args.no_margtimephi:
    n_marg_params = 2
elif check_for_arg_substring('margtime', li_args) or check_for_arg_substring('margphi', li_args):
    n_marg_params = 1
else:
    n_marg_params = 0

# Determine approximant used
approx = lalsim.GetApproximantFromString(args.approx)

spin_support = lalsim.SimInspiralGetSpinSupportFromApproximant(approx)

# User requests for spin parameter space have priority over waveform physics
# Check for no-spin flag
no_spin = check_for_arg_substring('noSpin', li_args) or \
          check_for_arg_substring('disable-spin', li_args)
if not no_spin:
    no_spin = spin_support == lalsim.SIM_INSPIRAL_SPINLESS

# Check for spin aligned flag
spin_aligned = check_for_arg_substring('spinAligned', li_args)
if not spin_aligned:
    spin_aligned = spin_support == lalsim.SIM_INSPIRAL_ALIGNEDSPIN

# Check for single spin flag
single_spin = check_for_arg_substring('singleSpin', li_args)
if not single_spin:
    single_spin = spin_support == lalsim.SIM_INSPIRAL_SINGLESPIN

# Figure out the possible number of parameters, without fixing or marginalizing
if no_spin:
    nPar = 9
elif spin_aligned:
    nPar = 11
elif single_spin:
    nPar = 12
else:
    nPar = 15

# Maximum sampling rate
srate_max = 16384

# Determine sampling rate, segment length, and SNR (--trigSNR takes precedence).
if args.amporder is None:
    amp_order = 0
else:
    try:
        amp_order = int(args.amporder)
    except ValueError:
        amp_order = lalsim.GetOrderFromString(args.amporder)

SNR = None
calcSNR = False if args.trigSNR else True

# Default to different starting frequencies for different eras:
if args.flow is not None:
    temp_flow = args.flow
elif args.era is 'advanced':
    temp_flow = 20.0
elif args.era is 'initial':
    temp_flow = 40.0
else:
    raise RuntimeError("No lower frequency bound provided.")

flow = temp_flow
# set fref to flow
fref = flow

# Default to different starting max distances for different eras:
if args.distance_max is not None:
    distance_max = args.distance_max
elif args.era is 'advanced':
    distance_max = 6000
elif args.era is 'initial':
    distance_max = 1000
else:
    raise RuntimeError("No lower frequency bound provided.")

# Keep track of fixed params
n_fixed_params = 0

# Determine trigger time (save precession read in)
trigtime = None
fhigh = None
if args.trigtime is not None:
    trigtime_as_string = args.trigtime
    trigtime = float(args.trigtime)

elif args.inj and args.event is not None:
    event = SimInspiralUtils.ReadSimInspiralFromFiles([args.inj])[args.event]
    trigtime = event.geocent_end_time + 1e-9*event.geocent_end_time_ns
    trigtime_as_string = str(trigtime)
    simidtmp = str(event.simulation_id)
    simid = simidtmp.split(':')[2]
    
    # Determine if fixed parameters have been asked for
    fixargs = ''
    if args.fix_rightascension:
		RA = event.longitude
		fixargs = fixargs + '  --fix-rightascension {} '.format(RA)
		n_fixed_params +=1

    if args.fix_declination:
		DEC = event.latitude
		fixargs = fixargs + '  --fix-declination {} '.format(DEC)
		n_fixed_params +=1

    if args.fix_distance:
		Dist = event.distance
		fixargs = fixargs + '  --fix-logdistance {} '.format(np.log(Dist))
		n_fixed_params +=1

    if args.fix_costheta_jn:
		commandline = './get_injection.py --inj {0} --event {1}'.format(args.inj,args.event)
		print('to get costheta_jn we must do a special command: {}'.format(commandline))
		proc = subprocess.Popen(shlex.split(commandline), stdout=subprocess.PIPE)
		output = proc.stdout.read()
		costheta_jn = output.split('costheta_jn')[1]
		costheta_jn = costheta_jn.split(':')[1]
		costheta_jn = costheta_jn.split('injected')[0]
		costheta_jn = costheta_jn.split('\n')[0]
		fixargs = fixargs + '  --fix-costheta_jn {} '.format(costheta_jn)
		n_fixed_params +=1

    # Determine upper frequency cutoff based on ISCO if requested
    if args.fhigh:
		f_isco = ISCO(event.mass1, event.mass2)
		fhigh = args.fhigh * f_isco

# Calculate any arguments not specified
if calcSNR:
    commandline = "lalapps_chirplen --flow {0} --m1 {1} --m2 {2}".format(flow,args.lowM1,args.lowM2) 
    print('Command that will determine some of our parameters = {}'.format(commandline))
    proc = subprocess.Popen(shlex.split(commandline), stdout=subprocess.PIPE)
    output = proc.stdout.read()
    dur = float(output.split()[19])
    fhighend  = float(output.split()[12])
    srate = power_log(2*fhighend)
    seglen = power_log(dur)

    print('F Low of waveform = {}'.format(float(output.split()[5])))
    print('F High of Waveform = {}'.format(fhighend))
    print('Sampling Rate based on FHigh = {}'.format(srate))
    print('Duration of Waveform = {}'.format(dur))
    print('Segment Length of waveform based on Duration = {}'.format(seglen))

    commandline = './pyburst_inj_snr {} --det-psd-func=H1=SimNoisePSDaLIGOZeroDetHighPower --det-psd-func=L1=SimNoisePSDaLIGOZeroDetHighPower --det-psd-func=V1=SimNoisePSDAdvVirgo  --low-frequency-cutoff={} --nyquist-frequency={} --waveform-length={}  --sim-id {} --skip-coherent-snr --no-print-single-ifo'.format(args.inj,flow,power_log(fhighend),seglen,simid)
    print('Command to find Network SNR so we can set temperature ladders: ' +commandline)
    
    proc = subprocess.Popen(shlex.split(commandline), stdout=subprocess.PIPE)
    output = proc.stdout.read()
    print(output)
    SNR = float(output.split(':')[1])

if args.trigSNR:
    SNR = args.trigSNR

if args.srate:
    srate = args.srate

if args.seglen:
    seglen = args.seglen

if srate > srate_max:
    srate = srate_max

# Use zero-noise realization by default
if not args.noisy:
    li_args.append('--0noise')
    li_args.append('--dataseed 0')


# Determine number of parameters
if args.nPar:
    nPar = args.nPar
else:
    nPar -= n_fixed_params
    nPar -= n_marg_params

    print "nPar: {}".format(nPar)

# Target likelihood value "sampled" by the hottest chain.
# Same as in lalinference_mcmc.
target_hot_like = nPar/2.

#tmin = 1.0
#tmax = SNR**2 / nPar
#print tmax

#delta_t = 1 + math.sqrt(2./nPar)
#print delta_t

#t, i = tmin, 0
#while t < tmax:
#    i += 1
#    t = tmin * delta_t ** i

#print i

# Ladder spacing flat in the log.  Analytic delta
if not args.temp_ladder_top_down:

    # Determine maximum temperature
    temp_max = None
    if args.tempMax:
        temp_max = args.tempMax

    elif SNR is not None:
        max_log_l = SNR*SNR/2
        temp_max = max_log_l / target_hot_like
    print "Max temperature of {} needed.".format(temp_max)

    # Determine spacing
    temp_delta = 1 + np.sqrt(2./nPar)

    n_chains = 1
    temp_min = args.tempMin
    temp = temp_min
    while temp < temp_max:
        n_chains += 1
        temp = temp_min * np.power(temp_delta, n_chains)
    print "{} temperatures needed.".format(n_chains)

    ppn = args.cores_per_node
    if n_chains > ppn and args.multiple_nodes:
        n_nodes = int(np.ceil(n_chains/ppn))
        n_cores = ppn
    else:
        n_nodes = 1
        n_cores = n_chains if n_chains < ppn else ppn

# Use malmquist prior by default
if not args.no_malmquist:
    li_args.append('--malmquistprior')

# Prepare lalinference_mcmc arguments
ifos = args.ifo
ifo_args = ['--ifo {}'.format(ifo) for ifo in ifos]
flow_args = ['--{}-flow {:g}'.format(ifo, flow) for ifo in ifos]

if fhigh:
    fhigh_args = ['--{}-fhigh {:g}'.format(ifo, fhigh) for ifo in ifos]
else:
    fhigh_args = None

# PSD-related args
psd_args = ''
 
psdlength = args.psdlength if args.psdlength else 32*seglen
psd_args += '--psdlength {}'.format(psdlength)

psdstart = args.psdstart if args.psdstart else trigtime-psdlength-seglen
psd_args += ' --psdstart {}'.format(psdstart)

if psd_files is not None:
    if not caches_specified:
        cache_args = ['--{}-cache interp:{}'.format(ifo, psd_files[ifo]) for ifo in ifos]
        caches_specified = True

    psd_args += ' --psd [{}]'.format(','.join([psd_files[ifo] for ifo in ifos]))

if not caches_specified:
    cache_args = \
        ['--{}-cache {}'.format(ifo, noise_psd_caches[ifo]) for ifo in ifos]

# Want to resume a run?
if args.resume:
	commandline = 'grep "seed" {0}'.format(glob.glob(submitFilePath + '.*')[0])
        proc = subprocess.Popen(shlex.split(commandline), stdout=subprocess.PIPE)
        output = proc.stdout.read()
	randomseed = output.split('seed:')[1]
	randomseed = randomseed.split('\n')[0]
	print('randomseed = {0}'.format(randomseed))
	resumeargs= '  --resume --randomseed {0}'.format(randomseed)

# Specify number of cores on the command line
runline = 'mpirun -n {} lalinference_mcmc'.format(n_chains)

with open(submitFilePath,'w') as outfile:
    outfile.write('#!/bin/bash\n')
    if on_quest:
        # MSUB directives
        outfile.write('#MSUB -A {}\n'.format(args.alloc))
        outfile.write('#MSUB -q {}\n'.format(args.queue))

        outfile.write('#MSUB -l walltime={}\n'.format(args.walltime))
        outfile.write('#MSUB -l nodes={}:ppn={}\n'.format(n_nodes,n_cores))

        if args.dep:
            outfile.write('#MSUB -l {}\n'.format(args.dep))

        if args.jobName:
            outfile.write('#MSUB -N {}\n'.format(args.jobName))

        # Give read permissions to screen output
        outfile.write('#MOAB -W umask=022\n')

        # Write stdout and stderr to the same file
        outfile.write('#MSUB -j oe\n')

        # Job working directory
        outfile.write('#MSUB -d {}\n'.format(out_dir))

        # Ensure core dump on failure
        outfile.write('ulimit -c unlimited\n')

        # Load modules
        outfile.writelines(['module load {}\n'.format(module) 
            for module in modules])
        outfile.write('\n')

        # Manually unload Active-python that is loaded by the intel MPI
        outfile.writelines(['module unload {}\n'.format(module) 
            for module in unload_modules])
        outfile.write('\n')
	
	# Add gcc 4.8.3 as required by quest as of November 2015
#	outfile.write('module load gcc/4.8.3')
#	outfile.write('\n')	


    # Load LALSuite environment
    outfile.writelines(['source {}\n'.format(rc) for rc in rcs])
    outfile.write('\n')

    # Add OpenMP environment variable
	
    outfile.write('export OMP_NUM_THREADS=1\n')


    # lalinference_mcmc command line
    outfile.write('{}\\\n'.format(runline))
    outfile.write('  {}\\\n'.format(' '.join(ifo_args)))
    if not caches_specified:
        outfile.write('  {}\\\n'.format(' '.join(cache_args)))
    outfile.write('  {}\\\n'.format(' '.join(flow_args)))
    if fhigh:
        outfile.write('  {}\\\n'.format(' '.join(fhigh_args)))

    if args.inj and args.event is not None:
        outfile.write('  --inj {} --event {}\\\n'.format(
            args.inj, args.event))

    if trigtime is not None:
        outfile.write('  --trigtime {}\\\n'.format(trigtime_as_string))

    outfile.write('  {}\\\n'.format(psd_args))
    outfile.write('  --srate {:g}\\\n'.format(srate))
    outfile.write('  --seglen {:g}\\\n'.format(seglen))
    outfile.write('  --approx {}\\\n'.format(args.approx))
    if args.resume:
	outfile.write('{}\\\n'.format(resumeargs))

    if amp_order is not None:
        outfile.write('  --amporder {}\\\n'.format(amp_order))

    if not args.no_margtimephi:
        outfile.write('  --margtimephi \\\n')
	outfile.write('  --no-detector-frame \\\n')

    outfile.write('  --distance-max {}\\\n'.format(distance_max))
    outfile.write('  --neff {}\\\n'.format(args.Neff))
    outfile.write('  --fref {0} --inj-fref {0}\\\n'.format(fref))
    if fixargs is not '':
        outfile.write('{}\\\n'.format(fixargs))
    outfile.write('  {}'.format('\\\n  '.join(li_args)))
    outfile.write('\n')
    outfile.write('\n')
    outfile.write('module load gcc/4.8.3\n')
    outfile.write('\n')
    # Post-processing command line

    outfile.write('cbcBayesPostProc.py --lalinfmcmc -i {} --event {} --outpath={} -d {}/PTMCMC.output.*.00 --dievidence --skyres=.5 --deltaLogL {}\n'.format(args.inj,args.event,webdir,out_dir,target_hot_like))
    outfile.write('\n')

# Create separate pp.sh in order to manually run PP when needed
ppFilePath = os.path.join(out_dir, 'pp.sh')
with open(ppFilePath,'w') as ppfile:
    	ppfile.write('#MSUB -A {}\n'.format(args.alloc))
    	ppfile.write('#MSUB -q {}\n'.format(args.queue))

    	ppfile.write('#MSUB -l walltime=00:08:00:00\n')
    	ppfile.write('#MSUB -l nodes=1:ppn=1\n')

    	ppfile.write('#MSUB -N fullPP\n')

    	# Give read permissions to screen output
    	ppfile.write('#MOAB -W umask=022\n')

    	# Write stdout and stderr to the same file
    	ppfile.write('#MSUB -j oe\n')

    	# Job working directory
    	ppfile.write('#MSUB -d {}\n'.format(args.homepath))

        ppfile.write('ulimit -c unlimited\n')
        ppfile.write('\n')

        ppfile.write('module load python\n')
        ppfile.write('module load mpi/openmpi-1.8.3-intel2015.0\n')
        ppfile.write('module load gcc/4.8.3\n')
        ppfile.write('\n')

        ppfile.write('source /projects/b1011/non-lsc/lscsoft-user-env.sh\n')
        ppfile.write('source /projects/b1011/ligo_project/lsc/o1_lalinference_20160402/etc/lscsoftrc\n')
        ppfile.write('\n')

	if args.plot_2d:
	        ppfile.write('cbcBayesPostProc.py --lalinfmcmc -i {} --event {} --outpath={} -d {}/PTMCMC.output.*.00 --dievidence  --skyres=.5 --deltaLogL {} --plot-2d\n'.format(args.inj,args.event,webdir,out_dir,target_hot_like))
	else:
		ppfile.write('cbcBayesPostProc.py --lalinfmcmc -i {} --event {} --outpath={} -d {}/PTMCMC.output.*.00 --dievidence  --skyres=.5 --deltaLogL {}\n'.format(args.inj,args.event,webdir,out_dir,target_hot_like))
	if args.ppall:
		for i in xrange(1,n_chains):
			ppfile.write('cbcBayesPostProc.py --lalinfmcmc -i {} --event {} --outpath={}{} -d {}/PTMCMC.output.*.0{} --dievidence  --skyres=.5 --deltaLogL {}\n'.format(args.inj,args.event,webdir,i,out_dir,i,target_hot_like))
	else:
		ppfile.close()

        system_call = 'chmod 755 {0}/pp.sh'.format(out_dir)
        os.system(system_call)

        ppsource = open('{0}/pp.sh'.format(args.homepath),"a+")
        ppsource.write('bash {0}/pp.sh\n'.format(out_dir))
        ppsource.close()

        ppmsub = open('{0}/ppmsub.sh'.format(args.homepath),"a+")
        ppmsub.write('msub {0}/pp.sh\n'.format(out_dir))
        ppmsub.close()

if args.compare:
    compFilePath = '/projects/b1011/spinning_runs/freezingparams_20160402/comp.sh'
    with open(compFilePath,'w') as compfile:
        compfile.write('#MSUB -A {}\n'.format(args.alloc))
        compfile.write('#MSUB -q {}\n'.format(args.queue))
        compfile.write('#MSUB -l walltime=00:08:00:00\n')
        compfile.write('#MSUB -l nodes=1:ppn=1\n')

        compfile.write('#MSUB -N fullcomp\n')

        # Give read permissions to screen output
        compfile.write('#MOAB -W umask=022\n')
        # Write stdout and stderr to the same file
        compfile.write('#MSUB -j oe\n')

        # Job working directory
        compfile.write('#MSUB -d {}\n'.format(args.homepath))

        compfile.write('\n')
        compfile.write('\n')
        compfile.write('ulimit -c unlimited\n')
        compfile.write('\n')

        compfile.write('module load python\n')
        compfile.write('module load mpi/openmpi-1.8.3-intel2015.0\n')
        compfile.write('module load gcc/4.8.3\n')
        compfile.write('\n')

        compfile.write('source /projects/b1011/non-lsc/lscsoft-user-env.sh\n')
        compfile.write('source /projects/b1011/ligo_project/lsc/o1_lalinference_20160402/etc/lscsoftrc\n')
        compfile.write('\n')


        compfile.write('cbcBayesCompPos.py -p https://ldas-jobs.ligo.caltech.edu/~leah.perri/freezingparams_20160402/{0}/none/post/posplots.html -n none -p https://ldas-jobs.ligo.caltech.edu/~leah.perri/freezingparams_20160402/{0}/skyloc/post/posplots.html -n skyloc -p https://ldas-jobs.ligo.caltech.edu/~leah.perri/freezingparams_20160402/{0}/skyloc_thetajn/post/posplots.html -n skyloc_thetajn -p https://ldas-jobs.ligo.caltech.edu/~leah.perri/freezingparams_20160402/{0}/skyloc_dist/post/posplots.html -n skyloc_dist -p https://ldas-jobs.ligo.caltech.edu/~leah.perri/freezingparams_20160402/{0}/skyloc_thetajn_dist/post/posplots.html -n skyloc_thetajn_dist -o ./comp -u scott.coughlin -x Balling1234 --ignore-missing-files --npixels-2d=200 --contour-dpi=200'.format(args.event))

        compmsub = open('{0}/compmsub.sh'.format(args.homepath),"a+")
        compmsub.write('msub {0}/comp.sh\n'.format(out_dir))
        compmsub.close()
# Create file of individual msub ../submit commands

outputLocation = ('{}/condor'.format(args.homepath))

if not os.path.isdir(outputLocation):
    os.makedirs(outputLocation)

dagFile = os.path.join(outputLocation,"msub.sh")

f = open(dagFile, "a+") # write mode

# The flag -a and -e will inform user when job has ended or been aborted

if args.emailyes:
        f.write('msub -m ae -M {} {}/submit\n'.format(args.email,out_dir))
else:
        f.write('msub {}/submit\n'.format(out_dir))
f.write('\n\n')
f.close()
# Make executable if not on quest
if not on_quest:
    st = os.stat(submitFilePath)
    os.chmod(submitFilePath, st.st_mode | stat.S_IEXEC)
