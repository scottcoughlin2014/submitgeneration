# submitgeneration
Example of how to run this code

usage: lalinference_mcmc_submit_new.py [-h] [--alloc ALLOC] [--queue QUEUE]
                                       [--jobName JOBNAME] [--dir DIR]
                                       [--dep DEP] [--name NAME]
                                       [--walltime WALLTIME]
                                       [--cores-per-node CORES_PER_NODE]
                                       [--multiple-nodes] [--nPar NPAR]
                                       [--branch BRANCH] [--rc RC]
                                       [--sim-quest] [--quest-iwd QUEST_IWD]
                                       [--era ERA] [--psd PSD [PSD ...]]
                                       [--ifo IFO [IFO ...]] [--inj INJ]
                                       [--event EVENT] [--trigtime TRIGTIME]
                                       [--approx APPROX] [--amporder AMPORDER]
                                       [--flow FLOW] [--lowM1 LOWM1]
                                       [--lowM2 LOWM2] [--fhigh FHIGH]
                                       [--srate SRATE] [--seglen SEGLEN]
                                       [--psdlength PSDLENGTH]
                                       [--psdstart PSDSTART]
                                       [--temp-ladder-top-down] [--noisy]
                                       [--distance-max DISTANCE_MAX]
                                       [--no-malmquist] [--no-margtimephi]
                                       [--trigSNR TRIGSNR] [--tempMin TEMPMIN]
                                       [--tempMax TEMPMAX] [--Neff NEFF]
                                       [--Niter NITER] [--fix-rightascension]
                                       [--fix-declination] [--fix-distance]
                                       [--fix-costheta_jn]

Generate a submit file for lalinference_mcmc on grail.

optional arguments:
  -h, --help            show this help message and exit

MSUB:
  --alloc ALLOC         Allocation to charge SUs to (default=b1011).
  --queue QUEUE         Queue for job (default=ligo).
  --jobName JOBNAME     Name of job, used for output file names and queue
                        listing.
  --dir DIR             Directory where submit script is written and executed
                        from (default=current directory).
  --dep DEP             Wait for dependent job ID to complete (default=None).
  --name NAME           Name of submit file (default=submit).
  --walltime WALLTIME   Walltime for job (default=2:00:00:00).
  --cores-per-node CORES_PER_NODE
                        Number of cores per node (default=16).
  --multiple-nodes      If nChains > 16 then use more than one node.
  --nPar NPAR           Number of dimensions for MCMC. Defaults for common
                        templates are set, assuming no PSD fitting.

env:
  --branch BRANCH       Branchname to use, assuming
                        /projects/p20251/USER/lsc/BRANCHNAME/etc/lscsoftrc
                        exists (default=master).
  --rc RC               Specify direct path to rc files to be sourced (e.g.
                        lscsoftrc). /projects/b1011/non-lsc/lscsoft-user-
                        env.sh added by default.
  --sim-quest           Act as if on Quest. Useful for setting up submit files
                        on local computer for uploading to Quest
  --quest-iwd QUEST_IWD
                        Working directory for Quest submit script.
lalinference_mcmc:
  --era ERA             Era ("initial" or "advanced") of detector PSD for SNR
                        calculation. If no cache arguments given, this will
                        add the appropriate analytical PSD arguments to the
                        submit file.
  --psd PSD [PSD ...]   Pre-computed PSD(s), either as a single xml, or list
                        of ascii files, one for each IFO in the order the IFOs
                        were specified. XMLs are converted using
                        lalinference_pipe_utils.
  --ifo IFO [IFO ...]   IFOs for the analysis.
  --inj INJ             Injection XML file.
  --event EVENT         Event number in XML to inject.
  --trigtime TRIGTIME   Trigger time of event. Automatically set when
                        injecting.
  --approx APPROX       Specify a template approximant (default SpinTaylorT4).
  --amporder AMPORDER   Specify amplitude order of template.
  --flow FLOW           Lower frequency bound for all detectors (default=40).
  --lowM1 LOWM1         Lowest reasonable M1 of the binary your are searching
                        for. Used to determine appropraite seglength.
  --lowM2 LOWM2         Lowest reasonable M2 of the binary your are searching
                        for. Used to determine appropraite seglength.
  --fhigh FHIGH         Upper frequency bound for all detectors, given as a
                        fraction of the injection's ISCO frequency.
  --srate SRATE         Sampling rate of the waveform. If not provided and an
                        injection is peformed, it is set to be sufficient for
                        the signal being injected. If no injection, it
                        defaults to a sufficient value for a 1.4-1.4 binary
                        coalescence (expensive!).
  --seglen SEGLEN       Length of data segment used for likelihood
                        compuatation. Same default behavior as "--srate".
  --psdlength PSDLENGTH
                        Length of data segment to use for PSD estimation.
                        Defaults to 32*seglen.
  --psdstart PSDSTART   GPS time to start PSD calculation. Defaults to
                        trigtime - psdlength - seglen
  --temp-ladder-top-down
                        Build the temperature ladder from the bottom up, using
                        an analytic prescription for the spacing that should
                        ensure communication between chains. Sets the number
                        of cores so that the hottest temperature should be
                        sampling the prior.
  --noisy               Use a non-zero noise realization.
  --distance-max DISTANCE_MAX
                        Hard outer prior boundary on distance (default=1 Gpc
                        iLIGO, 2 Gpc aLIGO).
  --no-malmquist        Do not use the Malmquist prior, which by default
                        approximates the selection effects imposed by the
                        detection processes as a cut in the SNR in the second
                        loudest detector.
  --no-margtimephi      Do not use the time and phase marginalized likelihood
                        function.
  --trigSNR TRIGSNR     SNR of the trigger (calculated automatically if
                        injecting).
  --tempMin TEMPMIN     Temperature of coldest chain (default=1.0).
  --tempMax TEMPMAX     Temperature of hotest chain. Determined automatically
                        if injecting, or trigSNR is given.
  --Neff NEFF           Requested number of independent samples.
  --Niter NITER         Maximum number of MCMC iterations to allow.
  --fix-rightascension  Fix RA
  --fix-declination     Fix DEC
  --fix-distance        Fix Distance
  --fix-costheta_jn     Fix costheta_jn

Help from LAL_MCMC

 ========== LALInference_MCMC ==========
    ----------------------------------------------
    --- Data Parameters --------------------------
    ----------------------------------------------
    --ifo IFO1 [--ifo IFO2 ...] IFOs can be H1,L1,V1
    --IFO1-cache cache1         Cache files 
    [--IFO2-cache2 cache2 ...]      lal PSDs: LAL{Ad}LIGO, LALVirgo
                                    lalsimuation PSDs: LALSim{Ad}LIGO, LALSim{Ad}Virgo
    --psdstart GPStime          GPS start time of PSD estimation data
    --psdlength length          Length of PSD estimation data in seconds
    --seglen length             Length of segments for PSD estimation and analysis in seconds
    (--dont-dump-extras)        If given, won't save PSD and SNR files
    (--trigtime GPStime)        GPS time of the trigger to analyse
                                    (optional when using --margtime or --margtimephi)
    (--segment-start)           GPS time of the start of the segment
                                     (optional with --trigtime,
                                      default: seglen-2 s before --trigtime)
    (--srate rate)              Downsample data to rate in Hz (4096.0,)
    (--padding PAD [sec]        Override default 0.4 seconds padding
    (--injectionsrate rate)     Downsample injection signal to rate in Hz (--srate)
    (--IFO1-flow freq1          Specify lower frequency cutoff for overlap integral (40.0)
     [--IFO2-flow freq2 ...])
    (--IFO1-fhigh freq1         Specify higher frequency cutoff for overlap integral (Nyquist
     [--IFO2-fhigh freq2 ...])      freq 0.5*srate)
    (--IFO1-channel chan1       Specify channel names when reading cache files
     [--IFO2-channel chan2 ...])
    (--IFO1-psd psd1-ascii.txt        Read in PSD from ascii file. This is not equivalent 
     [--IFO2-psd psd2-ascii.txt ...])     to using --IFO1-cache interp:file.txt since the former
                                          won't use the ascii psd to generate fake noise. 
    (--dataseed number)         Specify random seed to use when generating data
    (--lalinspiralinjection)    Enables injections via the LALInspiral package
    (--inj-fref)                Reference frequency of parameters in injection XML (default 100Hz)
    (--inj-lambda1)             value of lambda1 to be injected, LALSimulation only (0)
    (--inj-lambda2)             value of lambda2 to be injected, LALSimulation only (0)
    (--inj-lambdaT              value of lambdaT to be injected (0)
    (--inj-dlambdaT             value of dlambdaT to be injected (0)
    (--inj-spinOrder PNorder)   Specify twice the injection PN order (e.g. 5 <==> 2.5PN)
                                    of spin effects effects to use, only for LALSimulation
                                    (default: -1 <==> Use all spin effects).
    (--inj-tidalOrder PNorder)  Specify twice the injection PN order (e.g. 10 <==> 5PN)
                                    of tidal effects to use, only for LALSimulation
                                    (default: -1 <==> Use all tidal effects).
    (--inj-spin-frame FRAME     Specify injection spin frame: choice of total-j, orbital-l, view.
                                    (Default = OrbitalL).
    (--0noise)                  Sets the noise realisation to be identically zero
                                    (for the fake caches above only)
    

------------------------------------------------------------------------------------------------------------------
--- Calibration Errors Handling Arguments ------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------
(--AddCalibrationErrors) Adds calibration errors into the f domain datastream (that includes both noise and signal)
(--RandomCE) Add a random realization of phase and amplitude CE, using the S6/VSR2-3 error budget as an indication of the 1-sigma errors
(--ConstantCE) Assumes calibration errors are constant over the bandwidth (requires --IFO-constant_calamp and --IFO-constant_calpha for each ifo)
(--IFO-constant_calamp Constant amplitude CE for instrument IFO. 0.0 means no error, 0.1 means 10 percent
(--IFO-constant_calpha Constant phase CE for instrument IFO. 0.0 means no error, 5 means a  5 degree shift 
(--RandomLinearCE ) Assumes CE are given by a contant plateau plus a random jittering of a few percent.
		 After a given frequency f CE increase linearly with a given slope (requires --IFO-calamp_plateau, --IFO-calamp_knee, --IFO-calamp_slope AND similar with calamp<->calpha. Phase Errors in degree, Amp errors are relative (0.05=5%) 
(--IFO-calamp_plateau, --IFO-calamp_knee, --IFO-calamp_slope) Add on the i-th IFO's stream amplitude errors on the form (IFOi_c + jitter) for f<IFOi_f and (IFOi_c-f)*IFOi_slope for f>IFOi_f
(--IFO-calpha_plateau, --IFO-calpha_knee, --IFO-calpha_slope) Add on the i-th IFO's stream phase errors on the form (IFOi_c + jitter) for f<IFOi_f and (IFOi_c-f)*IFOi_slope for f>IFOi_f
 * Constant Calibration Model 
  (--MarginalizeConstantCalAmp ) If given, will add a constant value of Amplitude CE per each IFO on the top of the CBC parameters.
  (--MarginalizeConstantCalPha ) If given, will add a constant value of Phase CE per each IFO on the top of the CBC parameters.
  (--constcal_ampsigma ) If given, will use gaussian prior on the constant amplitude error with this sigma (e.g. 0.05=5%) .
  (--constcal_phasigma ) If given, will use gaussian prior on the constant phase error with this sigma (e.g. 5=5degs).
 * Spline Calibration Model 
  (--enable-spline-calibration)            Enable cubic-spline calibration error model.
  (--spcal-nodes N)           Set the number of spline nodes per detector (default 5)
  (--spcal-amp-uncertainty X) Set the prior on relative amplitude uncertainty (default 0.1)
  (--spcal-phase-uncertainty X) Set the prior on phase uncertanity in degrees (default 5)


    ----------------------------------------------
    --- MCMC Algorithm Parameters ----------------
    ----------------------------------------------
    (--nsteps n)        Maximum number of steps to take (1e7)
    (--neff N)          Number of independent samples to collect (nsteps)
    (--skip n)          Number of steps between writing samples to file (100)
    (--adapt-tau)       Adaptation decay power, results in adapt length of 10^tau (5)
    (--no-adapt)        Do not adapt run
    (--randomseed seed) Random seed of sampling distribution (random)

    ----------------------------------------------
    --- Parallel Tempering Algorithm Parameters --
    ----------------------------------------------
    (--temp-skip N)     Number of steps between temperature swap proposals (100)
    (--tempKill N)      Iteration number to stop temperature swapping (Niter)
    (--ntemps N)         Number of temperature chains in ladder (as many as needed)
    (--temp-min T)      Lowest temperature for parallel tempering (1.0)
    (--temp-max T)      Highest temperature for parallel tempering (50.0)
    (--anneal)          Anneal hot temperature linearly to T=1.0
    (--annealStart N)   Iteration number to start annealing (5*10^5)
    (--annealLength N)  Number of iterations to anneal all chains to T=1.0 (1*10^5)
    
    ----------------------------------------------
    --- Noise Model ------------------------------
    ----------------------------------------------
    (--psd-fit)         Run with PSD fitting
    (--psdNblock)       Number of noise parameters per IFO channel (8)
    (--psdFlatPrior)    Use flat prior on psd parameters (Gaussian)
    (--glitch-fit)       Run with glitch fitting
    (--glitchNmax)      Maximum number of glitch basis functions per IFO (20)
    
    ----------------------------------------------
    --- Output -----------------------------------
    ----------------------------------------------
    (--data-dump)       Output waveforms to file
    (--adapt-verbose)   Output parameter jump sizes and acceptance rate stats to file
    (--temp-verbose)    Output temperature swapping stats to file
    (--prop-verbose)    Output proposal stats to file
    (--prop-track)      Output proposal parameters
    (--outfile file)    Write output files <file>.<chain_number> 
                            (PTMCMC.output.<random_seed>.<mpi_thread>)
    
    ----------------------------------------------
    --- Injection Arguments ----------------------
    ----------------------------------------------
    (--inj injections.xml) Injection XML file to use
    (--event N)            Event number from Injection XML file to use

    ----------------------------------------------
    --- Template Arguments -----------------------
    ----------------------------------------------
    (--use-eta)            Jump in symmetric mass ratio eta, instead of q=m1/m2 (m1>m2)
    (--approx)             Specify a template approximant and phase order to use
                         (default TaylorF2threePointFivePN). Available approximants:
                         default modeldomain="time": GeneratePPN, TaylorT1, TaylorT2,
                                                       TaylorT3, TaylorT4, EOB, EOBNR,
                                                       EOBNRv2, EOBNRv2HM, SEOBNRv1,
                                                       SpinTaylor, SpinQuadTaylor, 
                                                       SpinTaylorFrameless, SpinTaylorT4,
                                                       PhenSpinTaylorRD, NumRel.
                         default modeldomain="frequency": TaylorF1, TaylorF2, TaylorF2RedSpin,
                                                       TaylorF2RedSpinTidal, IMRPhenomA,
                                                       IMRPhenomB, IMRPhenomP.
    (--amporder PNorder)            Specify a PN order in amplitude to use (defaults: LALSimulation: max available; LALInspiral: newtownian).
    (--fref f_ref)                  Specify a reference frequency at which parameters are defined (default 100).
    (--use-tidal)                   Enables tidal corrections, only with LALSimulation.
    (--use-tidalT)                  Enables reparmeterized tidal corrections, only with LALSimulation.
    (--spinOrder PNorder)           Specify twice the PN order (e.g. 5 <==> 2.5PN) of spin effects to use, only for LALSimulation (default: -1 <==> Use all spin effects).
    (--tidalOrder PNorder)          Specify twice the PN order (e.g. 10 <==> 5PN) of tidal effects to use, only for LALSimulation (default: -1 <==> Use all tidal effects).
    (--modeldomain)                 domain the waveform template will be computed in ("time" or "frequency"). If not given will use LALSim to decide
    (--spinAligned or --aligned-spin)  template will assume spins aligned with the orbital angular momentum.
    (--singleSpin)                  template will assume only the spin of the most massive binary component exists.
    (--noSpin, --disable-spin)      template will assume no spins (giving this will void spinOrder!=0) 
    (--no-detector-frame)              model will NOT use detector-centred coordinates and instead RA,dec

    ----------------------------------------------
    --- Starting Parameters ----------------------
    ----------------------------------------------
    You can generally have MCMC chains to start from a given parameter value by using --parname VALUE. Names currently known to the code are:
     time                         Waveform time (overrides random about trigtime).
     chirpmass                    Chirpmass
     eta                          Symmetric massratio (needs --use-eta)
     q                            Asymmetric massratio (a.k.a. q=m2/m1 with m1>m2)
     phase                        Coalescence phase.
     costheta_jn                  Cosine of angle between J and line of sight [rads]
     logdistance                  Log Distance (requires --use-logdistance)
     rightascension               Rightascensions
     declination                  Declination.
     polarisation                 Polarisation angle.
    * Spin Parameters:
     a_spin1                      Spin1 magnitude
     a_spin2                      Spin2 magnitude
     tilt_spin1                   Angle between spin1 and orbital angular momentum
     tilt_spin2                   Angle between spin2 and orbital angular momentum 
     phi_12                       Difference between spins' azimuthal angles 
     phi_jl                       Difference between total and orbital angular momentum azimuthal angles
    * Equation of State parameters (requires --use-tidal or --use-tidalT):
     lambda1                      lambda1.
     lambda2                      lambda2.
     lambdaT                      lambdaT.
     dLambdaT                     dLambdaT.

    ----------------------------------------------
    --- Prior Ranges -----------------------------
    ----------------------------------------------
    You can generally use --paramname-min MIN --paramname-max MAX to set the prior range for the parameter paramname
    The names known to the code are listed below.
    Component masses, total mass and time have dedicated options listed here:

    (--trigtime time)                       Center of the prior for the time variable.
    (--comp-min min)                        Minimum component mass (1.0).
    (--comp-max max)                        Maximum component mass (30.0).
    (--mass1-min min, --mass1-max max)      Min and max for mass1 (default: same as comp-min,comp-max, will over-ride these.
    (--mass2-min min, --mass2-max max)      Min and max for mass2 (default: same as comp-min,comp-max, will over-ride these.
    (--mtotal-min min)                      Minimum total mass (2.0).
    (--mtotal-max max)                      Maximum total mass (35.0).
    (--dt time)                             Width of time prior, centred around trigger (0.2s).

    (--varyFlow, --flowMin, --flowMax)       Allow the lower frequency bound of integration to vary in given range.
    (--pinparams)                            List of parameters to set to injected values [mchirp,asym_massratio,etc].
    ----------------------------------------------
    --- Fix Parameters ---------------------------
    ----------------------------------------------
    You can generally fix a parameter to be fixed to a given values by using --fix-paramname VALUE
    where the known names have been listed above.

    ----------------------------------------------
    --- Prior Arguments --------------------------
    ----------------------------------------------
    (--malmquistprior)               Impose selection effects on the prior (False)
    (--malmquist-loudest-snr)        Threshold SNR in the loudest detector (0.0)
    (--malmquist-second-loudest-snr) Threshold SNR in the second loudest detector (5.0)
    (--malmquist-network-snr)        Threshold network SNR (0.0)
    (--analyticnullprior)            Use analytic null prior
    (--nullprior)                    Use null prior in the sampled parameters

    
    ----------------------------------------------
    --- Likelihood Arguments ---------------------
    ----------------------------------------------
    (--zeroLogLike)                  Use flat, null likelihood
    (--studentTLikelihood)           Use the Student-T Likelihood that marginalizes over noise
    (--correlatedGaussianLikelihood) Use analytic, correlated Gaussian for Likelihood
    (--bimodalGaussianLikelihood)    Use analytic, bimodal correlated Gaussian for Likelihood
    (--rosenbrockLikelihood)         Use analytic, Rosenbrock banana for Likelihood
    (--noiseonly)                    Using noise-only likelihood
    (--margphi)                      Using marginalised phase likelihood
    (--margtime)                     Using marginalised time likelihood
    (--margtimephi)                  Using marginalised in time and phase likelihood

