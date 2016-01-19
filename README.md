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
                                       [--Niter NITER]
                                       [--fix-rightascension RIGHTASCENSION]

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
  --fix-rightascension RIGHTASCENSION
                        Fix RA

