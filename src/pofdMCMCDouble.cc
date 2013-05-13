#include<sstream>
#include<unistd.h>

#include "../include/pofdMCMCDouble.h"
#include "../include/affineExcept.h"

/*!
  \param[in] INITFILE Name of model initialization file; see initFileKnots
  \param[in] SPECFILE Name of fit specification file; see specFile
  \param[in] NWALKERS Number of walkers
  \param[in] NSAMPLES Number of samples (across all walkers) to do after burn;
                      realized to the closest larger multiple of nwalkers
  \param[in] INIT_STEPS Number of initialization steps, which are thrown away
                         even before starting burn-in process
  \param[in] MIN_BURN Minimum number of steps before burn-in
  \param[in] FIXED_BURN Do a fixed burn in of length MIN_BURN, not using
                         autocorrelation length to decide when burn in is 
                         finished
  \param[in] BURN_MULTIPLE This fraction of autocorrelation steps to add
                            before checking burn-in again
  \param[in] SCALEFAC Scale factor of Z distribution			    
 */
pofdMCMCDouble::pofdMCMCDouble(const std::string& INITFILE, 
			       const std::string& SPECFILE,
			       unsigned int NWALKERS, unsigned int NSAMPLES, 
			       unsigned int INIT_STEPS, unsigned int MIN_BURN, 
			       bool FIXED_BURN, float BURN_MULTIPLE, 
			       float SCALEFAC) :
  affineEnsemble(NWALKERS, 1, NSAMPLES, INIT_STEPS, MIN_BURN,
		 FIXED_BURN, BURN_MULTIPLE, SCALEFAC),
  initfile(INITFILE), specfile(SPECFILE) {
  //Note that we set NPARAMS to a bogus value (1) above, then 
  // have to change it later once we know how many model params we have
  // All of this is done in initChains
}

bool pofdMCMCDouble::initChainsMaster() {
  if (rank != 0)
    throw affineExcept("pofdMCMCDouble","initChainsMaster",
		       "Should not be called except on master node",1);

  //Model initialization file
  ifile.readFile(initfile,true,true);
  unsigned int ntot = ifile.getNTot();
  if (ntot == 0)
    throw affineExcept("pofdMCMCDouble","initChainsMaster",
		       "No model parameters read in",2);
  
  //Data/fit initialization file
  spec_info.readFile(specfile);
  if (spec_info.datafiles1.size() == 0)
      throw affineExcept("pofdMCMCDouble","initChainsMaster",
			 "No data files read",3);

  //Number of parameters is number of knots (ntot) + 2 for the
  // sigma in each band -- although some may be fixed
  unsigned int npar = ntot + 2;
  setNParams(npar);

  //Set up parameter names
  std::stringstream parname;
  unsigned int nknots = ifile.getNKnots();
  for (unsigned int i = 0; i < nknots; ++i) {
    parname.str("Knot");
    parname << i;
    setParamName(i, parname.str());
  }
  unsigned int nsigmas = ifile.getNSigmas();
  for (unsigned int i = 0; i < nsigmas; ++i) {
    parname.str("SigmaKnot");
    parname << i;
    setParamName(i + nknots, parname.str());
  }
  unsigned int noffsets = ifile.getNOffsets();
  for (unsigned int i = 0; i < noffsets; ++i) {
    parname.str("OffsetKnot");
    parname << i;
    setParamName(i + nknots + nsigmas, parname.str());
  }
  setParamName(ntot, "SigmaMult1");
  setParamName(ntot + 1, "SigmaMult2");

  // Set which parameters are fixed
  for (unsigned int i = 0; i < ntot; ++i)
    if (ifile.isKnotFixed(i)) this->fixParam(i);
  if (!spec_info.fit_sigma1) this->fixParam(ntot);
  if (!spec_info.fit_sigma2) this->fixParam(ntot + 1);

  //Initialize likelihood information -- priors, data, etc.
  likeSet.setPositions(ifile);
  likeSet.setFFTSize(spec_info.fftsize);
  if (spec_info.edge_fix) likeSet.setEdgeFix(); else likeSet.unSetEdgeFix();
  if (spec_info.edge_set) {
    likeSet.setEdgeInteg();
    likeSet.setNEdge(spec_info.nedge);
  }
  if (spec_info.bin_data) {
    likeSet.setBinData(); 
    likeSet.setNBins(spec_info.nbins);
  } else likeSet.unSetBinData();
  if (spec_info.has_wisdom_file) likeSet.addWisdom(spec_info.wisdom_file);

  //Set priors 
  if (spec_info.has_cfirbprior1)
    likeSet.setCFIRBPrior1(spec_info.cfirbprior_mean1,
			   spec_info.cfirbprior_stdev1);
  if (spec_info.has_cfirbprior2)
    likeSet.setCFIRBPrior2(spec_info.cfirbprior_mean2,
			   spec_info.cfirbprior_stdev2);
  if (spec_info.fit_sigma1 && spec_info.has_sigprior1)
    likeSet.setSigmaPrior1(spec_info.sigprior_stdev1);
  if (spec_info.fit_sigma2 && spec_info.has_sigprior2)
    likeSet.setSigmaPrior2(spec_info.sigprior_stdev2);

  //Read in data files
  if (spec_info.verbose || spec_info.ultraverbose)
      std::cout << "Reading in data files" << std::endl;
  likeSet.readDataFromFiles(spec_info.datafiles1, spec_info.datafiles2,
			    spec_info.psffiles1, spec_info.psffiles2, 
			    spec_info.sigmas1, spec_info.sigmas2,
			    spec_info.like_norm, spec_info.ignore_mask, 
			    spec_info.mean_sub, spec_info.beam_histogram);

  // Verbosity
  if (spec_info.verbose) setVerbose();
  if (spec_info.ultraverbose) setUltraVerbose();

  //Now, copy that information over to slaves
  int nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if (ultraverbose)
    std::cout << "Initializing " << nproc - 1 << " slave nodes from master"
	      << std::endl;

  int nnotinitialized;
  std::vector<bool> initialized(nproc, false);
  initialized[0] = true; //Master is initialized
  nnotinitialized = std::count(initialized.begin(), initialized.end(),
			       false);

  MPI_Status Info;
  int jnk, ismsg;
  while (nnotinitialized > 0) {
    //See if a message is available, wait for one if it isn't
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &ismsg, &Info);
    while (ismsg == 0) {
      usleep(10000); //Sleep for 1/100th of a second
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &ismsg, &Info);
    }

    MPI_Recv(&jnk, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, 
	     MPI_COMM_WORLD, &Info);
    int this_tag = Info.MPI_TAG;
    int this_rank = Info.MPI_SOURCE;
    if (this_tag == mcmc_affine::ERROR) {
      //An error was encountered -- the slave should have
      // notified all others of failure
      return false; 
    } else if (this_tag == pofd_mcmc::PMCMCSENDINIT) {
      MPI_Send(&jnk, 1, MPI_INT, this_rank, pofd_mcmc::PMCMCSENDINGINIT,
	       MPI_COMM_WORLD);
      MPI_Send(&npar, 1, MPI_UNSIGNED, this_rank, pofd_mcmc::PMCMCSENDNPAR,
	       MPI_COMM_WORLD);
      ifile.sendSelf(MPI_COMM_WORLD, this_rank);
      likeSet.sendSelf(MPI_COMM_WORLD, this_rank);
      
    } else if (this_tag == pofd_mcmc::PMCMCISREADY) {
      if (ultraverbose)
	std::cout << " Slave node " << this_rank << " initialized"
		  << std::endl;

      initialized[this_rank] = true;
    }
    nnotinitialized = std::count(initialized.begin(), initialized.end(),
				 false);
  }
  if (ultraverbose)
    std::cout << "All slave nodes initialized" << std::endl;


  //We can free all the data storage, since master doesn't need
  // any of it
  likeSet.freeData();

  // Generate initial positions
  if (ultraverbose)
    std::cout << "Setting up initial parameters" << std::endl;
  paramSet p(npar); //Generated parameter
  ifile.getParams(p); //Get central values from initialization file
  generateInitialPosition(p);

  is_init = true;

  //That's it!
  if (ultraverbose)
    std::cout << "Initialization completed" << std::endl;
  return true;
}

void pofdMCMCDouble::generateInitialPosition(const paramSet& p) {
  
  //Generate initial parameters for each walker
  if (rank != 0) return;

  unsigned int ntot = ifile.getNTot();
  if (ntot == 0)
    throw affineExcept("pofdMCMCDouble", "generateInitialPosition",
		       "No model info read in", 1);

  unsigned int npar = getNParams();
  if (p.getNParams() < npar)
    throw affineExcept("pofdMCMCDouble", "generateInitialPosition",
		       "Wrong number of params in input", 2);

  chains.clear();
  chains.addChunk(1);

  //Rather than generate actual likelihoods -- which would require some
  // complicated stuff -- we will assign them -infinite likelihood so 
  // that the first step is always taken
  paramSet pnew(npar); //Generated parameter

  //We may have to do trials to set up initial parameters, 
  const unsigned int maxtrials = 1000; //!< Maximum number of trials
  double trialval; //Test trial value

  //If we are fitting for sigma we need to generate a range of values
  // If there is a prior, we can use that as information for how to
  // distribute the multiplier, if not we need a default.
  double sm_stdev1, sm_stdev2;
  if (spec_info.has_sigprior1) 
    sm_stdev1 = spec_info.sigprior_stdev1;
  else
    sm_stdev1 = 0.1;
  if (spec_info.has_sigprior2) 
    sm_stdev2 = spec_info.sigprior_stdev2;
  else
    sm_stdev2 = 0.1;

  unsigned int nwalk = getNWalkers();
  for (unsigned int i = 0; i < nwalk; ++i) {
    //Fill in knot values.  Note that fixed parameters (sigma==0)
    // all get the same values for all walkers, which means all linear
    // combinations get the same value, so the parameter stays fixed
    ifile.generateRandomKnotValues(pnew, p);

    //Sigma multiplier value, one copy for each band
    if (spec_info.fit_sigma1) {
      unsigned int j;
      for (j = 0; j < maxtrials; ++j) {
	trialval = 1.0 + sm_stdev1 * rangen.gauss(); 
	if (trialval > 0) break;
      }
      if (j == maxtrials)
	throw affineExcept("pofdMCMCDouble", "initChainsMaster",
			   "Couldn't generate sigma multiplier1 value", 3);
      pnew[ntot] = trialval;
    } else pnew[ntot] = 1.0;

    if (spec_info.fit_sigma2) {
      unsigned int j;
      for (j = 0; j < maxtrials; ++j) {
	trialval = 1.0 + sm_stdev2 * rangen.gauss(); 
	if (trialval > 0) break;
      }
      if (j == maxtrials)
	throw affineExcept("pofdMCMCDouble", "initChainsMaster",
			   "Couldn't generate sigma multiplier2 value", 4);
      pnew[ntot + 1] = trialval;
    } else pnew[ntot + 1] = 1.0;
    
    //Add to chain
    //It's possible that -infinity will not be available on some
    // architectures -- the c++ standard is oddly quiet about this.
    //It seems to work on g++-4.7 though
    chains.addNewStep(i, pnew, -std::numeric_limits<double>::infinity());
    naccept[i] = 0;
  }
  chains.setSkipFirst();
}

bool pofdMCMCDouble::initChainsSlave() {
  if (rank == 0)
    throw affineExcept("pofdMCMCDouble","initChainsSlave",
		       "Should not be called on master node", 1);
  //Send message to master saying we are ready
  int jnk;
  MPI_Send(&jnk, 1, MPI_INT, 0, pofd_mcmc::PMCMCSENDINIT, MPI_COMM_WORLD);

  //And wait...
  MPI_Status Info;
  int ismsg;
  MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &ismsg, &Info);
  while (ismsg == 0) {
    usleep(10000); //Sleep for 1/100th of a second
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &ismsg, &Info);
  }

  //Got a message, process it
  MPI_Recv(&jnk, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, 
	   &Info);
  int this_tag = Info.MPI_TAG;
  int this_rank = Info.MPI_SOURCE;

  if (this_rank != 0) {
    //Got a message from another slave.  This is not allowed!
    MPI_Send(&jnk, 1, MPI_INT, 0, mcmc_affine::ERROR, MPI_COMM_WORLD);
    return false;
  } else {
    if (this_tag == mcmc_affine::STOP)
      return false;
    else if (this_tag == pofd_mcmc::PMCMCSENDINGINIT) {
      unsigned int npar;
      MPI_Recv(&npar, 1, MPI_UNSIGNED, 0, pofd_mcmc::PMCMCSENDNPAR,
	       MPI_COMM_WORLD, &Info);
      this->setNParams(npar);
      ifile.recieveCopy(MPI_COMM_WORLD, 0);
      likeSet.recieveCopy(MPI_COMM_WORLD, 0);

      MPI_Send(&jnk, 1, MPI_INT, 0, pofd_mcmc::PMCMCISREADY,
	       MPI_COMM_WORLD);
    }
  }
  
  is_init = true;
  return true;
}

void pofdMCMCDouble::initChains() {
  if (rank == 0) initChainsMaster(); else initChainsSlave();
  MPI_Barrier(MPI_COMM_WORLD);
}

/*!
  \param[in] p Parameters to evaluate model for
 */
double pofdMCMCDouble::getLogLike(const paramSet& p) {
  if (!is_init)
    throw affineExcept("pofdMCMC", "getLogLike",
		       "Called on unitialized object", 1);
  return likeSet.getLogLike(p);
}
