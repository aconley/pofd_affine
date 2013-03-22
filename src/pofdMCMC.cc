#include<sstream>
#include<unistd.h>

#include "../include/pofdMCMC.h"
#include "../include/specFile.h"
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
  \param[in] BURN_MULTIPLE This fraction of autocorrelation steps to add
                            before checking burn-in again
  \param[in] SCALEFAC Scale factor of Z distribution			    
 */
pofdMCMC::pofdMCMC(const std::string& INITFILE, const std::string& SPECFILE,
		   unsigned int NWALKERS, unsigned int NSAMPLES, 
		   unsigned int INIT_STEPS, unsigned int MIN_BURN, 
		   double BURN_MULTIPLE, double SCALEFAC) :
  affineEnsemble(NWALKERS, 1, NSAMPLES, INIT_STEPS, MIN_BURN,
		 BURN_MULTIPLE, SCALEFAC) {
  //Note that we set NPARAMS to a bogus value (1) above, then 
  // have to change it later once we know how many model params we have
  // All of this is done in initChains
}

bool pofdMCMC::initChainsMaster() {
  if (rank != 0)
    throw affineExcept("pofdMCMC","initChainsMaster",
		       "Should not be called except on master node",1);

  //Model initialization file
  ifile.readFile(initfile,true,true);
  unsigned int nknots = ifile.getNKnots();
  if (nknots == 0)
    throw affineExcept("pofdMCMC","initChainsMaster",
		       "No model knots read in",2);
  
  //Data/fit initialization file
  specFile spec_info(specfile);
  if (spec_info.datafiles.size() == 0)
      throw affineExcept("pofdMCMC","initChainsMaster",
			 "No data files read",3);


  //Number of parameters is number of knots + 1 for the
  // sigma -- although that may be fixed
  unsigned int npar = nknots+1;
  this->setNParams(npar);

  //Set up parameter names
  std::stringstream parname;
  for (unsigned int i = 0; i < nknots; ++i) {
    parname.str("Knot");
    parname << i;
    this->setParamName(i, parname.str());
  }
  this->setParamName(nknots, "SigmaMult");

  //For now just ignore where fixed
  for (unsigned int i = 0; i < nknots; ++i)
    if (ifile.isKnotFixed(i)) this->ignoreParam(i);
  //Ignore sigma if we aren't fitting it
  if (! spec_info.fit_sigma) this->ignoreParam(nknots);

  //Initialize likelihood information -- priors, data, etc.
  likeSet.setFFTSize(spec_info.fftsize);
  if (spec_info.edge_fix) likeSet.setEdgeFix(); else likeSet.unSetEdgeFix();
  likeSet.setNInterp(spec_info.ninterp);
  if (spec_info.bin_data) likeSet.setBinData(); else likeSet.unSetBinData();
  likeSet.setNBins(spec_info.nbins);
  if (spec_info.has_wisdom_file) likeSet.addWisdom(spec_info.wisdom_file);

  //Set priors 
  if (spec_info.has_cfirbprior)
    likeSet.setCFIRBPrior(spec_info.cfirbprior_mean,
			  spec_info.cfirbprior_stdev);
  if (spec_info.fit_sigma && spec_info.has_sigprior)
    likeSet.setSigmaPrior(spec_info.sigprior_stdev);

  //Read in data files
  if (spec_info.verbose || spec_info.ultraverbose)
      std::cout << "Reading in data files" << std::endl;
  likeSet.readDataFromFiles(spec_info.datafiles, spec_info.psffiles, 
			    spec_info.sigmas, spec_info.like_norm,
			    spec_info.ignore_mask, spec_info.mean_sub, 
			    spec_info.beam_histogram);


  //Now, copy that information over to slaves
  unsigned int nproc = MPI::COMM_WORLD.Get_size();

  unsigned int nnotinitialized;
  std::vector<bool> initialized(false,nproc);
  initialized[0] = true; //Master is initialized
  nnotinitialized = nproc-1;

  MPI::Status Info;
  int jnk;
  bool ismsg;
  while (nnotinitialized > 0) {
    //See if a message is available, wait for one if it isn't
    ismsg = MPI::COMM_WORLD.Iprobe(MPI::ANY_SOURCE, MPI::ANY_TAG,
				   Info);
    while (ismsg == false) {
      usleep(10000); //Sleep for 1/100th of a second
      ismsg = MPI::COMM_WORLD.Iprobe(MPI::ANY_SOURCE, MPI::ANY_TAG,
				     Info);
    }

    //We have a message -- process it
    MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,MPI::ANY_SOURCE,
			 MPI::ANY_TAG,Info);
    int this_tag = Info.Get_tag();
    int this_rank = Info.Get_source();
    if (this_tag == mcmc_affine::ERROR) {
      //An error was encountered -- the slave should have
      // notified all others of failure
      return false; 
    } else if (this_tag == pofd_mcmc::PMCMCSENDINIT) {
      MPI::COMM_WORLD.Send(&jnk, 1, MPI::INT, this_rank, 
			   pofd_mcmc::PMCMCSENDINGINIT);
      MPI::COMM_WORLD.Send(&npar, 1, MPI::UNSIGNED, this_rank, 
			   pofd_mcmc::PMCMCSENDNPAR);
      ifile.sendSelf(MPI::COMM_WORLD, this_rank);
      likeSet.sendSelf(MPI::COMM_WORLD, this_rank);
      
    } else if (this_tag == pofd_mcmc::PMCMCISREADY) {
      initialized[this_rank] = true;
    }
    nnotinitialized = std::count(initialized.begin(), initialized.end(),
				 false);
  }

  //We can free all the data storage, since master doesn't need
  // any of it
  likeSet.freeData();
  
  //Generate initial parameters for each walker
  //First, set up space to hold that first parameter
  chains.clear();
  chains.addChunk(1);
  chains.setSkipFirst();

  //Rather than generate actual likelihoods -- which would require some
  // complicated stuff -- we will assign them -infinite likelihood so 
  // that the first step is always taken
  paramSet p(npar); //Generated parameter

  //We may have to do trials to set up initial parameters, 
  const unsigned int maxtrials = 1000; //!< Maximum number of trials
  double trialval; //Test trial value

  //If we are fitting for sigma we need to generate a range of values
  // If there is a prior, we can use that as information for how to
  // distribute the multiplier, if not we need a default.
  double sm_stdev;
  if (spec_info.has_sigprior) 
    sm_stdev = spec_info.sigprior_stdev;
  else
    sm_stdev = 0.1;

  unsigned int nwalk = getNWalkers();
  for (unsigned int i = 0; i < nwalk; ++i) {
    //Fill in knot values.  Note that fixed parameters (sigma==0)
    // all get the same values for all walkers, which means all linear
    // combinations get the same value, so the parameter stays fixed
    ifile.generateRandomKnotValues(p);

    //Sigma multiplier value
    if (spec_info.fit_sigma) {
      unsigned int j;
      for (j = 0; j < maxtrials; ++j) {
	trialval = 1.0 + sm_stdev * rangen.gauss(); 
	if (trialval > 0) break;
      }
      if (j == maxtrials) {
	for (unsigned int k = 1; k < nproc; ++k) //Tell slaves to stop
	  MPI::COMM_WORLD.Send(&jnk, 1, MPI::INT, k, mcmc_affine::STOP);
	throw affineExcept("pofdMCMC", "initChainsMaster",
			   "Couldn't generate sigma multiplier value", 4);
      }
      p[nknots] = trialval;
    } else p[nknots] = 1.0;
    
    //Add to chain
    //It's possible that -infinity will not be available on some
    // architectures -- the c++ standard is oddly quiet about this.
    //It seems to work on g++-4.7 though
    chains.addNewStep(i, p, -std::numeric_limits<double>::infinity());
  }
  
  //That's it!
  return true;
}

bool pofdMCMC::initChainsSlave() {
  if (rank == 0)
    throw affineExcept("pofdMCMC","initChainsSlave",
		       "Should not be called on master node", 1);
  //Send message to master saying we are ready
  int jnk;
  MPI::COMM_WORLD.Send(&jnk, 1, MPI::INT, 0, pofd_mcmc::PMCMCSENDINIT);

  //And wait...
  MPI::Status Info;

  bool ismsg;
  //See if a message is available, wait for one if it isn't
  ismsg = MPI::COMM_WORLD.Iprobe(MPI::ANY_SOURCE, MPI::ANY_TAG, Info);
  while (ismsg == false) {
    usleep(10000); //Sleep for 1/100th of a second
    ismsg = MPI::COMM_WORLD.Iprobe(MPI::ANY_SOURCE, MPI::ANY_TAG,
				   Info);
  }

  //Got a message, process it
  MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,MPI::ANY_SOURCE,
		       MPI::ANY_TAG,Info);
  int this_tag = Info.Get_tag();
  int this_rank = Info.Get_source();

  if (this_rank != 0) {
    //Got a message from another slave.  This is not allowed!
    MPI::COMM_WORLD.Send(&jnk, 1, MPI::INT, 0, mcmc_affine::ERROR);
    return false;
  } else {
    if (this_tag == mcmc_affine::STOP)
      return false;
    else if (this_tag == pofd_mcmc::PMCMCSENDINGINIT) {
      unsigned int npar;
      MPI::COMM_WORLD.Recv(&npar, 1, MPI::UNSIGNED, 0, 
			   pofd_mcmc::PMCMCSENDNPAR);
      this->setNParams(npar);
      ifile.recieveCopy(MPI::COMM_WORLD, 0);
      likeSet.recieveCopy(MPI::COMM_WORLD, 0);

      MPI::COMM_WORLD.Send(&jnk, 1, MPI::INT, 0, pofd_mcmc::PMCMCISREADY);
    }
  }

  return true;
}

void pofdMCMC::initChains() {
  if (rank == 0) initChainsMaster(); else initChainsSlave();
}

/*!
  \param[in] p Parameters to evaluate model for
 */
double pofdMCMC::getLogLike(const paramSet& p) {
  return likeSet.getLogLike(p);
}
