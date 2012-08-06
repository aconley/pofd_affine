#include<pofdMCMC.h>
#include<affineExcept.h>

/*!
  \param[in] INITFILE Name of model initialization file; see initFileKnots
  \param[in] SPECFILE Name of fit specification file; see specFile
  \param[in] NWALKERS Number of walkers
  \param[in] NPARAMS  Number of parameters
  \param[in] NSAMPLES Number of samples (across all walkers) to do after burn;
                      realized to the closest larger multiple of nwalkers
  \param[in] INIT_STEPS Number of initialization steps, which are thrown away
                         even before starting burn-in process
  \param[in] MIN_BURN Minimum number of steps before burn-in
  \param[in] BURN_MULTIPLE This fraction of autocorrelation steps to add
                            before checking burn-in again
  \param[in] SCALEFAC Scale factor of Z distribution			    
 */
pofdMCMC(const std::string& INITFILE, const std::string& SPECFILE,
	 unsigned int NWALKERS, unsigned int NSAMPLES, 
	 unsigned int INIT_STEPS, unsigned int MIN_BURN, 
	 double BURN_MULTIPLE, double SCALEFAC) :
  affineEnsemble(NWALKERS, 1, NSAMPLES, INIT_STEPS, MIN_BURN,
		 BURN_MULTIPLE, SCALEFAC) {
  //Note that we set NPARAMS to a bogus value (1) above, then 
  // have to change it later once we know how many model params we have
  // All of this is done in initChains
}

void pofdMCMC::initChainsMaster() {
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

  //Use that to initialize likelihood 
  likeSet.setFFTSize( spec_info.fftsize );
  if (spec_info.edge_fix) likeSet.setEdgeFix(); else likeSet.unSetEdgeFix();
  likeSet.setNInterp( spec_info.ninterp );
  if (spec_info.bin_data) likeSet.setBinData(); else likeSet.unSetBinData();
  likeSet.setNBins( spec_info.nbins );
  if (spec_info.has_wisdom_file) likeSet.addWisdom(spec_info.wisdom_file);

  //Set priors 
  if (spec_info.has_cfirbprior)
    likeSet.setCFIRBPrior( spec_info.cfirbprior_mean,
			   spec_info.cfirbprior_stdev );
  if (spec_info.fit_sigma && spec_info.has_sigprior)
    likeSet.setSigmaPrior(spec_info.sigprior_stdev);

  //Read in data files
  if (spec_info.verbose || spec_info.ultraverbose)
      std::cout << "Reading in data files" << std::endl;
  likeSet.readDataFromFiles( spec_info.datafiles, spec_info.psffiles, 
			     spec_info.sigmas, spec_info.like_norm,
			     spec_info.ignore_mask, spec_info.mean_sub, 
			     spec_info.beam_histogram );

  //Now, copy that information over to slaves
}

void pofdMCMC::initChainsSlave() {
  //Send message to master saying we are ready

  //Copy info
}

void pofdMCMC::initChains() {
  if (rank == 0) initChainsMaster(); else initChainsSlave();
}
