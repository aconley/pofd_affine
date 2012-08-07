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


  //Number of parameters is number of knots + 1 for the
  // sigma -- although that may be fixed
  unsigned int nknots = ifile.getNKnots();
  nparams = nknots + 1;

  //We need to resize a few things now
  pstep.setNParams(nparams);
  chains.setNParams(nparams);

  //Set up parameter names
  has_name.resize(nparams);
  for (unsigned int i = 0; i < nparams; ++i) has_name[i] = false;
  parnames.resize(nparams);

  //For now just ignore where fixed
  ignore_params.resize(nparams);
  for (unsigned int i = 0; i < nknots-1; ++i)
    ignore_params[i] = ifile.isKnotFixed(i);
  //Ignore sigma if we aren't fitting it
  ignore_params[nknots] = (! spec_file.fit_sigma);
2
  //Initialize likelihood information -- priors, data, etc.
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
  unsigned int nproc = MPI::COMM_WORLD.Get_size();

  unsigned int ninitialized;
  std::vector<bool> initialized(false,nproc);
  initialized[0] = true; //Master is initialized
  ninitialized = 1;

  MPI::Status Info;
  int jnk;
  while (ninitialized < nproc) {
    //Wait for a message from a slave sayign they are ready
    MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,MPI::ANY_SOURCE,
			 MPI::ANY_SOURCE,Info);
    int this_tag = Info.Get_tag();
    int this_rank = Info.Get_source();
    if (this_tag == pofd_affine::DONE) {
      //An error was encountered
      return; //Should maybe send DONE message to all others?
    } else if (this_tag == pofd_affine::SENDINITFILE) {
      ifile.sendSelf(MPI::COMM_WORLD, this_rank);
    } else if (this_tag == pofd_affine::SENDLIKESET) {
      likeSet.sendSelf(MPI::COMM_WORLD, this_rank);
    } else if (this_tag == pofd_affine::ISREADY) {
      initialized[this_rank] = true;
    }
    ninitialized = std::count(initialized.begin(), initialized.end(),
			      true);
  }

  //We can free all the data storage, since master doesn't need
  // any of it
  

  //Generate initial parameters for each walker
  
}

void pofdMCMC::initChainsSlave() {
  //Send message to master saying we are ready

  //Copy info
}

void pofdMCMC::initChains() {
  if (rank == 0) initChainsMaster(); else initChainsSlave();
}
