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

  if (rank == 0) {
    //Master node
    ifile.readFile(initfile);
    specFileKnots spec(
  }

}
