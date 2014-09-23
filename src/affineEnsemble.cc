#include<ctime>
#include<cmath>
#include<limits>
#include<sstream>
#include<unistd.h>
#include<mpi.h>

#include "../include/global_settings.h"
#include "../include/affineEnsemble.h"
#include "../include/affineExcept.h"
#include "../include/hdf5utils.h"
#include "../include/hashbar.h"


/////////////////////////////////////

/*!
  \param[in] NWALKERS Number of walkers
  \param[in] NPARAMS  Number of parameters
  \param[in] NSAMPLES Number of samples (across all walkers) to do after burn;
                      realized to the closest larger multiple of nwalkers
  \param[in] INIT_STEPS Number of initialization steps, which are thrown away
                         even before starting burn-in process
  \param[in] INIT_TEMP Temperature factor used during initial steps
  \param[in] MIN_BURN Minimum number of steps before burn-in
  \param[in] FIXED_BURN Do a fixed burn in of length MIN_BURN, not using
                         autocorrelation length to decide when burn in is 
                         finished.  If this is set, INIT_STEPS is ignored.
  \param[in] BURN_MULTIPLE This fraction of autocorrelation steps to add
                            before checking burn-in again
  \param[in] SCALEFAC Scale factor of Z distribution
*/
affineEnsemble::affineEnsemble(unsigned int NWALKERS, unsigned int NPARAMS,
			       unsigned int NSAMPLES, unsigned int INIT_STEPS,
			       double INIT_TEMP, unsigned int MIN_BURN, 
			       bool FIXED_BURN, float BURN_MULTIPLE, 
			       float SCALEFAC) :
  nwalkers(NWALKERS), nparams(NPARAMS), has_any_names(false), 
  scalefac(SCALEFAC), init_steps(INIT_STEPS), init_temp(INIT_TEMP),
  min_burn(MIN_BURN), fixed_burn(FIXED_BURN), burn_multiple(BURN_MULTIPLE), 
  pstep(NPARAMS), is_init(false), has_initStep(false), initStep(NPARAMS),
  has_regenFirstStep(false), regenFirstStep(NPARAMS),
  chains(NWALKERS,NPARAMS), verbosity(0) {
  
  has_name.resize(nparams);
  has_name.assign(nparams, false);
  parnames.resize(nparams);

  //Set number of steps per walker, which won't quite match nsamples
  nsteps = static_cast<unsigned int>(NSAMPLES / static_cast<float>(nwalkers)
				     + 0.999999999999);

  acor_set = false;

  int nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if (nproc < 2)
    throw affineExcept("affineEnsemble", "affineEnsemble",
		       "Must have at least 2 nodes");

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    //Reserve room for queues
    procqueue.setCapacity(nproc);
    stepqueue.setCapacity(nwalkers / 2 + 1);  //Only store one half at a time

    //Master node; slaves don't need this
    naccept.resize(nwalkers);
    naccept.assign(nwalkers, 0);
    nfixed = 0;
    nignore = 0;
    nbonus = 0;
    
    param_state.resize(nparams);
    param_state.assign(nparams, 0);

    //Set random number generator seed from time
    unsigned long long int seed;
    seed = static_cast<unsigned long long int>(time(NULL));
    seed += static_cast<unsigned long long int>(clock());
    rangen.setSeed(seed);
  }
    
}

affineEnsemble::~affineEnsemble() {}

/*!
  \returns True if the ensemble is in a valid state, false otherwise

  Only meaningful for the master node, always return true for slaves
*/
bool affineEnsemble::isValid() const {
  if (rank != 0) return true;
  if (!is_init) return false;
  if (nwalkers < 2) return false; //Need at least 2!
  if (nparams == 0) return false;
  if (nsteps == 0) return false;
  if (min_burn < 10) return false;
  if ((!fixed_burn) && (burn_multiple <= 1.0)) return false;
  if (scalefac <= 0.0) return false;
  if (getNFitParams() == 0) return false; //!< All params fixed
  return true;

}

/*!
  \returns A new z value
*/
float affineEnsemble::generateZ() const {
  //Don't check isValid or if master for speed
  //Generate a random number satisfying g \propto 1/sqrt(z)
  // for z \elem [1/a,a]
  //The inverse cumulative distribution is
  //  G^{-1}(x) = [(a-1)*x+1]^2/a
  //a is scalefac in the nomenclature of this class
  double val;
  val = (scalefac - 1.0) * rangen.flt() + 1.0;
  return val * val / scalefac;
}

/*!
  \param[in] n New number of walkers
  
  Does not preserve contents unless new value is the same as current one.
*/
void affineEnsemble::setNWalkers(unsigned int n) {
  if (nwalkers == n) return;
  if (n == 0)
    throw affineExcept("affineEnsemble", "setNWalkers",
		       "n must be positive");

  chains.setNWalkers(n);
  
  unsigned int nsamples = nsteps * nwalkers;
  nsteps = static_cast<unsigned int>(nsamples / static_cast<float>(n)
				     + 0.999999999999);

  if (rank == 0) {
    stepqueue.setCapacity(n/2+1);  //Only store one half at a time

    naccept.resize(n);
    naccept.assign(n, 0);
  }
  
  nwalkers = n;
}

/*!
  \param[in] n New number of parameters
  
  Does not preserve contents unless new value is the same as current one.
*/
void affineEnsemble::setNParams(unsigned int n) {
  if (nparams == n) return;
  if (n == 0)
    throw affineExcept("affineEnsemble", "setNParams",
		       "n must be positive");

  has_initStep = false;
  initStep.setNParams(n);
  has_regenFirstStep = false;
  regenFirstStep.setNParams(n);

  chains.setNParams(n);
  nfixed = 0;
  nignore = 0;
  nbonus = 0;

  has_name.resize(n);
  has_name.assign(n, false);
  parnames.resize(n);

  if (rank == 0) {
    param_state.resize(n);
    clearParamState();
  }

  nparams = n;
}

/*!
  \returns Maximum likelihood of all steps taken
*/
double affineEnsemble::getMaxLogLike() const {
  if (rank != 0) return std::numeric_limits<double>::quiet_NaN();
  return chains.getMaxLogLike();
}

/*!
  \param[out] val Maximum likelihood of all steps taken
  \param[out] p Parameters corresponding to that likelihood
*/
void affineEnsemble::getMaxLogLikeParam(double& val, paramSet& p) const {
  if (rank != 0) {
    val = std::numeric_limits<double>::quiet_NaN();
    p.setNParams(nparams);
    for (unsigned int i = 0; i < nparams; ++i)
      p[i] = std::numeric_limits<double>::quiet_NaN();
    return;
  }
  chains.getMaxLogLikeParam(val,p);
}

/*!
  \param[in] flag Flag to compare param state to
  \returns Number of params with that flag set
*/
unsigned int affineEnsemble::countParamState(int flag) const {
  unsigned int n = 0;
  for (unsigned int i = 0; i < nparams; ++i)
    if (param_state[i] & flag) ++n;
  return n;
}

void affineEnsemble::clearParamState() {
  if (rank != 0) return;
  param_state.assign(param_state.size(), 0);
  nfixed = nignore = nbonus = 0;
}

/*!
  \returns Number of parameters being fit (nparams minus number of 
    fixed and bonus ones)
*/
unsigned int affineEnsemble::getNFitParams() const {
  return nparams - nfixed - nbonus;
}

/*!
  \returns Number of parameters being used in autocorrelation computation
*/
unsigned int affineEnsemble::getNAcorParams() const {
  return nparams - nignore;
}

/*!
  \returns Number of bonus parameters
*/
unsigned int affineEnsemble::getNBonusParams() const {
  return nbonus;
}

/*!
  \param[in] idx Parameter index of parameter to fix
*/
void affineEnsemble::fixParam(unsigned int idx) {
  if (rank != 0) return;
  param_state[idx] |= mcmc_affine::FIXED;
  param_state[idx] |= mcmc_affine::ACIGNORE;

  nfixed = countParamState(mcmc_affine::FIXED);
  nignore = countParamState(mcmc_affine::ACIGNORE);
}

/*!
  \param[in] idx Parameter index
  \returns True if that parameter is fixed, false otherwise
*/
bool affineEnsemble::isParamFixed(unsigned int idx) const {
  if (rank != 0) return false;
  return param_state[idx] & mcmc_affine::FIXED;
}

/*!
  \param[in] idx Index of parameter to ignore in autocorrelation computation
*/
void affineEnsemble::ignoreParamAcor(unsigned int idx) {
  if (rank != 0) return;
  param_state[idx] |= mcmc_affine::ACIGNORE;
  nignore = countParamState(mcmc_affine::ACIGNORE);
}

/*!
  \param[in] idx Parameter index
  \returns True if that parameter is being ignored in autocorrelation 
    computation, false otherwise
*/
bool affineEnsemble::isParamIgnoredAcor(unsigned int idx) const {
  if (rank != 0) return false;
  return param_state[idx] & mcmc_affine::ACIGNORE;
}

/*!
  \param[in] idx Index of parameter to set as Bonus

  Note that the parameter is automatically set to be ignored
  in the autocorrelation as well.
*/
void affineEnsemble::setParamBonus(unsigned int idx) {
  if (rank != 0) return;
  param_state[idx] |= mcmc_affine::ACIGNORE;
  param_state[idx] |= mcmc_affine::BONUS;

  nignore = countParamState(mcmc_affine::ACIGNORE);
  nbonus = countParamState(mcmc_affine::BONUS);
}

/*!
  \param[in] idx Parameter index
  \returns True if that parameter is a bonus parameter, false otherwise
*/
bool affineEnsemble::isParamBonus(unsigned int idx) const {
  if (rank != 0) return false;
  return param_state[idx] & mcmc_affine::BONUS;
}

/*!
  \param[in] idx Parameter index
  \param[in] parnm Name of parameter
*/
void affineEnsemble::setParamName(unsigned int idx, const std::string& parnm) {
  if (idx >= nparams) 
    throw affineExcept("affineEnsemble", "setParamName", "Invalid index");
  has_any_names = true;
  has_name[idx] = true;
  parnames[idx] = parnm;
}

/*!
  \param[in] idx Parameter index to remove name for
*/
void affineEnsemble::unsetParamName(unsigned int idx) {
  if (idx >= nparams) 
    throw affineExcept("affineEnsemble", "unsetParamName", "Invalid index");
  if (!has_name[idx]) return; //Wasn't set anyways
  has_name[idx] = false;
  has_any_names = has_name[idx];
  for (unsigned int i = 1; i < nparams; ++i)
    has_any_names |= has_name[i];
}

/*!
  \param[in] idx Parameter index
  \returns Name of parameter, or empty string if there is none
*/
std::string affineEnsemble::getParamName(unsigned int idx) const {
  if (idx >= nparams) 
    throw affineExcept("affineEnsemble", "getParamName", "Invalid index");
  if (!has_name[idx]) return std::string();
  return parnames[idx];
}

/*!
  \returns True if the autocorrelation was computed
*/
bool affineEnsemble::computeAcor() const {
  if (rank != 0) return false;
  if (!isValid())
    throw affineExcept("affineEnsemble", "computeAcor",
		       "Sampler not in valid state");
  
  if (acor.size() < nparams) acor.resize(nparams);
  bool success;
  success = chains.getAcorVector(acor, param_state);
  if (success) acor_set = true;
  return success;
}

/*!
  \param[out] retval Autocorrelation length for each parameter.  Resized
    to nparams.  Some entries may not be set (if ignored)
  \return True if the computation was successful
  
  Does not recompute the autocorrelation if already available
*/
bool affineEnsemble::getAcor(std::vector<float>& retval) const {
  if (rank != 0) 
    throw affineExcept("affineEnsemble", "getAcor", "Don't call on slave");
  bool success;
  if (!acor_set) success = computeAcor(); else success=true;
  if (!success) return success;
  retval.resize(nparams);
  for (unsigned int i = 0; i < nparams; ++i)
    retval[i] = acor[i];
  return success;
}

/*!
  \returns True if at least one step has been accepted on any walker
*/
bool affineEnsemble::hasOneStepBeenAccepted() const {
  if (rank != 0) 
    throw affineExcept("affineEnsemble", "hasOneStepBeenAccepted", 
		       "Don't call on slave");
  if (!isValid())
    throw affineExcept("affineEnsemble", "hasOneStepBeenAccepted",
		       "Sampler not in valid state");
  for (unsigned int i = 0; i < nwalkers; ++i)
    if (naccept[i] > 0) return true;
  return false;
}

/*!
  \returns retval Vector of acceptance fraction per walker
*/
void affineEnsemble::getAcceptanceFrac(std::vector<float>& retval) const {
  if (rank != 0) 
    throw affineExcept("affineEnsemble", "getAcceptanceFrac",
		       "Don't call on slave");
  if (!isValid())
    throw affineExcept("affineEnsemble", "computeAcceptanceFrac",
		       "Sampler not in valid state");
  
  retval.resize(nwalkers);
  for (unsigned int i = 0; i < nwalkers; ++i) {
    unsigned int ncurrsteps = chains.getNIters(i);
    if (ncurrsteps == 0)
      retval[i] = std::numeric_limits<float>::quiet_NaN();
    else 
      retval[i] = static_cast<float>(naccept[i]) /
	static_cast<float>(ncurrsteps);
  }
}

/*!
  \param[in] paridx Parameter index
  \returns Mean value of that parameter across all walkers and steps
*/
float affineEnsemble::getParamMean(unsigned int paridx) const {
  return chains.getParamMean(paridx);
}

/*!
  \param[in]  paridx Index of parameter
  \param[out] mean Mean value of parameter
  \param[out] var  Variance of parameter
  \param[out] lowlimit Lower limit of parameter. Returns Nan if non-computable
  \param[out] uplimit  Upper limit of parameter. Returns Nan if non-computable
  \param[in]  conflevel Confidence limit to use (def: 0.683, i.e. 1 sigma)

  Requires temporary storage equal to the total number of steps
*/
void affineEnsemble::getParamStats(unsigned int paridx, float& mean,
				   float& var, float& lowlimit,
				   float& uplimit, float conflevel) const {
  chains.getParamStats(paridx, mean, var, lowlimit, uplimit, conflevel);
}

/*!
  \param[in]  conflevel Confidence limit to use (def: 0.683, i.e. 1 sigma)
  \param[inout] os Stream to write to (def: std::cout)
*/
void affineEnsemble::printStatistics(float conflevel, 
				     std::ostream& os) const {
  if (rank != 0) 
    throw affineExcept("affineEnsemble", "printStatistics",
		       "Don't call on slave");
  if (!isValid())
    throw affineExcept("affineEnsemble", "printStatistics",
		       "Sampler not in valid state");
  std::string parname;
  float mn, var, low, up;
  for (unsigned int i = 0; i < nparams; ++i) {
    if (!has_name[i]) parname = "Unknown"; else
      parname = parnames[i];
    if ((param_state[i] & mcmc_affine::FIXED) != 0) {
      mn = chains.getParamMean(i);
      os << "Parameter: " << parname << " Fixed at value: " << mn << std::endl;
    } else {
      chains.getParamStats(i, mn, var, low, up, conflevel);
      os << "Parameter: " << parname << " Mean: " << mn << " Stdev: "
	 << sqrt(var) << std::endl;
      os << "  lower limit: " << low << " upper limit: "
	 << up << " (" << conflevel * 100.0 << "% limit)";
      if (param_state[i] & mcmc_affine::BONUS)
	os << " (bonus param)";
      os << std::endl;
    }
  }

}

/*!
  Manages the sampling.  Should be called for both master and slave
  processes.
  It proceeds as follows:
    1) Initialize the chains by calling initChains if not already done
    2) Possibly do initsteps steps in each walker, and discard
    3) Then do min_burn steps in each walker
    4) If fixed_burn is not set:
      4a) Compute the autocorrelation length, see if it is acceptable
      4b) Do enough additional steps so that burn_multiple*acor steps have
       been done
      4c) Compute acor again and make sure we have done more than that 
         many steps; if so, finish burn in.
      4d) Otherwise do more steps and try again
    5) Discard the previous steps, using the last ones as the first
       of the new set
    6) Then do nsteps additional steps
  Note that this generally should not be called more than once -- this
  is a driver that tries to do a full fit.  For more fine grained/low level 
  control, use doSteps.
*/
void affineEnsemble::sample() {
  //This routine is where all of the parallel bits happen
  // Different things happen if this is the master node or a slave one;
  //  slave ones just compute likelihoods, the master node suggests new
  //  steps and decides whether to accept them or not
  if (!is_init) initChains();
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) masterSample(); else slaveSample();
}

/*!
  \param[in] nsteps  Number of steps to do in each walker
  \param[in] initsteps Number of initial steps to do then discard.
                       If this is non-zero, the parameters are re-initialized
		       by calling generateInitialPosition with the best
		       step from this set.

  A quick and dirty routine to do a fixed number of steps.
  Assumes that initChains is already called.

*/
void affineEnsemble::doSteps(unsigned int nsteps, unsigned int initsteps) {
  if (!is_init) initChains();

  if (rank == 0) {
    if (!isValid())
      throw affineExcept("affineEnsemble", "doSteps",
			 "Calling with invalid model setup");
    int jnk;

    // Make sure our previous likelihood is valid
    if (verbosity >= 2)
      std::cout << "Computing likelihoods of initializing steps" << std::endl;
    calcLastLikelihood();

    //Do initial steps
    if (initsteps > 0) {
      if (verbosity >= 1) {
	if (verbosity >= 2) 
	  std::cout << "**********************************************"
		    << std::endl;
	std::cout << "Doing " << init_steps << " initial steps per walker"
		  << std::endl;
	if (verbosity >= 2) 
	  std::cout << "**********************************************"
		    << std::endl;
      }
      chains.addChunk(initsteps);
      for (unsigned int i = 0; i < initsteps; ++i) {
	if (verbosity >= 2)
	  std::cout << " Done " << i+1 << " of " << initsteps << " steps"
		    << std::endl;
	doMasterStep(init_temp);
      }

      // Make sure at least one step was accepted, or else something
      // must have gone wrong, like all the parameters being rejected
      // It is possible this could happen legitimately, by accident,
      // but catching bugs in the acceptance is worth the (small) risk
      if (!hasOneStepBeenAccepted())
	throw affineExcept("affineEnsemble", "doSteps",
			   "Failed to accept any initial steps");

      // Get best step, regenerate from that
      paramSet p(nparams);
      double llike;
      chains.getMaxLogLikeParam(llike, p);
      if (verbosity >= 2) {
	std::cout << "**********************************************"
		  << std::endl;
	std::cout << " Best likelihood from initial steps: " << llike
		  << std::endl;
	std::cout << "  For: " << p << std::endl;
	std::cout << " Generating new initial position from that" << std::endl;
      }
      generateInitialPosition(p);

      // Store new initial position
      has_regenFirstStep = true;
      regenFirstStep = p;
    }

    // Follow with main step loop
    if (verbosity >= 1) {
      std::cout << "Doing " << nsteps << " primary steps per walker"
		<< std::endl;
      if (verbosity >= 2)
	std::cout << "**********************************************"
		  << std::endl;
    }
    chains.addChunk(nsteps);
    for (unsigned int i = 0; i < nsteps; ++i) {
      if (verbosity >= 2)
	std::cout << " Done " << i+1 << " of " << nsteps << " steps"
		  << std::endl;
      doMasterStep();
    }
    
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    for (int i = 1; i < nproc; ++i)
      MPI_Send(&jnk, 1, MPI_INT, i, mcmc_affine::STOP, MPI_COMM_WORLD);
  } else
    slaveSample();
}

void affineEnsemble::masterSample() {
  const unsigned int maxhash = 60;

  if (rank != 0) 
    throw affineExcept("affineEnsemble", "masterSample", 
		       "Don't call on slave");

  if (!isValid())
    throw affineExcept("affineEnsemble", "masterSample",
		       "Calling with invalid model setup");
  
  //Make sure we have valid likelihoods
  calcLastLikelihood();

  // Progress bar variable, allocated if needed
  hashbar *progbar = NULL;

  //Do initial steps
  if (init_steps > 0) {
    if (verbosity >= 1) {
      if (verbosity >= 2) 
	std::cout << "**********************************************"
		  << std::endl;
      std::cout << "Doing " << init_steps << " initial steps per walker"
		<< std::endl;
      if (verbosity >= 2) 
	std::cout << "**********************************************"
		  << std::endl;
      if (verbosity == 1) {
	progbar = new hashbar(maxhash, init_steps, '=');
	progbar->initialize(std::cout);
      }
    }
    chains.addChunk(init_steps);
    for (unsigned int i = 0; i < init_steps; ++i) {
      doMasterStep(init_temp);
      if (verbosity == 1)
	progbar->update(i + 1, std::cout);
      else if (verbosity >= 2)
	std::cout << " Done " << i+1 << " of " << init_steps << " steps"
		  << std::endl;
    }
    if (verbosity == 1) {
      progbar->fill(std::cout);
      delete progbar;
    }

    // Make sure at least one step was accepted, or else something
    // must have gone wrong, like all the parameters being rejected
    // It is possible this could happen legitimately, by accident,
    // but catching bugs in the acceptance is worth the (small) risk
    if (!hasOneStepBeenAccepted())
      throw affineExcept("affineEnsemble", "masterSample",
			 "Failed to accept any initial steps");
    
    // Get best step, regenerate from that
    paramSet p(nparams);
    double llike;
    chains.getMaxLogLikeParam(llike, p);
    if (verbosity >= 2) {
      std::cout << "**********************************************"
		<< std::endl;
      std::cout << " Best likelihood from initial steps: " << llike
		<< std::endl;
      std::cout << "  For: " << p << std::endl;
      std::cout << " Generating new initial position from that" << std::endl;
    }
    generateInitialPosition(p);

    // Store new initial position
    has_regenFirstStep = true;
    regenFirstStep = p;
    
    // Get the likelihoods
    calcLastLikelihood();
  }

  // Do burn in.  No progress bar because we don't know how
  // many steps we are going to take
  doBurnIn();

  // Now do the main loop.  If the verbosity level is 1,
  // then output a hash progress bar.  If it's higher, don't,
  // because there will be other output interfering with it.
  // Then do extra steps
  if (verbosity >= 1) {
    std::cout << "Doing " << nsteps << " additional steps per"
	      << " walker, for " << nsteps*nwalkers
	      << " total steps" << std::endl;
    if (verbosity >= 2) 
      std::cout << "**********************************************"
		<< std::endl;
    if (verbosity == 1) {
      progbar = new hashbar(maxhash, nsteps, '=');
      progbar->initialize(std::cout);
    }
  }
  chains.addChunk(nsteps);
  for (unsigned int i = 0; i < nsteps; ++i) {
    // Actual work
    doMasterStep();

    // Update output
    if (verbosity == 1) progbar->update(i + 1, std::cout);
    else if (verbosity >= 2)
      std::cout << " Done " << i+1 << " of " << nsteps << " steps"
		<< std::endl;
  }
  if (verbosity == 1) {
    progbar->fill(std::cout);
    delete progbar;
  } else if (verbosity >= 2)
    std::cout << "Done all steps" << std::endl;

  // This would be a definite problem
  if (!hasOneStepBeenAccepted())
    throw affineExcept("affineEnsemble", "masterSample",
		       "Failed to accept any steps from main loop");
  
  //Tell slaves we are done
  int jnk, nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  for (int i = 1; i < nproc; ++i)
    MPI_Send(&jnk, 1, MPI_INT, i, mcmc_affine::STOP, MPI_COMM_WORLD);
}

/*!
  \returns Maximimum autocorrelation
*/
float affineEnsemble::getMaxAcor() const {
  if (!acor_set) 
    throw affineExcept("affineEnsemble", "getMaxAcor",
		       "Called without acor available");
  float maxval;
  unsigned int start_idx;
  for (start_idx = 0; start_idx < nparams; ++start_idx)
    if ((param_state[start_idx] & mcmc_affine::ACIGNORE) == 0) break;
  if (start_idx == nparams)
    throw affineExcept("affineEnsemble", "getMaxAcor", "All params ignored");
  maxval = acor[start_idx];
  for (unsigned int i = start_idx+1; i < nparams; ++i)
    if (((param_state[i] & mcmc_affine::ACIGNORE) == 0) && 
	acor[i] > maxval) maxval = acor[i];
  return maxval;
}

/*!
 This does the burn in.  Only called from master node
 1) Possibly do init_steps, taking the best one to re-seed the
     initial position for burn-in
 2) Do the burn in 
    a) if fixed_burn is set, just do that many steps
    b) if not, do steps until the number of steps is greater than
        burn_mult * autocorrelation 
 3) Throw away all steps, keeping the last step as the first one of
     the main loop (but don't count it in any statistics, etc.)
*/
void affineEnsemble::doBurnIn() throw(affineExcept) {

  const unsigned int max_acor_iters = 20;
  const unsigned int max_burn_iters = 30;

  if (rank != 0) 
    throw affineExcept("affineEnsemble", "doBurnIn", "Don't call on slave");

  if (verbosity >= 1) {
    if (verbosity >= 2)
      std::cout << "**********************************************"
		<< std::endl;
    std::cout << "Starting burn-in process" << std::endl;
    if (verbosity >= 2)
      std::cout << "**********************************************"
		<< std::endl;
  }

  //First do min_burn steps in each walker
  chains.addChunk(min_burn);
  for (unsigned int i = 0; i < min_burn; ++i) {
    if (verbosity >= 2)
      std::cout << " Done " << i+1 << " of " << min_burn << " steps"
		<< std::endl;
    doMasterStep();
  }

  // This would be a problem
  if (!hasOneStepBeenAccepted())
    throw affineExcept("affineEnsemble", "doBurnIn",
		       "Failed to accept any steps from initial burn in");

  if (fixed_burn) {
    if (verbosity >= 1) {
      if (verbosity >= 2)
	std::cout << "**********************************************"
		  << std::endl;
      std::cout << "Did fixed burn in of " << min_burn << " steps per"
		<< " walker." << std::endl;
      //Try to compute the acor and report it, but don't worry about
      // it if we can't
      bool acor_success = computeAcor();
      if (acor_success) 
	std::cout << " Maximum autocorrelation is: "
		  << getMaxAcor() << std::endl;
      if (verbosity >= 2)
	std::cout << "**********************************************"
		  << std::endl;
    }
  } else {
    // Autocorrelation based burn-in test
    // Compute the autocorrelation -- we may need to add more steps
    // to get this to work
    bool acor_success = computeAcor();
  
    //If it failed, we need more steps.  Add 25% of min burn
    // and try again.  Refuse to do this more than max_acor_iters times.
    unsigned int nextra = 10;
    if (!acor_success) {
      nextra = static_cast<unsigned int>(min_burn * 0.25);
      if (nextra < 10) nextra = 10; // Always do at least 10 steps
      for (unsigned int i = 0; i < max_acor_iters; ++i) {
	if (verbosity >= 1) {
	  if (verbosity >= 2)
	    std::cout << "**********************************************"
		      << std::endl;
	  std::cout << " Failed to compute acor after " << chains.getMinNIters()
		    << " steps.  Adding " << nextra << " more" << std::endl;
	  if (verbosity >= 2)
	    std::cout << "**********************************************"
		      << std::endl;
	}
	chains.addChunk(nextra);
	for (unsigned int i = 0; i < nextra; ++i) {
	  if (verbosity >= 2)
	    std::cout << " Done " << i+1 << " of " << nextra << " steps"
		      << std::endl;
	  doMasterStep();
	}
	acor_success = computeAcor();
	if (acor_success) break;
      }
      // We were completely unable to get the autocorrelation length
      if (!acor_success)
	throw affineExcept("affineEnsemble", "doBurnIn",
			   "Can't compute acor; increase min_burn");
    }
  
    //Okay, we have an acor estimate of some sort, even though
    // we probably aren't burned in yet.
    float max_acor = getMaxAcor();
    if (verbosity >= 1) {
      if (verbosity >= 2)
	std::cout << "**********************************************"
		  << std::endl;
      std::cout << " After " << chains.getMinNIters() 
		<< " steps, maximum autocorrelation is: "
		<< max_acor << std::endl;
      if (verbosity >= 2)
	std::cout << "**********************************************"
		  << std::endl;
    }
  
    unsigned int nsteps = chains.getMinNIters();
    unsigned int nminsteps = 
      static_cast<unsigned int>(burn_multiple*max_acor + 0.999999999999999);
    bool burned_in = (nsteps > nminsteps);
    unsigned int nburn_iters = 0;
    while (!burned_in) {
      if (nburn_iters > max_burn_iters)
	throw affineExcept("affineEnsemble", "doBurnIn",
			   "Failed to converge in burn in");

      //Figure out how many more steps to do.  Do half of the number
      // estimated
      unsigned int nmore = (nminsteps - nsteps) / 2 + 1;
      if (verbosity >= 1) {
	if (verbosity >= 2)
	  std::cout << "**********************************************"
		    << std::endl;
	std::cout << " Doing " << nmore << " additional steps" << std::endl;
	if (verbosity >= 2)
	  std::cout << "**********************************************"
		    << std::endl;
      }
      chains.addChunk(nmore);
      for (unsigned int i = 0; i < nmore; ++i) {
	if (verbosity >= 2)
	  std::cout << " Done " << i+1 << " of " << nmore << " steps"
		    << std::endl;
	doMasterStep();
      }
      nsteps += nmore;
      
      //Update acor, same as before
      acor_success = computeAcor();
      
      //Again, add more steps if we must to get acor; even though
      // we got it successfully before, it may not converge this time,
      // so we may need to do additional steps.
      nextra = static_cast<unsigned int>(min_burn * 0.25);
      if (nextra < 10) nextra = 10;
      for (unsigned int i = 0; i < max_acor_iters; ++i) {
	if (verbosity >= 1) {
	  if (verbosity >= 2)
	    std::cout << "**********************************************"
		      << std::endl;
	  std::cout << " Failed to compute acor after " << chains.getMinNIters()
		    << " steps.  Adding " << nextra << " more" << std::endl;
	  if (verbosity >= 2)
	    std::cout << "**********************************************"
		      << std::endl;
	}
	chains.addChunk(nextra);
	for (unsigned int i = 0; i < nextra; ++i) {
	  if (verbosity >= 2)
	    std::cout << " Done " << i+1 << " of " << nextra << " steps"
		      << std::endl;
	  doMasterStep();
	}
	acor_success = computeAcor();
	if (acor_success) break;
	nsteps += nextra;
      }
      if (!acor_success)
	throw affineExcept("affineEnsemble", "doBurnIn",
			   "Can't compute acor; increase min_burn");
      
      // Ok, we have an acor.  Now check to see if we have enough!
      max_acor = getMaxAcor();
      nminsteps = 
	static_cast<unsigned int>(burn_multiple*max_acor + 0.999999999999999);
      burned_in = (nsteps > nminsteps);
      nburn_iters += 1;
    }
  }

  if (verbosity >= 1) {
    if (verbosity >= 2)
      std::cout << "**********************************************"
		<< std::endl;
    std::cout << "Burned in after: " << chains.getMinNIters() 
	      << " steps" << std::endl;
    if (verbosity >= 2)
      std::cout << "**********************************************"
		<< std::endl;
  }

  //Throw away burn in, keeping last step as first step of new one
  // We will not count that last step as part of our statistics
  chains.clearPreserveLast(); 
  for (unsigned int i = 0; i < naccept.size(); ++i)
    naccept[i] = 0;
  chains.setSkipFirst();
}

/*!
  Redos the likelihood compuatation for the most recent step.
  The idea behind this is that sometimes -- like when initializing
  the chains -- you need to compute the likelihood of a step, but it's
  computationally difficult.  This provides a mechanism for doing that 
  which is called by the internal sampling routines to make sure
  they are working from valid likelihoods.

  This can not be called from slave nodes.

  The slave nodes should be in the middle of slaveSample for this
  to work.
*/
void affineEnsemble::calcLastLikelihood() {

  if (rank != 0)
    throw affineExcept("affineEnsemble", "calcLastLikelihood",
		       "Should only be called from master node");

  MPI_Status Info;
  bool parrej;
  int jnk, this_rank, this_tag, ismsg;
  
  // Set up variable for passing stuff around
  if (pstep.oldStep.getNParams() < nparams) 
    pstep.oldStep.setNParams(nparams);
  if (pstep.newStep.getNParams() < nparams) 
    pstep.newStep.setNParams(nparams);
  
  // Set up a queue of things that need updating
  affineQueue<unsigned int> needCalc(nwalkers);
  for (unsigned int i = 0; i < nwalkers; ++i)
    needCalc.push(i);
  
  int nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  
  unsigned int ndone = 0;
  unsigned int this_update;
  // This is pretty much a slightly modified version of what happens
  // in emptyMasterQueue
  while (ndone < nwalkers) {
    // If there are available procs, send them a step to do
    if (!procqueue.empty())
      try {
	for (unsigned int i = 0; i < needCalc.size(); ++i) {
	  if (procqueue.empty()) break; // No more procs available
	  this_rank = procqueue.pop();
	  this_update = needCalc.pop();
	  chains.getLastStep(this_update, pstep.newStep, pstep.newLogLike);
	  pstep.update_idx = this_update;
	  MPI_Send(&jnk, 1, MPI_INT, this_rank, mcmc_affine::SENDINGPOINT,
		   MPI_COMM_WORLD);
	  pstep.sendSelf(MPI_COMM_WORLD, this_rank);
	}
      } catch (const affineExcept& ex) {
	for (int i = 1; i < nproc; ++i)
	  MPI_Send(&jnk, 1, MPI_INT, i, mcmc_affine::STOP, MPI_COMM_WORLD);
	throw ex;
      }
    
    //No available procs, so wait for some sort of message
    // We use Iprobe to see if a message is available, and if not
    // sleep for 1/100th of a second and try again
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &ismsg, &Info);
    while (ismsg == 0) {
      usleep(mcmc_affine::usleeplen); //Sleep for 1/100th of a second
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &ismsg, &Info);
    }
    
    //There is a message, grab it and figure out what to do
    MPI_Recv(&jnk, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, 
	     MPI_COMM_WORLD, &Info);
    this_tag = Info.MPI_TAG;
    this_rank = Info.MPI_SOURCE;
    
    if (this_tag == mcmc_affine::ERROR) {
      //Slave had a problem, so send stop message to all slaves
      for (int i = 1; i < nproc; ++i)
	MPI_Send(&jnk, 1, MPI_INT, i, mcmc_affine::STOP, MPI_COMM_WORLD);
      throw affineExcept("affineEnsemble", "calcLastLikelihood",
			 "Problem encountered in slave");
    } else if (this_tag == mcmc_affine::SENDINGRESULT) {
      //Slave finished a computation, is sending us the result
      // First, find out if the parameter set was rejected
      MPI_Recv(&parrej, 1, MPI::BOOL, this_rank, mcmc_affine::SENDINGPARREJ,
	       MPI_COMM_WORLD, &Info);
      // And get the actual step
      pstep.receiveCopy(MPI_COMM_WORLD, this_rank);
      
      if (parrej) {
	// If step was rejected, treat this as an error and crash
	for (int i = 1; i < nproc; ++i)
	  MPI_Send(&jnk, 1, MPI_INT, i, mcmc_affine::STOP, MPI_COMM_WORLD);
	std::stringstream errstr;
	errstr << "Unable to compute likelihood of requested step: "
	       << std::endl << " Walker: " << pstep.update_idx << std::endl
	       << pstep.newStep;
	throw affineExcept("affineEnsemble", "calcLastLikelihood",
			   errstr.str());
      }
      
      // Update likelihood
      chains.replaceLastStep(pstep.update_idx, pstep.newStep, 
			     pstep.newLogLike);
      
      ndone += 1;
    } else if (this_tag == mcmc_affine::REQUESTPOINT) {
      procqueue.push(this_rank);
    } else {
      std::stringstream sstream;
      sstream << "Master got unexpected message from " << this_rank
	      << " with code: " << this_tag;
    }
  }

}

/*!
  \param[in] p Parameters to compute log likelihood for
  \returns The log likelihood, or nan if the parameters were rejected

  This is a convenience function -- generally the other version should
  be used as it can differentiate between a failed computation (which
  may return nan), and a parameter set that was simply rejected.
*/
double affineEnsemble::getLogLike(const paramSet& p) {
  bool parrej;
  double loglike = getLogLike(p, parrej);
  if (parrej)
    return std::numeric_limits<double>::quiet_NaN();
  else
    return loglike;
}


/*!
  \param[in] idx1  Walker that we are suggesting a new step for
  \param[in] idx2  Walker that we are combining with idx1
  \param[out] prstep Newly proposed step for idx1

  This generates a new step, checking the parameter limits to make sure
  they are obeyed, and obeying fixed parameters, and setting bonus
  parameters to NaN.
*/
void affineEnsemble::generateNewStep(unsigned int idx1, unsigned int idx2,
				     proposedStep& prstep) const 
  throw (affineExcept) {

  // Maximum number of times we will try to generate a new step
  const unsigned int maxiters = 40;

  bool is_valid;
  unsigned int iter;
  float val, omz;
  double tmp;

  // We need storage for three steps -- the old one, the one we are combining,
  // and the proposed new one.  We can use the oldStep/newStep fields of
  // prstep for the other two, but need internal storage for the middle
  if (prstep.oldStep.getNParams() < nparams) prstep.oldStep.setNParams(nparams);
  if (prstep.newStep.getNParams() < nparams) prstep.newStep.setNParams(nparams);
  if (params_tmp.getNParams() < nparams) params_tmp.setNParams(nparams);

  // Get steps to combine
  chains.getLastStep(idx1, prstep.oldStep, prstep.oldLogLike);
  chains.getLastStep(idx2, params_tmp, tmp);
  prstep.update_idx = idx1;
  prstep.newLogLike = std::numeric_limits<double>::quiet_NaN();

  // Start trying to generate new step
  prstep.z = generateZ();
  omz = 1.0 - prstep.z;
  for (unsigned int i = 0; i < nparams; ++i) 
    if (param_state[i] & mcmc_affine::FIXED) {
      //Fixed parameter, keep previous
      prstep.newStep.setParamValue(i, prstep.oldStep[i]);
    } else if (param_state[i] & mcmc_affine::BONUS) {
      prstep.newStep.setParamValue(i, std::numeric_limits<double>::quiet_NaN());
    } else {
      val = prstep.z * prstep.oldStep[i] + omz * params_tmp[i];
      prstep.newStep.setParamValue(i, val);
    }

  is_valid = areParamsValid(prstep.newStep);
  iter = 0;
  while (!is_valid) {
    if (iter > maxiters) {
      // We were unable to generate a new step after maxiters iters.
      // This shouldn't happen unless there's a bug, but check to see
      // if the input steps weren't valid.
      std::stringstream errstr;
      if (!areParamsValid(prstep.oldStep)) {
	errstr << "Unable to generate new step; input step from idx1 was"
	       << " also not valid" << std::endl;
	errstr << " " << prstep.oldStep;
	throw affineExcept("affineEnsemble", "generateNewStep",
			   errstr.str());
      }
      if (!areParamsValid(params_tmp)) {
	errstr << "Unable to generate new step; input step from idx2 was"
	       << " also not valid" << std::endl;
	errstr << " " << params_tmp;
	throw affineExcept("affineEnsemble", "generateNewStep",
			   errstr.str());
      }
      errstr << "Unable to generate new step after " << maxiters << " tries";
      throw affineExcept("affineEnsemble", "generateNewStep",
			 errstr.str());
    }

    // Try a new one
    prstep.z = generateZ();
    omz = 1.0 - prstep.z;
    for (unsigned int i = 0; i < nparams; ++i) 
      if (param_state[i] & mcmc_affine::FIXED) {
	//Fixed parameter, keep previous
	prstep.newStep.setParamValue(i, prstep.oldStep[i]);
      } else if (param_state[i] & mcmc_affine::BONUS) {
	// Don't have to do anything; should still be NaN
      } else {
	val = prstep.z * prstep.oldStep[i] + omz * params_tmp[i];
	prstep.newStep.setParamValue(i, val);
      }
    is_valid = areParamsValid(prstep.newStep);

    ++iter;
  }

}

/*
  \param[in] temp Temperature to use for step acceptance
*/
void affineEnsemble::doMasterStep(double temp) throw (affineExcept) {
  //We use a pushdown stack to keep track of which things to do
  //To make this parallel, we must split the number of walkers in
  // half, and update each half from the opposite one
  //Fill in queue
  if (rank != 0)
    throw affineExcept("affineEnsemble", "doMasterStep",
		       "Should only be called from master node");
  if (!stepqueue.empty())
    throw affineExcept("affineEnsemble", "doMasterStep",
		       "step queue should be empty at start");
  if (temp <= 0)
    throw affineExcept("affineEnsemble", "doMasterStep", "Invalid temperature");

  unsigned int minidx, maxidx;
  std::pair<unsigned int, unsigned int> pr;
  //Do the first half
  minidx = nwalkers/2; //Generate from second half
  maxidx = nwalkers;
  if (verbosity >= 3)
    std::cout << "  Generating new steps for 0:" << nwalkers/2
	      << " from " << minidx << ":" << maxidx << std::endl;
  for (unsigned int i = 0; i < nwalkers / 2; ++i) {
    pr.first  = i;
    pr.second = rangen.selectFromRange(minidx, maxidx);
    stepqueue.push(pr);
  }

  //Now run that
  if (verbosity >= 3)
    std::cout << "  Evaluating likelihoods" << std::endl;
  emptyMasterQueue(temp);

  //Then the second
  minidx = 0;
  maxidx = nwalkers/2;
  if (verbosity >= 3)
    std::cout << "  Generating new steps for " << nwalkers/2 << ":" << nwalkers
	      << " from " << minidx << ":" << maxidx << std::endl;
  for (unsigned int i = nwalkers/2; i < nwalkers; ++i) {
    pr.first  = i;
    pr.second = rangen.selectFromRange(minidx, maxidx);
    stepqueue.push(pr);
  }
  if (verbosity >= 3)
    std::cout << "  Evaluating likelihoods" << std::endl;
  emptyMasterQueue(temp);

  // The fact that we took another step means that acor is now invalid
  acor_set = false;
}

/*
  \param[in] temp Temperature to use for step acceptance

  Does all steps currently in the queue
*/
void affineEnsemble::emptyMasterQueue(double temp) throw (affineExcept) {
  //Step loop
  MPI_Status Info;
  int jnk, this_rank, this_tag, nproc, ismsg;
  unsigned int ndone, nsteps;
  std::pair<unsigned int, unsigned int> pr;
  double prob; // Probability of accepting a step
  double rval; // Random comparison value
  bool parrej; // Was a parameter set rejected in the likelihood call?

  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  nsteps = stepqueue.size();
  ndone = 0;

  double inv_temp = 1.0 / temp;

  if (nsteps == 0) return; //Nothing to do...

  while (ndone < nsteps) {
    //First, if there are available procs, send them a step if we can
    if (!procqueue.empty()) 
      try {
	for (unsigned int i = 0; i < stepqueue.size(); ++i) {
	  if (procqueue.empty()) break; //No more available procs
	  this_rank = procqueue.pop();
	  
	  //Figure out next thing to update
	  pr = stepqueue.pop();
	  //Generate actual value into pstep
	  generateNewStep(pr.first, pr.second, pstep);
	  
	  if (verbosity >= 3)
	    std::cout << "  Evaluating new step for walker: " << pr.first
		      << " using slave: " << this_rank << std::endl;
	  MPI_Send(&jnk, 1, MPI_INT, this_rank, mcmc_affine::SENDINGPOINT,
		   MPI_COMM_WORLD);
	  pstep.sendSelf(MPI_COMM_WORLD, this_rank);
	}
      } catch (const affineExcept& ex) {
	for (int i = 1; i < nproc; ++i)
	  MPI_Send(&jnk, 1, MPI_INT, i, mcmc_affine::STOP, MPI_COMM_WORLD);
	throw ex;
      }

    //No available procs, so wait for some sort of message
    // We use Iprobe to see if a message is available, and if not
    // sleep for 1/100th of a second and try again
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &ismsg, &Info);
    while (ismsg == 0) {
      usleep(mcmc_affine::usleeplen); //Sleep and then check again
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &ismsg, &Info);
    }

    //There is a message, grab it and figure out what to do
    MPI_Recv(&jnk, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, 
	     MPI_COMM_WORLD, &Info);
    this_tag = Info.MPI_TAG;
    this_rank = Info.MPI_SOURCE;
    if (this_tag == mcmc_affine::ERROR) {
      //Slave had a problem
      //Send stop message to all slaves
      for (int i = 1; i < nproc; ++i)
	MPI_Send(&jnk, 1, MPI_INT, i, mcmc_affine::STOP, MPI_COMM_WORLD);
      throw affineExcept("affineEnsemble", "emptyMasterQueue",
			 "Problem encountered in slave");
    } else if (this_tag == mcmc_affine::SENDINGRESULT) {
      //Slave finished a computation, is sending us the result
      //Note we wait for a REQUESTPOINT to actually add it to the
      // list of available procs, which makes initialization easier

      // First, find out if the parameter set was rejected
      MPI_Recv(&parrej, 1, MPI::BOOL, this_rank, mcmc_affine::SENDINGPARREJ,
	       MPI_COMM_WORLD, &Info);
      // And get the actual step
      pstep.receiveCopy(MPI_COMM_WORLD, this_rank);
      
      if (parrej) {
	// Parameters were invalid in some way, so automatically reject
	// the step.  
	if (verbosity >= 3)
	    std::cout << "  Rejecting new step as likelihood rejected params"
		      << std::endl;
	chains.addNewStep(pstep.update_idx, pstep.oldStep,
			  pstep.oldLogLike);
      } else {
	// It liked the params, so we have to figure out if we want to
	// keep the results
	// First, make sure the values are valid
	if (std::isinf(pstep.newLogLike) ||
	    std::isnan(pstep.newLogLike)) {
	  //That's not good -- exit
	  for (int i = 1; i < nproc; ++i)
	    MPI_Send(&jnk, 1, MPI_INT, i, mcmc_affine::STOP, MPI_COMM_WORLD);
	  throw affineExcept("affineEnsemble", "emptyMasterQueue",
			     "Invalid likelihood");
	}

	//Decide whether to accept the step or reject it, add to
	// chunk.  Note that the acceptance rule has a dependence
	// on z, the strech step, so isn't just the same as for MH MCMC.
	// In particular, we don't always accept a step with better logLike!
	//The probability of acceptance is
	// min(1, z^(n-1) P(new) / P(old)
	//We also allow this to happen at a different temperature.
	prob = exp((nparams - nfixed - 1) * log(pstep.z) + 
		   inv_temp * (pstep.newLogLike - pstep.oldLogLike));
	if (verbosity >= 3) {
	  std::cout << "  Got new step for " << pstep.update_idx 
		    << " from slave " << this_rank << std::endl;
	  std::cout << "   " << pstep << std::endl;
	  std::cout << "   Delta likelihood: "
		    << pstep.newLogLike - pstep.oldLogLike << std::endl;
	}

	if (prob >= 1) {
	  //Will always be accepted; this is the min(1, part
	  // This saves us a call to rangen.doub
	  if (verbosity >= 3)
	    std::cout << "  Accepting new step" << std::endl;
	  chains.addNewStep(pstep.update_idx, pstep.newStep,
			    pstep.newLogLike);
	  naccept[pstep.update_idx] += 1;
	} else {
	  rval = rangen.doub();
	  if (rval < prob) {
	    //Accept!
	    if (verbosity >= 3)
	      std::cout << "  Accepting new step (prob: "
			<< prob << " rval: " << rval <<")" << std::endl;
	    chains.addNewStep(pstep.update_idx, pstep.newStep,
			      pstep.newLogLike);
	    naccept[pstep.update_idx] += 1;
	  } else {
	    //reject, keep old step
	    if (verbosity >= 3)
	      std::cout << "  Rejecting new step (prob: "
			<< prob << " rval: " << rval <<")" << std::endl;
	    chains.addNewStep(pstep.update_idx, pstep.oldStep,
			      pstep.oldLogLike);
	  }
	}
      }
      ndone += 1;
    } else if (this_tag == mcmc_affine::REQUESTPOINT) {
      procqueue.push(this_rank);
    } else {
      std::stringstream sstream;
      sstream << "Master got unexpected message from " << this_rank
	      << " with code: " << this_tag;
      //Tell slaves to stop
      for (int i = 1; i < nproc; ++i)
	MPI_Send(&jnk, 1, MPI_INT, i, mcmc_affine::STOP, MPI_COMM_WORLD);
      throw affineExcept("affineEnsemble", "emptyMasterQueue", sstream.str());
    }
  }

}

void affineEnsemble::slaveSample() {
  //This one just sits here and computes likelihoods
  if (rank == 0)
    throw affineExcept("affineEnsemble", "slaveSample",
		       "Should not be called from master node");
  
  int jnk;
  MPI_Status Info;
  bool parrej; // Parameter set rejected

  try {
    while (true) { //Runs until we get a STOP message
      //Ask for a new point
      MPI_Send(&jnk, 1, MPI_INT, 0, mcmc_affine::REQUESTPOINT,
	       MPI_COMM_WORLD);

      //Now wait for one
      MPI_Recv(&jnk, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Info);

      int tag = Info.MPI_TAG;
      if (tag == mcmc_affine::STOP) return; //Master says -- we're done!
      else if (tag == mcmc_affine::SENDINGPOINT) {
	pstep.receiveCopy(MPI_COMM_WORLD, 0);
	
	// Compute likelihood
	pstep.newLogLike = getLogLike(pstep.newStep, parrej);

	// Fill bonus params (if any)
	fillBonusParams(pstep.newStep, parrej);

	//Send it back
	MPI_Send(&jnk, 1, MPI_INT, 0, mcmc_affine::SENDINGRESULT,
		 MPI_COMM_WORLD);
	MPI_Send(const_cast<bool*>(&parrej), 1, MPI::BOOL, 0, 
		 mcmc_affine::SENDINGPARREJ, MPI_COMM_WORLD);
	pstep.sendSelf(MPI_COMM_WORLD, 0);
      } else {
	std::cerr << "Unexpected message from master: "
                  << tag << " in slave: " << rank << std::endl;
	MPI_Send(&jnk, 1, MPI_INT, 0, mcmc_affine::ERROR,
		 MPI_COMM_WORLD);
	return;
      }
    }
  } catch (const affineExcept& ex) {
    std::cerr << "Error encountered for process: " << rank << std::endl;
    std::cerr << ex << std::endl;
    MPI_Send(&jnk, 1, MPI_INT, 0, mcmc_affine::ERROR, MPI_COMM_WORLD);
    return;
  } catch (const std::bad_alloc& ba) {
    std::cerr << "Allocation error encountered for process: " << rank 
	      << std::endl;
    std::cerr << ba.what() << std::endl;
    MPI_Send(&jnk, 1, MPI_INT, 0, mcmc_affine::ERROR, MPI_COMM_WORLD);
    return;
  }
}

/*!
  \param[inout] os Stream to write to
*/
void affineEnsemble::writeToStream(std::ostream& os) const {
  if (rank != 0) return;

  os << "Number of walkers: " << nwalkers << std::endl;
  os << "Number of parameters: " << nparams;
  if (verbosity >= 2) {
    if (nfixed > 0) 
      os << std::endl << " Number of fixed parameters: " << nfixed;
    if (nignore > 0) 
      os << std::endl << " Number of ignored parameters: " << nignore;
    if (nbonus > 0) 
      os << std::endl << " Number of bonus parameters: " << nbonus;
    if (has_any_names) {
      os << std::endl << "Parameter names: ";
      for (unsigned int i = 0; i < nparams; ++i)
	if (has_name[i]) os << std::endl << " " << i << ": " 
			    << parnames[i];
    }
  }

  if (init_steps > 0) {
    os << std::endl << "Will take initial steps to generate new starting"
       << " position";
    os << std::endl << " Number of initial steps: " << init_steps;
    if (init_temp != 1)
      os << std::endl << " Initial step temperature: " << init_temp;
  }
  if (min_burn > 0) {
    if (fixed_burn) {
      os << std::endl << "Will do fixed burn in of " << min_burn 
	 << " steps per walker";
    } else {
      os << std::endl << "Will do autocorrelation based burn-in";
      os << std::endl << " Minimum number of burn in steps: " << min_burn 
	 << " per walker";
      if (verbosity >= 2)
	os << std::endl << " Burn multiple: " << burn_multiple;
    }
  }

  if (nsteps > 0)
    os << std::endl << "Number of main loop steps per walker: " << nsteps;
}

/*!
  \param[in] filename File to write to as text
*/
void affineEnsemble::writeToFile(const std::string& filename) const {
  if (rank != 0) 
    throw affineExcept("affineEnsemble", "writeToFile",
		       "Should only be called from master node");  
  chains.writeToFile(filename);
}

/*!
  \param[in] objid HDF5 group to write to
*/
void affineEnsemble::writeToHDF5Handle(hid_t objid) const {
  if (rank != 0) 
    throw affineExcept("affineEnsemble", "writeToHDF5Handle",
		       "Should only be called from master node");  
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("affineEnsemble", "writeToHDF5Handle",
		       "Input handle is not valid");

  // This will get subgrouped in a few ways;
  //  AffineSettings has stuff related to how the affineMCMC worked
  //  ParamInfo has things like the names of the params, which were fixed,
  //    initial values, and limits
  //  The actual chains are written to the Chains group

  // Write some attributes
  hsize_t adims;
  hid_t mems_id, att_id;
  hbool_t btmp;

  ////////////////
  // Start with affineSettings
  hid_t groupid;
  groupid = H5Gcreate(objid, "AffineSettings", H5P_DEFAULT, H5P_DEFAULT, 
		      H5P_DEFAULT);
  if (H5Iget_ref(groupid) < 0)
    throw affineExcept("affineEnsemble", "writeToHDF5Handle",
		       "Failed to create AffineSettings HDF5 group");

  // First, simple 1 element objects
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate2(groupid, "ScaleFactor", H5T_NATIVE_FLOAT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_FLOAT, &scalefac);
  H5Aclose(att_id);
  att_id = H5Acreate2(groupid, "NInitSteps", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &init_steps);
  H5Aclose(att_id);
  att_id = H5Acreate2(groupid, "MinBurn", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &min_burn);
  H5Aclose(att_id);
  btmp = static_cast<hbool_t>(fixed_burn);
  att_id = H5Acreate2(groupid, "FixedBurn", H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, &btmp);
  H5Aclose(att_id);
  att_id = H5Acreate2(groupid, "BurnMultiple", H5T_NATIVE_FLOAT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_FLOAT, &burn_multiple);
  H5Aclose(att_id);
  att_id = H5Acreate2(groupid, "NSteps", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &nsteps);
  H5Aclose(att_id);
  H5Sclose(mems_id);

  H5Gclose(groupid);

  ////////////////
  // now ParamInfo
  groupid = H5Gcreate(objid, "ParamInfo", H5P_DEFAULT, H5P_DEFAULT, 
		      H5P_DEFAULT);
  if (H5Iget_ref(groupid) < 0)
    throw affineExcept("affineEnsemble", "writeToHDF5Handle",
		       "Failed to create ParamInfo HDF5 group");
  // Again, start with single value stuff
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, NULL);

  att_id = H5Acreate2(groupid, "NFixed", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &nfixed);
  H5Aclose(att_id);
  att_id = H5Acreate2(groupid, "NIgnore", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &nignore);
  H5Aclose(att_id);
  att_id = H5Acreate2(groupid, "NBonus", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &nbonus);
  H5Aclose(att_id);
  H5Sclose(mems_id);

  // Parameter state (fixed, ignored, bonus) information
  adims = nparams;
  mems_id = H5Screate_simple(1, &adims, NULL);
  hbool_t *batmp;
  batmp = new hbool_t[nparams];
  // Fixed first
  for (unsigned int i = 0; i < nparams; ++i) 
    batmp[i] = (param_state[i] & mcmc_affine::FIXED) != 0;
  att_id = H5Acreate2(groupid, "IsParamFixed", H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, batmp);
  H5Aclose(att_id);
  // Ignored params
  for (unsigned int i = 0; i < nparams; ++i) 
    batmp[i] = (param_state[i] & mcmc_affine::ACIGNORE) != 0;
  att_id = H5Acreate2(groupid, "IsParamAutocorrIgnore", H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, batmp);
  // Bonus params
  for (unsigned int i = 0; i < nparams; ++i) 
    batmp[i] = (param_state[i] & mcmc_affine::BONUS) != 0;
  att_id = H5Acreate2(groupid, "IsParamBonus", H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, batmp);
  delete[] batmp;
  H5Aclose(att_id);
  H5Sclose(mems_id);
 
  // Initial values
  if (has_initStep) {
    float *ftmp;
    ftmp = new float[nparams];
    for (unsigned int i = 0; i < nparams; ++i) 
      ftmp[i] = initStep[i];
    hdf5utils::writeDataFloats(groupid, "InitialPosition",
			       nparams, ftmp);
    delete[] ftmp;
  }
  if (has_regenFirstStep) {
    float *ftmp;
    ftmp = new float[nparams];
    for (unsigned int i = 0; i < nparams; ++i) 
      ftmp[i] = regenFirstStep[i];
    hdf5utils::writeDataFloats(groupid, "RegeneratedFirstStep",
			       nparams, ftmp);
    delete[] ftmp;
  }

  // Names of parameters
  if (has_any_names) 
    hdf5utils::writeAttStrings(groupid, "ParamNames", parnames);

  H5Gclose(groupid);

  // And the actual chains and likelihoods to the Chains group
  groupid = H5Gcreate(objid, "Chains", H5P_DEFAULT, H5P_DEFAULT, 
		      H5P_DEFAULT);
  if (H5Iget_ref(groupid) < 0)
    throw affineExcept("affineEnsemble", "writeToHDF5Handle",
		       "Failed to create Chains HDF5 group");

  // Number of accepted steps
  hdf5utils::writeDataUnsignedInts(groupid, "NAccept", naccept);

  // Chains
  chains.writeToHDF5Handle(groupid);

  H5Gclose(groupid);
}

/*!
  \param[in] filename File to write to as HDF5
*/
void affineEnsemble::writeToHDF5(const std::string& filename) const {

  if (rank != 0) 
    throw affineExcept("affineEnsemble", "writeToHDF5",
		       "Should only be called from master node");  

  // Create
  hid_t file_id;
  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
		      H5P_DEFAULT);

  if (H5Iget_ref(file_id) < 0) {
    H5Fclose(file_id);
    throw affineExcept("affineEnsemble", "writeToHDF5",
		       "Failed to open HDF5 file to write");
  }

  //Write
  writeToHDF5Handle(file_id);

  // All done
  H5Fclose(file_id);
}

/*!
  \param[inout] os Stream to write to
  \param[in] a Ensemble to write
*/
std::ostream& operator<<(std::ostream& os, const affineEnsemble& a) {
  a.writeToStream(os);
  return os;
}
