#include<ctime>
#include<cmath>
#include<limits>
#include<sstream>
#include<mpi.h>

#include<global_settings.h>
#include<affineEnsemble.h>
#include<affineExcept.h>

/////////////////////////////////////

/*!
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
affineEnsemble::affineEnsemble( unsigned int NWALKERS, unsigned int NPARAMS,
				unsigned int NSAMPLES, unsigned int INIT_STEPS,
				unsigned int MIN_BURN,
				double BURN_MULTIPLE, double SCALEFAC ) :
  nwalkers(NWALKERS), nparams(NPARAMS), has_any_names(false), 
  scalefac(SCALEFAC), init_steps(INIT_STEPS), min_burn(MIN_BURN), 
  burn_multiple( BURN_MULTIPLE ), pstep(NPARAMS), chains(NWALKERS,NPARAMS) {
  
  has_name.resize(nparams);
  has_name.assign(nparams, false);
  parnames.resize(nparams);

  nignore = 0;

  //Set number of steps per walker, which won't quite match nsamples
  nsteps = static_cast<unsigned int>( NSAMPLES/static_cast<double>(nwalkers)
				      + 0.999999999999 );


  acor_set = false;

  unsigned int nproc = MPI::COMM_WORLD.Get_size();
  if (nproc < 2)
    throw affineExcept("affineEnsemble","affineEnsemble",
		       "Must have at least 2 nodes",1);

  rank = MPI::COMM_WORLD.Get_rank();  
  if (rank == 0) {
    //Reserve room for queues
    procqueue.setCapacity(nproc);
    stepqueue.setCapacity(nwalkers/2+1);  //Only store one half at a time

    //Master node; slaves don't need this
    naccept.resize(nwalkers);
    naccept.assign(nwalkers, 0);
    ignore_params.resize(nparams);
    ignore_params.assign(nparams, false);

    //Set random number generator seed from time
    unsigned long long int seed;
    seed = static_cast<unsigned long long int>(time(NULL));
    seed += static_cast<unsigned long long int>(clock());
    rangen.setSeed(seed);
  }
    
}

affineEnsemble::~affineEnsemble() {}

/*!
  Only meaningful for the master node, always return true for slaves
 */
bool affineEnsemble::isValid() const {
  if (rank != 0) return true;
  if (nwalkers < 2) return false; //Need at least 2!
  if (nparams == 0) return false;
  if (nsteps == 0) return false;
  if (min_burn < 10) return false;
  if (burn_multiple <= 1.0) return false;
  if (scalefac <= 0.0) return false;
  if (nignore == nparams) return false;
  return true;

}

double affineEnsemble::generateZ() const {
  //Don't check isValid or if master for speed
  //Generate a random number satisfying g \propto 1/sqrt(z)
  // for z \elem [1/a,a]
  //The inverse cumulative distribution is
  //  G^{-1}(x) = [(a-1)*x+1]^2/a
  //a is scalefac in the nomenclature of this class
  double val;
  val = (scalefac-1.0)*rangen.doub()+1.0;
  return val*val/scalefac;
}

/*!
  \param[in] n New number of walkers
  
  Does not preserve contents unless new value is the same as current one.
 */
void affineEnsemble::setNWalkers(unsigned int n) {
  if (nwalkers == n) return;
  if (n == 0)
    throw affineExcept("affineEnsemble", "setNWalkers",
		       "n must be positive", 1);

  chains.setNWalkers(n);
  
  unsigned int nsamples = nsteps * nwalkers;
  nsteps = static_cast<unsigned int>( NSAMPLES/static_cast<double>(n)
				      + 0.999999999999 );

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
		       "n must be positive", 1);

  chains.setNParams(n);
  nignore = 0;

  has_name.resize(n);
  has_name.assign(n, false);
  parnames.resize(n);

  if (rank == 0) {
    ignore_params.resize(n);
    ignore_params.assign(n, false);
  }

  nparams = n;
}

double affineEnsemble::getMaxLogLike() const {
  if (rank != 0) return std::numeric_limits<double>::quiet_NaN();
  return chains.getMaxLogLike();
}

void affineEnsemble::getMaxLogLikeParam( double& val, paramSet& p ) const {
  if (rank != 0) {
    val = std::numeric_limits<double>::quiet_NaN();
    p.setNParams(nparams);
    for (unsigned int i = 0; i < nparams; ++i)
      p[i] = std::numeric_limits<double>::quiet_NaN();
    return;
  }
  chains.getMaxLogLikeParam(val,p);
}

void affineEnsemble::attentionAllParams() {
  if (rank != 0) return;
  ignore_params.assign( ignore_params.size(), false );
  nignore = 0;
}

/*!
  \returns True if anything was actually done, False if the parameter
           was already ignored
 */
bool affineEnsemble::ignoreParam(unsigned int idx) {
  if (rank != 0) return false;
  if (ignore_params[idx]) return false;
  ignore_params[idx] = true;
  nignore += 1;
  return true;
}

/*!
  \returns True if anything was actually done, False if the parameter
           was already being paid attention to
 */
bool affineEnsemble::attentionParam(unsigned int idx) {
  if (rank != 0) return false;
  if (!ignore_params[idx]) return false;
  ignore_params[idx] = false;
  nignore -= 1;
  return true;
}

void affineEnsemble::setParamName(unsigned int idx, const std::string& parnm) {
  if (idx >= nparams) throw affineExcept("affineEnsemble","setParamName",
					 "Invalid index",1);
  has_any_names = true;
  has_name[idx] = true;
  parnames[idx] = parnm;
}

void affineEnsemble::unsetParamName(unsigned int idx) {
  if (idx >= nparams) throw affineExcept("affineEnsemble","unsetParamName",
					 "Invalid index",1);
  if (!has_name[idx]) return; //Wasn't set anyways
  has_name[idx] = false;
  has_any_names = has_name[idx];
  for (unsigned int i = 1; i < nparams; ++i)
    has_any_names |= has_name[i];
}

std::string affineEnsemble::getParamName(unsigned int idx) const {
  if (idx >= nparams) throw affineExcept("affineEnsemble","getParamName",
					 "Invalid index",1);
  if (!has_name[idx]) return std::string();
  return parnames[idx];
}

bool affineEnsemble::computeAcor() const {
  if (rank != 0) return false;
  if (!isValid())
    throw affineExcept("affineEnsemble","computeAcor",
		       "Sampler not in valid state",1);
  
  if (acor.size() < nparams) acor.resize(nparams);
  bool success;
  success = chains.getAcorVector(acor, ignore_params);
  if (success) acor_set = true;
  return success;
}

bool affineEnsemble::getAcor(std::vector<double>& retval) const {
  if (rank != 0) 
    throw affineExcept("affineEnsemble","getAcor","Don't call on slave",1);
  bool success;
  if (! acor_set ) success = computeAcor(); else success=true;
  if (! success ) return success;
  retval.resize(nparams);
  for (unsigned int i = 0; i < nparams; ++i)
    retval[i] = acor[i];
  return success;
}

void affineEnsemble::getAcceptanceFrac(std::vector<double>& retval) const {
  if (rank != 0) 
    throw affineExcept("affineEnsemble","getAcceptanceFrac",
		       "Don't call on slave",1);
  if (!isValid())
    throw affineExcept("affineEnsemble","computeAcceptanceFrac",
		       "Sampler not in valid state",1);
  
  retval.resize(nwalkers);
  for (unsigned int i = 0; i < nwalkers; ++i) {
    unsigned int ncurrsteps = chains.getNIters(i);
    if (ncurrsteps == 0)
      retval[i] = std::numeric_limits<double>::quiet_NaN();
    else 
      retval[i] = static_cast<double>(naccept[i]) /
	static_cast<double>(ncurrsteps);
  }
}

void affineEnsemble::printStatistics(double conflevel, std::ostream& os) const {
  if (rank != 0) 
    throw affineExcept("affineEnsemble","printStatistics",
		       "Don't call on slave",1);
  if (!isValid())
    throw affineExcept("affineEnsemble","printStatistics",
		       "Sampler not in valid state",1);
  std::string parname;
  double mn, var, low, up;
  for (unsigned int i = 0; i < nparams; ++i) {
    if (!has_name[i]) parname = "Unknown"; else
      parname = parnames[i];
    chains.getParamStats(i,mn,var,low,up,conflevel);
    os << "Parameter: " << parname << " Mean: " << mn << " Stdev: "
       << sqrt(var) << std::endl;
    os << "  lower limit: " << low << " upper limit: "
       << up << " (" << conflevel*100.0 << "% limit)" << std::endl;
  }

}

/*!
  Manages the sampling.  Should be called for both master and slave
  processes.
  //It proceeds as follows:
  //  1) First, possibly do init_burn steps in each walker, and discard
  //  2) Then do min_burn steps in each walker
  //  2) Compute the autocorrelation length
  //  4) Do enough additional steps so that burn_multiple*acor steps have
  //     been done
  //  5) Compute acor again and make sure we still satisfy that criterion;
  //      if not do more steps, if so discard the previous steps
  //  6) Then do nsteps additional steps
 */
void affineEnsemble::sample() {
  //This routine is where all of the parallel bits happen
  // Different things happen if this is the master node or a slave one;
  //  slave ones just compute likelihoods, the master node suggests new
  //  steps and decides whether to accept them or not

  if (rank == 0) masterSample(); else slaveSample();
}

/*!
  A quick and dirty routine to do a fixed number of steps.
  Assumes that initChains is already called.

  \param[in] nsteps  Number of steps to do in each walker
  \param[in] initsteps Number of initial steps to do then discard
 */
void affineEnsemble::doSteps(unsigned int nsteps, unsigned int initsteps) {
  if (rank == 0) {
    if (!isValid())
      throw affineExcept("affineEnsemble","doSteps",
			 "Calling with invalid model setup",1);
    int jnk;

    if (initsteps > 0) {
      chains.addChunk(initsteps);
      for (unsigned int i = 0; i < initsteps; ++i)
	doMasterStep();
      chains.clearPreserveLast();
      for (unsigned int i = 0; i < naccept.size(); ++i)
	naccept[i] = 0;
      chains.setSkipFirst();
    }

    chains.addChunk(nsteps);
    for (unsigned int i = 0; i < nsteps; ++i)
      doMasterStep();

    unsigned int nproc = MPI::COMM_WORLD.Get_size();
    for (unsigned int i = 1; i < nproc; ++i)
      MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,mcmc_affine::STOP);
  } else {
    slaveSample();
  }
}

void affineEnsemble::masterSample() {
  if (!isValid())
    throw affineExcept("affineEnsemble","masterSample",
		       "Calling with invalid model setup",1);
  
  //Initialize chains; should also reserve space and clear old
  initChains();

  //Do burn in
  doBurnIn();
  if (verbose)
    std::cout << "Burned in after: " << chains.getMinNIters() 
	      << " steps" << std::endl;

  for (unsigned int i = 0; i < nwalkers; ++i)
    naccept[i] = 1; //Guess the last one was accepted...

  //Then do extra steps
  if (verbose) std::cout << "Doing " << nsteps << " additional steps per"
			 << " walker, for " << nsteps*nwalkers
			 << " total steps" << std::endl;
  chains.addChunk(nsteps);
  for (unsigned int i = 0; i < nsteps; ++i)
    doMasterStep();  
  
  //Tell slaves we are done
  int jnk;
  unsigned int nproc = MPI::COMM_WORLD.Get_size();
  for (unsigned int i = 1; i < nproc; ++i)
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,mcmc_affine::STOP);

}

double affineEnsemble::getMaxAcor() const {
  if (! acor_set ) 
    throw affineExcept("affineEnsemble","getMaxAcor",
		       "Called without acor available",1);
  double maxval;
  unsigned int start_idx;
  for (start_idx = 0; start_idx < nparams; ++start_idx)
    if ( ! ignore_params[start_idx] ) break;
  if (start_idx == nparams)
    throw affineExcept("affineEnsemble","getMaxAcor",
		       "All params ignored",2);
  maxval = acor[start_idx];
  for (unsigned int i = start_idx+1; i < nparams; ++i)
    if ( (!ignore_params[i]) && acor[i] > maxval ) maxval = acor[i];
  return maxval;
}

void affineEnsemble::doBurnIn() throw(affineExcept) {
  if (init_steps > 0) {
    chains.addChunk(init_steps);
    for (unsigned int i = 0; i < init_steps; ++i)
      doMasterStep();
    //Throw away burn in, keeping last step as first step of new one
    chains.clearPreserveLast(); 
    for (unsigned int i = 0; i < naccept.size(); ++i)
      naccept[i] = 0;
  }

  //First do min_burn steps in each walker
  chains.addChunk(min_burn);
  for (unsigned int i = 0; i < min_burn; ++i)
    doMasterStep();

  //See if we are burned in
  bool acor_success = computeAcor();
  
  //If it failed, we need more steps.  Add 50% of min burn
  // up to 50 times
  const unsigned int max_burn_iter = 50;
  unsigned int nextra = 10;
  if (! acor_success) {
    nextra = static_cast<unsigned int>(min_burn * 0.5);
    if (nextra < 10) nextra = 10;
    for (unsigned int i = 0; i < max_burn_iter; ++i) {
      if (verbose)
	std::cout << "Failed to compute acor after " << chains.getMinNIters()
		  << " steps.  Adding " << nextra << " more" << std::endl;
      chains.addChunk(nextra);
      for (unsigned int i = 0; i < nextra; ++i)
	doMasterStep();
      acor_success = computeAcor();
      if (acor_success) break;
    }
    if (!acor_success)
      throw affineExcept("affineEnsemble","doBurnIn",
			 "Can't compute acor; increase min_burn",1);
  }
  
  //Okay, we have an acor estimate of some sort
  double max_acor = getMaxAcor();
  if (verbose)
    std::cout << "After " << chains.getMinNIters() 
	      << " steps, maximum autocorrelation is: "
	      << max_acor << std::endl;
  
  unsigned int nsteps = min_burn;
  unsigned int nminsteps = 
    static_cast<unsigned int>(burn_multiple*max_acor + 0.999999999999999);
  bool burned_in = (nsteps > nminsteps);
  while (! burned_in ) {
    //Figure out how many more steps to do.  Do half of the number
    // estimated
    unsigned int nmore = (nminsteps-nsteps)/2 + 1;
    if (verbose) std::cout << "Doing " << nmore << " additional steps"
			   << std::endl;
    chains.addChunk(nmore);
    for (unsigned int i = 0; i < nmore; ++i)
      doMasterStep();
    nsteps += nmore;

    //Update acor, same as before
    acor_success = computeAcor();

    //Again, add more steps if we must
    nextra = static_cast<unsigned int>(min_burn * 0.5);
    if (nextra < 10) nextra = 10;
    for (unsigned int i = 0; i < max_burn_iter; ++i) {
      if (verbose)
	std::cout << "Failed to compute acor after " << chains.getMinNIters()
		  << " steps.  Adding " << nextra << " more" << std::endl;
      chains.addChunk(nextra);
      for (unsigned int i = 0; i < nextra; ++i)
	doMasterStep();
      acor_success = computeAcor();
      if (acor_success) break;
    }
    if (!acor_success)
      throw affineExcept("affineEnsemble","doBurnIn",
			 "Can't compute acor; increase min_burn",1);

    max_acor = getMaxAcor();
    nminsteps = 
      static_cast<unsigned int>(burn_multiple*max_acor + 0.999999999999999);
    burned_in = (nsteps > nminsteps);
  }

  //Throw away burn in, keeping last step as first step of new one
  chains.clearPreserveLast(); 
  for (unsigned int i = 0; i < naccept.size(); ++i)
    naccept[i] = 0;
  chains.setSkipFirst();
}

void affineEnsemble::generateNewStep( unsigned int idx1, unsigned int idx2,
				      proposedStep& prstep ) const {
  //Assume prstep is the right size, don't check for speed
  double zval = generateZ();
  //Get previous step for this walker
  chains.generateNewStep(zval,idx1,idx2,prstep.oldStep,
			 prstep.oldLogLike,prstep.newStep);
  prstep.update_idx = idx1;
  prstep.newLogLike=std::numeric_limits<double>::quiet_NaN();
}

void affineEnsemble::doMasterStep() throw (affineExcept) {
  //We use a pushdown stack to keep track of which things to do
  //To make this parallel, we must split the number of walkers in
  // half, and update each half from the opposite one
  //Fill in queue
  if (rank != 0)
    throw affineExcept("affineEnsemble","doMasterStep",
		       "Should only be called from master node",1);
  if (! stepqueue.empty() )
    throw affineExcept("affineEnsemble","doMasterStep",
		       "step queue should be empty at start",4);

  unsigned int minidx, maxidx;
  std::pair<unsigned int, unsigned int> pr;
  //Do the first half
  minidx = nwalkers/2; //Generate from second half
  maxidx = nwalkers;
  for (unsigned int i = 0; i < nwalkers/2; ++i) {
    pr.first  = i;
    pr.second = rangen.selectFromRange(minidx,maxidx);
    stepqueue.push( pr );
  }

  //Now run that
  emptyMasterQueue();

  //Then the second
  minidx = 0;
  maxidx = nwalkers/2;
  for (unsigned int i = nwalkers/2; i < nwalkers; ++i) {
    pr.first  = i;
    pr.second = rangen.selectFromRange(minidx,maxidx);
    stepqueue.push( pr );
  }
  emptyMasterQueue();

}

void affineEnsemble::emptyMasterQueue() throw (affineExcept) {
  //Step loop
  MPI::Status Info;
  int jnk, this_rank, this_tag;
  unsigned int ndone, nproc, nsteps;
  std::pair<unsigned int, unsigned int> pr;
  double prob;
  
  nproc = MPI::COMM_WORLD.Get_size();
  nsteps = stepqueue.size();
  ndone = 0;

  if (nsteps == 0) return; //Nothing to do...

  while (ndone < nsteps) {
    //First, if there are available procs, send them a step if we can
    if (!procqueue.empty()) 
      for (unsigned int i = 0; i < stepqueue.size(); ++i) {
	if (procqueue.empty()) break; //No more available procs
	this_rank = procqueue.pop();

	//Figure out next thing to update
	pr = stepqueue.pop();
	//Generate actual value into pstep
	generateNewStep( pr.first, pr.second, pstep );

	MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,this_rank,
			     mcmc_affine::SENDINGPOINT);
	pstep.sendSelf(MPI::COMM_WORLD,this_rank);
      }

    //No available procs, so wait for some sort of message
    MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,MPI::ANY_SOURCE,
			 MPI::ANY_SOURCE,Info);
    this_tag = Info.Get_tag();
    this_rank = Info.Get_source();
    if ( this_tag == mcmc_affine::ERROR ) {
      //Slave had a problem
      //Send stop message to all slaves
      for (unsigned int i = 1; i < nproc; ++i)
	MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,mcmc_affine::STOP);
      throw affineExcept("affineEnsemble","emptyMasterQueue",
			 "Problem encountered in slave",1);
    } else if ( this_tag == mcmc_affine::SENDINGRESULT ) {
      //Slave finished a computation, is sending us the result
      //Note we wait for a REQUESTPOINT to actually add it to the
      // list of available procs, which makes initialization easier
      pstep.recieveCopy(MPI::COMM_WORLD,this_rank);
      
      //Make sure it's okay
      if ( std::isinf( pstep.newLogLike ) ||
	   std::isnan( pstep.newLogLike ) ) {
	//That's not good -- exit
	for (unsigned int i = 1; i < nproc; ++i)
	  MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,mcmc_affine::STOP);
	throw affineExcept("affineEnsemble","emptyMasterQueue",
			   "Problem encountered in slave",2);
      }

      //Decide whether to accept the step or reject it, add to
      // chunk
      if (pstep.newLogLike > pstep.oldLogLike) {
	//Accept
	chains.addNewStep( pstep.update_idx, pstep.newStep,
			   pstep.newLogLike );
	naccept[pstep.update_idx] += 1;
      } else {
	prob = exp( pstep.newLogLike - pstep.oldLogLike );
	if (rangen.doub() < prob) {
	  //Also accept
	  chains.addNewStep( pstep.update_idx, pstep.newStep,
			     pstep.newLogLike );
	  naccept[pstep.update_idx] += 1;
	} else {
	  //reject, keep old step
	  chains.addNewStep( pstep.update_idx, pstep.oldStep,
			     pstep.oldLogLike );
	}
      }
      ndone += 1;
    } else if ( this_tag == mcmc_affine::REQUESTPOINT ) {
      procqueue.push(this_rank);
    } else {
      std::stringstream sstream;
      sstream << "Master got unexpected message from " << this_rank
	      << " with code: " << this_tag;
      //Tell slaves to stop
      for (unsigned int i = 1; i < nproc; ++i)
	MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,i,mcmc_affine::STOP);
      throw affineExcept("affineEnsemble","emptyMasterQueue",
			 sstream.str(),4);
    }
  }

}

void affineEnsemble::slaveSample() {
  //This one just sits here and computes likelihoods
  if (rank == 0)
    throw affineExcept("affineEnsemble","slaveSample",
		       "Should not be called from master node",1);
  
  initChains(); //Slave version

  int jnk;
  MPI::Status Info;

  try {
    while (true) { //Runs until we get a STOP message
      //Ask for a new point
      MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,mcmc_affine::REQUESTPOINT);

      //Now wait for one
      MPI::COMM_WORLD.Recv(&jnk,1,MPI::INT,0,MPI::ANY_TAG,Info);

      int tag = Info.Get_tag();
      if ( tag == mcmc_affine::STOP ) return; //Master says -- we're done!
      else if ( tag == mcmc_affine::SENDINGPOINT ) {
	pstep.recieveCopy(MPI::COMM_WORLD,0);
	
	//Compute likelihood
	pstep.newLogLike = getLogLike(pstep.newStep);

	//Send it back
	MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,mcmc_affine::SENDINGRESULT);
	pstep.sendSelf(MPI::COMM_WORLD,0);
      } else {
	std::cerr << "Unexpected message from master: "
                  << tag << " in slave: " << rank << std::endl;
	MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,mcmc_affine::ERROR);
	return;
      }
    }
  } catch ( const affineExcept& ex ) {
    std::cerr << "Error encountered for process: " << rank << std::endl;
    std::cerr << ex << std::endl;
    MPI::COMM_WORLD.Send(&jnk,1,MPI::INT,0,mcmc_affine::ERROR);
    return;
  }
}

/*!
  Doesn't output the first chunk (which should just have the initial step
  in it).
 */
void affineEnsemble::writeToFile( const std::string& filename ) const {
  if (rank != 0) 
    throw affineExcept("affineEnsemble","writeToFile",
		       "Should only be called from master node",1);  
  chains.writeToFile(filename);
}
