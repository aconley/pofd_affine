//affineEnsemble.h

#ifndef __affineEnsemble__
#define __affineEnsemble__

#include<utility>

#include "affineChainSet.h"
#include "affineQueue.h"
#include "proposedStep.h"
#include "paramSet.h"
#include "ran.h"
#include "global_settings.h"

/*!
  \brief Affine MCMC ensemble sampler.  See Foreman-Mackey et al. 2012
*/

class affineEnsemble {
private:
  unsigned int nwalkers; //!< Number of walkers
  unsigned int nparams; //!< Number of params
  
  //Names of parameters
  bool has_any_names; //!< Have we set any parameter names?
  std::vector<bool> has_name; //!< Does a given param have a name
  std::vector< std::string > parnames; //!< Optional names of params

  float scalefac; //!< Scaling parameter for proposal unsigned

  // Keep track of special parameter states
  unsigned int nfixed; //!< Number of fixed params
  unsigned int nignore; //!< Number ignored in acor
  std::vector<int> param_state; //!< Set to param_statue enum values

  unsigned int init_steps; //!< Number of initial steps to do before starting burn-in check
  double init_temp; //!< Temperature used during initial steps
  unsigned int min_burn; //!< Minimum number of burn steps to do
  bool fixed_burn; //!< Do fixed burn in rather than using autocorrelation
  float burn_multiple; //!< Multiple of autocorrelation steps to do for burn
  unsigned int nsteps; //!< Number of steps to do after burn per walker

  mutable std::vector<float> acor; //!< Holds autocorrelation
  mutable bool acor_set; //!< Has acor been computed?

  /*! \brief For keeping track of available procs */
  affineQueue<int> procqueue;
  /*! \brief For keeping track of what to do next */
  affineQueue< std::pair<int, int> > stepqueue;
  mutable proposedStep pstep; //!< Convenience variable for new steps

  mutable paramSet params_tmp; //!< Temporary internal step storage

  float getMaxAcor() const; //!< Return maximum autocorrelation length

  //Primary sampling routines
  void masterSample(); //!< Master node sampler routine
  void slaveSample(); //!< Slave node sampler routine

  //sampling sub-routines
  void doBurnIn() throw (affineExcept); //!< Master node burn in routine
  void doMasterStep(double=1.0) throw (affineExcept); //!< Does a step for all walkers, master node
  void emptyMasterQueue(double=1.0) throw (affineExcept); //!< Runs all steps in stepqueue as master node
  void calcLastLikelihood(); //!< Compute the likelihoods of the last step

protected:
  bool is_init; //!< Have the chains been initialized?

  bool has_initStep; //!< Is initial position recorded?
  paramSet initStep; //!< Initial position

  bool has_regenFirstStep; //!< Has post initial steps regenerated postion
  paramSet regenFirstStep; //!< Regenerated initial position (after init_steps)

  affineChainSet chains; //!< Holds actual steps
  std::vector<unsigned int> naccept; //!< Number of accepted steps

  mutable ran rangen; //!< Random number generator

  int rank; //!< Which node is this; if 0 master, otherwise slave

  // 0 is non-verbose, 1 is verbose, 2 is ultra-verbose, 3 is ultra-ultra, etc.
  unsigned int verbosity; //!< Verbosity level

public:
  /* \brief Constructor */
  affineEnsemble(unsigned int, unsigned int, unsigned int,
		 unsigned int=0, double=2.0, unsigned int=50, 
		 bool=false, float=5, float=2);
  ~affineEnsemble(); //!< Destructor

  bool isValid() const; //!< Are params valid?

  void setQuiet() { setVerbosity(0); } //!< Turn off verbose mode
  void setVerbose() { setVerbosity(1); } //!< Set verbose mode
  void setUltraVerbose() { setVerbosity(2); } //!< Set ultra-verbose mode
  virtual void setVerbosity(unsigned int V) { verbosity = V; } //!< Set verbosity level
  unsigned int getVerbosity() const { return verbosity; } //!< Get verbosity level

  unsigned int getNWalkers() const { return nwalkers; } //!< Get number of walkers
  void setNWalkers(unsigned int); //!< Set the number of walkers

  unsigned int getNParams() const { return nparams; } //!< Get number of params
  void setNParams(unsigned int); //!< Set the number of parameters

  unsigned int getNChunks() const { return chains.getNChunks(); } //!< Get number of chunks
  unsigned int getMinNIters() const { return chains.getMinNIters(); } //!< Get minimum number of iterations across all walkers
  unsigned int getNSteps() const { return nsteps; } //!< Get number of steps
  float getScalefac() const { return scalefac; } //!< Get scalefac
  void setScalefac(float val) { scalefac = val; } //!< Set scalefac

  double getMaxLogLike() const; //!< Get Maximum recorded Log Likelihood
  void getMaxLogLikeParam(double&, paramSet&) const; //!< Get parameters corresponding to maximum recorded log likelihood

  // Interact with parameter states
  // Fixing a parameter will also make it ignored
  void clearParamState(); //!< Fit all params and use all params in autocorrelation
  unsigned int getNFitParams() const; //!< Return number of params being fit
  unsigned int getNAcorParams() const; //!< Return number of parameters used in autocorrelation (burn in) check
  void fixParam(unsigned int); //!< Treat parameter as fixed
  bool isParamFixed(unsigned int) const; //!< Is a particular parameter fixed?
  void ignoreParamAcor(unsigned int); //!< Ignore parameter in acor computation
  bool isParamIgnoredAcor(unsigned int) const; //!< Is a parameter being ignored?

  // Parameter names
  void setParamName(unsigned int, const std::string&); //!< Set param names
  void unsetParamName(unsigned int); //!< Unset parameter name
  bool hasName(unsigned int idx) const { return has_name[idx]; } //!< Does a particular parameter have a name?
  std::string getParamName(unsigned int) const; //!< Get param name for particular parameter

  /*! \brief Sets random number generator seed */
  void setSeed(unsigned long long int seed) const { rangen.setSeed(seed); }
  float generateZ() const; //!< Generate a Z value

  bool computeAcor() const; //!< Computes autocorrelation

  bool getAcor(std::vector<float>&) const; //!< Returns acor, or computes it if not set
  bool hasOneStepBeenAccepted() const; //!< Has at least one step been accepted
  void getAcceptanceFrac(std::vector<float>&) const; //!< Returns acceptance fraction for each walker

  float getParamMean(unsigned int) const;
  void getParamStats(unsigned int, float&, float&, float&,
		     float&, float=0.683) const;
  virtual void printStatistics(float=0.683, std::ostream& = std::cout) const; //!< Output statistics for run

  //User must subclass these for their use.
  // Note that initChains should set up some sort of initial position as well
  virtual void initChains() = 0; //!< Set up information in each chain.  Must set is_init to true
  /*! \brief Generate initial position based on input parameter set

    The likelihoods of this initial step will be computed elsewhere.
    The initial parameters must be considered valid by the likelihood function.
  */
  virtual void generateInitialPosition(const paramSet&) = 0; //!< 

  /*! \brief Computes the log likelihood. */
  virtual double getLogLike(const paramSet&, bool& params_rejected) = 0;

  double getLogLike(const paramSet&); //!< Computes log likelihood

  /*! \brief Tests whether a given parameter set is valid */
  virtual bool areParamsValid(const paramSet&) const { return true; }

  /*! \brief Generate a new step */
  void generateNewStep(unsigned int, unsigned int, proposedStep& prstep) 
    const throw (affineExcept);

  void sample(); //!< Do burn in, get samples  
  void doSteps(unsigned int, unsigned int=0); //!< Do a fixed number of steps

  virtual void writeToStream(std::ostream&) const; //!< Write summary of fit parameters

  void writeToFile(const std::string&) const; //!< Write results to file
  virtual void writeToHDF5Handle(hid_t) const; //!< Write results to HDF5 handle
  void writeToHDF5(const std::string&) const; //!< Write results to HDF5 file
};

/*! \brief Write to stream */
std::ostream& operator<<(std::ostream& os, const affineEnsemble&);

#endif
