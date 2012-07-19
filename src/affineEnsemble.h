//affineEnsemble.h

#ifndef __affineEnsemble__
#define __affineEnsemble__

#include<utility>

#include<affineChainSet.h>
#include<affineQueue.h>
#include<proposedStep.h>
#include<paramSet.h>
#include<ran.h>

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

  double scalefac; //!< Scaling parameter for proposal density

  unsigned int nignore; //!< Number of params to ignore in burn-in computation
  std::vector<bool> ignore_params; //!< Which params to ignore in burn-in

  unsigned int init_steps; //!< Number of initial steps to do before starting burn-in check
  unsigned int min_burn; //!< Minimum number of burn steps to do
  unsigned int burn_multiple; //!< Multiple of autocorrelation steps to do for burn
  unsigned int nsteps; //!< Number of steps to do after burn per walker

  mutable std::vector<double> acor; //!< Holds autocorrelation
  mutable bool acor_set; //!< Has acor been computed?

  /*! \brief For keeping track of available procs */
  affineQueue<int> procqueue;

  /*! \brief For keeping track of what to do next */
  affineQueue< std::pair< unsigned int, unsigned int > > stepqueue;
  mutable proposedStep pstep; //!< Convenience variable for new steps

  double getMaxAcor() const; //!< Return maximum autocorrelation length

  //Primary sampling routines
  void masterSample(); //!< Master node sampler routine
  void slaveSample(); //!< Slave node sampler routine

  //sampling sub-routines
  void doBurnIn() throw (affineExcept); //!< Master node burn in routine
  void doMasterStep() throw (affineExcept); //!< Does a step for all walkers, master node
  void emptyMasterQueue() throw (affineExcept); //!< Runs all steps in stepqueue as master node

  /*! \brief Generate a new proposed Step */
  void generateNewStep( unsigned int, unsigned int, proposedStep&) const;


protected:

  affineChainSet chains; //!< Holds actual steps
  std::vector<unsigned int> naccept; //!< Number of accepted steps

  mutable ran rangen; //!< Random number generator

  unsigned int rank; //!< Which node is this; if 0 master, otherwise slave

  bool verbose; //!< Ouput information messages as we run

public:
  /* \brief Constructor */
  affineEnsemble(unsigned int, unsigned int, unsigned int,
		 unsigned int=50, unsigned int=50, double=5, double=2);
  ~affineEnsemble(); //!< Destructor

  bool isValid() const; //!< Are params valid?

  void setVerbose() { verbose = true; } //!< Set verbose mode
  void unsetVerbose() { verbose = false; } //!< Turn off verbose mode

  unsigned int getNWalkers() const { return nwalkers; } //!< Get number of walkers
  unsigned int getNParams() const { return nparams; } //!< Get number of params
  unsigned int getNChunks() const { return chains.getNChunks(); } //!< Get number of chunks
  unsigned int getMinNIters() const { return chains.getMinNIters(); } //!< Get minimum number of iterations across all walkers
  unsigned int getNSteps() const { return nsteps; } //!< Get number of steps
  double getScalefac() const { return scalefac; } //!< Get scalefac
  void setScalefac(double val) { scalefac = val; } //!< Return scalefac

  double getMaxLogLike() const; //!< Get Maximum recorded Log Likelihood
  void getMaxLogLikeParam(double& val, paramSet& p) const; //!< Get parameters corresponding to maximum recorded log likelihood

  //Dealing with ignorable params
  unsigned int getNIgnoreParams() const { return nignore; } //!< Return number of parameters being ignored for burn-in check
  void attentionAllParams(); //!< Don't ignore any parameters in burn-in check
  bool ignoreParam( unsigned int ); //!< Ignore this parameter in burn-in check
  bool attentionParam( unsigned int ); //!< Don't ignore this parameter is burn-in check
  /*! \brief Is a particular parameter being ignored in burn-in */
  bool isParamIgnored(unsigned int idx) const { return ignore_params[idx]; }
  
  void setParamName(unsigned int, const std::string&); //!< Set param names
  void unsetParamName(unsigned int); //!< Unset parameter name
  bool hasName(unsigned int idx) const { return has_name[idx]; } //!< Does a particular parameter have a name?
  std::string getParamName(unsigned int) const; //!< Get param name for particular parameter

  /*! \brief Sets random number generator seed */
  void setSeed( unsigned long long int seed ) const { rangen.setSeed(seed); }
  double generateZ() const; //!< Generate a Z value

  bool computeAcor() const; //!< Computes autocorrelation

  bool getAcor(std::vector<double>&) const; //!< Returns acor, or computes it if not set
  void getAcceptanceFrac(std::vector<double>&) const; //!< Returns acceptance fraction for each walker

  void printStatistics(double=0.683, std::ostream& = std::cout) const; //!< Output statistics for run

  //User must subclass these for their likelihood function
  virtual void initChains() = 0; //!< Initialize first step in chains
  virtual double getLogLike(const paramSet&) = 0; //!< Computes log likelihood

  void sample(); //!< Get samples
  
  void doSteps(unsigned int,unsigned int=0); //!< Do a fixed number of steps

  void writeToFile( const std::string& ) const; //!< Write results to file

};

#endif
