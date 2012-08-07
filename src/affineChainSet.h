//affineChainSet.h

#ifndef __affineChainSet__
#define __affineChainSet__

#include<paramSet.h>
#include<affineExcept.h>

/*!
  \brief Set of params and likelihoods
*/
struct affineStepChunk {
  unsigned int nwalkers; //!< Number of walkers
  unsigned int niters; //!< Number of steps in this chunk
  unsigned int nparams; //!< Number of params
  double* steps; //!< holds steps: nwalkers by niter by nparams
  double* logLike; //!< Log likelihoods, nwalkers by niters
  unsigned int *nsteps; //!< Number of filled elements for each walker (size nwalkers)

  affineStepChunk(unsigned int, unsigned int, unsigned int); //!< Constructor
  affineStepChunk(const affineStepChunk&); //!< Copy constructor
  ~affineStepChunk(); //!< Destructor

  /*! \brief Fills a chunk with the last step of another chunk, resizing as
  needed */
  void fillFromLastStep(const affineStepChunk&) throw (affineExcept);

  void clear(); //!< Frees memory

  /*! \brief Gets number of steps for specified walker */
  unsigned int getNSteps(unsigned int i) const { return nsteps[i]; }
  unsigned int getMinNSteps() const; //!< Get minimum number of steps accross all walkers

  affineStepChunk& operator=(const affineStepChunk&); //!< Copy
 
  double getMaxLogLike() const; //!< Return maximum log likelihood in this chunk
  void getMaxLogLikeParam(double&, paramSet&) const; //!< Return best Log Likelihood step and it's likelihood in this chunk

  /*! \brief Adds a new step, values in paramSet */
  bool addStep(unsigned int, const paramSet&, double);

  /*! \brief Get specified step */
  bool getStep(unsigned int, unsigned int, paramSet&, double&) const;

  /*! \brief get most resent step for given walker */
  bool getLastStep(unsigned int, paramSet&, double& ) const;
  
  /*! \brief Get pointer to particular set of parameters */
  double* getParamPointer(unsigned int, unsigned int); 
  /*! \brief Get pointer to particular set of parameters */
  double* const getParamPointer(unsigned int, unsigned int) const;
  /*! \brief Get pointer to last set of parameters for a given walker*/
  double* getLastParamPointer(unsigned int);
  /*! \brief Get pointer to last set of parameters for a given walker*/
  double* const getLastParamPointer(unsigned int) const;

  /*! \brief Get log likelihood for a particular walker and step */
  double getLogLike(unsigned int i1, unsigned int i2) const {
    return logLike[i1*niters+i2]; }

  void writeToStream( std::ostream& os ) const; //!< Write to stream
};

std::ostream& operator<<( std::ostream& os, const affineStepChunk& );

/*!
  \brief Holds chains
*/
class affineChainSet {
 private:
  static const unsigned int taumax;  //!< Acor computation constant
  static const unsigned int winmult; //!< Acor computation constant
  static const unsigned int maxlag;  //!< Acor computation constant
  static const unsigned int minfac;  //!< Acor computation constant

  //The ideal storage scheme here would be a nwalker x niter x nparams
  // array.  The problem is that we don't know in advance how
  // many iterations we will do...
  //So we try to approximate this by keeping a bunch of affineStepChunks
  // The idea is that the user will generally be doing big chunks of steps, 
  // rather than one at a time.  We are also taking advantage of the fact
  // that, due to the way affine MCMC works, the walkers are updated
  // sequentially, and hence all have about the same number of steps
  unsigned int nwalkers; //!< Number of walkers
  unsigned int nparams; //!< Number of params
  std::vector< affineStepChunk* > steps; //!< Holds pointers to nwalker by niter by nparams arrays

  /*! \brief Autocorrelation engine */
  int acor(double&,double&,double&,std::vector<double>&,
	   unsigned int) const throw (affineExcept);

  bool skipfirst; //!< True if the first chunk is just the initialization, and should be ignored in stats, etc.

 public:
  affineChainSet(unsigned int, unsigned int); //!< Constructor
  ~affineChainSet(); //!< Destructor

  void clear(); //!< Clears the chain
  void clearPreserveLast() throw (affineExcept); //!< Clears the chain, preserving the last step of the old one as the first of the new

  unsigned int getNWalkers() const { return nwalkers; } //!< Get number of walkers
  void setNWalkers(unsigned int) throw (affineExcept); //!< Change the number of walkers

  unsigned int getNParams() const { return nparams; } //!< Get number of params
  void setNParams(unsigned int) throw (affineExcept); //!< Change the number of parameters
  unsigned int getNIters() const; //!< Total number of iterations present across all walkers
  unsigned int getNIters(unsigned int) const throw (affineExcept); //!< Total number of iterations for given walker
  unsigned int getMinNIters() const; //!< Minimum number of iterations across all walkers

  unsigned int getNChunks() const { return steps.size(); } //!< Get number of chunks

  bool doSkipFirst() const { return skipfirst; } //!< Are we skipping the first tep in any results?
  void setSkipFirst() { skipfirst = true; } //!< Do skip the first step
  void unsetSkipFirst() { skipfirst = false; } //!< Don't skip the first step

  double getMaxLogLike() const; //!< Return maximum log likelihood over all walkers
  void getMaxLogLikeParam(double&,paramSet&) const; //!< Return best Log Likelihood step and it's likelihood over all walkers

  void addChunk(unsigned int); //!< Adds a new chunk of specified size
  bool addNewStep(unsigned int, const paramSet&, double); //!< Adds a new step

  /*! \brief Create a new step */
  void generateNewStep(double,unsigned int, unsigned int, paramSet&,
		       double&, paramSet&) const throw (affineExcept);

  /*! \brief Returns the most recent step for the specified walker */
  void getLastStep( unsigned int, paramSet&, double& ) const
    throw (affineExcept);

  /*! \brief Get a specified step */
  void getStep(unsigned int, unsigned int, unsigned int, paramSet&,
	       double&) const throw (affineExcept);

  /*! \brief Fills vector with param info */
  void getParamVector(unsigned int, unsigned int,
		      std::vector<double>& pvec) const throw (affineExcept);
  /*! \brief Fills vector with param info averaged over walkers*/
  void getAverageParamVector(unsigned int, std::vector<double>&) const 
    throw (affineExcept);

  /*! \brief Autocorrelation length for a particular parameter */
  double getAcor(unsigned int, double&, double&, 
		 bool&) const throw (affineExcept);
  bool getAcorVector(std::vector<double>&) const throw (affineExcept); //!< Autocorrelation length for all params
  bool getAcorVector(std::vector<double>&, const std::vector<bool> ignore) 
    const throw(affineExcept); //!< Autocorrelation length for all params
  
  affineChainSet& operator=(const affineChainSet&); //!< Replace this chain set with a copy of another
  affineChainSet& operator+=( const affineChainSet& ) throw (affineExcept);  //!<Merge another chain set onto this

  double getParamMean(unsigned int) const throw (affineExcept); //!< Get mean value of param across all walkers
  /* \brief Get statistics on a given parameter across all walkers */
  void getParamStats(unsigned int, double&, double&, double&,
		     double&, double=0.683) const
    throw (affineExcept); 
  void makeCovMatrix(double** covmat) const; //!< Get covariance of all params

  void writeToFile( const std::string&) const throw (affineExcept); //!< Write to file

};

#endif
