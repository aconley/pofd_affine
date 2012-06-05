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

  affineStepChunk(unsigned int, unsigned int, unsigned int);
  affineStepChunk(const affineStepChunk&);
  ~affineStepChunk();

  void fillFromLastStep(const affineStepChunk&) throw (affineExcept);

  void clear(); //!< Frees memory

  unsigned int getNSteps(unsigned int i) const { return nsteps[i]; }
  unsigned int getMinNSteps() const;

  affineStepChunk& operator=(const affineStepChunk&); //!< Copy
 
  double getMaxLogLike() const; //!< Return maximum log likelihood in this chunk
  void getMaxLogLikeParam(double&, paramSet&) const; //!< Return best Log Likelihood step and it's likelihood in this chunk

  /*! \brief Adds a new step, values in paramSet */
  bool addStep(unsigned int, const paramSet&, double);

  /*! \brief Get specified step */
  bool getStep(unsigned int, unsigned int, paramSet&, double&) const;

  /*! \brief get most resent step for given walker */
  bool getLastStep(unsigned int, paramSet&, double& ) const;

  double* getParamPointer(unsigned int, unsigned int);
  double* const getParamPointer(unsigned int, unsigned int) const;
  double* getLastParamPointer(unsigned int);
  double* const getLastParamPointer(unsigned int) const;

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
  //Acor computation parameters
  static const unsigned int taumax, winmult, maxlag, minfac;

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
  affineChainSet(unsigned int, unsigned int);
  ~affineChainSet();

  void clear(); //!< Clears the chain
  void clearPreserveLast() throw (affineExcept); //!< Clears the chain, preserving the last step of the old one as the first of the new

  unsigned int getNWalkers() const { return nwalkers; }
  unsigned int getNIters() const; //!< Total number of iterations present across all walkers
  unsigned int getNIters(unsigned int) const throw (affineExcept); //!< Total number of iterations for given walker
  unsigned int getMinNIters() const; //!< Minimum number of iterations across all walkers
  unsigned int getNParams() const { return nparams; }
  unsigned int getNChunks() const { return steps.size(); }

  bool doSkipFirst() const { return skipfirst; }
  void setSkipFirst() { skipfirst = true; }
  void unsetSkipFirst() { skipfirst = false; }

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
  void getParamStats(unsigned int, double&, double&, double&,
		     double&, double=0.683) const
    throw (affineExcept); //!< Get statistics on a given parameter across all walkers
  void makeCovMatrix(double** covmat) const; //!< Get covariance of all params

  void writeToFile( const std::string&) const throw (affineExcept); //!< Write to file

};

#endif
