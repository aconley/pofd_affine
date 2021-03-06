//affineChainSet.h

#ifndef __affineChainSet__
#define __affineChainSet__

#include "hdf5.h"

#include "paramSet.h"
#include "affineExcept.h"
#include "global_settings.h"

/*!
  \brief Set of params and likelihoods

  This stores a chunk of steps
*/
struct affineStepChunk {
  // The parameters are kept as floats, but the likelihoods as doubles
  // because they may have large constant factors (since normalization
  // is hard).
  unsigned int nwalkers; //!< Number of walkers
  unsigned int niters; //!< Number of steps in this chunk
  unsigned int nparams; //!< Number of params
  float* steps; //!< holds steps: nwalkers by niter by nparams
  double* logLike; //!< Log likelihoods, nwalkers by niters
  unsigned int *nsteps; //!< Number of filled elements for each walker (size nwalkers)

  affineStepChunk(unsigned int, unsigned int, unsigned int); //!< Constructor
  affineStepChunk(const affineStepChunk&); //!< Copy constructor
  affineStepChunk(affineStepChunk&&)=delete;
  ~affineStepChunk(); //!< Destructor

  /*! \brief Fills a chunk with the last step of another chunk, resizing as
  needed */
  void fillFromLastStep(const affineStepChunk&);

  void clear(); //!< Frees memory

  /*! \brief Gets number of steps for specified walker */
  unsigned int getNSteps(unsigned int i) const noexcept { return nsteps[i]; }
  unsigned int getMinNSteps() const noexcept; //!< Get minimum number of steps accross all walkers

  affineStepChunk& operator=(const affineStepChunk&); //!< Copy

  double getMaxLogLike() const noexcept; //!< Return maximum log likelihood in this chunk
  void getMaxLogLikeParam(double&, paramSet&) const; //!< Return best Log Likelihood step and it's likelihood in this chunk

  /*! \brief Adds a new step, values in paramSet */
  bool addStep(unsigned int, const paramSet&, double);

  /*! \brief Get specified step */
  bool getStep(unsigned int, unsigned int, paramSet&, double&) const;

  /*! \brief Get most resent step for given walker */
  bool getLastStep(unsigned int, paramSet&, double&) const;
  
  /*! \brief Replace last step */
  bool replaceLastStep(unsigned int, const paramSet&, double);

  /*! \brief Get pointer to particular set of parameters */
  float* getParamPointer(unsigned int, unsigned int); 
  /*! \brief Get pointer to particular set of parameters */
  float* const getParamPointer(unsigned int, unsigned int) const;
  /*! \brief Get pointer to last set of parameters for a given walker*/
  float* getLastParamPointer(unsigned int);
  /*! \brief Get pointer to last set of parameters for a given walker*/
  float* const getLastParamPointer(unsigned int) const;

  /*! \brief Get log likelihood for a particular walker and step */
  double getLogLike(unsigned int i1, unsigned int i2) const {
    return logLike[i1*niters+i2]; }

  void writeToStream(std::ostream& os) const; //!< Write to stream
};

std::ostream& operator<<(std::ostream& os, const affineStepChunk&);

/*!
  \brief Holds chains

  Chains are comprised of a set of affineStepChunks.
  This also contains routines for computing the autocorrelation of
  the chains, as well as serialization routines.
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
  int acor(double&, double&, double&, std::vector<float>&,
           unsigned int) const throw (affineExcept);

  bool skipfirst; //!< True if the first chunk is just the initialization, and should be ignored in stats, etc.

 public:
  affineChainSet(unsigned int, unsigned int); //!< Constructor
  affineChainSet(const affineChainSet&)=delete;
  affineChainSet(affineChainSet&&)=delete;
  ~affineChainSet(); //!< Destructor

  void clear(); //!< Clears the chain
  void clearPreserveLast() throw (affineExcept); //!< Clears the chain, preserving the last step of the old one as the first of the new

  unsigned int getNWalkers() const noexcept { return nwalkers; } //!< Get number of walkers
  void setNWalkers(unsigned int) throw (affineExcept); //!< Change the number of walkers

  unsigned int getNParams() const noexcept { return nparams; } //!< Get number of params
  void setNParams(unsigned int) throw (affineExcept); //!< Change the number of parameters
  unsigned int getNIters() const noexcept; //!< Total number of iterations present across all walkers
  unsigned int getNIters(unsigned int) const throw (affineExcept); //!< Total number of iterations for given walker
  unsigned int getMinNIters() const noexcept; //!< Minimum number of iterations across all walkers

  unsigned int getNChunks() const noexcept { return steps.size(); } //!< Get number of chunks

  bool doSkipFirst() const noexcept { return skipfirst; } //!< Are we skipping the first tep in any results?
  void setSkipFirst() noexcept { skipfirst = true; } //!< Do skip the first step
  void unsetSkipFirst() noexcept { skipfirst = false; } //!< Don't skip the first step

  double getMaxLogLike() const noexcept; //!< Return maximum log likelihood over all walkers
  void getMaxLogLikeParam(double&, paramSet&) const; //!< Return best Log Likelihood step and it's likelihood over all walkers

  void addChunk(unsigned int); //!< Adds a new chunk of specified size
  bool addNewStep(unsigned int, const paramSet&, double); //!< Adds a new step

  /*! \brief Returns the most recent step for the specified walker */
  void getLastStep(unsigned int, paramSet&, double&) const
    throw (affineExcept);

  /*! \brief Replace last step */
  void replaceLastStep(unsigned int, const paramSet&, double);

  /*! \brief Get a specified step */
  void getStep(unsigned int, unsigned int, unsigned int, paramSet&,
               double&) const throw (affineExcept);

  /*! \brief Fills vector with param info */
  void getParamVector(unsigned int, unsigned int,
                      std::vector<float>& pvec) const throw (affineExcept);
  /*! \brief Fills vector with param info averaged over walkers*/
  void getAverageParamVector(unsigned int, std::vector<float>&) const 
    throw (affineExcept);

  /*! \brief Autocorrelation length for a particular parameter */
  double getAcor(unsigned int, double&, double&, 
                 bool&) const throw (affineExcept);
  bool getAcorVector(std::vector<float>&) const throw (affineExcept); //!< Autocorrelation length for all params
  bool getAcorVector(std::vector<float>&, 
                     const std::vector<int> param_state) 
    const throw(affineExcept); //!< Autocorrelation length for all params
  
  affineChainSet& operator=(const affineChainSet&); //!< Replace this chain set with a copy of another
  affineChainSet& operator+=(const affineChainSet&) throw (affineExcept);  //!<Merge another chain set onto this

  float getParamMean(unsigned int) const throw (affineExcept); //!< Get mean value of param across all walkers
  /*! \brief Get statistics on a given parameter across all walkers */
  void getParamStats(unsigned int, float&, float&, float&,
                     float&, float=0.683) const
    throw (affineExcept); 
  void makeCovMatrix(float** covmat) const; //!< Get covariance of all params

  void writeToFile(const std::string&) const throw (affineExcept); //!< Write to file

  void writeToHDF5Handle(hid_t) const; //!< Write to HDF5 handle
  void writeToHDF5(const std::string&) const; //!< Write to HDF5 file
};

#endif
