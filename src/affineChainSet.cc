#include<limits>
#include<cmath>
#include<fstream>
#include<algorithm>
#include<iomanip>
#include<sstream>

#include "../include/affineChainSet.h"

//Parameters for acor computation
const unsigned int affineChainSet::taumax = 2; 
//Autocovariance up to lag = winmult*tau
const unsigned int affineChainSet::winmult = 5; 
const unsigned int affineChainSet::maxlag = taumax * winmult;
const unsigned int affineChainSet::minfac = 5;

/*!
  \param[in] NWALKERS Number of walkers
  \param[in] NITERS Maximum number of steps allowed in this chunk
  \param[in] NPARAMS Number of parameters in this chunk
*/
affineStepChunk::affineStepChunk(unsigned int NWALKERS,
				 unsigned int NITERS,
				 unsigned int NPARAMS) {
  nwalkers = NWALKERS;
  niters = NITERS;
  nparams = NPARAMS;

  //Allocate in large contiguous block for memory access efficiency
  unsigned int sz = nwalkers * niters * nparams;
  if (sz > 0) {
    steps = new float[sz];
    logLike = new double[nwalkers * niters];
    nsteps = new unsigned int[nwalkers];
    for (unsigned int i = 0; i < nwalkers; ++i)
      nsteps[i] = 0;
  } else {
    steps = NULL;
    logLike = NULL;
    nsteps = NULL;
  }
}

/*!
  \param[in] other affineStepChunk to copy
*/
affineStepChunk::affineStepChunk(const affineStepChunk& other) {
  nwalkers = other.nwalkers;
  niters = other.niters;
  nparams = other.nparams;

  unsigned int sz = nwalkers * niters * nparams;
  if (sz > 0) {
    steps = new float[sz];
    for (unsigned int i = 0; i < sz; ++i)
      steps[i] = other.steps[i];
    logLike = new double[nwalkers * niters];
    for (unsigned int i = 0; i < nwalkers*niters; ++i) 
      logLike[i] = other.logLike[i];
    nsteps = new unsigned int[nwalkers];
    for (unsigned int i = 0; i < nwalkers; ++i)
      nsteps[i] = other.nsteps[i];
  } else {
    steps = NULL;
    logLike = NULL;
    nsteps = NULL;
  }
}


affineStepChunk::~affineStepChunk() {
  clear();
}

/*!
  \param[in] other Other chunk to fill from

  Takes this chunk, throws out the contents, makes it one iteration long,
  and sets the only entry for each walker to the last step from the other
  walker. 
*/
void affineStepChunk::fillFromLastStep(const affineStepChunk& other)
  throw (affineExcept) {
  unsigned int sz;
  sz = other.nwalkers * other.niters * other.nparams;
  if (sz == 0)
    throw affineExcept("affineStepChunk","fillFromLastStep",
		       "Other chunk has zero size",1);
  if (other.getMinNSteps() < 1)
    throw affineExcept("affineStepChunk","fillFromLastStep",
		       "Other chunk doesn't have steps for all walkers",2);
  //Resize as needed
  if (nwalkers != other.nwalkers || nparams != other.nparams || niters != 1) {
    clear();
    nwalkers = other.nwalkers;
    nparams = other.nparams;
    niters = 1;
    sz = nwalkers * niters * nparams;
    if (sz > 0) {
      steps = new float[sz];
      logLike = new double[nwalkers*niters];
      nsteps = new unsigned int[nwalkers];
      for (unsigned int i = 0; i < nwalkers; ++i)
	nsteps[i] = 0;
    } else {
      steps = NULL;
      logLike = NULL;
      nsteps = NULL;
      return;
    }
  } 

  //Now, for each walker, set the last param
  unsigned int laststep;
  float *ptro, *ptrn;
  for (unsigned int i = 0; i < nwalkers; ++i) {
    laststep = other.nsteps[i]-1;
    ptro = other.getParamPointer(i, laststep);
    ptrn = getParamPointer(i, 0);
    for (unsigned int j = 0; j < nparams; ++j)
      ptrn[j] = ptro[j];
    //Recall logLike is nwalkers * niters, and niters is one
    logLike[i] = other.getLogLike(i, laststep);
    nsteps[i] = 1;
  }
}

//Frees memory
void affineStepChunk::clear() {
  if (steps != NULL) {
    delete[] steps;
    steps = NULL;
  }
  if (logLike != NULL) {
    delete[] logLike;
    logLike = NULL;
  }
  if (nsteps != NULL) {
    delete[] nsteps;
    nsteps = NULL;
  }
  nwalkers = niters = nparams = 0;
}

/*!
  \returns Minimum number of steps taken by any walker
*/
unsigned int affineStepChunk::getMinNSteps() const {
  if (nwalkers == 0) return 0;
  unsigned int minstep = nsteps[0];
  unsigned int currval;
  for (unsigned int i = 1; i < nwalkers; ++i) {
    currval = nsteps[i];
    if (currval < minstep) minstep = currval;
  }
  return minstep;
}

/*!
  \param[in] other affineStepChunk to copy
  \returns Reference to updated chunk
*/
affineStepChunk& affineStepChunk::operator=(const affineStepChunk& other) {
  if (this == &other) return *this; //Self copy
  
  //Unless they are exactly the same size, have to reallocate
  unsigned int sz = nwalkers * niters * nparams;
  if ((nwalkers != other.nwalkers) ||
      (niters != other.niters)   ||
      (nparams != other.nparams)) {

    clear();
    nwalkers = other.nwalkers;
    niters = other.niters;
    nparams = other.nparams;
    sz = nwalkers * niters * nparams;

    if (sz > 0) {
      steps = new float[sz];
      logLike = new double[nwalkers*niters];
      nsteps = new unsigned int[nwalkers];
    } 
  }

  //Now copy

  if (sz > 0) {
    for (unsigned int i = 0; i < sz; ++i)
      steps[i] = other.steps[i];
    for (unsigned int i = 0; i < nwalkers * niters; ++i)
      logLike[i] = other.logLike[i];
    for (unsigned int i = 0; i < nwalkers; ++i)
      nsteps[i] = other.nsteps[i];
  }
  return *this;
}

/*!
  \returns Maximum log likelihood of any step in this chunk, or NaN
    if no steps available
*/
double affineStepChunk::getMaxLogLike() const {
  if (nwalkers == 0 || niters == 0 || nparams == 0) 
    return std::numeric_limits<double>::quiet_NaN();
  double *sgl_ptr;
  double val;
  double maxval;
  unsigned int curr_nsteps;

  //Find the first walker with anything in it
  unsigned int startidx;
  for (startidx = 0; startidx < nwalkers; ++startidx)
    if (nsteps[startidx] > 0) break;
  if (startidx == nwalkers) //Didn't find anything
    return std::numeric_limits<double>::quiet_NaN();
  maxval = logLike[startidx*niters];

  for (unsigned int i = startidx; i < nwalkers; ++i) {
    sgl_ptr = logLike + i * niters;
    curr_nsteps = nsteps[i];
    for (unsigned int j = 0; j < curr_nsteps; ++j) {
      val = sgl_ptr[j];
      if (val > maxval) maxval = val;
    }
  }
  return maxval;
}

/*!
  \param[out] maxval Maximum log likelihood in this chunk, 
    or NaN if no steps taken
  \param[out] p parameter values at maximum likelihood.  Unmodified
    if no steps in this chunk
*/
void affineStepChunk::getMaxLogLikeParam(double& maxval, paramSet& p) const {
  if (nwalkers == 0 || niters == 0 || nparams == 0) {
    maxval = std::numeric_limits<double>::quiet_NaN();
    return;
  }
  p.setNParams(nparams);

  unsigned int iloc,jloc,curr_nsteps;
  double *sgl_ptr;
  double val;
  unsigned int startidx;
  for (startidx = 0; startidx < nwalkers; ++startidx)
    if (nsteps[startidx] > 0) break;
  if (startidx == nwalkers) {
    //Didn't find anything
    maxval = std::numeric_limits<double>::quiet_NaN();
    return;
  }
  maxval = logLike[startidx * niters];
  iloc = startidx;
  jloc = 0;
  for (unsigned int i = startidx; i < nwalkers; ++i) {
    sgl_ptr = logLike + i * niters;
    curr_nsteps = nsteps[i];
    for (unsigned int j = 0; j < curr_nsteps; ++j) {
      val = sgl_ptr[j];
      if (val > maxval) {
	iloc = i; jloc = j;
	maxval = val;
      }
    }
  }
  p.setParamValues(nparams, steps + (iloc * niters + jloc) * nparams);
  return;
}

/*!
  \param[in] walker_idx Walker to add step to
  \param[in] pars Parameters of step to be added
  \param[in] lglike Log-likelihood of step

  \returns true if the parameter was set, false otherwise
*/
bool affineStepChunk::addStep(unsigned int walker_idx, 
			      const paramSet& pars, double lglike) {
  if (walker_idx >= nwalkers) return false;
  if (pars.getNParams() != nparams) return false;
  
  //Now make sure we aren't out of capacity!
  unsigned int iter = nsteps[walker_idx];
  if (iter >= niters) return false; //No more room

  float *ptr;
  ptr = getParamPointer(walker_idx, iter);
  for (unsigned int i = 0; i < nparams; ++i) ptr[i] = pars[i];
  nsteps[walker_idx] += 1;
  logLike[walker_idx * niters + iter] = lglike;  
  return true;
}

/*!
  \param[in] walker_idx Walker to get last step from
  \param[in] pars Parameters of last step
  \param[in] lglike Log-likelihood of step
  \returns true if it got the step, false otherwise
 */
bool affineStepChunk::getLastStep(unsigned int walker_idx,
				  paramSet& pars, double& lglike) const {
  if (walker_idx >= nwalkers) return false;
  if (pars.getNParams() != nparams) return false;

  unsigned int csteps = nsteps[walker_idx];
  if (csteps == 0) return false;
  pars.setParamValues(nparams, getParamPointer(walker_idx, csteps-1));
  lglike = getLogLike(walker_idx, csteps - 1);
  return true;
}


/*!
  \param[in] walker_idx Walker to replace last step with
  \param[in] pars New parameters for last step
  \param[in] lglike Log-likelihood of new step
  \returns true if it was able to do the replacement, false otherwise

  This should be used with care!
*/
bool affineStepChunk::replaceLastStep(unsigned int walker_idx,
				      const paramSet& pars, double lglike) {
  if (walker_idx >= nwalkers) return false;
  if (pars.getNParams() > nparams) return false;

  unsigned int csteps = nsteps[walker_idx];
  if (csteps == 0) return false;
  float *ptr;
  ptr = getParamPointer(walker_idx, csteps-1);
  for (unsigned int i = 0; i < nparams; ++i)
    ptr[i] = pars[i];
  logLike[walker_idx * niters + csteps - 1] = lglike;
  return true;
}

/*!
  \param[in] walker_idx Walker to replace last step with
  \param[in] iter_idx Step to get from that walker
  \param[out] pars parameters of that step
  \param[out] lglike Log-likelihood of that step
  \returns true if it got the step, false otherwise
*/
bool affineStepChunk::getStep(unsigned int walker_idx, unsigned int iter_idx,
			      paramSet& pars, double& lglike) const {
  if (walker_idx >= nwalkers) return false;
  if (pars.getNParams() != nparams) return false;

  unsigned int csteps = nsteps[walker_idx];
  if (iter_idx >= csteps) return false;
  pars.setParamValues(nparams, getParamPointer(walker_idx, iter_idx));
  lglike = getLogLike(walker_idx, iter_idx);
  return true;
}

/*!
  \param[in] walker_idx Walker to replace last step with
  \param[in] iter_idx Step to get from that walker
  \returns Pointer to that set of parameters, NULL if invalid step
*/
float* affineStepChunk::getParamPointer(unsigned int walker_idx,
					unsigned int iter_idx) {
  if (walker_idx >= nwalkers) return NULL;
  if (iter_idx >= niters) return NULL;
  return steps + (niters * walker_idx + iter_idx) * nparams;
}

/*!
  \param[in] walker_idx Walker to replace last step with
  \param[in] iter_idx Step to get from that walker
  \returns Pointer to const of that set of parameters, NULL if invalid step
*/
float* const affineStepChunk::getParamPointer(unsigned int walker_idx,
					      unsigned int iter_idx) const {
  if (walker_idx >= nwalkers) return NULL;
  if (iter_idx >= niters) return NULL;
  return steps + (niters * walker_idx + iter_idx) * nparams;
}

/*!
  \param[in] walker_idx Walker to replace last step with
  \returns Pointer to last set of parameters, NULL if no steps in this chunk
*/
float* affineStepChunk::getLastParamPointer(unsigned int walker_idx) {
  if (walker_idx >= nwalkers) return NULL;
  unsigned int iter_idxp1 = nsteps[walker_idx];
  if (iter_idxp1 == 0) return NULL;
  return steps + (niters * walker_idx + iter_idxp1 - 1) * nparams;
}

/*!
  \param[in] walker_idx Walker to replace last step with
  \returns Pointer to const of last set of parameters, NULL 
     if no steps in this chunk
*/
float* const affineStepChunk::getLastParamPointer(unsigned int walker_idx) 
  const {
  if (walker_idx >= nwalkers) return NULL;
  unsigned int iter_idxp1 = nsteps[walker_idx];
  if (iter_idxp1 == 0) return NULL;
  return steps + (niters * walker_idx + iter_idxp1 - 1) * nparams;
}

/*!
  \param[inout] os Stream to write to
*/
void affineStepChunk::writeToStream(std::ostream& os) const {
  //Do each walker in turn
  if (nparams == 0 || nwalkers == 0 || niters == 0) return;
  unsigned int curr_nsteps;
  float *ptr;
  os << std::left << std::scientific << std::showpoint;
  for (unsigned int i = 0; i < nwalkers; ++i) {
    curr_nsteps = nsteps[i];
    for (unsigned int j = 0; j < curr_nsteps; ++j) {
      ptr = getParamPointer(i, j);
      os << std::setw(13) << std::setprecision(6);
      os << ptr[0];
      for (unsigned int k = 1; k < nparams; ++k)
	os << " " << ptr[k];
      os << " " << std::setw(15) << std::setprecision(8)
	 << logLike[i * niters + j] << " " << i << std::endl;
    }
  }
}

/////////////////////////////////////////////////////

/*!
  \param[in] NWALKERS Number of walkers
  \param[out] NPARAMS Number of parameters
*/
affineChainSet::affineChainSet(unsigned int NWALKERS,
			       unsigned int NPARAMS) {
  nwalkers = NWALKERS;
  nparams = NPARAMS;
  skipfirst = false;
}

affineChainSet::~affineChainSet() {
  clear();
}

void affineChainSet::clear() {
  for (unsigned int i = 0; i < steps.size(); ++i)
    if (steps[i] != NULL) delete steps[i];
  steps.clear();
  skipfirst = false;
}

/*!
  Clear out all values, perserving the last step from the last chunk.
  Throws and exception if the chain is empty
*/
void affineChainSet::clearPreserveLast() throw (affineExcept) {
  if (steps.size() == 0)
    throw affineExcept("affineChainSet","clearPreserveLast",
		       "No previous entries to preserve",1);
  affineStepChunk *newchunk = new affineStepChunk(nwalkers,1,nparams);
  newchunk->fillFromLastStep(*steps.back());
  for (unsigned int i = 0; i < steps.size(); ++i)
    if (steps[i] != NULL) delete steps[i];
  steps.clear();
  steps.push_back(newchunk);
}

/*!
  \param[in] n New number of walkers.  Must be positive.
  
  This generally will not preserve the contents unless the number
  of walkers doesn't change
 */
void affineChainSet::setNWalkers(unsigned int n) throw (affineExcept) {
  if (n == nwalkers) return;

  if (n == 0)
    throw affineExcept("affineChainSet", "setNWalkers",
		       "n must be positive", 1);

  //We have to delete all previous steps
  clear();

  nwalkers = n;
}

/*!
  \param[in] n New number of parameters.  Must be positive.
  
  This generally will not preserve the contents unless the number
  of parameters doesn't change
 */
void affineChainSet::setNParams(unsigned int n) throw (affineExcept) {
  if (n == nparams) return;

  if (n == 0)
    throw affineExcept("affineChainSet", "setNParams",
		       "n must be positive", 1);

  //We have to delete all previous steps
  clear();

  nparams = n;
}

/*!
  \param[in] walker_idx Which walker to add the step to
  \param[in] p          Parameters to add
  \param[in] logLike    Log Likelihood of step
  \returns true if the step was added
*/
bool affineChainSet::addNewStep(unsigned int walker_idx, const paramSet& p, 
				double logLike) {
  return steps.back()->addStep(walker_idx, p, logLike);
}

/*!
  \param[in] sz Size (in number of steps) of new chunk to add
*/
void affineChainSet::addChunk(unsigned int sz) {
  affineStepChunk *chnkptr;
  chnkptr = new affineStepChunk(nwalkers, sz, nparams);
  steps.push_back(chnkptr);
}

/*!
  \param[in] walker_idx Which walker look for last step of
  \param[in] p          Parameters of last step
  \param[in] lgLike     Log Likelihood of last step

  Throws an exception if no steps taken
*/
void affineChainSet::getLastStep(unsigned int walker_idx,
				 paramSet& par, double& lglike) const 
  throw (affineExcept) {

  // The complication here is that the last step may not be in the current
  // chunk, because it's possible we just added a new empty one.
  int nchunks = static_cast<int>(steps.size());
  if (nchunks == 0) 
    throw affineExcept("affineChainSet", "getLastStep",
		       "No steps taken",1);
  if (walker_idx >= nwalkers)
    throw affineExcept("affineChainSet", "getLastStep",
		       "Walker index invalid", 2);

  if (par.getNParams() != nparams)
    par.setNParams(nparams);

  bool succ;
  unsigned int csteps;
  const affineStepChunk* chunk_ptr;
  succ = false;
  for (int i = nchunks-1; i >= 0; --i) {
    chunk_ptr = steps[i];
    csteps = chunk_ptr->nsteps[walker_idx];
    if (csteps == 0) continue; // Empty chunk, move on

    // Found one with a step -- copy it over and end
    par.setParamValues(nparams, 
		       chunk_ptr->getParamPointer(walker_idx, csteps-1));
    lglike = chunk_ptr->getLogLike(walker_idx, csteps-1);
    succ = true;
    break;
  }
  if (!succ)
    throw affineExcept("affineChainSet", "getLastStep",
		       "Error getting last step",3);
}

/*!
  \param[in] walker_idx Which walker to replace last step of
  \param[in] p          Parameters to replace values with
  \param[in] lgLike     Log Likelihood of replacement step

  This should be used with extreme care.  Throws an exception if there
  is no last step to replace.
*/
void affineChainSet::replaceLastStep(unsigned int walker_idx,
				     const paramSet& pars, double lglike) {
  int nchunks = static_cast<int>(steps.size());
  if (nchunks == 0) 
    throw affineExcept("affineChainSet", "replaceLastStep",
		       "No steps taken",1);
  if (walker_idx >= nwalkers)
    throw affineExcept("affineChainSet", "replaceLastStep",
		       "Walker index invalid", 2);

  bool succ;
  unsigned int csteps;
  affineStepChunk* chunk_ptr;
  succ = false;
  for (int i = nchunks-1; i >= 0; --i) {
    chunk_ptr = steps[i];
    csteps = chunk_ptr->nsteps[walker_idx];
    if (csteps == 0) continue; // Empty chunk, move on
    
    // Found one with a step -- replace
    succ = chunk_ptr->replaceLastStep(walker_idx, pars, lglike);
    break;
  }
  if (!succ)
    throw affineExcept("affineChainSet", "replaceLastStep",
		       "Error replacing last step", 3);
}

/*!
  \param[in] chunkidx   Which chunk to get the step from
  \param[in] walker_idx Which walker 
  \param[in] iter_idx   Which step
  \param[out] p         Parameters of that step
  \param[out] lgLike    Log Likelihood of that step

  Throws an exception if no steps taken
*/
void affineChainSet::getStep(unsigned int chunkidx, 
			     unsigned int walker_idx, 
			     unsigned int iter_idx,
			     paramSet& par, double& lglike) 
  const throw (affineExcept) {
  if (chunkidx >= steps.size()) 
    throw affineExcept("affineChainSet", "getStep",
		       "Chunk idx invalid", 1);
  bool succ = steps[chunkidx]->getStep(walker_idx, iter_idx, par, lglike);
  if (!succ)
    throw affineExcept("affineChainSet", "getStep",
		       "Error getting specified step", 2);
}

/*!
  \param[in] a Other affineChainSet to copy
  \returns Reference to lhs
*/
affineChainSet& affineChainSet::operator=(const affineChainSet& a) {
  if (this == &a) return *this; //Self copy
  clear(); //Throw away current contents
  nwalkers = a.nwalkers;
  nparams = a.nparams;
  skipfirst = a.skipfirst;
  
  unsigned int nchunk = a.steps.size();
  if (nchunk == 0) return *this; //No contents...
  steps.resize(nchunk);
  for (unsigned int i = 0; i < nchunk; ++i)
    *steps[i] = *a.steps[i];
  return *this;
}

/*!
  \param[in] other affineChainSet to append to this one
  \returns Reference to lhs
*/
affineChainSet& affineChainSet::operator+=(const affineChainSet& other) 
  throw (affineExcept) {
  //Now things like the number of params and walkers have to be the same
  if (nwalkers != other.nwalkers)
    throw affineExcept("affineChainSet","operator+=",
		       "nwalkers doesn't match",1);
  if (nparams != other.nparams)
    throw affineExcept("affineChainSet","operator+=",
		       "nparams doesn't match",2);

  unsigned int nchunk1 = steps.size();
  unsigned int nchunk2 = other.steps.size();
  
  if (nchunk2 == 0) return *this; //Nothing to append
  if (nchunk1 == 0) return (*this = other); //Nothing in this one, so copy
  steps.reserve(nchunk1+nchunk2);
  for (unsigned int i = 0; i < nchunk2; ++i)
    steps.push_back(new affineStepChunk(*other.steps[i]));
  return *this;
}

/*!
  \returns Total number of iterations across all chunks
*/
unsigned int affineChainSet::getNIters() const {
  unsigned int ntotiter = 0;
  const affineStepChunk* chunk;
  unsigned int firststep = 0;
  if (skipfirst) firststep = 1;
  for (unsigned int i = firststep; i < steps.size(); ++i) {
    chunk = steps[i];
    for (unsigned int j = 0; j < nwalkers; ++j)
      ntotiter += chunk->nsteps[j];
  }
  return ntotiter;
}

/*!
  \param[in] walker Which walker
  \returns Total number of iterations across all chunks for specified walker
*/
unsigned int affineChainSet::getNIters(unsigned int walker) const 
  throw (affineExcept) {
  unsigned int ntotiter = 0;
 if (walker >= nwalkers) 
     throw affineExcept("affineChainSet","getNIters",
			"Specified walker exceeds number available",1);
 unsigned int firststep = 0;
 if (skipfirst) firststep = 1;
 for (unsigned int i = firststep; i < steps.size(); ++i)
   ntotiter += steps[i]->nsteps[walker];
 return ntotiter;
}

/*!
  \returns Number of iterations for walker with the smallest number

  Requires auxilliary storage equal to the number of walkers
*/
unsigned int affineChainSet::getMinNIters() const {
  unsigned int *iterarr;
  if (steps.size() == 0) return 0;
  iterarr = new unsigned int[nwalkers];
  const affineStepChunk *chunkptr;
  unsigned int firststep = 0;
  if (skipfirst) firststep = 1;

  chunkptr = steps[firststep];
  for (unsigned int i = 0; i < nwalkers; ++i)
    iterarr[i] = chunkptr->nsteps[i];
  for (unsigned int i = firststep+1; i < steps.size(); ++i) {
    chunkptr = steps[i];
    for (unsigned int j = 0; j < nwalkers; ++j)
      iterarr[j] += chunkptr->nsteps[j];
  }

  unsigned int currval, miniter;
  miniter = iterarr[0];
  for (unsigned int i = 1; i < nwalkers; ++i) {
    currval = iterarr[i];
    if (currval < miniter) miniter = currval;
  }
  delete[] iterarr;
  return miniter;
}

/*!
  \returns Maximum log likelihood of all steps taken
*/
double affineChainSet::getMaxLogLike() const {
  if (steps.size() == 0) return std::numeric_limits<double>::quiet_NaN();
  unsigned int firststep = 0;
  if (skipfirst) firststep = 1;
  double maxval = steps[firststep]->getMaxLogLike();
  double currval;
  for (unsigned int i = firststep+1; i < steps.size(); ++i) {
    currval = steps[i]->getMaxLogLike();
    if (currval > maxval) maxval = currval;
  }
  return maxval;
}

/*!
  \param[out] maxval Maximum log likelihood of all steps
  \param[out] p paramSet corresponding to that likelihood
*/
void affineChainSet::getMaxLogLikeParam(double& maxval, paramSet& p) const {
  p.setNParams(nparams);

  if (steps.size() == 0) {
    maxval = std::numeric_limits<double>::quiet_NaN();
    return;
  }

  double currval;
  unsigned int firststep = 0;
  if (skipfirst) firststep = 1;
  steps[firststep]->getMaxLogLikeParam(maxval, p);

  paramSet cp(nparams);
  for (unsigned int i = firststep + 1; i < steps.size(); ++i) {
    steps[i]->getMaxLogLikeParam(currval, cp);
    if (currval > maxval) {
      maxval = currval;
      p = cp;
    }
  }
}
 
/*!
  \param[in] walkeridx Walker index
  \param[in] paramidx Parameter index
  \param[out] pvec Contains all the steps for that parameter for specified
     walker
*/
void affineChainSet::getParamVector(unsigned int walkeridx,
				    unsigned int paramidx,
				    std::vector<float>& pvec) const 
  throw (affineExcept) {
  if (walkeridx >= nwalkers) 
    throw affineExcept("affineChainSet","getParamVector",
		       "Specified walkeridx exceeds number available",1);
  if (paramidx >= nparams) 
    throw affineExcept("affineChainSet","getParamVector",
		       "Specified paramidx exceeds number available",2);
  
  unsigned int sz = getNIters(walkeridx);
  pvec.resize(sz);
  if (sz == 0) return;
  
  unsigned int chunksz, ctr = 0, nciters, ncpars;
  const affineStepChunk* chunkptr;
  float *sgl_ptr;
  unsigned int firststep = 0;
  if (skipfirst) firststep = 1;
  for (unsigned int i = firststep; i < steps.size(); ++i) {
    chunkptr = steps[i];
    nciters = chunkptr->niters;
    ncpars  = chunkptr->nparams;
    chunksz = chunkptr->nsteps[walkeridx];
    sgl_ptr = chunkptr->steps+(walkeridx*nciters*ncpars);
    for (unsigned int j = 0; j < chunksz; ++j)
      pvec[ctr++] = sgl_ptr[j * ncpars + paramidx];
  }
}

/*!
  \param[in] paramidx parameter index
  \param[out] pvec Parameter values averaged across walkers.  Resized
    to the minimum number of steps taken for any walker.

  Only includes steps where all walkers have that step.
*/
void affineChainSet::getAverageParamVector(unsigned int paramidx,
					   std::vector<float>& pvec) const 
  throw (affineExcept) {
  if (paramidx >= nparams) 
    throw affineExcept("affineChainSet","getAverageParamVector",
		       "Specified paramidx exceeds number available",1);
  
  pvec.resize(getMinNIters());

  if (nwalkers == 1) {
    getParamVector(0, paramidx, pvec);
    return;
  }

  unsigned int ctr, nsteps, currsteps, nciters, ncpars;
  double cumval; //Do in double precision, then convert back
  double norm = 1.0 / static_cast<double>(nwalkers);
  ctr = 0;
  const affineStepChunk* chunkptr;
  unsigned int firststep = 0;
  if (skipfirst) firststep = 1;
  for (unsigned int i = firststep; i < steps.size(); ++i) {
    chunkptr = steps[i];
    nciters = chunkptr->niters;
    ncpars  = chunkptr->nparams;

    //First find the minimum number of steps for this param in this
    // chunk, store in nsteps
    nsteps = chunkptr->nsteps[0];
    for (unsigned int j = 1; j < nwalkers; ++j) {
      currsteps = chunkptr->nsteps[i];
      if (currsteps < nsteps) nsteps = currsteps;
    }
    if (nsteps == 0) continue; //None

    for (unsigned int j = 0; j < nsteps; ++j) {
      cumval = chunkptr->steps[j * ncpars + paramidx]; //Walker 0
      for (unsigned int k = 1; k < nwalkers; ++k)
	cumval += chunkptr->steps[(k * nciters + j) * ncpars + paramidx];
      pvec[ctr++] = static_cast<float>(norm * cumval);
    }
  }
  pvec.resize(ctr); //Should be non-destructive
}

/*!
  \param[out] mean   Mean of vector
  \param[out] sigma  Sigma of vector
  \param[out] tau    Autocorrelation length estimate of vector
  \param[in]  X      Data vector
  \param[in]  L      Number of elements in data vector

  Copied from the ACOR package by Jonathan Goodman:
  goodman@cims.nyu.edu, http://www.cims.nyu.edu/faculty/goodman

  Calls itself recursively, destroying X in the process.
*/
//  The way this handles errors is a bit funny, but is based on the orignial
//  code.  Basically, it's a good idea to do the first L < minfac*maxlag
//  test before calling this.
int affineChainSet::acor(double& mean, double& sigma,
			 double& tau, std::vector<float>& X,
			 unsigned int L) 
  const throw (affineExcept) {

  // Compute tau directly only if tau < taumax, otherwise use the pairwise sum
  double invL = 1.0 / static_cast<double>(L);

  // Compute the mean of X ... 
  mean = X[0];
  for (unsigned int i = 1; i < L; ++i) mean += X[i];
  mean *= invL;
  //  ... and subtract it away.
  for (unsigned int i = 0; i <  L; ++i) X[i] -= mean; 

  //If chunk is too small, stop recursion and just leave initial values alone
  if (L < minfac * maxlag) return 1;
   
  double C[maxlag + 1]; 
  // Here, s=0 is the variance, s = maxlag is the last one computed.
  for (unsigned int s = 0; s <= maxlag; ++s)  C[s] = 0.;  
  
  // Compute the autocovariance function . . . 
  unsigned int iMax = L - maxlag;
  double invImax = 1.0 / static_cast<double>(iMax);
  for (unsigned int i = 0; i < iMax; ++i) 
    // ...  first the inner products ...
    for (unsigned int s = 0; s <= maxlag; ++s)
      C[s] += X[i] * X[i + s];                              
  // ...  then the normalization.
  for (unsigned int s = 0; s <= maxlag; ++s) C[s] *= invImax;
  
  // The "diffusion coefficient" is the sum of the autocovariances
  double D = C[0];   
  // The rest of the C[s] are double counted since C[-s] = C[s].
  for (unsigned int s = 1; s <= maxlag; ++s) D += 2 * C[s];   
  // The standard error bar formula, if D were the complete sum.
  sigma = std::sqrt(D * invL);
  // A provisional estimate, since D is only part of the complete sum.
  tau   = D / C[0]; 
  
  // Stop if the D sum includes the given multiple of tau.
  // This is the self consistent window approach.
  if (tau * winmult < maxlag) return 0;    
  
  // If the provisional tau is so large that we don't think tau
  // is accurate, apply the acor procedure to the pairwase sums
  // of X.
  // The pairwise sequence is half the length (if L is even)
  unsigned int Lh = L/2;                          
  // The mean of the new sequence, to throw away.         
  double newMean;                                 
  unsigned int j1 = 0;
  unsigned int j2 = 1;
  for (unsigned int i = 0; i < Lh; ++i) {
    X[i] = X[j1] + X[j2];
    j1 += 2;
    j2 += 2; 
  }
  //Recursive call; note we ignore the return value!
  acor(newMean, sigma, tau, X, Lh); 
  
  // Reconstruct the fine time series numbers from the coarse series numbers.
  D     = 0.25 * sigma * sigma * static_cast<double>(L);    
  tau   = D / C[0];                 // As before, but with a corrected D.
  sigma = std::sqrt(D * invL);  // As before, again.
  return 0;
}

/*!
  \param[in] paramidx   Index of paramter
  \param[out] mean      Mean of parameter
  \param[out] sigma     Sigma of parameter
  \param[out] succ      True if the autocorrelation was computed
  \returns \f$\tau\f$, the estimated autocorrelation length
*/
double affineChainSet::getAcor(unsigned int paramidx, double& mean,
			       double& sigma, bool& succ) 
  const throw (affineExcept) {
  if (paramidx >= nparams) 
    throw affineExcept("affineChainSet","getAcor",
		       "Specified paramidx exceeds number available",1);

  double tau;
  sigma = std::numeric_limits<double>::quiet_NaN();
  tau = std::numeric_limits<double>::quiet_NaN();

  //Get the averaged parameters
  std::vector<float> pvec;
  getAverageParamVector(paramidx, pvec);

  //Make sure initial vector is long enough to allow for computation
  if (pvec.size() < minfac*maxlag) {
    succ = false;
    return std::numeric_limits<double>::quiet_NaN();
  }

  acor(mean, sigma, tau, pvec, pvec.size());
  succ = true;
  return tau;
}

/*!
  \param[out] tau The autocorrelation vector
  \returns True if the autocorrelation vector was computed for all params
*/
bool affineChainSet::getAcorVector(std::vector<float>& tau) const 
  throw (affineExcept) {
  double mn, sigma; //Throw away
  tau.resize(nparams);
  bool succ, indiv_succ;
  succ = true;
  for (unsigned int i = 0; i < nparams; ++i) {
    tau[i] = static_cast<float>(getAcor(i, mn, sigma, indiv_succ));
    succ &= indiv_succ;
  }
  return succ;
}

/*!
  \param[out] tau Autocorrelation length for each parameter
  \param[in] ignore Set to mcmc_affine::ACIGNORE where to ignore params.  
                    tau is set to NaN where this is true
  \returns True if the autocorrelation vector was computed for all params
            that were not ignored
*/
bool affineChainSet::getAcorVector(std::vector<float>& tau,
				   const std::vector<int> param_state) 
  const 
  throw (affineExcept) {
  if (param_state.size() < nparams)
    throw affineExcept("affineChainSet","getAcorVector",
		       "param_state is shorter than number of params",1);
  double mn, sigma; //Throw away
  tau.resize(nparams);
  bool succ, indiv_succ;
  succ = true;
  for (unsigned int i = 0; i < nparams; ++i) {
    if (param_state[i] & mcmc_affine::ACIGNORE)
      tau[i] = std::numeric_limits<float>::quiet_NaN();
    else {
      tau[i] = static_cast<float>(getAcor(i, mn, sigma, indiv_succ));
      succ &= indiv_succ;
    }
  }
  return succ;
}

/*
  \params[in] paridx Parameter index
  \returns Mean of parameter
*/
float affineChainSet::getParamMean(unsigned int paridx) const 
  throw (affineExcept) {
  if (paridx >= nparams)
    throw affineExcept("affineChainSet","getParamMean",
		       "Invalid parameter index",1);
  double mean = 0; //Do computation internally in double
  unsigned int chunksz, nciters, ncpars, ctr = 0;
  const affineStepChunk* chunkptr;
  float *ptr;
  for (unsigned int i = 0; i < steps.size(); ++i) {
    chunkptr = steps[i];
    nciters = chunkptr->niters;
    ncpars  = chunkptr->nparams;
    for (unsigned int j = 0; j < nwalkers; ++j) {
      chunksz = chunkptr->nsteps[j];
      ptr = chunkptr->steps + j * nciters * ncpars;
      for (unsigned int k = 0; k < chunksz; ++k) {
	mean += ptr[j * ncpars + paridx];
	++ctr;
      }
    }
  }
  if (ctr == 0) return std::numeric_limits<float>::quiet_NaN();
  else {
    mean /= static_cast<double>(ctr);
    return static_cast<float>(mean);
  }
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
void affineChainSet::getParamStats(unsigned int paridx, float& mean,
				   float& var, float& lowlimit,
				   float& uplimit, float conflevel) const 
  throw (affineExcept) {
  if (paridx >= nparams)
    throw affineExcept("affineChainSet", "getParamStats",
		       "Invalid parameter index", 1);

  mean = std::numeric_limits<float>::quiet_NaN();
  var = std::numeric_limits<float>::quiet_NaN();
  lowlimit = std::numeric_limits<float>::quiet_NaN();
  uplimit = std::numeric_limits<float>::quiet_NaN();

  //We will copy the full params into a vector for ease of sorting,
  // computation
  unsigned int sz = getNIters();
  if (sz == 0) return;
  if (skipfirst && steps.size() == 1) return;

  std::vector<float> pvec;
  pvec.resize(sz);
  
  //Stick in vector and add to mean as we go.  Use double internally
  float val, *ptr;
  double imean;
  unsigned int chunksz, firstidx, ctr, nciters, ncpars;
  const affineStepChunk* chunkptr;
  imean = 0.0;
  ctr = firstidx = 0;
  if (skipfirst) firstidx = 1;
  for (unsigned int stepidx = firstidx; stepidx < steps.size(); ++stepidx) {
    //For each chunk
    chunkptr = steps[stepidx];
    nciters = chunkptr->niters;
    ncpars = chunkptr->nparams;
    for (unsigned int i = 0; i < nwalkers; ++i) {
      //For each walker
      chunksz = chunkptr->nsteps[i];
      ptr = chunkptr->steps + i * nciters * ncpars;
      for (unsigned int j = 0; j < chunksz; ++j) {
	//For each step in each chunk
	//ptr = chunkptr->getParamPointer(i,j);
	//val = ptr[ paridx ];
	val = ptr[j * ncpars + paridx];
	imean += val;
	pvec[ctr++] = val;
      }
    }
  }
  if (ctr != sz) 
    throw affineExcept("affineChainSet", "getParamStats",
		       "Didn't get the number of expected steps",2);

  double norm = 1.0 / static_cast<double>(sz);
  imean *= norm;

  mean = static_cast<float>(imean);

  //Now compute variance using two-pass algorithm
  // Again, work in double precision internally
  double delta, sumdelta, ivar;
  ivar = sumdelta = 0.0;
  for (unsigned int i = 0; i < sz; ++i) {
    delta = pvec[i] - imean;
    sumdelta += delta;
    ivar += delta * delta;
  }
  ivar -= sumdelta * norm;
  ivar /= static_cast<double>(sz) - 1.0;
  var = static_cast<float>(ivar);
  if (var == 0.0) return; //Can't compute limits
  
  //Compute limits
  if (conflevel < 0)
    throw affineExcept("affineChainSet", "getParamStats",
		       "Invalid (negative) conf level",3 );
  if (conflevel > 1)
    throw affineExcept("affineChainSet", "getParamStats",
		       "Invalid (>1) conf level", 4);

  std::sort(pvec.begin(), pvec.end());
  //Now do the lower limit, rounding to get that index
  //I could add interpolation here, but not yet
  double lowfrac = 0.5 * (1.0 - conflevel);
  int lowidx = static_cast<int>(lowfrac * sz);
  if (lowidx != 0) lowlimit = pvec[lowidx];
  //Upper limit is the same
  double upfrac = 1.0 - 0.5 * conflevel;
  int upidx = static_cast<int>(upfrac * sz);
  if (upidx != 0) uplimit = pvec[upidx];
}

/*!
  \param[inout] covmat Covariance matrix.  Must be pre-allocated
    to nparams by nparams.

  Requires temporary storage equal to the full size of the chain
*/
void affineChainSet::makeCovMatrix(float** covmat) const {

  //Do internal computations in double precision

  if (nparams == 0 || nwalkers == 0) return;
  unsigned int sz = getMinNIters();
  if (sz < 2) {
    for (unsigned int i = 0; i < nparams; ++i)
      for (unsigned int j = 0; j < nparams; ++j)
	covmat[i][j] = std::numeric_limits<float>::quiet_NaN();
    return;
  }

  //We copy everything into an array, nparams x sz
  float *tmpdata;
  tmpdata = new float[nparams * sz];
  unsigned int chunksz, ctr, ctrbase, nciters, ncpars;
  const affineStepChunk* chunkptr;
  float *ptr;
  //Take at most sz from each
  for (unsigned int i = 0; i < nparams; ++i) {
    ctrbase = sz * i;
    ctr = 0;
    for (unsigned int j = 0; j < steps.size(); ++j) if (ctr <= sz) {
	chunkptr = steps[j];
	nciters = chunkptr->niters;
	ncpars = chunkptr->nparams;
	for (unsigned int k = 0; k < nwalkers; ++k) if (ctr <= sz) {
	    chunksz = chunkptr->nsteps[k];
	    ptr = chunkptr->steps + k * nciters * ncpars;
	    for (unsigned int m = 0; (m < chunksz) && (ctr <= sz); ++m) {
	      tmpdata[ctrbase + ctr] = ptr[j * ncpars + i];
	      ++ctr;
	    }
	  }
      }
  }

  //Now, mean subtract each parameter value
  double mean, norm;
  norm = 1.0 / static_cast<double>(sz);
  for (unsigned int i = 0; i < nparams; ++i) {
    ptr = tmpdata + sz * i;
    mean = ptr[0];
    for (unsigned int j = 1; j < sz; ++j)
      mean += ptr[j];
    mean *= norm;
    for (unsigned int j = 0; j < sz; ++j)
      ptr[j] -= static_cast<float>(mean);
  }

  //And compute covariances
  float *ptr2;
  double val, sumsq, norm2;
  norm2 = 1.0 / static_cast<double>(sz-1);
  for (unsigned int i = 0; i < nparams; ++i) {
    ptr = tmpdata + sz * i;
    //Diagonal
    val = static_cast<double>(ptr[0]);
    sumsq = val * val;
    for (unsigned int j = 1; j < sz; ++j) {
      val = static_cast<double>(ptr[j]);
      sumsq += val * val;
    }
    covmat[i][i] = static_cast<float>(sumsq * norm2);

    //Off diag, take advantage of symmetry
    for (unsigned int j = 0; j < i; ++j) {
      ptr2 = tmpdata + sz * j;
      sumsq = static_cast<double>(ptr[0] * ptr2[0]);
      for (unsigned int k = 1; k < sz; ++k)
	sumsq += static_cast<double>(ptr[k] * ptr2[k]);
      covmat[i][j] = covmat[j][i] = static_cast<float>(sumsq * norm2);
    }
  }

  delete[] tmpdata;
}

/*!
  \param[in] outfile File to write chains to
*/
void affineChainSet::writeToFile(const std::string& outfile) const 
  throw (affineExcept) {
  if (steps.size() == 0) return;
  std::ofstream ofs;
  ofs.open(outfile.c_str());
  if (! ofs) {
    ofs.close();
    throw affineExcept("affineChainSet","writeToFile",
		       "Failed to open file",1);
  }
  unsigned int startidx = 0;
  if (skipfirst) startidx = 1;
  for (unsigned int i = startidx; i < steps.size(); ++i)
    steps[i]->writeToStream(ofs);
  ofs.close();
}

/*!
  \param[in] filename HDF5 File to write to.  Will be overwritten if it exists.

  Only works if all walkers have the same number of steps.
*/
void affineChainSet::writeToHDF5(const std::string& filename) const {
  if (steps.size() == 0) return;

  // Open
  hid_t file_id;
  herr_t status;
  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
		      H5P_DEFAULT);
  // Write
  writeToHDF5(file_id);

  // All done
  status = H5Fclose(file_id);
}


/*!
  \param[in] objid HDF5 group to write to

  Only works if all walkers have the same number of steps.
  objid can be either a file id or a group id.  If it's a file id,
  everything is written to the root of that file.
*/
void affineChainSet::writeToHDF5(hid_t objid) const {
  if (steps.size() == 0) return;

  // Make sure all the walkers have the same number of iterations
  unsigned int nit = getNIters(0);
  unsigned int nc;
  for (unsigned int i = 1; i < nwalkers; ++i) {
    nc = getNIters(i);
    if (nc != nit) {
      std::stringstream errstr;
      errstr << "Walker " << i << " has different number of steps ("
	     << nc << ") than first walker (" << nit << ")";
      throw affineExcept("affineChainSet", "writeToHDF5", 
			 errstr.str(), 1);
    }
  }

  // Two datasets to write -- the steps (nwalkers * nsteps * nparams)
  //  and the likelihoods (nwalkers * nsteps).  
  // Start with the chain steps, writing in chunks
  //  of size nsteps * nparams, copying into a local variable
  hsize_t dims_steps[3];
  dims_steps[0] = static_cast<hsize_t>(nwalkers);
  dims_steps[1] = static_cast<hsize_t>(nit);
  dims_steps[2] = static_cast<hsize_t>(nparams);
  hsize_t offset_steps[3] = {0, 0, 0};
  hsize_t memsdim_steps[3]; // Working space for chunks we are writing
  memsdim_steps[0] = 1;
  memsdim_steps[1] = static_cast<hsize_t>(nit);
  memsdim_steps[2] = static_cast<hsize_t>(nparams);
  hsize_t count_steps[3] = {1, 1, 1}; // We write one of these blocks at a time

  // Data space
  hid_t dataspace_idsteps;
  herr_t status;
  dataspace_idsteps = H5Screate_simple(3, dims_steps, NULL);
  
  // Data set
  hid_t dataset_idsteps;
  dataset_idsteps = H5Dcreate2(objid, "chain", H5T_NATIVE_FLOAT, 
			       dataspace_idsteps, H5P_DEFAULT, H5P_DEFAULT, 
			       H5P_DEFAULT);
  
  // Write some attributes for that dataset
  hsize_t adims;
  adims = 1;
  hid_t mems_id, att_id;
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate2(dataset_idsteps, "nwalkers", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(att_id, H5T_NATIVE_UINT, &nwalkers);
  status = H5Aclose(att_id);
  att_id = H5Acreate2(dataset_idsteps, "nsteps", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(att_id, H5T_NATIVE_UINT, &nit);
  status = H5Aclose(att_id);
  att_id = H5Acreate2(dataset_idsteps, "nparams", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(att_id, H5T_NATIVE_UINT, &nparams);
  status = H5Aclose(att_id);
  status = H5Sclose(mems_id);

  // Write the data in chunks
  unsigned int cntr, firststep, nstp;
  float *fwrkarr, *sptr;
  const affineStepChunk* chunk;
  fwrkarr = new float[nit * nparams];
  if (skipfirst) firststep = 1; else firststep = 0;
  mems_id = H5Screate_simple(3, memsdim_steps, NULL);
  for (unsigned int walkidx = 0; walkidx < nwalkers; ++walkidx) {
    // Copy all steps for this walker into fwrkarr
    cntr = 0;
    for (unsigned int chunkidx = firststep; 
	 chunkidx < steps.size(); ++chunkidx) {
      chunk = steps[chunkidx];
      nstp = chunk->nsteps[chunkidx];
      for (unsigned int iteridx = 0; iteridx < nstp; ++iteridx) {
	sptr = chunk->getParamPointer(walkidx, iteridx);
	for (unsigned int pidx = 0; pidx < nparams; ++pidx)
	  fwrkarr[cntr++] = sptr[pidx];
      }
    }

    // Now write chunk
    offset_steps[0] = static_cast<hsize_t>(walkidx);
    status = H5Sselect_hyperslab(dataspace_idsteps, H5S_SELECT_SET, 
				 offset_steps, NULL, count_steps, 
				 memsdim_steps);
    status = H5Dwrite(dataset_idsteps, H5T_NATIVE_FLOAT, mems_id,
		      dataspace_idsteps, H5P_DEFAULT, fwrkarr);
  }
  delete[] fwrkarr;

  // Close up step variables
  status = H5Sclose(mems_id);
  status = H5Dclose(dataset_idsteps);
  status = H5Sclose(dataspace_idsteps);

  // Now do the same thing for the likelihood, writing in nwalker
  // chunks of size nsteps
  hsize_t dims_like[2];
  dims_like[0] = static_cast<hsize_t>(nwalkers);
  dims_like[1] = static_cast<hsize_t>(nit);
  hsize_t offset_like[2] = {0, 0};
  hsize_t memsdim_like[2];
  memsdim_like[0] = 1;
  memsdim_like[1] = static_cast<hsize_t>(nit);
  hsize_t count_like[2] = {1, 1};
  hid_t dataspace_idlike, dataset_idlike;
  dataspace_idlike = H5Screate_simple(2, dims_like, NULL);
  dataset_idlike = H5Dcreate2(objid, "likelihood", H5T_NATIVE_DOUBLE, 
			      dataspace_idlike, H5P_DEFAULT, H5P_DEFAULT, 
			      H5P_DEFAULT);
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate2(dataset_idlike, "nwalkers", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(att_id, H5T_NATIVE_UINT, &nwalkers);
  status = H5Aclose(att_id);
  att_id = H5Acreate2(dataset_idlike, "nsteps", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(att_id, H5T_NATIVE_UINT, &nit);
  status = H5Aclose(att_id);
  status = H5Sclose(mems_id);
  double *dwrkarr;
  dwrkarr = new double[nit];
  mems_id = H5Screate_simple(2, memsdim_like, NULL);
  for (unsigned int walkidx = 0; walkidx < nwalkers; ++walkidx) {
    cntr = 0;
    for (unsigned int chunkidx = firststep; 
	 chunkidx < steps.size(); ++chunkidx) {
      chunk = steps[chunkidx];
      nstp = chunk->nsteps[chunkidx];
      for (unsigned int iteridx = 0; iteridx < nstp; ++iteridx)
	dwrkarr[cntr++] = chunk->getLogLike(walkidx, iteridx);
    }

    // Now write chunk
    offset_like[0] = static_cast<hsize_t>(walkidx);
    status = H5Sselect_hyperslab(dataspace_idlike, H5S_SELECT_SET, 
				 offset_like, NULL, count_like, 
				 memsdim_like);
    status = H5Dwrite(dataset_idlike, H5T_NATIVE_DOUBLE, mems_id,
		      dataspace_idlike, H5P_DEFAULT, dwrkarr);
  }
  delete[] dwrkarr;
  status = H5Sclose(mems_id);
  status = H5Dclose(dataset_idlike);
  status = H5Sclose(dataspace_idlike);

}
