#include<limits>
#include<cmath>
#include<fstream>
#include<algorithm>

#include <affineChainSet.h>

//Parameters for acor computation
const unsigned int affineChainSet::taumax = 2; 
//Autocovariance up to lag = winmult*tau
const unsigned int affineChainSet::winmult = 5; 
const unsigned int affineChainSet::maxlag = taumax*winmult;
const unsigned int affineChainSet::minfac = 5;


affineStepChunk::affineStepChunk(unsigned int NWALKERS,
				 unsigned int NITERS,
				 unsigned int NPARAMS) {
  nwalkers = NWALKERS;
  niters   = NITERS;
  nparams  = NPARAMS;

  //Allocate in large contiguous block for memory access efficiency
  unsigned int sz = nwalkers*niters*nparams;
  if (sz > 0) {
    steps = new double[sz];
    logLike = new double[nwalkers*niters];
    nsteps = new unsigned int[nwalkers];
    for (unsigned int i = 0; i < nwalkers; ++i)
      nsteps[i] = 0;
  } else {
    steps = NULL;
    logLike = NULL;
    nsteps  = NULL;
  }
}

//copy constructor
affineStepChunk::affineStepChunk(const affineStepChunk& other) {
  nwalkers = other.nwalkers;
  niters   = other.niters;
  nparams  = other.nparams;

  unsigned int sz = nwalkers*niters*nparams;
  if (sz > 0) {
    steps = new double[sz];
    for (unsigned int i = 0; i < sz; ++i)
      steps[i] = other.steps[i];
    logLike = new double[nwalkers*niters];
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
  Fills a chunk with the last step of another chunk, resizing as
  needed
*/
void affineStepChunk::fillFromLastStep(const affineStepChunk& other)
  throw (affineExcept) {
  unsigned int sz;
  sz = other.nwalkers*other.niters*other.nparams;
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
    nparams  = other.nparams;
    niters   = 1;
    sz = nwalkers*niters*nparams;
    if (sz > 0) {
      steps = new double[sz];
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
  } else {
    sz = nwalkers*niters*nparams;
  }
  //Now, for each walker, set the last param
  unsigned int laststep;
  double *ptro, *ptrn;
  for (unsigned int i = 0; i < nwalkers; ++i) {
    laststep = other.nsteps[i]-1;
    ptro = other.getParamPointer(i,laststep);
    ptrn = getParamPointer(i,0);
    for (unsigned int j = 0; j < nparams; ++j)
      ptrn[j] = ptro[j];
    //Recall logLike is nwalkers * niters, and niters is one
    logLike[i] = other.getLogLike(i,laststep);
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

affineStepChunk& affineStepChunk::operator=(const affineStepChunk& other) {
  if (this == &other) return *this; //Self copy
  
  //Unless they are exactly the same size, have to reallocate
  if ( (nwalkers != other.nwalkers) ||
       (niters   != other.niters)   ||
       (nparams  != other.nparams) ) {
    clear();
    nwalkers = other.nwalkers;
    niters   = other.niters;
    nparams  = other.nparams;

    unsigned int sz = nwalkers*niters*nparams;
    if (sz > 0) {
      steps = new double[sz];
      logLike = new double[nwalkers*niters];
      nsteps = new unsigned int[nwalkers];
    } else {
      steps = NULL;
      logLike = NULL;
      nsteps = NULL;
    }
  }

  //Now copy
  unsigned int sz = nwalkers*niters*nparams;
  if (sz > 0) {
    for (unsigned int i = 0; i < sz; ++i)
      steps[i] = other.steps[i];
    for (unsigned int i = 0; i < nwalkers*niters; ++i)
      logLike[i] = other.logLike[i];
    for (unsigned int i = 0; i < nwalkers; ++i)
      nsteps[i] = other.nsteps[i];
  }
  return *this;
}

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
    sgl_ptr = logLike+i*niters;
    curr_nsteps = nsteps[i];
    for (unsigned int j = 0; j < curr_nsteps; ++j) {
      val = sgl_ptr[j];
      if (val > maxval) maxval = val;
    }
  }
  return maxval;
}

void affineStepChunk::getMaxLogLikeParam(double& maxval,paramSet& p) const {
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
    val = std::numeric_limits<double>::quiet_NaN();
    return;
  }
  maxval = logLike[startidx*niters];
  iloc = startidx;
  jloc = 0;
  for (unsigned int i = startidx; i < nwalkers; ++i) {
    sgl_ptr = logLike + i*niters;
    curr_nsteps = nsteps[i];
    for (unsigned int j = 0; j < curr_nsteps; ++j) {
      val = sgl_ptr[j];
      if (val > maxval) {
	iloc = i; jloc = j;
	maxval = val;
      }
    }
  }
  p.setParamValues(nparams,steps+(iloc*niters+jloc)*nparams);
  return;
}

/*!
  returns true if the parameter was set, false otherwise
 */
bool affineStepChunk::addStep(unsigned int walker_idx, 
			      const paramSet& pars, double lglike) {
  if (walker_idx >= nwalkers) return false;
  if (pars.getNParams() != nparams) return false;
  
  //Now make sure we aren't out of capacity!
  unsigned int iter = nsteps[walker_idx];
  if (iter >= niters) return false; //No more room

  logLike[walker_idx*niters+iter] = lglike;  
  double *ptr;
  ptr = getParamPointer(walker_idx,iter);
  for (unsigned int i = 0; i < nparams; ++i)
    ptr[i] = pars[i];
  nsteps[walker_idx] += 1;
  return true;
}

/*!
  \returns true if it got the step, false otherwise
 */
bool affineStepChunk::getLastStep(unsigned int walker_idx,
				  paramSet& pars, double& lglike) const {
  if (walker_idx >= nwalkers) return false;
  if (pars.getNParams() != nparams) return false;

  unsigned int csteps = nsteps[walker_idx];
  if (csteps == 0) return false;
  pars.setParamValues( nparams, getParamPointer(walker_idx,csteps-1) );
  lglike = getLogLike(walker_idx,csteps-1);
  return true;
}

/*!
  \returns true if it got the step, false otherwise
 */
bool affineStepChunk::getStep(unsigned int walker_idx, unsigned int iter_idx,
			      paramSet& pars, double& lglike) const {
  if (walker_idx >= nwalkers) return false;
  if (pars.getNParams() != nparams) return false;

  unsigned int csteps = nsteps[walker_idx];
  if (iter_idx >= csteps) return false;
  pars.setParamValues( nparams, getParamPointer(walker_idx,iter_idx) );
  lglike = getLogLike(walker_idx,iter_idx);
  return true;
}


double* affineStepChunk::getParamPointer( unsigned int walker_idx,
					  unsigned int iter_idx ) {
  if (walker_idx >= nwalkers) return NULL;
  if (iter_idx >= niters) return NULL;
  return steps + (niters*walker_idx + iter_idx)*nparams;
}

double* const affineStepChunk::getParamPointer( unsigned int walker_idx,
						unsigned int iter_idx ) const {
  if (walker_idx >= nwalkers) return NULL;
  if (iter_idx >= niters) return NULL;
  return steps + (niters*walker_idx + iter_idx)*nparams;
}

double* affineStepChunk::getLastParamPointer( unsigned int walker_idx ) {
  if (walker_idx >= nwalkers) return NULL;
  unsigned int iter_idxp1 = nsteps[walker_idx];
  if (iter_idxp1 == 0) return NULL;
  return steps + (niters*walker_idx + iter_idxp1 - 1)*nparams;
}

double* const affineStepChunk::getLastParamPointer( unsigned int walker_idx ) 
  const {
  if (walker_idx >= nwalkers) return NULL;
  unsigned int iter_idxp1 = nsteps[walker_idx];
  if (iter_idxp1 == 0) return NULL;
  return steps + (niters*walker_idx + iter_idxp1 - 1)*nparams;
}

void affineStepChunk::writeToStream( std::ostream& os ) const {
  //Do each walker in turn
  if (nparams == 0 || nwalkers == 0 || niters == 0) return;
  unsigned int curr_nsteps;
  double *ptr;
  for (unsigned int i = 0; i < nwalkers; ++i) {
    curr_nsteps = nsteps[i];
    for (unsigned int j = 0; j < curr_nsteps; ++j) {
      ptr = getParamPointer(i,j);
      os << ptr[0];
      for (unsigned int k = 1; k < nparams; ++k)
	os << " " << ptr[k];
      os << " " << logLike[i*niters+j] << " " << i << std::endl;
    }
  }
}

/////////////////////////////////////////////////////


affineChainSet::affineChainSet(unsigned int NWALKERS,
			       unsigned int NPARAMS) {
  nwalkers = NWALKERS;
  nparams  = NPARAMS;
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

void affineChainSet::clearPreserveLast() throw (affineExcept) {
  if (steps.size() == 0)
    throw affineExcept("affineChainSet","clearPreserveLast",
		       "No previous entries to preserve",1);
  affineStepChunk *newchunk = new affineStepChunk(nwalkers,1,nparams);
  newchunk->fillFromLastStep( *steps.back() );
  for (unsigned int i = 0; i < steps.size(); ++i)
    if (steps[i] != NULL) delete steps[i];
  steps.clear();
  steps.push_back(newchunk);
}

/*!
  \param[in] walker_idx Which walker to add the step to
  \param[in] p          Parameters to add
  \param[in] logLike    Log Likelihood of step
  \returns true if the step was added
 */
bool affineChainSet::addNewStep(unsigned int walker_idx, const paramSet& p, 
				double logLike) {
  return steps.back()->addStep(walker_idx,p,logLike);
}

void affineChainSet::addChunk(unsigned int sz) {
  affineStepChunk *chnkptr;
  chnkptr = new affineStepChunk(nwalkers,sz,nparams);
  steps.push_back( chnkptr );
}

/*!
  Implements stretch step

  \param[in] zval        Value of Z
  \param[in] idx1        Index of walker we are updating
  \param[in] idx2        Index of walker we are updating from
  \param[out] oldStep    Previous step from idx1
  \param[out] oldLogLike Previous Log Likelihood of idx1 walker
  \param[out] newStep    Proposed new step

  Takes the most recent available step for the specified walkers
 */

void affineChainSet::generateNewStep(double zval, unsigned int idx1,
				     unsigned int idx2, paramSet& oldStep,
				     double& oldLogLike, 
				     paramSet& newStep) const 
  throw (affineExcept) {
  unsigned int nsteps = steps.size();
  if (nsteps == 0) 
    throw affineExcept("affineChainSet","generateNewStep",
		       "No steps taken to generate from",1);
  if (idx1 >= nwalkers) 
    throw affineExcept("affineChainSet","generateNewStep",
		       "Invalid walker index for first walker",2);
  if (idx2 >= nwalkers) 
    throw affineExcept("affineChainSet","generateNewStep",
		       "Invalid walker index for second walker",4);

  if (oldStep.getNParams() != nparams) oldStep.setNParams(nparams);
  if (newStep.getNParams() != nparams) newStep.setNParams(nparams);

  //We have to figure out which Chunk to ask -- this is -not-
  // necessarily the one on the end, which could have no
  // steps in it!  And, to make things worse, we don't assume
  // that idx1 and idx2 have to come from the same chunk
  int chunknum1;
  for (chunknum1 = static_cast<int>(nsteps)-1; chunknum1 >= 0; --chunknum1)
    if (steps[chunknum1]->nsteps[idx1] > 0) break;
  if (chunknum1 < 0)     
    throw affineExcept("affineChainSet","generateNewStep",
		       "No steps taken to generate from for walker 1",2);
  int chunknum2;
  for (chunknum2 = static_cast<int>(nsteps)-1; chunknum2 >= 0; --chunknum2)
    if (steps[chunknum2]->nsteps[idx2] > 0) break;
  if (chunknum2 < 0)     
    throw affineExcept("affineChainSet","generateNewStep",
		       "No steps taken to generate from for walker 2",4);

  const affineStepChunk* const cptr1 = steps[chunknum1];
  const affineStepChunk* const cptr2 = steps[chunknum2];
  unsigned int csteps1 = cptr1->nsteps[idx1];

  //Get pointers to the parameter sets we are playing with
  //Old step (to update)
  double *ptr1 = cptr1->getLastParamPointer( idx1 );
  //Step to combine with old step
  double *ptr2 = cptr2->getLastParamPointer( idx2 );

  oldStep.setParamValues( nparams, ptr1 );
  oldLogLike = cptr1->getLogLike( idx1, csteps1-1 );

  //Compute stretch move
  double val, omz;
  omz = 1.0 - zval;
  for (unsigned int i = 0; i < nparams; ++i) {
    val = zval * ptr1[i] + omz*ptr2[i];
    newStep.setParamValue(i,val);
  }
}

void affineChainSet::getLastStep( unsigned int walker_idx,
				  paramSet& par, double& lglike ) const 
  throw (affineExcept) {
  if (steps.size() == 0) 
    throw affineExcept("affineChainSet","getLastStep",
		       "No steps taken",1);
  bool succ = steps.back()->getLastStep(walker_idx,par,lglike);
  if (!succ)
    throw affineExcept("affineChainSet","getLastStep",
		       "Error getting last step",2);
}

void affineChainSet::getStep(unsigned int chunkidx, 
			     unsigned int walker_idx, 
			     unsigned int iter_idx,
			     paramSet& par, double& lglike) 
  const throw (affineExcept) {
  if (chunkidx >= steps.size()) 
    throw affineExcept("affineChainSet","getStep",
		       "Chunk idx invalid",1);
  bool succ = steps[chunkidx]->getStep(walker_idx,iter_idx,par,lglike);
  if (!succ)
    throw affineExcept("affineChainSet","getStep",
		       "Error getting specified step",2);
}

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
    steps.push_back( new affineStepChunk( *other.steps[i] ) );
  return *this;
}
 
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

//Requires auxilliary storage equal to the number of walkers
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


void affineChainSet::getMaxLogLikeParam(double& maxval, paramSet& p) const {
  p.setNParams(nparams);
  if (steps.size() == 0) {
    maxval = std::numeric_limits<double>::quiet_NaN();
  }
  double currval;
  unsigned int firststep = 0;
  if (skipfirst) firststep = 1;
  steps[firststep]->getMaxLogLikeParam(maxval,p);

  paramSet cp(nparams);
  for (unsigned int i = firststep+1; i < steps.size(); ++i) {
    steps[i]->getMaxLogLikeParam(currval,cp);
    if (currval > maxval) {
      maxval = currval;
      p = cp;
    }
  }
}
 


void affineChainSet::getParamVector(unsigned int walkeridx,
				    unsigned int paramidx,
				    std::vector<double>& pvec) const 
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
  double *sgl_ptr;
  unsigned int firststep = 0;
  if (skipfirst) firststep = 1;
  for (unsigned int i = firststep; i < steps.size(); ++i) {
    chunkptr = steps[i];
    nciters = chunkptr->niters;
    ncpars  = chunkptr->nparams;
    chunksz = chunkptr->nsteps[walkeridx];
    sgl_ptr = chunkptr->steps+(walkeridx*nciters*ncpars);
    for (unsigned int j = 0; j < chunksz; ++j)
      pvec[ctr++] = sgl_ptr[j*ncpars+paramidx];
  }
}

/*!
  Returns parameter list averaged over all walkers.
  Only includes steps where all walkers have that step.
 */

void affineChainSet::getAverageParamVector(unsigned int paramidx,
					   std::vector<double>& pvec) const 
  throw (affineExcept) {
  if (paramidx >= nparams) 
    throw affineExcept("affineChainSet","getAverageParamVector",
		       "Specified paramidx exceeds number available",1);
  
  pvec.resize(getMinNIters());

  if (nwalkers == 1) {
    getParamVector(0,paramidx,pvec);
    return;
  }

  unsigned int ctr, nsteps, currsteps, nciters, ncpars;
  double cumval;
  double norm = 1.0/static_cast<double>(nwalkers);
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
    currsteps = 0;
    for (unsigned int j = 1; j < nwalkers; ++j) {
      currsteps = chunkptr->nsteps[i];
      if (currsteps < nsteps) nsteps = currsteps;
    }
    if (nsteps == 0) continue; //None

    for (unsigned int j = 0; j < nsteps; ++j) {
      cumval = chunkptr->steps[ j*ncpars + paramidx ]; //Walker 0
      for (unsigned int k = 1; k < nwalkers; ++k)
	cumval += chunkptr->steps[ (k*nciters+j)*ncpars + paramidx ];
      pvec[ctr++] = norm*cumval;
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
			 double& tau, std::vector<double>& X,
			 unsigned int L) 
  const throw (affineExcept) {

  // Compute tau directly only if tau < taumax, otherwise use the pairwise sum
  double invL = 1.0 / static_cast<double>(L);

  // Compute the mean of X ... 
  mean = X[0];
  for ( unsigned int i = 1; i < L; ++i) mean += X[i];
  mean *= invL;
  //  ... and subtract it away.
  for ( unsigned int i = 0; i <  L; ++i ) X[i] -= mean; 

  //If chunk is too small, stop recursion and just leave initial values alone
  if ( L < minfac*maxlag ) return 1;
   
  double C[maxlag+1]; 
  // Here, s=0 is the variance, s = maxlag is the last one computed.
  for ( unsigned int s = 0; s <= maxlag; ++s )  C[s] = 0.;  
  
  // Compute the autocovariance function . . . 
  unsigned int iMax = L - maxlag;
  double invImax = 1.0 / static_cast<double>(iMax);
  for ( unsigned int i = 0; i < iMax; ++i ) 
    // ...  first the inner products ...
    for ( unsigned int s = 0; s <= maxlag; ++s )
      C[s] += X[i]*X[i+s];                              
  // ...  then the normalization.
  for ( unsigned int s = 0; s <= maxlag; ++s ) C[s] *= invImax;
  
  // The "diffusion coefficient" is the sum of the autocovariances
  double D = C[0];   
  // The rest of the C[s] are double counted since C[-s] = C[s].
  for ( unsigned int s = 1; s <= maxlag; ++s ) D += 2*C[s];   
  // The standard error bar formula, if D were the complete sum.
  sigma = std::sqrt( D * invL );
  // A provisional estimate, since D is only part of the complete sum.
  tau   = D / C[0]; 
  
  // Stop if the D sum includes the given multiple of tau.
  // This is the self consistent window approach.
  if ( tau*winmult < maxlag ) return 0;    
  
  // If the provisional tau is so large that we don't think tau
  // is accurate, apply the acor procedure to the pairwase sums
  // of X.
  // The pairwise sequence is half the length (if L is even)
  unsigned int Lh = L/2;                          
  // The mean of the new sequence, to throw away.         
  double newMean;                                 
  unsigned int j1 = 0;
  unsigned int j2 = 1;
  for ( unsigned int i = 0; i < Lh; ++i ) {
    X[i] = X[j1] + X[j2];
    j1  += 2;
    j2  += 2; 
  }
  //Recursive call; note we ignore the return value!
  acor( newMean, sigma, tau, X, Lh); 
  
  // Reconstruct the fine time series numbers from the coarse series numbers.
  D     = 0.25*(sigma) * (sigma) * static_cast<double>(L);    
  tau   = D/C[0];                 // As before, but with a corrected D.
  sigma = std::sqrt( D * invL );  // As before, again.
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

  sigma = std::numeric_limits<double>::quiet_NaN();

  //Get the averaged parameters
  std::vector<double> pvec;
  getAverageParamVector(paramidx,pvec);

  //Make sure initial vector is long enough to allow for computation
  if ( pvec.size() < minfac*maxlag ) {
    succ = false;
    return std::numeric_limits<double>::quiet_NaN();
  }

  double tau;
  acor(mean,sigma,tau,pvec,pvec.size());
  succ = true;
  return tau;
}

/*!
  \returns True if the autocorrelation vector was computed for all params
 */
bool affineChainSet::getAcorVector(std::vector<double>& tau) const 
  throw (affineExcept) {
  double mn,sigma; //Throw away
  tau.resize(nparams);
  bool succ, indiv_succ;
  succ = true;
  for (unsigned int i = 0; i < nparams; ++i) {
    tau[i] = getAcor(i,mn,sigma,indiv_succ);
    succ &= indiv_succ;
  }
  return succ;
}


/*!
  \param[in] ignore Set to true where to ignore params.  tau is set to NaN
                    where this is true
  \returns True if the autocorrelation vector was computed for all params
            that were not ignored
 */
bool affineChainSet::getAcorVector(std::vector<double>& tau,
				   const std::vector<bool> ignore) const 
  throw (affineExcept) {
  if (ignore.size() < nparams)
    throw affineExcept("affineChainSet","getAcorVector",
		       "ignore is shorter than number of params",1);
  double mn,sigma; //Throw away
  tau.resize(nparams);
  bool succ, indiv_succ;
  succ = true;
  for (unsigned int i = 0; i < nparams; ++i) {
    if (! ignore[i]) {
      tau[i] = getAcor(i,mn,sigma,indiv_succ);
      succ &= indiv_succ;
    } else
      tau[i] = std::numeric_limits<double>::quiet_NaN();
  }
  return succ;
}

double affineChainSet::getParamMean(unsigned int paridx) const 
  throw (affineExcept) {
  if (paridx >= nparams)
    throw affineExcept("affineChainSet","getParamMean",
		       "Invalid parameter index",1);
  double mean = 0;
  unsigned int chunksz, nciters, ncpars, ctr = 0;
  const affineStepChunk* chunkptr;
  double *ptr;
  for (unsigned int i = 0; i < steps.size(); ++i) {
    chunkptr = steps[i];
    nciters = chunkptr->niters;
    ncpars  = chunkptr->nparams;
    for (unsigned int j = 0; j < nwalkers; ++j) {
      chunksz = chunkptr->nsteps[j];
      ptr = chunkptr->steps + j*nciters*ncpars;
      for (unsigned int k = 0; k < chunksz; ++k) {
	mean += ptr[j*ncpars+paridx];
	++ctr;
      }
    }
  }
  if (ctr == 0) return std::numeric_limits<double>::quiet_NaN();
  else return mean /= static_cast<double>(ctr);
}

/*!
  \param[in] idx Index of parameter
  \param[out] mean Mean value of parameter
  \param[out] var  Variance of parameter
  \param[out] lowlimit Lower limit of parameter. Returns Nan if non-computable
  \param[out] uplimit  Upper limit of parameter. Returns Nan if non-computable
  \param[in]  conflevel Confidence limit to use (def: 0.683, i.e. 1 sigma)

  Requires temporary storage equal to the total number of steps
*/
void affineChainSet::getParamStats(unsigned int paridx, double& mean,
				   double& var, double& lowlimit,
				   double& uplimit, double conflevel) const 
  throw (affineExcept) {
  if (paridx >= nparams)
    throw affineExcept("affineChainSet","getParamStats",
		       "Invalid parameter index",1);

  mean     = std::numeric_limits<double>::quiet_NaN();
  var      = std::numeric_limits<double>::quiet_NaN();
  lowlimit = std::numeric_limits<double>::quiet_NaN();
  uplimit  = std::numeric_limits<double>::quiet_NaN();

  //We will copy the full params into a vector for ease of sorting,
  // computation
  unsigned int sz = getNIters();
  if (sz == 0) return;
  if (skipfirst && steps.size() == 1) return;

  std::vector< double > pvec;
  pvec.resize(sz);
  
  //Stick in vector and add to mean as we go
  double val;
  unsigned int chunksz, ctr = 0, nciters, ncpars;
  const affineStepChunk* chunkptr;
  double *ptr;
  mean = 0.0;
  unsigned int firstidx = 0;
  if (skipfirst) firstidx = 1;
  for (unsigned int stepidx = firstidx; stepidx < steps.size(); ++stepidx) {
    //For each chunk
    chunkptr = steps[stepidx];
    nciters = chunkptr->niters;
    ncpars  = chunkptr->nparams;
    for (unsigned int i = 0; i < nwalkers; ++i) {
      //For each walker
      chunksz = chunkptr->nsteps[i];
      ptr     = chunkptr->steps + i*nciters*ncpars;
      for (unsigned int j = 0; j < chunksz; ++j) {
	//For each step in each chunk
	//ptr = chunkptr->getParamPointer(i,j);
	//val = ptr[ paridx ];
	val = ptr[ j*ncpars + paridx ];
	mean += val;
	pvec[ctr++] = val;
      }
    }
  }
  if (ctr != sz) 
    throw affineExcept("affineChainSet","getParamStats",
		       "Didn't get the number of expected steps",2);

  double norm = 1.0/static_cast<double>(sz);
  mean *= norm;

  //Now compute variance using two-pass algorithm
  double delta, sumdelta;
  var = sumdelta = 0.0;
  for (unsigned int i = 0; i < sz; ++i) {
    delta = pvec[i] - mean;
    sumdelta += delta;
    var += delta*delta;
  }
  var -= sumdelta*norm;
  var /= static_cast<double>(sz)-1.0;
  if (var == 0.0) return; //Can't compute limits
  
  //Compute limits
  if (conflevel < 0)
    throw affineExcept("affineChainSet","getParamStats",
		       "Invalid (negative) conf level",4);
  if (conflevel > 1)
    throw affineExcept("affineChainSet","getParamStats",
		       "Invalid (>1) conf level",8);

  std::sort( pvec.begin(), pvec.end() );
  //Now do the lower limit, rounding to get that index
  //I could add interpolation here, but not yet
  double lowfrac = 0.5*(1.0 - conflevel);
  int lowidx = static_cast<int>( lowfrac * sz );
  if (lowidx != 0) lowlimit = pvec[lowidx];
  //Upper limit is the same
  double upfrac = 1.0 - 0.5*conflevel;
  int upidx = static_cast<int>( upfrac * sz );
  if (upidx != 0) uplimit = pvec[upidx];
}

/*!
  Requires temporary storage equal to the full size of the chain
*/
void affineChainSet::makeCovMatrix(double** covmat) const {
  if (nparams == 0 || nwalkers == 0) return;
  unsigned int sz = getMinNIters();
  if (sz < 2) {
    for (unsigned int i = 0; i < nparams; ++i)
      for (unsigned int j = 0; j < nparams; ++j)
	covmat[i][j] = std::numeric_limits<double>::quiet_NaN();
    return;
  }

  //We copy everything into an array, nparams x sz
  double *tmpdata;
  tmpdata = new double[nparams*sz];
  unsigned int chunksz, ctr, ctrbase, nciters, ncpars;
  const affineStepChunk* chunkptr;
  double *ptr;
  //Take at most sz from each
  for (unsigned int i = 0; i < nparams; ++i) {
    ctrbase = sz*i;
    ctr = 0;
    for (unsigned int j = 0; j < steps.size(); ++j) if (ctr <= sz) {
	chunkptr = steps[j];
	nciters = chunkptr->niters;
	ncpars  = chunkptr->nparams;
	for (unsigned int k = 0; k < nwalkers; ++k) if (ctr <= sz) {
	    chunksz = chunkptr->nsteps[k];
	    ptr = chunkptr->steps + k*nciters*ncpars;
	    for (unsigned int m = 0; (m < chunksz) && (ctr <= sz); ++m) {
	      tmpdata[ctrbase+ctr] = ptr[j*ncpars+i];
	      ++ctr;
	    }
	  }
      }
  }

  //Now, mean subtract each
  double mean, norm;
  norm = 1.0 / static_cast<double>(sz);
  for (unsigned int i = 0; i < nparams; ++i) {
    ptr = tmpdata + sz*i;
    mean = ptr[0];
    for (unsigned int j = 1; j < sz; ++j)
      mean += ptr[j];
    mean *= norm;
    for (unsigned int j = 0; j < sz; ++j)
      ptr[j] -= mean;
  }

  //And compute covariances
  double norm2 = 1.0 / static_cast<double>(sz-1);
  double val, sumsq, *ptr2;
  for (unsigned int i = 0; i < nparams; ++i) {
    ptr = tmpdata + sz*i;
    //Diagonal
    val = ptr[0];
    sumsq = val*val;
    for (unsigned int j = 1; j < sz; ++j) {
      val = ptr[j];
      sumsq += val*val;
    }
    covmat[i][i] = sumsq * norm2;

    //Off diag, take advantage of symmetry
    for (unsigned int j = 0; j < i; ++j) {
      ptr2 = tmpdata + sz*j;
      sumsq = ptr[0]*ptr2[0];
      for (unsigned int k = 1; k < sz; ++k)
	sumsq += ptr[k]*ptr2[k];
      covmat[i][j] = covmat[j][i] = sumsq*norm2;
    }
  }

  delete[] tmpdata;
}


void affineChainSet::writeToFile( const std::string& outfile ) const 
  throw (affineExcept) {
  if (steps.size() == 0) return;
  std::ofstream ofs;
  ofs.open( outfile.c_str() );
  if (! ofs ) {
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
