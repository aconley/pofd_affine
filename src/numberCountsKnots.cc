//numberCountsKnots.cc
#include<iostream>
#include<cmath>
#include<iomanip>
#include<limits>
#include<cstdlib>
#include<ctime>
#include<sstream>

#include<global_settings.h>
#include<numberCountsKnots.h>

#include<utility.h>
#include<affineExcept.h>

numberCountsKnots::numberCountsKnots() {
  nknots = 0;
  knots = NULL;
  logknotvals = NULL;
  knotvals_loaded = false;
}

numberCountsKnots::numberCountsKnots( unsigned int NKNOTS ) :
  nknots(NKNOTS) {
  if (nknots > 0) {
    knots = new double[nknots];
    logknotvals = new double[nknots];
  } else
    knots = logknotvals = NULL;
  knotvals_loaded = false;
}  

numberCountsKnots::numberCountsKnots( const std::vector<double>& S ) {
  nknots = 0;
  knots = NULL;
  logknotvals = NULL;
  setKnotPositions(S);
  knotvals_loaded = false;
}

numberCountsKnots::numberCountsKnots( unsigned int n, const double* const S) {
  nknots = 0;
  knots = NULL;
  logknotvals = NULL;
  setKnotPositions(n,S);
  knotvals_loaded = false;
}

numberCountsKnots::numberCountsKnots( const numberCountsKnots& other ) {
  if ( this == &other ) return; //Self-copy
  nknots = 0;
  knots = NULL;
  logknotvals = NULL;
  setNKnots(other.nknots);
  for (unsigned int i = 0; i < nknots; ++i)
    knots[i] = other.knots[i];
  if (other.knotvals_loaded) {
    for (unsigned int i = 0; i < nknots; ++i)
      knots[i] = other.logknotvals[i];
  }
  knotvals_loaded = other.knotvals_loaded;
}

numberCountsKnots::~numberCountsKnots() {
  if (knots != NULL) delete[] knots;
  if (logknotvals != NULL) delete[] logknotvals;
}

void numberCountsKnots::setNKnots(unsigned int n) {
  if ( nknots == n ) return;
  if ( knots != NULL ) delete[] knots;
  if ( logknotvals != NULL ) delete[] logknotvals;
  if ( n > 0 ) {
    knots = new double[n];
    logknotvals = new double[n];
  } else {
    knots = logknotvals = NULL;
  }
  nknots = n;
  knotvals_loaded = false;
}

/*!
  \param[in] S Input knot positions
*/
void numberCountsKnots::setKnotPositions(const std::vector<double>& S) {
  unsigned int n = S.size();
  if (n != nknots) setNKnots(n);
  for (unsigned int i = 0; i < nknots; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsKnots","setKnots",
			 "Negative knot positions not allowed",1);
  for (unsigned int i = 0; i < nknots; ++i)
    knots[i] = S[i];
}

/*!
  \param[in] n Number of knots
  \param[in] S Input knot positions
*/
void numberCountsKnots::setKnotPositions(unsigned int n, 
					 const double* const S) {
  if (n != nknots) setNKnots(n);
  for (unsigned int i = 0; i < nknots; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsKnots","setKnots",
			 "Negative knot positions not allowed",1);
  for (unsigned int i = 0; i < nknots; ++i)
    knots[i] = S[i];
}

/*!
  \param[in] F Parameters to set in model
 */
void numberCountsKnots::setParams(const paramSet& F) {
  if (F.getNParams() <= 2)
    throw affineExcept("numberCountsKnots","setKnots",
		       "No knots present",1);
  if (nknots != (F.getNParams())) 
    throw affineExcept("numberCountsKnots","setKnots",
		       "Number of knot values different than expected",2);
  for (unsigned int i = 0; i < nknots; ++i)
    logknotvals[i] = pofd_mcmc::logfac*F[i]; //convert to base 2
  knotvals_loaded = true;
}

/*!
  \returns True if the model parameters are valid
 */
bool numberCountsKnots::isValid() const {
  if (nknots == 0) return false;
  if (!knotvals_loaded) return false;
  for (unsigned int i = 0; i < nknots; ++i)
    if ( std::isnan(knots[i]) ) return false;
  if ( knots[0] <= 0.0 ) return false;
  for (unsigned int i = 1; i < nknots; ++i)
    if (knots[i] <= knots[i-1] ) return false;
  for (unsigned int i = 0; i < nknots; ++i)
    if ( std::isnan(logknotvals[i]) ) return false;
  return true;
}

double numberCountsKnots::getMinFlux() const {
  if (nknots == 0) return std::numeric_limits<double>::quiet_NaN();
  return knots[0];
}

double numberCountsKnots::getMaxFlux() const {
  if (nknots == 0) return std::numeric_limits<double>::quiet_NaN();
  return knots[nknots-1];
}

void numberCountsKnots::sendSelf(MPI::Comm& comm, int dest) const {
  comm.Send(&nknots,1,MPI::UNSIGNED,dest,pofd_mcmc::NCKSENDNKNOTS);
  if (nknots != 0) {
    comm.Send(knots,nknots,MPI::DOUBLE,dest,pofd_mcmc::NCKSENDKNOTS);
    comm.Send(&knotvals_loaded,1,MPI::BOOL,dest,pofd_mcmc::NCKSENDKNOTSLOADED);
    if (knotvals_loaded)
      comm.Send(logknotvals,nknots,MPI::DOUBLE,dest,
		pofd_mcmc::NCKSENDLOGKNOTVALS);
  }
}

void numberCountsKnots::recieveCopy(MPI::Comm& comm, int src) {
  unsigned int n;
  comm.Recv(&n,1,MPI::UNSIGNED,src,pofd_mcmc::NCKSENDNKNOTS);
  if (n != 0) {
    if (n != nknots) setNKnots(n);
    comm.Recv(knots,nknots,MPI::DOUBLE,src,pofd_mcmc::NCKSENDKNOTS);
    comm.Recv(&knotvals_loaded,1,MPI::BOOL,src,pofd_mcmc::NCKSENDKNOTSLOADED);
    if (knotvals_loaded)
      comm.Recv(logknotvals,nknots,MPI::DOUBLE,src,
		pofd_mcmc::NCKSENDLOGKNOTVALS);
  }
}

bool numberCountsKnots::writeToStream(std::ostream& os) const {
  os << "Model parameters: " << std::endl;
  if (knotvals_loaded) {
    os << " " << std::left << std::setw(13) << "#Flux knot" << "  "
       << std::setw(13) << "Knot value" << std::endl;
    //Convert to log10 for output
    for (unsigned int i = 0; i < nknots; ++i)
      os << " " << std::left << std::setw(13) << knots[i] << "  "
	 << std::setw(13) << pofd_mcmc::ilogfac * logknotvals[i] << std::endl; 
  } else
    os << "Number of knots: " << nknots << std::endl;
  return true;
}

std::ostream& operator<<(std::ostream& os, const numberCountsKnots& b) {
  b.writeToStream(os);
  return os;
}

////////////////////////////////////////

initFileKnots::initFileKnots() : nknots(0), has_sigma(false), 
				 has_lower_limits(false),
				 has_upper_limits(false) {

  knotpos = knotval = sigma = lowlim = uplim = NULL;
  has_lowlim = has_uplim = NULL;

  //Set random number generator seed from time
  unsigned long long int seed;
  seed = static_cast<unsigned long long int>(time(NULL));
  seed += static_cast<unsigned long long int>(clock());
  rangen.setSeed(seed);
}

/*
  \param[in] flname File to read from
  \param[in] read_sigma      Read in (and require) knot sigmas
  \param[in] read_limits     Try to read limits; this will turn on require_sigma
  
  see initFileKnots::readFile for more details of the file format
*/
initFileKnots::initFileKnots(const std::string& flname, 
			     bool read_sigma, bool read_limits) :
  nknots(0), has_sigma(false), has_lower_limits(false), 
  has_upper_limits(false) {
  
  knotpos = knotval = sigma = lowlim = uplim = NULL;
  has_lowlim = has_uplim = NULL;

  //Set random number generator seed from time
  unsigned long long int seed;
  seed = static_cast<unsigned long long int>(time(NULL));
  seed += static_cast<unsigned long long int>(clock());
  rangen.setSeed(seed);

  readFile(flname, read_sigma, read_limits);
}

initFileKnots::~initFileKnots() {
  if (knotpos != NULL) delete[] knotpos;
  if (knotval != NULL) delete[] knotval;
  if (sigma != NULL) delete[] sigma;
  if (has_lowlim != NULL) delete[] has_lowlim;
  if (lowlim != NULL) delete[] lowlim;
  if (has_uplim != NULL) delete[] has_uplim;
  if (uplim != NULL) delete[] uplim;
}

/*!
  \param[in] flname File to read from
  \param[in] read_sigma      Read in (and require) knot sigmas
  \param[in] read_limits     Try to read limits; this will turn on require_sigma

  The file format is a bunch of lines of the form
  knotpos   knotval   [sigma [ lowlim [ uplim ]]
  So knotpos, knotval are always required
  sigma is optionally required if read_sigma is set
  lowlim and uplim may be present, and are looked for if read_limits is set.
    sigma must also be present, and the first element found is lowlim.
    If another is also found, it is interpreted as uplim
 */
void initFileKnots::readFile(const std::string& flname, 
			     bool read_sigma, bool read_limits) {
  if (read_limits) read_sigma = true;

  //Clear old data
  nknots = 0;
  if (knotpos != NULL) delete[] knotpos;
  if (knotval != NULL) delete[] knotval;
  if (sigma != NULL) delete[] sigma;
  if (has_lowlim != NULL) delete[] has_lowlim;
  if (lowlim != NULL) delete[] lowlim;
  if (has_uplim != NULL) delete[] has_uplim;
  if (uplim != NULL) delete[] uplim;
  knotpos = knotval = sigma = lowlim = uplim = NULL;
  has_lowlim = has_uplim = NULL;
  has_sigma = has_lower_limits = has_upper_limits = false;

  //Figure out how many elements we require
  unsigned int nreq = 2; //Pos, value
  if (read_sigma) nreq += 1; //Sigma -- read limits not required ever

  std::ifstream initfs( flname.c_str() );
  if (!initfs) {
    initfs.close();
    std::stringstream errmsg;
    errmsg << "Unable to open file:" << flname << std::endl;
    throw affineExcept("initFileKnots","readFile",errmsg.str(),1);
  }

  //Do the read into temporary vectors, then copy
  std::string line;
  std::vector<std::string> words;
  std::stringstream str;
  double currval;
  std::vector<double> kp, kv, ks, kl, ku;
  std::vector<bool> hl, hu;

  while (!initfs.eof()) {
    std::getline(initfs,line);
    if (line[0] == '#') continue; //Skip comments

    //Parse into words, stipping spaces
    utility::stringwords(line,words);
    if (words.size() == 0) continue; //Nothing on line (with spaces removed)
    if (words[0][0] == '#') continue; //Comment line
    if (words.size() < nreq) continue; //Has wrong number of entries
    str.str(words[0]); str.clear(); str >> currval;
    kp.push_back(currval);
    str.str(words[1]); str.clear(); str >> currval;
    kv.push_back(currval);
    if (read_sigma) {
      has_sigma = true;
      str.str(words[2]); str.clear(); str >> currval;
      ks.push_back(currval);
    }
    if (read_limits) { 
      if (words.size() > 3) {
	has_lower_limits = true;
	str.str(words[3]); str.clear(); str >> currval;
	hl.push_back(true);
	kl.push_back(currval);
	if (words.size() > 4) {
	  has_upper_limits = true;
	  str.str(words[4]); str.clear(); str >> currval;
	  hu.push_back(true);
	  ku.push_back(currval);
	} else {
	  hu.push_back(false);
	  ku.push_back(std::numeric_limits<double>::quiet_NaN());
	}
      } else {
	hl.push_back(false);
	kl.push_back( std::numeric_limits<double>::quiet_NaN() );
	hu.push_back(false);
	ku.push_back( std::numeric_limits<double>::quiet_NaN() );
      }
    }
  }
  initfs.close();
  
  nknots = kp.size();
  if (nknots == 0) {
    std::stringstream errstr;
    errstr << "No knot positions or values found in " << flname;
    throw affineExcept("initFileKnots","readFiles",errstr.str(),2);
  }

  //Copy into vectors
  knotpos = new double[nknots];
  for (unsigned int i = 0; i < nknots; ++i) knotpos[i] = kp[i];
  knotval = new double[nknots];
  for (unsigned int i = 0; i < nknots; ++i) knotval[i] = kv[i];
  if (has_sigma) {
    sigma = new double[nknots];
    for (unsigned int i = 0; i < nknots; ++i) sigma[i] = ks[i];
  }
  if (has_lower_limits) {
    has_lowlim = new bool[nknots];
    lowlim = new double[nknots];
    for (unsigned int i = 0; i < nknots; ++i) has_lowlim[i] = hl[i];
    for (unsigned int i = 0; i < nknots; ++i) lowlim[i] = kl[i];
  }
  if (has_upper_limits) {
    has_uplim = new bool[nknots];
    uplim = new double[nknots];
    for (unsigned int i = 0; i < nknots; ++i) has_uplim[i] = hu[i];
    for (unsigned int i = 0; i < nknots; ++i) uplim[i] = ku[i];
  }

  //Make sure lower/upper limits don't cross
  if (has_lower_limits && has_upper_limits)
    for (unsigned int i = 0; i < nknots; ++i) {
      if ( !(has_uplim[i] && has_lowlim[i]) ) continue;
      if (uplim[i] < lowlim[i]) {
	std::stringstream errstr;
	errstr << "Lower/Upper limits cross at index: " << i << std::endl;
	errstr << " Lower limit: " << lowlim[i] 
	       << " Upper limit: " << uplim[i];
	throw affineExcept("initFileKnots", "readFiles", errstr.str(), 3);
      }
      if ( (sigma[i] > 0.) && (uplim[i] == lowlim[i]) ) {
	std::stringstream errstr;
	errstr << "Lower/Upper limits meet at index: " << i 
	       << " but sigma is not zero" << std::endl;
	errstr << " Lower limit: " << lowlim[i] << " Upper limit: " << uplim[i]
	       << " sigma: " << sigma[i];
	throw affineExcept("initFileKnots", "readFiles", errstr.str(), 4);
      }
    }

  //Make sure that if sigmas is 0 then the mean value falls within
  // the range of any limits
  if (has_sigma && (has_lower_limits || has_upper_limits)) {
    for (unsigned int i = 0; i < nknots; ++i)
      if (sigma[i] > 0) {
	if (has_lower_limits && has_lowlim[i] && (knotval[i] < lowlim[i])) {
	  std::stringstream errstr;
	  errstr << "At knot " << i << " sigma is zero but mean value "
		 << knotval[i] << std::endl << " lies below lower limit "
		 << lowlim[i];
	  throw affineExcept("initFileKnots", "readFiles", errstr.str(), 5);
	}
	if (has_upper_limits && has_uplim[i] && (knotval[i] > uplim[i])) {
	  std::stringstream errstr;
	  errstr << "At knot " << i << " sigma is zero but mean value "
		 << knotval[i] << std::endl << " lies above upper limit "
		 << uplim[i];
	  throw affineExcept("initFileKnots", "readFiles", errstr.str(), 6);
	}
      }
  }
}

std::pair<double,double> initFileKnots::getKnot(unsigned int idx) const {
  if (nknots == 0)
    throw affineExcept("initFileKnots","getKnot",
		       "No knot information read in",1);
  if (idx >= nknots)
    throw affineExcept("initFileKnots","getKnot",
		       "Invalid knot index",2);
  return std::make_pair(knotpos[idx],knotval[idx]);
}

/*
  \param[out] kp Set to knot positions on output
 */
void initFileKnots::getKnotPos(std::vector<double>& kp) const {
  if (nknots == 0)
    throw affineExcept("initFileKnots","getKnotPos",
		       "No knot information read in",1);
  kp.resize(nknots);
  for (unsigned int i = 0; i < nknots; ++i)
    kp[i] = knotpos[i];
}

/*
  \param[in] idx Index of knot to get
  \returns Knot position at that index
*/
double initFileKnots::getKnotPos(unsigned int idx) const {
  if (nknots == 0)
    throw affineExcept("initFileKnots","getKnotPos",
		       "No knot information read in",1);
  if (idx >= nknots)
    throw affineExcept("initFileKnots","getKnotPos",
		       "Invalid knot index",2);
  return knotpos[idx];
}

/*
  \param[inout] model Modified on output; knot positions are set

  This will change the number of knots in the model if they
  don't match.
 */
void initFileKnots::getKnotPos(numberCountsKnots& model) const {
  if (nknots == 0)
    throw affineExcept("initFileKnots","getKnotPos",
		       "No knot information read in",1);
  model.setKnotPositions(nknots, knotpos);
}

/*
  \param[out] kp Set to knot positions on output
 */
void initFileKnots::getKnotVals(std::vector<double>& kv) const {
  if (nknots == 0)
    throw affineExcept("initFileKnots","getKnotVals",
		       "No knot information read in",1);
  kv.resize(nknots);
  for (unsigned int i = 0; i < nknots; ++i)
    kv[i] = knotval[i];
}

/*
  \param[in] idx Index of knot to get
  \returns Knot value at that index
*/
double initFileKnots::getKnotValue(unsigned int idx) const {
  if (nknots == 0)
    throw affineExcept("initFileKnots","getKnotValue",
		       "No knot information read in",1);
  if (idx >= nknots)
    throw affineExcept("initFileKnots","getKnotValue",
		       "Invalid knot index",2);
  return knotval[idx];
}


/*
  This only fills the first nknots parameters
 */
void initFileKnots::getParams(paramSet& p) const {
  if (nknots == 0)
    throw affineExcept("initFileKnots","getParams",
		       "No information loaded",1);
  if (p.getNParams() < nknots)
    throw affineExcept("initFileKnots","getParams",
		       "Not enough space in provided paramSet",2);
  for (unsigned int i = 0; i < nknots; ++i)
    p[i] = knotval[i];
}

/*
  This only fills the first nknots parameters
 */
void initFileKnots::generateRandomKnotValues(paramSet& p) const {
  const unsigned int maxiters = 1000; //Maximum number of generation attempts
  if (nknots == 0)
    throw affineExcept("initFileKnots","generateRandomKnotValues",
		       "No knot information read in",1);
    
  //Make sure p is big enough; don't resize, complain
  if (p.getNParams() < nknots)
    throw affineExcept("initFileKnots","generateRandomKnotValues",
		       "Not enough space in provided paramSet",2);

  //Deal with simple case -- everything fixed
  //So just return central values
  if (!has_sigma) {
    for (unsigned int i = 0; i < nknots; ++i)
      p[i] = knotval[i];
    return;
  }

  //Now we have at least some sigmas
  //The simple case is if there are no limits.  If there are, we will
  // have to do trials.
  if (!(has_lower_limits || has_upper_limits)) {
    for (unsigned int i = 0; i < nknots; ++i)
      if (sigma[i] > 0)
	p[i] = rangen.gauss() * sigma[i] + knotval[i];
      else
	p[i] = knotval[i];
  } else {
    //Both sigmas and limits
    bool goodval;
    double trialval;
    unsigned int iters;
    for (unsigned int i = 0; i < nknots; ++i) {
      if (sigma[i] > 0) {
	//Some sanity checks
	if (has_lowlim[i] && (lowlim[i] > knotval[i]+sigma[i]*4.0)) {
	  std::stringstream errstr;
	  errstr << "Lower limit is too far above central value; will take too"
		 << " long to " << std::endl << "generate value for param idx: "
		 << i;
	  throw affineExcept("initFileKnots","generateRandomKnotValues",
			     errstr.str(),4);
	}
	if (has_uplim[i] && (uplim[i] < knotval[i]-sigma[i]*4.0)) {
	  std::stringstream errstr;
	  errstr << "Upper limit is too far below central value; will take too"
		 << " long to " << std::endl << "generate value for param idx: "
		 << i;
	  throw affineExcept("initFileKnots","generateRandomKnotValues",
			     errstr.str(),8);
	}
	
	if (has_lowlim[i] && has_uplim[i]) {
	  //If range between them is much smaller than
	  // sigma, just chose uniformly to save time
	  double rng = uplim[i] - lowlim[i];
	  if (rng < 1e-4*sigma[i]) {
	    p[i] = lowlim[i] + rng*rangen.doub();
	  } else {
	    //Trial
	    goodval = false;
	    iters = 0;
	    while (!goodval) {
	      if (iters >= maxiters) {
		std::stringstream errstr;
		errstr << "Failed to generate acceptable value for param "
		       << i << " after " << iters << " attempts";
		throw affineExcept("initFileKnots","generateRandomKnotValues",
				   errstr.str(),16);
	      }
	      trialval = rangen.gauss() * sigma[i] + knotval[i];
	      if ( (trialval >= lowlim[i]) && (trialval <= uplim[i]) ) 
		goodval = true;
	      ++iters;
	    }
	    p[i] = trialval;
	  }
	} else if (has_lowlim[i]) {
	  //Lower limit only
	  goodval = false;
	  iters = 0;
	  while (!goodval) {
	    if (iters >= maxiters) {
	      std::stringstream errstr;
	      errstr << "Failed to generate acceptable value for param "
		     << i << " after " << iters << " attempts";
	      throw affineExcept("initFileKnots","generateRandomKnotValues",
				 errstr.str(),16);
	    }
	    trialval = rangen.gauss() * sigma[i] + knotval[i];
	    if (trialval >= lowlim[i]) goodval = true;
	    ++iters;
	  }
	  p[i] = trialval;
	} else if (has_uplim[i]) {
	  //Upper limit only
	  goodval = false;
	  iters = 0;
	  while (!goodval) {
	    if (iters >= maxiters) {
	      std::stringstream errstr;
	      errstr << "Failed to generate acceptable value for param "
		     << i << " after " << iters << " attempts";
	      throw affineExcept("initFileKnots","generateRandomKnotValues",
				 errstr.str(),16);
	    }
	    trialval = rangen.gauss() * sigma[i] + knotval[i];
	    if (trialval <= uplim[i]) goodval = true;
	    ++iters;
	  }
	  p[i] = trialval;
	} else {
	  //No limit, easy cakes
	  p[i] = rangen.gauss() * sigma[i] + knotval[i];
	}
      } else {
	//Sigma is 0.  The read operation makes sure that, in this case,
	// the limits include the mean.
	p[i] = knotval[i];
      }
    }
  }

}

double initFileKnots::getKnotSigma(unsigned int idx) const {
  if (idx >= nknots)
    throw affineExcept("initFileKnots","getKnotSigma","Invalid knot index",1);
  if (!has_sigma) return std::numeric_limits<double>::quiet_NaN();
  return sigma[idx];
}

/*
  \param[in] idx Knot index
  \returns Whether knot is fixed
 */
bool initFileKnots::isKnotFixed(unsigned int idx) const {
  if (idx >= nknots)
    throw affineExcept("initFileKnots","isKnotFixed","Invalid knot index",1);
  if (!has_sigma) return false;
  if (sigma[idx] == 0) return true;
  return false;
}

/*
  \param[in] idx Knot index
  \returns Whether knot has a lower limit
 */
bool initFileKnots::knotHasLowerLimit(unsigned int idx) const {
  if (!has_lower_limits) return false;
  if (idx >= nknots)
    throw affineExcept("initFileKnots","knotHasLowerLimit",
		       "Invalid knot index",1);
  return has_lowlim[idx];
}

/*
  \param[in] idx Knot index
  \returns Lower limit on knot, or NaN if none
 */
double initFileKnots::getLowerLimit(unsigned int idx) const {
  if (!has_lower_limits) return std::numeric_limits<double>::quiet_NaN();
  if (idx >= nknots)
    throw affineExcept("initFileKnots","getLowerLimit",
		       "Invalid knot index",1);
  if (!has_lowlim[idx]) return std::numeric_limits<double>::quiet_NaN();
  return lowlim[idx];
}

/*
  \param[in] idx Knot index
  \returns Whether knot has a upper limit
 */
bool initFileKnots::knotHasUpperLimit(unsigned int idx) const {
  if (!has_upper_limits) return false;
  if (idx >= nknots)
    throw affineExcept("initFileKnots","knotHasUpperLimit",
		       "Invalid knot index",1);
  return has_uplim[idx];
}


/*
  \param[in] idx Knot index
  \returns Upper limit on knot, or NaN if none
 */
double initFileKnots::getUpperLimit(unsigned int idx) const {
  if (!has_upper_limits) return std::numeric_limits<double>::quiet_NaN();
  if (idx >= nknots)
    throw affineExcept("initFileKnots","getUpperLimit",
		       "Invalid knot index",1);
  if (!has_uplim[idx]) return std::numeric_limits<double>::quiet_NaN();
  return uplim[idx];
}


bool initFileKnots::isValid(const paramSet& p) const {
  if (! (has_lower_limits || has_upper_limits) ) return true;
  if (p.getNParams() < nknots)
    throw affineExcept("initFileKnots","isValid",
		       "Not enough params in paramSet to test validity",1);
  double val;
  for (unsigned int i = 0; i < nknots; ++i) {
    val = p[i];
    if (has_lowlim[i] && (val < lowlim[i])) return false;
    if (has_uplim[i] && (val < uplim[i])) return false;
  }
  return true;
}

void initFileKnots::sendSelf(MPI::Comm& comm, int dest) const {
  comm.Send(&nknots,1,MPI::UNSIGNED,dest,pofd_mcmc::IFKSENDNKNOTS);
  if (nknots != 0) {
    comm.Send(knotpos,nknots,MPI::DOUBLE,dest,pofd_mcmc::IFKSENDKNOTPOS);
    comm.Send(knotval,nknots,MPI::DOUBLE,dest,pofd_mcmc::IFKSENDKNOTVAL);
    comm.Send(&has_sigma,1,MPI::BOOL,dest,pofd_mcmc::IFKHASSIGMA);
    if (has_sigma)
      comm.Send(sigma,nknots,MPI::DOUBLE,dest,pofd_mcmc::IFKSENDSIGMA);
    comm.Send(&has_lower_limits,1,MPI::BOOL,dest,pofd_mcmc::IFKHASLOWERLIMITS);
    if (has_lower_limits) {
      comm.Send(has_lowlim,nknots,MPI::BOOL,dest,pofd_mcmc::IFKSENDHASLOWLIM);
      comm.Send(lowlim,nknots,MPI::DOUBLE,dest,pofd_mcmc::IFKSENDLOWLIM);
    }
    comm.Send(&has_upper_limits,1,MPI::BOOL,dest,pofd_mcmc::IFKHASUPPERLIMITS);
    if (has_lower_limits) {
      comm.Send(has_lowlim,nknots,MPI::BOOL,dest,pofd_mcmc::IFKSENDHASLOWLIM);
      comm.Send(lowlim,nknots,MPI::DOUBLE,dest,pofd_mcmc::IFKSENDLOWLIM);
    }
  }
}

void initFileKnots::recieveCopy(MPI::Comm& comm, int src) {
  //Delete everything for simplicity
  if (knotpos != NULL) delete[] knotpos;
  if (knotval != NULL) delete[] knotval;
  if (sigma != NULL) delete[] sigma;
  if (has_lowlim != NULL) delete[] has_lowlim;
  if (lowlim != NULL) delete[] lowlim;
  if (has_uplim != NULL) delete[] has_uplim;
  if (uplim != NULL) delete[] uplim;
  knotpos = knotval = sigma = lowlim = uplim = NULL;
  has_lowlim = has_uplim = NULL;

  comm.Recv(&nknots,1,MPI::UNSIGNED,src,pofd_mcmc::IFKSENDNKNOTS);
  if (nknots > 0) {
    knotpos = new double[nknots];
    comm.Recv(knotpos,nknots,MPI::DOUBLE,src,pofd_mcmc::IFKSENDKNOTPOS);
    knotval = new double[nknots];
    comm.Recv(knotval,nknots,MPI::DOUBLE,src,pofd_mcmc::IFKSENDKNOTVAL);
    comm.Recv(&has_sigma,1,MPI::BOOL,src,pofd_mcmc::IFKHASSIGMA);
    if (has_sigma) {
      sigma = new double[nknots];
      comm.Recv(sigma,nknots,MPI::DOUBLE,src,pofd_mcmc::IFKSENDSIGMA);
    }
    comm.Recv(&has_lower_limits,1,MPI::BOOL,src,pofd_mcmc::IFKHASLOWERLIMITS);
    if (has_lower_limits) {
      has_lowlim = new bool[nknots];
      comm.Recv(has_lowlim,nknots,MPI::BOOL,src,pofd_mcmc::IFKSENDHASLOWLIM);
      lowlim = new double[nknots];
      comm.Recv(lowlim,nknots,MPI::DOUBLE,src,pofd_mcmc::IFKSENDLOWLIM);
    }
    comm.Recv(&has_upper_limits,1,MPI::BOOL,src,pofd_mcmc::IFKHASUPPERLIMITS);
    if (has_upper_limits) {
      has_uplim = new bool[nknots];
      comm.Recv(has_uplim,nknots,MPI::BOOL,src,pofd_mcmc::IFKSENDHASUPLIM);
      uplim = new double[nknots];
      comm.Recv(uplim,nknots,MPI::DOUBLE,src,pofd_mcmc::IFKSENDUPLIM);
    }
  }
}
