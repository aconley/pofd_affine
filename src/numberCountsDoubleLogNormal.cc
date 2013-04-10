#include<cmath>
#include<iomanip>
#include<sstream> 

#include "../include/numberCountsDoubleLogNormal.h"
#include "../include/global_settings.h"
#include "../include/affineExcept.h"

//Function to pass to GSL integrator
/*! \brief Evaluates flux1^power1 * exp(const1*mu + const2*sigma^2) dN/dS1 */
static double evalPowfNDoubleLogNormal(double,void*); 

const unsigned int numberCountsDoubleLogNormal::nvarr = 17;

numberCountsDoubleLogNormal::numberCountsDoubleLogNormal() : 
  nknots(0), knots(NULL), logknots(NULL), logknotvals(NULL), splinelog(NULL), 
  nsigmaknots(0), sigmaknots(NULL), sigmavals(NULL), sigmainterp(NULL), 
  noffsetknots(0), offsetknots(NULL), offsetvals(NULL), offsetinterp(NULL),
  knots_valid(false), sigmas_valid(false), offsets_valid(false),
  knotpos_loaded(false), sigmapos_loaded(false), offsetpos_loaded(false),
  knotvals_loaded(false), sigmavals_loaded(false), offsetvals_loaded(false),
  nRWork(0), RWorkValid(NULL), RWork1(NULL), RWork2(NULL), RWork3(NULL) {
  
  acc = gsl_interp_accel_alloc();
  accsigma = gsl_interp_accel_alloc();
  accoffset = gsl_interp_accel_alloc();
  gsl_work = gsl_integration_workspace_alloc(1000);
  varr = new void*[nvarr];
}

/*!
  \param[in] NKNOTS Number of model knots for band 1 counts model
  \param[in] NSIGMA Number of model knots for band 2 color model sigma
  \param[in] NOFFSET Number of model knots for band 2 color model offset
*/
numberCountsDoubleLogNormal::numberCountsDoubleLogNormal(unsigned int NKNOTS,
							 unsigned int NSIGMA,
							 unsigned int NOFFSET) :
  nknots(NKNOTS), nsigmaknots(NSIGMA), noffsetknots(NOFFSET), 
  knots_valid(false), sigmas_valid(false), offsets_valid(false),
  knotpos_loaded(false), sigmapos_loaded(false), offsetpos_loaded(false),
  knotvals_loaded(false), sigmavals_loaded(false), offsetvals_loaded(false),
  nRWork(0), RWorkValid(NULL), RWork1(NULL), RWork2(NULL), RWork3(NULL) {

  if (NKNOTS > 0) {
    knots = new double[NKNOTS];
    logknots = new double[NKNOTS]; 
    logknotvals = new double[NKNOTS];
    if (NKNOTS <= 2)
      throw affineExcept("numberCountsDoubleLogNormal","constructor",
			 "nknots must be > 2",1);
    splinelog = gsl_spline_alloc(gsl_interp_cspline,
				 static_cast<size_t>(NKNOTS));
    for (unsigned int i = 0; i < NKNOTS; ++i)
      knots[i] = std::numeric_limits<double>::quiet_NaN();
    for (unsigned int i = 0; i < NKNOTS; ++i)
      logknots[i] = std::numeric_limits<double>::quiet_NaN();
    for (unsigned int i = 0; i < NKNOTS; ++i)
      logknotvals[i] = std::numeric_limits<double>::quiet_NaN();
  } else {
    knots = logknots = logknotvals = NULL;
    splinelog = NULL;
  }
  acc = gsl_interp_accel_alloc();

  sigmainterp = NULL;
  if (NSIGMA > 0) {
    sigmaknots = new double[NSIGMA];
    sigmavals = new double[NSIGMA];
    for (unsigned int i = 0; i < NSIGMA; ++i)
      sigmaknots[i] = std::numeric_limits<double>::quiet_NaN();
    for (unsigned int i = 0; i < NSIGMA; ++i)
      sigmavals[i] = std::numeric_limits<double>::quiet_NaN();
    if (NSIGMA == 2)
      sigmainterp = gsl_interp_alloc(gsl_interp_linear,
				     static_cast<size_t>(NSIGMA));
    else if (NSIGMA > 2)
      sigmainterp = gsl_interp_alloc(gsl_interp_cspline,
				     static_cast<size_t>(NSIGMA));
  } else sigmaknots = sigmavals = NULL;
  accsigma = gsl_interp_accel_alloc();

  offsetinterp = NULL;
  if (NOFFSET > 0) {
    offsetknots = new double[NOFFSET];
    offsetvals = new double[NOFFSET];
    for (unsigned int i = 0; i < NOFFSET; ++i)
      offsetknots[i] =  std::numeric_limits<double>::quiet_NaN();
    for (unsigned int i = 0; i < NOFFSET; ++i)
      offsetvals[i] =  std::numeric_limits<double>::quiet_NaN();
    if (NOFFSET == 2)
	offsetinterp = gsl_interp_alloc(gsl_interp_linear,
					static_cast<size_t>(NOFFSET));
    else if (NOFFSET > 2)
      offsetinterp = gsl_interp_alloc(gsl_interp_cspline,
				      static_cast<size_t>(NOFFSET));
  } else offsetknots = offsetvals = NULL;
  accoffset = gsl_interp_accel_alloc();

  gsl_work = gsl_integration_workspace_alloc(1000);
  varr = new void*[nvarr];
}
  
/*!
  \param[in] KNOTS Model knots positions for band 1 counts model
  \param[in] SIGMAS Model knots positions for band 2 color model sigma
  \param[in] OFFSETS Model knots positions for band 2 color model offset
*/
numberCountsDoubleLogNormal::
numberCountsDoubleLogNormal(const std::vector<double>& KNOTS,
			    const std::vector<double>& SIGMAS,
			    const std::vector<double>& OFFSETS) :
  knots_valid(false), sigmas_valid(false), offsets_valid(false),
  knotpos_loaded(false), sigmapos_loaded(false), offsetpos_loaded(false),
  knotvals_loaded(false), sigmavals_loaded(false), offsetvals_loaded(false),
  nRWork(0), RWorkValid(NULL), RWork1(NULL), RWork2(NULL), RWork3(NULL) {

  knots = 0; knots = logknots = logknotvals = NULL; 
  splinelog = NULL; 
  setKnotPositions(KNOTS);
  acc = gsl_interp_accel_alloc();

  nsigmaknots = 0; sigmaknots = sigmavals = NULL; sigmainterp = NULL;
  setSigmaPositions(SIGMAS);
  accsigma = gsl_interp_accel_alloc();

  noffsetknots = 0; offsetknots = offsetvals = NULL; offsetinterp = NULL;
  setOffsetPositions(OFFSETS);
  accoffset = gsl_interp_accel_alloc();

  gsl_work = gsl_integration_workspace_alloc(1000);
  varr = new void*[nvarr];
}

/*!
  \param[in] nk Number of elements in K
  \param[in] K Model knots positions for band 1 counts model
  \param[in] ns Number of elements in S
  \param[in] S Model knots positions for band 2 color model sigma  
  \param[in] no Number of elements in O
  \param[in] O Model knots positions for band 2 color model offset
*/
numberCountsDoubleLogNormal::
numberCountsDoubleLogNormal(unsigned int nk, const double* const K,
			    unsigned int ns, const double* const S,
			    unsigned int no, const double* const O) :
  knots_valid(false), sigmas_valid(false), offsets_valid(false),
  knotpos_loaded(false), sigmapos_loaded(false), offsetpos_loaded(false),
  knotvals_loaded(false), sigmavals_loaded(false), offsetvals_loaded(false),
  nRWork(0), RWorkValid(NULL), RWork1(NULL), RWork2(NULL), RWork3(NULL) {

  knots = 0; knots = logknots = logknotvals = NULL; splinelog = NULL;
  setKnotPositions(nk,K);
  acc = gsl_interp_accel_alloc();

  nsigmaknots = 0; sigmaknots = sigmavals = NULL; sigmainterp = NULL;
  setSigmaPositions(ns,S);
  accsigma = gsl_interp_accel_alloc();

  noffsetknots = 0; offsetknots = offsetvals = NULL; offsetinterp = NULL;
  setOffsetPositions(no,O);
  accoffset = gsl_interp_accel_alloc();

  gsl_work = gsl_integration_workspace_alloc(1000);
  varr = new void*[nvarr];
}

/*!
  \param[in] other Instance to copy from
*/
numberCountsDoubleLogNormal::
numberCountsDoubleLogNormal(const numberCountsDoubleLogNormal& other) {
  if (this == &other) return; //Self-copy
  nknots = nsigmaknots = noffsetknots = 0;
  knots = logknots = logknotvals = sigmaknots = 
    sigmavals = offsetknots = offsetvals = NULL;
  acc = accsigma = accoffset = NULL;
  splinelog = NULL;
  sigmainterp = offsetinterp = NULL;
  nRWork = 0;
  RWorkValid = NULL;
  RWork1 = RWork2 = RWork3 = NULL;
  knots_valid = sigmas_valid = offsets_valid = false;
  knotpos_loaded = sigmapos_loaded = offsetpos_loaded = false;
  knotvals_loaded = sigmavals_loaded = offsetvals_loaded = false;

  acc = gsl_interp_accel_alloc();
  accsigma = gsl_interp_accel_alloc();
  accoffset = gsl_interp_accel_alloc();

  //Knot positions
  setNKnots(other.nknots);
  if (other.knotpos_loaded) {
    for (unsigned int i = 0; i < nknots; ++i)
      knots[i] = other.knots[i];
    for (unsigned int i = 0; i < nknots; ++i)
      logknots[i] = other.logknots[i];
  }
  knotpos_loaded = other.knotpos_loaded;
  setNSigmas(other.nsigmaknots);
  if (other.sigmapos_loaded)
    for (unsigned int i = 0; i < nsigmaknots; ++i)
      sigmaknots[i] = other.sigmaknots[i];
  sigmapos_loaded = other.sigmapos_loaded;
  setNOffsets(other.noffsetknots);
  if (other.offsetpos_loaded)
    for (unsigned int i = 0; i < noffsetknots; ++i)
      offsetknots[i] = other.offsetknots[i];
  offsetpos_loaded = other.offsetpos_loaded;

  //Knot values
  if (other.knotvals_loaded) {
    for (unsigned int i = 0; i < nknots; ++i)
      logknotvals[i] = other.logknotvals[i];
    knotvals_loaded = other.knotvals_loaded;
    if (nknots > 2)
      gsl_spline_init(splinelog, logknots, logknotvals,
		      static_cast<size_t>(other.nknots));
    checkKnotsValid();
  }
  if (other.sigmavals_loaded) {
    for (unsigned int i = 0; i < nsigmaknots; ++i)
      sigmavals[i] = other.sigmavals[i];
    sigmavals_loaded = other.sigmavals_loaded;
    if (nsigmaknots > 1)
      gsl_interp_init(sigmainterp, sigmaknots, sigmavals,
		      static_cast<size_t>(other.nsigmaknots));
    checkSigmasValid();
  }
  if (other.offsetvals_loaded) {
    for (unsigned int i = 0; i < noffsetknots; ++i)
      offsetvals[i] = other.offsetvals[i];
    offsetvals_loaded = other.offsetvals_loaded;
    if (noffsetknots > 1)
      gsl_interp_init(offsetinterp, offsetknots, offsetvals,
		      static_cast<size_t>(other.noffsetknots));
    checkOffsetsValid();
  }

  gsl_work = gsl_integration_workspace_alloc(1000);
  varr = new void*[nvarr];
}

numberCountsDoubleLogNormal::~numberCountsDoubleLogNormal() {
  if (knots != NULL) delete[] knots;
  if (logknots != NULL) delete[] logknots;
  if (logknotvals != NULL) delete[] logknotvals;
  if (splinelog != NULL) gsl_spline_free(splinelog);
  gsl_interp_accel_free(acc);

  if (sigmaknots != NULL) delete[] sigmaknots;
  if (sigmavals != NULL) delete[] sigmavals;
  if (sigmainterp != NULL) gsl_interp_free(sigmainterp);
  gsl_interp_accel_free(accsigma);

  if (offsetknots != NULL) delete[] offsetknots;
  if (offsetvals != NULL) delete[] offsetvals;
  if (offsetinterp != NULL) gsl_interp_free(offsetinterp);
  gsl_interp_accel_free(accoffset);

  gsl_integration_workspace_free(gsl_work);
  delete[] varr;

  if (RWork1 != NULL) delete[] RWork1;
  if (RWork2 != NULL) delete[] RWork2;
  if (RWork3 != NULL) delete[] RWork3;
  if (RWorkValid != NULL) delete[] RWorkValid;
}

/*!
  \param[in] other Instance to copy from
*/
numberCountsDoubleLogNormal& numberCountsDoubleLogNormal::
operator=(const numberCountsDoubleLogNormal& other) {
  if (this == &other) return *this; //Self-copy

  knots_valid = sigmas_valid = offsets_valid = false;
  knotpos_loaded = sigmapos_loaded = offsetpos_loaded = false;
  knotvals_loaded = sigmavals_loaded = offsetvals_loaded = false;

  //Allocate space; will also delete if other is empty
  setNKnots(other.nknots);
  if (other.knotpos_loaded) {
    for (unsigned int i = 0; i < nknots; ++i)
      knots[i] = other.knots[i];
    for (unsigned int i = 0; i < nknots; ++i)
      logknots[i] = other.logknots[i];
    knotpos_loaded = other.knotpos_loaded;
  }
  setNSigmas(other.nsigmaknots);
  if (other.sigmapos_loaded) {
    for (unsigned int i = 0; i < nsigmaknots; ++i)
      sigmaknots[i] = other.sigmaknots[i];
    sigmapos_loaded = other.sigmapos_loaded;
  }
  setNOffsets(other.noffsetknots);
  if (other.offsetpos_loaded) {
    for (unsigned int i = 0; i < noffsetknots; ++i)
      offsetknots[i] = other.offsetknots[i];
    offsetpos_loaded = other.offsetpos_loaded;
  }

  //Knot values
  if (other.knotvals_loaded) {
    for (unsigned int i = 0; i < nknots; ++i)
      logknotvals[i] = other.logknotvals[i];
    knotvals_loaded = other.knotvals_loaded;
    if (nknots > 2)
      gsl_spline_init(splinelog, logknots, logknotvals,
		      static_cast<size_t>(other.nknots));
    checkKnotsValid();
  }
  if (other.sigmavals_loaded) {
    for (unsigned int i = 0; i < nsigmaknots; ++i)
      sigmavals[i] = other.sigmavals[i];
    sigmavals_loaded = other.sigmavals_loaded;
    if (nsigmaknots > 1)
      gsl_interp_init(sigmainterp, sigmaknots, sigmavals,
		      static_cast<size_t>(other.nsigmaknots));
    checkSigmasValid();
  }
  if (other.offsetvals_loaded) {
    for (unsigned int i = 0; i < noffsetknots; ++i)
      offsetvals[i] = other.offsetvals[i];
    offsetvals_loaded = other.offsetvals_loaded;
    if (noffsetknots > 1)
      gsl_interp_init(offsetinterp, offsetknots, offsetvals,
		      static_cast<size_t>(other.noffsetknots));
    checkOffsetsValid();
  }

  return *this;
}

/*!
  \param[in] n Number of knots to set
*/
void numberCountsDoubleLogNormal::setNKnots(unsigned int n) {
  if (nknots == n) return;

  if (knots != NULL) delete[] knots;
  if (logknotvals != NULL) delete[] logknotvals;
  if (logknots != NULL) delete[] logknots;
  if (splinelog != NULL) gsl_spline_free(splinelog);
  knots_valid = knotpos_loaded = knotvals_loaded = false;

  if (n > 0) {
    knots = new double[n];
    logknots = new double[n];
    logknotvals = new double[n];
  } else {
    knots = logknots = logknotvals = NULL;
  }
  if (n > 2)
    splinelog=gsl_spline_alloc(gsl_interp_cspline,
			       static_cast<size_t>(n));
  else 
    splinelog=NULL;
  nknots = n;
}

/*!
  \param[in] n Number of sigmas
*/
void numberCountsDoubleLogNormal::setNSigmas(unsigned int n) {
  if (nsigmaknots == n) return;

  if (sigmaknots != NULL) delete[] sigmaknots;
  if (sigmavals != NULL) delete[] sigmavals;
  if (sigmainterp != NULL) gsl_interp_free(sigmainterp);
  sigmas_valid = sigmapos_loaded = sigmavals_loaded = false;

  if (n > 0) {
    sigmaknots = new double[n];
    sigmavals = new double[n];
  } else {
    sigmaknots = sigmavals = NULL;
  }
  sigmainterp = NULL;
  if (n == 2)
      sigmainterp=gsl_interp_alloc(gsl_interp_linear,
				   static_cast<size_t>(n));
  else if (n > 2)
    sigmainterp=gsl_interp_alloc(gsl_interp_cspline,
				 static_cast<size_t>(n));
  nsigmaknots = n;
}

/*!
  \param[in] n Number of offsets
*/
void numberCountsDoubleLogNormal::setNOffsets(unsigned int n) {
  if (noffsetknots == n) return;

  if (offsetknots != NULL) delete[] offsetknots;
  if (offsetvals != NULL) delete[] offsetvals;
  if (offsetinterp != NULL) gsl_interp_free(offsetinterp);
  offsets_valid = offsetpos_loaded = offsetvals_loaded = false;

  if (n > 0) {
    offsetknots = new double[n];
    offsetvals = new double[n];
  } else {
    offsetknots = offsetvals = NULL;
  }
  offsetinterp = NULL;
  if (n == 2) 
    offsetinterp=gsl_interp_alloc(gsl_interp_linear,
				  static_cast<size_t>(n));
  else if (n > 2)
    offsetinterp=gsl_interp_alloc(gsl_interp_cspline,
				  static_cast<size_t>(n));
  noffsetknots = n;
}

/*!
  \param[in] S Input knot positions for band 1 model
*/
void numberCountsDoubleLogNormal::
setKnotPositions(const std::vector<double>& S) {
  unsigned int n = S.size();
  setNKnots(n); //Also invalidates everything
  for (unsigned int i = 0; i < nknots; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsDoubleLogNormal","setKnotPositions",
			 "Non-positive knot positions not allowed",1);
  for (unsigned int i = 0; i < nknots; ++i)
    knots[i] = S[i];
  for (unsigned int i = 0; i < nknots; ++i)
    logknots[i] = log2(knots[i]);
  if (nknots > 0) knotpos_loaded = true;
}

/*!
  \param[in] n Number of knots in band 1 model
  \param[in] S Input knot positions
*/
void numberCountsDoubleLogNormal::setKnotPositions(unsigned int n, 
						   const double* const S) {
  setNKnots(n); //Also Invalidates everything
  for (unsigned int i = 0; i < nknots; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsDoubleLogNormal","setKnotPositions",
			 "Non-positive knot positions not allowed",1);
  for (unsigned int i = 0; i < nknots; ++i)
    knots[i] = S[i];
  for (unsigned int i = 0; i < nknots; ++i)
    logknots[i] = log2(knots[i]);
  if (nknots > 0) knotpos_loaded = true;
}

/*!
  \param[in] idx Index into knot positions
  \returns Position of idx-th knot

  Access is not checked.
*/
double numberCountsDoubleLogNormal::
getKnotPosition(const unsigned int i) const {
  return knots[i];
}


/*!
  \param[in] S Input knot positions for band 2 color model sigma
*/
void numberCountsDoubleLogNormal::
setSigmaPositions(const std::vector<double>& S) {
  unsigned int n = S.size();
  setNSigmas(n); //Also invalidates
  for (unsigned int i = 0; i < n; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsDoubleLogNormal","setSigmaPositions",
			 "Non-positive sigma knot positions not allowed",1);
  for (unsigned int i = 0; i < n; ++i)
    sigmaknots[i] = S[i];
  if (nsigmaknots > 0) sigmapos_loaded = true;
}

/*!
  \param[in] n Number of knots for band 2 color model sigma
  \param[in] S Input knot positions
*/
void numberCountsDoubleLogNormal::setSigmaPositions(unsigned int n, 
						    const double* const S) {
  setNSigmas(n); //Invalidates
  for (unsigned int i = 0; i < n; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsDoubleLogNormal","setKnots",
			 "Non-positive sigma knot positions not allowed",1);
  for (unsigned int i = 0; i < n; ++i)
    sigmaknots[i] = S[i];
  if (nsigmaknots > 0) sigmapos_loaded = true;
}

/*!
  \param[in] idx Index into sigma knot positions
  \returns Position of idx-th sigma knot

  Access is not checked.
*/
double numberCountsDoubleLogNormal::
getSigmaPosition(const unsigned int i) const {
  return sigmaknots[i];
}

/*!
  \param[in] S Input knot positions for band 2 color model offset
*/
void numberCountsDoubleLogNormal::
setOffsetPositions(const std::vector<double>& S) {
  unsigned int n = S.size();
  setNOffsets(n);
  for (unsigned int i = 0; i < n; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsDoubleLogNormal","setOffsetPositions",
			 "Non-positive offset knot positions not allowed",1);
  for (unsigned int i = 0; i < n; ++i)
    offsetknots[i] = S[i];
  if (noffsetknots > 0) offsetpos_loaded = true;
}

/*!
  \param[in] idx Index into offset knot positions
  \returns Position of idx-th offset knot

  Access is not checked.
*/
double numberCountsDoubleLogNormal::
getOffsetPosition(const unsigned int i) const {
  return offsetknots[i];
}

/*!
  \param[in] n Number of knots for band 2 color model offset
  \param[in] S Input knot positions
*/
void numberCountsDoubleLogNormal::setOffsetPositions(unsigned int n, 
						     const double* const S) {
  setNOffsets(n);
  for (unsigned int i = 0; i < n; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsDoubleLogNormal","setKnots",
			 "Non-positive offset knot positions not allowed",1);
  for (unsigned int i = 0; i < n; ++i)
    offsetknots[i] = S[i];
  if (noffsetknots > 0) offsetpos_loaded = true;
}

/*!
  \param[out] K Band 1 knot positions vector
  \param[out] S Band 2 color model sigma knot positions vector 
  \param[out] O Band 2 color model offset knot positions vector 
*/
void numberCountsDoubleLogNormal::
getPositions(std::vector<double>& K, std::vector<double>& S,
	     std::vector<double>& O) const {
  if (!knotpos_loaded) 
    throw affineExcept("numberCountsDoubleLogNormal", "getPositions",
		       "Knot positions not set", 1);
  if (!sigmapos_loaded) 
    throw affineExcept("numberCountsDoubleLogNormal", "getPositions",
		       "Sigma knot positions not set", 2);
  if (!offsetpos_loaded) 
    throw affineExcept("numberCountsDoubleLogNormal", "getPositions",
		       "Offset knot positions not set", 3);

  K.resize(nknots);
  if (nknots > 0)
    for (unsigned int i = 0; i < nknots; ++i)
      K[i] = knots[i];
  S.resize(nsigmaknots);
  if (nsigmaknots > 0)
    for (unsigned int i = 0; i < nsigmaknots; ++i)
      S[i] = sigmaknots[i];
  O.resize(noffsetknots);
  if (noffsetknots > 0)
    for (unsigned int i = 0; i < noffsetknots; ++i)
      O[i] = offsetknots[i];
}

/*!
  \param[in] K Band 1 knot positions vector
  \param[in] S Band 2 color model sigma knot positions vector 
  \param[in] O Band 2 color model offset knot positions vector 
*/
void numberCountsDoubleLogNormal::
setPositions(const std::vector<double>& K, const std::vector<double>& S,
	     const std::vector<double>& O) {
  setKnotPositions(K);
  setSigmaPositions(S);
  setOffsetPositions(O);
}

/*! 
  \param[in] sz New size of R work arrays

  Allocates R work arrays.  Only upsizes 
*/
void numberCountsDoubleLogNormal::setRWorkSize(unsigned int sz) const {
  if (sz <= nRWork) return;
  if (RWork1 != NULL) delete[] RWork1;
  if (RWork2 != NULL) delete[] RWork2;
  if (RWork3 != NULL) delete[] RWork3;
  if (RWorkValid != NULL) delete[] RWorkValid;
  if (sz > 0) {
    RWork1 = new double[sz];
    RWork2 = new double[sz];
    RWork3 = new double[sz];
    RWorkValid = new bool[sz];
  } else {
    RWork1 = RWork2 = RWork3 = NULL;
    RWorkValid = NULL;
  }
  nRWork = sz;
}

//////////////////////////////////////////////////////////////////////////
//The stuff below here is specific to this particular model, the stuff above
// may be useful for an abstract base class of band1 spline plus color models.
// Currently this isn't implemented because we only really have this one model

/*!
  \param[out] F Parameters to get from model
  
  Will set the first nknot + nsigma + noffset parameters, ignoring any others.
  If the knot values aren't loaded, sets them to NaN.
 */
void numberCountsDoubleLogNormal::getParams(paramSet& F) const {
  unsigned int nneeded = nknots + nsigmaknots + noffsetknots;
  if (F.getNParams() < nneeded)
    throw affineExcept("numberCountsDoubleLogNormal", "getKnots",
		       "Not enough space in output variable", 1);
  if (knotvals_loaded && offsetvals_loaded && sigmavals_loaded) {
    for (unsigned int i = 0; i < nknots; ++i)
      F[i] = pofd_mcmc::ilogfac * logknotvals[i]; //Convert to base 10
    for (unsigned int i = 0; i < nsigmaknots; ++i)
      F[i + nknots] = sigmavals[i];
    for (unsigned int i = 0; i < noffsetknots; ++i)
      F[i + nknots + nsigmaknots] = offsetvals[i];
  } else 
    for (unsigned int i = 0; i < nneeded; ++i)
      F[i] = std::numeric_limits<double>::quiet_NaN();
}

/*!
  \param[in] F Parameters to set in model

  This will ignore any parameters beyond those it needs.
 */
void numberCountsDoubleLogNormal::setParams(const paramSet& F) {
  unsigned int nneeded = nknots + nsigmaknots + noffsetknots;
  if (nneeded > F.getNParams())
    throw affineExcept("numberCountsDoubleLogNormal","setKnots",
		       "Not enough parameters present to set",1);
  if (!(knotpos_loaded && sigmapos_loaded && offsetpos_loaded))
    throw affineExcept("numberCountsDoubleLogNormal","setKnots",
		       "Some positions not set", 2);

  for (unsigned int i = 0; i < nknots; ++i)
    logknotvals[i] = pofd_mcmc::logfac*F[i]; //Convert to base 2
  if (nknots > 2) 
    gsl_spline_init(splinelog, logknots, logknotvals,
		    static_cast<size_t>(nknots));
  else throw affineExcept("numberCountsDoubleLogNormal", "setParams",
			  "Need at least 3 knots", 3);
  knotvals_loaded = true;
  checkKnotsValid();

  for (unsigned int i = nknots; i < nknots+nsigmaknots; ++i)
    sigmavals[i-nknots] = F[i];
  if (nsigmaknots > 1)
    gsl_interp_init(sigmainterp, sigmaknots, sigmavals,
		    static_cast<size_t>(nsigmaknots));
  sigmavals_loaded = true;
  checkSigmasValid();
  
  unsigned int istart = nknots+nsigmaknots;
  for (unsigned int i = istart; i < istart+noffsetknots; ++i)
    offsetvals[i-istart] = F[i];
  if (noffsetknots > 1)
    gsl_interp_init(offsetinterp, offsetknots, offsetvals,
		    static_cast<size_t>(noffsetknots));
  offsetvals_loaded = true;
  checkOffsetsValid();
}


/*!
  This is an expensive check, so we use this by storing
  the results whenever the knots are set.
*/
void numberCountsDoubleLogNormal::checkKnotsValid() const {
  knots_valid = false;
  if (!knotpos_loaded) return;
  if (!knotvals_loaded) return;
  if (nknots < 3) return;
  for (unsigned int i = 0; i < nknots; ++i)
    if (std::isnan(knots[i]) || std::isinf(knots[i])) return;
  if (knots[0] <= 0.0) return;
  for (unsigned int i = 1; i < nknots; ++i)
    if (knots[i] <= knots[i-1]) return;
  for (unsigned int i = 0; i < nknots; ++i)
    if (std::isnan(logknotvals[i])) return;
  knots_valid = true;
}

/*!
  This is an expensive check, so we use this by storing
  the results whenever the sigmas are set.
*/
void numberCountsDoubleLogNormal::checkSigmasValid() const {
  sigmas_valid = false;
  if (!sigmapos_loaded) return;
  if (!sigmavals_loaded) return;
  if (nsigmaknots == 0) return;
  for (unsigned int i = 0; i < nsigmaknots; ++i)
    if (std::isnan(sigmaknots[i]) || std::isinf(sigmaknots[i])) return;
  if (sigmaknots[0] <= 0.0) return;
  for (unsigned int i = 1; i < nsigmaknots; ++i)
    if (sigmaknots[i] <= sigmaknots[i-1]) return;
  for (unsigned int i = 0; i < nsigmaknots; ++i)
    if (std::isnan(sigmavals[i]) || std::isinf(sigmavals[i])) return;
  sigmas_valid = true;
}

/*!
  \returns True if a valid set of offset knots has been loaded
  
  This is an expensive check, so we use this by storing
  the results whenever the offsets are set.
*/
void numberCountsDoubleLogNormal::checkOffsetsValid() const {
  offsets_valid = false;
  if (!offsetpos_loaded) return;
  if (!offsetvals_loaded) return;
  if (noffsetknots == 0) return;
  for (unsigned int i = 0; i < noffsetknots; ++i)
    if (std::isnan(offsetknots[i]) || std::isinf(offsetknots[i])) return;
  if (offsetknots[0] <= 0.0) return;
  for (unsigned int i = 1; i < noffsetknots; ++i)
    if (offsetknots[i] <= offsetknots[i-1]) return;
  for (unsigned int i = 0; i < noffsetknots; ++i)
    if (std::isnan(offsetvals[i])) return;
  offsets_valid = true;
}

/*!
  \returns True if the model parameters are valid
 */
bool numberCountsDoubleLogNormal::isValid() const {
  if (!(knots_valid && sigmas_valid && offsets_valid)) return false;
  return true;
}

/*!
  \param[in] f1 Band 1 flux value to evaluate sigma for
  \returns Value of band 2 color model sigma at f1

  Assumes validity already checked
*/
double numberCountsDoubleLogNormal::getSigmaInner(double f1) const {
  if (nsigmaknots == 1) return sigmavals[0];
  if (f1 <= sigmaknots[0]) return sigmavals[0];
  if (f1 >= sigmaknots[nsigmaknots-1]) return sigmavals[nsigmaknots-1];
  return gsl_interp_eval(sigmainterp, sigmaknots, sigmavals, 
			 f1, accsigma);
}

/*!
  \param[in] f1 Band 1 flux value to evaluate offset for
  \returns Value of band 2 color model offset at f1

  Assumes validity already checked
*/
double numberCountsDoubleLogNormal::getOffsetInner(double f1) const {
  if (noffsetknots == 1) return offsetvals[0];
  if (f1 <= offsetknots[0]) return offsetvals[0];
  if (f1 >= offsetknots[noffsetknots-1]) return offsetvals[noffsetknots-1];
  return gsl_interp_eval(offsetinterp, offsetknots, offsetvals,
			 f1, accoffset);
}

/*!
  \param[in] f1 Band 1 flux value to evaluate sigma for
  \returns Value of band 2 color model sigma at f1

  Does validity checks on model state.
*/
double numberCountsDoubleLogNormal::getSigma(double f1) const {
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  if (nsigmaknots < 1) return std::numeric_limits<double>::quiet_NaN();
  if (std::isnan(f1) || std::isinf(f1)) 
    return std::numeric_limits<double>::quiet_NaN();
  return getSigmaInner(f1);
}

/*!
  \param[in] f1 Band 1 flux value to evaluate offset for
  \returns Value of band 2 color model offset at f1

  Does validity checks on model state.
*/
double numberCountsDoubleLogNormal::getOffset(double f1) const {
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  if (noffsetknots < 1) return std::numeric_limits<double>::quiet_NaN();
  if (std::isnan(f1) || std::isinf(f1)) 
    return std::numeric_limits<double>::quiet_NaN();
  return getOffsetInner(f1);
}

/*!
  \param[in] f1 Flux density in band 1
  \param[in] f2 Flux density in band 2
  \returns Number counts at f1, f2

  Assumes model validity already checked
*/
double numberCountsDoubleLogNormal::
getNumberCountsInner(double f1, double f2) const {
  const double normfac = 1.0 / sqrt(2 * M_PI);
  if (f1 < knots[0] || f1 >= knots[nknots-1] || f2 <= 0.0) 
    return 0.0; //Out of range

  //This is the n_1 bit
  double cnts = exp2(gsl_spline_eval(splinelog, log2(f1), acc)); 

  //Counts in band 2, Log Normal in f2/f1, multiply them onto n_1
  double if1 = 1.0 / f1;
  double isigma = 1.0/getSigmaInner(f1);
  double tfac = (log(f2 * if1) - getOffsetInner(f1)) * isigma;
  cnts *= normfac * isigma * exp(-0.5 * tfac * tfac) / f2; //yes, it's 1/f2 here
  return cnts;
}

/*!
  \param[in] f1 Flux density in band 1
  \param[in] f2 Flux density in band 2
  \returns Number counts at f1, f2

  Does validity checks on input.
*/
double numberCountsDoubleLogNormal::getNumberCounts(double f1, double f2) 
  const {
  if ((nknots < 2) || (nsigmaknots < 1) || (noffsetknots < 1))
    return std::numeric_limits<double>::quiet_NaN();
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  if (std::isnan(f1) || std::isinf(f1)) 
    return std::numeric_limits<double>::quiet_NaN();
  if (std::isnan(f2) || std::isinf(f2)) 
    return std::numeric_limits<double>::quiet_NaN();
  return getNumberCountsInner(f1,f2);
}

/*
  \param[in] band Which band to select; 0 for band 1, 1 for band 2
  \returns The minimum flux supported by the model.

  Fluxes must be positive, so return smallest non-zero double value.
 */
double numberCountsDoubleLogNormal::getMinFlux(unsigned int band) const {
  if (nknots == 0) return std::numeric_limits<double>::quiet_NaN();
  if (band == 0) return std::numeric_limits<double>::min();
  else if (band == 1) return std::numeric_limits<double>::min();
  else throw affineExcept("numberCountsDoubleLogNormal","getMinFlux",
			  "Invalid band (must be 0 or 1)",1);
  return std::numeric_limits<double>::quiet_NaN();
}

/*
  \param[in] band Which band to select; 0 for band 1, 1 for band 2
  \returns The maximum flux supported by the model

  This is not entirely defined for band 2, so return a value a few sigma
  above the top value.  It is well defined for band 1.
 */
double numberCountsDoubleLogNormal::getMaxFlux(unsigned int band) const {
  const double sigmult = 4.0; //Number of sigma above
  if (nknots == 0) return std::numeric_limits<double>::quiet_NaN();
  if (band == 0) return knots[nknots-1];
  else if (band == 1) {
    double kv = knots[nknots-1];
    //Get the expectation and variance in S2/S1 at the top S1 knot
    double sg = getSigmaInner(kv);
    double mu = getOffsetInner(kv);
    double sg2 = sg * sg;
    double mn = exp(mu + 0.5 * sg2);
    double var = exp(2 * mu + sg2)*(exp(sg2) - 1.0);
    return kv * (mn + sigmult * sqrt(var));
  } else throw affineExcept("numberCountsDoubleLogNormal","getMaxFlux",
			    "Invalid band (must be 0 or 1)",1);
}

/*!
  Compute
  \f[
    \int dS_1 \int dS_2 \, S_1^{\alpha} S_2^{\beta} \frac{dN}{dS_1 dS_2} =
      \int dS_1 \, S_1^{\alpha+\beta} n_1\left(S_1\right) \exp \left[
       \beta \mu \left(S_1\right) + \frac{1}{2} \beta^2 \sigma^2
       \left( S_1 \right) \right]
  \f]
  where \f$n_1\left(S_1\right)\f$ is the number counts in band 1.
  This simplification takes advantage of the rather special form
  of the number counts model, which was basically designed to keep
  this simple and avoid actually having to integrate over \f$S_2\f$.  

  Note that to compute the average flux of each object, you have to 
  divide by the total number (i.e., with \f$\alpha = \beta = 0\f$).
  If you want the mean flux per area, however, you don't divide.

  \param[in] alpha   Power of flux in band 1
  \param[in] beta    Power of flux in band 2
  \returns Integral

  The function evaluation is done in evalPowfNDoubleLogNormal, 
  which is the thing that is passed to the GSL integration routines.
 */
double numberCountsDoubleLogNormal::splineInt(double alpha, double beta) const {
  if (nknots < 2) return std::numeric_limits<double>::quiet_NaN();
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  double result, error;
  void *params;
  gsl_function F;
  double minknot = knots[0];
  double maxknot = knots[nknots-1];
  unsigned int noff = noffsetknots;
  unsigned int nsig = nsigmaknots;
  
  //What we actually pass to evalPowfNDoubleLogNormal is
  // power = alpha+beta
  // const1 = beta
  // const2 = 1/2 beta^2
  //So it evaluates 
  // S_1^power1 n_1(S_1) exp( const1*mu(S_1) + const2*sigma^2(S_2))
  double power = alpha + beta;
  double const1 = beta;
  double const2 = 0.5 * beta * beta;

  //There are a -ton- of other things to set though, so that
  // evalfN knows what to do in detail (minima, maxima, etc.)

  //Stuff we always need
  varr[0] = static_cast<void*>(&power);
  varr[1] = static_cast<void*>(&const1);
  varr[2] = static_cast<void*>(&const2);
  varr[3] = static_cast<void*>(splinelog);
  varr[4] = static_cast<void*>(acc);
  varr[5] = static_cast<void*>(&minknot);
  varr[6] = static_cast<void*>(&maxknot);
  if (const1 != 0) {
    //We will be evaluating the mu term, so need to pass these
    varr[7]  = static_cast<void*>(offsetinterp);
    varr[8]  = static_cast<void*>(accoffset);
    varr[9]  = static_cast<void*>(&noff);
    varr[10] = static_cast<void*>(offsetknots);
    varr[11] = static_cast<void*>(offsetvals);
  }
  if (const2 != 0) {
    //We will be evaluating the sigma term, so need to pass these
    varr[12] = static_cast<void*>(sigmainterp);
    varr[13] = static_cast<void*>(accsigma);
    varr[14] = static_cast<void*>(&nsig);
    varr[15] = static_cast<void*>(sigmaknots);
    varr[16] = static_cast<void*>(sigmavals);
  }

  params = static_cast<void*>(varr);

  F.function = &evalPowfNDoubleLogNormal;
  F.params = params;

  gsl_integration_qag(&F, minknot, maxknot, 0, 1e-5, 1000,
		      GSL_INTEG_GAUSS41, gsl_work, &result, &error); 
  return result;
}

/*!
  \returns The Total number of sources per area, however area is defined
           by the model.
*/
double numberCountsDoubleLogNormal::getNS() const {
  return splineInt(0.0,0.0);
}

/*!
  \param[in] band Which band to select; 0 for band 1, 1 for band 2
  \returns The flux per area -- however the model defines area and
           flux
*/
double numberCountsDoubleLogNormal::
getFluxPerArea(unsigned int band) const {
  if (band == 0) return splineInt(1.0, 0.0);
  else if (band == 1) return splineInt(0.0, 1.0);
  else throw affineExcept("numberCountsDoubleLogNormal",
			  "getFluxPerArea","Band must be 0 or 1",1);
}

/*!
  \param[in] band Which band to select; 0 for band 1, 1 for band 2
  \returns The flux density^2 per area -- however the model defines area and
           flux
*/
double numberCountsDoubleLogNormal::
getFluxSqPerArea(unsigned int band) const {
  if (band == 0) return splineInt(2.0, 0.0);
  else if (band == 1) return splineInt(0.0, 2.0);
  else throw affineExcept("numberCountsDoubleLogNormal",
			  "getFluxSqPerArea","Band must be 0 or 1",1);
}

/*!
  \param[in] p1 Power of band 1 flux density
  \param[in] p2 Power of band 2 flux density
  \returns The flux density (band1)^p1 * flux density (band2)^p2 per area
 */
double numberCountsDoubleLogNormal::
getFluxPowerPerArea(double p1, double p2) const {
  return splineInt(p1, p2);
}

/*!
  \param[in] f1   Flux 1
  \param[in] f2   Flux 2
  \param[in] bm   Beam
  \param[in] sgn  Sign index [pp,pn,np,nn]
  \returns R but without the pixel area factor

  Does not check inputs for validity, and does not include area prefactor
*/
double numberCountsDoubleLogNormal::getRInternal(double f1,double f2,
						 const doublebeam& bm, 
						 unsigned int sgn) const {
  double ieta1, ieta2, retval;
  const double* iparr1;
  const double* iparr2;
  unsigned int npsf = bm.getNPix(sgn);
  iparr1 = bm.getInvPixArr1(sgn);
  iparr2 = bm.getInvPixArr2(sgn);
  retval = 0.0;
  if (bm.hasWeights(sgn)) {
    const double* warr = bm.getWeights(sgn);
    for (unsigned int i = 0; i < npsf; ++i) {
      ieta1 = iparr1[i];
      ieta2 = iparr2[i];
      retval += warr[i] * ieta1 * ieta2 * 
	getNumberCountsInner(f1 * ieta1, f2 * ieta2);
    } 
  } else {
    for (unsigned int i = 0; i < npsf; ++i) {
      ieta1 = iparr1[i];
      ieta2 = iparr2[i];
      retval += ieta1 * ieta2 * getNumberCountsInner(f1 * ieta1, f2 * ieta2);
    } 
  }
  return retval;
}

/*!
  \param[in] n1   Number of fluxes, band 1
  \param[in] f1   Flux in band 1, length n1
  \param[in] n2   Number of fluxes, band 2
  \param[in] f2   Flux 2, length n2
  \param[in] bm   Beam
  \param[in] sgn  Sign index [pp,pn,np,nn] 0-3
  \param[out] R   R value is added to this, dimension n1*n2.  
                  Does not include pixel area factor

  It is up to the caller to initialize R (to zeros, for example)
  Does not check inputs for validity, and does not include area prefactor.
*/
void numberCountsDoubleLogNormal::
getRInternal(unsigned int n1, const double* const f1,
	     unsigned int n2, const double* const f2,
	     const doublebeam& bm, unsigned int sgn,
	     double* R) const {

  double minknot = knots[0];
  double maxknot = knots[nknots-1];

  double maxf1 = maxknot / bm.getMax1(sgn); //!< Max non-zero flux, band 1

  unsigned int npsf = bm.getNPix(sgn);
  setRWorkSize(npsf);

  const double normfac = 1.0/sqrt(2*M_PI);
  double *rowptr;
  const double *iparr1;
  const double *iparr2;
  double ieta1, workval, f1val, f1prod, f2val, f2prod, isigma, tfac, if2;

  iparr1 = bm.getInvPixArr1(sgn);
  iparr2 = bm.getInvPixArr2(sgn);
  if (bm.hasWeights(sgn)) {
    const double* warr = bm.getWeights(sgn);
    for (unsigned int i = 0; i < n1; ++i) {
      rowptr = R + i * n2;
      f1val = f1[i];
      if (f1val <= 0 || f1val >= maxf1) { 
	//R always zero at these boundaries or outside them
	//Do nothing
      } else {
	//Precompute f1 related values
	for (unsigned int j = 0; j < npsf; ++j) {
	  ieta1 = iparr1[j];
	  f1prod = f1val * ieta1;
	  if (f1prod >= minknot && f1prod < maxknot) {
	    RWork1[j] = log(f1prod) + getOffsetInner(f1prod);
	    isigma = 1.0 / getSigmaInner(f1prod);
	    RWork3[j] = -0.5 * isigma * isigma;
	    RWork2[j] = warr[j] * normfac * ieta1 * isigma *
	      exp2(gsl_spline_eval(splinelog, log2(f1prod), acc));
	    RWorkValid[j] = true;
	  } else RWorkValid[j] = false;
	}
	//Now loop over flux 2 values
	for (unsigned int j = 0; j < n2; ++j) {
	  f2val = f2[j];
	  workval = 0;
	  if (f2val > 0) {
	    if2 = 1.0 / f2val;
	    for (unsigned int k = 0; k < npsf; ++k)
	      if (RWorkValid[k]) {
		f2prod = f2val * iparr2[k];
		tfac = log(f2prod) - RWork1[k];
		workval += RWork2[k] * if2 * exp(tfac * tfac * RWork3[k]);
	      }
	  }
	  rowptr[j] += workval;
	}
      }
    }
  } else {
    for (unsigned int i = 0; i < n1; ++i) {
      rowptr = R + i * n2;
      f1val = f1[i];
      if (f1val <= 0.0 || f1val >= maxf1) { 
	//Do nothing
      } else {
	for (unsigned int j = 0; j < npsf; ++j) {
	  ieta1 = iparr1[j];
	  f1prod = f1val * ieta1;
	  if (f1prod >= minknot && f1prod < maxknot) {
	    RWork1[j] = log(f1prod) + getOffsetInner(f1prod);
	    isigma = 1.0 / getSigmaInner(f1prod);
	    RWork3[j] = -0.5 * isigma * isigma;
	    RWork2[j] = normfac * ieta1 * isigma *
	      exp2(gsl_spline_eval(splinelog, log2(f1prod), acc));
	    RWorkValid[j] = true;
	  } else RWorkValid[j] = false;
	}
	//Now loop over flux 2 values
	for (unsigned int j = 0; j < n2; ++j) {
	  f2val = f2[j];
	  workval = 0;
	  if (f2val > 0) {
	    if2 = 1.0/f2val;
	    for (unsigned int k = 0; k < npsf; ++k)
	      if (RWorkValid[k]) {
		f2prod = f2val * iparr2[k];
		tfac = log(f2prod) - RWork1[k];
		workval += RWork2[k] * if2 * exp(tfac * tfac * RWork3[k]);
	      }
	  }
	  rowptr[j] += workval;
	}
      }
    }
  }
}


// Pure brute force
/*
  \param[in] fluxdensity  Position to get R at
  \param[in] params Parameters
  \param[in] beam Beam
  \param[in] bmtype What type of R value to get (pos-pos, pos-neg, neg-pos,
                   neg-neg or the sum of all types that are present)
  \returns R value, including the beam pixel size area factor  
*/
double numberCountsDoubleLogNormal::getR(double f1, double f2, 
					 const doublebeam& bm,
					 const rtype bmtype) const {

  if (!isValid())
    return std::numeric_limits<double>::quiet_NaN();

  if ((f1 <= 0.0) || (f2 <= 0.0)) return 0.0;
  if (f1 > knots[nknots-1]) {
    //Since max(beam) = 1
    return 0.0;
  }

  double Rval;

  switch (bmtype) {
  case BEAMPOS :
    if (! bm.hasSign(0))
      return std::numeric_limits<double>::quiet_NaN();
    Rval = getRInternal(f1,f2,bm,0);
    break;
  case BEAMPOSNEG :
    if (! bm.hasSign(1))
      return std::numeric_limits<double>::quiet_NaN();
    Rval = getRInternal(f1,f2,bm,1);
    break;
  case BEAMNEGPOS :
    if (! bm.hasSign(2))
      return std::numeric_limits<double>::quiet_NaN();
    Rval = getRInternal(f1,f2,bm,2);
    break;
  case BEAMNEG :
    if (! bm.hasSign(3))
      return std::numeric_limits<double>::quiet_NaN();
    Rval = getRInternal(f1,f2,bm,3);
    break;
  case BEAMALL :
    Rval = 0.0;
    for (unsigned int i = 0; i < 4; ++i) {
      if (bm.hasSign(i))
	Rval += getRInternal(f1,f2,bm,i);
    }
    break;
  default :
    throw affineExcept("numberCountsDoubleLogNormal","getR",
		       "Invalid bmtype",1);
    break;
  }

  double prefac;
  prefac = bm.getPixSize()/3600.0;  //To sq deg
  
  return prefac * prefac * Rval;

}

/*!
  Array version.

  \param[in]  n1 Number of fluxes along dimension 1
  \param[in]  f1 Fluxes along dimension 1
  \param[in]  n2 Number of fluxes along dimension 2
  \param[in]  f2 Fluxes along dimension 2
  \param[in]  bm Holds beam information
  \param[out] vals Row major array (n1 x n2) giving values of R, including
                    the beam pixel size prefactor
  \param[in]  bmtype What type of R value to get (pos-pos, pos-neg, neg-pos,
                   neg-neg or the sum of all types)

  This does zero out vals, but it is up to the caller to allocate memory
  for it.
 */ 
void numberCountsDoubleLogNormal::getR(unsigned int n1, const double* const f1,
				       unsigned int n2, const double* const f2, 
				       const doublebeam& bm, double* vals,
				       const rtype bmtype) 
  const {

  if ((!isValid()) || (n1 == 0) || (n2 == 0)) {
    for (unsigned int i = 0; i < n1*n2; ++i)
	vals[i] = std::numeric_limits<double>::quiet_NaN();
    return;
  }

  //Zero out R
  for (unsigned int i = 0; i < n1*n2; ++i)
    vals[i] = 0.0;

  //Compute R
  switch (bmtype) {
  case BEAMPOS :
    if (! bm.hasSign(0)) {
      for (unsigned int i = 0; i < n1*n2; ++i)
	vals[i] = std::numeric_limits<double>::quiet_NaN();
      return;
    }
    getRInternal(n1,f1,n2,f2,bm,0,vals);
    break;
  case BEAMPOSNEG :
    if (! bm.hasSign(1)) {
      for (unsigned int i = 0; i < n1*n2; ++i)
	vals[i] = std::numeric_limits<double>::quiet_NaN();
      return;
    }
    getRInternal(n1,f1,n2,f2,bm,1,vals);
    break;
  case BEAMNEGPOS :
   if (! bm.hasSign(2)) {
      for (unsigned int i = 0; i < n1*n2; ++i)
	vals[i] = std::numeric_limits<double>::quiet_NaN();
      return;
    }
    getRInternal(n1,f1,n2,f2,bm,2,vals);
    break;
  case BEAMNEG :
   if (! bm.hasSign(2)) {
      for (unsigned int i = 0; i < n1*n2; ++i)
	vals[i] = std::numeric_limits<double>::quiet_NaN();
      return;
    }
    getRInternal(n1,f1,n2,f2,bm,3,vals);
    break;
  case BEAMALL :
    for (unsigned int k = 0; k < 4; ++k) {
      if (bm.hasSign(k)) getRInternal(n1,f1,n2,f2,bm,k,vals);
    }
    break;
  default :
    throw affineExcept("numberCountsDoubleLogNormal","getR",
		       "Invalid bmtype",1);
    break;
  }

  //Put in prefac
  double prefac;
  prefac = bm.getPixSize() / 3600.0;  //To sq deg
  prefac = prefac * prefac;
  for (unsigned int i = 0; i < n1*n2; ++i)
    vals[i] *= prefac;
}


void numberCountsDoubleLogNormal::sendSelf(MPI_Comm comm, int dest) const {

  //Knots
  MPI_Send(const_cast<bool*>(&knotpos_loaded), 1, MPI::BOOL, dest, 
	   pofd_mcmc::NCDCSENDKPLOADED, comm);
  MPI_Send(const_cast<unsigned int*>(&nknots), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::NCDCSENDNKNOTS, comm);
  if (knotpos_loaded && nknots > 0) {
    MPI_Send(knots, nknots, MPI_DOUBLE, dest, pofd_mcmc::NCDCSENDKNOTS, comm);
    MPI_Send(logknots, nknots, MPI_DOUBLE, dest,
	     pofd_mcmc::NCDCSENDLOGKNOTS, comm);
  }
  MPI_Send(const_cast<bool*>(&knotvals_loaded), 1, MPI::BOOL, dest,
	   pofd_mcmc::NCDCSENDKVLOADED, comm);
  if (knotpos_loaded && knotvals_loaded && nknots > 0) 
    MPI_Send(logknotvals, nknots, MPI_DOUBLE, dest,
	     pofd_mcmc::NCDCSENDLOGKNOTVALS, comm);

  //Sigmas
  MPI_Send(const_cast<bool*>(&sigmapos_loaded), 1, MPI::BOOL, dest, 
	   pofd_mcmc::NCDCSENDSPLOADED, comm);
  MPI_Send(const_cast<unsigned int*>(&nsigmaknots), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::NCDCSENDNSIGMAKNOTS, comm);
  if (sigmapos_loaded && nsigmaknots > 0)
    MPI_Send(sigmaknots, nsigmaknots, MPI_DOUBLE, dest,
	     pofd_mcmc::NCDCSENDSIGMAKNOTS, comm);
  MPI_Send(const_cast<bool*>(&sigmavals_loaded), 1, MPI::BOOL, dest,
	   pofd_mcmc::NCDCSENDSVLOADED, comm);
  if (sigmapos_loaded && sigmavals_loaded && nsigmaknots > 0)
    MPI_Send(sigmavals, nsigmaknots, MPI_DOUBLE, dest,
	     pofd_mcmc::NCDCSENDSIGMAVALS, comm);

  //Offsets
  MPI_Send(const_cast<bool*>(&offsetpos_loaded), 1, MPI::BOOL, dest, 
	   pofd_mcmc::NCDCSENDOPLOADED, comm);
  MPI_Send(const_cast<unsigned int*>(&noffsetknots), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::NCDCSENDNOFFSETKNOTS, comm);
  if (offsetpos_loaded && noffsetknots > 0) 
    MPI_Send(offsetknots, noffsetknots, MPI_DOUBLE, dest,
	     pofd_mcmc::NCDCSENDOFFSETKNOTS, comm);
  MPI_Send(const_cast<bool*>(&offsetvals_loaded), 1, MPI::BOOL, dest,
	   pofd_mcmc::NCDCSENDOVLOADED, comm);
  if (sigmapos_loaded && offsetvals_loaded && noffsetknots > 0) 
    MPI_Send(offsetvals, noffsetknots, MPI_DOUBLE, dest,
	     pofd_mcmc::NCDCSENDOFFSETVALS, comm);

}

void numberCountsDoubleLogNormal::recieveCopy(MPI_Comm comm, int src) {
  unsigned int n;
  bool loaded;

  //Knots
  MPI_Status Info;
  MPI_Recv(&loaded, 1, MPI::BOOL, src, pofd_mcmc::NCDCSENDKPLOADED,
	   comm, &Info);
  MPI_Recv(&n, 1, MPI_UNSIGNED, src, pofd_mcmc::NCDCSENDNKNOTS, comm, &Info);
  setNKnots(n); //sets loaded values to false
  if (loaded && n > 0) {
    MPI_Recv(knots, n, MPI_DOUBLE, src, pofd_mcmc::NCDCSENDKNOTS,
	     comm, &Info);
    MPI_Recv(logknots, n, MPI_DOUBLE, src, pofd_mcmc::NCDCSENDLOGKNOTS,
	     comm, &Info);
    knotpos_loaded = loaded;
  }
  MPI_Recv(&loaded, 1, MPI::BOOL, src, pofd_mcmc::NCDCSENDKVLOADED,
	   comm, &Info);
  if (knotpos_loaded && loaded && nknots > 0) {
    MPI_Recv(logknotvals,nknots,MPI_DOUBLE,src,
	     pofd_mcmc::NCDCSENDLOGKNOTVALS, comm, &Info);
    if (nknots > 2)
      gsl_spline_init(splinelog, logknots, logknotvals,
		      static_cast<size_t>(nknots));
    knotvals_loaded = loaded;
  }
  checkKnotsValid();

  //Sigmas
  MPI_Recv(&loaded, 1, MPI::BOOL, src, pofd_mcmc::NCDCSENDSPLOADED,
	   comm, &Info);
  MPI_Recv(&n, 1, MPI_UNSIGNED, src, pofd_mcmc::NCDCSENDNSIGMAKNOTS,
	   comm, &Info);
  setNSigmas(n);
  if (loaded && nsigmaknots > 0) {
    MPI_Recv(sigmaknots, n, MPI_DOUBLE, src, pofd_mcmc::NCDCSENDSIGMAKNOTS,
	     comm, &Info);
    sigmapos_loaded = loaded;
  }
  MPI_Recv(&loaded, 1, MPI::BOOL, src, pofd_mcmc::NCDCSENDSVLOADED,
	   comm, &Info);
  if (sigmapos_loaded && loaded && nsigmaknots > 0) {
    MPI_Recv(sigmavals, n, MPI_DOUBLE, src, pofd_mcmc::NCDCSENDSIGMAVALS,
	     comm, &Info);
    if (n > 1)
      gsl_interp_init(sigmainterp, sigmaknots, sigmavals,
		      static_cast<size_t>(nsigmaknots));
    sigmavals_loaded = loaded;
  }
  checkSigmasValid();

  //Offsets
  MPI_Recv(&loaded, 1, MPI::BOOL, src, pofd_mcmc::NCDCSENDOPLOADED,
	   comm, &Info);
  MPI_Recv(&n, 1, MPI_UNSIGNED, src, pofd_mcmc::NCDCSENDNOFFSETKNOTS,
	   comm, &Info);
  setNOffsets(n);
  if (loaded && noffsetknots > 0) {
    MPI_Recv(offsetknots, n, MPI_DOUBLE, src, pofd_mcmc::NCDCSENDOFFSETKNOTS,
	     comm, &Info);
    offsetpos_loaded = loaded;
  }
  MPI_Recv(&loaded, 1, MPI::BOOL, src, pofd_mcmc::NCDCSENDOVLOADED,
	   comm, &Info);
  if (offsetpos_loaded & loaded && noffsetknots > 0) {
    MPI_Recv(offsetvals, n, MPI_DOUBLE, src, pofd_mcmc::NCDCSENDSIGMAVALS,
	     comm, &Info);
    if (n > 1)
      gsl_interp_init(sigmainterp, sigmaknots, sigmavals,
		      static_cast<size_t>(nsigmaknots));
    offsetvals_loaded = loaded;
  }
  checkOffsetsValid();

}

bool numberCountsDoubleLogNormal::writeToStream(std::ostream& os) const {
  os << "Model parameters: " << std::endl;
  if (knotvals_loaded) {
    os << " " << std::left << std::setw(13) << "#Flux knot" << "  "
       << std::setw(13) << "Knot value" << std::endl;
    //Convert to log10 for output
    for (unsigned int i = 0; i < nknots; ++i)
      os << " " << std::left << std::setw(13) << knots[i] << "  "
	 << std::setw(13) << pofd_mcmc::ilogfac * logknotvals[i] << std::endl; 

    os << " " << std::left << std::setw(13) << "#Sigma knot" << "  "
       << std::setw(13) << "Sigma value" << std::endl;
    for (unsigned int i = 0; i < nsigmaknots; ++i)
      os << " " << std::left << std::setw(13) << sigmaknots[i] << "  "
	 << std::setw(13) << sigmavals[i] << std::endl; 

    os << " " << std::left << std::setw(13) << "#Offset knot" << "  "
       << std::setw(13) << "Offset value" << std::endl;
    for (unsigned int i = 0; i < noffsetknots; ++i)
      os << " " << std::left << std::setw(13) << offsetknots[i] << "  "
	 << std::setw(13) << offsetvals[i] << std::endl; 
  } else {
    os << " Number of knots: " << nknots << " sigmas: " << nsigmaknots
       << " offsets: " << noffsetknots << std::endl;
  }
  return true;
}


/*!
  Internal function for use in model integrations.
 */
static double evalPowfNDoubleLogNormal(double s1, void* params) {
  //Model is: ( f^power * exp( const1*mu + const2*sigma^2 ) ) * n1(f)
  // where f is the flux in the first band and n1 is the number
  // counts in band 1 -- see splineInt.
  //Params are:
  // params[0]  power
  // params[1]  const1
  // params[2]  const2
  // params[3]  knot interpolant (log)
  // params[4]  knot accelerator
  // params[5]  minknot
  // params[6]  maxknot
  // params[7]  offset interpolant (log)
  // params[8]  offset accelerator
  // params[9]  noffsets
  // params[10] offset positions
  // params[11] offset values
  // params[12] sigma interpolant (log)
  // params[13] sigma accelerator
  // params[14] nsigmas
  // params[15] sigma positions
  // params[16] sigma values
  //But this really has to be an array of pointers to void to work
  void** vptr = static_cast<void**>(params);

  //First get min/max knot in band 1 for quick return if we are outside that
  double minknot = *static_cast<double*>(vptr[5]);
  double maxknot = *static_cast<double*>(vptr[6]);
  if (s1 < minknot || s1 >= maxknot) return 0.0;

  //Get coeffs
  double power  = *static_cast<double*>(vptr[0]);
  double const1 = *static_cast<double*>(vptr[1]);
  double const2 = *static_cast<double*>(vptr[2]);

  //Construct thing we multiply spline by  
  double prefac;
  //Construct s1^power part
  if (fabs(power) < 1e-6)
    prefac = 1.0;
  else if (fabs(power - 1.0) < 1e-6) 
    prefac = s1; 
  else if (fabs(power - 2.0) < 1e-6) 
    prefac = s1 * s1;
  else prefac = pow(s1, power);

  //Now exponential part
  if (const1 != 0 || const2 != 0) {
    double expbit; //Part to exponentiate
    if (const1 != 0) {
      //Evaluate offset at s1
      unsigned int noffsets = *static_cast<unsigned int*>(vptr[9]);
      double *offsetpos = static_cast<double*>(vptr[10]);
      double *offsetval = static_cast<double*>(vptr[11]);
      double mu;
      if (noffsets == 1) mu = offsetval[0];
      else if (s1 <= offsetpos[0]) mu = offsetval[0];
      else if (s1 >= offsetpos[noffsets-1]) mu = offsetval[noffsets-1];
      else {
	gsl_interp* ospl = static_cast<gsl_interp*>(vptr[7]);
	gsl_interp_accel* oacc = static_cast<gsl_interp_accel*>(vptr[8]);
	mu = gsl_interp_eval(ospl, offsetpos, offsetval, s1, oacc);
      }
      expbit = const1 * mu;
    } else expbit = 0.0;

    if (const2 != 0) {
      //Get sigma bit -> const2* sigma^2
      double sigma;
      unsigned int nsigmas = *static_cast<unsigned int*>(vptr[14]);
      double *sigmapos = static_cast<double*>(vptr[15]);
      double *sigmaval = static_cast<double*>(vptr[16]);
      if (nsigmas == 1) sigma = sigmaval[0];
      else if (s1 <= sigmapos[0]) sigma = sigmaval[0];
      else if (s1 >= sigmapos[nsigmas-1]) sigma = sigmaval[nsigmas-1];
      else {
	gsl_interp* sspl = static_cast<gsl_interp*>(vptr[12]);
	gsl_interp_accel* sacc = static_cast<gsl_interp_accel*>(vptr[13]);
	sigma = gsl_interp_eval(sspl, sigmapos, sigmaval, s1, sacc);
      }
      expbit += const2 * sigma * sigma;
    } 

    prefac *= exp(expbit);

  } //Otherwise exp bit is just 1

  //Now multiply in n(band1)
  gsl_spline* spl = static_cast<gsl_spline*>(vptr[3]);
  gsl_interp_accel* acc = static_cast<gsl_interp_accel*>(vptr[4]);
  double splval = exp2(gsl_spline_eval(spl,log2(s1),acc));
  return prefac * splval;
}

/////////////////////////////////////

initFileDoubleLogNormal::initFileDoubleLogNormal() : 
  nknots(0), nsigmas(0), noffsets(0), sigmaidx(0), offsetidx(0),
  has_sigma(false), has_lower_limits(false), has_upper_limits(false) {

  knotpos = knotval = sigma = lowlim = uplim = NULL;
  has_lowlim = has_uplim = NULL;

  //Set random number generator seed from time
  unsigned long long int seed;
  seed = static_cast<unsigned long long int>(time(NULL));
  seed += static_cast<unsigned long long int>(clock());
  rangen.setSeed(seed);
}

/*!
  \param[in] flname File to read from
  \param[in] read_sigma      Read in (and require) knot sigmas
  \param[in] read_limits     Try to read limits; this will turn on require_sigma

  See initFileDoubleLogNormal::readFile for details of the file format
*/
initFileDoubleLogNormal::initFileDoubleLogNormal(const std::string& flname, 
						 bool read_sigma, 
						 bool read_limits) : 
  nknots(0), nsigmas(0), noffsets(0), sigmaidx(0), offsetidx(0),
  has_sigma(false), has_lower_limits(false), has_upper_limits(false) {

  knotpos = knotval = sigma = lowlim = uplim = NULL;
  has_lowlim = has_uplim = NULL;

  //Set random number generator seed from time
  unsigned long long int seed;
  seed = static_cast<unsigned long long int>(time(NULL));
  seed += static_cast<unsigned long long int>(clock());
  rangen.setSeed(seed);

  readFile(flname, read_sigma, read_limits);

}

initFileDoubleLogNormal::~initFileDoubleLogNormal() {
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

  The file format is: 
  first, a line giving the n1 ns no
  where n1 is the number of knots in the band 1 model, ns is the number
  of knots in the color sigma model, and no is the number of knots in the
  color offset model.  This is followed by n1+ns+no lines of the form
  knotpos   knotval   [sigma [ lowlim [ uplim ]]
  So knotpos, knotval are always required
  sigma is optionally required if read_sigma is set
  lowlim and uplim may be present, and are looked for if read_limits is set.
    sigma must also be present, and the first element found is lowlim.
    If another is also found, it is interpreted as uplim
 */
void initFileDoubleLogNormal::readFile(const std::string& flname, 
				       bool read_sigma, bool read_limits) {
  if (read_limits) read_sigma = true;

  //Clear old data
  nknots = 0; 
  nsigmas = 0;
  noffsets = 0;
  sigmaidx = 0;
  offsetidx = 0;
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

  unsigned int nk, ns, no; //Number of knots, sigmas, offsets
  std::vector<double> wvec1, wvec2, wvec3, wvec4, wvec5;
  std::vector<bool> h4, h5;
  std::string line;
  std::vector<std::string> words;
  std::stringstream str;
  double currval;

  std::ifstream initfs(flname.c_str());
  if (!initfs) {
    initfs.close();
    std::stringstream errmsg;
    errmsg << "Unable to open file:" << flname << std::endl;
    throw affineExcept("initFileDoubleLogNormal","readFile",errmsg.str(),1);
  }

  //Read in number
  initfs >> nk >> ns >> no;
  if (nk < 2) {
    initfs.close();
    throw affineExcept("initFileDoubleLogNormal","readFile",
		       "Need at least 2 band 1 knots",2);
  }
  if (ns < 1) {
    initfs.close();
    throw affineExcept("initFileDoubleLogNormal","readFile",
		       "Need at least one sigma color model knot",3);

  }
  if (no < 1) {
    initfs.close();
    throw affineExcept("initFileDoubleLogNormal","readFile",
		       "Need at least one offset color model knot",4);
  }
  nknots = sigmaidx = nk; //These could be combined, but meh
  nsigmas = ns;
  noffsets = no;
  offsetidx = nk+ns;
  
  //Read in values
  while (!initfs.eof()) {
    std::getline(initfs, line);
    if (line[0] == '#') continue; //Comment
    utility::stringwords(line, words);
    if (words.size() == 0) continue; //Nothing on line (with spaces removed)
    if (words[0][0] == '#') continue; //Comment line
    if (words.size() < nreq) continue; //Has wrong number of entries
    str.str(words[0]); str.clear(); str >> currval;
    wvec1.push_back(currval);
    str.str(words[1]); str.clear(); str >> currval;
    wvec2.push_back(currval); 

    //Now, all the options
    if (read_sigma) {
      has_sigma = true;
      str.str(words[2]); str.clear(); str >> currval;
      wvec3.push_back(currval);
    }
    if (read_limits) { 
      if (words.size() > 3) {
	has_lower_limits = true;
	str.str(words[3]); str.clear(); str >> currval;
	h4.push_back(true);
	wvec4.push_back(currval);
	if (words.size() > 4) {
	  has_upper_limits = true;
	  str.str(words[4]); str.clear(); str >> currval;
	  h5.push_back(true);
	  wvec5.push_back(currval);
	} else {
	  h5.push_back(false);
	  wvec5.push_back(std::numeric_limits<double>::quiet_NaN());
	}
      } else {
	h4.push_back(false);
	wvec4.push_back(std::numeric_limits<double>::quiet_NaN());
	h5.push_back(false);
	wvec5.push_back(std::numeric_limits<double>::quiet_NaN());
      }
    }
  }
  initfs.close();

  unsigned int ntot = nk+ns+no;
  if (wvec1.size() != ntot) {
    std::stringstream errstr;
    errstr << "Expected " << ntot << " values, got: " 
	   << wvec1.size();
    throw affineExcept("initFileDoubleLogNormal","readFile",
		       errstr.str(),5);
  }

  //Copy into vectors
  knotpos = new double[ntot];
  for (unsigned int i = 0; i < ntot; ++i) knotpos[i] = wvec1[i];
  knotval = new double[ntot];
  for (unsigned int i = 0; i < ntot; ++i) knotval[i] = wvec2[i];
  if (has_sigma) {
    sigma = new double[ntot];
    for (unsigned int i = 0; i < ntot; ++i) sigma[i] = wvec3[i];
  }
  if (has_lower_limits) {
    has_lowlim = new bool[ntot];
    lowlim = new double[ntot];
    for (unsigned int i = 0; i < ntot; ++i) has_lowlim[i] = h4[i];
    for (unsigned int i = 0; i < ntot; ++i) lowlim[i] = wvec4[i];
  }
  if (has_upper_limits) {
    has_uplim = new bool[ntot];
    uplim = new double[ntot];
    for (unsigned int i = 0; i < ntot; ++i) has_uplim[i] = h5[i];
    for (unsigned int i = 0; i < ntot; ++i) uplim[i] = wvec5[i];
  }

  //Make sure lower/upper limits don't cross
  if (has_lower_limits && has_upper_limits)
    for (unsigned int i = 0; i < ntot; ++i) {
      if (!(has_uplim[i] && has_lowlim[i])) continue;
      if (uplim[i] < lowlim[i]) {
	std::stringstream errstr;
	errstr << "Lower/Upper limits cross at index: " << i << std::endl;
	errstr << " Lower limit: " << lowlim[i] 
	       << " Upper limit: " << uplim[i];
	throw affineExcept("initFileDoubleLogNormal", "readFiles", 
			   errstr.str(), 6);
      }
      if ((sigma[i] > 0.) && (uplim[i] == lowlim[i])) {
	std::stringstream errstr;
	errstr << "Lower/Upper limits meet at index: " << i 
	       << " but sigma is not zero" << std::endl;
	errstr << " Lower limit: " << lowlim[i] << " Upper limit: " << uplim[i]
	       << " sigma: " << sigma[i];
	throw affineExcept("initFileDoubleLogNormal", "readFiles", 
			   errstr.str(), 7);
      }
    }

  //Make sure that if sigmas is 0 then the mean value falls within
  // the range of any limits
  if (has_sigma && (has_lower_limits || has_upper_limits)) {
    for (unsigned int i = 0; i < ntot; ++i)
      if (sigma[i] > 0) {
	if (has_lower_limits && has_lowlim[i] && (knotval[i] < lowlim[i])) {
	  std::stringstream errstr;
	  errstr << "At knot " << i << " sigma is zero but mean value "
		 << knotval[i] << std::endl << " lies below lower limit "
		 << lowlim[i];
	  throw affineExcept("initFileDoubleLogNormal", "readFiles", 
			     errstr.str(), 8);
	}
	if (has_upper_limits && has_uplim[i] && (knotval[i] > uplim[i])) {
	  std::stringstream errstr;
	  errstr << "At knot " << i << " sigma is zero but mean value "
		 << knotval[i] << std::endl << " lies above upper limit "
		 << uplim[i];
	  throw affineExcept("initFileDoubleLogNormal", "readFiles", 
			     errstr.str(), 9);
	}
      }
  }
}

std::pair<double,double> 
initFileDoubleLogNormal::getKnot(unsigned int idx) const {
  if (nknots == 0)
    throw affineExcept("initFileDoubleLogNormal","getKnot",
		       "No knot information read in",1);
  if (idx >= nknots)
    throw affineExcept("initFileDoubleLogNormal","getKnot",
		       "Invalid knot index",2);
  return std::make_pair(knotpos[idx],knotval[idx]);
}

std::pair<double,double> 
initFileDoubleLogNormal::getSigma(unsigned int idx) const {
  if (nsigmas == 0)
    throw affineExcept("initFileDoubleLogNormal","getSigma",
		       "No sigma information read in",1);
  if (idx >= nsigmas)
    throw affineExcept("initFileDoubleLogNormal","getSigma",
		       "Invalid sigma index",2);
  return std::make_pair(knotpos[idx+sigmaidx],knotval[idx+sigmaidx]);
}

std::pair<double,double> 
initFileDoubleLogNormal::getOffset(unsigned int idx) const {
  if (noffsets == 0)
    throw affineExcept("initFileDoubleLogNormal","getOffset",
		       "No offset information read in",1);
  if (idx >= noffsets)
    throw affineExcept("initFileDoubleLogNormal","getOffset",
		       "Invalid offset index",2);
  return std::make_pair(knotpos[idx+offsetidx],knotval[idx+offsetidx]);
}


/*
  \param[inout] model Modified on output; knot positions are set for
                      band 1 model and color model

  This will change the number of knots in the model if they
  don't match.
 */
void 
initFileDoubleLogNormal::getModelPositions(numberCountsDoubleLogNormal& model) 
  const {
  if (nknots + nsigmas + noffsets == 0)
    throw affineExcept("initFileDoubleLogNormal","getModelPositions",
		       "No knot information read in",1);
  model.setKnotPositions(nknots, knotpos);
  model.setSigmaPositions(nsigmas, knotpos+sigmaidx);
  model.setOffsetPositions(noffsets, knotpos+offsetidx);
}

/*
  \param[out] kp Set to knot positions on output for band 1 number counts
 */
void initFileDoubleLogNormal::getKnotPos(std::vector<double>& kp) const {
  if (nknots == 0)
    throw affineExcept("initFileDoubleLogNormal","getKnotPos",
		       "No knot information read in",1);
  kp.resize(nknots);
  for (unsigned int i = 0; i < nknots; ++i)
    kp[i] = knotpos[i];

}

/*
  \param[out] kp Set to knot positions on output for band 1 number counts
 */
void initFileDoubleLogNormal::getKnotVals(std::vector<double>& kv) const {
  if (nknots == 0)
    throw affineExcept("initFileDoubleLogNormal","getKnotVals",
		       "No knot information read in",1);
  kv.resize(nknots);
  for (unsigned int i = 0; i < nknots; ++i)
    kv[i] = knotval[i];
}

/*
  \param[out] kp Set to knot positions on output for color model sigma
 */
void initFileDoubleLogNormal::getSigmaPos(std::vector<double>& kp) const {
  if (nsigmas == 0)
    throw affineExcept("initFileDoubleLogNormal","getSigmaPos",
		       "No sigma information read in",1);
  kp.resize(nsigmas);
  for (unsigned int i = 0; i < nsigmas; ++i)
    kp[i] = knotpos[i+sigmaidx];
}

/*
  \param[out] kp Set to knot values on output for color model sigma
 */
void initFileDoubleLogNormal::getSigmaVals(std::vector<double>& kv) const {
  if (nsigmas == 0)
    throw affineExcept("initFileDoubleLogNormal","getSigmaVals",
		       "No sigma information read in",1);
  kv.resize(nsigmas);
  for (unsigned int i = 0; i < nsigmas; ++i)
    kv[i] = knotval[i+sigmaidx];
}

/*
  \param[out] kp Set to knot positions on output for color model offset
 */
void initFileDoubleLogNormal::getOffsetPos(std::vector<double>& kp) const {
  if (noffsets == 0)
    throw affineExcept("initFileDoubleLogNormal","getOffsetPos",
		       "No offset information read in",1);
  kp.resize(noffsets);
  for (unsigned int i = 0; i < noffsets; ++i)
    kp[i] = knotpos[i+offsetidx];
}

/*
  \param[out] kp Set to knot values on output for color model offset
 */
void initFileDoubleLogNormal::getOffsetVals(std::vector<double>& kv) const {
  if (noffsets == 0)
    throw affineExcept("initFileDoubleLogNormal","getOffsetVals",
		       "No sigma information read in",1);
  kv.resize(noffsets);
  for (unsigned int i = 0; i < noffsets; ++i)
    kv[i] = knotval[i+offsetidx];
}

/*
  \param[out] p Filled on output with the mean knot values.

  This only fills the first nknots + nsigmas + noffsets parameters
 */
void initFileDoubleLogNormal::getParams(paramSet& p) const {
  if (nknots == 0)
    throw affineExcept("initFileDoubleLogNormal","getParams",
		       "No information loaded",1);
  if (p.getNParams() < nknots+noffsets+nsigmas)
    throw affineExcept("initFileDoubleLogNormal","getParams",
		       "Not enough space in provided paramSet",2);
  for (unsigned int i = 0; i < nknots+nsigmas+noffsets; ++i)
    p[i] = knotval[i];
}

/*
  \param[out] p Filled on output with the mean knot values.

  This only fills the first nknots+nsigmas+noffsets parameters
 */
void initFileDoubleLogNormal::generateRandomKnotValues(paramSet& p) const {
  const unsigned int maxiters = 1000; //Maximum number of generation attempts
  unsigned int ntot = nknots + noffsets + nsigmas;
  if (ntot == 0)
    throw affineExcept("initFileDoubleLogNormal", "generateRandomKnotValues",
		       "No knot information read in", 1);
    
  //Make sure p is big enough; don't resize, complain
  if (p.getNParams() < ntot)
    throw affineExcept("initFileDoubleLogNormal", "generateRandomKnotValues",
		       "Not enough space in provided paramSet", 2);

  //Deal with simple case -- everything fixed
  //So just return central values
  if (!has_sigma) {
    for (unsigned int i = 0; i < ntot; ++i)
      p[i] = knotval[i];
    return;
  }

  //Now we have at least some sigmas
  //The simple case is if there are no limits.  If there are, we will
  // have to do trials.
  if (!(has_lower_limits || has_upper_limits)) {
    for (unsigned int i = 0; i < ntot; ++i)
      p[i] = rangen.gauss() * sigma[i] + knotval[i];
  } else {
    //Both sigmas and limits
    bool goodval;
    double trialval;
    unsigned int iters;
    for (unsigned int i = 0; i < ntot; ++i) {
      if (sigma[i] > 0) {
	//Some sanity checks
	if (has_lowlim[i] && (lowlim[i] > knotval[i]+sigma[i]*4.0)) {
	  std::stringstream errstr;
	  errstr << "Lower limit is too far above central value; will take too"
		 << " long to " << std::endl << "generate value for param idx: "
		 << i;
	  throw affineExcept("initFileDoubleLogNormal",
			     "generateRandomKnotValues",
			     errstr.str(), 3);
	}
	if (has_uplim[i] && (uplim[i] < knotval[i]-sigma[i]*4.0)) {
	  std::stringstream errstr;
	  errstr << "Upper limit is too far below central value; will take too"
		 << " long to " << std::endl << "generate value for param idx: "
		 << i;
	  throw affineExcept("initFileDoubleLogNormal",
			     "generateRandomKnotValues",
			     errstr.str(), 4);
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
		throw affineExcept("initFileDoubleLogNormal",
				   "generateRandomKnotValues",
				   errstr.str(), 5);
	      }
	      trialval = rangen.gauss() * sigma[i] + knotval[i];
	      if ((trialval >= lowlim[i]) && (trialval <= uplim[i])) 
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
	      throw affineExcept("initFileDoubleLogNormal",
				 "generateRandomKnotValues",
				 errstr.str(), 6);
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
	      throw affineExcept("initFileDoubleLogNormal",
				 "generateRandomKnotValues",
				 errstr.str(), 7);
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

/*!
  \param[in] idx Which parameter index to get
  \returns Knot sigma at idx
  
  Note that idx ranges over all the parameters -- including
  the color model.  The valid range is thus [0,nknots + nsigmas + noffsets).
*/
double initFileDoubleLogNormal::getKnotSigma(unsigned int idx) const {
  if (idx >= nknots+nsigmas+noffsets)
    throw affineExcept("initFileDoubleLogNormal",
		       "getKnotSigma","Invalid knot index",1);
  if (!has_sigma) return std::numeric_limits<double>::quiet_NaN();
  return sigma[idx];
}

/*!
  \param[in] idx Which paraemter index to return value for
  
  \returns True if the knot is fixed, false if not.

  Note that idx ranges over all the parameters -- including
  the color model. The valid range is thus [0,nknots + nsigmas + noffsets).
*/

bool initFileDoubleLogNormal::isKnotFixed(unsigned int idx) const {
  if (idx >= nknots+nsigmas+noffsets)
    throw affineExcept("initFileDoubleLogNormal","isKnotFixed",
		       "Invalid knot index",1);
  if (!has_sigma) return false;
  if (sigma[idx] == 0) return true;
  return false;
}

/*
  \param[in] idx Knot index
  \returns Whether knot has a lower limit
 */
bool initFileDoubleLogNormal::knotHasLowerLimit(unsigned int idx) const {
  if (!has_lower_limits) return false;
  if (idx >= nknots+nsigmas+noffsets)
    throw affineExcept("initFileDoubleLogNormal","knotHasLowerLimit",
		       "Invalid knot index",1);
  return has_lowlim[idx];
}

/*
  \param[in] idx Knot index
  \returns Lower limit on knot, or NaN if none
 */
double initFileDoubleLogNormal::getLowerLimit(unsigned int idx) const {
  if (!has_lower_limits) return std::numeric_limits<double>::quiet_NaN();
  if (idx >= nknots+nsigmas+noffsets)
    throw affineExcept("initFileDoubleLogNormal","getLowerLimit",
		       "Invalid knot index",1);
  if (!has_lowlim[idx]) return std::numeric_limits<double>::quiet_NaN();
  return lowlim[idx];
}

/*
  \param[in] idx Knot index
  \returns Whether knot has a upper limit
 */
bool initFileDoubleLogNormal::knotHasUpperLimit(unsigned int idx) const {
  if (!has_upper_limits) return false;
  if (idx >= nknots+nsigmas+noffsets)
    throw affineExcept("initFileDoubleLogNormal","knotHasUpperLimit",
		       "Invalid knot index",1);
  return has_uplim[idx];
}


/*
  \param[in] idx Knot index
  \returns Upper limit on knot, or NaN if none
 */
double initFileDoubleLogNormal::getUpperLimit(unsigned int idx) const {
  if (!has_upper_limits) return std::numeric_limits<double>::quiet_NaN();
  if (idx >= nknots+nsigmas+noffsets)
    throw affineExcept("initFileDoubleLogNormal","getUpperLimit",
		       "Invalid knot index",1);
  if (!has_uplim[idx]) return std::numeric_limits<double>::quiet_NaN();
  return uplim[idx];
}


bool initFileDoubleLogNormal::isValid(const paramSet& p) const {
  if (! (has_lower_limits || has_upper_limits)) return true;
  unsigned int ntot = nknots + nsigmas + noffsets;
  if (p.getNParams() < ntot)
    throw affineExcept("initFileDoubleLogNormal","isValid",
		       "Not enough params in paramSet to test validity",1);
  double val;
  for (unsigned int i = 0; i < ntot; ++i) {
    val = p[i];
    if (has_lowlim[i] && (val < lowlim[i])) return false;
    if (has_uplim[i] && (val < uplim[i])) return false;
  }
  return true;
}

/*!
  \param[in] comm MPI communication handle
  \param[in] dest Destination for send operation
 */
void initFileDoubleLogNormal::sendSelf(MPI_Comm comm, int dest) const {
  MPI_Send(const_cast<unsigned int*>(&nknots), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::IFDLNSENDNKNOTS, comm);
  MPI_Send(const_cast<unsigned int*>(&nsigmas), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::IFDLNSENDNSIGMAS, comm);
  MPI_Send(const_cast<unsigned int*>(&noffsets), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::IFDLNSENDNOFFSETS, comm);
  unsigned int ntot = nknots+nsigmas+noffsets;
  if (ntot != 0) {
    MPI_Send(knotpos, ntot, MPI_DOUBLE, dest, 
	     pofd_mcmc::IFDLNSENDKNOTPOS, comm);
    MPI_Send(knotval, ntot, MPI_DOUBLE, dest, 
	     pofd_mcmc::IFDLNSENDKNOTVAL, comm);
    MPI_Send(const_cast<bool*>(&has_sigma), 1, MPI::BOOL, dest, 
	     pofd_mcmc::IFDLNHASSIGMA, comm);
    if (has_sigma)
      MPI_Send(sigma, ntot, MPI_DOUBLE, dest, pofd_mcmc::IFDLNSENDSIGMA, comm);
    MPI_Send(const_cast<bool*>(&has_lower_limits), 1, MPI::BOOL, dest,
	     pofd_mcmc::IFDLNHASLOWERLIMITS, comm);
    if (has_lower_limits) {
      MPI_Send(has_lowlim, ntot, MPI::BOOL, dest,
	       pofd_mcmc::IFDLNSENDHASLOWLIM, comm);
      MPI_Send(lowlim, ntot, MPI_DOUBLE, dest, 
	       pofd_mcmc::IFDLNSENDLOWLIM, comm);
    }
    MPI_Send(const_cast<bool*>(&has_upper_limits), 1, MPI::BOOL, dest,
	     pofd_mcmc::IFDLNHASUPPERLIMITS, comm);
    if (has_upper_limits) {
      MPI_Send(has_uplim, ntot, MPI::BOOL, dest,
	       pofd_mcmc::IFDLNSENDHASUPLIM, comm);
      MPI_Send(uplim, ntot, MPI_DOUBLE, dest, pofd_mcmc::IFDLNSENDUPLIM, comm);
    }
  }
}

void initFileDoubleLogNormal::recieveCopy(MPI_Comm comm, int src) {
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

  MPI_Status Info;
  MPI_Recv(&nknots, 1, MPI_UNSIGNED, src, pofd_mcmc::IFDLNSENDNKNOTS,
	   comm, &Info);
  MPI_Recv(&nsigmas, 1, MPI_UNSIGNED, src, pofd_mcmc::IFDLNSENDNSIGMAS,
	   comm, &Info);
  MPI_Recv(&noffsets, 1, MPI_UNSIGNED, src, pofd_mcmc::IFDLNSENDNOFFSETS,
	   comm, &Info);
  unsigned int ntot = nknots + nsigmas + noffsets;
  sigmaidx = nknots;
  offsetidx = nknots + nsigmas;
  if (ntot > 0) {
    knotpos = new double[ntot];
    MPI_Recv(knotpos, ntot, MPI_DOUBLE, src, pofd_mcmc::IFDLNSENDKNOTPOS,
	     comm, &Info);
    knotval = new double[ntot];
    MPI_Recv(knotval, ntot, MPI_DOUBLE, src, pofd_mcmc::IFDLNSENDKNOTVAL,
	     comm, &Info);
    MPI_Recv(&has_sigma, 1, MPI::BOOL, src, pofd_mcmc::IFDLNHASSIGMA,
	     comm, &Info);
    if (has_sigma) {
      sigma = new double[ntot];
      MPI_Recv(sigma, ntot, MPI_DOUBLE, src, pofd_mcmc::IFDLNSENDSIGMA,
	       comm, &Info);
    }
    MPI_Recv(&has_lower_limits, 1, MPI::BOOL, src,
	     pofd_mcmc::IFDLNHASLOWERLIMITS, comm, &Info);
    if (has_lower_limits) {
      has_lowlim = new bool[ntot];
      MPI_Recv(has_lowlim, ntot, MPI::BOOL, src, 
	       pofd_mcmc::IFDLNSENDHASLOWLIM, comm, &Info);
      lowlim = new double[ntot];
      MPI_Recv(lowlim, ntot, MPI_DOUBLE, src, pofd_mcmc::IFDLNSENDLOWLIM,
	       comm, &Info);
    }
    MPI_Recv(&has_upper_limits, 1, MPI::BOOL, src, 
	     pofd_mcmc::IFDLNHASUPPERLIMITS, comm, &Info);
    if (has_upper_limits) {
      has_uplim = new bool[ntot];
      MPI_Recv(has_uplim, ntot, MPI::BOOL, src, pofd_mcmc::IFDLNSENDHASUPLIM,
	       comm, &Info);
      uplim = new double[ntot];
      MPI_Recv(uplim, ntot, MPI_DOUBLE, src, pofd_mcmc::IFDLNSENDUPLIM,
	       comm, &Info);
    }
  }
}
