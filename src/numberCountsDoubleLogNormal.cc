#include<cmath>
#include<iomanip>
#include<sstream>
#include<limits>
#include<cstring>

#include<fftw3.h>

#include<gsl/gsl_errno.h>
#include<gsl/gsl_math.h>

#include "../include/numberCountsDoubleLogNormal.h"
#include "../include/global_settings.h"
#include "../include/hdf5utils.h"
#include "../include/affineExcept.h"

// Determines sign component of fluxes (pp, pn, np, nn)
inline unsigned int signComp(double x1, double x2) noexcept {
  if (x1 >= 0) {
    if (x2 >= 0) return 0;
    else return 1;
  } else {
    if (x2 >= 0) return 2;
    else return 3;
  }
}

//Functions to pass to GSL integrator
/*! \brief Evaluates flux1^power1 * exp(const1*mu + const2*sigma^2) dN/dS1 */
static double evalPowfNDoubleLogNormal(double, void*) noexcept; 

// Function for root solving to find maximum band 2 flux
struct lognorm_params {
  // Some repetition here in params for efficiency
  double mu, sig2, isig, isig2, target;
};
static double lognorm(double, void*);

/*! \brief Evaluates sqrt(2 pi) dN / dS1 dS2 with S1 as an argument and S2 a parameter*/
static double evalCounts(double, void*) noexcept; 

numberCountsDoubleLogNormal::numberCountsDoubleLogNormal() : 
  nknots(0), knots(nullptr), logknots(nullptr), logknotvals(nullptr),
  splinelog(nullptr), nsigmaknots(0), sigmaknots(nullptr), sigmavals(nullptr),
  min_sigma(0.0), sigmainterp(nullptr), noffsetknots(0), offsetknots(nullptr),
  offsetvals(nullptr), offsetinterp(nullptr),
  knots_valid(false), sigmas_valid(false), offsets_valid(false),
  knotpos_loaded(false), sigmapos_loaded(false), offsetpos_loaded(false),
  knotvals_loaded(false), sigmavals_loaded(false), offsetvals_loaded(false),
  fslv(nullptr) {
  
  nRWork = 0;
  RWork = nullptr;

  nf2work = 0;
  f2work_inv = f2work_log = nullptr;
  f2work_sgn = nullptr;

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
  nknots(NKNOTS), nsigmaknots(NSIGMA), min_sigma(0.0), noffsetknots(NOFFSET), 
  knots_valid(false), sigmas_valid(false), offsets_valid(false),
  knotpos_loaded(false), sigmapos_loaded(false), offsetpos_loaded(false),
  knotvals_loaded(false), sigmavals_loaded(false), offsetvals_loaded(false),
  fslv(nullptr) {

  nRWork = 0;
  RWork = nullptr;

  nf2work = 0;
  f2work_inv = f2work_log = nullptr;
  f2work_sgn = nullptr;

  if (NKNOTS > 0) {
    knots = new double[NKNOTS];
    logknots = new double[NKNOTS]; 
    logknotvals = new double[NKNOTS];
    if (NKNOTS <= 2)
      throw affineExcept("numberCountsDoubleLogNormal", "constructor",
                         "nknots must be > 2");
    splinelog = gsl_spline_alloc(gsl_interp_cspline,
                                 static_cast<size_t>(NKNOTS));
    for (unsigned int i = 0; i < NKNOTS; ++i)
      knots[i] = std::numeric_limits<double>::quiet_NaN();
    for (unsigned int i = 0; i < NKNOTS; ++i)
      logknots[i] = std::numeric_limits<double>::quiet_NaN();
    for (unsigned int i = 0; i < NKNOTS; ++i)
      logknotvals[i] = std::numeric_limits<double>::quiet_NaN();
  } else {
    knots = logknots = logknotvals = nullptr;
    splinelog = nullptr;
  }
  acc = gsl_interp_accel_alloc();

  sigmainterp = nullptr;
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
  } else sigmaknots = sigmavals = nullptr;
  accsigma = gsl_interp_accel_alloc();

  offsetinterp = nullptr;
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
  } else offsetknots = offsetvals = nullptr;
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
numberCountsDoubleLogNormal(const std::vector<float>& KNOTS,
                            const std::vector<float>& SIGMAS,
                            const std::vector<float>& OFFSETS) :
  knots_valid(false), sigmas_valid(false), offsets_valid(false),
  knotpos_loaded(false), sigmapos_loaded(false), offsetpos_loaded(false),
  knotvals_loaded(false), sigmavals_loaded(false), offsetvals_loaded(false),
  fslv(nullptr) {

  nRWork = 0;
  RWork = nullptr;

  nf2work = 0;
  f2work_inv = f2work_log = nullptr;
  f2work_sgn = nullptr;

  knots = 0; knots = logknots = logknotvals = nullptr; 
  splinelog = nullptr; 
  setKnotPositions(KNOTS);
  acc = gsl_interp_accel_alloc();

  nsigmaknots = 0; sigmaknots = sigmavals = nullptr; sigmainterp = nullptr;
  min_sigma = 0.0;
  setSigmaPositions(SIGMAS);
  accsigma = gsl_interp_accel_alloc();

  noffsetknots = 0; offsetknots = offsetvals = nullptr; offsetinterp = nullptr;
  setOffsetPositions(OFFSETS);
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
  fslv(nullptr) {

  nRWork = 0;
  RWork = nullptr;

  nf2work = 0;
  f2work_inv = f2work_log = nullptr;
  f2work_sgn = nullptr;

  knots = 0; knots = logknots = logknotvals = nullptr; 
  splinelog = nullptr; 
  setKnotPositions(KNOTS);
  acc = gsl_interp_accel_alloc();

  nsigmaknots = 0; sigmaknots = sigmavals = nullptr; sigmainterp = nullptr;
  min_sigma = 0.0;
  setSigmaPositions(SIGMAS);
  accsigma = gsl_interp_accel_alloc();

  noffsetknots = 0; offsetknots = offsetvals = nullptr; offsetinterp = nullptr;
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
numberCountsDoubleLogNormal(unsigned int nk, const float* const K,
                            unsigned int ns, const float* const S,
                            unsigned int no, const float* const O) :
  knots_valid(false), sigmas_valid(false), offsets_valid(false),
  knotpos_loaded(false), sigmapos_loaded(false), offsetpos_loaded(false),
  knotvals_loaded(false), sigmavals_loaded(false), offsetvals_loaded(false),
  fslv(nullptr) {

  nRWork = 0;
  RWork = nullptr;

  nf2work = 0;
  f2work_inv = f2work_log = nullptr;
  f2work_sgn = nullptr;

  knots = 0; knots = logknots = logknotvals = nullptr; splinelog = nullptr;
  setKnotPositions(nk,K);
  acc = gsl_interp_accel_alloc();

  nsigmaknots = 0; sigmaknots = sigmavals = nullptr; sigmainterp = nullptr;
  min_sigma = 0.0;
  setSigmaPositions(ns,S);
  accsigma = gsl_interp_accel_alloc();

  noffsetknots = 0; offsetknots = offsetvals = nullptr; offsetinterp = nullptr;
  setOffsetPositions(no,O);
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
  fslv(nullptr) {

  nRWork = 0;
  RWork = nullptr;

  nf2work = 0;
  f2work_inv = f2work_log = nullptr;
  f2work_sgn = nullptr;

  knots = 0; knots = logknots = logknotvals = nullptr; splinelog = nullptr;
  setKnotPositions(nk,K);
  acc = gsl_interp_accel_alloc();

  nsigmaknots = 0; sigmaknots = sigmavals = nullptr; sigmainterp = nullptr;
  min_sigma = 0.0;
  setSigmaPositions(ns,S);
  accsigma = gsl_interp_accel_alloc();

  noffsetknots = 0; offsetknots = offsetvals = nullptr; offsetinterp = nullptr;
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
  min_sigma = 0.0;
  knots = logknots = logknotvals = sigmaknots = 
    sigmavals = offsetknots = offsetvals = nullptr;
  acc = accsigma = accoffset = nullptr;
  splinelog = nullptr;
  sigmainterp = offsetinterp = nullptr;

  nRWork = 0;
  RWork = nullptr;

  nf2work = 0;
  f2work_inv = f2work_log = nullptr;
  f2work_sgn = nullptr;

  fslv = nullptr;

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
  min_sigma = other.min_sigma;
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
  if (knots != nullptr) delete[] knots;
  if (logknots != nullptr) delete[] logknots;
  if (logknotvals != nullptr) delete[] logknotvals;
  if (splinelog != nullptr) gsl_spline_free(splinelog);
  gsl_interp_accel_free(acc);

  if (sigmaknots != nullptr) delete[] sigmaknots;
  if (sigmavals != nullptr) delete[] sigmavals;
  if (sigmainterp != nullptr) gsl_interp_free(sigmainterp);
  gsl_interp_accel_free(accsigma);

  if (offsetknots != nullptr) delete[] offsetknots;
  if (offsetvals != nullptr) delete[] offsetvals;
  if (offsetinterp != nullptr) gsl_interp_free(offsetinterp);
  gsl_interp_accel_free(accoffset);

  gsl_integration_workspace_free(gsl_work);
  delete[] varr;

  if (fslv != nullptr) gsl_root_fsolver_free(fslv);
  if (RWork != nullptr) fftw_free(RWork);
  if (f2work_sgn != nullptr) fftw_free(f2work_sgn);
  if (f2work_inv != nullptr) fftw_free(f2work_inv);
  if (f2work_log != nullptr) fftw_free(f2work_log);
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
  min_sigma = other.min_sigma;
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

  If number of knots doesn't change, contents are preserved.
  Otherwise all values are destroyed (but not sigma or offset)
*/
void numberCountsDoubleLogNormal::setNKnots(unsigned int n) {
  if (nknots == n) return;

  if (knots != nullptr) delete[] knots;
  if (logknotvals != nullptr) delete[] logknotvals;
  if (logknots != nullptr) delete[] logknots;
  if (splinelog != nullptr) gsl_spline_free(splinelog);
  knots_valid = knotpos_loaded = knotvals_loaded = false;

  if (n > 0) {
    knots = new double[n];
    logknots = new double[n];
    logknotvals = new double[n];
  } else {
    knots = logknots = logknotvals = nullptr;
  }
  if (n > 2)
    splinelog=gsl_spline_alloc(gsl_interp_cspline,
                               static_cast<size_t>(n));
  else 
    splinelog=nullptr;
  nknots = n;
}

/*!
  \param[in] n Number of sigmas

  If number of sigmas doesn't change, contents are preserved.
  Otherwise all values are destroyed (but not knots or offset)
*/
void numberCountsDoubleLogNormal::setNSigmas(unsigned int n) {
  if (nsigmaknots == n) return;

  if (sigmaknots != nullptr) delete[] sigmaknots;
  if (sigmavals != nullptr) delete[] sigmavals;
  if (sigmainterp != nullptr) gsl_interp_free(sigmainterp);
  sigmas_valid = sigmapos_loaded = sigmavals_loaded = false;

  if (n > 0) {
    sigmaknots = new double[n];
    sigmavals = new double[n];
  } else {
    sigmaknots = sigmavals = nullptr;
  }
  sigmainterp = nullptr;
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

  If number of sigmas doesn't change, contents are preserved.
  Otherwise all values are destroyed (but not knots or sigmas)
*/
void numberCountsDoubleLogNormal::setNOffsets(unsigned int n) {
  if (noffsetknots == n) return;

  if (offsetknots != nullptr) delete[] offsetknots;
  if (offsetvals != nullptr) delete[] offsetvals;
  if (offsetinterp != nullptr) gsl_interp_free(offsetinterp);
  offsets_valid = offsetpos_loaded = offsetvals_loaded = false;

  if (n > 0) {
    offsetknots = new double[n];
    offsetvals = new double[n];
  } else {
    offsetknots = offsetvals = nullptr;
  }
  offsetinterp = nullptr;
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
setKnotPositions(const std::vector<float>& S) {
  unsigned int n = S.size();
  setNKnots(n); //Also invalidates everything
  for (unsigned int i = 0; i < nknots; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsDoubleLogNormal", "setKnotPositions",
                         "Non-positive knot positions not allowed");
  // Knot positions -must- increase, the GSL requires
  for (unsigned int i = 1; i < nknots; ++i)
    if (S[i] <= S[i-1]) {
      std::stringstream errstr;
      errstr << "Knot positions not monotonically increasing: S[" << i-1 
             << "]=" << S[i-1] << " S[" << i <<"]=" << S[i] << std::endl;
      throw affineExcept("numberCountsDoubleLogNormal", "setKnotPositions",
                         errstr.str());
    }
  for (unsigned int i = 0; i < nknots; ++i)
    knots[i] = static_cast<double>(S[i]);
  for (unsigned int i = 0; i < nknots; ++i)
    logknots[i] = log2(knots[i]);
  if (nknots > 0) knotpos_loaded = true;
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
      throw affineExcept("numberCountsDoubleLogNormal", "setKnotPositions",
                         "Non-positive knot positions not allowed");
  // Knot positions -must- increase, the GSL requires
  for (unsigned int i = 1; i < nknots; ++i)
    if (S[i] <= S[i-1]) {
      std::stringstream errstr;
      errstr << "Knot positions not monotonically increasing: S[" << i-1 
             << "]=" << S[i-1] << " S[" << i <<"]=" << S[i] << std::endl;
      throw affineExcept("numberCountsDoubleLogNormal", "setKnotPositions",
                         errstr.str());
    }
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
                                                   const float* const S) {
  setNKnots(n); //Also Invalidates everything
  for (unsigned int i = 0; i < nknots; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsDoubleLogNormal","setKnotPositions",
                         "Non-positive knot positions not allowed");
  for (unsigned int i = 1; i < nknots; ++i)
    if (S[i] <= S[i-1]) {
      std::stringstream errstr;
      errstr << "Knot positions not monotonically increasing: S[" << i-1 
             << "]=" << S[i-1] << " S[" << i <<"]=" << S[i] << std::endl;
      throw affineExcept("numberCountsDoubleLogNormal", "setKnotPositions",
                         errstr.str());
    }
  for (unsigned int i = 0; i < nknots; ++i)
    knots[i] = static_cast<double>(S[i]);
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
      throw affineExcept("numberCountsDoubleLogNormal", "setKnotPositions",
                         "Non-positive knot positions not allowed");
  for (unsigned int i = 1; i < nknots; ++i)
    if (S[i] <= S[i-1]) {
      std::stringstream errstr;
      errstr << "Knot positions not monotonically increasing: S[" << i-1 
             << "]=" << S[i-1] << " S[" << i <<"]=" << S[i] << std::endl;
      throw affineExcept("numberCountsDoubleLogNormal", "setKnotPositions",
                         errstr.str());
    }
  for (unsigned int i = 0; i < nknots; ++i)
    knots[i] = S[i];
  for (unsigned int i = 0; i < nknots; ++i)
    logknots[i] = log2(knots[i]);
  if (nknots > 0) knotpos_loaded = true;
}

/*!
  \param[in] i Index into knot positions
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
setSigmaPositions(const std::vector<float>& S) {
  unsigned int n = S.size();
  setNSigmas(n); //Also invalidates
  for (unsigned int i = 0; i < n; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsDoubleLogNormal", "setSigmaPositions",
                         "Non-positive sigma knot positions not allowed");
  for (unsigned int i = 1; i < nsigmaknots; ++i)
    if (S[i] <= S[i-1]) {
      std::stringstream errstr;
      errstr << "Sigma positions not monotonically increasing: S[" << i-1 
             << "]=" << S[i-1] << " S[" << i <<"]=" << S[i] << std::endl;
      throw affineExcept("numberCountsDoubleLogNormal", "setSigmaPositions",
                         errstr.str());
    }
  for (unsigned int i = 0; i < n; ++i)
    sigmaknots[i] = static_cast<double>(S[i]);
  if (nsigmaknots > 0) sigmapos_loaded = true;
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
      throw affineExcept("numberCountsDoubleLogNormal", "setSigmaPositions",
                         "Non-positive sigma knot positions not allowed");
  for (unsigned int i = 1; i < nsigmaknots; ++i)
    if (S[i] <= S[i-1]) {
      std::stringstream errstr;
      errstr << "Sigma positions not monotonically increasing: S[" << i-1 
             << "]=" << S[i-1] << " S[" << i <<"]=" << S[i] << std::endl;
      throw affineExcept("numberCountsDoubleLogNormal", "setSigmaPositions",
                         errstr.str());
    }
  for (unsigned int i = 0; i < n; ++i)
    sigmaknots[i] = S[i];
  if (nsigmaknots > 0) sigmapos_loaded = true;
}

/*!
  \param[in] n Number of knots for band 2 color model sigma
  \param[in] S Input knot positions
*/
void numberCountsDoubleLogNormal::setSigmaPositions(unsigned int n, 
                                                    const float* const S) {
  setNSigmas(n); //Invalidates
  for (unsigned int i = 0; i < n; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsDoubleLogNormal", 
                         "setSigmaPositions",
                                           "Non-positive sigma knot positions not allowed");
  for (unsigned int i = 1; i < nsigmaknots; ++i)
    if (S[i] <= S[i-1]) {
      std::stringstream errstr;
      errstr << "Sigma positions not monotonically increasing: S[" << i-1 
             << "]=" << S[i-1] << " S[" << i <<"]=" << S[i] << std::endl;
      throw affineExcept("numberCountsDoubleLogNormal", 
                         "setSigmaPositions", errstr.str());
    }
  for (unsigned int i = 0; i < n; ++i)
    sigmaknots[i] = static_cast<double>(S[i]);
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
      throw affineExcept("numberCountsDoubleLogNormal", "setKnots",
                         "Non-positive sigma knot positions not allowed");
  for (unsigned int i = 1; i < nsigmaknots; ++i)
    if (S[i] <= S[i-1]) {
      std::stringstream errstr;
      errstr << "Sigma positions not monotonically increasing: S[" << i-1 
             << "]=" << S[i-1] << " S[" << i <<"]=" << S[i] << std::endl;
      throw affineExcept("numberCountsDoubleLogNormal", "setSigmaPositions",
                         errstr.str());
    }
  for (unsigned int i = 0; i < n; ++i)
    sigmaknots[i] = S[i];
  if (nsigmaknots > 0) sigmapos_loaded = true;
}

/*!
  \param[in] mins New minimum sigma value
*/
void numberCountsDoubleLogNormal::setMinSigma(double mins) {
  if (mins < 0)
    throw affineExcept("numberCountsDoubleLogNormal", "setMinSigma",
                       "Invalid (negative) min Sigma");
  min_sigma = mins;
}

/*!
  \param[in] i Index into sigma knot positions
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
setOffsetPositions(const std::vector<float>& S) {
  unsigned int n = S.size();
  setNOffsets(n);
  for (unsigned int i = 0; i < n; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsDoubleLogNormal", "setOffsetPositions",
                         "Non-positive offset knot positions not allowed");
  for (unsigned int i = 1; i < noffsetknots; ++i)
    if (S[i] <= S[i-1]) {
      std::stringstream errstr;
      errstr << "Offset positions not monotonically increasing: S[" << i-1 
             << "]=" << S[i-1] << " S[" << i <<"]=" << S[i] << std::endl;
      throw affineExcept("numberCountsDoubleLogNormal", "setOffsetPositions",
                         errstr.str());
    }
  for (unsigned int i = 0; i < n; ++i)
    offsetknots[i] = static_cast<double>(S[i]);
  if (noffsetknots > 0) offsetpos_loaded = true;
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
      throw affineExcept("numberCountsDoubleLogNormal", "setOffsetPositions",
                         "Non-positive offset knot positions not allowed");
  for (unsigned int i = 1; i < noffsetknots; ++i)
    if (S[i] <= S[i-1]) {
      std::stringstream errstr;
      errstr << "Offset positions not monotonically increasing: S[" << i-1 
             << "]=" << S[i-1] << " S[" << i <<"]=" << S[i] << std::endl;
      throw affineExcept("numberCountsDoubleLogNormal", "setOffsetPositions",
                         errstr.str());
    }
  for (unsigned int i = 0; i < n; ++i)
    offsetknots[i] = S[i];
  if (noffsetknots > 0) offsetpos_loaded = true;
}

/*!
  \param[in] i Index into offset knot positions
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
                                                     const float* const S) {
  setNOffsets(n);
  for (unsigned int i = 0; i < n; ++i)
    if (S[i] <= 0.0)
      throw affineExcept("numberCountsDoubleLogNormal", "setKnots",
                         "Non-positive offset knot positions not allowed");
  for (unsigned int i = 1; i < noffsetknots; ++i)
    if (S[i] <= S[i-1]) {
      std::stringstream errstr;
      errstr << "Offset positions not monotonically increasing: S[" << i-1 
             << "]=" << S[i-1] << " S[" << i <<"]=" << S[i] << std::endl;
      throw affineExcept("numberCountsDoubleLogNormal", "setOffsetPositions",
                         errstr.str());
    }
  for (unsigned int i = 0; i < n; ++i)
    offsetknots[i] = static_cast<double>(S[i]);
  if (noffsetknots > 0) offsetpos_loaded = true;
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
      throw affineExcept("numberCountsDoubleLogNormal", "setKnots",
                         "Non-positive offset knot positions not allowed");
  for (unsigned int i = 1; i < noffsetknots; ++i)
    if (S[i] <= S[i-1]) {
      std::stringstream errstr;
      errstr << "Offset positions not monotonically increasing: S[" << i-1 
             << "]=" << S[i-1] << " S[" << i <<"]=" << S[i] << std::endl;
      throw affineExcept("numberCountsDoubleLogNormal", "setOffsetPositions",
                         errstr.str());
    }
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
                       "Knot positions not set");
  if (!sigmapos_loaded) 
    throw affineExcept("numberCountsDoubleLogNormal", "getPositions",
                       "Sigma knot positions not set");
  if (!offsetpos_loaded) 
    throw affineExcept("numberCountsDoubleLogNormal", "getPositions",
                       "Offset knot positions not set");

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
setPositions(const std::vector<float>& K, const std::vector<float>& S,
             const std::vector<float>& O) {
  setKnotPositions(K);
  setSigmaPositions(S);
  setOffsetPositions(O);
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

//////////////////////////////////////////////////////////////////////////
//The stuff below here is specific to this particular model, the stuff above
// may be useful for an abstract base class of band1 spline plus color models.
// Currently this isn't implemented because we only really have this one model,
// so adding the subclass infrastructure would just slow down the code.

/*!
  \param[out] F Parameters to get from model
  
  Will set the first nknot + nsigma + noffset parameters, ignoring any others.
  If the knot values aren't loaded, sets them to NaN.
*/
void numberCountsDoubleLogNormal::getParams(paramSet& F) const {
  unsigned int nneeded = nknots + nsigmaknots + noffsetknots;
  if (F.getNParams() < nneeded)
    throw affineExcept("numberCountsDoubleLogNormal", "getKnots",
                       "Not enough space in output variable");
  if (knotvals_loaded && offsetvals_loaded && sigmavals_loaded) {
    //External storage is base 10 (and float)
    for (unsigned int i = 0; i < nknots; ++i)
      F[i] = static_cast<float>(pofd_mcmc::ilogfac * logknotvals[i]); 
    for (unsigned int i = 0; i < nsigmaknots; ++i)
      F[i + nknots] = static_cast<float>(sigmavals[i]);
    for (unsigned int i = 0; i < noffsetknots; ++i)
      F[i + nknots + nsigmaknots] = static_cast<float>(offsetvals[i]);
  } else 
    for (unsigned int i = 0; i < nneeded; ++i)
      F[i] = std::numeric_limits<float>::quiet_NaN();
}

/*!
  \param[in] F Parameters to set in model

  This will ignore any parameters beyond those it needs.
*/
void numberCountsDoubleLogNormal::setParams(const paramSet& F) {
  unsigned int nneeded = nknots + nsigmaknots + noffsetknots;
  if (nneeded > F.getNParams())
    throw affineExcept("numberCountsDoubleLogNormal", "setParams",
                       "Not enough parameters present to set");
  if (!knotpos_loaded)
    throw affineExcept("numberCountsDoubleLogNormal", "setParams",
                       "Knot positions not set");
  if (!sigmapos_loaded)
    throw affineExcept("numberCountsDoubleLogNormal", "setParams",
                       "Sigma knot positions not set");
  if (!offsetpos_loaded)
    throw affineExcept("numberCountsDoubleLogNormal", "setParams",
                       "Offset knot positions not set");

  // Internal storage is log2 and double, inputs are log10 float
  for (unsigned int i = 0; i < nknots; ++i)
    logknotvals[i] = pofd_mcmc::logfac * static_cast<double>(F[i]); 
  if (nknots > 2) 
    gsl_spline_init(splinelog, logknots, logknotvals,
                    static_cast<size_t>(nknots));
  else throw affineExcept("numberCountsDoubleLogNormal", "setParams",
                          "Need at least 3 knots");
  knotvals_loaded = true;
  checkKnotsValid();

  for (unsigned int i = nknots; i < nknots+nsigmaknots; ++i)
    sigmavals[i-nknots] = static_cast<double>(F[i]);
  if (nsigmaknots > 1)
    gsl_interp_init(sigmainterp, sigmaknots, sigmavals,
                    static_cast<size_t>(nsigmaknots));
  sigmavals_loaded = true;
  checkSigmasValid();
  
  unsigned int istart = nknots+nsigmaknots;
  for (unsigned int i = istart; i < istart+noffsetknots; ++i)
    offsetvals[i-istart] = static_cast<double>(F[i]);
  if (noffsetknots > 1)
    gsl_interp_init(offsetinterp, offsetknots, offsetvals,
                    static_cast<size_t>(noffsetknots));
  offsetvals_loaded = true;
  checkOffsetsValid();
}


/*!
  This isn't a free thing to check, so we use this by storing
  the results whenever the knots are set rather than every
  time we call the model.
*/
void numberCountsDoubleLogNormal::checkKnotsValid() const noexcept {
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
  This isn't a free thing to check, so we use this by storing
  the results whenever the sigmas are set rather than every
  time we call the model.
*/
void numberCountsDoubleLogNormal::checkSigmasValid() const noexcept {
  sigmas_valid = false;
  if (!sigmapos_loaded) return;
  if (!sigmavals_loaded) return;
  if (nsigmaknots == 0) return;
  if (min_sigma < 0.0) return;
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
  
  This isn't a free thing to check, so we use this by storing
  the results whenever the offsets are set rather than every time
  we do anything.
*/
void numberCountsDoubleLogNormal::checkOffsetsValid() const noexcept {
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
bool numberCountsDoubleLogNormal::isValid() const noexcept {
  if (!(knots_valid && sigmas_valid && offsets_valid)) return false;
  return true;
}

/*!
  \param[in] p1 First parameter set
  \param[in] p2 Second parameter set
  \returns Relative distance between p1, p2 over entries this model cares 
     about.  This is defined as the sqrt of the sum differences
     between the two parameter sets divided by the sqrt of the
     length of the first parameter (or, if the lenght is 0, just
     the sqrt of the sum of the squared differences).
*/
float numberCountsDoubleLogNormal::paramRelativeDistance(const paramSet& p1, 
                                                         const paramSet& p2) 
  const throw(affineExcept) {

  unsigned int ntot = getNParams();

  if (ntot == 0) return std::numeric_limits<double>::quiet_NaN();
  if (p1.getNParams() < ntot)
    throw affineExcept("numberCountsDoubleLogNormal", 
                       "paramRelativeDistance",
                       "paramSet 1 doesn't have enough entries");
  if (p2.getNParams() < ntot)
    throw affineExcept("numberCountsDoubleLogNormal", 
                       "paramRelativeDistance",
                       "paramSet 2 doesn't have enough entries");
  if (&p1 == &p2) return 0.0; // Same object!
  float val = p1[0];
  float diffval = val - p2[0];
  float lensum = val * val;
  float diffsum = diffval * diffval;
  for (unsigned int i = 1; i < ntot; ++i) {
    val = p1[i];
    lensum += val * val;
    diffval = val - p2[i];
    diffsum += diffval * diffval;
  }
  if (lensum == 0)
    return sqrt(diffsum);
  return sqrt(diffsum / lensum);
}

/*!
  \param[in] f1 Band 1 flux value to evaluate sigma for
  \returns Value of band 2 color model sigma at f1

  Assumes validity already checked
*/
double numberCountsDoubleLogNormal::getSigmaInner(double f1) const noexcept {
  if (nsigmaknots == 1) return sigmavals[0];
  if (f1 <= sigmaknots[0]) return sigmavals[0];
  if (f1 >= sigmaknots[nsigmaknots-1]) return sigmavals[nsigmaknots-1];
  double val = gsl_interp_eval(sigmainterp, sigmaknots, sigmavals, 
                               f1, accsigma);
  return (val < min_sigma) ? min_sigma : val;
}

/*!
  \param[in] f1 Band 1 flux value to evaluate offset for
  \returns Value of band 2 color model offset at f1

  Assumes validity already checked
*/
double numberCountsDoubleLogNormal::getOffsetInner(double f1) const noexcept {
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
double numberCountsDoubleLogNormal::getSigma(double f1) const noexcept {
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
double numberCountsDoubleLogNormal::getOffset(double f1) const noexcept {
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
getNumberCountsInner(double f1, double f2) const noexcept {
  const double normfac = 1.0 / sqrt(2 * M_PI);
  const double prefac = -0.5 * pofd_mcmc::ilog2toe;
  if (f1 < knots[0] || f1 >= knots[nknots-1] || f2 <= 0.0) 
    return 0.0; //Out of range

  //This is the n_1 bit.  Not zero because of the above check.
  double cnts = exp2(gsl_spline_eval(splinelog, log2(f1), acc)); 

  //Counts in band 2, Log Normal in f2/f1, multiply them onto n_1
  double if1 = 1.0 / f1;
  double isigma = 1.0 / getSigmaInner(f1);
  double tfac = (log(f2 * if1) - getOffsetInner(f1)) * isigma;
  //yes, it's 1/f2 here
  cnts *= normfac * isigma * exp2(prefac * tfac * tfac) / f2; 
  return cnts;
}

/*!
  \param[in] f1 Flux density in band 1
  \param[in] f2 Flux density in band 2
  \returns Number counts at f1, f2

  Does validity checks on input.
*/
double numberCountsDoubleLogNormal::getNumberCounts(double f1, double f2) 
  const noexcept {
  if ((nknots < 2) || (nsigmaknots < 1) || (noffsetknots < 1))
    return std::numeric_limits<double>::quiet_NaN();
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  if (std::isnan(f1) || std::isinf(f1)) 
    return std::numeric_limits<double>::quiet_NaN();
  if (std::isnan(f2) || std::isinf(f2)) 
    return std::numeric_limits<double>::quiet_NaN();
  return getNumberCountsInner(f1, f2);
}

/*!
  \param[in] f1 Flux density in band 1
  \returns Number counts in band 1 at f1: \f$dN / dS_1 \left(f1\right) \f$
*/
double numberCountsDoubleLogNormal::getBand1NumberCounts(double f1) 
  const noexcept {
  if ((nknots < 2) || (nsigmaknots < 1) || (noffsetknots < 1))
    return std::numeric_limits<double>::quiet_NaN();
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  if (std::isnan(f1) || std::isinf(f1)) 
    return std::numeric_limits<double>::quiet_NaN();
  if (f1 < knots[0] || f1 >= knots[nknots-1]) 
    return 0.0; //Out of range
  
  // The model is set up so that dN / dS_1 is just n1
  return exp2(gsl_spline_eval(splinelog, log2(f1), acc)); 
}

/*!
  \param[in] f2 Flux density in band 2
  \returns Number counts in band 2 at f2: \f$dN / dS_2 \left(f2\right) \f$
*/
double numberCountsDoubleLogNormal::getBand2NumberCounts(double f2) 
  const noexcept {

  const double normfac = 1.0 / sqrt(2 * M_PI);

  if ((nknots < 2) || (nsigmaknots < 1) || (noffsetknots < 1))
    return std::numeric_limits<double>::quiet_NaN();
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  if (std::isnan(f2) || std::isinf(f2)) 
    return std::numeric_limits<double>::quiet_NaN();
  if (f2 <= 0.0) return 0.0; // Out of range

  // dN / dS_1 was simple, but dN / dS_2 is not.  We have to integrate numerically
  //  and that means the actual computation is done in another function (evalCounts)
  //  that can be passed to the GSL integration routines.  This is just the nice(r)
  //  user interface to that

  double result, error;
  void *params;
  gsl_function F;
  double minknot = knots[0];
  double maxknot = knots[nknots-1];
  unsigned int noff = noffsetknots; // Copy to avoid casting errors
  unsigned int nsig = nsigmaknots;
  double mins = min_sigma;
  double f2_internal = f2;
  //  We use the same varr as splineInt, but pack stuff in differently
  varr[0] = static_cast<void*>(&f2_internal);
  varr[1] = static_cast<void*>(splinelog);
  varr[2] = static_cast<void*>(acc);
  varr[3] = static_cast<void*>(&minknot);
  varr[4] = static_cast<void*>(&maxknot);
  varr[5]  = static_cast<void*>(offsetinterp);
  varr[6]  = static_cast<void*>(accoffset);
  varr[7]  = static_cast<void*>(&noff);
  varr[8] = static_cast<void*>(offsetknots);
  varr[9] = static_cast<void*>(offsetvals);
  varr[10] = static_cast<void*>(sigmainterp);
  varr[11] = static_cast<void*>(accsigma);
  varr[12] = static_cast<void*>(&nsig);
  varr[13] = static_cast<void*>(sigmaknots);
  varr[14] = static_cast<void*>(sigmavals);
  varr[15] = static_cast<void*>(&mins);

  params = static_cast<void*>(varr);

  F.function = &evalCounts;
  F.params = params;
  // This integrates dN / dS1 dS2 over S1 for a specific value of S2 
  //  (which is stored in params[0])
  gsl_integration_qag(&F, minknot, maxknot, 0, 1e-5, 1000,
                      GSL_INTEG_GAUSS41, gsl_work, &result, &error); 

  // The 1 / sqrt(2 pi) bit isn't included in evalCounts, so add it here
  return result * normfac;
}

/*
  \returns The minimum flux supported by the model in each band as a pair

  Fluxes must be positive, so return smallest non-zero double value.
 */
dblpair numberCountsDoubleLogNormal::getMinFlux() const noexcept {
  if (nknots == 0) 
    return std::make_pair(std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::quiet_NaN());
  return std::make_pair(knots[0], std::numeric_limits<double>::min());
}

/*
  \returns The maximum flux supported by the model in each band

  This is not entirely defined for band 2, so return a value a few sigma
  above the top value.  It is well defined for band 1.
 */
dblpair numberCountsDoubleLogNormal::getMaxFlux() const noexcept {
  const double nsig = 3.0; // Gaussian equivalent sigma for band 2 solver
  if (nknots == 0)     
    return std::make_pair(std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::quiet_NaN());
  double kv = knots[nknots - 1];

  //Get the S2/S1 at the top S1 knot that represents nsig down
  // in equivalent Gaussian terms
  double sg = getSigmaInner(kv);
  double mu = getOffsetInner(kv);
  double m2 = logNormalSolver(mu, sg, nsig);

  return std::make_pair(kv, kv * m2);
}

/*!
  \param[in] mu Mean parameter of Log-Normal distribution (not the actual mean)
  \param[in] sig Sigma parameter of Log-Normal distribution (not the actual 
                  sigma)
  \param[in] nsig Number of Gaussian equivalent sigma down from peak
  \returns Value of log Normal argument that is the required amount
            down from the peak.  On the positive side of the mode.
*/
double numberCountsDoubleLogNormal::logNormalSolver(double mu, double sig, 
                                                    double nsig) const {

  const unsigned int max_bracket_iters = 20;
  const double prefac = 1.0 / sqrt(mcmc_affine::two_pi);
  const double tol = 1e-3; // Bracket tolerance

  // Argument checks
  if (sig <= 0.0)
    throw affineExcept("numberCountsDoubleLogNormal", "logNormalSolver",
                       "Invalid (non-positive) sigma");
  if (nsig < 0.0)
    throw affineExcept("numberCountsDoubleLogNormal", "logNormalSolver",
                       "Invalid (negative) nsig");
  double sig2 = sig * sig;
    
  // We work in terms of u = log(x) to avoid issues of negative x
  // First, find the peak
  double mode_u = mu - sig2;
  // Quick return if asking for peak
  if (nsig == 0) return exp(mode_u);
  double isig = 1.0 / sig;
  double val_mode = prefac * isig * exp(0.5 * sig2 - mu);
  
  // Now find the target value, which is the equivalent of nsig
  //  sigma from the peak for a Gaussian
  double target = val_mode * exp(-0.5 * nsig * nsig);

  // Set up interface to GSL solver -- first, the GSL function
  gsl_function F;
  lognorm_params params = {mu, sig2, isig, isig * isig, target};
  F.function = &lognorm;
  F.params = &params;

  // Set up bracket.  Lower limit is the mode, upper limit has
  // to be iterated too.
  void *p = static_cast<void*>(&params);
  double u_low = mode_u;
  double delta_u = log(exp(mode_u) + 2 * sig) - mode_u;
  double fval = lognorm(u_low + delta_u, p);
  if (fval >= 0.0) {
    // Haven't bracketed yet, must increase search size
    for (unsigned int i = 0; i < max_bracket_iters; ++i) {
      delta_u *= 1.5;
      fval = lognorm(u_low + delta_u, p);
      if (fval < 0.0) break;
    }
    if (fval >= 0.0) {
      std::stringstream errstr;
      errstr << "Unable to bracket Log Normal with mu: " << mu
             << " sigma: " << sig << " mode: " << exp(mode_u)
             << " target value: " << target;
      throw affineExcept("numberCountsDoubleLogNormal", "logNormalSolver",
                         errstr.str());
    }
  }
  double u_high = u_low + delta_u;

  // Allocate solver if not previously allocated
  if (fslv == nullptr) fslv = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
  gsl_root_fsolver_set (fslv, &F, u_low, u_high);
  
  // Solve
  int status = GSL_CONTINUE;
  double u;
  for (unsigned int i = 0; i < max_bracket_iters; ++i) {
    gsl_root_fsolver_iterate(fslv);
    u = gsl_root_fsolver_root(fslv);
    u_low = gsl_root_fsolver_x_lower(fslv);
    u_high = gsl_root_fsolver_x_upper(fslv);
    status = gsl_root_test_interval(u_low, u_high, 0, tol);
    if (status == GSL_SUCCESS) break;
  }
  if (status != GSL_SUCCESS) {
    std::stringstream errstr;
    errstr << "Unable to brent bracket Log Normal with mu: " << mu
           << " sigma: " << sig << " mode: " << exp(mode_u)
           << " target value: " << target;
    throw affineExcept("numberCountsDoubleLogNormal", "logNormalSolver",
                       errstr.str());
  }

  return exp(u);
}

/*!
  \param[in] bm Beam to compute range for
  \returns A pair of min/max fluxes over which R is expected to be nonzero.

  Doesn't check for model validity.  Not cleanly defined for the second band.
  We also include flux=0, even though R is technically zero there.
*/
std::pair<dblpair, dblpair> 
numberCountsDoubleLogNormal::getRRangeInternal(const doublebeam& bm) const 
  throw(affineExcept) {

  double minflux1, maxflux1, minflux2, maxflux2;
  minflux1 = maxflux1 = minflux2 = maxflux2 = 0.0;
  
  dblpair maxf = getMaxFlux();
  double cval;
  if (bm.hasSign(0)) {
    // pos/pos bit.  Affects maximum values in bands 1 and 2
    cval = maxf.first * bm.getMinMax1(0).second;
    if (cval > maxflux1) maxflux1 = cval;
    cval = maxf.second * bm.getMinMax2(0).second;
    if (cval > maxflux2) maxflux2 = cval;
  }
  if (bm.hasSign(1)) {
    //pos/neg bit.  Affects maximum in band 1, minimum in band 2
    cval = maxf.first * bm.getMinMax1(1).second;
    if (cval > maxflux1) maxflux1 = cval;
    cval = - maxf.second * bm.getMinMax2(1).second;
    if (cval < minflux2) minflux2 = cval;
  }
  if (bm.hasSign(2)) {
    //neg/pos bit.  Affects minimum in band 1, maximum in band 2
    cval = - maxf.first * bm.getMinMax1(2).second;
    if (cval < minflux1) minflux1 = cval;
    cval = maxf.second * bm.getMinMax2(2).second;
    if (cval > maxflux2) maxflux2 = cval;
  }
  if (bm.hasSign(3)) {
    //neg/neg bit.  Affects minimum in bands 1 and 2
    cval = - maxf.first * bm.getMinMax1(3).second;
    if (cval < minflux1) minflux1 = cval;
    cval = - maxf.second * bm.getMinMax2(3).second;
    if (cval < minflux2) minflux2 = cval;
  }

  return std::make_pair(std::make_pair(minflux1, maxflux1),
                        std::make_pair(minflux2, maxflux2));
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
  If you want the mean flux to some power per area, however, you 
  don't divide.

  \param[in] alpha   Power of flux in band 1
  \param[in] beta    Power of flux in band 2
  \returns Integral

  The function evaluation is done in evalPowfNDoubleLogNormal, 
  which is the thing that is passed to the GSL integration routines.
*/
double numberCountsDoubleLogNormal::splineInt(double alpha, double beta) 
  const noexcept {
  if (nknots < 2) return std::numeric_limits<double>::quiet_NaN();
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  double result, error;
  void *params;
  gsl_function F;
  double minknot = knots[0];
  double maxknot = knots[nknots-1];
  unsigned int noff = noffsetknots; // Copy to avoid casting issues
  unsigned int nsig = nsigmaknots;
  double mins = min_sigma;
  
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
    varr[17] = static_cast<void*>(&mins);
  }

  params = static_cast<void*>(varr);

  F.function = &evalPowfNDoubleLogNormal;
  F.params = params;

  gsl_integration_qag(&F, minknot, maxknot, 0, 1e-5, 1000,
                      GSL_INTEG_GAUSS41, gsl_work, &result, &error); 
  return result;
}

  
/*!
  \param[in] bm Beam to compute range for
  \returns A pair of min/max fluxes over which R is expected to be nonzero.

  We also include flux=0, even though R is technically zero there.
  This is well defined in the first band, not as well in the second,
  but this is still useful.
*/
std::pair<dblpair, dblpair> 
numberCountsDoubleLogNormal::getRRange(const doublebeam& bm) const 
  throw(affineExcept) {

  if (!isValid())
    throw affineExcept("numberCountsDoubleLogNormal", "getRRange",
                       "Model is not valid");
  if (!bm.hasData())
    throw affineExcept("numberCountsDoubleLogNormal", "getRRange",
                       "Beam is empty");

  return getRRangeInternal(bm);
}

/*!
  \returns The Total number of sources per area, however area is defined
           by the model.
*/
double numberCountsDoubleLogNormal::getNS() const noexcept {
  return splineInt(0.0,0.0);
}

/*!
  \param[in] band Which band to select; 0 for band 1, 1 for band 2
  \returns The flux per area -- however the model defines area and
           flux
*/
double numberCountsDoubleLogNormal::
getFluxPerArea(unsigned int band) const noexcept {
  if (band == 0) return splineInt(1.0, 0.0);
  else if (band == 1) return splineInt(0.0, 1.0);
  else return std::numeric_limits<double>::quiet_NaN();
}

/*!
  \param[in] band Which band to select; 0 for band 1, 1 for band 2
  \returns The flux density^2 per area -- however the model defines area and
           flux
*/
double numberCountsDoubleLogNormal::
getFluxSqPerArea(unsigned int band) const noexcept {
  if (band == 0) return splineInt(2.0, 0.0);
  else if (band == 1) return splineInt(0.0, 2.0);
  else return std::numeric_limits<double>::quiet_NaN();
}

/*!
  \param[in] p1 Power of band 1 flux density
  \param[in] p2 Power of band 2 flux density
  \returns The flux density (band1)^p1 * flux density (band2)^p2 per area
*/
double numberCountsDoubleLogNormal::
getFluxPowerPerArea(double p1, double p2) const noexcept {
  return splineInt(p1, p2);
}

/*!
  \param[in] f1   Flux 1
  \param[in] f2   Flux 2
  \param[in] bm   Beam
  \returns R 
*/
double numberCountsDoubleLogNormal::getR(double f1, double f2,
                                         const doublebeam& bm) const {
  
  // Validity
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  if (!bm.hasData()) return std::numeric_limits<double>::quiet_NaN();

  // Quick returns
  if ((f1 == 0) || (f2 == 0)) return 0.0;
  unsigned int sgn = signComp(f1, f2); // Figure out which sign component
  if (!bm.hasSign(sgn)) return 0.0; 

  double ieta1, ieta2, retval, af1, af2;
  retval = 0.0;
  af1 = fabs(f1);
  af2 = fabs(f2);
  if (bm.isHistogrammed(sgn)) {
    unsigned int npsf = bm.getNHist(sgn);
    const double *__restrict__ warr = bm.getBinWeights(sgn);
    const double *__restrict__ iparr1 = bm.getBinVals1(sgn);
    const double *__restrict__ iparr2 = bm.getBinVals2(sgn);
    for (unsigned int i = 0; i < npsf; ++i) {
      ieta1 = iparr1[i];
      ieta2 = iparr2[i];
      retval += warr[i] * ieta1 * ieta2 * 
        getNumberCountsInner(af1 * ieta1, af2 * ieta2);
    } 
  } else {
    unsigned int npsf = bm.getNPix(sgn);
    const double *__restrict__ iparr1 = bm.getInvPixArr1(sgn);
    const double *__restrict__ iparr2 = bm.getInvPixArr2(sgn);
    for (unsigned int i = 0; i < npsf; ++i) {
      ieta1 = iparr1[i];
      ieta2 = iparr2[i];
      retval += ieta1 * ieta2 * getNumberCountsInner(af1 * ieta1, 
                                                     af2 * ieta2);
    } 
  }

  double prefac;
  prefac = bm.getPixSize() / 3600.0;  //To sq deg
  return prefac * prefac * retval;
}

/*! 
  \param[in] sz New size of R work arrays.  This should be enough
    to hold two beam components, not just one.

  Allocates R work arrays.  Only upsizes -- that is, the arrays are already
  sized and you ask for a smaller size, nothing is done.
*/
void numberCountsDoubleLogNormal::setRWorkSize(unsigned int sz) const {
  // We want to pack these in relatively efficiently.
  // So we do this logially as a sz by 3 array, where the second index
  //  are 
  //    0: n1(f1/eta1) / (eta1 * sigma(f1/eta1) or 0 if the counts are 0
  //    1: log(f1/eta1) + offset(f1/eta1)
  //    2: -0.5 / sigma(f1/eta1)**2
  // If the 0th component is 0 (e.g., n1 is 0), then the others are not set
  // Note we need to pack two beam components in: pp + pn, or np + nn.
  //  So make sure to make it big enough to hold the larger of those.

  if (sz <= nRWork) return;

  if (RWork != nullptr) fftw_free(RWork);
  if (sz > 0) {
    RWork = (double*) fftw_malloc(sizeof(double) * sz * 3);
  } else RWork = nullptr;
  nRWork = sz;
}

/*! 
  \param[in] sz New size of f2 work arrays.  This should hold the
      largest f2 size you plan on using.

  Allocates f2 work arrays.  Only upsizes -- that is, the arrays are already
  sized and you ask for a smaller size, nothing is done.
*/
void numberCountsDoubleLogNormal::setf2WorkSize(unsigned int sz) const {
  if (sz <= nf2work) return;

  if (f2work_sgn != nullptr) fftw_free(f2work_sgn);
  if (f2work_inv != nullptr) fftw_free(f2work_inv);
  if (f2work_log != nullptr) fftw_free(f2work_log);
  if (sz > 0) {
    f2work_sgn = (unsigned int*) fftw_malloc(sizeof(unsigned int) * sz);
    f2work_inv = (double*) fftw_malloc(sizeof(double) * sz);
    f2work_log = (double*) fftw_malloc(sizeof(double) * sz);
  } else {
    f2work_inv = f2work_log = nullptr;
    f2work_sgn = nullptr;
  }
  nf2work = sz;
}

/*!
  \param[in] n1   Number of fluxes, band 1
  \param[in] f1   Flux in band 1, length n1
  \param[in] n2   Number of fluxes, band 2
  \param[in] f2   Flux 2, length n2
  \param[in] bm   Beam
  \param[out] R   R value, dimension n1 * n2 pre-allocated by caller.
                  Does not include pixel area factor

  Does not check inputs for validity, and does not include area prefactor.
*/
void numberCountsDoubleLogNormal::getR(unsigned int n1, const double* const f1,
                                       unsigned int n2, const double* const f2,
                                       const doublebeam& bm, 
                                       double *__restrict__ R) const {

  // Quick return
  if ((!isValid()) || (!bm.hasData()) || (n1 == 0) || (n2 == 0)) {
    for (unsigned int i = 0; i < n1 * n2; ++i)
      R[i] = std::numeric_limits<double>::quiet_NaN();
    return;
  }

  double minknot = knots[0];
  double maxknot = knots[nknots-1];

  // This computation is done in a more optimized way than
  //  the scalar version above.  As a result it is much more complicated.
  // First, we are going to load all the beam components into local
  //  arrays to avoid function call overhead in the inner loop.
  bool hassign[4];
  for (unsigned int i = 0; i < 4; ++i) hassign[i] = bm.hasSign(i);
  bool ishist[4];
  for (unsigned int i = 0; i < 4; ++i) ishist[i] = bm.isHistogrammed(i);
  unsigned int npsf[4];
  for (unsigned int i = 0; i < 4; ++i)
    if (hassign[i])
      npsf[i] = ishist[i] ? bm.getNHist(i) : bm.getNPix(i);
    else
      npsf[i] = 0;
  const double* __restrict__ wtptr[4];
  for (unsigned int i = 0; i < 4; ++i) 
    wtptr[i] = ishist[i] ? bm.getBinWeights(i) : nullptr;
  const double* __restrict__ ibmptr1[4];
  for (unsigned int i = 0; i < 4; ++i)
    ibmptr1[i] = ishist[i] ? bm.getBinVals1(i) : bm.getInvPixArr1(i);
  const double* __restrict__ logratioptr[4]; // Pointer to log(eta1/eta2)
  for (unsigned int i = 0; i < 4; ++i)
    if (hassign[i])
      logratioptr[i] = ishist[i] ? bm.getBinLogRatio(i) : bm.getLogRatio(i);
  // If abs(f1) is bigger than this, R is zero for this component
  double maxf1[4];
  for (unsigned int i = 0; i < 4; ++i)
    maxf1[i] = hassign[i] ? maxknot * bm.getMinMax1(i).second : 0.0;

  // Make sure there's enough room to do the computation
  // We need to be able to hold the larger of pp + pn or np + nn
  unsigned int maxnpsf = std::max(npsf[0] + npsf[1],
                                  npsf[2] + npsf[3]);
  setRWorkSize(maxnpsf);

  // Pixel area factor and Lognorm factor
  const double normfac = 1.0 / sqrt(2 * M_PI);
  double pixsizedeg = bm.getPixSize() / 3600.0; // To sq deg
  double prefac = pixsizedeg * pixsizedeg * normfac;

  // Main computation
  // The basic framework for speeding this up is that many of the
  //  components of the model depend only on the flux / beam in the
  //  first band.  Thus, as we loop over the 2D grid of fluxes in ij, 
  //  for each i we try to precompute as many of the pieces used in
  //  the j loop as possible by storing them in the RWork arrays.
  //  However, there is a complication -- for each sign of f1, there
  //  are potentially two beam pieces.  For example, for f1 > 0, we
  //  have to compute both the parts from the pp and pn beam.
  //  In principle, this means that there may be some wasted effort here
  //  if, in fact, there are no f2 < 0 (pn) components in the fluxes.
  //  Rather than checking that, however, in practice this should be extremely
  //  rare because if the pn beam component is present, we are going to
  //  want to explore f2 < 0 space to compute the P(D).
  //   RWork[*,0] = n1(f1/eta1) / eta1 * sigma1(f1/eta1) * pixel size prefac, 
  //      possibly multiplied by the wts if the beam is histogrammed
  //   RWork[*,1] will hold log(f1) + offset(f1/eta1) - log(eta1 / eta2)
  //               noting that the last term is something doublebeam
  //               provides
  //   RWork[*,2] = -0.5 / (log(2) * sigma(f1/eta1)^2)
  //                Where the log(2) allows us to use exp2 rather than exp
  //                (since exp2 is faster)
  //  We only test for out of range (R always 0) on the first flux density
  //   because the model doesn't actually terminate sharply in the other band

  // Pointers into RWork array for sgn components
  double *__restrict__ RWptr, *__restrict__ cRW;
  double *__restrict__ rowptr; // Row pointer into output (R)
  const double *__restrict__ wptr; // Weights 
  const double *__restrict__ iparr1; // Inverse beam 1
  const double *__restrict__ plogratio; // log(eta1 / eta2);
  unsigned int sgn1, sgn, f2sgn; // Sign index for this component
  unsigned int curr_n; // Number of beam elements for current component
  double f1val, f2val, if2val, f1prod, ieta1, workval;
  double isigma, tfac, cts, logf1val, logf2val;

  bool hasNegX1 = hassign[2] || hassign[3]; // Any neg x1 beam components
  std::memset(R, 0, n1 * n2 * sizeof(double));

  // Pretabulated f2 values
  setf2WorkSize(n2);
  for (unsigned int j = 0; j < n2; ++j) {
    f2val = f2[j];

    if (f2val == 0) {
      f2work_inv[j] = 0.0; // Sign to skip
      continue;
    }
    
    if (f2val > 0) {
      f2work_sgn[j] = 0;
    } else if (f2val < 0) {
      f2work_sgn[j] = 1;
      f2val = fabs(f2val);
    } 

    f2work_inv[j] = 1.0  / f2val;
    f2work_log[j] = log(f2val);
  }

  for (unsigned int i = 0; i < n1; ++i) { // Loop over flux1
    f1val = f1[i];
    if ((f1val == 0) || (f1val < 0 && !hasNegX1)) continue; // R will be zero
    // Set sign and |f1| if needed
    if (f1val > 0) sgn1 = 0; else { sgn1 = 2; f1val = fabs(f1val); }
    logf1val = log(f1val);

    // Number in neg/pos bits, if present
    unsigned int np = hassign[sgn1] ? npsf[sgn1] : 0;
    unsigned int nn = hassign[sgn1 + 1] ? npsf[sgn1 + 1] : 0;

    // Loop over two sign bits (the sign on the second beam), pre-computing
    //  RWork information as described above.  Some index trickery here,
    //  so this isn't the clearest code.  But it's inner loop stuff, so
    //  clarity is of secondary importance
    for (unsigned int k = 0; k < 2; ++k) { // k = 0 is pos, k = 1 is neg
      sgn = sgn1 + k; // Current full sign (0 through 4, usual scheme)
      curr_n = (k == 0) ? np : nn;
      // Skip stuff out of range; note we set RWork[:, 0] to 0 below
      //  in the else, which we use to mark values to skip
      if (curr_n > 0 && f1val < maxf1[sgn]) {
          iparr1 = ibmptr1[sgn]; // Band 1 beam pixel array pointer
                plogratio = logratioptr[sgn]; // log (beam1 / beam2) pointer
                // &RWork[k * np, 0] -- piece of RWork for this sign
                RWptr = RWork + k * np * 3; 
                if (ishist[sgn]) {
                  // With beam hist
                  wptr = wtptr[sgn];
                  //Loop over beam pixels in band 1
                  for (unsigned int j = 0; j < curr_n; ++j) { 
                    ieta1 = iparr1[j];
                    f1prod = f1val * ieta1;
                    // cRW will be the RWork piece for this sign and beam pixel
                    cRW = RWptr + 3 * j; 
                    if (f1prod >= minknot && f1prod < maxknot) {
                      // Band 1 number counts
                      cts = exp2(gsl_spline_eval(splinelog, log2(f1prod), acc));
                      if (cts > 0) {
                              // Pieces we can pre-compute
                              isigma = 1.0 / getSigmaInner(f1prod);
                              cRW[0] = wptr[j] * ieta1 * isigma * cts * prefac;
                              cRW[1] = logf1val + getOffsetInner(f1prod) - plogratio[j];
                              cRW[2] = -0.5 * isigma * isigma * pofd_mcmc::ilog2toe;
                      } else cRW[0] = 0; // Mark as no-use
                    } else cRW[0] = 0; // Same
                  }
                } else {
                  // No beam hist; same basic plan, but no weights
                  for (unsigned int j = 0; j < curr_n; ++j) {
                    ieta1 = iparr1[j];
                    f1prod = f1val * ieta1;
                    cRW = RWptr + 3*j;
                    if (f1prod >= minknot && f1prod < maxknot) {
                      cts = exp2(gsl_spline_eval(splinelog, log2(f1prod), acc));
                      if (cts > 0) {
                              isigma = 1.0 / getSigmaInner(f1prod);
                              cRW[0] = cts * ieta1 * isigma * prefac;
                              cRW[1] = logf1val + getOffsetInner(f1prod) - plogratio[j];
                              cRW[2] = -0.5 * isigma * isigma * pofd_mcmc::ilog2toe;
                      } else cRW[0] = 0;
                    } else cRW[0] = 0;
                  }
                }
        } else {
                // Zero out Rwork[:,0] to mark as not valid, but only the
                //  part for this sign component
                unsigned int minidx = (k == 0) ? 0 : np;
                unsigned int maxidx = (k == 0) ? np : (np+nn);
                for (unsigned int j = minidx; j < maxidx; ++j)
            RWork[3 * j] = 0.0;
        }
      }

      // Now that we have pre-computed all the things that only
      //  depend on the flux in band 1, do the band 2 loop, writing
      //  the final values into R.
      rowptr = R + i * n2; // row pointer into R
      double cRW0;
      for (unsigned int j = 0; j < n2; ++j) {
        if2val = f2work_inv[j]; 
        // Testing on this is -much- faster, oddly, than testing on sgn
        if (if2val == 0) continue; // f2 = 0, R will be 0

        // sgn will be pp, pn, np, nn component
        f2sgn = f2work_sgn[j];
        sgn = sgn1 + f2sgn; // Recall sgn1 is from f1

        if (hassign[sgn]) {
                workval = 0;
                curr_n = npsf[sgn];
                RWptr = RWork + f2sgn * np * 3; //&RWork[sgnoff*sgnbreak, 0]
                logf2val = f2work_log[j];
                for (unsigned int k = 0; k < curr_n; ++k) {
                  cRW = RWptr + k * 3; //&RWork[sgnoff*sgnbreak + k, 0]
                  cRW0 = cRW[0]; //RWork[sgnoff*sgnbreak + k, 0]
                  if (cRW0 != 0) {
                  tfac = logf2val - cRW[1];
                  workval += cRW0 * exp2(tfac * tfac * cRW[2]);
                }
              }
              rowptr[j] = workval * f2work_inv[j];
      } // Recall that R was zeroed
    } // End loop over RFlux2
  } // End loop over RFlux1
}

/*!
  \param[inout] comm MPI communicator
  \param[in] dest Destination to send messages to
*/
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
  MPI_Send(const_cast<double*>(&min_sigma), 1, MPI_DOUBLE, dest,
           pofd_mcmc::NCDCSENDMINSIGMA, comm);

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
  if (offsetpos_loaded && offsetvals_loaded && noffsetknots > 0) 
    MPI_Send(offsetvals, noffsetknots, MPI_DOUBLE, dest,
             pofd_mcmc::NCDCSENDOFFSETVALS, comm);

}

/*!
  \param[inout] comm MPI communicator
  \param[in] src Where messages will come from
*/
void numberCountsDoubleLogNormal::receiveCopy(MPI_Comm comm, int src) {
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
  MPI_Recv(&min_sigma, 1, MPI_DOUBLE, src, pofd_mcmc::NCDCSENDMINSIGMA,
           comm, &Info);
  checkSigmasValid();

  //Offsets
  MPI_Recv(&loaded, 1, MPI::BOOL, src, pofd_mcmc::NCDCSENDOPLOADED,
           comm, &Info);
  MPI_Recv(&n, 1, MPI_UNSIGNED, src, pofd_mcmc::NCDCSENDNOFFSETKNOTS,
           comm, &Info);
  setNOffsets(n);
  if (loaded && noffsetknots > 0) {
    MPI_Recv(offsetknots, noffsetknots, MPI_DOUBLE, src, 
             pofd_mcmc::NCDCSENDOFFSETKNOTS, comm, &Info);
    offsetpos_loaded = loaded;
  }
  MPI_Recv(&loaded, 1, MPI::BOOL, src, pofd_mcmc::NCDCSENDOVLOADED,
           comm, &Info);
  if (offsetpos_loaded & loaded && noffsetknots > 0) {
    MPI_Recv(offsetvals, noffsetknots, MPI_DOUBLE, src, 
             pofd_mcmc::NCDCSENDOFFSETVALS, comm, &Info);
    if (n > 1)
      gsl_interp_init(sigmainterp, sigmaknots, sigmavals,
                      static_cast<size_t>(nsigmaknots));
    offsetvals_loaded = loaded;
  }
  checkOffsetsValid();
}

/*!
  \param[in] alpha Regularization multiplier.  Must be positive (not checked)
  \returns log Likelhood penalty (negative)

  Computes difference operator Tikhonov regularization penalty on
  model, where the derivative is taken in log/log space.  The penalty
  is on the sum of the squares of the differences between the slopes
  and the mean slope. It is only applied to the band 1 model.
*/
double numberCountsDoubleLogNormal::differenceRegularize(double alpha) const {
  if (!knots_valid) return std::numeric_limits<double>::quiet_NaN();
  if (nknots == 0 || nknots == 1) return 0.0;

  double log_penalty=0.0, mean_deriv=0.0;
  double delta_logknotpos, delta_logknotval, delta;
  for (unsigned int i = 1; i < nknots; ++i) {
    delta_logknotpos = logknots[i] - logknots[i - 1];
    delta_logknotval = logknotvals[i] - logknotvals[i - 1];
    mean_deriv += delta_logknotval / delta_logknotpos;
  }
  mean_deriv /= static_cast<double>(nknots-1);
  for (unsigned int i = 1; i < nknots; ++i) {
    delta_logknotpos = logknots[i] - logknots[i - 1];
    delta_logknotval = logknotvals[i] - logknotvals[i - 1];
    delta = delta_logknotval / delta_logknotpos - mean_deriv;
    log_penalty -= delta * delta;
  }
  return alpha * log_penalty;
}

/*!
  \param[inout] objid HDF5 handle to write to
  \param[in] writevals If true, write the knot values in addition to
      positions.  If the knots aren't loaded, this is ignored.
*/
void numberCountsDoubleLogNormal::writeToHDF5Handle(hid_t objid,
                                                    bool writevals) const {
  hsize_t adims;
  hid_t mems_id, att_id;

  if (H5Iget_ref(objid) < 0)
    throw affineExcept("numberCountsDoubleLogNormal", "writeToHDF5Handle",
                       "Input handle is not valid");

  // Name of model
  // Name of model
  const std::string modeltype("numberCountsDoubleLogNormal");
  hdf5utils::writeAttString(objid, "ModelType", modeltype);

  // Number of knots
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, nullptr);
  att_id = H5Acreate2(objid, "NKnots", H5T_NATIVE_UINT,
                      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &nknots);
  H5Aclose(att_id);
  att_id = H5Acreate2(objid, "NSigmaKnots", H5T_NATIVE_UINT,
                      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &nsigmaknots);
  H5Aclose(att_id);
  att_id = H5Acreate2(objid, "NOffsetKnots", H5T_NATIVE_UINT,
                      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &noffsetknots);
  H5Aclose(att_id);
  att_id = H5Acreate2(objid, "MinSigma", H5T_NATIVE_DOUBLE,
                      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &min_sigma);
  H5Aclose(att_id);
  H5Sclose(mems_id);

  // Knot positions as data
  if (knotpos_loaded)
    hdf5utils::writeDataDoubles(objid, "KnotPositions", nknots, knots);
  if (sigmapos_loaded)
    hdf5utils::writeDataDoubles(objid, "SigmaKnotPositions", nsigmaknots, 
                                sigmaknots);
  if (offsetpos_loaded)
    hdf5utils::writeDataDoubles(objid, "OffsetKnotPositions", noffsetknots, 
                                offsetknots);

  if (writevals && knotvals_loaded) {
    // Convert to log 10 for write
    double* kv;
    kv = new double[nknots];
    for (unsigned int i = 0; i < nknots; ++i)
      kv[i] = pofd_mcmc::ilogfac * logknotvals[i];
    hdf5utils::writeDataDoubles(objid, "Log10KnotValues", nknots, kv);
    delete[] kv;
  }

  if (writevals && sigmavals_loaded)
    hdf5utils::writeDataDoubles(objid, "SigmaKnotValues", nsigmaknots,
                                sigmavals);

  if (writevals && offsetvals_loaded)
    hdf5utils::writeDataDoubles(objid, "OffsetKnotValues", noffsetknots,
                                offsetvals);
}

/*!
  \param[in] objid HDF5 handle to read from

  Initializes a model from an HDF5 file, setting the knot positions
  but not values.  Checks to make sure the model type is correct
*/
void numberCountsDoubleLogNormal::readFromHDF5Handle(hid_t objid) {
  if (H5Iget_ref(objid) < 0)
    throw affineExcept("numberCountsDoubleLogNormal", 
                       "readFromHDF5Handle",
                       "Input handle is not valid");

  std::string mtype = hdf5utils::readAttString(objid, "ModelType");
  if (mtype != "numberCountsDoubleLogNormal") {
    std::stringstream errstr;
    errstr << "Unexpected model type; wanted numberCountsDoubleLogNormal"
           << " got " << mtype;
    throw affineExcept("numberCountsDoubleLogNormal", 
                       "readFromHDF5Handle", errstr.str());
  }

  unsigned int f_nknots = hdf5utils::readAttUnsignedInt(objid, "NKnots");
  setNKnots(f_nknots);
  if (f_nknots > 0) {
    // We could read directly into knots, but that could be a problem
    // for subclasses, which might want to do extra work.
    double *newknots = new double[f_nknots];
    hdf5utils::readDataDoubles(objid, "KnotPositions", f_nknots, newknots);
    setKnotPositions(f_nknots, newknots);
    delete[] newknots;
  }
  f_nknots = hdf5utils::readAttUnsignedInt(objid, "NSigmaKnots");
  setNSigmas(f_nknots);
  if (f_nknots > 0) {
    double *newknots = new double[f_nknots];
    hdf5utils::readDataDoubles(objid, "SigmaKnotPositions", 
                               f_nknots, newknots);
    setSigmaPositions(f_nknots, newknots);
    delete[] newknots;
  }
  // Note -- not reading MinSigma to keep backwards compatability
  f_nknots = hdf5utils::readAttUnsignedInt(objid, "NOffsetKnots");
  setNOffsets(f_nknots);
  if (f_nknots > 0) {
    double *newknots = new double[f_nknots];
    hdf5utils::readDataDoubles(objid, "OffsetKnotPositions", 
                               f_nknots, newknots);
    setOffsetPositions(f_nknots, newknots);
    delete[] newknots;
  }
}

/*!
  \param[inout] os Stream to write to.
  \returns True
*/
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
  \param[in] s1 Flux in band 1
  \param[in] params Model information jammed into something the GSL
     can work with.
  Internal function for use in band 2 number counts projection.
  This is just the number counts with some constants removed,
  and the idea is to integrate this over S_1
*/
static double evalCounts(double s1, void* params) noexcept {
  // Model is: n1(s1) / s1 * LogNormal(s2 / s1; mu(s1), sigma(s1))
  //  where s1 is the flux in the first band and n1 is the number
  // counts in band 1 -- see splineInt.  However, the constant 1 / sqrt(2 * pi)
  //  of LogNorm is left out, since it can be applied outside the integral
  //Params are:
  // params[0]  s_2
  // params[1]  knot interpolant (log)
  // params[2]  knot accelerator
  // params[3]  minknot
  // params[4]  maxknot
  // params[5]  offset interpolant (log)
  // params[6]  offset accelerator
  // params[7]  noffsets
  // params[8]  offset positions
  // params[9]  offset values
  // params[10] sigma interpolant (log)
  // params[11] sigma accelerator
  // params[12] nsigmas
  // params[13] sigma positions
  // params[14] sigma values
  // params[15] min_sigma
  //But this really has to be an array of pointers to void to work
  void** vptr = static_cast<void**>(params);

  // s2 quick check
  double s2 = *static_cast<double*>(vptr[0]);
  if (s2 <= 0.0) return 0.0;

  // min/max knot in band 1 for quick return if s1 is outside that
  double minknot = *static_cast<double*>(vptr[3]);
  double maxknot = *static_cast<double*>(vptr[4]);
  if (s1 < minknot || s1 >= maxknot) return 0.0;

  // Evaluate n_1(s1)
  gsl_spline* spl = static_cast<gsl_spline*>(vptr[1]);
  gsl_interp_accel* acc = static_cast<gsl_interp_accel*>(vptr[2]);
  double n1bit = exp2(gsl_spline_eval(spl, log2(s1), acc));

  // Log normal part
  // Offset
  unsigned int noffsets = *static_cast<unsigned int*>(vptr[7]);
  double *offsetpos = static_cast<double*>(vptr[8]);
  double *offsetval = static_cast<double*>(vptr[9]);
  double mu;
  if (noffsets == 1) mu = offsetval[0];
  else if (s1 <= offsetpos[0]) mu = offsetval[0];
  else if (s1 >= offsetpos[noffsets - 1]) mu = offsetval[noffsets - 1];
  else {
    gsl_interp* ospl = static_cast<gsl_interp*>(vptr[5]);
    gsl_interp_accel* oacc = static_cast<gsl_interp_accel*>(vptr[6]);
    mu = gsl_interp_eval(ospl, offsetpos, offsetval, s1, oacc);
  }

  // Sigma
  double sigma, isigma;
  unsigned int nsigmas = *static_cast<unsigned int*>(vptr[12]);
  double *sigmapos = static_cast<double*>(vptr[13]);
  double *sigmaval = static_cast<double*>(vptr[14]);
  double min_sigma = *static_cast<double*>(vptr[15]);
  if (nsigmas == 1) sigma = sigmaval[0];
  else if (s1 <= sigmapos[0]) sigma = sigmaval[0];
  else if (s1 >= sigmapos[nsigmas - 1]) sigma = sigmaval[nsigmas - 1];
  else {
    gsl_interp* sspl = static_cast<gsl_interp*>(vptr[10]);
    gsl_interp_accel* sacc = static_cast<gsl_interp_accel*>(vptr[11]);
    sigma = gsl_interp_eval(sspl, sigmapos, sigmaval, s1, sacc);
  }
  if (sigma < min_sigma) sigma = min_sigma;
  isigma = 1.0 / sigma;

  // Log normal evaluation; note that s2 and s1 are both required to be
  // positive because we explicitly check that above (for s1, that is because
  // minknot must be > 0 in this model)
  // using log(s2/s1) = - log(s1/s2) and the fact we sqr
  double s1overs2 = s1 / s2;
  double exparg = (log(s1overs2) + mu) * isigma;  
  double lnorm = isigma * s1overs2 * exp(-0.5 * exparg * exparg);
  
  return n1bit * lnorm / s1;
}

/*!
  \param[in] s1 Flux in band 1
  \param[in] params Model information jammed into something the GSL
     can work with.
  Internal function for use in model integrations.
*/
static double evalPowfNDoubleLogNormal(double s1, void* params) noexcept {
  //Model is: ( s1^power * exp( const1*mu + const2*sigma^2 ) ) * n1(s1)
  // where s1 is the flux in the first band and n1 is the number
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
  // params[17] min_sigma
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
      double min_sigma = *static_cast<double*>(vptr[17]);
      if (nsigmas == 1) sigma = sigmaval[0];
      else if (s1 <= sigmapos[0]) sigma = sigmaval[0];
      else if (s1 >= sigmapos[nsigmas-1]) sigma = sigmaval[nsigmas-1];
      else {
        gsl_interp* sspl = static_cast<gsl_interp*>(vptr[12]);
        gsl_interp_accel* sacc = static_cast<gsl_interp_accel*>(vptr[13]);
        sigma = gsl_interp_eval(sspl, sigmapos, sigmaval, s1, sacc);
      }
      if (sigma < min_sigma) sigma = min_sigma;
      expbit += const2 * sigma * sigma;
    } 

    prefac *= exp(expbit);

  } //Otherwise exp bit is just 1

  //Now multiply in n(band1)
  gsl_spline* spl = static_cast<gsl_spline*>(vptr[3]);
  gsl_interp_accel* acc = static_cast<gsl_interp_accel*>(vptr[4]);
  double splval = exp2(gsl_spline_eval(spl, log2(s1), acc));
  return prefac * splval;
}

/*!
  \param[in] u Log of log-Normal argument
  \param[in] params Parameters of log Normal
  
  This is used for root finding on a Log Normal distribution.
  It works in terms of u = log x instead of x to avoid any
  delicate issues with the LogNormal only being defined for positive x.

  The derivative is not supported because the GSL root finders that
  use derivatives are much less stable.
*/
static double lognorm(double u, void* params) {
  const double prefac = 1.0 / sqrt(2 * M_PI);

  struct lognorm_params *p = (struct lognorm_params *) params;

  double mu = p->mu;
  double isig = p->isig;
  double isig2 = p->isig2;
  double target = p->target;  // We are trying to solve L(x) == target

  double eu = exp(-u);
  double arg = u - mu;
  double L = prefac * exp(-0.5 * arg * arg * isig2) * isig * eu; // Log Norm
  return L - target;
}
