#include<sstream>
#include<cmath>
#include<cstring>
#include<limits>

#include "../include/global_settings.h"
#include "../include/PDFactory.h"
#include "../include/affineExcept.h"

const double PDFactory::subedgemult = 1e-5;

/*!
  \param[in] NINTERP number of interpolation points
*/
PDFactory::PDFactory(unsigned int NINTERP) { 
  init(NINTERP);
}

/*
  \param[in] wisfile Wisdom file filename
  \param[in] NINETERP number of interpolation points
*/
PDFactory::PDFactory(const std::string& wisfile, unsigned int NINTERP) {
  init(NINTERP);
  addWisdom(wisfile);
}

PDFactory::~PDFactory() {
  if (RinterpFlux_pos != NULL) delete[] RinterpFlux_pos;
  if (RinterpVals_pos != NULL) delete[] RinterpVals_pos;
  if (acc_pos != NULL)  gsl_interp_accel_free(acc_pos);
  if (spline_pos != NULL) gsl_spline_free(spline_pos);

  if (RinterpFlux_neg != NULL) delete[] RinterpFlux_neg;
  if (RinterpVals_neg != NULL) delete[] RinterpVals_neg;
  if (acc_neg != NULL)  gsl_interp_accel_free(acc_neg);
  if (spline_neg != NULL) gsl_spline_free(spline_neg);

  if (rvals != NULL) fftw_free(rvals);
  if (rtrans != NULL) fftw_free(rtrans);
  if (pofd != NULL) fftw_free(pofd);
  if (pval != NULL) fftw_free(pval);

  if (plan != NULL) fftw_destroy_plan(plan); 
  if (plan_inv != NULL) fftw_destroy_plan(plan_inv);
}

/*!
  \param[in] NINTERP number of interpolation points
*/
void PDFactory::init(unsigned int NINTERP) {
  currsize = 0;
  ninterp = NINTERP;

#ifdef TIMING
  resetTime();
#endif

  acc_pos = gsl_interp_accel_alloc();
  spline_pos = NULL;
  acc_neg = gsl_interp_accel_alloc();
  spline_neg = NULL;

  interpvars_allocated_pos = false;
  RinterpFlux_pos = NULL;
  RinterpVals_pos = NULL;

  interpvars_allocated_neg = false;
  RinterpFlux_neg = NULL;
  RinterpVals_neg = NULL;

  plan = plan_inv = NULL;
  rvars_allocated = false;
  rvals = NULL;
  rdflux = false;
  rtrans = NULL;
  pval = NULL;
  pofd = NULL;

  dflux = 0.0;

  verbose = false;
  has_wisdom = false;
  fftw_plan_style = FFTW_MEASURE;

  rinitialized = false;
  initialized = false;

}

#ifdef TIMING
void PDFactory::resetTime() {
  RTime = p0Time = fftTime = posTime = copyTime = normTime = 0;
  edgeTime = meanTime = logTime = 0;
}

/*!
  \param[in] nindent Number of indentation spaces
*/
void PDFactory::summarizeTime(unsigned int nindent) const {
  std::string prestring(nindent,' ');
  prestring = "  ";
    
  std::cout << "R time: " << prestring 
	    << 1.0*RTime/CLOCKS_PER_SEC << "s" << std::endl;
  std::cout << "p0 time: " << prestring 
	    << 1.0*p0Time/CLOCKS_PER_SEC << "s" << std::endl;
  std::cout << "fft time: " << prestring 
	    << 1.0*fftTime/CLOCKS_PER_SEC << "s" << std::endl;
  std::cout << "pos time: " << prestring 
	    << 1.0*posTime/CLOCKS_PER_SEC << "s" << std::endl;
  std::cout << "copy time: " << prestring 
	    << 1.0*copyTime/CLOCKS_PER_SEC << "s" << std::endl;
  std::cout << "norm time: " << prestring 
	    << 1.0*normTime/CLOCKS_PER_SEC << "s" << std::endl;
  std::cout << "edge time: " << prestring 
	    << 1.0*edgeTime/CLOCKS_PER_SEC << "s" << std::endl;
  std::cout << "mean time: " << prestring 
	    << 1.0*meanTime/CLOCKS_PER_SEC << "s" << std::endl;
  std::cout << "log time: " << prestring 
	    << 1.0*logTime/CLOCKS_PER_SEC << "s" << std::endl;
}
#endif

/*
  \param[in] NSIZE new size
  \returns True if a resize was needed
*/
//I don't check for 0 NSIZE because that should never happen
// in practice, and it isn't worth the idiot-proofing
// rtrans always gets nulled in here, then refilled when you
//  call initPD
bool PDFactory::resize(unsigned int NSIZE) {
  if (NSIZE == currsize) return false;
  freeRvars();
  currsize = NSIZE;
  if (plan != NULL) { fftw_destroy_plan(plan); plan = NULL; }
  if (plan_inv != NULL) { fftw_destroy_plan(plan_inv); plan_inv = NULL; }  
  initialized = false;
  return true;
}

// Doesn't necessarily invalidate the plans
void PDFactory::allocateRvars() {
  if (rvars_allocated) return;
  if (currsize == 0)
    throw affineExcept("PDFactory", "allocate_rvars", "Invalid (0) currsize");
  rvals = (double*) fftw_malloc(sizeof(double) * currsize);
  pofd = (double*) fftw_malloc(sizeof(double) * currsize);
  unsigned int fsize = currsize / 2 + 1;
  rtrans = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fsize);
  pval = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fsize);
  rvars_allocated = true;
}

// Doesn't necessarily invalidate the plans
void PDFactory::freeRvars() {
  if (rvals != NULL) { fftw_free(rvals); rvals=NULL; }
  if (rtrans != NULL) {fftw_free(rtrans); rtrans=NULL; }
  if (pval != NULL) { fftw_free(pval); pval = NULL; }
  if (pofd != NULL) { fftw_free(pofd); pofd = NULL; }
  rvars_allocated = false;
  initialized = false;
}

/*!
  \param[in] NINTERP New interpolation length
*/
void PDFactory::setNInterp(unsigned int NINTERP) {
  if (NINTERP == ninterp) return;
  ninterp = NINTERP;
  if (interpvars_allocated_pos) {
    freeInterpPos();
    allocateInterpPos();
  }
  if (interpvars_allocated_neg) {
    freeInterpNeg();
    allocateInterpNeg();
  }
}

void PDFactory::allocateInterpPos() {
  if (interpvars_allocated_pos) return;
  if (ninterp == 0)
    throw affineExcept("PDFactory", "allocateInterpPos",
		       "Invalid (0) ninterp");
  
  //Note that acc_pos is always allocated
  spline_pos = gsl_spline_alloc( gsl_interp_cspline, 
				 static_cast<size_t>(ninterp));
  RinterpFlux_pos = new double[ninterp];
  RinterpVals_pos = new double[ninterp];
  interpvars_allocated_pos = true;
  rinitialized = false; // Can't match rinterp any more
}


void PDFactory::freeInterpPos() {
  if (RinterpFlux_pos != NULL) { 
    delete[] RinterpFlux_pos; 
    RinterpFlux_pos = NULL; 
  }
  if (RinterpVals_pos != NULL) { 
    delete[] RinterpVals_pos; 
    RinterpVals_pos = NULL; 
  }
  interpvars_allocated_pos = false;
  rinitialized = false; // Can't match rinterp any more
  initialized = false;
}

void PDFactory::allocateInterpNeg() {
  if (interpvars_allocated_neg) return;
  if (ninterp == 0)
    throw affineExcept("PDFactory", "allocateInterpNeg",
		       "Invalid (0) ninterp");
  
  //Note that acc_neg is always allocated
  spline_neg = gsl_spline_alloc(gsl_interp_cspline, 
				static_cast<size_t>(ninterp));
  RinterpFlux_neg = new double[ninterp];
  RinterpVals_neg = new double[ninterp];
  interpvars_allocated_neg = true;
  rinitialized = false; // Can't match rinterp any more
}

void PDFactory::freeInterpNeg() {
  if (RinterpFlux_neg != NULL) { 
    delete[] RinterpFlux_neg; 
    RinterpFlux_neg = NULL; 
  }
  if (RinterpVals_neg != NULL) { 
    delete[] RinterpVals_neg; 
    RinterpVals_neg = NULL; 
  }
  interpvars_allocated_neg = false;
  initialized = false;
  rinitialized = false; // Can't match rinterp any more
}


/*!
  Frees all internal memory
*/
void PDFactory::free() {
  freeRvars();
  freeInterpPos();
  freeInterpNeg();
}

/*!
  \param[in] filename Name of wisdom file
*/
void PDFactory::addWisdom(const std::string& filename) {
  FILE *fp = NULL;
  fp = fopen(filename.c_str(), "r");
  if (fp == NULL) {
    std::stringstream str;
    str << "Error opening wisdom file: " << filename;
    throw affineExcept("PDFactory", "addWisdom", str.str());
  }
  if (fftw_import_wisdom_from_file(fp) == 0) {
    std::stringstream str;
    str << "Error reading wisdom file: " << filename;
    throw affineExcept("PDFactory", "addWisdom", str.str());
  }
  fclose(fp);
  fftw_plan_style = FFTW_WISDOM_ONLY;
  has_wisdom = true;
  wisdom_file = filename;
  if (plan != NULL) {
    fftw_destroy_plan(plan); 
    plan = NULL;
  }
  if (plan_inv != NULL) {
    fftw_destroy_plan(plan_inv);
    plan_inv = NULL;
  }

  initialized = false;
}

/*!
  \param[in] n New transform size.
  Sets up transforms, resizing if needed
*/
void PDFactory::setupTransforms(unsigned int n) {

  // The transform plan is that the forward transform
  //  takes rvals to rtrans.  We then do things to rtrans
  //  to turn it into pval, including shifting, adding noise,
  //  etc, and then transform from pval to pofd.
  // However, because we may move around the arrays,
  //  one should always call these plans in array-execute form

  if (n == 0)
    throw affineExcept("PDFactoryDouble", "setupTransforms",
		       "Invalid (0) transform size");

  // Make sure we have enough room
  resize(n);

  //We need the R variables allocated so we can plan to them
  if (!rvars_allocated) allocateRvars();

  int intn = static_cast<int>(n);
  if (plan == NULL) {
    if (plan != NULL) fftw_destroy_plan(plan);
    plan = fftw_plan_dft_r2c_1d(intn, rvals, rtrans,
				fftw_plan_style);
  }
  if (plan == NULL) {
    std::stringstream str;
    str << "Plan creation failed for forward transform of size: " << n;
    if (has_wisdom) str << std::endl << "Your wisdom file may not have"
			<< " that size";
    throw affineExcept("PDFactory", "setupTransforms", str.str());
  }

  if (plan_inv == NULL) {
    if (plan_inv != NULL) fftw_destroy_plan(plan_inv);
    plan_inv = fftw_plan_dft_c2r_1d(intn, pval, pofd,
				    fftw_plan_style);
  }
  if (plan_inv == NULL) {
    std::stringstream str;
    str << "Plan creation failed for inverse transform of size: " << n;
    if (has_wisdom) str << std::endl << "Your wisdom file may not have"
			<< " that size";
    throw affineExcept("PDFactory", "setupTransforms", str.str());
  }
}

/*!
  \param[in] n Number of elements
  \param[in] minflux Minimum flux to use in R
  \param[in] maxflux Maximum flux to use in R

  Sets up flux information for R.  minflux may not actually be
   achieved if it is negative, because we try to include flux=0 in that case.
   If maxflux is negative, it is treated as zero.  And if minflux is
   greater than zero -- it is treated as zero.

  We don't actually store an array of fluxes.  We just store the
   information to recreate it.  For purely positive fluxes, this is
   trivial -- the r flux array is i * dflux + minflux_R after this routine.
  
  If there are negative fluxes, this is more complicated.
  The positive fluxes go up to index wrapRidx.  So we have
   0 <= i <= wrapRidx  flux = i * dflux
   wrapRidx < i < n = (i - wrapRidx - 1) * dflux + minflux_R
  The latter doesn't quite get up to 0, since that was already included in
   the first part. 
*/
void PDFactory::initRFlux(unsigned int n, double minflux, double maxflux) {
  // Make sure there is room
  resize(n);

  if (n == 0)
    throw affineExcept("PDFactory", "initRFlux", "Invalid (0) n");
  if (maxflux < minflux) std::swap(minflux, maxflux);

  if (n == 1) {
    if (minflux < 0.0)
      throw affineExcept("PDFactory", "initRFlux", 
			 "n must be > 1 with negative fluxes");
    dflux = maxflux - minflux;
    minflux_R = minflux;
    wrapRidx = 0;
    return;
  }

  double inm1 = 1.0 / static_cast<double>(n - 1);
  if (maxflux > 0) {
    if (minflux > 0) minflux = 0.0; // Always need zero in there somehow
    if (minflux == 0) {
      // Simple case -- 0 to maxflux
      wrapRidx = n - 1;
      minflux_R = 0.0;
      dflux = maxflux * inm1;
    } else {
      // Most complicated case -- maxflux > 0, minflux < 0
      //  We would really like to have flux = 0 included.  Also, we 
      // are wrapping negative fluxes around to the top of the array.  
      // This doesn't work out for arbitrary min/max flux, 
      // so something has to give. The choice here is to tweak minflux 
      // by the smallest amount possible to make it work while keeping
      // maxflux the same.
      dflux = (maxflux - minflux) * inm1; // Initial value
      dflux = maxflux / (n - floor(-minflux / dflux) - 1.0); //Adjust
      
      // Figure out what index we go up to with positive fills
      wrapRidx = static_cast<unsigned int>(maxflux / dflux + 0.9999999999);
      minflux_R = (wrapRidx+1) * dflux - static_cast<double>(n) * dflux;
    }
  } else {
    // Maxflux is either 0 or negative.  We will therefore treat it
    // as zero, since we want to include it.
    if (minflux >= 0) // But then we need minflux to be negative
      throw affineExcept("PDFactory", "initRFlux", 
			 "No effective range in R fluxes");
    dflux = - minflux * inm1;
    wrapRidx = 0; // Includes 0, rest of the array is negative
    minflux_R = minflux;
  }
  rinitialized = false; // Any r likely no longer matches this
  initialized = false; // If there was a P(D), it's no longer valid
}

/*!
  \param[in] model number counts model to use for fill.  Params must be set
  \param[in] bm Beam 

  The R and interpolation variables must have already been allocated
  This fills in the interpolation variables with R, then interpolates
  out to the full range.

  Doesn't check validity of model or beam -- caller should do that  
*/
void PDFactory::initRInterp(const numberCounts& model, const beam& bm) {
  
  if (!rvars_allocated)
    throw affineExcept("PDFactory", "computeR",
		       "R variables must have been previously allocated");

  // Figure out which ones we need to fill
  bool has_pos = bm.hasPos();
  bool has_neg = bm.hasNeg();
  if (!(has_pos || has_neg))
    throw affineExcept("PDFactory", "computeR",
		       "Beam has neither positive nor negative bits");

  double modelmin, modelmax;
  modelmin = model.getMinFlux();
  modelmax = model.getMaxFlux();
  double ininterpm1 = 1.0 / static_cast<double>(ninterp-1);

  if (has_pos) {
    if (!interpvars_allocated_pos) allocateInterpPos();
    // Figure out limits -- we bascially just want to go out to
    //  where R is nonzero because we will be interpolating in log space.  
    //  That doesn't work on the lower edge (since R is nonzero for any
    //  value > 0), so we just use some constant times the lowest model flux.
    double mininterpflux = modelmin * subedgemult;
    double maxinterpflux = modelmax * bm.getMinMaxPos().second;

    // R[maxinterpflux] = 0.0, which we don't want to include in our
    //  log/log interpolation, so step it back slightly
    double dinterp = 0.001 * (maxinterpflux - mininterpflux) * ininterpm1;
    maxinterpflux -= dinterp;

    // Note that the interpolation is in log2 space, so therefore
    //  we want the interpolation flux to be spaced that way
    double dinterpflux = (log2(maxinterpflux / mininterpflux)) * ininterpm1;
    for (unsigned int i = 0; i < ninterp; ++i)
      RinterpFlux_pos[i] = mininterpflux * 
	exp2(static_cast<double>(i) * dinterpflux);
    
    // Now actually get R
#ifdef TIMING
    std::clock_t starttime = std::clock();
#endif
    model.getR(ninterp, RinterpFlux_pos, bm, RinterpVals_pos);
#ifdef TIMING
    RTime += std::clock() - starttime;
#endif
    
    //Load the values into the spline.  Note we take the log --
    // we are interpolating in log space 
    double val;
    for (unsigned int i = 0; i < ninterp; ++i) {
      val = RinterpVals_pos[i]; // We are paranoid about < 0 values
      if (val > 0) RinterpVals_pos[i] = log2(val); //log2 is faster than ln
      else RinterpVals_pos[i] = pofd_mcmc::smalllogval;
    }
    gsl_spline_init(spline_pos, RinterpFlux_pos, RinterpVals_pos, 
		    static_cast<size_t>(ninterp));

  }
  // Same, but for negative beam
  if (has_neg) {
    if (!interpvars_allocated_neg) allocateInterpNeg();

    double maxinterpflux = -modelmin * subedgemult;
    double mininterpflux = -modelmax * bm.getMinMaxNeg().second;

    // Step back from both edges slightly
    double dinterp = 0.001 * (maxinterpflux - mininterpflux) * ininterpm1;
    mininterpflux += dinterp; 
    maxinterpflux -= dinterp;

    // Note that the interpolation is in log2 space, so therefore
    //  we want the interpolation flux to be spaced that way
    // Except we want it spaced with more values at the top end now
    // Also recall that the GSL spline class must have values in increasing
    //  order
    double dinterpflux = (log2(mininterpflux / maxinterpflux)) * ininterpm1;
    for (unsigned int i = 0; i < ninterp; ++i)
      RinterpFlux_neg[i] = maxinterpflux * 
	exp2(static_cast<double>(ninterp - i - 1) * dinterpflux);
    
#ifdef TIMING
    std::clock_t starttime = std::clock();
#endif
    model.getR(ninterp, RinterpFlux_neg, bm, RinterpVals_neg);
#ifdef TIMING
    RTime += std::clock() - starttime;
#endif

    double val;
    for (unsigned int i = 0; i < ninterp; ++i) {
      val = RinterpVals_neg[i];
      if (val > 0) RinterpVals_neg[i] = log2(val); //log2 is faster than ln
      else RinterpVals_neg[i] = pofd_mcmc::smalllogval;
    }
    gsl_spline_init(spline_neg, RinterpFlux_neg, RinterpVals_neg, 
		    static_cast<size_t>(ninterp));
  }
  rinitialized = false; // Any r likely no longer matches this
  initialized = false; //!< any transform no longer matches Rinterp
}
  
/*!
  \param[in] n New size of transform.  Will resize if different than old
  \param[in] minflux Minimum flux desired in R fill
  \param[in] maxflux Maximum flux desired in R fill
  \param[in] model Number counts model
  \param[in] bm Beam
  \param[in] muldr Multiply R by dflux

  This is the interface to initRInterp and initRFlux, plus does the
  filling of the real R.  It does not check the model or beam for
  validity, so the caller (initPD) must do so.
*/
void PDFactory::initR(unsigned int n, double minflux, double maxflux,
		      const numberCounts& model, const beam& bm,
		      bool muldr) {

  // Set up interp vars, flux ranges
  initRInterp(model, bm);
  initRFlux(n, minflux, maxflux);

  // Now interpolate out R into rvals.

  //Now interpolate out; note r still needs to be multiplied by dflux
  // for use later.  We figure out which bins are covered by the
  // interpolation for efficiency
  bool haspos = bm.hasPos() && (maxflux > 0);
  bool hasneg = bm.hasNeg() && (minflux < 0) && (wrapRidx > 0);
  
  std::memset(rvals, 0, n * sizeof(double));
  if (haspos) {
    double minI = RinterpFlux_pos[0]; // Min interpolated value
    double maxI = RinterpFlux_pos[ninterp-1]; // Max interpolated value
    
    // This is the minimum R index which will be non-zero for the
    // positive bit
    unsigned int minitidx = static_cast<unsigned int>(minI / dflux + 
						      0.9999999999999999);
    // And this the maximum
    unsigned int maxitidx = static_cast<unsigned int>(maxI / dflux);
    if (maxitidx > wrapRidx) maxitidx = wrapRidx;

#ifdef TIMING
    std::clock_t starttime = std::clock();
#endif
    double cflux, splval;
    for (unsigned int i = minitidx; i <= maxitidx; ++i) {
      cflux = static_cast<double>(i) * dflux; //Min is always 0 in positive R
      splval = gsl_spline_eval(spline_pos, cflux, acc_pos);
      rvals[i] = exp2(splval); // Recall -- spline is log2
    }
#ifdef TIMING
    RTime += std::clock() - starttime;
#endif

  }
  if (hasneg) {
    double minI = RinterpFlux_neg[0]; // Min interpolated value
    double maxI = RinterpFlux_neg[ninterp-1]; // Max interpolated value
    
    // This is the minimum R index which will be non-zero for the
    // negative bit.  The computation is more complicated here
    unsigned int minitidx = 
      static_cast<unsigned int>((minI - minflux_R) / dflux + 
				wrapRidx + 1 + 0.9999999999999999);
    if (minitidx <= wrapRidx) minitidx = wrapRidx + 1; // Necessary?
    unsigned int maxitidx = 
      static_cast<unsigned int>((maxI - minflux_R) / dflux +
				wrapRidx + 1);
    if (maxitidx > n - 1) maxitidx = n - 1;

#ifdef TIMING
    std::clock_t starttime = std::clock();
#endif
    double cflux, splval;
    double mval = minflux_R - static_cast<double>(wrapRidx + 1)*dflux;
    for (unsigned int i = minitidx; i <= maxitidx; ++i) {
      cflux = static_cast<double>(i)* dflux + mval;
      splval = gsl_spline_eval(spline_neg, cflux, acc_neg);
      rvals[i] = exp2(splval);
    }
#ifdef TIMING
    RTime += std::clock() - starttime;
#endif
  }

  if (muldr) {
    for (unsigned int i = 0; i < n; ++i) rvals[i] *= dflux;
    rdflux = true;
  } else rdflux = false;

  rinitialized = true;
}

/*!
  The mean is stored in mn, the variance in var_noi
*/
void PDFactory::getMeanVarFromR() {

  if (!rinitialized)
    throw affineExcept("PDFactory", "getMeanVarFromR",
		       "R variables not computed");

  // We take advantage of the relations:
  //  mn = \int x R dx
  //  var = \int x^2 R dx
  // and use the trapezoid rule. The fact that we may have wrapped
  // makes this a bit complicated...  we deal with this by doing a
  // straigt sum (that is, ignore the 0.5 at the edges) then fixing
  // for that.  Well, actually, only the lower edge needs special treatment
  double cf, prod;
  mn = 0.0;
  var_noi = 0.0;
  // Pos fluxes
  // Can start at 1 because cf = 0 at i = 0, so no contribution
  for (unsigned int i = 1; i <= wrapRidx-1; ++i) {
    cf = static_cast<double>(i) * dflux;
    prod = cf * rvals[i];
    mn += prod;
    var_noi += cf * prod;
  }
  // Upper edge
  cf = static_cast<double>(wrapRidx) * dflux;
  prod = 0.5 * cf * rvals[wrapRidx];
  mn += prod;
  var_noi += cf * prod;
  // Neg fluxes
  for (unsigned int i = wrapRidx + 1; i < currsize; ++i) {
    cf = static_cast<double>(i - wrapRidx - 1) * dflux + minflux_R;
    prod = cf * rvals[i];
    mn += prod;
    var_noi += cf * prod;
  }
  // Correct lower limit
  if (wrapRidx < currsize - 1) { // Have negative fluxes, index is wrapRidx+1
    cf = minflux_R;
    prod = 0.5 * cf * rvals[wrapRidx - 1];
    mn -= prod;
    var_noi -= cf * prod;
  } else { // No neg fluxes - lower limit is first element.  But these are 0
    //cf = 0;
    //prod = 0.5 * cf * rvals[0]; 
    //mn -= prod;
    //var_noi -= cf * prod;
  }
  
  // Deal with presence/absence of dflux in R
  if (!rdflux) {
    mn *= dflux;
    var_noi *= dflux;
  }
}

/*!
  \returns Split point index

  This is a utility function for finding the split point to 'unwrap' a P(D).
  The details are very slightly different for 1D and 2D because of the
  way that flux densities are stored, so it is implemented
  both here and in PDFactoryDouble.

  Uses a bunch of internal variables that are assumed set.
*/
unsigned int PDFactory::findSplitPoint() const {
  // The logic here is a bit complex.  This code has caused a lot
  //  of issues in fits, so it's worth explaining in some detail.
  //
  // What we have is a P(D) (in pofd) that is set up so that the peak 
  //  (really, the mean, but in practice this is pretty close to the peak) 
  //  is near index zero.  Because we compute the P(D) with an FFT, the P(D) 
  //  wraps around the end of the array, so that the top part of the array
  //  is actually the P(D) at negative flux densities, but the P(D) of 
  //  positive flux densities is also coming up the array from 0.  The issue
  //  is that these overlap slightly, and we want to find a good place to
  //  split the array into negative and positive flux densities.
  // Here's what we are hoping for:
  //  1) We want the split point to be more than a few (n1) sigma away from
  //      the mean value in both negative and positive flux density.
  //  2) We want the value of the data array at the split to be less
  //      than some small fraction of the peak value.
  //  3) We want the split point to be modestly stable in that running
  //      slightly different models, or on different machines should
  //      not produce wildly different split points.
  // The third point isn't entirely necessary for the code to work, but
  //  it makes debugging much, much easier.  The problem is that adding
  //  that third requirement makes things rather more difficult.  Without
  //  that we just split on the smallest value that is far enough away.
  //  But that isn't numerically stable because we frequently have really
  //  tiny values in pd that are just numeric noise.
  // The steps used here to achieve these goals is:
  //  1) Make sure that the total flux density range is larger than 
  //      n1 * sigma.  If it isn't, we can never satisfy point 1,
  //      so throw an exception.
  //  2) Find the minimum location and value.
  //  3) Find the maximum (peak) value.
  //  4) If that location satisfies:
  //       a) value > f1 * peak
  //       b) value <= f2 * peak
  //       c) more than n1 * sigma away from edges
  //     then accept that as the split point.  b and c are
  //     just requirements 1 and 2, but a is about requirement 3.
  //     If the min is < f1 * peak (where f1 is some small number),
  //     then the min point is just noise, so we have to do something
  //     else.
  //  5) Now start at n1 * sigma from the top and go down to n1 * sigma
  //      from the bottom.  Return the first point encountered that
  //      is less than f2 * peak.
  //  6) If no point is found, throw an exception.
  // Here we are taking advantage of the assumption that the P(D) will
  //  have much more power in its positive tail than its negative tail.
  //  This is not true in general, but is true for the data we are
  //  planning on using this code on because the beam is mostly positive.
  //  For data sets where the beam is not like that, this code will
  //  be non-optimal.

  // These are are control constants for accepting points
  const double n1 = 3.0; // Minimum number of sigma away from mean
  const double f1 = 1e-10; // Values smaller than f1 * peak considered noise
  const double f2 = 1e-5; // Split must be less than f2 * peak

  if (sg <= 0)
    throw affineExcept("PDFactory", "findSplitPoint", 
		       "Invalid (non-positive) sigma");

  double fluxrange = currsize * dflux;
  double n1sig = n1 * sg;

  // Step 1
  if (fluxrange < 2 * n1sig) {
    std::stringstream errstr;
    errstr << "Not enough flux range; range is: "
	   << fluxrange << " sigma is: " << sg
	   << " so required flux range is: "
	   << 2 * n1sig;
    throw affineExcept("PDFactory", "findSplitPoint", 
		       errstr.str());
  }

  // Step 2 and 3
  unsigned int splitidx;
  double peak, minval, currval;
  splitidx = 0;
  peak = minval = pofd[0];
  for (unsigned int i = 1; i < currsize; ++i) {
    currval = pofd[i];
    if (currval > peak) peak = currval;
    if (currval < minval) {
      splitidx = i;
      minval = currval;
    }
  }

  // Step 4
  // First we check minval <= f2 * peak.  If this
  //  isn't true at the minimum, there is no way
  //  to split this.
  double targval = f2 * peak;
  if (minval > targval) {
    std::stringstream errstr;
    errstr << "Minimum in pd is too large relative to peak;"
	   << " can't find split point with minimum: "
	   << minval << " peak value: " << peak
	   << " Required split value < " << targval;
    throw affineExcept("PDFactory", "findSplitPoint", 
		       errstr.str());
  }
  // Now see if this point satisfies our criteria 1 and 2
  //  These are the flux densities of the proposed splitidx
  //  in postive and negative space (abs of the latter)
  double fwrap_pos = static_cast<double>(splitidx) * dflux;
  double fwrap_neg = 
    fabs(static_cast<double>(splitidx - wrapRidx - 1) * dflux +
	 minflux_R);
  if ((minval > f1 * peak) && (fwrap_pos >= n1sig) && (fwrap_neg >= n1sig))
    return splitidx;
  
  // Step 5
  // Figure out the search index range from -n1sig to +n1sig
  unsigned int topidx = 
    static_cast<unsigned int>((n1sig - minflux_R) / dflux) + wrapRidx;
  unsigned int botidx = static_cast<unsigned int>(n1sig / dflux);
  for (unsigned int i = topidx; i > botidx; --i) {
    currval = pofd[i];
    if (currval <= targval) return i; // Success!
  }

  // Step 6
  // Didn't find one -- throw exception
  std::stringstream errstr;
  errstr << "Unable to find split point in range "
	 << topidx << " down to " << botidx
	 << " with value less than " << targval;
  throw affineExcept("PDFactory", "findSplitPoint", 
		     errstr.str());
}



/*!
  \param[out] pd Holds P(D) on output, normalized, mean subtracted,
                 and with positivity enforced.
  
  This does the unwrapping of the internal P(D) into pd,
  as well as removing negative values, normalizing, fixing
  up the flux scale.
*/
// This should only ever be called by getPD, so the inputs aren't checked
//  Among other implications, we are assuming the P(D) is in pofd,
//  the sigma is in sg, etc.
void PDFactory::unwrapAndNormalizePD(PD& pd) const {
  // The idea is to find the point where the postive and negative bits
  //  cross each other and split there.  The details of that
  //  computation are in PDFactory::findSplitPoint.

  // First, Enforce positivity
#ifdef TIMING
  starttime = std::clock();
#endif
  for (unsigned int idx = 0; idx < currsize; ++idx)
    if (pofd[idx] < 0) pofd[idx] = 0.0;
#ifdef TIMING
  posTime += std::clock() - starttime;
#endif

  // Find the point we want to split at
  unsigned int splitidx = findSplitPoint();
  
  // Copy over.  Things above splitidx in pofd go into the bottom of
  // pd.pd_, then the stuff below that in pofd goes above that in pd.pd_
  // in the same order.
  pd.resize(currsize);
  double *ptr_curr, *ptr_out; // convenience vars
  ptr_curr = pofd + splitidx;
  ptr_out = pd.pd_;
  //for (unsigned int i = 0; i < currsize - splitidx; ++i)
  //  ptr_out[i] = ptr_curr[i];
  std::memcpy(ptr_out, ptr_curr, (currsize - splitidx) * sizeof(double));
  ptr_curr = pofd;
  ptr_out = pd.pd_ + currsize - splitidx;
  //for (unsigned int i = 0; i < splitidx; ++i)
  //  ptr_out[i] = ptr_curr[i];
  std::memcpy(ptr_out, ptr_curr, splitidx * sizeof(double));

  pd.logflat = false;
  pd.minflux = 0.0; pd.dflux = dflux;

#ifdef TIMING
  copyTime += std::clock() - starttime;
#endif

  //Normalize
#ifdef TIMING
  starttime = std::clock();
#endif
  pd.normalize();
#ifdef TIMING
  normTime += std::clock() - starttime;
#endif

  //Now mean subtract flux axis
#ifdef TIMING
  starttime = std::clock();
#endif
  double tmn; //True mean
  tmn = pd.getMean(false);
  if (std::isinf(tmn) || std::isnan(tmn)) {
    std::stringstream str;
    str << "Un-shift amounts not finite: " << tmn << " " << std::endl;
    str << "At length: " << currsize << " with noise: " << sg;
    throw affineExcept("PDFactory", "unwrapandNormalizePD", str.str());
  }
  pd.minflux = -tmn;
#ifdef TIMING
  meanTime += std::clock() - starttime;
#endif
}

/*!
  \param[in] n       Size of transform 
  \param[in] minflux Minimum flux generated in R
  \param[in] maxflux Maximum flux generated in R
  \param[in] model   number counts model to use for fill.  Params must be set
  \param[in] bm      Beam 
  \returns True if the P(D) could be initialized, false if something
           about the parameters prevented initialization.  Note that
	   a genuine error results in throwing an exception, not setting this
	   to false.

  Gets ready for P(D) computation by preparing R, including forward
  transforming it.  The idea is to compute everything that doesn't
  require knowing sigma here so we can call this multiple times on
  maps with the same beam but different noise values.
*/
bool PDFactory::initPD(unsigned int n, double minflux, double maxflux, 
		       const numberCounts& model, const beam& bm) {

  if (n == 0)
    throw affineExcept("PDFactory", "initPD", "Invalid (zero) n");  
  if (n == 1)
    throw affineExcept("PDFactory", "initPD", "Invalid (==0) n");  
  if (maxflux <= 0.0) {
    std::stringstream errstr;
    errstr << "Invalid (non-positive) maxflux (" << maxflux << ")";
    throw affineExcept("PDFactory", "initPD", errstr.str());
  }
  if (minflux == maxflux) {
    std::stringstream errstr;
    errstr << "Invalid (0) range between min/max flux with value "
	   << minflux;
    throw affineExcept("PDFactory", "initPD", errstr.str());
  }
  if (!model.isValid())
    throw affineExcept("PDFactory", "initPD", "Invalid model");
  if (!bm.hasData())
    throw affineExcept("PDFactory", "initPD", "Empty beam");

  // We are about to nuke these
  rinitialized = false;
  initialized = false;

  //Make the plans, or keep the old ones if possible
  // Note we have to do this before we fill R, as plan construction
  // may overwrite the values.  This allocates R and sets currfftlen
  setupTransforms(n);

  // Initialize R, multiplying it by dflux
  initR(n, minflux, maxflux, model, bm, true);

  //Now that we have R, use it to compute the mean and variance
  // (stored in mn and var_noi)
  getMeanVarFromR();

  // Set this to the value without instrument noise for now
  sg = sqrt(var_noi);

  //Decide if we will shift and pad, and if so by how much.
  // The idea is to shift the mean to zero -- but we only
  // do the shift if the sigma is larger than one actual step size
  // because otherwise we can't represent it well.
  doshift = (sg > dflux) && (fabs(mn) > dflux);
  if (doshift) shift = -mn; else shift = 0.0;

  if (verbose) {
    std::cout << " Initial mean estimate: " << mn << std::endl;
    std::cout << " Initial stdev estimate (no inst noise): " << sg << std::endl;
    if (doshift)
      std::cout << " Additional shift applied: " << shift << std::endl;
    else 
      std::cout << " Not applying additional shift" << std::endl;
  }

  //Compute forward transform of this r value, store in rtrans
  //Have to use argument version, since the address of rtrans, etc. can move
#ifdef TIMING
  starttime = std::clock();
#endif
  fftw_execute_dft_r2c(plan, rvals, rtrans); 
#ifdef TIMING
  fftTime += std::clock() - starttime;
#endif

  initialized = true;

  return true;
}

/*!
  Computes the P(D) for a specific value of the instrument noise,
  assuming that initPD has been called first.
 
  \param[in] sigma Instrument noise sigma
  \param[out] pd Holds P(D) on output
  \param[in] setLog If true, pd is log(P(D) on output; convenient
              for likelihood evaluation.

  You must call initPD first.
*/
void PDFactory::getPD(double sigma, PD& pd, bool setLog) {

  // The basic idea is to compute the P(D) from the previously filled
  // R values, adding in noise and all that fun stuff, filling pd
  // for output

  if (!initialized)
    throw affineExcept("PDFactory", "getPD", "Must call initPD first");

  //Output array from 2D FFT is n/2+1
  unsigned int n = currsize;
  unsigned int ncplx = n / 2 + 1;
      
  //Calculate p(omega) = exp( r(omega) - r(0) ),
  // which is what we will transform back into pofd.
  // The forward transform of R is stored in rtrans 
  // There are some complications because of shifts and all that.
  // The frequencies are:
  //  f = i/dflux*n  
  // We work in w instead of f (2 pi f)
  // and actually compute 
  //  exp( r(omega) - r(0) - i*shift*omega - 1/2 sigma^2 omega^2 )
  if (verbose) std::cout << "  Computing p(w)" << std::endl;

#ifdef TIMING
  std::clock_t starttime = std::clock();
#endif

  double r0, expfac, rval, ival;
  r0 = rtrans[0][0]; //r[0] is pure real
  double iflux = mcmc_affine::two_pi / (n * dflux);

  // Four cases -- with and without shifting, with and without sigma.
  //  Write them out explicitly since this is inner loop stuff
  if (doshift) {
    double w;
    if (sigma > 0) {
      double sigfac = 0.5 * sigma * sigma;
      for (unsigned int idx = 1; idx < ncplx; ++idx) {
	w = iflux * static_cast<double>(idx);
	rval = rtrans[idx][0] - r0 - sigfac * w * w;
	ival = rtrans[idx][1] - shift * w;
	expfac = exp(rval);
	pval[idx][0] = expfac * cos(ival);
	pval[idx][1] = expfac * sin(ival);
      } 
    } else {
      for (unsigned int idx = 1; idx < ncplx; ++idx) {
	w = iflux * static_cast<double>(idx);
	rval = rtrans[idx][0] - r0;
	ival = rtrans[idx][1] - shift * w;
	expfac = exp(rval);
	pval[idx][0] = expfac * cos(ival);
	pval[idx][1] = expfac * sin(ival);
      } 
    }
  } else {
    double expfac, ival;
    if (sigma > 0.0) {
      double w, rval;
      double sigfac = 0.5 * sigma * sigma;
      for (unsigned int idx = 1; idx < ncplx; ++idx) {
	w = iflux * static_cast<double>(idx);
	rval = rtrans[idx][0] - r0 - sigfac * w * w;
	ival = rtrans[idx][1];
	expfac = exp(rval);
	pval[idx][0] = expfac * cos(ival);
	pval[idx][1] = expfac * sin(ival);
      }
    } else {
      for (unsigned int idx = 1; idx < ncplx; ++idx) {
	expfac = exp(rtrans[idx][0]-r0);
	ival = rtrans[idx][1];
	pval[idx][0] = expfac * cos(ival);
	pval[idx][1] = expfac * sin(ival);
      }
    }
  }
  //p(0) is special
  pval[0][0] = 1.0;
  pval[0][1] = 0.0;

#ifdef TIMING
  p0Time += std::clock() - starttime;
#endif

  //Transform back
  if (verbose) std::cout << " Reverse transform" << std::endl;
#ifdef TIMING
  starttime = std::clock();
#endif
  // From pval into pofd
  fftw_execute_dft_c2r(plan_inv, pval, pofd);
#ifdef TIMING
  fftTime += std::clock() - starttime;
#endif

  // Copy into output variable, also normalizing, mean subtracting, 
  // making positive
  sg = sqrt(var_noi + sigma * sigma); // Used in unwrapPD
  unwrapAndNormalizePD(pd);

  //Turn PD to log for more efficient log computation of likelihood
#ifdef TIMING
  starttime = std::clock();
#endif
  if (setLog) pd.applyLog(false);
#ifdef TIMING
  logTime += std::clock() - starttime;
#endif
}

/*!
  \param[in] filename File to write to

  You must call initPD or initR first, or bad things will probably happen.
*/
void PDFactory::writeRToHDF5(const std::string& filename) const {
  if (!rinitialized )
    throw affineExcept("PDFactory", "writeRToHDF5",
		       "Must call initR or initPD first");

  hid_t file_id;
  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
		      H5P_DEFAULT);

  if (H5Iget_ref(file_id) < 0) {
    H5Fclose(file_id);
    throw affineExcept("PDFactory", "writeToHDF5",
		     "Failed to open HDF5 file to write");
  }

  // Write it as one dataset -- Rflux, R. 
  hsize_t adims;
  hid_t mems_id, att_id, dat_id;
  
  // Properties
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate2(file_id, "dflux", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dflux);
  H5Aclose(att_id);
  H5Sclose(mems_id);

  // Rflux.  Set up a temporary array to put RFlux in
  double* RFlux = new double[currsize];
  for (unsigned int i = 0; i <= wrapRidx; ++i)
    RFlux[i] = static_cast<double>(i) * dflux;
  double mval = minflux_R - static_cast<double>(wrapRidx + 1)*dflux;
  for (unsigned int i = wrapRidx + 1; i < currsize; ++i)
    RFlux[i] = static_cast<double>(i) * dflux - mval;
  adims = currsize;
  mems_id = H5Screate_simple(1, &adims, NULL);
  dat_id = H5Dcreate2(file_id, "RFlux", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
	   H5P_DEFAULT, RFlux);
  H5Dclose(dat_id);
  delete[] RFlux;

  // R -- which we may need to copy to remove the dflux
  dat_id = H5Dcreate2(file_id, "R", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (rdflux) {
    double* tmp = new double[currsize];
    double idflux = 1.0 / dflux;
    for (unsigned int i = 0; i < currsize; ++i) tmp[i] = rvals[i] * idflux;
    H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,  H5P_DEFAULT, tmp);
    delete[] tmp;
  } else
    H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,  H5P_DEFAULT, rvals);
  H5Dclose(dat_id);
  H5Sclose(mems_id);

  // Done
  H5Fclose(file_id);
}


void PDFactory::sendSelf(MPI_Comm comm, int dest) const {
  MPI_Send(const_cast<unsigned int*>(&fftw_plan_style), 1, MPI_UNSIGNED,
	   dest, pofd_mcmc::PDFSENDPLANSTYLE, comm);
  MPI_Send(const_cast<bool*>(&has_wisdom), 1, MPI::BOOL, dest, 
	   pofd_mcmc::PDFHASWISDOM, comm);
  if (has_wisdom) {
    //Send wisdom file name
    unsigned int nstr = wisdom_file.size()+1;
    char *cstr = new char[nstr];
    std::strncpy(cstr, wisdom_file.c_str(), nstr);
    MPI_Send(&nstr, 1, MPI_UNSIGNED, dest, pofd_mcmc::PDFWISLEN, comm);
    MPI_Send(cstr, nstr, MPI_CHAR, dest, pofd_mcmc::PDFWISNAME, comm);
    delete[] cstr;
  }
  MPI_Send(const_cast<bool*>(&verbose), 1, MPI::BOOL, dest,
	   pofd_mcmc::PDFVERBOSE, comm);
  MPI_Send(const_cast<unsigned int*>(&ninterp), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::PDFNINTERP, comm);
}

//Note this doesn't copy over interal variables
void PDFactory::receiveCopy(MPI_Comm comm, int src) {
  MPI_Status Info;
  MPI_Recv(&fftw_plan_style, 1, MPI_UNSIGNED, src, 
	   pofd_mcmc::PDFSENDPLANSTYLE, comm, &Info);
  MPI_Recv(&has_wisdom, 1, MPI::BOOL, src, pofd_mcmc::PDFHASWISDOM,
	   comm, &Info);
  if (has_wisdom) {
    //Receive wisdom file name
    unsigned int nstr;
    MPI_Recv(&nstr, 1, MPI_UNSIGNED, src, pofd_mcmc::PDFWISLEN, comm, &Info);
    char *cstr = new char[nstr];
    MPI_Recv(cstr, nstr, MPI_CHAR, src, pofd_mcmc::PDFWISNAME, comm, &Info);    
    wisdom_file = std::string(cstr);
    delete[] cstr;
    addWisdom(wisdom_file);
  }
  MPI_Recv(&verbose, 1, MPI::BOOL, src, pofd_mcmc::PDFVERBOSE, comm, &Info);
  unsigned int nnewinterp;
  MPI_Recv(&nnewinterp, 1, MPI_UNSIGNED, src, pofd_mcmc::PDFNINTERP,
	   comm, &Info);
  if (nnewinterp != ninterp) {
    freeInterpPos();
    freeInterpNeg();
    ninterp = nnewinterp;
  }
  rinitialized = false;
  initialized = false;
}
