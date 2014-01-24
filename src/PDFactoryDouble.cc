#include<sstream>
#include<cmath>
#include<cstring>
#include<limits>
#include<iomanip>

#include "../include/global_settings.h"
#include "../include/PDFactoryDouble.h"
#include "../include/affineExcept.h"

const double PDFactoryDouble::lowEdgeRMult=1e-9; //!< How far to go down on edge
//Control of how we do the edge integrals -- linear or log?
const bool PDFactoryDouble::use_edge_log_x = true; //!< edge log x spacing
const bool PDFactoryDouble::use_edge_log_y = false; //!< edge log y spacing


// Helper function for RFlux setup.  After this is called
//  rflux is filled with values from minflux_realized to maxflux 
//  in steps of dflux with wrapping to negative values occuring at wrapidx.
// That is, for all indices [0, wrapidx] the flux densities are positive,
//  above that they are negative (which may be the empty set).
void initRFluxInternal(unsigned int n, double minflux, double maxflux,
		       double* const rflux, double& dflux, 
		       double& minflux_realized, unsigned int& wrapidx) {
  // Assume n > 1 (checked by initRFlux)
  double inm1 = 1.0 / static_cast<double>(n - 1);
  dflux = (maxflux - minflux) * inm1;
  if (minflux >= 0.0) {
    // Easy case
    rflux[0] = minflux;
    for (unsigned int i = 1; i < n; ++i)
      rflux[i] = static_cast<double>(i) * dflux + minflux;
    minflux_realized = minflux;
    wrapidx = n - 1;
  } else {
    // Here the complication is that we would really like to have
    // RFlux = 0 included in the array.  Also, we are wrapping 
    // negative fluxes around to the top of the array.
    // We do this by tweaking minflux slightly
    dflux = maxflux / (n - floor(-minflux / dflux) - 2.0);
    // Figure out what index we go up to with positive fills
    wrapidx = static_cast<unsigned int>(maxflux / dflux);
    rflux[0] = 0.0;
    for (unsigned int i = 1; i <= wrapidx; ++i) // Pos Rflux
      rflux[i] = static_cast<double>(i) * dflux;

    // Figure out new minimum flux
    double wrapval = - static_cast<double>(n) * dflux;
    for (unsigned int i = wrapidx + 1; i < n; ++i) // Wrapped neg Rflux
      rflux[i] = static_cast<double>(i) * dflux + wrapval;
    minflux_realized = rflux[wrapidx + 1];
  }
}		       

/*!
  \param[in] nedge Size of edge integrals
*/
PDFactoryDouble::PDFactoryDouble(unsigned int nedge) {
  init(nedge);
}

/*
  \param[in] wisfile Wisdom file filename
  \param[in] nedge Size of edge integrals
*/
PDFactoryDouble::PDFactoryDouble(const std::string& wisfile,
				 unsigned int nedge) {
  init(nedge);
  addWisdom(wisfile);
}

PDFactoryDouble::~PDFactoryDouble() {
  if (RFlux1 != NULL) fftw_free(RFlux1);
  if (RFlux2 != NULL) fftw_free(RFlux2);

  if (rvals != NULL)  fftw_free(rvals);
  if (rtrans != NULL) fftw_free(rtrans);
  if (pofd != NULL)   fftw_free(pofd);
  if (pval != NULL)   fftw_free(pval);

  if (REdgeFlux1 != NULL) fftw_free(REdgeFlux1);
  if (REdgeFlux2 != NULL) fftw_free(REdgeFlux2);
  if (REdgeWork != NULL)  fftw_free(REdgeWork);

  if (plan != NULL) fftw_destroy_plan(plan); 
  if (plan_inv != NULL) fftw_destroy_plan(plan_inv);
}

/*!
  \param[in] NEDGE Number of edge integral steps
*/
void PDFactoryDouble::init(unsigned int NEDGE) {
  currsize = 0;

#ifdef TIMING
  resetTime();
#endif

  rvars_allocated = false;
  RFlux1 = RFlux2 = NULL;
  RwrapIdx1 = RwrapIdx2 = 0;
  rdflux = false;
  rvals = NULL;
  rtrans = NULL;
  pofd = NULL;
  pval = NULL;

  nedge = NEDGE;
  edgevars_allocated = false;
  REdgeFlux1 = NULL;
  REdgeFlux2 = NULL;
  REdgeWork  = NULL;

  dflux1 = dflux2 = 0.0;

  plan = plan_inv = NULL;

  verbose = false;
  has_wisdom = false;
  fftw_plan_style = FFTW_MEASURE;

  mn1 = mn2 = var_noi1 = var_noi2 = sg1 = sg2 = 
    std::numeric_limits<double>::quiet_NaN();
  rinitialized = false;
  initialized = false;
}

/*!
  \param[in] nedg Edge integration size
*/
void PDFactoryDouble::setNEdge(unsigned int nedg) {
  if (nedg == nedge) return;
  nedge = nedg;
  if (edgevars_allocated) {
    freeEdgevars();
    allocateEdgevars();
  }
}


#ifdef TIMING
void PDFactoryDouble::resetTime() {
  RTime = p0Time = fftTime = posTime = copyTime = normTime = 0;
  edgeTime = meanTime = logTime = 0;
}

/*!
  \param[in] nindent Number of indentation spaces
*/
void PDFactoryDouble::summarizeTime(unsigned int nindent) const {
  std::string prestring(nindent,' ');
    
  std::cout << "R time: " << prestring 
	    << 1.0*RTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "p0 time: " << prestring 
	    << 1.0*p0Time/CLOCKS_PER_SEC << std::endl;
  std::cout << "fft time: " << prestring 
	    << 1.0*fftTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "pos time: " << prestring 
	    << 1.0*posTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "copy time: " << prestring 
	    << 1.0*copyTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "norm time: " << prestring 
	    << 1.0*normTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "edge time: " << prestring 
	    << 1.0*edgeTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "mean time: " << prestring 
	    << 1.0*meanTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "log time: " << prestring 
	    << 1.0*logTime/CLOCKS_PER_SEC << std::endl;
}
#endif

/*
  \param[in] NSIZE new size
*/
void PDFactoryDouble::resize(unsigned int NSIZE) {
  if (NSIZE == currsize) return;
  freeRvars();
  currsize = NSIZE;
  if (plan != NULL) { fftw_destroy_plan(plan); plan = NULL; }
  if (plan_inv != NULL) { fftw_destroy_plan(plan_inv); plan_inv = NULL; }  
  rinitialized = false;
  initialized = false;
}

// Note this doesn't necessarily invalidate the plans, as long as
// they are always called in array execute mode
void PDFactoryDouble::allocateRvars() {
  if (rvars_allocated) return;
  if (currsize == 0)
    throw affineExcept("PDFactoryDouble", "allocate_rvars",
		       "Invalid (0) currsize");
  RFlux1 = (double*) fftw_malloc(sizeof(double) * currsize);
  RFlux2 = (double*) fftw_malloc(sizeof(double) * currsize);
  unsigned int fsize = currsize * currsize;
  rvals = (double*) fftw_malloc(sizeof(double) * fsize);
  pofd  = (double*) fftw_malloc(sizeof(double) * fsize);
  fsize = currsize * (currsize / 2 + 1);
  rtrans = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fsize);
  pval = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fsize);

  rvars_allocated = true;
  rinitialized = false;
  initialized = false;
}

// Also doesn't necessarily invalidate the plans
void PDFactoryDouble::freeRvars() {
  if (RFlux1 != NULL)     { fftw_free(RFlux1); RFlux1=NULL; }
  if (RFlux2 != NULL)     { fftw_free(RFlux2); RFlux2=NULL; }
  if (rvals != NULL)      { fftw_free(rvals); rvals=NULL; }
  if (rtrans != NULL)     { fftw_free(rtrans); rtrans=NULL; }
  if (pval != NULL)       { fftw_free(pval); pval = NULL; }
  if (pofd != NULL)       { fftw_free(pofd); pofd = NULL; }
  rvars_allocated = false;
  rinitialized = false;
  initialized = false;
}

// Note this doesn't affect the initialization state for the 
// R forward transform because the stuff done in the edge vars
// is copied over before the transform
void PDFactoryDouble::allocateEdgevars() {
  if (edgevars_allocated) return;
  if (nedge > 0) {
    REdgeFlux1 = (double*) fftw_malloc(sizeof(double)*nedge);
    REdgeFlux2 = (double*) fftw_malloc(sizeof(double)*nedge);
    REdgeWork  = (double*) fftw_malloc(sizeof(double)*nedge*nedge);
    edgevars_allocated = true;
  } else {
    REdgeWork = REdgeFlux1 = REdgeFlux2 = NULL;
    edgevars_allocated = false;
  }
}

void PDFactoryDouble::freeEdgevars() {
  if (REdgeFlux1 != NULL) { fftw_free(REdgeFlux1); REdgeFlux1 = NULL; }
  if (REdgeFlux2 != NULL) { fftw_free(REdgeFlux2); REdgeFlux2 = NULL; }
  if (REdgeWork != NULL)  { fftw_free(REdgeWork); REdgeWork = NULL; }
  edgevars_allocated = false;
}

/*!
  Frees all internal memory
*/
void PDFactoryDouble::free() {
  freeRvars();
  freeEdgevars();
}

/*!
  \param[in] filename Name of wisdom file
*/
void PDFactoryDouble::addWisdom(const std::string& filename) {
  FILE *fp = NULL;
  fp = fopen(filename.c_str(), "r");
  if (fp == NULL) {
    std::stringstream str;
    str << "Error opening wisdom file: " << filename;
    throw affineExcept("PDFactoryDouble", "addWisdom", str.str());
  }
  if (fftw_import_wisdom_from_file(fp) == 0) {
    std::stringstream str;
    str << "Error reading wisdom file: " << filename;
    throw affineExcept("PDFactoryDouble", "addWisdom", str.str());
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
  \param[in] n Number of elements
  \param[in] minflux1 Minimum flux band 1 to target in R
  \param[in] maxflux1 Maximum flux band 1 to use in R
  \param[in] minflux2 Minimum flux band 2 to target in R
  \param[in] maxflux2 Maximum flux band 2 to use in R

  Sets up RFlux1 and Rflux2.  They may not quite get to minflux in
  each band if negative, since we make small adjustments to make sure
  flux = 0 is included.
*/
void PDFactoryDouble::initRFlux(unsigned int n, double minflux1, 
				double maxflux1, double minflux2,
				double maxflux2) {
  // Make sure there is room
  resize(n);

  if (n == 0)
    throw affineExcept("PDFactoryDouble", "initR", "Invalid (0) n");
  if (n == 1) {
    dflux1 = dflux2 = 0.1;
    RFlux1[0] = minflux1;
    RFlux2[0] = minflux2;
    RwrapIdx1 = RwrapIdx2 = 1;
    return;
  }

  if (maxflux1 < minflux1) std::swap(minflux1, maxflux1);
  if (maxflux2 < minflux2) std::swap(minflux2, maxflux2);
  
  double st;
  initRFluxInternal(n, minflux1, maxflux1, RFlux1, dflux1, st, RwrapIdx1);
  initRFluxInternal(n, minflux2, maxflux2, RFlux2, dflux2, st, RwrapIdx2);

  rinitialized = false;
  initialized = false;
}


/*!
  \param[in] n Size of transforms
  \param[in] n        Size of transform 
  \param[in] minflux1 Minimum flux to use in R, band 1
  \param[in] maxflux1 Maximum flux to use in R, band 1
  \param[in] minflux2 Minimum flux to use in R, band 2
  \param[in] maxflux2 Maximum flux to use in R, band 2
  \param[in] model   number counts model to use for fill.  Params must be set
  \param[in] bm      Histogrammed inverse beam
  \param[in] setEdge Set the edge of R using edge integration
  \param[in] muldflux Multiply R by dflux1 * dflux2

  dflux1, dflux2 is also set, are are RFlux1, RFlux2 and the edge vars
*/
void PDFactoryDouble::initR(unsigned int n, double minflux1,
			    double maxflux1, double minflux2,
			    double maxflux2, const numberCountsDouble& model,
			    const doublebeam& bm, bool setEdge,
			    bool muldflux) {

  // This version is much more complex than the 1D case because of the
  // edge bits
  if (!rvars_allocated) allocateRvars();

  initRFlux(n, minflux1, maxflux1, minflux2, maxflux2);


  //Now fill in R.  The edges require special care (the part between
  // 0 and dflux in each dimension, which is the bottom edge of the
  // R array even with wrapping).  This
  // first call will fill nonsense values into the lower edges, which
  // we will later overwrite
  double *rptr;  

#ifdef TIMING
  std::clock_t starttime = std::clock();
#endif
  model.getR(n, RFlux1, n, RFlux2, bm, rvals);
#ifdef TIMING
  RTime += std::clock() - starttime;
#endif

  if (setEdge) {
    //Now fill in the lower edges by doing integrals
    // and setting to the mean value inside that.
    //Use the trapezoidal rule in either log or linear flux
    // space

    //Edge bits
    //Minimum values in integral; maximum are dflux1, dflux2
    double minedge1 = dflux1 *  PDFactoryDouble::lowEdgeRMult;
    double minedge2 = dflux2 *  PDFactoryDouble::lowEdgeRMult;
    double inedgem1 = 1.0 / static_cast<double>(nedge-1);
    double dinterpfluxedge1, dinterpfluxedge2;
    double iRxnorm = 0.0, iRynorm = 0.0, iR00norm = 0.0;
    
    if (setEdge) {
      if (!edgevars_allocated) allocateEdgevars();
      if (use_edge_log_x) {
	dinterpfluxedge1 = -log(PDFactoryDouble::lowEdgeRMult) * inedgem1;
	for (unsigned int i = 0; i < nedge; ++i)
	  REdgeFlux1[i] = minedge1 * 
	    exp(static_cast<double>(i) * dinterpfluxedge1);
      } else {
	dinterpfluxedge1 = (dflux1 - minedge1) * inedgem1;
	for (unsigned int i = 0; i < nedge; ++i)
	  REdgeFlux1[i] = minedge1 + static_cast<double>(i) * dinterpfluxedge1;
      }
      if (use_edge_log_y) {
	dinterpfluxedge2 = -log(PDFactoryDouble::lowEdgeRMult) * inedgem1;
	for (unsigned int i = 0; i < nedge; ++i)
	  REdgeFlux2[i] = minedge2 * 
	    exp(static_cast<double>(i) * dinterpfluxedge2);
      } else {
	dinterpfluxedge2 = (dflux2 - minedge2) * inedgem1;
	for (unsigned int i = 0; i < nedge; ++i)
	  REdgeFlux2[i] = minedge2 + static_cast<double>(i) * dinterpfluxedge2;
      }
      iRxnorm  = dinterpfluxedge1 / (dflux1 - minedge1);
      iRynorm  = dinterpfluxedge2 / (dflux2 - minedge2);
      iR00norm = dinterpfluxedge1 * dinterpfluxedge2 /
	((dflux1 - minedge1) * (dflux2 - minedge2));
    }

    //First, do r[0,0]
    double scriptr;
#ifdef TIMING
    starttime = std::clock();
#endif
    model.getR(nedge, REdgeFlux1, nedge, REdgeFlux2, bm, REdgeWork);
#ifdef TIMING
    RTime += std::clock() - starttime;
#endif

    //Do y integral first, store in REdgeWork[0,*]
    if (use_edge_log_y) {
      for (unsigned int i = 0; i < nedge; ++i) {
	rptr = REdgeWork + i * nedge; //row pointer
	scriptr = 0.5*(REdgeFlux2[0] * rptr[0] + 
		       REdgeFlux2[nedge-1] * rptr[nedge-1]);
	for (unsigned int j = 1; j < nedge-1; ++j)
	  scriptr += REdgeFlux2[j] * rptr[j];
	REdgeWork[i] = scriptr;
      }
    } else {
      for (unsigned int i = 0; i < nedge; ++i) {
	rptr = REdgeWork + i * nedge; 
	scriptr = 0.5*(rptr[0] + rptr[nedge-1]);
	for (unsigned int j = 1; j < nedge-1; ++j)
	  scriptr += rptr[j];
	REdgeWork[i] = scriptr;
      }
    }
    //Now X integral, put in integral step size and area
    // of bin, store in R[0,0]
    if (use_edge_log_x) {
      scriptr = 0.5*(REdgeFlux1[0] * REdgeWork[0]+
		     REdgeFlux1[nedge-1] * REdgeWork[nedge-1]);
      for (unsigned int i = 1; i < nedge-1; ++i)
	scriptr += REdgeFlux1[i] * REdgeWork[i];
      rvals[0] = scriptr * iR00norm;
    } else {
      scriptr = 0.5*(REdgeWork[0] + REdgeWork[nedge-1]);
      for (unsigned int i = 1; i < nedge-1; ++i)
	scriptr += REdgeWork[i];
      rvals[0] = scriptr * iR00norm;
    }
    
    //Now do Rx = R[0,y], integral along x
    double fixed_value;
    for (unsigned int j = 1; j < n; ++j) {
      fixed_value = RFlux2[j];
      //REdgeWork is more than big enough
#ifdef TIMING
      starttime = std::clock();
#endif
      model.getR(nedge, REdgeFlux1, 1, &fixed_value, bm, REdgeWork);
#ifdef TIMING
      RTime += std::clock() - starttime;
#endif
      if (use_edge_log_x) {
	scriptr = 0.5*(REdgeFlux1[0] * REdgeWork[0]+
		       REdgeFlux1[nedge-1] * REdgeWork[nedge-1]);
	for (unsigned int i = 1; i < nedge-1; ++i)
	  scriptr += REdgeFlux1[i] * REdgeWork[i];
      } else {
	scriptr = 0.5*(REdgeWork[0] + REdgeWork[nedge-1]);
	for (unsigned int i = 1; i < nedge-1; ++i)
	  scriptr += REdgeWork[i];
      }
      rvals[j] = scriptr * iRxnorm;
    }
    
    //And Ry = R[x,0]
    for (unsigned int i = 1; i < n; ++i) {
      fixed_value = RFlux1[i];
      //REdgeWork is more than big enough
#ifdef TIMING
      starttime = std::clock();
#endif
      model.getR(1, &fixed_value, nedge, REdgeFlux2, bm, REdgeWork);
#ifdef TIMING
      RTime += std::clock() - starttime;
#endif
      if (use_edge_log_y) {
	scriptr = 0.5*(REdgeFlux2[0] * REdgeWork[0]+
		       REdgeFlux2[nedge-1] * REdgeWork[nedge-1]);
	for (unsigned int j = 1; j < nedge-1; ++j)
	  scriptr += REdgeFlux2[j] * REdgeWork[j];
      } else {
	scriptr = 0.5*(REdgeWork[0] + REdgeWork[nedge-1]);
	for (unsigned int j = 1; j < nedge-1; ++j)
	  scriptr += REdgeWork[j];
      }
      rvals[i*n] = scriptr * iRynorm;
    }
  } else {
    //Just set edges to zero
    std::memset(rvals, 0, n * sizeof(double));
    for (unsigned int i = 1; i < n; ++i)
      rvals[i*n] = 0.0;
  }

  //Multiply R by dflux1 * dflux2
  if (muldflux) {
    double fluxfac = dflux1 * dflux2;
    for (unsigned int i = 0; i < n * n; ++i)
      rvals[i] *= fluxfac;
    rdflux = true;
  } else rdflux = false;
  
  initialized = false;
  rinitialized = true;
}

/*!
  Computes mean and variance of the P(D) using R along each dimension,
  then stores them in the internal state.
   
  You must have called initR before calling this.
*/
void PDFactoryDouble::getRStats() {

  if (!rinitialized)
    throw affineExcept("PDFactoryDouble", "getRMoment",
		       "R must be initialized before calling this");

  // Recall that the formulae for the mean and central 2nd moment
  // along each axis (not including instrumental noise) are
  //  <x> = \int x R dx dy
  //  <y> = \int y R dx dy
  //  <(x - <x>)^2> = \int x^2 R dx dy
  //  <(y - <y>)^2> = \int y^2 R dx dy
  // We ignore wrapping issues, since the contributions should be negligible
  //  there anyways, and that would signficantly increase the complexity
  //  of the code.

  // And... trap rule the integrals.
  // Always move along j since array access is
  //  faster that way (rvals is row-major order)
  // This is a simple calculation, somewhat tedious to write out.

  double cf1, cf2, cR, mean1, mean2, var1, var2, prod1, prod2;
  double *rowptr;
  
  // First row, has a 0.5 factor
  rowptr = rvals;
  cf1 = RFlux1[0];
  cR = 0.5 * rowptr[0]; // This handles the 0.5 factor
  prod1 = 0.5 * cf1 * cR; // Extra 0.5 for first col
  mean1 = prod1;
  var1 = cf1 * prod1;
  cf2 = RFlux2[0];
  prod2 = 0.5 * cf2 * cR;
  mean2 = prod2;
  var2 = cf2 * prod2;
  for (unsigned int j = 1; j < currsize-1; ++j) {
    cR = 0.5 * rowptr[j]; // Still 0.5 for first row
    prod1 = cf1 * cR;
    mean1 += prod1;
    var1 += cf1 * prod1;
    cf2 = RFlux2[j];
    prod2 = cf2 * cR;
    mean2 += prod2;
    var2 += cf2 * prod2;
  }
  cR = 0.5 * rowptr[currsize-1];
  prod1 = 0.5 * cf1 * cR;
  mean1 += prod1;
  var1 += cf1 * prod1;
  cf2 = RFlux2[currsize-1];
  prod2 = 0.5 * cf2 * cR;
  mean2 += prod2;
  var2 += cf2 * prod2;

  // Do middle rows
  for (unsigned int i = 1; i < currsize-1; ++i) {
    rowptr = rvals + i * currsize; 
    cf1 = RFlux1[i];
    cR = rowptr[0]; // No 0.5 any more -- we are in the middle
    prod1 = 0.5 * cf1 * cR; // 0.5 for first col
    mean1 += prod1;
    var1 += cf1 * prod1;
    cf2 = RFlux2[0];
    prod2 = 0.5 * cf2 * cR;
    mean2 += prod2;
    var2 += cf2 * prod2;
    for (unsigned int j = 1; j < currsize - 1; ++j) {
      cR = rowptr[j];
      prod1 = cf1 * cR;
      mean1 += prod1;
      var1 += cf1 * prod1;
      cf2 = RFlux2[j];
      prod2 = cf2 * cR;
      mean2 += prod2;
      var2 += cf2 * prod2;
    }
    cR = rowptr[currsize - 1];
    prod1 = 0.5 * cf1 * cR;
    mean1 += prod1;
    var1 += cf1 * prod1;
    cf2 = RFlux2[currsize - 1];
    prod2 = 0.5 * cf2 * cR;
    mean2 += prod2;
    var2 += cf2 * prod2;
  }

  // And final row, has an extra 0.5 factor again
  rowptr = rvals + (currsize - 1) * currsize;
  cf1 = RFlux1[currsize - 1];
  cR = 0.5 * rowptr[0]; 
  prod1 = 0.5 * cf1 * cR;
  mean1 += prod1;
  var1 += cf1 * prod1;
  cf2 = RFlux2[0];
  prod2 = 0.5 * cf2 * cR;
  mean2 += prod2;
  var2 += cf2 * prod2;
  for (unsigned int j = 1; j < currsize - 1; ++j) {
    cR = 0.5 * rowptr[j];
    prod1 = cf1 * cR;
    mean1 += prod1;
    var1 += cf1 * prod1;
    cf2 = RFlux2[j];
    prod2 = cf2 * cR;
    mean2 += prod2;
    var2 += cf2 * prod2;
  }
  cR = 0.5 * rowptr[currsize - 1];
  prod1 = 0.5 * cf1 * cR;
  mean1 += prod1;
  var1 += cf1 * prod1;
  cf2 = RFlux2[currsize - 1];
  prod2 = 0.5 * cf2 * cR;
  mean2 += prod2;
  var2 += cf2 * prod2;

  // Apply dflux factors if needed
  if (!rdflux) {
    double dfact = dflux1 * dflux2;
    mean1 *= dfact;
    mean2 *= dfact;
    var1 *= dfact;
    var2 *= dfact;
  }

  mn1 = mean1;
  mn2 = mean2;
  var_noi1 = var1;
  var_noi2 = var2;
}

/*!
  \param[in] n Transform size
*/
void PDFactoryDouble::setupTransforms(unsigned int n) {

  //The transform plan is that the forward transform
  // takes rvals to rtrans.  We then do things to rtrans
  // to turn it into pval, including shifting, adding noise,
  // etc, and then transform from pval to pofd.

  if (n == 0)
    throw affineExcept("PDFactoryDouble", "setupTransforms",
		       "Invalid (0) transform size");

  // Make sure we have enough room
  resize(n);

  //We need the R variables allocated so we can plan to them
  if (!rvars_allocated) allocateRvars();

  int intn = static_cast<int>(n);
  if (plan == NULL) {
    plan = fftw_plan_dft_r2c_2d(intn, intn, rvals, rtrans,
				fftw_plan_style);
  }
  if (plan == NULL) {
    std::stringstream str;
    str << "Plan creation failed for forward transform of size: " << n;
    if (has_wisdom) str << std::endl << "Your wisdom file may not have"
			<< " that size";
    throw affineExcept("PDFactoryDouble", "setupTransforms", str.str());
  }
  if (plan_inv != NULL) {
    plan_inv = fftw_plan_dft_c2r_2d(intn, intn, pval, pofd,
				    fftw_plan_style);
    if (plan_inv == NULL) {
      std::stringstream str;
      str << "Plan creation failed for inverse transform of size: " << intn
	  << " by " << intn;
      if (has_wisdom) str << std::endl << "Your wisdom file may not have"
			  << " that size";
      throw affineExcept("PDFactoryDouble", "setupTransforms", str.str());
    }
  }
}

/*!
  \param[out] pd Holds P(D) on output, normalized, mean subtracted,
                 and with positivity enforced.

  This does the unwrapping of the internal P(D) into pd.
*/
// This should only ever be called by getPD, so we don't really
// check the inputs
void PDFactoryDouble::unwrapPD(PDDouble& pd) const {
  // Our acceptance testing is a bit complicated.
  // If the minimum point is more than nsig1 away from the expected mean (0)
  //  then we just accept it.  If it is more than nsig2 away, we make sure
  //  that the min/max ratio of the P(D) along that axis is more than
  //  maxminratio.  If it is less than nsig2, we just flat out reject.
  const double nsig1 = 4.0; 
  const double nsig2 = 2.0;
  const double maxminratio = 1e5;

  //Enforce positivity
#ifdef TIMING
  starttime = std::clock();
#endif
  for (unsigned int idx = 0; idx < currsize * currsize; ++idx)
    if (pofd[idx] < 0) pofd[idx] = 0.0;
#ifdef TIMING
  posTime += std::clock() - starttime;
#endif

  // Figure out the indices of the minimum in each dimension
  // Rather than find the 2D index, which would be noisier, 
  //  intergrate out each axis to find the minimum.
  // Start from the top and move down -- if there is a tie for
  //  some reason we prefer the index be high because that is more
  //  likely to be right in practice
  // Start with min along the first index.  We also find the max, but
  // don't keep track of its index
#ifdef TIMING
  starttime = std::clock();
#endif
  int nm1 = static_cast<int>(currsize - 1);
  int mdx = nm1; // Curr min index
  double cval, minval, maxval, *rowptr; 
  // Initialize minval to last col
  rowptr = pofd + nm1 * currsize;
  minval = 0.5 * rowptr[0];
  for (unsigned int j = 1; j < currsize - 1; ++j) minval += rowptr[j];
  minval += 0.5 * rowptr[currsize - 1];
  maxval = minval;
  // Now do others
  for (int i = currsize - 2; i >= 0; --i) {
    rowptr = pofd + i * currsize;
    cval = 0.5 * rowptr[0];
    for (unsigned int j = 1; j < currsize - 1; ++j) cval += rowptr[j];
    cval += 0.5 * rowptr[currsize - 1];
    if (cval < minval) {
      minval = cval;
      mdx = i;
    } else if (cval > maxval) maxval = cval;
  }
  unsigned int minidx1 = static_cast<unsigned>(mdx);

  // Sanity check
  double fwrap_plus = RFlux1[minidx1]; // Wrap in pos flux
  double fwrap_minus = 
    static_cast<double>(currsize - minidx1) * dflux1; // Abs neg wrap
  double cs1, cs2;
  cs1 = nsig1 * sg1;
  cs2 = nsig2 * sg1;
  if ((fwrap_plus > cs1) || (fwrap_minus > cs1)) {
    // Worth further investigation
    if (fwrap_plus < cs2) {
      std::stringstream errstr;
      errstr << "Top wrapping problem dim 1; wrapping point at "
	     << fwrap_plus << " which is only " << fwrap_plus / sg1
	     << " sigma away from expected (0) mean with sigma " << sg1;
      throw affineExcept("PDFactoryDouble", "unwrapPD", errstr.str());
    }
    if (fwrap_minus < cs2) {
      std::stringstream errstr;
      errstr << "Bottom wrapping problem dim 1; wrapping point at "
	     << -fwrap_minus << " which is only " << fwrap_minus / sg1
	     << " sigma away from expected (0) mean, with sigma " << sg1;
      throw affineExcept("PDFactoryDouble", "unwrapPD", errstr.str());
    } 
    // Min/max ratio test
    if (maxval / minval < maxminratio) {
      std::stringstream errstr;
      errstr << "Dim 1 wrapping problem with wrapping fluxes: "
	     << fwrap_plus << " and " << -fwrap_minus << " with min/max ratio: "
	     << maxval / minval << " and sigma: " << sg1;
      throw affineExcept("PDFactoryDouble", "unwrapPD", errstr.str());
    }
  }
  
  // Now second dimension.  This one is slower due to stride issues.
  // That could be improved with an auxilliary array, but it doesn't
  // seem worth it, frankly
  mdx = nm1; // Index of current minimum
  minval = 0.5 * pofd[nm1];
  for (unsigned int i = 1; i < currsize - 1; ++i) 
    minval += pofd[i * currsize + nm1];
  minval += 0.5 * pofd[nm1 * currsize + nm1];
  maxval = minval;
  for (int j = currsize - 2; j >= 0; --j) {
    cval = 0.5 * pofd[j];
    for (unsigned int i = 1; i < currsize - 1; ++i) 
      cval += pofd[i * currsize + j];
    cval += 0.5 * pofd[nm1 * currsize + j];
    if (cval < minval) {
      minval = cval;
      mdx = j;
    } else if (cval > maxval) maxval = cval;
  }
  unsigned int minidx2 = static_cast<unsigned>(mdx);

  // Same sanity check
  fwrap_plus = RFlux2[minidx2]; 
  fwrap_minus = static_cast<double>(currsize - minidx2) * dflux2;
  cs1 = nsig1 * sg2;
  cs2 = nsig2 * sg2;
  if ((fwrap_plus > cs1) || (fwrap_minus > cs1)) {
    // Worth further investigation
    if (fwrap_plus < cs2) {
      std::stringstream errstr;
      errstr << "Top wrapping problem dim 2; wrapping point at "
	     << fwrap_plus << " which is only " << fwrap_plus / sg2
	     << " sigma away from expected (0) mean with sigma " << sg2;
      throw affineExcept("PDFactoryDouble", "unwrapPD", errstr.str());
    }
    if (fwrap_minus < cs2) {
      std::stringstream errstr;
      errstr << "Bottom wrapping problem dim 2; wrapping point at "
	     << -fwrap_minus << " which is only " << fwrap_minus / sg2
	     << " sigma away from expected (0) mean, with sigma " << sg2;
      throw affineExcept("PDFactoryDouble", "unwrapPD", errstr.str());
    } 
    // Min/max ratio test
    if (maxval / minval < maxminratio) {
      std::stringstream errstr;
      errstr << "Dim 2 wrapping problem with wrapping fluxes: "
	     << fwrap_plus << " and " << -fwrap_minus << " with min/max ratio: "
	     << maxval / minval << " and sigma: " << sg2;
      throw affineExcept("PDFactoryDouble", "unwrapPD", errstr.str());

    }
  }

  // Now the actual copying, which is an exercise in index gymnastics
  pd.resize(currsize, currsize);
  size_t size_double = sizeof(double); // Size of double
  std::memset(pd.pd_, 0, currsize * currsize * size_double);

  size_t rowsz;
  double *ptr_curr, *ptr_out; // convenience pointers
  // We start with the neg, neg bit -- that is stuff >= both
  // minidx1 and minidx2, which ends up in the 0, 0 part of the output
  ptr_out = pd.pd_;
  ptr_curr = pofd + minidx1 * currsize + minidx2;
  rowsz = (currsize - minidx2) * size_double;
  for (unsigned int i = 0; i < currsize - minidx1; ++i)
    std::memcpy(ptr_out + i * currsize, ptr_curr + i * currsize, rowsz);

  // Next pos neg
  ptr_out = pd.pd_ + (currsize - minidx1) * currsize;
  ptr_curr = pofd + minidx2;
  for (unsigned int i = 0; i < minidx1; ++i) 
    std::memcpy(ptr_out + i * currsize, ptr_curr + i * currsize, rowsz);

  // pos, pos
  ptr_out = pd.pd_ + (currsize - minidx1) * currsize + (currsize - minidx2);
  ptr_curr = pofd;
  rowsz = minidx2 * size_double;
  for (unsigned int i = 0; i < minidx1; ++i)
    std::memcpy(ptr_out + i * currsize, ptr_curr + i * currsize, rowsz);
  
  // neg, pos
  ptr_out = pd.pd_ + (currsize - minidx2);
  ptr_curr = pofd + minidx1 * currsize;
  for (unsigned int i = 0; i < currsize - minidx1; ++i)
    std::memcpy(ptr_out + i * currsize, ptr_curr + i * currsize, rowsz);

  pd.logflat = false;
  pd.minflux1 = 0.0; pd.dflux1 = dflux1;
  pd.minflux2 = 0.0; pd.dflux2 = dflux2;

#ifdef TIMING
  copyTime += std::clock() - starttime;
#endif

  // Normalize
#ifdef TIMING
  starttime = std::clock();
#endif
  pd.normalize();
#ifdef TIMING
  normTime += std::clock() - starttime;
#endif

  // Mean subtract axes
#ifdef TIMING
  starttime = std::clock();
#endif
  double tmn1, tmn2; //True means
  pd.getMeans(tmn1, tmn2, false);
  if (std::isinf(tmn1) || std::isnan(tmn1) || std::isinf(tmn2) ||
       std::isnan(tmn2)) {
    std::stringstream str;
    str << "Un-shift amounts not finite band1: " << tmn1 << " "
        << tmn2 << std::endl;
    str << "At length: " << currsize << " with sigma: " << sg1 << " "
        << sg2;
    throw affineExcept("PDFactoryDouble", "unwrapPD", str.str());
  }
  pd.minflux1 = -tmn1;
  pd.minflux2 = -tmn2;
#ifdef TIMING
  meanTime += std::clock() - starttime;
#endif

}


/*!
  \param[in] n Size of transform (square)
  \param[in] minflux1 Minimum flux generated for R, dimension 1
  \param[in] maxflux1 Maximum flux generated for R, dimension 1
  \param[in] minflux2 Minimum flux generated for R, dimension 2
  \param[in] maxflux2 Maximum flux generator for R, dimension 2
  \param[in] model    number counts model to use for fill.  Params must be set
  \param[in] bm       Beam
  \param[in] setEdge  Use integral of mean values at the edges

  \returns True if the P(D) could be initialized, false if something
           about the parameters prevented initialization.  Note that
	   a genuine error results in throwing an exception, not setting this
	   to false.

  Gets ready for P(D) computation by preparing R.
  The idea is to compute everything that doesn't require knowing the
  exact sigmas here so we can call this multiple times on
  maps with the same beam but different noise values.
*/
bool PDFactoryDouble::initPD(unsigned int n, double minflux1, double maxflux1, 
			     double minflux2, double maxflux2, 
			     const numberCountsDouble& model,
			     const doublebeam& bm, bool setEdge) {

  if (n == 0)
    throw affineExcept("PDFactoryDouble", "initPD", 
		       "Invalid (non-positive) n");
  if (maxflux1 <= 0.0)
    throw affineExcept("PDFactoryDouble", "initPD",
		       "Invalid (non-positive) maxflux1");
  if (maxflux2 <= 0.0)
    throw affineExcept("PDFactoryDouble", "initPD",
		       "Invalid (non-positive) maxflux2");
  
  initialized = false;

  //Make the plans, or keep the old ones if possible
  // Note we have to do this before we fill R, as plan construction
  // may overwrite the values.  This allocates R and resizes
  setupTransforms(n);

  // Compute R, including multiplying it by dflux
  initR(n, minflux1, maxflux2, minflux1, maxflux2, model, bm, 
	setEdge, true);

  // Get mean and variance estimates from R
  getRStats();

  // Set the sigmas to the no-noise values for now
  sg1 = sqrt(var_noi1);
  sg2 = sqrt(var_noi2);

  //Decide if we will shift and pad, and if so by how much.
  // The idea is to shift the mean to zero -- but we only
  // do the shift if the sigma is larger than one actual step size
  // because otherwise we can't represent it well.
  doshift1 = (sg1 > dflux1) && (fabs(mn1) > dflux1);
  if (doshift1) shift1 = -mn1; else shift1 = 0.0;
  doshift2 = (sg2 > dflux2) && (fabs(mn2) > dflux2);
  if (doshift2) shift2 = -mn2; else shift2 = 0.0;

  if (verbose) {
    std::cout << " Initial mean estimate band1: " << mn1 << " band2: "
	      << mn2 << std::endl;
    std::cout << " Initial stdev estimate band1: " << sg1 << " band2: "
	      << sg2 << std::endl;
    if (doshift1)
      std::cout << " Additional shift in band 1: " << shift1 << std::endl;
    else std::cout << " Not applying additional shift in band 1" << std::endl;
    if (doshift2)
      std::cout << " Additional shift in band 2: " << shift2 << std::endl;
    else std::cout << " Not applying additional shift in band 2" << std::endl;
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
  Calculates P(D) for a dataset using direct fills
 
  \param[in] sigma1 Sigma in dimension 1
  \param[in] sigma2 Sigma in dimension 2
  \param[out] pd Holds P(D1,D2) on output
  \param[in] setLog If true, pd is log(P(D1,D2)) on output; convenient
              for likelihood evaluation.

  You must call initPD first.
*/
void PDFactoryDouble::getPD(double sigma1, double sigma2,  
			    PDDouble& pd, bool setLog) {

  // The basic idea is to compute the P(D) from the previously filled
  // R values, adding in noise and all that fun stuff, filling pd
  // for output

  if (!initialized )
    throw affineExcept("PDFactoryDouble", "getPD",
		       "Must call initPD first");

  //Output array from 2D FFT is n * (n/2+1)
  unsigned int n = currsize;
  unsigned int ncplx = n / 2 + 1;
      
  //Calculate p(omega) = exp( r(omega1,omega2) - r(0,0) ),
  // which is what we will transform back into pofd.
  // There are some complications because of shifts and all that.
  //The 2D output real FFT format makes this a bit tricky.  The output
  // array is n by (n/2+1).  The first dimension has both
  // positive and negative frequencies, the second dimension has
  // only positive dimensions.
  //In particular, if i is the index over the first dimension,
  // the frequencies along the first dimension are:
  //  f1 = i/delta1*n        for i = [0,n/2]
  //     = - (n-i)/delta1*n  for i = [n/2+1,n-1]
  // where delta1 is dflux1.
  // And if j is the second dimension index, then
  //  f2 = j/delta2*n  
  // and delta2 = dflux2.  We work in w instead of f (2 pi f)
  // and actually compute 
  //  exp( r(omega1,omega2) - r(0,0) - i*shift1*omega1 - i*shift2*omega2
  //       - 1/2 sigma1^2 omega1^2 - 1/2 sigma2^2 * omega2^2)
  
  double r0, expfac, rval, ival;
  fftw_complex *row_current_out; //Output out variable
  fftw_complex *r_input_rowptr; //Row pointer into out_part (i.e., input)
  r0 = rtrans[0][0]; //r[0,0] is pure real
  
  double iflux1 = mcmc_affine::two_pi / (n * dflux1);
  double iflux2 = mcmc_affine::two_pi / (n * dflux2);
  
  if (verbose) std::cout << "  Computing p(w1,w2)" << std::endl;

#ifdef TIMING
  std::clock_t starttime = std::clock();
#endif

  // Now we have 16 potential cases -- both shift bits, both with
  //  and without sigma in each.  However, in practice it is unlikely
  //  that we are shifting in one band and not the other, or have noise
  //  in one band and not the other, so we will pay the small inefficiency
  //  price for those cases by only exploring the 4 cases
  //  (any shift vs. noshift) by (any sigma vs. no sigma).
  if (doshift1 || doshift2) {
    double meanprod1, didx2, w1, w2;
    if (sigma1 > 0 || sigma2 > 0) {
      double sigfac1 = 0.5 * sigma1 * sigma1;
      double sigfac2 = 0.5 * sigma2 * sigma2;
      double sigprod1;
      //First, Pos freq
      for (unsigned int idx1 = 0; idx1 < ncplx; ++idx1) {
	r_input_rowptr = rtrans + idx1 * ncplx; //Input
	row_current_out = pval + idx1 * ncplx; //Output
	w1 = iflux1 * static_cast<double>(idx1);
	meanprod1 = shift1 * w1;
	sigprod1  = sigfac1 * w1 * w1;
	for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
	  didx2 = static_cast<double>(idx2);
	  w2 = iflux2 * didx2;
	  rval = r_input_rowptr[idx2][0] - r0 - sigprod1 - sigfac2 * w2 * w2;
	  ival = r_input_rowptr[idx2][1] - meanprod1 - shift2 * w2;
	  expfac = exp(rval);
	  row_current_out[idx2][0] = expfac * cos(ival);
	  row_current_out[idx2][1] = expfac * sin(ival);
	}
      }
      //Now, Neg freq
      for (unsigned int idx1 = ncplx; idx1 < n; ++idx1) {
	r_input_rowptr = rtrans + idx1 * ncplx; //Input
	row_current_out = pval + idx1 * ncplx; //Output
	w1 = - iflux1 * static_cast<double>(n - idx1);
	meanprod1 = shift1 * w1;
	sigprod1  = sigfac1 * w1 * w1;
	for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
	  didx2 = static_cast<double>(idx2);
	  w2 = iflux2 * didx2;
	  rval = r_input_rowptr[idx2][0] - r0 - sigprod1 - sigfac2 * w2 * w2;
	  ival = r_input_rowptr[idx2][1] - meanprod1 - shift2 * w2;
	  expfac = exp(rval);
	  row_current_out[idx2][0] = expfac * cos(ival);
	  row_current_out[idx2][1] = expfac * sin(ival);
	}
      }
    } else {
      // Shifts, but no sigmas
      for (unsigned int idx1 = 0; idx1 < ncplx; ++idx1) {
	r_input_rowptr = rtrans + idx1 * ncplx; //Input
	row_current_out = pval + idx1 * ncplx; //Output
	w1 = iflux1 * static_cast<double>(idx1);
	meanprod1 = shift1 * w1;
	for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
	  didx2 = static_cast<double>(idx2);
	  w2 = iflux2 * didx2;
	  rval = r_input_rowptr[idx2][0] - r0;
	  ival = r_input_rowptr[idx2][1] - meanprod1 - shift2 * w2;
	  expfac = exp(rval);
	  row_current_out[idx2][0] = expfac * cos(ival);
	  row_current_out[idx2][1] = expfac * sin(ival);
	}
      }
      //Now, Neg freq
      for (unsigned int idx1 = ncplx; idx1 < n; ++idx1) {
	r_input_rowptr = rtrans + idx1 * ncplx; //Input
	row_current_out = pval + idx1 * ncplx; //Output
	w1 = - iflux1 * static_cast<double>(n - idx1);
	meanprod1 = shift1 * w1;
	for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
	  didx2 = static_cast<double>(idx2);
	  w2 = iflux2 * didx2;
	  rval = r_input_rowptr[idx2][0] - r0;
	  ival = r_input_rowptr[idx2][1] - meanprod1 - shift2 * w2;
	  expfac = exp(rval);
	  row_current_out[idx2][0] = expfac * cos(ival);
	  row_current_out[idx2][1] = expfac * sin(ival);
	}
      }
    }
  } else {
    // No shifts, but there may be sigmas
    // TODO -- is there a symmetry here which simplifies this case,
    //  since the sigmas are quadratic?
    if (sigma1 > 0 || sigma2 > 0) {
      double w1, w2, didx2;
      double sigfac1 = 0.5 * sigma1 * sigma1;
      double sigfac2 = 0.5 * sigma2 * sigma2;
      double sigprod1;
      //First, Pos freq
      for (unsigned int idx1 = 0; idx1 < ncplx; ++idx1) {
	r_input_rowptr = rtrans + idx1 * ncplx; //Input
	row_current_out = pval + idx1 * ncplx; //Output
	w1 = iflux1 * static_cast<double>(idx1);
	sigprod1  = sigfac1 * w1 * w1;
	for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
	  didx2 = static_cast<double>(idx2);
	  w2 = iflux2 * didx2;
	  rval = r_input_rowptr[idx2][0] - r0 - sigprod1 - sigfac2 * w2 * w2;
	  ival = r_input_rowptr[idx2][1];
	  expfac = exp(rval);
	  row_current_out[idx2][0] = expfac * cos(ival);
	  row_current_out[idx2][1] = expfac * sin(ival);
	}
      }
      //Now, Neg freq
      for (unsigned int idx1 = ncplx; idx1 < n; ++idx1) {
	r_input_rowptr = rtrans + idx1 * ncplx; //Input
	row_current_out = pval + idx1 * ncplx; //Output
	w1 = - iflux1 * static_cast<double>(n - idx1);
	for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
	  didx2 = static_cast<double>(idx2);
	  w2 = iflux2 * didx2;
	  rval = r_input_rowptr[idx2][0] - r0 - sigprod1 - sigfac2 * w2 * w2;
	  ival = r_input_rowptr[idx2][1];
	  expfac = exp(rval);
	  row_current_out[idx2][0] = expfac * cos(ival);
	  row_current_out[idx2][1] = expfac * sin(ival);
	}
      }
    } else {
      //No shifts, sigmas
      for (unsigned int idx1 = 0; idx1 < n; ++idx1) {
	r_input_rowptr = rtrans + idx1 * ncplx;
	row_current_out = pval + idx1 * ncplx;
	for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
	  rval = r_input_rowptr[idx2][0] - r0;
	  ival = r_input_rowptr[idx2][1];
	  expfac = exp(rval);
	  row_current_out[idx2][0] = expfac * cos(ival);
	  row_current_out[idx2][1] = expfac * sin(ival);
	}
      }
    }
  }

  //p(0,0) is special, but shared by all cases.
  pval[0][0] = 1.0;
  pval[0][1] = 0.0;

#ifdef TIMING
  p0Time += std::clock() - starttime;
#endif

  //Now transform back
  if (verbose) std::cout << " Reverse transform" << std::endl;
#ifdef TIMING
  starttime = std::clock();
#endif
  // From pval into pofd
  fftw_execute_dft_c2r(plan, pval, pofd);
#ifdef TIMING
  fftTime += std::clock() - starttime;
#endif

  // Copy into output variable, also normalizing, mean subtracting, 
  // making positive
  sg1 = sqrt(var_noi1 + sigma1 * sigma1); // Used in unwrapPD
  sg2 = sqrt(var_noi2 + sigma2 * sigma2);
  unwrapPD(pd);

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
  \param[in] comm Communicator
  \param[in] dest Destination of messages
*/
void PDFactoryDouble::sendSelf(MPI_Comm comm, int dest) const {
  MPI_Send(const_cast<unsigned int*>(&fftw_plan_style), 1, MPI_UNSIGNED,
	   dest, pofd_mcmc::PDFDSENDPLANSTYLE, comm);
  MPI_Send(const_cast<bool*>(&has_wisdom), 1, MPI::BOOL, dest,
	   pofd_mcmc::PDFDHASWISDOM, comm);
  if (has_wisdom) {
    //Send wisdom file name
    unsigned int nstr = wisdom_file.size()+1;
    char *cstr = new char[nstr];
    std::strncpy( cstr, wisdom_file.c_str(), nstr );
    MPI_Send(&nstr,1,MPI_UNSIGNED,dest,pofd_mcmc::PDFDWISLEN, comm);
    MPI_Send(cstr, nstr, MPI_CHAR, dest,pofd_mcmc::PDFDWISNAME, comm);
    delete[] cstr;
  }
  MPI_Send(const_cast<bool*>(&verbose), 1, MPI::BOOL, dest,
	   pofd_mcmc::PDFDVERBOSE, comm);
  MPI_Send(const_cast<unsigned int*>(&nedge), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::PDFDNEDGE, comm);
}

/*!
  \param[in] comm Communicator
  \param[in] src Source of messages

  Doesn't copy over interal variables
*/
void PDFactoryDouble::recieveCopy(MPI_Comm comm, int src) {
  MPI_Status Info;
  MPI_Recv(&fftw_plan_style, 1, MPI_UNSIGNED, src, 
	   pofd_mcmc::PDFDSENDPLANSTYLE, comm, &Info);
  MPI_Recv(&has_wisdom, 1, MPI::BOOL, src, pofd_mcmc::PDFDHASWISDOM,
	   comm, &Info);
  if (has_wisdom) {
    //Recieve wisdom file name
    unsigned int nstr;
    MPI_Recv(&nstr, 1, MPI_UNSIGNED, src, pofd_mcmc::PDFDWISLEN, comm, &Info);
    char *cstr = new char[nstr];
    MPI_Recv(cstr, nstr, MPI_CHAR, src, pofd_mcmc::PDFDWISNAME, comm, &Info); 
    wisdom_file = std::string(cstr);
    delete[] cstr;
    addWisdom(wisdom_file);
  }
  MPI_Recv(&verbose, 1, MPI::BOOL, src, pofd_mcmc::PDFDVERBOSE, comm, &Info);
  unsigned int newnedge;
  MPI_Recv(&newnedge, 1, MPI_UNSIGNED, src, pofd_mcmc::PDFDNEDGE, comm, &Info);
  if (newnedge != nedge) {
    freeEdgevars();
    nedge   = newnedge;
  }
  rinitialized = false;
  initialized = false;
}
