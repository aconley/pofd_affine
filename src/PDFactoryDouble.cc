#include<sstream>
#include<cstdlib>
#include<cmath>
#include<cstring>
#include<limits>
#include<iomanip>

#include "../include/global_settings.h"
#include "../include/PDFactoryDouble.h"
#include "../include/hdf5utils.h"
#include "../include/affineExcept.h"

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
    dflux = maxflux / (n - floor(-minflux / dflux) - 1.0);
    // Figure out what index we go up to with positive fills
    wrapidx = static_cast<unsigned int>(maxflux / dflux + 0.999999999);
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
  if (RFlux1 != nullptr) fftw_free(RFlux1);
  if (RFlux2 != nullptr) fftw_free(RFlux2);

  if (rvals != nullptr) fftw_free(rvals);
  if (rsum != nullptr) fftw_free(rsum);
  if (prtrans != nullptr) fftw_free(prtrans);
  if (pofd != nullptr) fftw_free(pofd);
  if (pval != nullptr) fftw_free(pval);

  if (REdgeFlux1 != nullptr) fftw_free(REdgeFlux1);
  if (REdgeFlux2 != nullptr) fftw_free(REdgeFlux2);
  if (RxEdgeWork != nullptr) fftw_free(RxEdgeWork);
  if (RyEdgeWork != nullptr) fftw_free(RyEdgeWork);

  if (plan != nullptr) fftw_destroy_plan(plan); 
  if (plan_inv != nullptr) fftw_destroy_plan(plan_inv);
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
  RFlux1 = RFlux2 = nullptr;
  RwrapIdx1 = RwrapIdx2 = 0;
  rdflux = false;
  rvals = nullptr;
  rsum = nullptr;
  prtrans = nullptr;
  pofd = nullptr;
  pval = nullptr;

  nedge = NEDGE;
  edgevars_allocated = false;
  nedgework = 0;
  REdgeFlux1 = nullptr;
  REdgeFlux2 = nullptr;
  RxEdgeWork = nullptr;
  RyEdgeWork = nullptr;

  dflux1 = dflux2 = 0.0;

  plan = plan_inv = nullptr;

  verbose = false;
  has_wisdom = false;
  fftw_plan_style = FFTW_MEASURE;

  mn1 = mn2 = var_noi1 = var_noi2 = std::numeric_limits<double>::quiet_NaN();
  sg1pos = sg1neg = sg2pos = sg2neg = std::numeric_limits<double>::quiet_NaN();
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
  RTime = p0Time = fftTime = posTime = copyTime = splitTime = normTime = 0;
  edgeTime = meanTime = logTime = 0;
}

/*!
  \param[in] nindent Number of indentation spaces
*/
void PDFactoryDouble::summarizeTime(unsigned int nindent) const {
  std::string prestring(nindent, ' ');
    
  std::cout << "R time: " << prestring 
            << 1.0*RTime/CLOCKS_PER_SEC << "s" << std::endl;
  std::cout << "p0 time: " << prestring 
            << 1.0*p0Time/CLOCKS_PER_SEC << "s" << std::endl;
  std::cout << "fft time: " << prestring 
            << 1.0*fftTime/CLOCKS_PER_SEC << "s" << std::endl;
  std::cout << "pos time: " << prestring 
            << 1.0*posTime/CLOCKS_PER_SEC << "s" << std::endl;
  std::cout << "split time: " << prestring 
            << 1.0*splitTime/CLOCKS_PER_SEC << "s" << std::endl;
  std::cout << "copy time: " << prestring 
            << 1.0*copyTime/CLOCKS_PER_SEC << "s" << std::endl;
  std::cout << "norm time: " << prestring 
            << 1.0*normTime/CLOCKS_PER_SEC << "s" << std::endl;
  std::cout << "mean time: " << prestring 
            << 1.0*meanTime/CLOCKS_PER_SEC << "s" << std::endl;
  std::cout << "log time: " << prestring 
            << 1.0*logTime/CLOCKS_PER_SEC << "s" << std::endl;
}
#endif

/*
  \param[in] NSIZE new size
*/
void PDFactoryDouble::resize(unsigned int NSIZE) {
  if (NSIZE == currsize) return;
  freeRvars();
  currsize = NSIZE;
  if (plan != nullptr) { fftw_destroy_plan(plan); plan = nullptr; }
  if (plan_inv != nullptr) {
    fftw_destroy_plan(plan_inv);
    plan_inv = nullptr;
  }  
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
  rsum = (double*) fftw_malloc(sizeof(double) * currsize);
  unsigned int fsize = currsize * currsize;
  rvals = (double*) fftw_malloc(sizeof(double) * fsize);
  pofd  = (double*) fftw_malloc(sizeof(double) * fsize);
  fsize = currsize * (currsize / 2 + 1);
  prtrans = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fsize);
  pval = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fsize);

  rvars_allocated = true;
  rinitialized = false;
  initialized = false;
}

// Also doesn't necessarily invalidate the plans
void PDFactoryDouble::freeRvars() {
  if (RFlux1 != nullptr) { fftw_free(RFlux1); RFlux1=nullptr; }
  if (RFlux2 != nullptr) { fftw_free(RFlux2); RFlux2=nullptr; }
  if (rvals != nullptr)  { fftw_free(rvals); rvals=nullptr; }
  if (rsum != nullptr)   { fftw_free(rsum); rsum=nullptr; }
  if (prtrans != nullptr) { fftw_free(prtrans); prtrans=nullptr; }
  if (pval != nullptr)   { fftw_free(pval); pval = nullptr; }
  if (pofd != nullptr)   { fftw_free(pofd); pofd = nullptr; }
  rvars_allocated = false;
  rinitialized = false;
  initialized = false;
}

// Note this doesn't affect the initialization state for the 
// R forward transform because the stuff done in the edge vars
// is copied over before the transform
void PDFactoryDouble::allocateEdgevars() {
  // Deal with edge variables.  Two types: REdgeFlux, 
  if (!edgevars_allocated) {
    if (nedge > 0) {
      REdgeFlux1 = (double*) fftw_malloc(sizeof(double) * nedge);
      REdgeFlux2 = (double*) fftw_malloc(sizeof(double) * nedge);
      edgevars_allocated = true;
    } else {
      REdgeFlux1 = REdgeFlux2 = nullptr;
      edgevars_allocated = false;
    }
  }
  if (currsize > nedgework) { // Don't downsize
    // Resize (or allocate) REdgeWork[12]
    if (RxEdgeWork != nullptr) fftw_free(RxEdgeWork);
    if (RyEdgeWork != nullptr) fftw_free(RyEdgeWork);
    if (currsize > 0) {
      RxEdgeWork = (double*) fftw_malloc(sizeof(double) * currsize);
      RyEdgeWork = (double*) fftw_malloc(sizeof(double) * currsize);
    } else RxEdgeWork = RyEdgeWork = nullptr;
    nedgework = currsize;
  }
}

void PDFactoryDouble::freeEdgevars() {
  if (REdgeFlux1 != nullptr) { fftw_free(REdgeFlux1); REdgeFlux1 = nullptr; }
  if (REdgeFlux2 != nullptr) { fftw_free(REdgeFlux2); REdgeFlux2 = nullptr; }
  if (RxEdgeWork != nullptr) { fftw_free(RxEdgeWork); RxEdgeWork = nullptr; }
  if (RyEdgeWork != nullptr) { fftw_free(RyEdgeWork); RyEdgeWork = nullptr; }
  edgevars_allocated = false;
  nedgework = 0;
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
  FILE *fp = nullptr;
  fp = fopen(filename.c_str(), "r");
  if (fp == nullptr) {
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
  if (plan != nullptr) {
    fftw_destroy_plan(plan); 
    plan = nullptr;
  }
  if (plan_inv != nullptr) {
    fftw_destroy_plan(plan_inv);
    plan_inv = nullptr;
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
  resize(n); // Set currsize if needed
  if (!rvars_allocated) allocateRvars();

  initRFlux(n, minflux1, maxflux1, minflux2, maxflux2);

  // Now fill in R.  The edges require special care (the part between
  // 0 and dflux in each dimension, which is the bottom edge of the
  // R array even with wrapping).  We do the edges first, then
  // the main array.  We use rvals as working space for the edges,
  // which means we have to store the edge values somewhere to later
  // dump into the main array (since getR on the main array will fill
  // the edges with junk).

  double r00val = 0.0; // Will hold R[0, 0] if computed
  if (setEdge) {
    // Since we will use rvals for temporary storage, it better have
    //  enough room!
    if (n < nedge)
      throw affineExcept("PDFactoryDouble", "initR",
                         "Nedge is less than R fill size");

    //Now fill in the lower edges by doing integrals
    // and setting to the mean value inside that.
    //Use the trapezoidal rule in either log or linear flux
    // space
    double *rptr;

    //Edge bits
    //Minimum values in integral; maximum are dflux1, dflux2
    double minedge1 = dflux1 *  PDFactoryDouble::lowEdgeRMult;
    double minedge2 = dflux2 *  PDFactoryDouble::lowEdgeRMult;
    double inedgem1 = 1.0 / static_cast<double>(nedge-1);
    double dinterpfluxedge1, dinterpfluxedge2;
    double iRxnorm = 0.0, iRynorm = 0.0, iR00norm = 0.0;
    
    // Allocate R vars, resize REdgeWorks[12] if needed
    allocateEdgevars();
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

    //First, do r[0,0]
    double scriptr;
#ifdef TIMING
    starttime = std::clock();
#endif
    // Stored into rvals
    model.getR(nedge, REdgeFlux1, nedge, REdgeFlux2, bm, rvals);
#ifdef TIMING
    RTime += std::clock() - starttime;
#endif

    //Do y integral first, store in RxEdgeWork as temporary
    if (use_edge_log_y) {
      for (unsigned int i = 0; i < nedge; ++i) {
        rptr = rvals + i * nedge; //row pointer
        scriptr = 0.5 * REdgeFlux2[0] * rptr[0];
        for (unsigned int j = 1; j < nedge - 1; ++j)
          scriptr += REdgeFlux2[j] * rptr[j];
        scriptr += 0.5 * REdgeFlux2[nedge-1] * rptr[nedge-1];
        RxEdgeWork[i] = scriptr;
      }
    } else {
      for (unsigned int i = 0; i < nedge; ++i) {
        rptr = rvals + i * nedge;
        scriptr = 0.5 * rptr[0];
        for (unsigned int j = 1; j < nedge-1; ++j)
          scriptr += rptr[j];
        scriptr += 0.5 * rptr[nedge - 1];
        RxEdgeWork[i] = scriptr;
      }
    }
    // Now X integral along RxEdgeWork, put in integral step size and area
    // of bin, store in r00val
    if (use_edge_log_x) {
      scriptr = 0.5 * REdgeFlux1[0] * RxEdgeWork[0];
      for (unsigned int i = 1; i < nedge - 1; ++i)
        scriptr += REdgeFlux1[i] * RxEdgeWork[i];
      scriptr += 0.5 * REdgeFlux1[nedge - 1] * RxEdgeWork[nedge - 1];
      r00val = scriptr * iR00norm;
    } else {
      scriptr = 0.5*RxEdgeWork[0];
      for (unsigned int i = 1; i < nedge-1; ++i)
        scriptr += RxEdgeWork[i];
      scriptr += 0.5 * RxEdgeWork[nedge - 1];
      r00val = scriptr * iR00norm;
    }

    // Now do Rx = R[0, y] (RxEdgeWork), integral along x
    //  First, compute R over the Edge x values for the real y, store in rvals
#ifdef TIMING
    starttime = std::clock();
#endif
    model.getR(nedge, REdgeFlux1, n, RFlux2, bm, rvals);
#ifdef TIMING
    RTime += std::clock() - starttime;
#endif

    // Now we integrate along each x and store the results
    // in RxEdgeWork.
    if (use_edge_log_x) {
      double rf;
      // Do i = 0
      rf = REdgeFlux1[0];
      rptr = rvals;
      for (unsigned int j = 0; j < n; ++j)
        RxEdgeWork[j] = 0.5 * rf * rptr[j];
      // i = 1 to nedge - 1
      for (unsigned int i = 1; i < nedge - 1; ++i) {
        rf = REdgeFlux1[i];
        rptr = rvals + i * n;
        for (unsigned int j = 0; j < n; ++j)
          RxEdgeWork[j] += rf * rptr[j];
      }
      // And i = nedge - 1
      rf = REdgeFlux1[nedge - 1];
      rptr = rvals + (nedge - 1) * n;
      for (unsigned int j = 0; j < n; ++j)
        RxEdgeWork[j] += 0.5 * rf * rptr[j];
    } else {
      for (unsigned int j = 0; j < n; ++j) // i = 0
        RxEdgeWork[j] = 0.5 * rvals[j];
      for (unsigned int i = 1; i < nedge - 1; ++i) {
        rptr = rvals + i * n;
        for (unsigned int j = 0; j < n; ++j)
          RxEdgeWork[j] += rptr[j];
      }
      // And i = nedge - 1
      rptr = rvals + (nedge - 1) * n;
      for (unsigned int j = 0; j < n; ++j)
        RxEdgeWork[j] += 0.5 * rptr[j];
    }    
    for (unsigned int j = 0; j < n; ++j)
      RxEdgeWork[j] *= iRxnorm;

    //And Ry = R[x, 0] (RyEdgeWork)
#ifdef TIMING
    starttime = std::clock();
#endif
    model.getR(n, RFlux1, nedge, REdgeFlux2, bm, rvals);
#ifdef TIMING
    RTime += std::clock() - starttime;
#endif    

    // Now we integrate along each x and store the results
    // in RyEdgeWork.
    if (use_edge_log_y) {
      for (unsigned int i = 0; i < n; ++i) {
        rptr = rvals + i * nedge;
        scriptr = 0.5 * REdgeFlux2[0] * rptr[0]; // j = 0
        for (unsigned int j = 1; j < nedge - 1; ++j)
          scriptr += REdgeFlux2[j] * rptr[j];
        scriptr += 0.5 * REdgeFlux2[nedge - 1] * rptr[nedge - 1]; // j = top
        RyEdgeWork[i] = scriptr;
      }
    } else {
      for (unsigned int i = 0; i < n; ++i) {
        rptr = rvals + i * nedge;
        scriptr = 0.5 * rptr[0]; // j = 0
        for (unsigned int j = 1; j < nedge - 1; ++j)
          scriptr += rptr[j];
        scriptr += 0.5 * rptr[nedge - 1]; // j = top
        RyEdgeWork[i] = scriptr;
      }
    }
    for (unsigned int i = 0; i < n; ++i)
      RyEdgeWork[i] *= iRynorm;
  }

  // Fill main array.  This will set the edges to nonsense
#ifdef TIMING
  std::clock_t starttime = std::clock();
#endif
  model.getR(n, RFlux1, n, RFlux2, bm, rvals);
#ifdef TIMING
  RTime += std::clock() - starttime;
#endif
  
  if (setEdge) {
    // Copy over from REdgeWork arrays
    std::memcpy(rvals + 1, RxEdgeWork + 1, (n - 1) * sizeof(double));
    for (unsigned int i = 1; i < n; ++i)
      rvals[i * n] = RyEdgeWork[i];
    rvals[0] = r00val;
  } else {
    // Set edges to zero
    std::memset(rvals, 0, n * sizeof(double)); // includes 0, 0
    for (unsigned int i = 1; i < n; ++i)
      rvals[i * n] = 0.0;
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
  // takes rvals to prtrans.  We then do things to prtrans
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
  if (plan == nullptr) {
    plan = fftw_plan_dft_r2c_2d(intn, intn, rvals, prtrans,
                                fftw_plan_style);
  }
  if (plan == nullptr) {
    std::stringstream str;
    str << "Plan creation failed for forward transform of size: " << n;
    if (has_wisdom) str << std::endl << "Your wisdom file may not have"
                        << " that size";
    throw affineExcept("PDFactoryDouble", "setupTransforms", str.str());
  }
  if (plan_inv == nullptr) {
    plan_inv = fftw_plan_dft_c2r_2d(intn, intn, pval, pofd,
                                    fftw_plan_style);
    if (plan_inv == nullptr) {
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
  \param[in] data Marginalized P(D) to find split point at.  Must be 
  purly non-negative.
  \param[in] dflux Flux step
  \param[in] sigmapos Estimated sigma for marginalized data for positive
  flux densities
  \param[in] sigmaneg Estimated sigma for marginalized data for negative
  flux densities
  \returns Split point index

  This is a utility function for finding the split point to 'unwrap' a P(D).
  The user will call this once for each dimension, and is responsible
  for doing the marginalization.
*/
unsigned int PDFactoryDouble::findSplitPoint(const double* const data,
                                             double dflux, double sigmapos,
                                             double sigmaneg) const {
  // This is basically the same thing as PDFactory::findSplitPoint,
  //  but because some of the internal details of how things
  //  are stored differ, it is implemented both here and in PDFactory.
  //  See that for a much more detailed explanation of what this
  //  is trying to do.

  // These are are control constants for accepting points
  //  We use slightly more forgiving values here than in PDFactory
  //  because 2D is hard, and not as well sampled in practice.
  const double n1 = 2.25;     // Minimum number of sigma away from mean
  const double f1 = 5e-14;   // Jitter estimate #1 factor
  const double f2 = 1.5e-5;  // Peak relative estimate
  const double f3 = 5;       // Jitter estimate #2

  if (sigmapos <= 0)
    throw affineExcept("PDFactory", "findSplitPoint", 
                       "Invalid (non-positive) sigma pos");
  if (sigmaneg <= 0)
    throw affineExcept("PDFactory", "findSplitPoint", 
                       "Invalid (non-positive) sigma neg");

  double fluxrange = currsize * dflux;
  double n1sigpos = n1 * sigmapos;
  double n1signeg = n1 * sigmaneg;

  // Step 1
  if (fluxrange < (n1sigpos + n1signeg)) {
    std::stringstream errstr;
    errstr << "Not enough flux range; range is: "
           << fluxrange << " sigmas are: " << sigmapos
           << " and " << sigmaneg
           << " so required flux range is: "
           << (sigmapos + sigmaneg);
    throw affineExcept("PDFactoryDouble", "findSplitPoint", 
                       errstr.str());
  }

  // Step 2
  double peak, currval;
  peak = data[0];
  for (unsigned int i = 1; i < currsize; ++i) {
    currval = data[i];
    if (currval > peak) peak = currval;
  }
  if (peak <= 0.0)
    throw affineExcept("PDFactoryDouble", "findSplitPoint",
                       "Peak P(D) value was zero");

  // Step 3
  // Figure out the search index range from -n1signeg to +n1sigpos
  int topidx = currsize - static_cast<int>(n1signeg / dflux) - 1;
  int botidx = static_cast<int>(n1sigpos / dflux);
  if (topidx >= currsize || botidx < 0 || botidx > topidx) {
    std::stringstream errstr;
    errstr << "Logic error in topidx/botidx with values "
           << topidx << " and " << botidx << " and size "
           << currsize << " with dflux " << dflux
           << " n1signeg " << n1signeg << " and n1sigpos "
           << n1sigpos;
    throw affineExcept("PDFactoryDouble", "findSplitPoint", 
                       errstr.str());   
  }
  // Find minimum value
  int minidx;
  double minval;
  minidx = botidx;
  minval = data[minidx];
  for (int i = botidx + 1; i <= topidx; ++i) {
    currval = data[i];
    if (currval < minval) {
      minidx = i;
      minval = currval;
    }
  }

  // Step 4
  double targval = f2 * peak;
  if (minval > targval) {
    std::stringstream errstr;
    errstr << "Minimum in marginalized pd is too large "
           << " relative to peak; can't find split point "
           << "with minimum: " << minval << " peak value: " 
           << peak << " Required split value < " << targval
           << " in index range " << botidx << " to " << topidx
           << " in P(D) of size: " << currsize;
    throw affineExcept("PDFactoryDouble", "findSplitPoint", 
                       errstr.str());
  }

  // Step 5
  const int pre_jitter_n = 17; // Elements before jitter to check
  double jitter = 0.0; // Estimate
  if (minidx > pre_jitter_n) {
    // Check if these values are decreasing
    bool monotonic_decrease = true;
    int initidx = minidx - pre_jitter_n;

    double currval, prevval;
    jitter = prevval = data[initidx];
    for (int i = initidx + 1; i < minidx; ++i) {
      currval = data[i];
      if (currval >= prevval)
        monotonic_decrease = false;
      jitter += currval;
      prevval = currval;
    }
    if (monotonic_decrease) // Success
      return static_cast<unsigned int>(minidx);
    else {
      // Seems to be in the jitter; use 1.1x the mean as a jitter estimator
      // estimate
      jitter += minval;
      jitter = 1.1 * jitter / static_cast<double>(pre_jitter_n + 1);
    }
  } else {
    // We don't have room!  This is odd, but possible
    // Just form an estimate from f1
    jitter = f1 * peak;
    if (minval > jitter) // Success!
      return static_cast<unsigned int>(minidx);
  }

  // Step 6
  // Update jitter estimate.  Beware of minimum being zero!
  if (minval <= 0.0) {
    // Try to use some neighboring points to get nicer estimate
    const int nminup = 3;
    int bi, ti;
    if (minidx < botidx) bi = botidx; 
    else bi = (minidx >= nminup) ? (minidx - nminup) : minidx;
    if (minidx > topidx - nminup ) ti = topidx - nminup; 
    else ti = minidx + nminup;
    ti = (ti >= currsize) ? (currsize - 1) : ti;
    for (int i = bi; i < ti; ++i) {
      currval = data[i];
      if (currval > minval) minval = currval;
    }
  }
  if (minval > 0)
    jitter = (f3 * minval) < jitter ? (f3 * minval) : jitter;

  // Step 7
  // Make sure there's a range to search
  if ((botidx >= topidx) || (abs(topidx - botidx) < 1)) {
    std::stringstream errstr;
    errstr << "Topidx (" << topidx << ") vs. botidx ("
           << botidx << ") search issue -- can't find range"
           << " to search";
    throw affineExcept("PDFactoryDouble", "findSplitPoint", 
                       errstr.str());
  }
  // Recall we search down.  We want to find the first
  //  point below the jitter estimate and return the point above
  //  that (which must be above the jitter value).  First make
  //  sure that's true for the first point.
  if (data[topidx] < jitter) {
    std::stringstream errstr;
    errstr << "Topidx (" << topidx << ") is already below jitter"
           << " estimate (" << jitter << ") with value "
           << data[topidx] << " minval: " << minval
           << " and peak: " << peak << "; can't find split";
    throw affineExcept("PDFactoryDouble", "findSplitPoint", 
                       errstr.str());
  }
  // Now search
  for (int i = topidx - 1; i >= botidx; --i)
    if (data[i] < jitter) return static_cast<unsigned int>(i + 1); // Success

  // Step 8
  // Didn't find one -- throw exception.  But first find minimum
  //  so that error message is more useful
  minval = data[topidx];
  for (int i = topidx - 1; i >= botidx; --i) 
    if (data[i] < minval) minval = data[i];
  std::stringstream errstr;
  errstr << "Unable to find split point in range "
         << topidx << " down to " << botidx
         << " with value less than " << jitter
         << "; minimum value found: " << minval;
  throw affineExcept("PDFactoryDouble", "findSplitPoint", 
                     errstr.str());
}

/*!
  \param[out] pd Holds P(D) on output, normalized, mean subtracted,
  and with positivity enforced.

  This does the unwrapping of the internal P(D) into pd.
*/
// This should only ever be called by getPD, so we don't really
// check the inputs
void PDFactoryDouble::unwrapAndNormalizePD(PDDouble& pd) const {

  //Enforce positivity
#ifdef TIMING
  starttime = std::clock();
#endif
  for (unsigned int idx = 0; idx < currsize * currsize; ++idx)
    if (pofd[idx] < 0) pofd[idx] = 0.0;
#ifdef TIMING
  posTime += std::clock() - starttime;
#endif

  // Figure out the indices of the split points
  // Do index 1 first
#ifdef TIMING
  starttime = std::clock();
#endif
  // Sum into rsum, making this a 1D problem
  double cval, *rowptr;
  for (unsigned int i = 0; i < currsize; ++i) {
    rowptr = pofd + i * currsize;
    cval = 0.5 * rowptr[0];
    for (unsigned int j = 1; j < currsize - 1; ++j)
      cval += rowptr[j];
    cval += 0.5 * rowptr[currsize-1];
    rsum[i] = cval;
  }

  unsigned int splitidx1;
  try {
    splitidx1 = findSplitPoint(rsum, dflux1, sg1pos, sg1neg);
  } catch (const affineExcept& ex) {
    // Add information to the exception about which band we were in.
    std::stringstream newerrstr("");
    newerrstr << "In dimension one, problem with split point: " 
              << ex.getErrStr();
    throw affineExcept(ex.getErrClass(), ex.getErrMethod(), 
                       newerrstr.str());
  }

  // Now second dimension.  This one is slower due to stride issues,
  //  but otherwise is the same. Doing the sum in this order is much faster
  for (unsigned int j = 0; j < currsize; ++j) rsum[j] = 0.5 * pofd[j]; // i=0
  for (unsigned int i = 1; i < currsize - 1; ++i) { // center
    rowptr = pofd + i * currsize;
    for (unsigned int j = 0; j < currsize; ++j)
      rsum[j] += rowptr[j];
  }
  rowptr = pofd + (currsize - 1) * currsize; // i = currsize-1
  for (unsigned int j = 0; j < currsize; ++j)
    rsum[j] += 0.5 * rowptr[j];

  unsigned int splitidx2;
  try {
    splitidx2 = findSplitPoint(rsum, dflux2, sg2pos, sg2neg);
  } catch (const affineExcept& ex) {
    std::stringstream newerrstr("");
    newerrstr << "In dimension two, problem with split point: " 
              << ex.getErrStr();
    throw affineExcept(ex.getErrClass(), ex.getErrMethod(), 
                       newerrstr.str());
  }
#ifdef TIMING
  splitTime += std::clock() - starttime;
#endif

  // Now the actual copying, which is an exercise in index gymnastics
#ifdef TIMING
  starttime = std::clock();
#endif
  pd.resize(currsize, currsize);
  size_t size_double = sizeof(double); // Size of double
  std::memset(pd.pd_, 0, currsize * currsize * size_double);

  size_t rowsz;
  double *ptr_curr, *ptr_out; // convenience pointers
  // We start with the neg, neg bit -- that is stuff >= both
  // splitidx1 and splitidx2, which ends up in the 0, 0 part of the output
  ptr_out = pd.pd_;
  ptr_curr = pofd + splitidx1 * currsize + splitidx2;
  rowsz = (currsize - splitidx2) * size_double;
  for (unsigned int i = 0; i < currsize - splitidx1; ++i)
    std::memcpy(ptr_out + i * currsize, ptr_curr + i * currsize, rowsz);

  // Next pos neg
  ptr_out = pd.pd_ + (currsize - splitidx1) * currsize;
  ptr_curr = pofd + splitidx2;
  for (unsigned int i = 0; i < splitidx1; ++i) 
    std::memcpy(ptr_out + i * currsize, ptr_curr + i * currsize, rowsz);

  // pos, pos
  ptr_out = pd.pd_ + (currsize - splitidx1) * currsize + (currsize - splitidx2);
  ptr_curr = pofd;
  rowsz = splitidx2 * size_double;
  for (unsigned int i = 0; i < splitidx1; ++i)
    std::memcpy(ptr_out + i * currsize, ptr_curr + i * currsize, rowsz);
  
  // neg, pos
  ptr_out = pd.pd_ + (currsize - splitidx2);
  ptr_curr = pofd + splitidx1 * currsize;
  for (unsigned int i = 0; i < currsize - splitidx1; ++i)
    std::memcpy(ptr_out + i * currsize, ptr_curr + i * currsize, rowsz);
#ifdef TIMING
  copyTime += std::clock() - starttime;
#endif

  pd.logflat = false;
  pd.minflux1 = 0.0; pd.dflux1 = dflux1;
  pd.minflux2 = 0.0; pd.dflux2 = dflux2;

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
  dblpair tmn; // True mean from P(D) rather than estimate
  tmn = pd.getMeans(false);
  if (std::isinf(tmn.first) || std::isnan(tmn.first) || 
      std::isinf(tmn.second) || std::isnan(tmn.second)) {
    std::stringstream str;
    str << "Un-shift amounts not finite band1: " << tmn.first << " band2: "
        << tmn.second << std::endl;
    str << "At length: " << currsize << " with sigmas: " << sg1pos << " "
        << sg1neg << " and " << sg2pos << " " << sg2neg;
    throw affineExcept("PDFactoryDouble", "unwrapAndNormalizePD", 
                       str.str());
  }
  pd.minflux1 = -tmn.first;
  pd.minflux2 = -tmn.second;
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

  After this, prtrans holds exp(rtrans - rtrans[0] - shift*omega)
  where rtrans is the FFTed R and the shift is optional.  Then getPD
  only needs to add the sigma bits to get p.
*/
bool PDFactoryDouble::initPD(unsigned int n, double minflux1, double maxflux1, 
                             double minflux2, double maxflux2, 
                             const numberCountsDouble& model,
                             const doublebeam& bm, bool setEdge) {
  if (n == 0)
    throw affineExcept("PDFactoryDouble", "initPD", "Invalid (zero) n");
  if (n == 1)
    throw affineExcept("PDFactoryDouble", "initPD", "Invalid n (==1)");
  if (maxflux1 <= 0.0) {
    std::stringstream errstr;
    errstr << "Invalid (non-positive) maxflux1 (" << maxflux1 << ")";
    throw affineExcept("PDFactoryDouble", "initPD", errstr.str());
  }
  if (maxflux2 <= 0.0) {
    std::stringstream errstr;
    errstr << "Invalid (non-positive) maxflux2 (" << maxflux2 << ")";
    throw affineExcept("PDFactoryDouble", "initPD", errstr.str());
  }
  if (maxflux1 == minflux1) {
    std::stringstream errstr;
    errstr << "Invalid range (0) between min/max flux in band 1 with value: "
           << maxflux1;
    throw affineExcept("PDFactoryDouble", "initPD", errstr.str());
  }
  if (maxflux2 == minflux2) {
    std::stringstream errstr;
    errstr << "Invalid range (0) between min/max flux in band 2 with value: "
           << maxflux2;
    throw affineExcept("PDFactoryDouble", "initPD", errstr.str());
  }
  
  initialized = false;

  //Make the plans, or keep the old ones if possible
  // Note we have to do this before we fill R, as plan construction
  // may overwrite the values.  This allocates R and resizes
  setupTransforms(n);

  // Compute R, including multiplying it by dflux
  initR(n, minflux1, maxflux1, minflux2, maxflux2, model, bm, 
        setEdge, true);

  // Get mean and variance estimates from R
  getRStats();

  // Set the sigmas to the no-noise values for now
  sg1pos = sqrt(var_noi1);
  sg2pos = sqrt(var_noi2);

  //Decide if we will shift and pad, and if so by how much.
  // The idea is to shift the mean to zero -- but we only
  // do the shift if the sigma is larger than one actual step size
  // because otherwise we can't represent it well.
  doshift1 = (sg1pos > dflux1) && (fabs(mn1) > dflux1);
  if (doshift1) shift1 = -mn1; else shift1 = 0.0;
  doshift2 = (sg2pos > dflux2) && (fabs(mn2) > dflux2);
  if (doshift2) shift2 = -mn2; else shift2 = 0.0;

  if (verbose) {
    std::cout << " Initial mean estimate band1: " << mn1 << " band2: "
              << mn2 << std::endl;
    std::cout << " Initial stdev estimate band1 (no inst noise): " 
              << sg1pos << " band2: "
              << sg2pos << std::endl;
    if (doshift1)
      std::cout << " Additional shift in band 1: " << shift1 << std::endl;
    else std::cout << " Not applying additional shift in band 1" << std::endl;
    if (doshift2)
      std::cout << " Additional shift in band 2: " << shift2 << std::endl;
    else std::cout << " Not applying additional shift in band 2" << std::endl;
  }

  //Compute forward transform of this r value, store in prtrans
  //Have to use argument version, since the address of prtrans, etc. can move
#ifdef TIMING
  starttime = std::clock();
#endif
  fftw_execute_dft_r2c(plan, rvals, prtrans);
#ifdef TIMING
  fftTime += std::clock() - starttime;
#endif
  initialized = true;

  // Now replace FFT(R) with most of p in prtrans in place.
  // The 2D output real FFT format makes this a bit tricky.  The output
  //  array is n by (n/2+1).  The first dimension has both
  //  positive and negative frequencies, the second dimension has
  //  only positive dimensions.
  // In particular, if i is the index over the first dimension,
  //  the frequencies along the first dimension are:
  //   f1 = i/delta1*n        for i = [0,n/2]
  //      = - (n-i)/delta1*n  for i = [n/2+1,n-1]
  //  where delta1 is dflux1.
  // And if j is the second dimension index, then
  //   f2 = j/delta2*n  
  //  and delta2 = dflux2.
  // We work in w instead of f (2 pi f) and actually compute 
  //   exp(r(omega1,omega2) - r(0,0) - i*shift1*omega1 - i*shift2*omega2)
  //  where r is FFT(R)
  unsigned int ncplx = n / 2 + 1;
  double r0, expfac, rval, ival;
  fftw_complex *rowptr; // Row pointer into prtrans
  double iflux1 = mcmc_affine::two_pi / (n * dflux1);
  double iflux2 = mcmc_affine::two_pi / (n * dflux2);

#ifdef TIMING
  std::clock_t starttime = std::clock();
#endif

  r0 = prtrans[0][0]; // r[0,0] is pure real

  if (doshift1 && doshift2) {
    // Should be the most common case
    double shiftfac1 = shift1 * iflux1;
    double shiftfac2 = shift2 * iflux2;
    double shiftprod1;
    // First, positive frequencies
    for (unsigned int idx1 = 0; idx1 < ncplx; ++idx1) {
      rowptr = prtrans + idx1 * ncplx;
      shiftprod1 = shiftfac1 * static_cast<double>(idx1);
      for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
        rval = rowptr[idx2][0] - r0;
        ival = rowptr[idx2][1] - shiftprod1 -
          shiftfac2 * static_cast<double>(idx2);
        expfac = exp(rval);
        rowptr[idx2][0] = expfac * cos(ival);
        rowptr[idx2][1] = expfac * sin(ival);
      }
    }
    // Now negative frequencies
    for (unsigned int idx1 = ncplx; idx1 < n; ++idx1) {
      rowptr = prtrans + idx1 * ncplx;
      shiftprod1 = - shiftfac1 * static_cast<double>(n - idx1);
      for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
        rval = rowptr[idx2][0] - r0;
        ival = rowptr[idx2][1] - shiftprod1 -
          shiftfac2 * static_cast<double>(idx2);
        expfac = exp(rval);
        rowptr[idx2][0] = expfac * cos(ival);
        rowptr[idx2][1] = expfac * sin(ival);
      }
    }
  } else if (doshift1) {
    // Only a shift in band 1.  Probably pretty uncommon
    double shiftfac1 = shift1 * iflux1;
    double shiftprod1;
    for (unsigned int idx1 = 0; idx1 < ncplx; ++idx1) {
      rowptr = prtrans + idx1 * ncplx;
      shiftprod1 = shiftfac1 * static_cast<double>(idx1);
      for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
        rval = rowptr[idx2][0] - r0;
        ival = rowptr[idx2][1] - shiftprod1;
        expfac = exp(rval);
        rowptr[idx2][0] = expfac * cos(ival);
        rowptr[idx2][1] = expfac * sin(ival);
      }
    }
    for (unsigned int idx1 = ncplx; idx1 < n; ++idx1) {
      rowptr = prtrans + idx1 * ncplx;
      shiftprod1 = - shiftfac1 * static_cast<double>(n - idx1);
      for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
        rval = rowptr[idx2][0] - r0;
        ival = rowptr[idx2][1] - shiftprod1;
        expfac = exp(rval);
        rowptr[idx2][0] = expfac * cos(ival);
        rowptr[idx2][1] = expfac * sin(ival);
      }
    }
  } else if (doshift2) {
    // Only shift in band 2 -- also probably pretty uncommon
    double shiftfac2 = shift2 * iflux2;
    // Note we no longer have to split into neg/pos frequency bits
    for (unsigned int idx1 = 0; idx1 < n; ++idx1) {
      rowptr = prtrans + idx1 * ncplx;
      for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
        rval = rowptr[idx2][0] - r0;
        ival = rowptr[idx2][1] - shiftfac2 * static_cast<double>(idx2);
        expfac = exp(rval);
        rowptr[idx2][0] = expfac * cos(ival);
        rowptr[idx2][1] = expfac * sin(ival);
      }
    }
  } else {
    // No shifts -- very easy, more common than the previous 2,
    //  although only in testing
    for (unsigned int idx1 = 0; idx1 < n; ++idx1) {
      rowptr = prtrans + idx1 * ncplx;
      for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
        rval = rowptr[idx2][0] - r0;
        ival = rowptr[idx2][1];
        expfac = exp(rval);
        rowptr[idx2][0] = expfac * cos(ival);
        rowptr[idx2][1] = expfac * sin(ival);
      }
    }
  }
  prtrans[0][0] = 1.0;
  prtrans[0][1] = 0.0;

#ifdef TIMING
  p0Time += std::clock() - starttime;
#endif
  
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

  // If set, we assume that the P(D) will tend to cut off much more
  // sharply at negative flux densities, which arises when the beam
  // is mostly positive.  This is used in the split point computation,
  // and is a good approximation for Herschel data but perhaps not in
  // general.
  const bool assume_negative_sharp = true;

  // The basic idea is to compute the P(D) from the previously filled
  // R values, adding in noise filling pd for output.
  // Note that pval will not be filled if the sigmas are zero
  if (!initialized )
    throw affineExcept("PDFactoryDouble", "getPD",
                       "Must call initPD first");
  if (!rvars_allocated) 
    throw affineExcept("PDFactoryDouble", "getPD",
                       "Should not have been able to get here without allocating R variables");

  //Output array from 2D FFT is n * (n/2+1)
  unsigned int n = currsize;
  unsigned int ncplx = n / 2 + 1;

  // Finish computing P and transform.
  // For simplicity, assume if there is
  // noise in either band there will be noise in both, which should
  // always be true in practice.
  if (sigma1 > 0 || sigma2 > 0) {
    // Real data case (that is, with noise)

#ifdef TIMING
    std::clock_t starttime = std::clock();
#endif

    fftw_complex *r_input_rowptr;
    fftw_complex *row_current_out;
    
    double iflux1 = mcmc_affine::two_pi / (n * dflux1);
    double iflux2 = mcmc_affine::two_pi / (n * dflux2);
    double sigfac1 = -0.5 * sigma1 * sigma1 * iflux1 * iflux1;
    double sigfac2 = -0.5 * sigma2 * sigma2 * iflux2 * iflux2;
    double expfac, sigprod1, didx1, didx2;
    //First, Pos freq
    for (unsigned int idx1 = 0; idx1 < ncplx; ++idx1) {
      r_input_rowptr = prtrans + idx1 * ncplx;
      row_current_out = pval + idx1 * ncplx;
      didx1 = static_cast<double>(idx1);
      sigprod1  = sigfac1 * didx1 * didx1; // -1/2 sigma1^2 w1^2
      for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
        didx2 = static_cast<double>(idx2);
        //  exp(-1/2 sigma1^2 w1^2 - 1/2 sigma2^2 w2^2)
        expfac = exp(sigprod1 + sigfac2 * didx2 * didx2);
        row_current_out[idx2][0] = expfac * r_input_rowptr[idx2][0];
        row_current_out[idx2][1] = expfac * r_input_rowptr[idx2][1];
      }
    }
    // Negative frequency bit; needs to be a separate loop because
    //  the direction that w1 increases flips
    for (unsigned int idx1 = ncplx; idx1 < n; ++idx1) {
      r_input_rowptr = prtrans + idx1 * ncplx;
      row_current_out = pval + idx1 * ncplx;
      didx1 = static_cast<double>(n - idx1); // Drop - sgn since we only ^2
      sigprod1 = sigfac1 * didx1 * didx1; // -1/2 sigma1^2 w1^2
      for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
        didx2 = static_cast<double>(idx2);
        //  exp(-1/2 sigma1^2 w1^2 - 1/2 sigma2^2 w2^2)
        expfac = exp(sigprod1 + sigfac2 * didx2 * didx2);
        row_current_out[idx2][0] = expfac * r_input_rowptr[idx2][0];
        row_current_out[idx2][1] = expfac * r_input_rowptr[idx2][1];
      }
    }
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
    fftw_execute_dft_c2r(plan_inv, pval, pofd);
#ifdef TIMING
    fftTime += std::clock() - starttime;
#endif
  } else {
    // Noiseless, probably a test case.  We can just directly
    // transform prtrans now.
    if (verbose) std::cout << " Reverse transform" << std::endl;
#ifdef TIMING
    starttime = std::clock();
#endif
    // From pval into pofd
    fftw_execute_dft_c2r(plan_inv, prtrans, pofd);
#ifdef TIMING
    fftTime += std::clock() - starttime;
#endif
  }
  
  // Copy into output variable, also normalizing, mean subtracting, 
  // making positive.  Need sigma estimates for unwrapAndNormalize
  // As is the case for PDFactory, we are going to assume that there 
  // isn't much negative beam so that the sigma in the negative direction
  // tends to be dominated by instrument noise.
  sg1pos = sqrt(var_noi1 + sigma1 * sigma1);
  sg2pos = sqrt(var_noi2 + sigma2 * sigma2);
  if (assume_negative_sharp) {
    if (9 * sigma1 * sigma1 > var_noi1) sg1neg = sigma1; else sg1neg = sg1pos;
    if (9 * sigma2 * sigma2 > var_noi2) sg2neg = sigma2; else sg2neg = sg2pos;
  } else {
    sg1neg = sg1pos;
    sg2neg = sg2pos;
  }
  
  try {
    unwrapAndNormalizePD(pd);
  } catch  (const affineExcept& ex) {
    // Add information to the exception about actual sigma values
    std::stringstream newerrstr("");
    newerrstr << ex.getErrStr() << std::endl
              << " for sigma values: " << sigma1 << " " << sigma2;
    throw affineExcept(ex.getErrClass(), ex.getErrMethod(), 
                       newerrstr.str());
  }

  //Turn PD to log for more efficient log computation of log likelihood
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

  You must call initR or initPD first.
*/
void PDFactoryDouble::writeRToHDF5(const std::string& filename) const {
  if (!rinitialized)
    throw affineExcept("PDFactoryDouble", "writeRToHDF5",
                       "Must call initPD or initR first");

  hid_t file_id;
  file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                      H5P_DEFAULT);

  if (H5Iget_ref(file_id) < 0) {
    H5Fclose(file_id);
    throw affineExcept("PDFactoryDouble", "writeToHDF5",
                       "Failed to open HDF5 file to write");
  }

  // Write it as one dataset -- Rflux1, Rflux2, R. 
  hsize_t adims;
  hid_t mems_id, att_id;
  
  // First, some properties
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, nullptr);
  att_id = H5Acreate2(file_id, "dFlux1", H5T_NATIVE_DOUBLE,
                      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dflux1);
  H5Aclose(att_id);
  att_id = H5Acreate2(file_id, "dFlux2", H5T_NATIVE_DOUBLE,
                      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dflux2);
  H5Aclose(att_id);
  H5Sclose(mems_id);

  // Rfluxes
  hdf5utils::writeDataDoubles(file_id, "RFlux1", currsize, RFlux1);
  hdf5utils::writeDataDoubles(file_id, "RFlux2", currsize, RFlux2);

  // R -- which we may need to copy to remove the dflux
  if (rdflux) {
    double* tmp = new double[currsize * currsize];
    double idflux = 1.0 / (dflux1 * dflux2);
    for (unsigned int i = 0; i < currsize * currsize; ++i) 
      tmp[i] = rvals[i] * idflux;
    hdf5utils::writeData2DDoubles(file_id, "R", currsize, currsize, tmp);
    delete[] tmp;
  } else
    hdf5utils::writeData2DDoubles(file_id, "R", currsize, currsize, rvals);

  // Done
  H5Fclose(file_id);
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
void PDFactoryDouble::receiveCopy(MPI_Comm comm, int src) {
  MPI_Status Info;
  MPI_Recv(&fftw_plan_style, 1, MPI_UNSIGNED, src, 
           pofd_mcmc::PDFDSENDPLANSTYLE, comm, &Info);
  MPI_Recv(&has_wisdom, 1, MPI::BOOL, src, pofd_mcmc::PDFDHASWISDOM,
           comm, &Info);
  if (has_wisdom) {
    //Receive wisdom file name
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
