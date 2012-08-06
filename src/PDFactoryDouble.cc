#include<sstream>
#include<cmath>
#include<cstring>
#include<limits>
#include<iomanip>

#include<global_settings.h>
#include<PDFactoryDouble.h>
#include<affineExcept.h>

const double PDFactoryDouble::lowEdgeRMult=1e-9;
//Control of how we do the edge integrals -- linear or log?
const bool PDFactoryDouble::use_edge_log_x = true;
const bool PDFactoryDouble::use_edge_log_y = false;

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
  lastfftlen = 0;
  currsize = 0;

#ifdef TIMING
  resetTime();
#endif

  rvars_allocated = false;
  RFlux1 = RFlux2 = NULL;
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
  fftw_plan_style = FFTW_ESTIMATE;

  max_sigma1 = max_sigma2 = std::numeric_limits<double>::quiet_NaN();
  mn1 = mn2 = var_noi1 = var_noi2 = sg1 = sg2 = 
    std::numeric_limits<double>::quiet_NaN();
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
  RTime = RStatsTime = p0Time = fftTime = posTime = copyTime = normTime = 0;
  edgeTime = meanTime = logTime = 0;
}

void PDFactoryDouble::summarizeTime(unsigned int nindent) const {
  std::string prestring(nindent,' ');
    
  std::cout << "R time: " << prestring 
	    << 1.0*RTime/CLOCKS_PER_SEC << std::endl;
  std::cout << "R stats time: " << prestring 
	    << 1.0*RStatsTime/CLOCKS_PER_SEC << std::endl;
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
  Doesn't force resize if there is already enough room available

  \param[in] NSIZE new size
  \returns True if a resize was needed
 */
bool PDFactoryDouble::resize(unsigned int NSIZE) {
  if (NSIZE > currsize) {
    strict_resize(NSIZE);
    return true;
  } else return false;
}

/*!
  Forces a resize

  \param[in] NSIZE new size (must be > 0)
 */
//I don't check for 0 NSIZE because that should never happen
// in practice, and it isn't worth the idiot-proofing
// rtrans always gets nulled in here, then refilled when you
//  call initPD
void PDFactoryDouble::strict_resize(unsigned int NSIZE) {
  if (NSIZE == currsize) return;
  freeRvars();
  currsize = NSIZE;
  initialized = false;
}

void PDFactoryDouble::allocateRvars() {
  if (rvars_allocated) return;
  if (currsize == 0)
    throw affineExcept("PDFactory","allocate_rvars",
		       "Invalid (0) currsize",1);
  RFlux1 = (double*) fftw_malloc(sizeof(double)*currsize);
  RFlux2 = (double*) fftw_malloc(sizeof(double)*currsize);
  unsigned int fsize = currsize*currsize;
  rvals = (double*) fftw_malloc(sizeof(double)*fsize);
  pofd  = (double*) fftw_malloc(sizeof(double)*fsize);
  fsize = currsize*(currsize/2+1);
  rtrans = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*fsize);
  pval = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*fsize);

  rvars_allocated = true;
  initialized = false;
}

void PDFactoryDouble::freeRvars() {
  if (RFlux1 != NULL)     { fftw_free(RFlux1); RFlux1=NULL; }
  if (RFlux2 != NULL)     { fftw_free(RFlux2); RFlux2=NULL; }
  if (rvals != NULL)      { fftw_free(rvals); rvals=NULL; }
  if (rtrans != NULL)     { fftw_free(rtrans); rtrans=NULL; }
  if (pval != NULL)       { fftw_free(pval); pval = NULL; }
  if (pofd != NULL)       { fftw_free(pofd); pofd = NULL; }
  rvars_allocated = false;
  initialized = false;
}

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
  initialized = false;
}

void PDFactoryDouble::freeEdgevars() {
  if (REdgeFlux1 != NULL) { fftw_free(REdgeFlux1); REdgeFlux1 = NULL; }
  if (REdgeFlux2 != NULL) { fftw_free(REdgeFlux2); REdgeFlux2 = NULL; }
  if (REdgeWork != NULL)  { fftw_free(REdgeWork); REdgeWork = NULL; }
  edgevars_allocated = false;
  initialized = false;
}

/*!
  \param[in] filename Name of wisdom file
*/
bool PDFactoryDouble::addWisdom(const std::string& filename) {
  FILE *fp = NULL;
  fp = fopen( filename.c_str(), "r" );
  if (!fp) {
    std::stringstream str;
    str << "Error opening wisdom file: " << filename;
    throw affineExcept("PDFactoryDouble","addWisdom",str.str(),1);
  }
  if (fftw_import_wisdom_from_file(fp) == 0) {
    std::stringstream str;
    str << "Error reading wisdom file: " << wisdom_file;
    throw affineExcept("PDFactoryDouble","addWisdom",str.str(),1);
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

  return true;
}

/*!
  This assumes R, Rflux, etc. is set and doesn't check, so be careful
  It also ignores the edges.  The instrument noise is not included
  in the returned variances.
*/
void PDFactoryDouble::getRStats(unsigned int n, double& mn1, double& mn2,
				double& var1, double& var2) const {

  //  <x> = \int x R dx dy
  //  <y> = \int y R dx dy
  //  Var[x] = \int x^2 R dx dy
  //  Var[y] = \int y^2 R dx dy

  //Use trap rule for all computations
  //Always move along j since array access is
  // faster that way (rvals is row-major order)
  //This is a simple calculation, somewhat tedious to write out
  double cf1, cf2, crv;
  register double tmp;
  double *rowptr;

  //i=0
  cf1   = RFlux1[0];
  cf2   = RFlux2[0];
  crv   = rvals[0];
  mn1   = 0.5*cf1*crv;
  mn2   = 0.5*cf2*crv;
  var1  = mn1*cf1;
  var2  = mn2*cf2;
  for (unsigned int j = 1; j < n-1; ++j) {
    crv    = rvals[j];
    cf2    = RFlux2[j];
    tmp    = cf1*crv;
    mn1   += tmp;
    var1  += tmp*cf1;
    tmp    = cf2*crv;
    mn2   += tmp;
    var2  += tmp*cf2;
  }
  cf2 = RFlux2[n-1];
  crv = rvals[n-1];
  tmp    = 0.5*cf1*crv;
  mn1   += tmp;
  var1  += tmp*cf1;
  tmp    = 0.5*cf2*crv;
  mn2   += tmp;
  var2  += tmp*cf2;

  //That was the i=0 case, which is the first step and hence must
  // be multiplied by 0.5
  mn1  *= 0.5;
  mn2  *= 0.5;
  var1 *= 0.5;
  var2 *= 0.5;

  //Now main body, multiplier is one
  for (unsigned int i = 1; i < n-1; ++i) {
    cf1 = RFlux1[i];
    cf2 = RFlux2[0];
    rowptr = rvals + i*n;
    crv    = rowptr[0];
    tmp    = 0.5*cf1*crv;
    mn1   += tmp;
    var1  += tmp*cf1;
    tmp    = 0.5*cf2*crv;
    mn2   += tmp;
    var2  += tmp*cf2;
    for (unsigned int j = 1; j < n-1; ++j) {
      crv = rowptr[j];
      cf2 = RFlux2[j];
      tmp    = cf1*crv;
      mn1   += tmp;
      var1  += tmp*cf1;
      tmp    = cf2*crv;
      mn2   += tmp;
      var2  += tmp*cf2;
    }
    cf2 = RFlux2[n-1];
    crv = rowptr[n-1];
    tmp    = 0.5*cf1*crv;
    mn1   += tmp;
    var1  += tmp*cf1;
    tmp    = 0.5*cf2*crv;
    mn2   += tmp;
    var2  += tmp*cf2;
  }

  //i=n-1, 0.5 multiplier, store in tmp vars
  double tmn1, tmn2, tvar1, tvar2;
  cf1 = RFlux1[n-1];
  cf2 = RFlux2[0];
  rowptr = rvals + (n-1)*n;
  crv    = rowptr[0];
  tmn1   = 0.5*cf1*crv;
  tvar1  = tmn1 * cf1;
  tmn2   = 0.5*cf2*crv;
  tvar2  = tmn2 * cf2;
  for (unsigned int j = 1; j < n-1; ++j) {
    crv     = rowptr[j];
    cf2     = RFlux2[j];
    tmp     = cf1*crv;
    tmn1   += tmp;
    tvar1  += tmp*cf1;
    tmp     = cf2*crv;
    tmn2   += tmp;
    tvar2  += tmp*cf2;
  }
  cf2     = RFlux2[n-1];
  crv     = rowptr[n-1];
  tmp     = 0.5*cf1*crv;
  tmn1   += tmp;
  tvar1  += tmp*cf1;
  tmp     = 0.5*cf2*crv;
  tmn2   += tmp;
  tvar2  += tmp*cf2;
   
  //Now apply the i = n-1 multiplier (0.5)
  mn1 += 0.5*tmn1;
  mn2 += 0.5*tmn2;
  var1 += 0.5*tvar1;
  var2 += 0.5*tvar2;
 
  //Grid size factor in integral
  double fluxfac = dflux1*dflux2;
  mn1   *= fluxfac;
  mn2   *= fluxfac;
  var1  *= fluxfac;
  var2  *= fluxfac;
  
}
  

/*!
  Gets ready for P(D) computation by preparing R
 
  \param[in] n Size of transform (square)
  \param[in] sigma1 Maximum allowed sigma in dimension 1
  \param[in] sigma2 Maximum allwed sigma in dimension 2
  \param[in] maxflux1 Maximum flux generated for R, dimension 1
  \param[in] maxflux2 Maximum flux generator for R, dimension 2
  \param[in] model    number counts model to use for fill.  Params must be set
  \param[in] bm       Beam
  \param[in] setEdge  Use integral of mean values at the edges

  Note that n is the transform size; the output array will generally
  be smaller because of padding.  Furthermore, because of mean shifting,
  the maximum flux often won't quite match the target values.
 */
void PDFactoryDouble::initPD(unsigned int n, double sigma1,
			     double sigma2, double maxflux1, 
			     double maxflux2, 
			     const numberCountsDouble& model,
			     const doublebeam& bm,
			     bool setEdge) {

  if (n == 0)
    throw affineExcept("PDFactoryDouble","initPD",
		       "Invalid (non-positive) n",1);  

  if (sigma1 < 0.0)
    throw affineExcept("PDFactoryDouble","initPD",
		       "Invalid (negative) sigma1",2);
  if (sigma2 < 0.0)
    throw affineExcept("PDFactoryDouble","initPD",
		       "Invalid (negative) sigma2",4);
  if (maxflux1 <= 0.0)
    throw affineExcept("PDFactoryDouble","initPD",
		       "Invalid (non-positive) maxflux1",8);
  if (maxflux2 <= 0.0)
    throw affineExcept("PDFactoryDouble","initPD",
		       "Invalid (non-positive) maxflux2",16);

  //Make sure we have enough room
  bool did_resize = resize(n);

  //Initial guess at dflux, will be changed later
  double inm1 = 1.0 / static_cast<double>(n-1);
  dflux1 = maxflux1 * inm1;
  dflux2 = maxflux2 * inm1;

  //Now fill R.  Some setup is required first
  if (!rvars_allocated) allocateRvars();
  //Fill in flux values for main and edge pieces
  //Main bit
  for (unsigned int i = 0; i < n; ++i)
    RFlux1[i] = static_cast<double>(i)*dflux1;
  for (unsigned int i = 0; i < n; ++i)
    RFlux2[i] = static_cast<double>(i)*dflux2;

  //Now fill in R.  The edges require special care.  This
  // first call will fill nonsense values into the lower edges, which
  // we will later overwrite
  double *rptr;  

#ifdef TIMING
  std::clock_t starttime = std::clock();
#endif
  model.getR(n,RFlux1,n,RFlux2,bm,rvals,numberCountsDouble::BEAMALL);
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
    double inedgem1 = 1.0/static_cast<double>(nedge-1);
    double dinterpfluxedge1, dinterpfluxedge2;
    double iRxnorm = 0.0, iRynorm = 0.0, iR00norm = 0.0;
    
    if (setEdge) {
      if (!edgevars_allocated) allocateEdgevars();
      if (use_edge_log_x) {
	dinterpfluxedge1 = -log(PDFactoryDouble::lowEdgeRMult)*inedgem1;
	for (unsigned int i = 0; i < nedge; ++i)
	  REdgeFlux1[i] = minedge1*exp(static_cast<double>(i)*dinterpfluxedge1);
      } else {
	dinterpfluxedge1 = (dflux1-minedge1)*inedgem1;
	for (unsigned int i = 0; i < nedge; ++i)
	  REdgeFlux1[i] = minedge1 + static_cast<double>(i)*dinterpfluxedge1;
      }
      if (use_edge_log_y) {
	dinterpfluxedge2 = -log(PDFactoryDouble::lowEdgeRMult)*inedgem1;
	for (unsigned int i = 0; i < nedge; ++i)
	  REdgeFlux2[i] = minedge2*exp(static_cast<double>(i)*dinterpfluxedge2);
      } else {
	dinterpfluxedge2 = (dflux2-minedge2)*inedgem1;
	for (unsigned int i = 0; i < nedge; ++i)
	  REdgeFlux2[i] = minedge2 + static_cast<double>(i)*dinterpfluxedge2;
      }
      iRxnorm  = dinterpfluxedge1/(dflux1-minedge1);
      iRynorm  = dinterpfluxedge2/(dflux2-minedge2);
      iR00norm = dinterpfluxedge1*dinterpfluxedge2/
	( (dflux1-minedge1) * (dflux2-minedge2) );
    }

    //First, do r[0,0]
    double scriptr;
#ifdef TIMING
    starttime = std::clock();
#endif
    model.getR(nedge,REdgeFlux1,nedge,REdgeFlux2,bm,REdgeWork,
	       numberCountsDouble::BEAMALL);
#ifdef TIMING
    RTime += std::clock() - starttime;
#endif

    //Do y integral first, store in REdgeWork[0,*]
    if (use_edge_log_y) {
      for (unsigned int i = 0; i < nedge; ++i) {
	rptr = REdgeWork + i*nedge; //row pointer
	scriptr = 0.5*(REdgeFlux2[0]*rptr[0] + 
		       REdgeFlux2[nedge-1]*rptr[nedge-1]);
	for (unsigned int j = 1; j < nedge-1; ++j)
	  scriptr += REdgeFlux2[j]*rptr[j];
	REdgeWork[i] = scriptr;
      }
    } else {
      for (unsigned int i = 0; i < nedge; ++i) {
	rptr = REdgeWork + i*nedge; 
	scriptr = 0.5*(rptr[0] + rptr[nedge-1]);
	for (unsigned int j = 1; j < nedge-1; ++j)
	  scriptr += rptr[j];
	REdgeWork[i] = scriptr;
      }
    }
    //Now X integral, put in integral step size and area
    // of bin, store in R[0,0]
    if (use_edge_log_x) {
      scriptr = 0.5*(REdgeFlux1[0]*REdgeWork[0]+
		     REdgeFlux1[nedge-1]*REdgeWork[nedge-1]);
      for (unsigned int i = 1; i < nedge-1; ++i)
	scriptr += REdgeFlux1[i]*REdgeWork[i];
      rvals[0] = scriptr*iR00norm;
    } else {
      scriptr = 0.5*(REdgeWork[0]+REdgeWork[nedge-1]);
      for (unsigned int i = 1; i < nedge-1; ++i)
	scriptr += REdgeWork[i];
      rvals[0] = scriptr*iR00norm;
    }
    
    //Now do Rx = R[0,y], integral along x
    double fixed_value;
    for (unsigned int j = 1; j < n; ++j) {
      fixed_value = RFlux2[j];
      //REdgeWork is more than big enough
#ifdef TIMING
      starttime = std::clock();
#endif
      model.getR(nedge,REdgeFlux1,1,&fixed_value,bm,REdgeWork,
		 numberCountsDouble::BEAMALL);
#ifdef TIMING
      RTime += std::clock() - starttime;
#endif
      if (use_edge_log_x) {
	scriptr = 0.5*(REdgeFlux1[0]*REdgeWork[0]+
		       REdgeFlux1[nedge-1]*REdgeWork[nedge-1]);
	for (unsigned int i = 1; i < nedge-1; ++i)
	  scriptr += REdgeFlux1[i]*REdgeWork[i];
      } else {
	scriptr = 0.5*(REdgeWork[0]+REdgeWork[nedge-1]);
	for (unsigned int i = 1; i < nedge-1; ++i)
	  scriptr += REdgeWork[i];
      }
      rvals[j] = scriptr*iRxnorm;
    }
    
    //And Ry = R[x,0]
    for (unsigned int i = 1; i < n; ++i) {
      fixed_value = RFlux1[i];
      //REdgeWork is more than big enough
#ifdef TIMING
      starttime = std::clock();
#endif
      model.getR(1,&fixed_value,nedge,REdgeFlux2,bm,REdgeWork,
		 numberCountsDouble::BEAMALL);
#ifdef TIMING
      RTime += std::clock() - starttime;
#endif
      if (use_edge_log_y) {
	scriptr = 0.5*(REdgeFlux2[0]*REdgeWork[0]+
		       REdgeFlux2[nedge-1]*REdgeWork[nedge-1]);
	for (unsigned int j = 1; j < nedge-1; ++j)
	  scriptr += REdgeFlux2[j]*REdgeWork[j];
      } else {
	scriptr = 0.5*(REdgeWork[0]+REdgeWork[nedge-1]);
	for (unsigned int j = 1; j < nedge-1; ++j)
	  scriptr += REdgeWork[j];
      }
      rvals[i*n] = scriptr*iRynorm;
    }
  } else {
    //Just set edges to zero
    for (unsigned int i = 0; i < n; ++i)
      rvals[i] = 0.0;
    for (unsigned int i = 1; i < n; ++i)
      rvals[i*n] = 0.0;
  }

  //Whew!  R is now filled.  Now estimate the mean and standard
  // deviation of the resulting P(D)
#ifdef TIMING
  starttime = std::clock();
#endif
  getRStats(n,mn1,mn2,var_noi1,var_noi2);
#ifdef TIMING
  RStatsTime += std::clock() - starttime;
#endif
  //Add the instrument noise into the variances
  sg1 = sqrt( var_noi1 + sigma1*sigma1 );
  sg2 = sqrt( var_noi2 + sigma2*sigma2 );


  //Decide if we will shift and pad, and if so by how much
  //Only do shifts if the noise is larger than one actual step size
  // Otherwise we can't represent it well.
  bool dopad1 = (sigma1 > dflux1);
  bool dopad2 = (sigma2 > dflux2);
  doshift1 = ( dopad1 && ( mn1 < pofd_mcmc::n_sigma_shift2d*sg1) );
  doshift2 = ( dopad2 && ( mn2 < pofd_mcmc::n_sigma_shift2d*sg2 ) );
  if ( doshift1 ) shift1 = pofd_mcmc::n_sigma_shift2d*sg1 - mn1; else
    shift1=0.0;
  if ( doshift2 ) shift2 = pofd_mcmc::n_sigma_shift2d*sg2 - mn2; else
    shift2=0.0;

  if (verbose) {
    std::cout << " Initial mean estimate band1: " << mn1 << " band2: "
	      << mn2 << std::endl;
    std::cout << " Initial stdev estimate band1: " << sg1 << " band2: "
	      << sg2 << std::endl;
    std::cout << " Additional shift applied band1: " << shift1 << " band2: "
	      << shift2 << std::endl;
  }

  //Make sure that maxflux is large enough that we don't get
  // bad aliasing wrap from the top around into the lower P(D) values.
  if (maxflux1 <= pofd_mcmc::n_sigma_pad*sg1)
    throw affineExcept("PDFactoryDouble","initPD",
		       "Top wrap problem, dimension 1",128);
  if (maxflux2 <= pofd_mcmc::n_sigma_pad*sg2)
    throw affineExcept("PDFactoryDouble","initPD",
		       "Top wrap problem, dimension 2",
		       256);

  //The other side of the equation is that we want to zero-pad the
  // top, and later discard that stuff.  The idea is as follows:
  // the 'target mean' of the calculated P(D) will lie at mn+shift
  // in each dimension.  We assume that anything within n_sigma_pad2d*sg
  // is 'contaminated'.  That means that, if n_sigma_pad2d*sg >
  // mn+shift, there will be some wrapping around the bottom of the P(D)
  // to contaminate the top by an amount n_sigma_pad2d*sg - (mn+shift).
  // We therefore zero pad and discard anything above
  // maxflux - (n_sigma_pad2d*sg - (mn+shift))
  if (dopad1) {
    double contam = pofd_mcmc::n_sigma_pad2d*sg1 - (mn1+shift1);
    if (contam <= 0) maxidx1 = n; else {
      double topflux = maxflux1 - contam;
      if (topflux < 0)
	throw affineExcept("PDFactoryDouble","initPD",
			 "Padding problem, dimension 1",512);
      maxidx1 = static_cast< unsigned int>( topflux/dflux1 );
      if (maxidx1 > n)
	throw affineExcept("PDFactoryDouble","initPD",
			 "Padding problem, dimension 1",
			 1024);
      //Now pad!
      for (unsigned int i = maxidx1; i < n; ++i) {
	rptr = rvals + n*i;
	for (unsigned int j = 0; j < n; ++j)
	  rptr[j] = 0.0;
      }
    }
  } else maxidx1 = n;

  if (dopad2) {
    double contam = pofd_mcmc::n_sigma_pad2d*sg2 - (mn2+shift2);
      if (contam <= 0) maxidx2 = n; else {
      double topflux = maxflux2 - contam;
      if (topflux < 0)
	throw affineExcept("PDFactoryDouble","initPD",
			 "Padding problem, dimension 2",512);
      maxidx2 = static_cast< unsigned int>( topflux/dflux2 );
      if (maxidx2 > n)
	throw affineExcept("PDFactoryDouble","initPD",
			 "Padding problem, dimension 2",1024);
      for (unsigned int i = 0; i < maxidx1; ++i) {
	rptr = rvals + n*i;
	for (unsigned int j = maxidx2; j < n; ++j)
	  rptr[j] = 0.0;
      }
      for (unsigned int i = maxidx2; i < n; ++i) {
	rptr = rvals + n*i;
	for (unsigned int j = 0; j < n; ++j)
	  rptr[j] = 0.0;
      }
    }
  } else maxidx2 = n;

  //Multiply r by dflux factor to represent the actual
  // number of sources in each bin.  Note we do this
  // after calling getRStats, since the thing we have
  // is no longer strictly R after this.
  double dfactor = dflux1*dflux2;
  for (unsigned int i = 0; i < maxidx1; ++i) {
    rptr = rvals + n*i;
    for (unsigned int j = 0; j < maxidx2; ++j)
      rptr[j] *= dfactor;
  }
  
  //Make the plans, or keep the old ones if possible
  //Note that the forward transform dumps into out_part,
  // but the backwards one comes from pval.  The idea
  // is that out_part holds the working bit.  These are
  // convolved together into pval.  This is inefficient
  // if there is only one sign present, but the special
  // case doesn't seem worth the effort
  //If we resized, we must make the new plans because the
  // addresses changed
  //We will have to use the advanced interfact to point
  // specifically at the rvals subindex we are using on the
  // forward plan, but the backwards plan is fine
  int intn = static_cast<int>(n);
  if ( did_resize || (lastfftlen != n) || (plan == NULL) ) {
    if (plan != NULL) fftw_destroy_plan(plan);
    plan = fftw_plan_dft_r2c_2d(intn, intn, rvals, rtrans,
				fftw_plan_style);
  }
  if ( did_resize || (lastfftlen != n) || (plan_inv == NULL) ) {
    if (plan_inv != NULL) fftw_destroy_plan(plan_inv);
    plan_inv = fftw_plan_dft_c2r_2d(intn, intn, pval, pofd,
				    fftw_plan_style);
  }
  if (plan == NULL) {
    std::stringstream str;
    str << "Plan creation failed for forward transform of size: " << 
      n << std::endl;
    throw affineExcept("PDFactoryDouble","initPD",str.str(),32);
  }
  if (plan_inv == NULL) {
    std::stringstream str;
    str << "Plan creation failed for inverse transform of size: " << 
      n << std::endl;
    throw affineExcept("PDFactoryDouble","initPD",str.str(),64);
  }

  //Compute forward transform of this r value, store in rtrans
#ifdef TIMING
  starttime = std::clock();
#endif
  fftw_execute_dft_r2c(plan,rvals,rtrans);
#ifdef TIMING
  fftTime += std::clock() - starttime;
#endif
  
  lastfftlen = n;
  max_sigma1 = sigma1;
  max_sigma2 = sigma2;
  initialized = true;
}

/*!
  Calculates P(D) for a dataset using direct fills
 
  \param[in] sigma1 Sigma in dimension 1
  \param[in] sigma2 Sigma in dimension 2
  \param[out] pd Holds P(D1,D2) on output
  \param[in] setLog If true, pd is log(P(D1,D2)) on output; convenient
              for likelihood evaluation.
  \param[in] edgeFix  Apply a fix to the lower edges to minimize wrapping
                       effects using a Gaussian to each row/col

  Note that n is the transform size; the output array will generally
  be smaller because of padding.  Furthermore, because of mean shifting,
  the maximum flux often won't quite match the target values.

  You must call initPD first, or bad things will probably happen.
*/
void PDFactoryDouble::getPD( double sigma1, double sigma2,  
			     PDDouble& pd, bool setLog, bool edgeFix) {

  // The basic idea is to compute the P(D) from the previously filled
  // R values, adding in noise and all that fun stuff, filling pd
  // for output

  if (! initialized )
    throw affineExcept("PDFactoryDouble","getPD",
		     "Must call initPD first",1);
  if (sigma1 > max_sigma1) {
    std::stringstream errstr("");
    errstr << "Sigma 1 value " << sigma1 
	   << " larger than maximum prepared value " << max_sigma1
	   << std::endl;
    errstr << "initPD should have been called with at least " << sigma1;
    throw affineExcept("PDFactoryDouble","getPD",errstr.str(),2);
  }
  if (sigma2 > max_sigma2) {
    std::stringstream errstr("");
    errstr << "Sigma 1 value " << sigma2 
	   << " larger than maximum prepared value " << max_sigma2
	   << std::endl;
    errstr << "initPD should have been called with at least " << sigma2;
    throw affineExcept("PDFactoryDouble","getPD",errstr.str(),2);
  }

  //Output array from 2D FFT is n * (n/2+1)
  unsigned int n = lastfftlen;
  unsigned int ncplx = n/2 + 1;
      
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

  if (doshift1 && doshift2) {
    //Both shifted -- this should be the most common case,
    // and corresponds to noise in both bands
    double sigfac1 = 0.5*sigma1*sigma1;
    double sigfac2 = 0.5*sigma2*sigma2;
    
    double sigprod1, meanprod1, didx2, w1, w2;
    
    //First, Pos freq
    for (unsigned int idx1 = 0; idx1 < ncplx; ++idx1) {
      r_input_rowptr = rtrans + idx1*ncplx; //Input
      row_current_out = pval + idx1*ncplx; //Output
      w1 = iflux1 * static_cast<double>(idx1);
      meanprod1 = shift1*w1;
      sigprod1  = sigfac1*w1*w1;
      for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
	didx2 = static_cast<double>(idx2);
	w2 = iflux2 * didx2;
	rval = r_input_rowptr[idx2][0] - r0 - sigprod1 - sigfac2*w2*w2;
	ival = r_input_rowptr[idx2][1] - meanprod1 - shift2*w2;
	expfac = exp( rval );
	row_current_out[idx2][0] = expfac*cos(ival);
	row_current_out[idx2][1] = expfac*sin(ival);
      }
    }
    //Now, Neg freq
    for (unsigned int idx1 = ncplx; idx1 < n; ++idx1) {
      r_input_rowptr = rtrans + idx1*ncplx; //Input
      row_current_out = pval + idx1*ncplx; //Output
      w1 = - iflux1 * static_cast<double>(n - idx1);
      meanprod1 = shift1*w1;
      sigprod1  = sigfac1*w1*w1;
      for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
	didx2 = static_cast<double>(idx2);
	w2 = iflux2 * didx2;
	rval = r_input_rowptr[idx2][0] - r0 - sigprod1 - sigfac2*w2*w2;
	ival = r_input_rowptr[idx2][1] - meanprod1 - shift2*w2;
	expfac = exp( rval );
	row_current_out[idx2][0] = expfac*cos(ival);
	row_current_out[idx2][1] = expfac*sin(ival);
      }
    }
    //p(0,0) is special
    pval[0][0] = 1.0;
    pval[0][1] = 0.0;
  } else if (doshift1) {
    //Only shift in band 1
    double sigfac1 = 0.5*sigma1*sigma1;
    double sigprod1, meanprod1, w1;
    for (unsigned int idx1 = 0; idx1 < ncplx; ++idx1) {
      r_input_rowptr = rtrans + idx1*ncplx; //Input
      row_current_out = pval + idx1*ncplx; //Output
      w1 = iflux1 * static_cast<double>(idx1);
      meanprod1 = shift1*w1;
      sigprod1  = sigfac1*w1*w1;
      for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
	rval = r_input_rowptr[idx2][0] - r0 - sigprod1;
	ival = r_input_rowptr[idx2][1] - meanprod1;
	expfac = exp( rval );
	row_current_out[idx2][0] = expfac*cos(ival);
	row_current_out[idx2][1] = expfac*sin(ival);
      }
    }
    for (unsigned int idx1 = ncplx; idx1 < n; ++idx1) {
      r_input_rowptr = rtrans + idx1*ncplx; //Input
      row_current_out = pval + idx1*ncplx; //Output
      w1 = - iflux1 * static_cast<double>(n - idx1);
      meanprod1 = shift1*w1;
      sigprod1  = sigfac1*w1*w1;
      for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
	rval = r_input_rowptr[idx2][0] - r0 - sigprod1;
	ival = r_input_rowptr[idx2][1] - meanprod1;
	expfac = exp( rval );
	row_current_out[idx2][0] = expfac*cos(ival);
	row_current_out[idx2][1] = expfac*sin(ival);
      }
    }
    pval[0][0] = 1.0;
    pval[0][1] = 0.0;
  } else if (doshift2) {
    //And only shift band 2
    //Only flux 2 shifted, has sigma
    //Can do all of band 1 in one loop, since we ignore
    // freq1
    double sigfac2 = 0.5*sigma2*sigma2;
    double didx2, w2;
    for (unsigned int idx1 = 0; idx1 < n; ++idx1) {
      r_input_rowptr = rtrans + idx1*ncplx; //Input
      row_current_out = pval + idx1*ncplx; //Output
      for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
	didx2 = static_cast<double>(idx2);
	w2 = iflux2 * didx2;
	rval = r_input_rowptr[idx2][0] - r0 - sigfac2*w2*w2;
	ival = r_input_rowptr[idx2][1] - shift2*w2;
	expfac = exp( rval );
	row_current_out[idx2][0] = expfac*cos(ival);
	row_current_out[idx2][1] = expfac*sin(ival);
      }
    }
    pval[0][0] = 1.0;
    pval[0][1] = 0.0;
  } else {
    //No shifts, sigmas
    for (unsigned int idx1 = 0; idx1 < n; ++idx1) {
      r_input_rowptr = rtrans + idx1*ncplx;
      row_current_out = pval + idx1*ncplx;
      for (unsigned int idx2 = 0; idx2 < ncplx; ++idx2) {
	rval = r_input_rowptr[idx2][0] - r0;
	ival = r_input_rowptr[idx2][1];
	expfac = exp(rval);
	row_current_out[idx2][0] = expfac*cos(ival);
	row_current_out[idx2][1] = expfac*sin(ival);
      }
    }
    pval[0][0] = 1.0;
    pval[0][1] = 0.0;
  }
#ifdef TIMING
  p0Time += std::clock() - starttime;
#endif

  
  //Now transform back
  if (verbose) std::cout << " Reverse transform" << std::endl;
#ifdef TIMING
  starttime = std::clock();
#endif
  fftw_execute(plan_inv); //overwrites pofd with reverse transform of pval
#ifdef TIMING
  fftTime += std::clock() - starttime;
#endif

  //Enforce positivity
#ifdef TIMING
  starttime = std::clock();
#endif
  for (unsigned int idx = 0; idx < n*n; ++idx)
    if (pofd[idx] < 0) pofd[idx] = 0.0;
#ifdef TIMING
  posTime += std::clock() - starttime;
#endif

  //Copy into output variable
  pd.resize(maxidx1,maxidx2);
  double *pdptr, *pofdptr;
#ifdef TIMING
  starttime = std::clock();
#endif
  for (unsigned int i = 0; i < maxidx1; ++i) {
    pdptr = pd.pd_ + maxidx2*i; //Row pointer for output variable
    pofdptr = pofd + n*i; //Row pointer for variable we are copying
    for (unsigned int j = 0; j < maxidx2; ++j)
      pdptr[j] = pofdptr[j];
  }
  pd.logflat = false;
  pd.minflux1 = 0.0; pd.dflux1 = dflux1;
  pd.minflux2 = 0.0; pd.dflux2 = dflux2;
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

  //Edge fix
#ifdef TIMING
  starttime = std::clock();
#endif
  if (edgeFix) pd.edgeFix();
#ifdef TIMING
  edgeTime += std::clock() - starttime;
#endif
  

  //Now mean subtract each axis
#ifdef TIMING
  starttime = std::clock();
#endif
  double tmn1, tmn2; //True means
  pd.getMeans(tmn1,tmn2,false);
  if ( std::isinf(tmn1) || std::isnan(tmn1) || std::isinf(tmn2) ||
       std::isnan(tmn2) ) {
    std::stringstream str;
    str << "Un-shift amounts not finite band1: " << tmn1 << " "
	<< tmn2 << std::endl;
    str << "At length: " << n << " with noise: " << sigma1 << " "
	<< sigma2;
    throw affineExcept("PDFactoryDouble","getPD",str.str(),8);
  }
  if (verbose) {
    std::cout << " Expected mean band1: " << std::fixed 
	      << std::setprecision(6) << shift1+mn1 << " band2: "
	      << shift2 + mn2 << std::endl;
    std::cout << " Realized mean band1: " << std::fixed 
	      << std::setprecision(6) << tmn1 << " band2: "
	      << tmn2 << std::endl;
  }
  pd.minflux1 = -tmn1;
  pd.minflux2 = -tmn2;
#ifdef TIMING
  meanTime += std::clock() - starttime;
#endif

  //Turn PD to log for more efficient log computation of likelihood
#ifdef TIMING
  starttime = std::clock();
#endif
  if (setLog) pd.applyLog(false);
#ifdef TIMING
  logTime += std::clock() - starttime;
#endif
}
 
void PDFactoryDouble::SendSelf(MPI::Comm& comm, int dest) const {
  comm.Send(&fftw_plan_style,1,MPI::UNSIGNED,dest,pofd_mcmc::PDFDSENDPLANSTYLE);
  comm.Send(&has_wisdom,1,MPI::UNSIGNED,dest,pofd_mcmc::PDFDHASWISDOM);
  if (has_wisdom) {
    //Send wisdom file name
    unsigned int nstr = wisdom_file.size()+1;
    char *cstr = new char[nstr];
    std::strncpy( cstr, wisdom_file.c_str(), nstr );
    comm.Send(&nstr,1,MPI::UNSIGNED,dest,pofd_mcmc::PDFDWISLEN);
    comm.Send(cstr, nstr, MPI::CHAR, dest,pofd_mcmc::PDFDWISNAME);
    delete[] cstr;
  }
  comm.Send(&verbose,1,MPI::BOOL,dest,pofd_mcmc::PDFDVERBOSE);
  comm.Send(&nedge,1,MPI::UNSIGNED,dest,pofd_mcmc::PDFDNEDGE);
}

//Note this doesn't copy over interal variables
void PDFactoryDouble::RecieveCopy(MPI::Comm& comm, int src) {
  comm.Recv(&fftw_plan_style,1,MPI::UNSIGNED,src,pofd_mcmc::PDFDSENDPLANSTYLE);
  comm.Recv(&has_wisdom,1,MPI::UNSIGNED,src,pofd_mcmc::PDFDHASWISDOM);
  if (has_wisdom) {
    //Recieve wisdom file name
    unsigned int nstr;
    comm.Recv(&nstr,1,MPI::UNSIGNED,src,pofd_mcmc::PDFDWISLEN);
    char *cstr = new char[nstr];
    comm.Recv(cstr, nstr, MPI::CHAR, src,pofd_mcmc::PDFDWISNAME);    
    wisdom_file = std::string(cstr);
    delete[] cstr;
    addWisdom(wisdom_file);
  }
  comm.Recv(&verbose,1,MPI::BOOL,src,pofd_mcmc::PDFDVERBOSE);
  unsigned int newnedge;
  comm.Recv(&newnedge,1,MPI::UNSIGNED,src,pofd_mcmc::PDFDNEDGE);
  if (newnedge != nedge) {
    freeEdgevars();
    nedge   = newnedge;
  }
  initialized = false;
}
