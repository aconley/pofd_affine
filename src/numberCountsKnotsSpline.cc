#include<cmath>
#include<sstream>

#include<numberCountsKnotsSpline.h>
#include<global_settings.h>
#include<affineExcept.h>

//Function to pass to GSL integrator
static double evalPowfNKnotsSpline(double,void*); //!< Evaluates f^pow dN/dS


numberCountsKnotsSpline::numberCountsKnotsSpline() : numberCountsKnots() {
  logknots = NULL;
  acc = gsl_interp_accel_alloc();
  splinelog = NULL;
  gsl_work = gsl_integration_workspace_alloc(1000);
  varr = new void*[5];
}

numberCountsKnotsSpline::numberCountsKnotsSpline(unsigned int NKNOTS) :
  numberCountsKnots(NKNOTS) {
  if (NKNOTS > 0) logknots = new double[NKNOTS]; else logknots = NULL;
  acc = gsl_interp_accel_alloc();
  if (NKNOTS > 0)
    splinelog = gsl_spline_alloc( gsl_interp_cspline,
				  static_cast<size_t>(NKNOTS));
  else {
    splinelog=NULL;
  }
  gsl_work = gsl_integration_workspace_alloc(1000);
  varr = new void*[5];
}

numberCountsKnotsSpline::numberCountsKnotsSpline(const std::vector<double>& v):
  numberCountsKnots(v) {
  if (nknots > 0) {
    logknots = new double[nknots]; 
    for (unsigned int i = 0; i < nknots; ++i)
      logknots[i] = log2( knots[i] );
  } else logknots = NULL;
  acc = gsl_interp_accel_alloc();
  if (v.size() > 0)
    splinelog = gsl_spline_alloc( gsl_interp_cspline,
				  static_cast<size_t>(v.size()));
  else
    splinelog = NULL;
  gsl_work = gsl_integration_workspace_alloc(1000);
  varr = new void*[5];
}
numberCountsKnotsSpline::numberCountsKnotsSpline( unsigned int n, 
						  const double* const S) :
  numberCountsKnots(n,S) {
  if (nknots > 0) {
    logknots = new double[nknots]; 
    for (unsigned int i = 0; i < nknots; ++i)
      logknots[i] = log2( knots[i] );
  } else logknots = NULL;
  acc = gsl_interp_accel_alloc();
  splinelog = gsl_spline_alloc( gsl_interp_cspline,
				static_cast<size_t>(n));
  gsl_work = gsl_integration_workspace_alloc(1000);
  varr = new void*[5];
}

numberCountsKnotsSpline::numberCountsKnotsSpline( const numberCountsKnotsSpline& 
						  other ){
  if ( this == &other ) return; //Self-copy
  nknots = 0;
  knots = logknots = logknotvals = NULL;
  acc = NULL;
  splinelog=NULL;

  setNKnots(other.nknots);
  for (unsigned int i = 0; i < nknots; ++i)
    knots[i] = other.knots[i];
  for (unsigned int i = 0; i < nknots; ++i)
    logknots[i] = other.logknots[i];
  if (other.knotvals_loaded) {
    for (unsigned int i = 0; i < nknots; ++i)
      logknotvals[i] = other.logknotvals[i];
    gsl_spline_init( splinelog, logknots, logknotvals,
		     static_cast<size_t>(other.nknots) );
  }
  knotvals_loaded = other.knotvals_loaded;
  gsl_work = gsl_integration_workspace_alloc(1000);
  varr = new void*[5];
}

numberCountsKnotsSpline::~numberCountsKnotsSpline() {
  gsl_interp_accel_free(acc);
  if (splinelog != NULL) gsl_spline_free(splinelog);
  gsl_integration_workspace_free(gsl_work);
  delete[] varr;
}

void numberCountsKnotsSpline::setNKnots(unsigned int n) {
  if ( nknots == n ) return;
  if ( knots != NULL ) delete[] knots;
  if ( logknotvals != NULL ) delete[] logknotvals;
  if ( logknots != NULL) delete[] logknots;
  if (splinelog != NULL) gsl_spline_free(splinelog);

  if ( n > 0 ) {
    knots = new double[n];
    logknots = new double[n];
    logknotvals = new double[n];
  } else {
    knots = logknots = logknotvals = NULL;
  }
  if ( n > 1 )
    splinelog=gsl_spline_alloc( gsl_interp_cspline,
				static_cast<size_t>(n));
  else 
    splinelog=NULL;
  nknots = n;
  knotvals_loaded = false;
}

numberCountsKnotsSpline& 
numberCountsKnotsSpline::operator=(const numberCountsKnotsSpline& other) {
  if ( this == &other ) return *this; //Self-copy
  if (other.nknots == 0) {
    if (nknots != 0) {
      if (knots != NULL) delete[] knots;
      if (logknots != NULL) delete[] logknots;
      if (logknotvals != NULL) delete[] logknotvals;
      if (splinelog != NULL) gsl_spline_free(splinelog);
      knots = logknots = logknotvals = NULL;
      splinelog = NULL;
    }
    knotvals_loaded = false;
  } else {
    setNKnots(other.nknots);
    for (unsigned int i = 0; i < nknots; ++i)
      knots[i] = other.knots[i];
    for (unsigned int i = 0; i < nknots; ++i)
      logknots[i] = other.logknots[i];
    if (other.knotvals_loaded) {
      for (unsigned int i = 0; i < nknots; ++i)
	logknotvals[i] = other.logknotvals[i];
      gsl_spline_init( splinelog, logknots, logknotvals,
		       static_cast<size_t>(other.nknots) );
    }
    knotvals_loaded = other.knotvals_loaded;
  }
  return *this;
}
 
void numberCountsKnotsSpline::setKnotPositions(const std::vector<double>& vec) {
  numberCountsKnots::setKnotPositions(vec);
  for (unsigned int i = 0; i < nknots; ++i)
    logknots[i] = log2(knots[i]);
}

void numberCountsKnotsSpline::setKnotPositions(unsigned int n, 
					       const double* const vals) {
  numberCountsKnots::setKnotPositions(n,vals);
  for (unsigned int i = 0; i < nknots; ++i)
    logknots[i] = log2(knots[i]);
}


/*!
  \param[in] F Parameters.  Will ignore any after
                    the first nknots
 */
void numberCountsKnotsSpline::setParams(const paramSet& F) {
  if (nknots > F.getNParams())
    throw affineExcept("numberCountsKnotsSpline","setKnots",
		       "Number of knot values different than expected",2);
  for (unsigned int i = 0; i < nknots; ++i)
    logknotvals[i] = pofd_mcmc::logfac*F[i];
  gsl_spline_init( splinelog, logknots, logknotvals,
		   static_cast<size_t>(nknots) );
  knotvals_loaded = true;
}


/*!
  Note that this has some overhead, so within inner computation loops
  shouldn't be called directly.
 */
double numberCountsKnotsSpline::getNumberCounts(double value) const {
  if (nknots < 2) return std::numeric_limits<double>::quiet_NaN();
  if (! isValid() ) return std::numeric_limits<double>::quiet_NaN();
  if (value < knots[0] || value >= knots[nknots-1]) return 0.0; //Out of range
  return exp2( gsl_spline_eval( splinelog, log2(value), acc ) ); 
}

/*!
  Computes
  \f[
   \int dS\, S^{\alpha} \frac{dN}{dS}
  \f]

  \param[in] alpha  Power of flux
  \returns Integral
 */
double numberCountsKnotsSpline::splineInt(double alpha) const {
  if (nknots < 2) return std::numeric_limits<double>::quiet_NaN();
  if (! isValid() ) return std::numeric_limits<double>::quiet_NaN();
  double result, error;
  void *params;
  gsl_function F;
  double minknot = knots[0];
  double maxknot = knots[nknots-1];

  varr[0] = static_cast<void*>(&alpha);
  varr[1] = static_cast<void*>(splinelog);
  varr[2] = static_cast<void*>(acc);
  varr[3] = static_cast<void*>(&minknot);
  varr[4] = static_cast<void*>(&maxknot);
  params = static_cast<void*>(varr);

  F.function = &evalPowfNKnotsSpline;
  F.params = params;

  gsl_integration_qag(&F, minknot, maxknot, 0.0, 1e-5, 1000,
  		      GSL_INTEG_GAUSS41, gsl_work, &result, &error); 
  return result;
}

double numberCountsKnotsSpline::getNS() const {
  return splineInt(0);
}

double numberCountsKnotsSpline::getMeanFluxPerArea() const {
  return splineInt(1);
}

double numberCountsKnotsSpline::getMeanFluxSqPerArea() const {
  return splineInt(2);
}
  
// Pure brute force
/*
  \params[in] fluxdensity  Position to get R at
  \params[in] params Parameters
  \params[in] beam Beam
  \returns R(fluxdensity) for the positive beam
 */
double numberCountsKnotsSpline::getRPos(double fluxdensity,
					const beam& bm) const {

  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  if (! bm.hasPos()) return std::numeric_limits<double>::quiet_NaN();

  if (fluxdensity <= 0.0) return 0.0;
  double minknot = knots[0];
  double maxknot = knots[nknots-1];
  if (fluxdensity > maxknot) return 0.0; //Since max(beam) = 1

  double prefac;
  prefac = bm.getPixSize()/3600.0;  //To sq deg
  prefac = prefac*prefac;

  unsigned int npsf = bm.getNPos();
  
  double currarg, ieta, retval;
  retval = 0.0;
  const double *ipixarr = bm.getPosInvPixArr();
  if (bm.hasPosWeights()) {
    const double* warr = bm.getPosWeights();
    for (unsigned int i = 0; i < npsf; ++i) {
      ieta = ipixarr[i];
      currarg = fluxdensity * ieta;
      if (currarg < minknot || currarg >= maxknot) continue;
      retval += warr[i]*
	exp2( gsl_spline_eval( splinelog, log2(currarg), acc ) )*ieta;
    } 
  } else {
    for (unsigned int i = 0; i < npsf; ++i) {
      ieta = ipixarr[i];
      currarg = fluxdensity * ieta;
      if (currarg < minknot || currarg >= maxknot) continue;
      retval += exp2( gsl_spline_eval( splinelog, log2(currarg), acc ) )*ieta;
    } 
  }

  return prefac*retval;
}

/*!
  \param[in] n Number of elements to get R for (length of flux)
  \param[in] flux Fluxes to get R for
  \param[in] bm   Beam
  \param[in,out] R Values of R from the positive beam are added onto this 
                 (len n)

  Note you need to zero your input array first if you just want R+
 */
void numberCountsKnotsSpline::getRPos(unsigned int n,const double* const flux,
				      const beam& bm,double* R) const {
  if (n == 0) return;
  if ( (!isValid()) || (!bm.hasPos()) ) {
    for (unsigned int i = 0; i < n; ++i)
      R[i] = std::numeric_limits<double>::quiet_NaN();
    return;
  }

  double minknot = knots[0];
  double maxknot = knots[nknots-1];

  double prefac;
  prefac = bm.getPixSize()/3600.0;  //To sq deg
  prefac = prefac*prefac;

  unsigned int npsf = bm.getNPos();
  
  double currarg, ieta, workval, fluxdensity;
  const double *ipixarr = bm.getPosInvPixArr();
  if (bm.hasPosWeights()) { 
    //Binned beam
    const double* warr = bm.getPosWeights();
    for (unsigned int j = 0; j < n; ++j) {
      fluxdensity = flux[j];
      workval = 0.0;
      for (unsigned int i = 0; i < npsf; ++i) {
	ieta = ipixarr[i];
	currarg = fluxdensity * ieta;
	if (currarg < minknot || currarg >= maxknot) continue;
	workval += warr[i]*
	  exp2( gsl_spline_eval( splinelog, log2(currarg), acc ) )*ieta;
      } 
      R[j] += prefac*workval;
    }
  } else {
    //Unbinned beam
    for (unsigned int j = 0; j < n; ++j) {
      fluxdensity = flux[j];
      workval = 0.0;
      for (unsigned int i = 0; i < npsf; ++i) {
	ieta = ipixarr[i];
	currarg = fluxdensity * ieta;
	if (currarg < minknot || currarg >= maxknot) continue;
	workval += exp2( gsl_spline_eval( splinelog, log2(currarg), acc ) ) * 
	  ieta;
      } 
      R[j] += prefac*workval;
    }
  }
}

/*
  \params[in] fluxdensity  Position to get R at
  \params[in] params Parameters
  \params[in] bm Beam
  \returns R(fluxdensity) for the negative beam
 */
double numberCountsKnotsSpline::getRNeg(double fluxdensity,
					const beam& bm) const {

  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  if (!bm.hasNeg()) return std::numeric_limits<double>::quiet_NaN();

  if (fluxdensity <= 0.0) return 0.0;
  double minknot = knots[0];
  double maxknot = knots[nknots-1];
  if (fluxdensity > maxknot) return 0.0; //Since max(beam) = 1

  double prefac;
  prefac = bm.getPixSize()/3600.0;  //To sq deg
  prefac = prefac*prefac;

  unsigned int npsf = bm.getNNeg();
  
  double currarg, ieta, retval;
  retval = 0.0;
  const double *ipixarr = bm.getNegInvPixArr();
  if (bm.hasNegWeights()) {
    const double* warr = bm.getPosWeights();
    for (unsigned int i = 0; i < npsf; ++i) {
      ieta = ipixarr[i];
      currarg = fluxdensity * ieta;
      if (currarg < minknot || currarg >= maxknot) continue;
      retval += warr[i]*
	exp2( gsl_spline_eval( splinelog, log2(currarg), acc ) )*ieta;
    } 
  } else {
    for (unsigned int i = 0; i < npsf; ++i) {
      ieta = ipixarr[i];
      currarg = fluxdensity * ieta;
      if (currarg < minknot || currarg >= maxknot) continue;
      retval += exp2( gsl_spline_eval( splinelog, log2(currarg), acc ) )*ieta;
    } 
  }
  return prefac*retval;
}

/*!
  \param[in] n Number of elements to get R for (length of flux)
  \param[in] flux Fluxes to get R for
  \param[in] bm   Beam
  \param[in,out] R Values of R from the negative beam are added onto this 
          (len n)
 */
void numberCountsKnotsSpline::getRNeg(unsigned int n,const double* const flux,
				      const beam& bm,double* R) const {
  if (n == 0) return;
  if ( (!isValid()) || (!bm.hasNeg()) ) {
    for (unsigned int i = 0; i < n; ++i)
      R[i] = std::numeric_limits<double>::quiet_NaN();
    return;
  }

  double minknot = knots[0];
  double maxknot = knots[nknots-1];

  double prefac;
  prefac = bm.getPixSize()/3600.0;  //To sq deg
  prefac = prefac*prefac;

  unsigned int npsf = bm.getNNeg();
  double currarg, ieta, workval, fluxdensity;
  const double *ipixarr = bm.getNegInvPixArr();
  if (bm.hasNegWeights()) {
    const double* warr = bm.getPosWeights();
    for (unsigned int j = 0; j < n; ++j) {
      fluxdensity = flux[j];
      workval = 0.0;
      for (unsigned int i = 0; i < npsf; ++i) {
	ieta = ipixarr[i];
	currarg = fluxdensity * ieta;
	if (currarg < minknot || currarg >= maxknot) continue;
	workval += warr[i]*
	  exp2( gsl_spline_eval( splinelog, log2(currarg), acc ) )*ieta;
      } 
      R[j] += prefac*workval;
    }
  } else {
    for (unsigned int j = 0; j < n; ++j) {
      fluxdensity = flux[j];
      workval = 0.0;
      for (unsigned int i = 0; i < npsf; ++i) {
	ieta = ipixarr[i];
	  currarg = fluxdensity * ieta;
	if (currarg < minknot || currarg >= maxknot) continue;
	workval += exp2( gsl_spline_eval( splinelog, log2(currarg), acc ) )
	  * ieta;
      } 
      R[j] += prefac*workval;
    }
  }
}

/*
  \param[in] fluxdensity  Position to get R at
  \param[in] params Parameters
  \param[in] bm Beam
  \param[in] rt Controls whether one gets the positive, negative, or sum
       of both
  \returns R(fluxdensity)
*/
double numberCountsKnotsSpline::getR(double fluxdensity,const beam& bm,
				     rtype rt) const {
  switch (rt) {
  case BEAMBOTH :
    if (bm.hasPos()) {
      if (bm.hasNeg()) return getRPos(fluxdensity,bm)+getRNeg(fluxdensity,bm);
      else return getRPos(fluxdensity,bm);
    }
    break;
  case BEAMPOS :
    return getRPos(fluxdensity,bm);
    break;
  case BEAMNEG :
    return getRNeg(fluxdensity,bm);
    break;
  default :
    std::stringstream errstr;
    errstr << "Unknown R beam specification: " << rt;
    throw affineExcept("numberCountsKnotsSpline","getR",
		       errstr.str(),1);
  }
  return std::numeric_limits<double>::quiet_NaN();
}

/*!
  \param[in] n Number of elements to get R for (length of flux)
  \param[in] flux Fluxes to get R for
  \param[in] bm Beam
  \param[in,out] R Values of R from the negative beam are added onto 
                 this (len n)
  \param[in] rt Type of R desired (pos, neg, sum of both)
 */
void numberCountsKnotsSpline::getR(unsigned int n,const double* const flux,
				   const beam& bm,double* R, rtype rt) const {
  if (n == 0) return;
  for (unsigned int i = 0; i < n; ++i) R[i] = 0.0;

  switch (rt) {
  case BEAMBOTH :
    if (bm.hasPos()) getRPos(n,flux,bm,R);
    if (bm.hasNeg()) getRNeg(n,flux,bm,R);
    break;
  case BEAMPOS :
    getRPos(n,flux,bm,R);
    break;
  case BEAMNEG :
    getRNeg(n,flux,bm,R);
    break;
  default :
    std::stringstream errstr;
    errstr << "Unknown R beam specification: " << rt;
    throw affineExcept("numberCountsKnotsSpline","getR",
		       errstr.str(),1);
  }
  return;
}

void numberCountsKnotsSpline::SendSelf(MPI::Comm& comm, int dest) const {
  numberCountsKnots::SendSelf(comm,dest);
  if (knotvals_loaded) {
    comm.Send(logknots,nknots,MPI::DOUBLE,dest,
	      pofd_mcmc::NCKSSENDLOGKNOTS);
  }
}

void numberCountsKnotsSpline::RecieveCopy(MPI::Comm& comm, int src) {
  unsigned int oldnknots = nknots;
  numberCountsKnots::RecieveCopy(comm,src);

  if (nknots != oldnknots) {
    if (splinelog != NULL) gsl_spline_free(splinelog);
    if (nknots > 0) splinelog = gsl_spline_alloc( gsl_interp_cspline,
						  static_cast<size_t>(nknots));
  }
  if (knotvals_loaded) {
    comm.Recv(logknots,nknots,MPI::DOUBLE,src,pofd_mcmc::NCKSSENDLOGKNOTS);
    gsl_spline_init( splinelog, logknots, logknotvals,
		     static_cast<size_t>(nknots) );
  }
}

static double evalPowfNKnotsSpline(double x, void* params) {
  //Params are:
  // parmas[0] power
  // params[1] spline (log2)
  // params[2] accelerator
  // params[3] minknot
  // params[4] maxknot
  //But this really has to be an array of pointers to void to work
  void** vptr = static_cast<void**>(params);

  double power = *static_cast<double*>(vptr[0]);
  gsl_spline* spl = static_cast<gsl_spline*>(vptr[1]);
  gsl_interp_accel* acc = static_cast<gsl_interp_accel*>(vptr[2]);
  double minknot = *static_cast<double*>(vptr[3]);
  double maxknot = *static_cast<double*>(vptr[4]);
  
  if (x < minknot || x >= maxknot) return 0.0;

  double splval = exp2( gsl_spline_eval(spl,log2(x),acc) );

  if ( fabs(power) < 1e-2 ) return splval;
  if ( fabs(power-1.0) < 1e-2) return x*splval;
  if ( fabs(power-2.0) < 1e-2) return x*x*splval;
  return pow(x,power)*splval;

}
