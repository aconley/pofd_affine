#include<cmath>
#include<sstream>

#include "../include/numberCountsKnotsSpline.h"
#include "../include/global_settings.h"
#include "../include/affineExcept.h"

//Function to pass to GSL integrator
static double evalPowfNKnotsSpline(double, void*); //!< Evaluates f^pow dN/dS

numberCountsKnotsSpline::numberCountsKnotsSpline() : numberCountsKnots() {
  logknots = NULL;
  acc = gsl_interp_accel_alloc();
  splinelog = NULL;
  gsl_work = gsl_integration_workspace_alloc(1000);
  varr = new void*[5];
}

/*!
  \param[in] NKNOTS Number of knots
*/
numberCountsKnotsSpline::numberCountsKnotsSpline(unsigned int NKNOTS) :
  numberCountsKnots(NKNOTS) {
  if (NKNOTS > 0) logknots = new double[NKNOTS]; else logknots = NULL;
  acc = gsl_interp_accel_alloc();
  if (NKNOTS > 0)
    splinelog = gsl_spline_alloc(gsl_interp_cspline,
				 static_cast<size_t>(NKNOTS));
  else {
    splinelog=NULL;
  }
  gsl_work = gsl_integration_workspace_alloc(1000);
  varr = new void*[5];
}

/*!
  \param[in] v Vector of knot positions
*/
numberCountsKnotsSpline::numberCountsKnotsSpline(const std::vector<float>& v):
  numberCountsKnots(v) {
  if (nknots > 0) {
    logknots = new double[nknots]; 
    for (unsigned int i = 0; i < nknots; ++i)
      logknots[i] = log2(static_cast<double>(knots[i]));
  } else logknots = NULL;
  acc = gsl_interp_accel_alloc();
  if (v.size() > 0)
    splinelog = gsl_spline_alloc(gsl_interp_cspline,
				 static_cast<size_t>(v.size()));
  else
    splinelog = NULL;
  gsl_work = gsl_integration_workspace_alloc(1000);
  varr = new void*[5];
}

/*!
  \param[in] v Vector of knot positions
*/
numberCountsKnotsSpline::numberCountsKnotsSpline(const std::vector<double>& v):
  numberCountsKnots(v) {
  if (nknots > 0) {
    logknots = new double[nknots]; 
    for (unsigned int i = 0; i < nknots; ++i)
      logknots[i] = log2(knots[i]);
  } else logknots = NULL;
  acc = gsl_interp_accel_alloc();
  if (v.size() > 0)
    splinelog = gsl_spline_alloc(gsl_interp_cspline,
				 static_cast<size_t>(v.size()));
  else
    splinelog = NULL;
  gsl_work = gsl_integration_workspace_alloc(1000);
  varr = new void*[5];
}

/*!
  \param[in] n Number of knots
  \param[in] S Array of knot positions.  Must be at least length n.
*/
numberCountsKnotsSpline::numberCountsKnotsSpline(unsigned int n, 
						 const float* const S) :
  numberCountsKnots(n, S) {
  if (nknots > 0) {
    logknots = new double[nknots]; 
    for (unsigned int i = 0; i < nknots; ++i)
      logknots[i] = log2(static_cast<double>(knots[i]));
  } else logknots = NULL;
  acc = gsl_interp_accel_alloc();
  splinelog = gsl_spline_alloc(gsl_interp_cspline,
			       static_cast<size_t>(n));
  gsl_work = gsl_integration_workspace_alloc(1000);
  varr = new void*[5];
}

/*!
  \param[in] n Number of knots
  \param[in] S Array of knot positions.  Must be at least length n.
*/
numberCountsKnotsSpline::numberCountsKnotsSpline(unsigned int n, 
						 const double* const S) :
  numberCountsKnots(n, S) {
  if (nknots > 0) {
    logknots = new double[nknots]; 
    for (unsigned int i = 0; i < nknots; ++i)
      logknots[i] = log2(knots[i]);
  } else logknots = NULL;
  acc = gsl_interp_accel_alloc();
  splinelog = gsl_spline_alloc(gsl_interp_cspline,
			       static_cast<size_t>(n));
  gsl_work = gsl_integration_workspace_alloc(1000);
  varr = new void*[5];
}

/*!
  \param[in] other numberCountsKnotsSpline to copy
*/
numberCountsKnotsSpline::
numberCountsKnotsSpline(const numberCountsKnotsSpline& other){
  if (this == &other) return; //Self-copy
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
    gsl_spline_init(splinelog, logknots, logknotvals,
		    static_cast<size_t>(other.nknots));
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

/*!
  \param[in] n New number of knots

  Does nothing if same as the current number.  Otherwise contents
  are destroyed.
*/
void numberCountsKnotsSpline::setNKnots(unsigned int n) {
  if (nknots == n) return;
  if (knots != NULL) delete[] knots;
  if (logknotvals != NULL) delete[] logknotvals;
  if (logknots != NULL) delete[] logknots;
  if (splinelog != NULL) gsl_spline_free(splinelog);

  if (n > 0) {
    knots = new double[n];
    logknots = new double[n];
    logknotvals = new double[n];
  } else {
    knots = logknots = logknotvals = NULL;
  }
  if (n > 1)
    splinelog=gsl_spline_alloc(gsl_interp_cspline,
			       static_cast<size_t>(n));
  else 
    splinelog=NULL;
  nknots = n;
  knotvals_loaded = false;
}

/*!
  \param[in] other numberCountsKnotsSpline instance to copy
  \returns Reference to this object
*/
numberCountsKnotsSpline& 
numberCountsKnotsSpline::operator=(const numberCountsKnotsSpline& other) {
  if (this == &other) return *this; //Self-copy
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
      gsl_spline_init(splinelog, logknots, logknotvals,
		      static_cast<size_t>(other.nknots));
    }
    knotvals_loaded = other.knotvals_loaded;
  }
  return *this;
}
 
/*!
  \param[in] vec Vector of knot positions.  

  This will change the number of knots if needed.
*/
void numberCountsKnotsSpline::setKnotPositions(const std::vector<float>& vec) {
  numberCountsKnots::setKnotPositions(vec);
  for (unsigned int i = 0; i < nknots; ++i)
    logknots[i] = log2(knots[i]);
}

/*!
  \param[in] vec Vector of knot positions.  

  This will change the number of knots if needed.
*/
void numberCountsKnotsSpline::setKnotPositions(const std::vector<double>& vec) {
  numberCountsKnots::setKnotPositions(vec);
  for (unsigned int i = 0; i < nknots; ++i)
    logknots[i] = log2(knots[i]);
}

/*!
  \param[in] n Number of elements in vals
  \param[in] vals Array of knot positions.  

  This will change the number of knots if needed.
*/
void numberCountsKnotsSpline::setKnotPositions(unsigned int n, 
					       const float* const vals) {
  numberCountsKnots::setKnotPositions(n, vals);
  for (unsigned int i = 0; i < nknots; ++i)
    logknots[i] = log2(knots[i]);
}

/*!
  \param[in] n Number of elements in vals
  \param[in] vals Array of knot positions.  

  This will change the number of knots if needed.
*/
void numberCountsKnotsSpline::setKnotPositions(unsigned int n, 
					       const double* const vals) {
  numberCountsKnots::setKnotPositions(n, vals);
  for (unsigned int i = 0; i < nknots; ++i)
    logknots[i] = log2(knots[i]);
}


/*!
  \param[in] F Parameters.  Will ignore any after
                    the first nknots
*/
void numberCountsKnotsSpline::setParams(const paramSet& F) {
  if (nknots > F.getNParams())
    throw affineExcept("numberCountsKnotsSpline", "setKnots",
		       "Number of knot values different than expected");
  for (unsigned int i = 0; i < nknots; ++i)
    logknotvals[i] = pofd_mcmc::logfac * static_cast<double>(F[i]);
  gsl_spline_init(splinelog, logknots, logknotvals,
		  static_cast<size_t>(nknots));
  knotvals_loaded = true;
}


/*!
  \param[in] fdens Flux density to evaluate number counts at
  \returns Number counts evaluated at fdens, or NaN if model is
    not valid
  
  Does input validity checks.
*/
double numberCountsKnotsSpline::getNumberCounts(double fdens) const {
  if (nknots < 2) return std::numeric_limits<double>::quiet_NaN();
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  if (fdens < knots[0] || fdens >= knots[nknots-1]) return 0.0; //Out of range
  return exp2(gsl_spline_eval(splinelog, log2(fdens), acc)); 
}

/*!
  \param[in] alpha  Power of flux
  \returns Integral over number counts

  Computes
  \f[
   \int dS\, S^{\alpha} \frac{dN}{dS}
  \f]
*/
double numberCountsKnotsSpline::splineInt(double alpha) const {
  if (nknots < 2) return std::numeric_limits<double>::quiet_NaN();
  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
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

/*!
  \returns Number of sources per unit area
*/
double numberCountsKnotsSpline::getNS() const {
  return splineInt(0);
}

/*!
  \returns Total flux density per unit area
*/
double numberCountsKnotsSpline::getFluxPerArea() const {
  return splineInt(1);
}

/*!
  \returns Total flux density squared per unit area
*/
double numberCountsKnotsSpline::getFluxSqPerArea() const {
  return splineInt(2);
}
  
// Pure brute force R computation
/*
  \params[in] fluxdensity  Position to get R at
  \params[in] params Parameters
  \params[in] beam Beam
  \returns R(fluxdensity)
*/
double numberCountsKnotsSpline::getR(double fluxdensity, const beam& bm) const {

  if (!isValid()) return std::numeric_limits<double>::quiet_NaN();
  if (!bm.hasData()) 
    throw affineExcept("numberCountsKnotsSpline", "getR",
		       "Beam has no data");

  // Quick return
  if (fluxdensity > 0 && !bm.hasPos()) return 0.0;
  else if (fluxdensity < 0 && !bm.hasNeg()) return 0.0;
  else if (fluxdensity == 0) return 0.0;

  double minknot = knots[0];
  double maxknot = knots[nknots - 1];

  // Find range where R is nonzero, if we are outside it easy
  dblpair Rrange = getRRangeInternal(bm);
  if (fluxdensity >= Rrange.second) return 0.0;
  if (fluxdensity <= Rrange.first) return 0.0;

  double prefac;
  prefac = bm.getPixSize() / 3600.0;  //To sq deg
  prefac = prefac * prefac;
  
  double currarg, ieta, retval, cts;
  retval = 0.0;
  if (fluxdensity > 0) { // Note that R(0) = 0, so we just test > <
    if (bm.isPosHist()) {
      // Histogrammed beam
      unsigned int nposhist = bm.getNHistPos();
      const double* warr = bm.getPosHistWeights();
      const double *ipixarr = bm.getPosHist();
      for (unsigned int i = 0; i < nposhist; ++i) {
	ieta = ipixarr[i];
	currarg = fluxdensity * ieta;
	if (currarg < minknot || currarg >= maxknot) continue;
	cts = exp2(gsl_spline_eval(splinelog, log2(currarg), acc));
	retval += warr[i] * ieta * cts;
      }
    } else {
      // Raw beam
      unsigned int npsf = bm.getNPos();
      const double *ipixarr = bm.getPosInvPixArr();
      for (unsigned int i = 0; i < npsf; ++i) {
	ieta = ipixarr[i];
	currarg = fluxdensity * ieta;
	if (currarg < minknot || currarg >= maxknot) continue;
	retval += exp2(gsl_spline_eval(splinelog, log2(currarg), acc)) * ieta;
      }
    }
  } else if (fluxdensity < 0) {
    if (bm.isNegHist()) {
      unsigned int nneghist = bm.getNHistNeg();
      const double* warr = bm.getNegHistWeights();
      const double *ipixarr = bm.getNegHist();
      for (unsigned int i = 0; i < nneghist; ++i) {
	ieta = ipixarr[i];
	currarg = - fluxdensity * ieta; // Note the -
	if (currarg < minknot || currarg >= maxknot) continue;
	cts = exp2(gsl_spline_eval(splinelog, log2(currarg), acc));
	retval += warr[i] * ieta * cts;
      }
    } else {
      unsigned int npsf = bm.getNNeg();
      const double *ipixarr = bm.getNegInvPixArr();
      for (unsigned int i = 0; i < npsf; ++i) {
	ieta = ipixarr[i];
	currarg = - fluxdensity * ieta;
	if (currarg < minknot || currarg >= maxknot) continue;
	retval += exp2(gsl_spline_eval(splinelog, log2(currarg), acc)) * ieta;
      }
    }
  }

  return prefac * retval;
}

/*!
  \param[in] n    Number of elements to get R for (length of flux and R)
  \param[in] flux Fluxes to get R for.  Note we don't assume these are monotonic
  \param[in] bm   Beam
  \param[inout]   R Computed value of R.  Must be pre-allocated to length n
                   (or more) by caller.  Only the first n values are set.
*/
void numberCountsKnotsSpline::getR(unsigned int n, const double* const flux,
				   const beam& bm, double* R) const {
  if (n == 0) return;
  if (!isValid()) {
    for (unsigned int i = 0; i < n; ++i)
      R[i] = std::numeric_limits<double>::quiet_NaN();
    return;
  }

  // Find range where R is nonzero, if we are outside it easy
  dblpair Rrange = getRRangeInternal(bm);

  double minknot = knots[0];
  double maxknot = knots[nknots-1];

  double prefac;
  prefac = bm.getPixSize() / 3600.0;  //To sq deg
  prefac = prefac * prefac;

  // Precalculate and store local refs to all the arrays we may need
  unsigned int npos = 0, nneg = 0;
  bool haspos = bm.hasPos();
  bool hasneg = bm.hasNeg();
  bool poshist = false, neghist = false;
  const double* wtptr_pos = NULL;
  const double* wtptr_neg = NULL;
  const double* ibmptr_pos = NULL;
  const double* ibmptr_neg = NULL;
  if (haspos) {
    poshist = bm.isPosHist();
    npos = poshist ? bm.getNHistPos() : bm.getNPos();
    if (poshist) wtptr_pos = bm.getPosHistWeights();
    ibmptr_pos = poshist ? bm.getPosHist() : bm.getPosInvPixArr();
  }
  if (hasneg) {
    neghist = bm.isNegHist();
    nneg = neghist ? bm.getNHistNeg() : bm.getNNeg();
    if (neghist) wtptr_neg = bm.getNegHistWeights();
    ibmptr_neg = neghist ? bm.getNegHist() : bm.getNegInvPixArr();
  }

  double currarg, ieta, workval, fluxdensity, cts;
  for (unsigned int i = 0; i < n; ++i) {
    fluxdensity = flux[i];
    workval = 0.0;
    // Again -- R[0] = 0, so we can just test > < 
    // There is a three part test here -- first on the sign of the
    //  fluxdensity, then if the corresponding beam is available,
    //  and then, if those are satisfied, if this fluxdensity is
    //  outside the range where the contribution is non-zero
    if (fluxdensity > 0 && haspos && fluxdensity <= Rrange.second) { 
      if (poshist) {
	for (unsigned int j = 0; j < npos; ++j) {
	  ieta = ibmptr_pos[j];
	  currarg = fluxdensity * ieta;
	  if (currarg < minknot || currarg >= maxknot) continue;
	  cts = exp2(gsl_spline_eval(splinelog, log2(currarg), acc));
	  if (cts > 0) workval += wtptr_pos[j] * cts * ieta;
	}
      } else {
	for (unsigned int j = 0; j < npos; ++j) {
	  ieta = ibmptr_pos[j];
	  currarg = fluxdensity * ieta;
	  if (currarg < minknot || currarg >= maxknot) continue;
	  cts = exp2(gsl_spline_eval(splinelog, log2(currarg), acc));
	  if (cts > 0) workval += cts * ieta;
	}
      }
    } else if (hasneg && fluxdensity < 0 && fluxdensity >= Rrange.first) {
      // Note the order is flopped on the first two tests because
      //  not having the negative beam is a common case, which is not
      //  true for not having a positive beam.  
      if (neghist) {
	for (unsigned int j = 0; j < nneg; ++j) {
	  ieta = ibmptr_neg[j];
	  currarg = - fluxdensity * ieta; // Note the - sign
	  if (currarg < minknot || currarg >= maxknot) continue;
	  cts = exp2(gsl_spline_eval(splinelog, log2(currarg), acc));
	  if (cts > 0) workval += wtptr_neg[j] * cts * ieta;
	}
      } else {
	for (unsigned int j = 0; j < nneg; ++j) {
	  ieta = ibmptr_neg[j];
	  currarg = - fluxdensity * ieta;
	  if (currarg < minknot || currarg >= maxknot) continue;
	  cts = exp2(gsl_spline_eval(splinelog, log2(currarg), acc));
	  if (cts > 0) workval += cts * ieta;
	}
      }
    }
    R[i] = prefac * workval;
  }
}

/*!
  \param[in] alpha Regularization multiplier.  Must be positive (not checked)
  \returns log Likelhood penalty

  Computes difference operator Tikhonov regularization penalty
  on model, where the derivative is taken in log/log space
 */
double numberCountsKnotsSpline::differenceRegularize(double alpha) const {
  if (!knotvals_loaded) return std::numeric_limits<double>::quiet_NaN();
  if (nknots == 0 || nknots == 1) return 0.0;

  double log_penalty = 0.0;
  double deriv, delta_logknotpos, delta_logknotval;
  for (unsigned int i = 1; i < nknots; ++i) {
    delta_logknotpos = logknots[i] - logknots[i - 1];
    delta_logknotval = logknotvals[i] - logknotvals[i - 1];
    deriv = delta_logknotval / delta_logknotpos;
    log_penalty -= deriv * deriv;
  }
  return alpha * log_penalty;
}

/*!
  \param[inout] objid HDF5 handle to write to
*/
void numberCountsKnotsSpline::writeToHDF5Handle(hid_t objid) const {
  hsize_t adims;
  hid_t mems_id, att_id;

  if (H5Iget_ref(objid) < 0)
    throw affineExcept("numberCountsKnotsSpline", "writeToHDF5Handle",
		       "Input handle is not valid");

  // Name of model
  const char modeltype[] = "numberCountsKnotsSpline";
  hid_t datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, strlen(modeltype)); 
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate1(objid, "model_type", datatype,
		      mems_id, H5P_DEFAULT);
  H5Awrite(att_id, datatype, modeltype);
  H5Aclose(att_id);
  H5Sclose(mems_id);
  
  // Other writes
  numberCountsKnots::writeToHDF5Handle(objid);
}

/*!
  \param[in] comm Communicator
  \param[in] dest Destination of messages
*/
void numberCountsKnotsSpline::sendSelf(MPI_Comm comm, int dest) const {
  numberCountsKnots::sendSelf(comm, dest);
  if (nknots > 0)
    MPI_Send(logknots, nknots, MPI_DOUBLE, dest, 
	     pofd_mcmc::NCKSSENDLOGKNOTS, comm);
}

/*!
  \param[in] comm Communicator
  \param[in] src Source of messages
*/
void numberCountsKnotsSpline::recieveCopy(MPI_Comm comm, int src) {
  MPI_Status Info;
  unsigned int oldnknots = nknots;
  numberCountsKnots::recieveCopy(comm, src);

  if (nknots != oldnknots) {
    if (splinelog != NULL) gsl_spline_free(splinelog);
    if (nknots > 0) 
      splinelog = gsl_spline_alloc(gsl_interp_cspline,
				   static_cast<size_t>(nknots));
    else
      splinelog = NULL;
  }
  if (nknots > 0) {
    MPI_Recv(logknots, nknots, MPI_DOUBLE, src, pofd_mcmc::NCKSSENDLOGKNOTS,
	     comm, &Info);
    if (knotvals_loaded)
      gsl_spline_init(splinelog, logknots, logknotvals,
		      static_cast<size_t>(nknots));
  }
}

/*!
  \param[in] x Flux density to evaluate at
  \param[in] params Array containing information about spline
  \returns x^power times the number counts evaluated at x
*/
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

  double splval = exp2(gsl_spline_eval(spl,log2(x),acc));

  if (fabs(power) < 1e-2) return splval;
  if (fabs(power-1.0) < 1e-2) return x*splval;
  if (fabs(power-2.0) < 1e-2) return x*x*splval;
  return pow(x, power) * splval;

}
