//PD.cc
#include<limits>
#include<sstream>

#include<fitsio.h>
#include<fftw3.h>
#include<hdf5.h>

#include "../include/global_settings.h"
#include "../include/PD.h"
#include "../include/affineExcept.h"

const double PD::lowsigval = 2.5;

/*!
  \param[in] N Number of elements
  \param[in] MINFLUX minimum flux value
  \param[in] DFLUX Delta flux value
*/
PD::PD(unsigned int N, double MINFLUX, double DFLUX) :
  n(N), capacity(N), logflat(false), minflux(MINFLUX), dflux(DFLUX) {
  if (capacity == 0) pd_ = NULL; else
    pd_ = (double *) fftw_malloc(sizeof(double) * capacity);
}

PD::~PD() {
  if (pd_ != NULL) fftw_free(pd_);
}

/*!
  \param[in] N New number of elements

  Does nothing if new size is smaller than old capacity.
  Otherwise doesn't preserve data.
*/
void PD::resize(unsigned int N) {
  //Doesn't actually resize arrays if it can avoid it
  unsigned int newcap = N;
  if (newcap > capacity) {
    if (pd_ != NULL) fftw_free(pd_);
    if (newcap > 0) pd_ = (double *) fftw_malloc(sizeof(double) * newcap);
    else pd_ = NULL;
    capacity = newcap;
  }
  n = N;
}

/*!
  Tries to preserve data
*/
void PD::shrink() {
  unsigned int newcap = n;
  if (newcap < capacity) {
    if (newcap > 0) {
      double* tmp = (double*) fftw_malloc(sizeof(double) * newcap);
      for (unsigned int i = 0; i < newcap; ++i)
	tmp[i] = pd_[i];
      if (pd_ != NULL) fftw_free(pd_);
      pd_ = tmp;
    } else {
      if (pd_ != NULL) fftw_free(pd_);
      pd_ = NULL;
    }
    capacity = newcap;
  }
}

/*!
  \param[in] N New number of elements

  Doesn't preserve data unless new size is the same as old
*/
void PD::strict_resize(unsigned int N) {
  if (N != capacity) {
    if (pd_ != NULL) fftw_free(pd_);
    if (N > 0) pd_ = (double*) fftw_malloc(sizeof(double) * N);
    else pd_ = NULL;
    capacity = N;
  }
  n = N;
}

/*!
  \returns Total of all elements
*/
double PD::getTotal() const {
  if (n == 0)
    return std::numeric_limits<double>::quiet_NaN();
  double retval;
  if (logflat) {
    retval = exp2(pd_[0]);
    for (unsigned int i = 1; i < n; ++i)
      retval += exp2(pd_[i]);
  } else {
    retval = pd_[0];
    for (unsigned int i = 1; i < n; ++i)
      retval += pd_[i];
  }
  return retval;
}

/*!
  \returns Integral of P(D)

  Uses the trapezoidal rule
*/
double PD::getIntegral() const {
  if (n == 0)
    return std::numeric_limits<double>::quiet_NaN();
  
  double tot;
  if (logflat) {
    tot = 0.5*exp2(pd_[0]);
    for (unsigned int i = 1; i < n-1; ++i)
      tot += exp2(pd_[i]);
    tot += 0.5*exp2(pd_[n-1]);
  } else {
    tot = 0.5*pd_[0];
    for (unsigned int i = 1; i < n-1; ++i)
      tot += pd_[i];
    tot += 0.5*pd_[n-1];
  }
  return tot*dflux;
}

/*!
  Normalize the P(D), using the trapezoidal rule
  to integrate
*/
void PD::normalize() {
  if (n == 0)
    throw affineExcept("PD", "normalize", 
		       "No information present to normalize");
  double tot = getIntegral();
  unsigned int sz = n;
  if (logflat) {
    double lgtot = log2(tot);
    for (unsigned int i = 0; i < sz; ++i)
      pd_[i] -= lgtot;
  } else {
    double itot = 1.0/tot;
    for (unsigned int i = 0; i < sz; ++i)
      pd_[i] *= itot;
  }
}

/*!
  \param[in] nocheck Don't check whether or not values are positive

  If doing check and value is less than or equal to zero, set to 
  pofd_mcmc::smalllogval
*/
void PD::applyLog(bool nocheck) {
  if (logflat) return;
  unsigned int sz = n;
  if (nocheck) {
    for (unsigned int i = 0; i < sz; ++i) 
      pd_[i] = log2(pd_[i]);
  } else {
    double val;
    for (unsigned int i = 0; i < sz; ++i) {
      val = pd_[i];
      if (val > 0.0) pd_[i] = log2(val);
      else pd_[i] = pofd_mcmc::smalllogval;
    }
  }
  logflat = true;
}

void PD::deLog() {
  if (!logflat) return;
  unsigned int sz = n;
  for (unsigned int i = 0; i < sz; ++i)
    pd_[i] = exp2(pd_[i]);
  logflat = false;
}

/*
  \param[in] donorm Do not assume P(D) is normalized
*/
void PD::edgeFix(bool donorm) {
  //Compute mean and stdev
  if (n < 3) return; //No point

  if (logflat)
    throw affineExcept("PD", "edgeFix", "Not supported for logged PDs");

  //Get mean and vars
  dblpair mnvar = getMeanAndVar(donorm);
  if (std::isnan(mnvar.first) || std::isinf(mnvar.first) ||
      std::isnan(mnvar.second) || std::isinf(mnvar.second)) {
    std::stringstream errstr;
    errstr << "Problem with mean/var of P(D): " << std::endl;
    if (std::isnan(mnvar.first)) errstr << std::endl << "Mean is NaN";
    if (std::isinf(mnvar.first)) errstr << std::endl << "Mean is Inf";
    if (std::isnan(mnvar.second)) errstr << std::endl << "Var is NaN";
    if (std::isinf(mnvar.second)) errstr << std::endl << "Var is Inf";
    throw affineExcept("PD", "edgeFix", errstr.str());
  }
  double istdev = 1.0 / sqrt(mnvar.second);

  //Figure out what indexes these represent in x and y
  double maxfluxfix;
  int maxidx;
  maxfluxfix = mnvar.first - PD::lowsigval*sqrt(mnvar.second);
  maxidx = static_cast<int>((maxfluxfix - minflux) / dflux);
  maxfluxfix = minflux + maxidx * dflux;
  
  double pdval, tval, preconst, stepfac, subfac;
  if (maxidx > 1) {
    pdval = pd_[maxidx];
    tval = (maxfluxfix - mnvar.first) * istdev;
    preconst = pdval * exp(0.5 * tval * tval);
    subfac = (minflux - mnvar.first) * istdev;
    stepfac = dflux * istdev;
    for (int i = 0; i < maxidx; ++i) {
      tval = subfac + i * stepfac;
      pd_[i] = preconst * exp(-0.5 * tval * tval);
    }
  }
}


/*!
  \param[in] donorm Do not assume that P(D) is normalized.

  \returns Mean of the P(D).
*/
double PD::getMean(bool donorm) const {

  if (n == 0) 
    return std::numeric_limits<double>::quiet_NaN();
  
  //We use the trapezoidal rule here for the integrals
  // so it isn't quite a simple sum
  double mean = 0.0;

  if (logflat) {
    //i=0 has weight 0
    for (unsigned int i = 1; i < n-1; ++i)
      mean += static_cast<double>(i) * exp2(pd_[i]);
    mean += 0.5 * static_cast<double>(n-1) * exp2(pd_[n-1]);
  } else {
    for (unsigned int i = 1; i < n-1; ++i)
      mean += static_cast<double>(i) * pd_[i];
    mean += 0.5 * static_cast<double>(n-1) * pd_[n-1];
  }
  //Add on step sizes for each integral,
  // which is both area and step size in x,y
  //dflux twice -- once for the conversion from i to flux, once for
  // the step size
  mean *= dflux*dflux;

  if (donorm) mean /= getIntegral();

  mean += minflux;
  return mean;
}

/*!
  \param[in] donorm Do not assume that P(D) is normalized.
  \returns A pair containing the mean and variance.
*/
dblpair PD::getMeanAndVar(bool donorm) const {
  if (n == 0) 
    return std::make_pair(std::numeric_limits<double>::quiet_NaN(),
			  std::numeric_limits<double>::quiet_NaN());
  
  double normfac = 1.0;
  if (donorm) normfac = 1.0 / getIntegral();
  
  //It is more efficient to do this as a <x^2> - <x>^2 calculation
  //The problem is that that way results on fine cancellations.
  //So we do this in two steps -- get the mean, then compute the
  // variance.  This is a lot more expensive, especially if the
  // PD is stored as a log, but is more accurate.
  //Why not just call getMeans?  To avoid calling getIntegral
  // twice.  After this, mean1 and mean2 will be the actual
  // means/dflux - minflux.
  double mean = 0.0;
  if (logflat) {
    //i=0 has weight 0
    for (unsigned int i = 1; i < n-1; ++i)
      mean += static_cast<double>(i) * exp2(pd_[i]);
    //0.5 factor since last step of trap
    mean += 0.5 * static_cast<double>(n-1) * exp2(pd_[n-1]); 
  } else {
    for (unsigned int i = 1; i < n-1; ++i)
      mean += static_cast<double>(i) * pd_[i];
    mean += 0.5 * static_cast<double>(n-1) * pd_[n-1]; //0.5 since last trap
  }

  //We need two factors of dflux -- one for xval=dflux*i and one for the
  // step in the integration
  mean   *= dflux * dflux;

  //Apply normalization
  if (donorm) mean *= normfac;

  //Now, compute the variance
  double var = 0.0;
  double xval;
  if (logflat) {
    //i=0 has weight 0
    for (unsigned int i = 1; i < n-1; ++i) {
      xval = static_cast<double>(i)*dflux - mean;
      var += xval * xval * exp2(pd_[i]);
    }
    xval = static_cast<double>(n-1)*dflux - mean;
    var += 0.5 * xval * xval * exp2(pd_[n-1]);
  } else {
    for (unsigned int i = 1; i < n-1; ++i) {
      xval = static_cast<double>(i)*dflux - mean;
      var += xval * xval * pd_[i];
    }
    xval = static_cast<double>(n-1)*dflux - mean;
    var += 0.5 * xval * xval * pd_[n-1];
  }

  //This time only one factor of dflux, for the integral step
  var *= dflux;
  if (donorm) var *= normfac;

  //Correct mean for offset
  mean += minflux;

  return std::make_pair(mean, var);
}

/*!
  \param[in] other PD to copy from
*/
PD& PD::operator=(const PD& other) {
  if (this == &other) return *this; //Self-copy
  resize(other.n);
  minflux = other.minflux;
  dflux   = other.dflux;
  if (n > 0) memcpy(pd_, other.pd_, n * sizeof(double));
  logflat = other.logflat;
  return *this;
}

/*!
  \param[in] N New number of elements
  \param[in] MINFLUX New minimum flux density
  \param[in] DFLUX New delta flux density
  \param[in] PD New values
  \param[in] LOG Is the input PD log2?
*/
void PD::fill(unsigned int N, double MINFLUX, double DFLUX,
	      const double* const PD, bool LOG) {
  logflat = LOG;
  resize(N);
  minflux = MINFLUX;
  dflux = DFLUX;
  if (n > 0) memcpy(pd_, PD, n * sizeof(double));
}

/*!
  \param[in] x Flux density to evaluate P(D) at
  \returns Interpolated P(D) value.  Will be log2 if P(D) is
    stored as a log.

  Uses linear interpolation.  Note that this will work differently
  if the P(D) has been converted to log values.
*/
double PD::getPDVal(double x) const {
  if (pd_ == NULL) return std::numeric_limits<double>::quiet_NaN();

  //look up the effective indexes
  int idx = static_cast<int>((x - minflux) / dflux);

  double maxfluxlm = minflux + static_cast<double>(n-2)*dflux;

  if (x < minflux) return pd_[0];
  if (x > maxfluxlm) return pd_[n-1]; 
  double p0 = pd_[idx];
  double dx = x - (minflux + static_cast<double>(idx) * dflux);
  return p0 + dx / dflux * (pd_[idx + 1] - p0);
}

/*!
  \param[inout] os Stream to output values to
*/
std::ostream& PD::writeToStream(std::ostream& os) const {
  os << n << " " << minflux << " " << dflux << std::endl;
  if (n > 0) {
    os << pd_[0];
    for (unsigned int i = 1; i < n-1; ++i)
      os << " " << pd_[i];
    os << std::endl;
  }
  return os;
}

/*!
  \param[in] outputfile File to write to in FITS format
  \returns 0 on success, an error code (!=0) for anything else
*/
int PD::writeToFits(const std::string& outputfile) const {

  //Make the fits file
  int status = 0;
  fitsfile *fp;

  fits_create_file(&fp, outputfile.c_str(), &status);

  if (status) {
    fits_report_error(stderr,status);
    throw affineExcept("PD", "writeToFits", "Error creating FITS output file");
  }

  long axissize[1];
  axissize[0] = static_cast<long>(n);
  fits_create_img(fp, DOUBLE_IMG, 1, axissize, &status);
  
  //Add "WCS" info to hdr
  float crpix = 1;
  double tmpval;
  fits_write_key(fp, TSTRING, const_cast<char*>("CTYPE1"),
		 const_cast<char*>("FLUX"),
		 const_cast<char*>("Type of Data axis 1"),&status);
  fits_write_key(fp, TFLOAT, const_cast<char*>("CRPIX1"), &crpix, 
		 const_cast<char*>("Ref pix of axis 1"), &status);
  tmpval = minflux;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CRVAL1"), &tmpval, 
		 const_cast<char*>("val at ref pix"), &status);
  tmpval = dflux;
  fits_write_key(fp, TDOUBLE, const_cast<char*>("CDELT1"), &tmpval,
		 const_cast<char*>("delta along axis 1"), &status);

  int lg = static_cast<int>(logflat);
  fits_write_key(fp, TLOGICAL, const_cast<char*>("LOG"),&lg,
		 const_cast<char*>("Is log P(D) stored?"), &status);

  //Do data writing.  
  long fpixel[1] = { 1 };
  fits_write_pix(fp, TDOUBLE, fpixel, n, pd_, &status);

  fits_close_file(fp, &status);

  if (status) {
    fits_report_error(stderr, status);
    throw affineExcept("PD", "writeToFits", "Error doing FITS write");
  }
  return status;
}

/*!
  \param[in] outputfile File to write to
*/
void PD::writeToHDF5(const std::string& outputfile) const {
  hid_t file_id;
  file_id = H5Fcreate(outputfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
		      H5P_DEFAULT);
  if (H5Iget_ref(file_id) < 0) {
    H5Fclose(file_id);
    throw affineExcept("PD", "writeToHDF5", 
		       "Failed to open HDF5 file to write");
  }

  hsize_t adims;
  hid_t mems_id, att_id, dat_id;
  
  // Properties
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate2(file_id, "isLog", H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, &logflat);
  H5Aclose(att_id);
  att_id = H5Acreate2(file_id, "dflux", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &dflux);
  H5Aclose(att_id);
  att_id = H5Acreate2(file_id, "minflux", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_DOUBLE, &minflux);
  H5Aclose(att_id);
  H5Sclose(mems_id);
  
  // Rflux -- by making temporary array
  double *flux = new double[n];
  for (unsigned int i = 0; i < n; ++i) 
    flux[i] = static_cast<double>(i) * dflux + minflux;
  adims = n;
  mems_id = H5Screate_simple(1, &adims, NULL);
  dat_id = H5Dcreate2(file_id, "flux", H5T_NATIVE_DOUBLE,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
	   H5P_DEFAULT, flux);
  H5Dclose(dat_id);
  delete[] flux;

  dat_id = H5Dcreate2(file_id, "PD", H5T_NATIVE_DOUBLE, mems_id,
		      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dat_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
	   H5P_DEFAULT, pd_);
  H5Dclose(dat_id);
  H5Sclose(mems_id);

  H5Fclose(file_id);
}

/*!
  \param[in] data Data to evaluate the log Likelihood of
  \returns Log likelihood
*/
double PD::getLogLike(const fitsData& data) const {
  if (pd_ == NULL) throw affineExcept("PD", "getLogLike",
				      "pd not filled before likelihood calc");
  unsigned int ndata = data.getN();
  if (ndata == 0) throw affineExcept("PD", "getLogLike",
				     "No data present");

  if (data.isBinned()) return getLogLikeBinned(data);
  else return getLogLikeUnbinned(data);
}

/*!
  \param[in] data Data to evaluate the log Likelihood of
  \returns Log likelihood

  Because this uses interpolation, slightly different values
  will result if the P(D) is stored in log form or not.
*/
double PD::getLogLikeBinned(const fitsData& data) const {

  if (!data.isBinned())
    throw affineExcept("PD", "getLogLikeBinned", "Data is not binned");

  //Quantities for edge test
  double maxflux = minflux + static_cast<double>(n-1)*dflux;

  int idx; //!< Index look up
  const unsigned int* bins;
  unsigned int nbins, ninbin;
  double cflux, bincent0, bindelta, loglike;
  double t, delt, maxfluxval;
  loglike = 0.0;
  bins = data.getBinnedData();
  bincent0 = data.getBinCent0();
  bindelta = data.getBinDelta();
  nbins = data.getNBins();
  double idflux = 1.0 / dflux;
  maxfluxval = maxflux - dflux;
  if (logflat) {
    for (unsigned int i = 0; i < nbins; ++i) {
      ninbin = bins[i];
      if (ninbin == 0) continue;
      cflux = bincent0 + static_cast<double>(i) * bindelta;
      //Get effective index
      if (cflux < minflux) loglike += static_cast<double>(ninbin) * pd_[0];
      else if (cflux > maxfluxval) 
        loglike += static_cast<double>(ninbin)*pd_[n-1];
      else {
        //Not off edge
        delt = (cflux - minflux) * idflux;
        idx = static_cast<int>(delt);
        t   = delt - static_cast<double>(idx);
        loglike += ((1.0-t) * pd_[idx] + t * pd_[idx + 1]) *
          static_cast<double>(ninbin);
      }
    }
  } else {
    //Not stored as log -- inefficient, but supported
    //Note that it would be insane to do this multiplicatively,
    // then take the log.  Also, it's better to interpolate
    // in log space than interpolate, then log
    for (unsigned int i = 0; i < nbins; ++i) {
      ninbin = bins[i];
      if (ninbin == 0) continue;
      cflux = bincent0 + static_cast<double>(i)*bindelta;
      delt = (cflux - minflux) * idflux;
      idx = static_cast<int>(delt);
      if (cflux < minflux) 
        loglike += static_cast<double>(ninbin) * log2(pd_[0]);
      else if (cflux > maxfluxval) 
        loglike += static_cast<double>(ninbin) * log2(pd_[n - 1]);
      else {
        t   = delt - static_cast<double>(idx);
        loglike += ((1.0 - t) * log2(pd_[idx]) + t * log2(pd_[idx + 1])) *
          static_cast<double>(ninbin);
      }
    }
  }
  //This has been base 2 -- convert back to base e
  return pofd_mcmc::log2toe * loglike;
}

/*!
  \param[in] data Data to evaluate the log Likelihood of
  \returns Log likelihood

  Because this uses interpolation, slightly different values
  will result if the P(D) is stored in log form or not.
*/
double PD::getLogLikeUnbinned(const fitsData& data) const {

  unsigned int ndata = data.getN();

  //Quantities for edge test
  double maxflux = minflux + static_cast<double>(n-1)*dflux;

  int idx; //!< Index look up
  const double* flux;
  double cflux, loglike, t, delt, maxfluxval;
  loglike = 0.0;
  flux = data.getData();
  double idflux = 1.0/dflux;
  maxfluxval = maxflux - dflux;

  if (logflat) {
    //Stored as log2 P(D)
    for (unsigned int i = 0; i < ndata; ++i) {
      cflux = flux[i];
      //Get effective index
      if (cflux <= minflux) loglike += pd_[0];
      else if (cflux >= maxfluxval) loglike += pd_[n - 1];
      else {
        //Not off edge
        delt = (cflux - minflux)*idflux;
        idx = static_cast<int>(delt);
        t   = delt - static_cast<double>(idx);
        loglike += (1.0 - t) * pd_[idx] + t * pd_[idx + 1];
      }
    }
  } else {
    //Not stored as log -- inefficient, but supported
    //Note that it would be insane to do this multiplicatively,
    // then take the log.  Also, it's better to interpolate
    // in log space than interpolate, then log
    for (unsigned int i = 0; i < ndata; ++i) {
      cflux = flux[i]; 
      if (cflux <= minflux) loglike += log2(pd_[0]);
      else if (cflux >= maxfluxval) loglike += log2(pd_[n - 1]);
      else {
        delt = (cflux-minflux)*idflux;
        idx = static_cast<int>(delt);
        t   = delt - static_cast<double>(idx);
        loglike += (1.0-t)*log2(pd_[idx]) + t*log2(pd_[idx+1]);
      }
    }
  }
  //This has been base 2 -- convert back to base e
  return pofd_mcmc::log2toe * loglike;
}

/*!
  \param[inout] os Stream to write to
  \param[in] b PD to output
*/
std::ostream& operator<<(std::ostream& os, const PD& b) {
  b.writeToStream(os);
  return os;
}
