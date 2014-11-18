//beam.cc

#include<limits>
#include<iostream>
#include<cmath>
#include<sstream>
#include<cstring>

#include<fitsio.h>

#include "../include/beam.h"
#include "../include/global_settings.h"
#include "../include/affineExcept.h"

const unsigned int beam::histothresh = 15;

beam::beam() {
  nneg = npos = 0; 
  totneg = totpos = 0; 
  totnegsq = totpossq = 0.0;
  pixsize=0.0; 
  pospixarr = negpixarr = nullptr;
  posinvpixarr = neginvpixarr = nullptr;
  haspos = hasneg = false;
  minval = 0.0;
  nbins = 0;
  is_pos_histogrammed = is_neg_histogrammed = false;
  posnbins = negnbins = 0;
  posweights = negweights = nullptr;
  poshistval = neghistval = nullptr;
  
  if (minval < 0.0)
    throw affineExcept("beam", "beam", "Invalid (non-positive) minval");
}

/*!
  \param[in] filename   File to read from
  \param[in] histogram  Do beam histogramming
  \param[in] NBINS Number of histogram bins, if histogramming is applied
  \param[in] MINVAL Only beam components with an absoute value larger than
                    this are kept.
*/
beam::beam(const std::string& filename, bool histogram, 
	   unsigned int NBINS, double MINVAL) {
  nneg = npos = 0; 
  totneg = totpos = 0; 
  totnegsq = totpossq = 0.0;
  pixsize=0.0; 
  pospixarr = negpixarr = nullptr;
  posinvpixarr = neginvpixarr = nullptr;
  haspos = hasneg = false;
  is_pos_histogrammed = is_neg_histogrammed = false;
  posnbins = negnbins = 0;
  posweights = negweights = nullptr;
  poshistval = neghistval = nullptr;

  if (minval < 0.0)
    throw affineExcept("beam", "beam", "Invalid (non-positive) minval");

  readFile(filename, MINVAL);
  if (histogram) makeHistogram(NBINS);
}

/*!
  \param[in] inp Beam to copy
*/
beam::beam(const beam& inp) {
  negpixarr = pospixarr = nullptr;
  posinvpixarr = neginvpixarr = nullptr;
  posweights = negweights = nullptr;
  poshistval = neghistval = nullptr;
  npos = inp.npos;
  nneg = inp.nneg;
  haspos = inp.haspos;
  hasneg = inp.hasneg;
  is_pos_histogrammed = inp.is_pos_histogrammed;
  is_neg_histogrammed = inp.is_neg_histogrammed;
  nbins = inp.nbins;
  posnbins = inp.posnbins;
  negnbins = inp.negnbins;
  pixsize = inp.pixsize;
  totpos = inp.totpos;
  totneg = inp.totneg;
  totpossq = inp.totpossq;
  totnegsq = inp.totnegsq;
  minval = inp.minval;
  if (haspos) {
    pospixarr = new double[npos];
    for (unsigned int i = 0; i < npos; ++i) pospixarr[i]=inp.pospixarr[i];
    posinvpixarr = new double[npos];
    for (unsigned int i = 0; i < npos; ++i) 
      posinvpixarr[i]=inp.posinvpixarr[i];
    if (is_pos_histogrammed) {
      posweights = new double[posnbins];
      for (unsigned int i = 0; i < posnbins; ++i)
	posweights[i] = inp.posweights[i];
      poshistval = new double[posnbins];
      for (unsigned int i = 0; i < posnbins; ++i)
	poshistval[i] = inp.poshistval[i];
    }
  }
  if (hasneg) {
    negpixarr = new double[nneg];
    for (unsigned int i = 0; i < nneg; ++i) negpixarr[i]=inp.negpixarr[i];
    neginvpixarr = new double[nneg];
    for (unsigned int i = 0; i < nneg; ++i) 
      neginvpixarr[i]=inp.neginvpixarr[i];
    if (is_neg_histogrammed) {
      negweights = new double[negnbins];
      for (unsigned int i = 0; i < negnbins; ++i)
	negweights[i] = inp.negweights[i];
      neghistval = new double[negnbins];
      for (unsigned int i = 0; i < negnbins; ++i)
	neghistval[i] = inp.neghistval[i];
    }
  }
}

/*!
  \param[in] other beam to copy
*/
beam& beam::operator=(const beam& other) {
  if (this == &other) return *this;
  cleanup();
  npos = other.npos;
  nneg = other.nneg;
  haspos = other.haspos;
  hasneg = other.hasneg;
  nbins = other.nbins;
  is_pos_histogrammed = other.is_pos_histogrammed;
  is_neg_histogrammed = other.is_neg_histogrammed;
  posnbins = other.posnbins;
  negnbins = other.negnbins;
  pixsize = other.pixsize;
  totpos = other.totpos;
  totneg = other.totneg;
  totpossq = other.totpossq;
  totnegsq = other.totnegsq;
  minval = other.minval;
  if (haspos) {
    pospixarr = new double[npos];
    for (unsigned int i = 0; i < npos; ++i) pospixarr[i]=other.pospixarr[i];
    posinvpixarr = new double[npos];
    for (unsigned int i = 0; i < npos; ++i) 
      posinvpixarr[i]=other.posinvpixarr[i];
    if (is_pos_histogrammed) {
      posweights = new double[posnbins];
      for (unsigned int i = 0; i < posnbins; ++i) 
	posweights[i]=other.posweights[i];
      poshistval = new double[posnbins];
      for (unsigned int i = 0; i < posnbins; ++i) 
	poshistval[i]=other.poshistval[i];
    }
  }
  if (hasneg) {
    negpixarr = new double[nneg];
    for (unsigned int i = 0; i < nneg; ++i) negpixarr[i]=other.negpixarr[i];
    neginvpixarr = new double[nneg];
    for (unsigned int i = 0; i < nneg; ++i) 
      neginvpixarr[i]=other.neginvpixarr[i];
    if (is_neg_histogrammed) {
      negweights = new double[negnbins];
      for (unsigned int i = 0; i < negnbins; ++i) 
	negweights[i]=other.negweights[i];
      neghistval = new double[negnbins];
      for (unsigned int i = 0; i < negnbins; ++i) 
	neghistval[i]=other.neghistval[i];
    }
  }

  return *this;
}

void beam::cleanup() {
  if (pospixarr != nullptr) { delete[] pospixarr; pospixarr = nullptr; }
  if (negpixarr != nullptr) { delete[] negpixarr; negpixarr = nullptr; }
  if (posinvpixarr != nullptr) {
    delete[] posinvpixarr;
    posinvpixarr = nullptr;
  }
  if (neginvpixarr != nullptr) {
    delete[] neginvpixarr;
    neginvpixarr = nullptr;
  }
  if (posweights != nullptr) { delete[] posweights; posweights = nullptr; }
  if (negweights != nullptr) { delete[] negweights; negweights = nullptr; }
  if (poshistval != nullptr) { delete[] posweights; poshistval = nullptr; }
  if (neghistval != nullptr) { delete[] negweights; neghistval = nullptr; }
  haspos = hasneg = false;
  is_pos_histogrammed = is_neg_histogrammed = false;
  npos = nneg = 0;
  posnbins = negnbins = 0;
  pixsize = 0.0;
  totpos = totneg = 0.0;
  totpossq = totnegsq = 0.0;
}

bool beam::revSort(const double& d1, const double& d2) const {
  return d1 > d2;
}

/*!
  \param[in] filename   File to read from
  \param[in] MINVAL     Minimum beam value to use
*/
void beam::readFile(const std::string& filename, double MINVAL) {

  int status;
  fitsfile *fptr;
  char card[FLEN_CARD];
  int nkeys;

  cleanup(); //Free arrays, etc

  if (MINVAL <= 0.0)
    minval = 0.0;
  else
    minval = MINVAL;
  double negminval = -minval;

  // Do the read
  status = 0;
  fptr = nullptr;

  fits_open_file(&fptr, filename.c_str(), READONLY, &status);
  if (fptr == nullptr) {
    fits_close_file(fptr, &status);
    std::stringstream errstr;
    errstr << "Error opening input file: " << filename << std::endl;
    throw affineExcept("beam", "readFile", errstr.str());
  }
  if (status) {
    fits_report_error(stderr, status);
    fits_close_file(fptr, &status);
    std::stringstream errstr;
    errstr << "Error opening input file: " << filename << std::endl;
    throw affineExcept("beam", "readFile", errstr.str());
  }

  fits_get_hdrspace(fptr, &nkeys, nullptr, &status);
  if (status) {
    fits_report_error(stderr, status);
    fits_close_file(fptr, &status);
    std::stringstream errstr;
    errstr << "Error getting header space for: " << filename;
    throw affineExcept("beam", "readFile", errstr.str());
  }

  //Get the pixel scale
  fits_read_key(fptr, TDOUBLE, const_cast<char*>("PIXSIZE"),
		&pixsize,card, &status);
  if (status) {
    //Try pixscale
    status=0;
    fits_read_key(fptr, TDOUBLE, const_cast<char*>("PIXSCALE"),
		  &pixsize, card, &status);
    if (status) {
      //Ok, give up
      fits_report_error(stderr, status);
      fits_close_file(fptr, &status);
      std::stringstream errstr;
      errstr << "Unable to find pixel scale information in "
	     << filename;
      throw affineExcept("beam", "readFile", errstr.str());
    }
  }

  //Read the actual data
  int naxis;
  long naxes[2];
  fits_get_img_dim(fptr, &naxis, &status);
  fits_get_img_size(fptr, 2, naxes, &status);
  if (status || naxis != 2) {
    fits_close_file(fptr, &status);
    throw affineExcept("beam", "readFile", "Input BEAM is not 2D");
  }

  unsigned int n = naxes[0] * naxes[1];
  double* pixarr = new double[n];
  long fpixel[2];
  fpixel[0] = 1;
  for (unsigned int i = 1; i <= static_cast<unsigned int>(naxes[1]); ++i) {
    fpixel[1] = i;
    if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], 0, 
		      pixarr + naxes[0]*(i-1), 0, &status)) 
      break;   /* jump out of loop on error */
  }

  //Clean up
  fits_close_file(fptr,&status);
  if (status) {
    fits_report_error(stderr,status);
    delete[] pixarr;
    throw affineExcept("beam", "readFile", "Error closing beam file");
  } 

  //Now, sort in to positive and negative
  npos = nneg = 0;
  for (unsigned int i = 0; i < n; ++i) {
    if (std::isnan(pixarr[i]) || std::isinf(pixarr[i])) continue;
    if (pixarr[i] > minval) npos++;
    if (pixarr[i] < negminval) nneg++;
  }

  if ((npos == 0) && (nneg == 0))
    throw affineExcept("beam", "readFile", "Beam file has only zero values");

  //Get pos and abs(neg) into arrays -- don't store values closer to 0
  // then minval
  unsigned int ctr;
  if (npos > 0) {
    ctr = 0;
    pospixarr = new double[npos];
    for (unsigned int i = 0; i < n; ++i) {
      if (std::isnan(pixarr[i]) || std::isinf(pixarr[i])) continue;
      if (pixarr[i] > minval) pospixarr[ctr++] = pixarr[i];
    }
    haspos = true;
  } else haspos = false;

  if (nneg > 0) {
    ctr = 0;
    negpixarr = new double[nneg];
    for (unsigned int i = 0; i < n; ++i) {
      if (std::isnan(pixarr[i]) || std::isinf(pixarr[i])) continue;
      if (pixarr[i] < negminval) negpixarr[ctr++] = fabs(pixarr[i]);
    }
    hasneg = true;
  } else hasneg = false;
  delete[] pixarr;

  //Store pixels in sorted order
  if (haspos) std::sort(pospixarr, pospixarr + npos);
  if (hasneg) std::sort(negpixarr, negpixarr + nneg);

  //Get totals
  double val;
  if (haspos) {
    totpos = pospixarr[0];
    totpossq = totpos * totpos;
    for (unsigned int i = 1; i < npos; ++i) {
      val = pospixarr[i];
      totpos += val;
      totpossq += val * val;
    }
  }

  if (hasneg) {
    totneg = negpixarr[0];
    totnegsq = totneg * totneg;
    for (unsigned int i = 1; i < nneg; ++i) {
      val = negpixarr[i];
      totneg += val;
      totnegsq += val * val;
    }
  }

  // Compute inverse pixel arrays
  if (haspos) {
    posinvpixarr = new double[npos];
    for (unsigned int i = 0; i < npos; ++i)
      posinvpixarr[i] = 1.0 / pospixarr[i];
  }
  if (hasneg) {
    neginvpixarr = new double[nneg];
    for (unsigned int i = 0; i < nneg; ++i)
      neginvpixarr[i] = 1.0 / negpixarr[i];
  }
}


/*!
  \param[in] NBINS Number of bins to use.

  Empty bins are discarded.  This actually stores the inverse pixel
  histogram
*/
void beam::makeHistogram(unsigned int NBINS) {
  if (NBINS == 0)
    throw affineExcept("beam", "makeHistogram",
		       "Invalid (non-positive) nbins");
  nbins = NBINS;

  // Clean up hist values
  is_pos_histogrammed = is_neg_histogrammed = false;
  posnbins = negnbins = 0;
  if (posweights != nullptr) { delete[] posweights; posweights = nullptr; }
  if (negweights != nullptr) { delete[] negweights; negweights = nullptr; }
  if (poshistval != nullptr) { delete[] poshistval; poshistval = nullptr; }
  if (neghistval != nullptr) { delete[] neghistval; neghistval = nullptr; }

  bool do_pos = haspos && (npos >= histothresh);
  bool do_neg = hasneg && (nneg >= histothresh);
  if (!(do_pos || do_neg)) return; // Simple

  // Temporary variables
  unsigned int* tmp_wt;
  double* tmp_val;
  tmp_wt = new unsigned int[nbins];
  tmp_val = new double[nbins];

  // Note that we histogram in the beam, -then- invert it.  The other
  // way doesn't seem to work as well
  const double log2outscale = 0.0014419741739063218; // log2(1.001)
  unsigned int idx;
  double curr_val;
  if (do_pos) {
    // Figure out bin size.  We work in log2 space
    double minlog2val = log2(pospixarr[0]) - log2outscale;
    double maxlog2val = log2(pospixarr[npos - 1]) + log2outscale;
    double histstep = (maxlog2val - minlog2val) / static_cast<double>(nbins);
    std::memset(tmp_wt, 0, nbins * sizeof(unsigned int));
    std::memset(tmp_val, 0, nbins * sizeof(double));
    for (unsigned int i = 0; i < npos; ++i) {
      curr_val = pospixarr[i];
      idx = static_cast<unsigned int>((log2(curr_val) - minlog2val) / histstep);
      tmp_wt[idx] += 1;
      tmp_val[idx] += curr_val;
    }
    for (unsigned int i = 0; i < nbins; ++i)
      if (tmp_wt[i] > 0) ++posnbins;
    idx = 0;
    posweights = new double[posnbins];
    poshistval = new double[posnbins];
    for (unsigned int i = 0; i < nbins; ++i)
      if (tmp_wt[i] > 0) {
	posweights[idx] = static_cast<double>(tmp_wt[i]);
	poshistval[idx] = tmp_val[i] / posweights[idx];
	++idx;
      }
    for (unsigned int i = 0; i < posnbins; ++i)
      poshistval[i] = 1.0 / poshistval[i];
    is_pos_histogrammed = true;

  }

  if (do_neg) {
    double minlog2val = log2(negpixarr[0]) - log2outscale;
    double maxlog2val = log2(negpixarr[nneg - 1]) + log2outscale;
    double histstep = (maxlog2val - minlog2val) / static_cast<double>(nbins);
    std::memset(tmp_wt, 0, nbins * sizeof(unsigned int));
    std::memset(tmp_val, 0, nbins * sizeof(double));
    for (unsigned int i = 0; i < nneg; ++i) {
      curr_val = negpixarr[i]; // Alread fabsed
      idx = static_cast<unsigned int>((log2(curr_val) - minlog2val) / histstep);
      tmp_wt[idx] += 1;
      tmp_val[idx] += curr_val;
    }
    for (unsigned int i = 0; i < nbins; ++i)
      if (tmp_wt[i] > 0) ++negnbins;
    idx = 0;
    negweights = new double[negnbins];
    neghistval = new double[negnbins];
    for (unsigned int i = 0; i < nbins; ++i)
      if (tmp_wt[i] > 0) {
	negweights[idx] = static_cast<double>(tmp_wt[i]);
	neghistval[idx] = tmp_val[i] / negweights[idx];
	++idx;
      }
    for (unsigned int i = 0; i < negnbins; ++i)
      neghistval[i] = 1.0 / neghistval[i];
    is_neg_histogrammed = true;
  }

  delete[] tmp_wt;
  delete[] tmp_val;
}

/*!
  \returns Minimum and maximum value of positive beam
*/
dblpair beam::getMinMaxPos() const {
  if (!haspos) 
    return std::make_pair(std::numeric_limits<double>::quiet_NaN(),
			  std::numeric_limits<double>::quiet_NaN());
  return std::make_pair(pospixarr[0], pospixarr[npos-1]);
}

/*!
  \returns Minimum and maximum absolute value of negative beam
*/
dblpair beam::getMinMaxNeg() const {
  if (!hasneg) 
    return std::make_pair(std::numeric_limits<double>::quiet_NaN(),
			  std::numeric_limits<double>::quiet_NaN());
  return std::make_pair(negpixarr[0], negpixarr[nneg-1]);
}

/*!
  \param[in] exponent Power to raise beam to
  \param[out] array Loaded with positive beam raised to exponent

  Size of array is unchecked and assumed to be at least large enough
  to hold all of pixels.
*/
void beam::powerPos(double exponent, double* array) const {
  if (!haspos) return;
  if (fabs(exponent) < 1e-4) {
    for (unsigned int i = 0; i < npos; ++i) array[i]=1.0;
    return;
  }
  for (unsigned int i = 0; i < npos; ++i)
    array[i] = pow(pospixarr[i], exponent);
}

/*!
  \param[in] exponent Power to raise beam to
  \param[out] array Loaded with beam raised to exponent

  Size of array is unchecked and assumed to be at least large enough
  to hold all of pixels.
*/
void beam::powerNeg(double exponent, double* array) const {
  if (! hasneg) return;
  if (fabs(exponent) < 1e-4) {
    for (unsigned int i = 0; i < nneg; ++i) array[i] = 1.0;
    return;
  }
  for (unsigned int i = 0; i < nneg; ++i)
    array[i] = pow(negpixarr[i], exponent);
}

/*!
  \returns Effective area of beam, in square degrees
*/
double beam::getEffectiveArea() const {
  if (!(haspos || hasneg)) return 0.0;
  double convfac = pixsize / 3600.0;
  double area = 0.0;
  if (haspos) area = totpos;
  if (hasneg) area += totneg;
  return area * convfac * convfac;
}

/*!
  \returns Effective area of positive beam, in square degrees
*/
double beam::getEffectiveAreaPos() const {
  if (!haspos) return 0.0;
  double convfac = pixsize / 3600.0;
  double area = totpos * convfac * convfac;
  return area;
}

/*!
  \returns Effective area of negative beam, in square degrees
*/
double beam::getEffectiveAreaNeg() const {
  if (!hasneg) return 0.0;
  double convfac = pixsize / 3600.0;
  double area = totneg * convfac * convfac;
  return area;
}

/*!
  \returns Effective area of beam, in pixels
*/
double beam::getEffectiveAreaPix() const {
  if (!(haspos || hasneg)) return 0.0;
  double area = 0.0;
  if (haspos) area = totpos;
  if (hasneg) area += totneg;
  return area;
}

/*!
  \returns Effective area of squared beam, in square degrees
*/
double beam::getEffectiveAreaSq() const {
  if (!(haspos || hasneg)) return 0.0;
  double convfac = pixsize / 3600.0;
  double area = 0.0;
  if (haspos) area = totpossq * convfac * convfac;
  if (hasneg) area += totnegsq * convfac*convfac;
  return area;
}

/*!
  \param[in] comm MPI communicator
  \param[in] dest Destination for messages
*/
void beam::sendSelf(MPI_Comm comm, int dest) const {
  MPI_Send(const_cast<double*>(&pixsize), 1, MPI_DOUBLE, dest,
	   pofd_mcmc::BEAMSENDPIXSIZE, comm);
  MPI_Send(const_cast<double*>(&minval), 1, MPI_DOUBLE, dest,
	   pofd_mcmc::BEAMSENDMINVAL, comm);
  MPI_Send(const_cast<unsigned int*>(&nbins), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::BEAMSENDNBINS, comm);

  // Pos beam
  MPI_Send(const_cast<unsigned int*>(&npos), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::BEAMSENDNPOS, comm);
  if (npos > 0) {
    MPI_Send(pospixarr, npos, MPI_DOUBLE, dest, 
	     pofd_mcmc::BEAMSENDPOSPIXARR, comm);
    MPI_Send(posinvpixarr, npos, MPI_DOUBLE, dest,
	     pofd_mcmc::BEAMSENDINVPOSPIXARR, comm);
    MPI_Send(const_cast<double*>(&totpos), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::BEAMSENDTOTPOS, comm);
    MPI_Send(const_cast<double*>(&totpossq), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::BEAMSENDTOTSQPOS, comm);
    MPI_Send(const_cast<bool*>(&is_pos_histogrammed), 1, MPI::BOOL, dest,
	     pofd_mcmc::BEAMSENDISPOSHIST, comm);
    if (is_pos_histogrammed) {
      MPI_Send(const_cast<unsigned int*>(&posnbins), 1, MPI_UNSIGNED, dest,
	       pofd_mcmc::BEAMSENDPOSNBINS, comm);
      MPI_Send(posweights, posnbins, MPI_DOUBLE, dest,
	       pofd_mcmc::BEAMSENDPOSWEIGHTS, comm);
      MPI_Send(poshistval, posnbins, MPI_DOUBLE, dest,
	       pofd_mcmc::BEAMSENDPOSHISTVAL, comm);
    }
  }

  // Neg beam
  MPI_Send(const_cast<unsigned int*>(&nneg), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::BEAMSENDNNEG, comm);
  if (nneg > 0) {
    MPI_Send(negpixarr, nneg, MPI_DOUBLE, dest, 
	     pofd_mcmc::BEAMSENDNEGPIXARR, comm);
    MPI_Send(neginvpixarr, nneg, MPI_DOUBLE, dest,
	     pofd_mcmc::BEAMSENDINVNEGPIXARR, comm);
    MPI_Send(const_cast<double*>(&totneg), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::BEAMSENDTOTNEG, comm);
    MPI_Send(const_cast<double*>(&totnegsq), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::BEAMSENDTOTSQNEG, comm);
    MPI_Send(const_cast<bool*>(&is_neg_histogrammed), 1, MPI::BOOL, dest,
	     pofd_mcmc::BEAMSENDISNEGHIST, comm);
    if (is_neg_histogrammed) {
      MPI_Send(const_cast<unsigned int*>(&negnbins), 1, MPI_UNSIGNED, dest,
	       pofd_mcmc::BEAMSENDNEGNBINS, comm);
      MPI_Send(negweights, negnbins, MPI_DOUBLE, dest,
	       pofd_mcmc::BEAMSENDNEGWEIGHTS, comm);
      MPI_Send(neghistval, negnbins, MPI_DOUBLE, dest,
	       pofd_mcmc::BEAMSENDNEGHISTVAL, comm);
    }
  }
}

/*!
  \param[in] comm MPI communicator
  \param[in] src Source for messages
*/
void beam::receiveCopy(MPI_Comm comm, int src) {
  MPI_Status Info;
  cleanup(); // Just easier to clear it all

  MPI_Recv(&pixsize, 1, MPI_DOUBLE, src, pofd_mcmc::BEAMSENDPIXSIZE,
	   comm, &Info);
  MPI_Recv(&minval, 1, MPI_DOUBLE, src, pofd_mcmc::BEAMSENDMINVAL, comm, &Info);
  MPI_Recv(&nbins, 1, MPI_UNSIGNED, src, pofd_mcmc::BEAMSENDNBINS, comm, &Info);

  //Pos side
  MPI_Recv(&npos, 1, MPI_UNSIGNED, src, pofd_mcmc::BEAMSENDNPOS,
	   comm, &Info);
  if (npos > 0) {
    pospixarr = new double[npos];
    posinvpixarr = new double[npos];
    MPI_Recv(pospixarr, npos, MPI_DOUBLE, src, pofd_mcmc::BEAMSENDPOSPIXARR,
	     comm, &Info);
    MPI_Recv(posinvpixarr, npos, MPI_DOUBLE, src,
	     pofd_mcmc::BEAMSENDINVPOSPIXARR, comm, &Info);
    MPI_Recv(&totpos, 1, MPI_DOUBLE, src, pofd_mcmc::BEAMSENDTOTPOS,
	     comm, &Info);
    MPI_Recv(&totpossq, 1, MPI_DOUBLE, src, pofd_mcmc::BEAMSENDTOTSQPOS,
	     comm, &Info);
    MPI_Recv(&is_pos_histogrammed, 1, MPI::BOOL, src,
	     pofd_mcmc::BEAMSENDISPOSHIST, comm, &Info);
    if (is_pos_histogrammed) {
      MPI_Recv(&posnbins, 1, MPI_UNSIGNED, src,
	       pofd_mcmc::BEAMSENDPOSNBINS, comm, &Info);
      if (posnbins > 0) {
	posweights = new double[posnbins];
	MPI_Recv(posweights, posnbins, MPI_DOUBLE, src,
		 pofd_mcmc::BEAMSENDPOSWEIGHTS, comm, &Info);
	poshistval = new double[posnbins];
	MPI_Recv(poshistval, posnbins, MPI_DOUBLE, src,
		 pofd_mcmc::BEAMSENDPOSHISTVAL, comm, &Info);
      }
    }
    haspos = true;
  } else haspos = false;

  //Neg side
  MPI_Recv(&nneg, 1, MPI_UNSIGNED, src, pofd_mcmc::BEAMSENDNNEG,
	   comm, &Info);
  if (nneg > 0) {
    negpixarr = new double[nneg];
    neginvpixarr = new double[nneg];
    MPI_Recv(negpixarr, nneg, MPI_DOUBLE, src, pofd_mcmc::BEAMSENDNEGPIXARR,
	     comm, &Info);
    MPI_Recv(neginvpixarr, nneg, MPI_DOUBLE, src,
	     pofd_mcmc::BEAMSENDINVNEGPIXARR, comm, &Info);
    MPI_Recv(&totneg, 1, MPI_DOUBLE, src, pofd_mcmc::BEAMSENDTOTNEG,
	     comm, &Info);
    MPI_Recv(&totnegsq, 1, MPI_DOUBLE, src, pofd_mcmc::BEAMSENDTOTSQNEG,
	     comm, &Info);
    MPI_Recv(&is_neg_histogrammed, 1, MPI::BOOL, src,
	     pofd_mcmc::BEAMSENDISNEGHIST, comm, &Info);
    if (is_neg_histogrammed) {
      MPI_Recv(&negnbins, 1, MPI_UNSIGNED, src,
	       pofd_mcmc::BEAMSENDNEGNBINS, comm, &Info);
      if (negnbins > 0) {
	negweights = new double[negnbins];
	MPI_Recv(negweights, negnbins, MPI_DOUBLE, src,
		 pofd_mcmc::BEAMSENDNEGWEIGHTS, comm, &Info);
	neghistval = new double[negnbins];
	MPI_Recv(neghistval, negnbins, MPI_DOUBLE, src,
		 pofd_mcmc::BEAMSENDNEGHISTVAL, comm, &Info);
      }
    }
    hasneg = true;
  } else hasneg = false;
}

