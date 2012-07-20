//beam.cc

#include <limits>
#include <iostream>
#include <cmath>

#include<fitsio.h>
#include<beam.h>
#include<global_settings.h>
#include<affineExcept.h>

const unsigned int beam::histothresh = 15;

beam::beam() {
  nneg = npos = 0; 
  totneg = totpos = 0; 
  totnegsq = totpossq = 0.0;
  pixsize=0.0; 
  pospixarr = negpixarr = NULL;
  posinvpixarr = neginvpixarr = NULL;
  haspos = hasneg = false;
  hasposweights = hasnegweights = false;
  posweights = negweights = NULL;
}

/*!
  \param[in] filename   File to read from
  \param[in] histogram  Do beam histogramming
  \param[in] histogramlogstep Log step for histogram bins
*/
beam::beam( const std::string& filename, bool histogram, 
	    double histogramlogstep ) {
  nneg = npos = 0; 
  totneg = totpos = 0; 
  totnegsq = totpossq = 0.0;
  pixsize=0.0; 
  pospixarr = negpixarr = NULL;
  posinvpixarr = neginvpixarr = NULL;
  haspos = hasneg = false;
  hasposweights = hasnegweights = false;
  posweights = negweights = NULL;
  readFile(filename, histogram, histogramlogstep);
}

beam::beam( const beam& inp ) {
  negpixarr = pospixarr = NULL;
  cleanup();
  npos = inp.npos;
  nneg = inp.nneg;
  haspos = inp.haspos;
  hasneg = inp.hasneg;
  hasposweights = inp.hasposweights;
  hasnegweights = inp.hasnegweights;
  pixsize = inp.pixsize;
  totpos = inp.totpos;
  totneg = inp.totneg;
  totpossq = inp.totpossq;
  totnegsq = inp.totnegsq;
  if (haspos) {
    pospixarr = new double[npos];
    for (unsigned int i = 0; i < npos; ++i) pospixarr[i]=inp.pospixarr[i];
    posinvpixarr = new double[npos];
    for (unsigned int i = 0; i < npos; ++i) 
      posinvpixarr[i]=inp.posinvpixarr[i];
  }
  if (hasposweights) {
    posweights = new double[npos];
    for (unsigned int i = 0; i < npos; ++i) posweights[i]=inp.posweights[i];
  }
  if (hasneg) {
    negpixarr = new double[nneg];
    for (unsigned int i = 0; i < nneg; ++i) negpixarr[i]=inp.negpixarr[i];
    neginvpixarr = new double[nneg];
    for (unsigned int i = 0; i < nneg; ++i) 
      neginvpixarr[i]=inp.neginvpixarr[i];
  }
  if (hasnegweights) {
    negweights = new double[nneg];
    for (unsigned int i = 0; i < nneg; ++i) negweights[i]=inp.negweights[i];
  }
}

beam& beam::operator=(const beam& other) {
  if (this == &other) return *this;
  cleanup();
  npos = other.npos;
  nneg = other.nneg;
  haspos = other.haspos;
  hasneg = other.hasneg;
  hasposweights = other.hasposweights;
  hasnegweights = other.hasnegweights;
  pixsize = other.pixsize;
  totpos = other.totpos;
  totneg = other.totneg;
  totpossq = other.totpossq;
  totnegsq = other.totnegsq;
  if (haspos) {
    pospixarr = new double[npos];
    for (unsigned int i = 0; i < npos; ++i) pospixarr[i]=other.pospixarr[i];
    posinvpixarr = new double[npos];
    for (unsigned int i = 0; i < npos; ++i) 
      posinvpixarr[i]=other.posinvpixarr[i];
  }
  if (hasposweights) {
    posweights = new double[npos];
    for (unsigned int i = 0; i < npos; ++i) posweights[i]=other.posweights[i];
  }
  if (hasneg) {
    negpixarr = new double[nneg];
    for (unsigned int i = 0; i < nneg; ++i) negpixarr[i]=other.negpixarr[i];
    neginvpixarr = new double[nneg];
    for (unsigned int i = 0; i < nneg; ++i) 
      neginvpixarr[i]=other.neginvpixarr[i];
  }
  if (hasnegweights) {
    negweights = new double[nneg];
    for (unsigned int i = 0; i < nneg; ++i) negweights[i]=other.negweights[i];
  }
  return *this;
}

void beam::cleanup() {
  if (pospixarr != NULL) delete[] pospixarr;
  pospixarr = NULL;
  if (negpixarr != NULL) delete[] negpixarr;
  negpixarr = NULL;
  if (posinvpixarr != NULL) delete[] posinvpixarr;
  posinvpixarr = NULL;
  if (neginvpixarr != NULL) delete[] neginvpixarr;
  neginvpixarr = NULL;
  if (posweights != NULL) delete[] posweights;
  posweights = NULL;
  if (negweights != NULL) delete[] negweights;
  negweights = NULL;
  haspos = hasneg = false;
  hasposweights = hasnegweights = false;
  npos = nneg = 0;
  pixsize = 0.0;
  totpos = totneg = 0.0;
  totpossq = totnegsq = 0.0;
}

bool beam::revSort(const double& d1, const double& d2) const {
  return d1 > d2;
}

/*!
  \param[in] filename   File to read from
  \param[in] histogram Do beam histogramming?
  \param[in] histogramlogstep Log step if histogramming
*/
bool beam::readFile(const std::string& filename,bool histogram, 
		    double histogramlogstep) {

  int status;
  fitsfile *fptr;
  char card[FLEN_CARD];
  int nkeys;

  cleanup(); //Free arrays, etc
  status = 0;
  fptr = NULL;

  fits_open_file(&fptr, filename.c_str(), READONLY, &status);
  if (fptr == NULL) {
    std::cerr << "Error opening input file: " << filename << std::endl;
    fits_close_file(fptr,&status);
    return false;
  }
  if (status) {
    std::cerr << "Error opening input file: " << filename << std::endl;
    fits_report_error(stderr,status);
    fits_close_file(fptr,&status);
    return false;
  }

  fits_get_hdrspace(fptr,&nkeys,NULL,&status);
  if (status) {
    std::cerr << "Error getting header space for: " << filename << std::endl;
    fits_report_error(stderr,status);
    fits_close_file(fptr,&status);
    return false;
  }

  //Get the pixel scale
  fits_read_key(fptr, TDOUBLE, const_cast<char*>("PIXSIZE"),
		&pixsize,card,&status);
  if (status) {
    //Try pixscale
    status=0;
    fits_read_key(fptr, TDOUBLE, const_cast<char*>("PIXSCALE"),
		  &pixsize,card,&status);
    if (status) {
      //Ok, give up
      std::cerr << "Unable to find pixel scale information in "
		<< filename << std::endl;
      fits_report_error(stderr,status);
      fits_close_file(fptr,&status);
      return false;
    }
  }

  //Read the actual data
  int naxis;
  long naxes[2];
  fits_get_img_dim(fptr, &naxis, &status);
  fits_get_img_size(fptr, 2, naxes, &status);
  if (status || naxis != 2) {
    std::cerr << "ERROR: input BEAM is not 2D" << std::endl;
    fits_close_file(fptr,&status);
    return false;
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
    std::cerr << "Error closing BEAM file" << std::endl;
    fits_report_error(stderr,status);
    delete[] pixarr;
    return false;
  } 

  //Now, sort in to positive and negative
  npos = nneg = 0;
  for (unsigned int i = 0; i < n; ++i) {
    if ( std::isnan(pixarr[i]) || std::isinf(pixarr[i] ) ) continue;
    if (pixarr[i] > 0.0) npos++;
    if (pixarr[i] < 0.0) nneg++;
  }

  if ( (npos == 0) && (hasneg == 0) )
    throw affineExcept("beam","readFile",
		       "Beam file has only zero values",1);

  //Get pos and abs(neg) into arrays -- don't store 0 values
  unsigned int ctr;
  if (npos > 0) {
    ctr = 0;
    pospixarr = new double[npos];
    for (unsigned int i = 0; i < n; ++i) {
      if ( std::isnan(pixarr[i]) || std::isinf(pixarr[i] ) ) continue;
      if (pixarr[i] > 0.0) pospixarr[ctr++] = pixarr[i];
    }
    haspos = true;
  }

  if (nneg > 0) {
    ctr = 0;
    negpixarr = new double[nneg];
    for (unsigned int i = 0; i < n; ++i) {
      if ( std::isnan(pixarr[i]) || std::isinf(pixarr[i] ) ) continue;
      if (pixarr[i] < 0.0) negpixarr[ctr++] = fabs(pixarr[i]);
    }
    hasneg = true;
  }
  delete[] pixarr;


  //Store pixels in sorted order, get inv pixarr
  if (haspos) std::sort(pospixarr,pospixarr+npos);
  if (hasneg) std::sort(negpixarr,negpixarr+nneg);

  //Get totals
  double val;
  if (haspos) {
    totpos = pospixarr[0];
    totpossq = totpos*totpos;
    for (unsigned int i = 1; i < npos; ++i) {
      val = pospixarr[i];
      totpos += val;
      totpossq += val*val;
    }
  }

  if (hasneg) {
    totneg = negpixarr[0];
    totnegsq = totneg*totneg;
    for (unsigned int i = 1; i < nneg; ++i) {
      val = negpixarr[i];
      totneg = val;
      totnegsq = val*val;
    }
  }

  //Histogram if desired
  if (histogram) {
    //The strategy is to bin in log space in each bin,
    // averaging over the beam values in each bin with a weight.
    // For some reason, averaging over the 1/beam values works
    // quite badly, even though that in principle is something
    // we care about more.
    if (histogramlogstep <= 0)
      throw affineExcept("beam","readFile",
			 "Invalid (non-positive) histogramlogstep",2);
    if (haspos && (npos > histothresh)) {
      //Find min values, now that zeros have
      // been removed and absolute values taken
      double minval, maxval, val;
      unsigned int idx;
      minval = log(pospixarr[0]);
      maxval = log(pospixarr[npos-1]);
      unsigned int nbins = 
	static_cast<unsigned int>( (maxval-minval)/histogramlogstep )+2;
      unsigned int *ninbin = new unsigned int[nbins];
      double *meaninbin = new double[nbins];
      for (unsigned int i = 0; i < nbins; ++i) ninbin[i]=0;
      for (unsigned int i = 0; i < nbins; ++i) meaninbin[i]=0.0;
      for (unsigned int i = 0; i < npos; ++i) {
	val = pospixarr[i];
	idx = static_cast<unsigned int>( (log(val)-minval)/histogramlogstep );
	ninbin[idx] += 1;
	meaninbin[idx] += val;
      }
      //New arrays
      unsigned int newnpix = 0;
      for (unsigned int i = 0; i < nbins; ++i)
	if (ninbin[i] > 0) ++newnpix;
      double* newpixarr  = new double[newnpix];
      posweights = new double[newnpix];
      unsigned int ctr = 0;
      for (unsigned int i = 0; i < nbins; ++i) 
	if (ninbin[i] > 0) {
	  val = static_cast<double>(ninbin[i]);
	  posweights[ctr] = val;
	  newpixarr[ctr] = meaninbin[i]/val;
	  ++ctr;
	}
      delete[] pospixarr;
      pospixarr = newpixarr;
      npos = newnpix;
      hasposweights = true;
      delete[] ninbin;
      delete[] meaninbin;
    }
    if (hasneg && (nneg > histothresh)) {
      double minval, maxval, val;
      unsigned int idx;
      minval = log(negpixarr[0]);
      maxval = log(negpixarr[nneg-1]);
      unsigned int nbins = 
	static_cast<unsigned int>( (maxval-minval)/histogramlogstep )+2;
      unsigned int *ninbin = new unsigned int[nbins];
      double *meaninbin = new double[nbins];
      for (unsigned int i = 0; i < nbins; ++i) ninbin[i]=0;
      for (unsigned int i = 0; i < nbins; ++i) meaninbin[i]=0.0;
      for (unsigned int i = 0; i < nneg; ++i) {
	val = negpixarr[i];
	idx = static_cast<unsigned int>( (log(val)-minval)/histogramlogstep );
	ninbin[idx] += 1;
	meaninbin[idx] += val;
      }
      unsigned int newnpix = 0;
      for (unsigned int i = 0; i < nbins; ++i)
	if (ninbin[i] > 0) ++newnpix;
      double* newpixarr = new double[newnpix];
      negweights   = new double[newnpix];
      unsigned int ctr = 0;
      for (unsigned int i = 0; i < nbins; ++i) 
	if (ninbin[i] > 0) {
	  val = static_cast<double>(ninbin[i]);
	  negweights[ctr] = val;
	  newpixarr[ctr] = meaninbin[i]/val;
	  ++ctr;
	}
      delete[] negpixarr;
      negpixarr = newpixarr;
      nneg = newnpix;
      hasnegweights = true;
      delete[] ninbin;
      delete[] meaninbin;
    }
  }

  if (haspos) {
    posinvpixarr = new double[npos];
    for (unsigned int i = 0; i < npos; ++i)
      posinvpixarr[i] = 1.0/pospixarr[i];
  }
  if (hasneg) {
    neginvpixarr = new double[nneg];
    for (unsigned int i = 0; i < nneg; ++i)
      neginvpixarr[i] = 1.0/negpixarr[i];
  }
  return true;
}

double beam::getMinPos() const {
  if (! haspos ) return std::numeric_limits<double>::quiet_NaN();
  return pospixarr[0];
}

double beam::getMaxPos() const {
  if (! haspos ) return std::numeric_limits<double>::quiet_NaN();
  return pospixarr[npos-1];
}

double beam::getMinAbsNeg() const {
  if (! hasneg ) return std::numeric_limits<double>::quiet_NaN();
  return negpixarr[0];
}

double beam::getMaxAbsNeg() const {
  if (! hasneg ) return std::numeric_limits<double>::quiet_NaN();
  return negpixarr[nneg-1];
}

std::pair<double,double> beam::getRangePos() const {
  if (!haspos) {
    double val = std::numeric_limits<double>::quiet_NaN();
    return std::pair< double, double > ( val, val );
  }
  double lowval = *std::min_element( pospixarr, pospixarr+npos );
  double highval = *std::max_element( pospixarr, pospixarr+npos );
  return std::pair< double, double > ( lowval, highval );
}

std::pair<double,double> beam::getRangeNeg() const {
  if (!hasneg) {
    double val = std::numeric_limits<double>::quiet_NaN();
    return std::pair< double, double > ( val, val );
  }
  double lowval = *std::min_element( negpixarr, negpixarr+nneg );
  double highval = *std::max_element( negpixarr, negpixarr+nneg );
  return std::pair< double, double > ( lowval, highval );
}

/*!
  Size of array is unchecked and assumed to be at least large enough
  to hold all of pixels.
 */
void beam::powerPos( double exponent, double* array ) const {
  if (!haspos) return;
  if ( fabs(exponent) < 1e-4 ) {
    for (unsigned int i = 0; i < npos; ++i) array[i]=1.0;
    return;
  }
  for (unsigned int i = 0; i < npos; ++i)
    array[i] = pow( pospixarr[i], exponent );
}

/*!
  Size of array is unchecked and assumed to be at least large enough
  to hold all of pixels.
 */
void beam::powerNeg( double exponent, double* array ) const {
  if ( ! hasneg ) return;
  if ( fabs(exponent) < 1e-4 ) {
    for (unsigned int i = 0; i < nneg; ++i) array[i]=1.0;
    return;
  }
  for (unsigned int i = 0; i < nneg; ++i)
    array[i] = pow( negpixarr[i], exponent );
}

double beam::getEffectiveArea() const {
  if (! (haspos || hasneg) )
    return std::numeric_limits<double>::quiet_NaN();
  double convfac = pixsize / 3600.0;
  double area = 0.0;
  if (haspos) area = totpos * convfac * convfac;
  if (hasneg) area += totneg * convfac*convfac;
  return area;
}

double beam::getEffectiveAreaPos() const {
  if (!haspos) return std::numeric_limits<double>::quiet_NaN();
  double convfac = pixsize / 3600.0;
  double area = totpos * convfac * convfac;
  return area;
}

double beam::getEffectiveAreaNeg() const {
  if (!hasneg) return std::numeric_limits<double>::quiet_NaN();
  double convfac = pixsize / 3600.0;
  double area = totneg * convfac * convfac;
  return area;
}

double beam::getEffectiveAreaPix() const {
  if (! (haspos || hasneg) )
    return std::numeric_limits<double>::quiet_NaN();
  double area = 0.0;
  if (haspos) area = totpos;
  if (hasneg) area += totneg;
  return area;
}

double beam::getEffectiveAreaSq() const {
  if (! (haspos || hasneg) )
    return std::numeric_limits<double>::quiet_NaN();
  double convfac = pixsize / 3600.0;
  double area = 0.0;
  if (haspos) area = totpossq * convfac * convfac;
  if (hasneg) area += totnegsq * convfac*convfac;
  return area;
}

void beam::SendSelf(MPI::Comm& comm, int dest) const {
  comm.Send(&pixsize,1,MPI::DOUBLE,dest,pofd_mcmc::BEAMSENDPIXSIZE);
  comm.Send(&npos,1,MPI::UNSIGNED,dest,pofd_mcmc::BEAMSENDNPOS);

  if (npos > 0) {
    comm.Send(&hasposweights,1,MPI::BOOL,dest,
	      pofd_mcmc::BEAMSENDHASPOSWEIGHTS);
    comm.Send(pospixarr,npos,MPI::DOUBLE,dest,
	      pofd_mcmc::BEAMSENDPOSPIXARR);
    comm.Send(posinvpixarr,npos,MPI::DOUBLE,dest,
	      pofd_mcmc::BEAMSENDINVPOSPIXARR);
    if (hasposweights)
      comm.Send(posweights,npos,MPI::DOUBLE,dest,
		pofd_mcmc::BEAMSENDPOSWEIGHTS);
    comm.Send(&totpos,1,MPI::DOUBLE,dest,pofd_mcmc::BEAMSENDTOTPOS);
    comm.Send(&totpossq,1,MPI::DOUBLE,dest,pofd_mcmc::BEAMSENDTOTSQPOS);
  }
  comm.Send(&nneg,1,MPI::UNSIGNED,dest,pofd_mcmc::BEAMSENDNNEG);
  if (nneg > 0) {
    comm.Send(&hasnegweights,1,MPI::BOOL,dest,
	      pofd_mcmc::BEAMSENDHASNEGWEIGHTS);
    comm.Send(negpixarr,nneg,MPI::DOUBLE,dest,
	      pofd_mcmc::BEAMSENDNEGPIXARR);
    comm.Send(neginvpixarr,nneg,MPI::DOUBLE,dest,
	      pofd_mcmc::BEAMSENDINVNEGPIXARR);
    if (hasnegweights)
      comm.Send(negweights,npos,MPI::DOUBLE,dest,
		pofd_mcmc::BEAMSENDNEGWEIGHTS);
    comm.Send(&totneg,1,MPI::DOUBLE,dest,pofd_mcmc::BEAMSENDTOTNEG);
    comm.Send(&totnegsq,1,MPI::DOUBLE,dest,pofd_mcmc::BEAMSENDTOTSQNEG);
  }
}

void beam::RecieveCopy(MPI::Comm& comm, int src) {
  unsigned int new_n;
  comm.Recv(&pixsize,1,MPI::DOUBLE,src,pofd_mcmc::BEAMSENDPIXSIZE);

  //Pos side
  comm.Recv(&new_n,1,MPI::UNSIGNED,src,pofd_mcmc::BEAMSENDNPOS);
  if (new_n != npos) {
    if (pospixarr != NULL) delete[] pospixarr;
    if (posinvpixarr != NULL) delete[] posinvpixarr;
    if (posweights != NULL) delete[] posweights;
    posweights = NULL;
    if (new_n > 0) pospixarr = new double[new_n]; else pospixarr = NULL;
    if (new_n > 0) posinvpixarr = new double[new_n]; else posinvpixarr = NULL;
    npos = new_n;
  }
  if (npos > 0) {
    comm.Recv(&hasposweights,1,MPI::BOOL,src,
	      pofd_mcmc::BEAMSENDHASPOSWEIGHTS);
    comm.Recv(pospixarr,npos,MPI::DOUBLE,src,pofd_mcmc::BEAMSENDPOSPIXARR);
    comm.Recv(posinvpixarr,npos,MPI::DOUBLE,src,
	      pofd_mcmc::BEAMSENDINVPOSPIXARR);
    if (hasposweights) {
      if (posweights == NULL) posweights = new double[npos];
      comm.Recv(posweights,npos,MPI::DOUBLE,src,
		pofd_mcmc::BEAMSENDPOSWEIGHTS);
    }
    comm.Recv(&totpos,1,MPI::DOUBLE,src,pofd_mcmc::BEAMSENDTOTPOS);
    comm.Recv(&totpossq,1,MPI::DOUBLE,src,pofd_mcmc::BEAMSENDTOTSQPOS);
    haspos = true;
  } else haspos = false;

  //Neg side
  comm.Recv(&new_n,1,MPI::UNSIGNED,src,pofd_mcmc::BEAMSENDNNEG);
  if (new_n != nneg) {
    if (negpixarr != NULL) delete[] negpixarr;
    if (neginvpixarr != NULL) delete[] neginvpixarr;
    if (negweights != NULL) delete[] negweights;
    negweights = NULL;
    if (new_n > 0) negpixarr = new double[new_n]; else negpixarr = NULL;
    if (new_n > 0) neginvpixarr = new double[new_n]; else neginvpixarr = NULL;
    nneg = new_n;
  }
  if (nneg > 0) {
    comm.Recv(&hasnegweights,1,MPI::BOOL,src,
	      pofd_mcmc::BEAMSENDHASNEGWEIGHTS);
    comm.Recv(negpixarr,nneg,MPI::DOUBLE,src,pofd_mcmc::BEAMSENDNEGPIXARR);
    comm.Recv(neginvpixarr,nneg,MPI::DOUBLE,src,
	      pofd_mcmc::BEAMSENDINVNEGPIXARR);
    if (hasnegweights) {
      if (negweights == NULL) negweights = new double[nneg];
      comm.Recv(negweights,nneg,MPI::DOUBLE,src,
		pofd_mcmc::BEAMSENDNEGWEIGHTS);
    }
    comm.Recv(&totneg,1,MPI::DOUBLE,src,pofd_mcmc::BEAMSENDTOTNEG);
    comm.Recv(&totnegsq,1,MPI::DOUBLE,src,pofd_mcmc::BEAMSENDTOTSQNEG);
    hasneg = true;
  } else hasneg = false;

}
