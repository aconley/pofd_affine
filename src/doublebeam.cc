//doublebeam.cc

#include <limits>
#include <iostream>
#include <cmath>
#include <sstream>
#include <map>

#include<fitsio.h>

#include "../include/doublebeam.h"
#include "../include/global_settings.h"
#include "../include/affineExcept.h"

const unsigned int doublebeam::histothresh = 15;

doublebeam::doublebeam() {
  hassign[0] = hassign[1] = hassign[2] = hassign[3] = false;
  npix[0] = npix[1] = npix[2] = npix[3] = 0;
  tot1[0] =tot1[1] = tot1[2] = tot1[3] = 0.0;
  tot2[0] =tot2[1] = tot2[2] = tot2[3] = 0.0;
  totsq1[0] = totsq1[1] = totsq1[2] = totsq1[3] = 0.0;
  totsq2[0] = totsq2[1] = totsq2[2] = totsq2[3] = 0.0;
  totsm1 = totsm2 = 0.0;
  pixsize = 0.0; 
  pixarr1[0] = pixarr1[1] = pixarr1[2] = pixarr1[3] = NULL;
  pixarr2[0] = pixarr2[1] = pixarr2[2] = pixarr2[3] = NULL;
  ipixarr1[0] = ipixarr1[1] = ipixarr1[2] = ipixarr1[3] = NULL;
  ipixarr2[0] = ipixarr2[1] = ipixarr2[2] = ipixarr2[3] = NULL;
  has_weights[0] = has_weights[1] = has_weights[2] = has_weights[3] = false;
  weights[0] = weights[1] = weights[2] = weights[3] = NULL;
}

/*!
  \param[in] filename1   File to read from, beam 1
  \param[in] filename2   File to read from, beam 2
  \param[in] histogram Do beam histogramming?
  \param[in] histogramlogstep Log step if histogramming
*/
doublebeam::doublebeam(const std::string& filename1,
		       const std::string& filename2,
		       bool histogram, double histogramlogstep ) {
  hassign[0] = hassign[1] = hassign[2] = hassign[3] = false;
  npix[0] = npix[1] = npix[2] = npix[3] = 0;
  tot1[0] = tot1[1] = tot1[2] = tot1[3] = 0.0;
  tot2[0] = tot2[1] = tot2[2] = tot2[3] = 0.0;
  totsq1[0] = totsq1[1] = totsq1[2] = totsq1[3] = 0.0;
  totsq2[0] = totsq2[1] = totsq2[2] = totsq2[3] = 0.0;
  totsm1 = totsm2 = 0.0;
  pixsize = 0.0; 
  pixarr1[0] = pixarr1[1] = pixarr1[2] = pixarr1[3] = NULL;
  pixarr2[0] = pixarr2[1] = pixarr2[2] = pixarr2[3] = NULL;
  ipixarr1[0] = ipixarr1[1] = ipixarr1[2] = ipixarr1[3] = NULL;
  ipixarr2[0] = ipixarr2[1] = ipixarr2[2] = ipixarr2[3] = NULL;
  has_weights[0] = has_weights[1] = has_weights[2] = has_weights[3] = false;
  weights[0] = weights[1] = weights[2] = weights[3] = NULL;
  readFiles(filename1, filename2, histogram, histogramlogstep);
}

doublebeam::doublebeam(const doublebeam& inp) {
  pixarr1[0] = pixarr1[1] = pixarr1[2] = pixarr1[3] = NULL;
  pixarr2[0] = pixarr2[1] = pixarr2[2] = pixarr2[3] = NULL;
  ipixarr1[0] = ipixarr1[1] = ipixarr1[2] = ipixarr1[3] = NULL;
  ipixarr2[0] = ipixarr2[1] = ipixarr2[2] = ipixarr2[3] = NULL;
  weights[0] = weights[1] = weights[2] = weights[3] = NULL;
  for (unsigned int i = 0; i < 4; ++i) npix[i] = inp.npix[i];
  for (unsigned int i = 0; i < 4; ++i) hassign[i] = inp.hassign[i];
  for (unsigned int i = 0; i < 4; ++i) tot1[i] = inp.tot1[i];
  for (unsigned int i = 0; i < 4; ++i) tot2[i] = inp.tot2[i];
  for (unsigned int i = 0; i < 4; ++i) totsq1[i] = inp.totsq1[i];
  for (unsigned int i = 0; i < 4; ++i) totsq2[i] = inp.totsq2[i];
  totsm1 = inp.totsm1;
  totsm2 = inp.totsm2;
  pixsize = inp.pixsize;
  for (unsigned int  i = 0; i < 4; ++i)
    if (hassign[i]) {
      pixarr1[i] = new double[npix[i]];
      pixarr2[i] = new double[npix[i]];
      ipixarr1[i] = new double[npix[i]];
      ipixarr2[i] = new double[npix[i]];
      for (unsigned int j = 0; j < npix[i]; ++j)
	pixarr1[i][j] = inp.pixarr1[i][j];
      for (unsigned int j = 0; j < npix[i]; ++j)
	pixarr2[i][j] = inp.pixarr2[i][j];
      for (unsigned int j = 0; j < npix[i]; ++j)
	ipixarr1[i][j] = inp.ipixarr1[i][j];
      for (unsigned int j = 0; j < npix[i]; ++j)
	ipixarr2[i][j] = inp.ipixarr2[i][j];
      has_weights[i] = inp.has_weights[i];
      if (has_weights[i]) {
	weights[i] = new double[npix[i]];
	for (unsigned int j = 0; j < npix[i]; ++j)
	  weights[i][j] = inp.weights[i][j];
      }
    }
}

/*!
  \param[in] other doublebeam to copy
*/
doublebeam& doublebeam::operator=(const doublebeam& other) {
  if (this == &other) return *this;
  cleanup();
  for (unsigned int i = 0; i < 4; ++i) npix[i] = other.npix[i];
  for (unsigned int i = 0; i < 4; ++i) hassign[i] = other.hassign[i];
  for (unsigned int i = 0; i < 4; ++i) tot1[i] = other.tot1[i];
  for (unsigned int i = 0; i < 4; ++i) tot2[i] = other.tot2[i];
  for (unsigned int i = 0; i < 4; ++i) totsq1[i] = other.totsq1[i];
  for (unsigned int i = 0; i < 4; ++i) totsq2[i] = other.totsq2[i];
  totsm1 = other.totsm1;
  totsm2 = other.totsm2;
  pixsize = other.pixsize;
  for (unsigned int i = 0; i < 4; ++i)
    if (hassign[i]) {
      pixarr1[i] = new double[npix[i]];
      pixarr2[i] = new double[npix[i]];
      ipixarr1[i] = new double[npix[i]];
      ipixarr2[i] = new double[npix[i]];
      for (unsigned int j = 0; j < npix[i]; ++j)
	pixarr1[i][j] = other.pixarr1[i][j];
      for (unsigned int j = 0; j < npix[i]; ++j)
	pixarr2[i][j] = other.pixarr2[i][j];
      for (unsigned int j = 0; j < npix[i]; ++j)
	ipixarr1[i][j] = other.ipixarr1[i][j];
      for (unsigned int j = 0; j < npix[i]; ++j)
	ipixarr2[i][j] = other.ipixarr2[i][j];
      has_weights[i] = other.has_weights[i];
      if (has_weights[i]) {
	weights[i] = new double[npix[i]];
	for (unsigned int j = 0; j < npix[i]; ++j)
	  weights[i][j] = other.weights[i][j];
      }
    }
  return *this;
}

unsigned int doublebeam::getTotalNPix() const {
  unsigned int n = 0;
  for (unsigned int i = 0; i < 4; ++i)
    if (hassign[i]) n += npix[i];
  return n;
}

void doublebeam::cleanup() {
  for (unsigned int i = 0; i < 4; ++i) {
    if (pixarr1[i] != NULL) delete[] pixarr1[i];
    if (pixarr2[i] != NULL) delete[] pixarr2[i];
    if (ipixarr1[i] != NULL) delete[] ipixarr1[i];
    if (ipixarr2[i] != NULL) delete[] ipixarr2[i];
    if (weights[i] != NULL) delete[] weights[i];
    pixarr1[i] = pixarr2[i] = NULL;
    ipixarr1[i] = ipixarr2[i] = NULL;
    weights[i] = NULL;
  }
  for (unsigned int i = 0; i < 4; ++i) npix[i] = 0;
  for (unsigned int i = 0; i < 4; ++i) hassign[i] = false;
  for (unsigned int i = 0; i < 4; ++i) tot1[i] = 0.0;
  for (unsigned int i = 0; i < 4; ++i) tot2[i] = 0.0;
  for (unsigned int i = 0; i < 4; ++i) totsq1[i] = 0.0;
  for (unsigned int i = 0; i < 4; ++i) totsq2[i] = 0.0;
  totsm1 = totsm2 = 0.0;
  for (unsigned int i = 0; i < 4; ++i) has_weights[i] = false;
  pixsize = 0.0;
}

/*!
  \param[in] n Size of beam arrays
  \param[in] bm1 Beam 1 array
  \param[in] bm2 Beam 2 array
  \param[in] pixscale Pixel scale in arcsec per pixel
  \param[in] histogram Do beam histogramming?
  \param[in] histogramlogstep Log step if histogramming

  Will ignore pixels that are zero or with an absolute value
  larger than one.  Also may histogram beams, depending on settings
*/
void doublebeam::setBeams(unsigned int n, const double* const bm1,
			  const double* const bm2, double pixscale,
			  bool histogram, double histogramlogstep ) {
  cleanup();  //Free arrays, etc.
  if (n == 0) return;

  //Now, we sort into pp, pn, np, nn bits
  //First step -- count how many
  double bmval1, bmval2;
  for (unsigned int i = 0; i < n; ++i) {
    bmval1 = bm1[i]; 
    bmval2 = bm2[i]; 
    if (std::isnan(bmval1) || std::isinf(bmval1)) continue;
    if (std::isnan(bmval2) || std::isinf(bmval2 )) continue;
    if ((bmval1 == 0) || (bmval2 == 0) ) continue;
    if ((fabs(bmval1) > 1) || (fabs(bmval2) > 1)) continue;
    if (bmval1 > 0.0)
      if (bmval2 > 0.0) ++npix[0]; else ++npix[1];
    else
      if (bmval2 > 0.0) ++npix[2]; else ++npix[3];
  }

  if ((npix[0]+npix[1]+npix[2]+npix[3]) == 0)
    throw affineExcept("doublebeam","setBeams",
		       "Doublebeam file has only zero values",1);

  //Get pos and abs(neg) into arrays -- don't store 0 values
  double abmval1, abmval2;
  unsigned int ppctr, pnctr, npctr, nnctr;
  for (unsigned int i = 0; i < 4; ++i) 
    if (npix[i] > 0) {
      pixarr1[i] = new double[npix[i]];
      pixarr2[i] = new double[npix[i]];
      ipixarr1[i] = new double[npix[i]];
      ipixarr2[i] = new double[npix[i]];
    }
  ppctr = pnctr = npctr = nnctr = 0;
  for (unsigned int i = 0; i < n; ++i) {
    bmval1 = bm1[i]; 
    bmval2 = bm2[i]; 
    if (std::isnan(bmval1) || std::isinf(bmval1)) continue;
    if (std::isnan(bmval2) || std::isinf(bmval2 )) continue;
    if ((bmval1 == 0) || (bmval2 == 0)) continue;
    abmval1 = fabs(bmval1);
    abmval2 = fabs(bmval2);
    if  ((abmval1 > 1) || (abmval2 > 1)) continue;
    if (bmval1 > 0.0)
      if (bmval2 > 0.0) {
	pixarr1[0][ppctr]=bmval1; 
	pixarr2[0][ppctr++]=bmval2;
      } else {
	pixarr1[1][pnctr]=bmval1; 
	pixarr2[1][pnctr++]=abmval2;
      }
    else
      if (bmval2 > 0.0) {
	pixarr1[2][npctr]=abmval1; 
	pixarr2[2][npctr++]=bmval2;
      } else {
	pixarr1[3][nnctr]=abmval1;
	pixarr2[3][nnctr++]=abmval2;
      }
  }
  for (unsigned int i = 0; i < 4; ++i)
    if (npix[i] > 0) {
      for (unsigned int j = 0; j < npix[i]; ++j)
	ipixarr1[i][j] = 1.0/pixarr1[i][j];
      for (unsigned int j = 0; j < npix[i]; ++j)
	ipixarr2[i][j] = 1.0/pixarr2[i][j];
      hassign[i] = true;
    }
  
  //Get totals
  double val;
  double *pptr;
  totsm1 = totsm2 = 0.0;
  for (unsigned int i = 0; i < 4; ++i) 
    if (hassign[i]) {
      pptr = pixarr1[i];
      val = pptr[0];
      tot1[i] = val;
      totsq1[i] = val * val;
      for (unsigned int j = 1; j < npix[i]; ++j) {
	val = pptr[j];
	tot1[i] += val;
	totsq1[i] += val * val;
      }
      totsm1 += tot1[i];
      pptr = pixarr2[i];
      val = pptr[0];
      tot2[i] = val;
      totsq2[i] = val * val;
      for (unsigned int j = 1; j < npix[i]; ++j) {
	val = pptr[j];
	tot2[i] += val;
	totsq2[i] += val * val;
      }
      totsm2 += tot2[i];
    }

  //Histogram if desired
  if (histogram) {
    //The strategy is to bin in log space in each bin,
    // averaging over the beam values in each bin with a weight.
    // For some reason, averaging over the 1/beam values works
    // quite badly, even though that in principle is something
    // we care about more.
    if (histogramlogstep <= 0)
      throw affineExcept("doublebeam","setBeams",
			 "Invalid (non-positive) histogramlogstep",2);
    for (unsigned int sgn = 0; sgn < 4; ++sgn)
      if (npix[sgn] > histothresh) {
	//Find min values for both bands, now that zeros have
	// been removed and absolute values taken
	double minval1, minval2, maxval2, val;
	minval1 = pixarr1[sgn][0];
	for (unsigned int i = 1; i < npix[sgn]; ++i) {
	  val = pixarr1[sgn][i];
	  if (val < minval1) minval1 = val;
	}
	minval2 = maxval2 = pixarr2[sgn][0];
	for (unsigned int i = 1; i < npix[sgn]; ++i) {
	  val = pixarr2[sgn][i];
	  if (val < minval2) minval2 = val;
	  if (val > maxval2) maxval2 = val;
	}
	minval1 = log(minval1); //Binning in log space
	minval2 = log(minval2);
	maxval2 = log(maxval2);

	unsigned int nbins2 = 
	  static_cast<unsigned int>((maxval2 - minval2) / histogramlogstep) + 2;

	//Now build up number of bins and mean values
	//This will be sparse, so use a map
	std::map<unsigned int, bmhist> hist;
	std::map<unsigned int, bmhist>::iterator histpos;
	bmhist histelem;
	unsigned int idx1, idx2, fullidx;
	double val1, val2;
	for (unsigned int i = 0; i < npix[sgn]; ++i) {
	  //Get index value, histogramming on beam value
	  //Since we are in log space, histogramming on beam or inverse
	  // beam is irrelevant
	  val1 = log(pixarr1[sgn][i]);
	  val2 = log(pixarr2[sgn][i]);
	  idx1 = static_cast<unsigned int>((val1 - minval1) / histogramlogstep);
	  idx2 = static_cast<unsigned int>((val2 - minval2) / histogramlogstep);

	  fullidx = idx1 * nbins2 + idx2;
	  
	  //See if already present
	  histpos = hist.find(fullidx);
	  if (histpos == hist.end()) {
	    //New
	    histelem.cnt = 1;
	    histelem.tot1 = pixarr1[sgn][i];
	    histelem.tot2 = pixarr2[sgn][i];
	    hist[fullidx] = histelem;
	  } else {
	    //Not new
	    histpos->second.cnt += 1;
	    histpos->second.tot1 += pixarr1[sgn][i];
	    histpos->second.tot2 += pixarr2[sgn][i];
	  }
	}

	//Copy new stuff over
	unsigned int newnpix = hist.size();
	delete[] pixarr1[sgn];
	pixarr1[sgn] = new double[newnpix];
	delete[] pixarr2[sgn];
	pixarr2[sgn] = new double[newnpix];
	delete[] ipixarr1[sgn];
	ipixarr1[sgn] = new double[newnpix];
	delete[] ipixarr2[sgn];
	ipixarr2[sgn] = new double[newnpix];
	weights[sgn] = new double[newnpix];
	unsigned int ctr=0, npres;
	for (histpos = hist.begin(); histpos != hist.end(); ++histpos) {
	  npres = histpos->second.cnt;
	  weights[sgn][ctr] = npres;
	  pixarr1[sgn][ctr] = histpos->second.tot1/npres;
	  pixarr2[sgn][ctr] = histpos->second.tot2/npres;
	  ipixarr1[sgn][ctr] = 1.0/pixarr1[sgn][ctr];
	  ipixarr2[sgn][ctr] = 1.0/pixarr2[sgn][ctr];
	  ++ctr;
	}
	npix[sgn] = newnpix;

	has_weights[sgn] = true;
      } else {
	has_weights[sgn] = false;
	weights[sgn] = NULL;
      }
  } 
  pixsize = pixscale;
}

/*!
  \param[in] filename1   File to read from, beam 1
  \param[in] filename2   File to read from, beam 2
  \param[in] histogram Do beam histogramming?
  \param[in] histogramlogstep Log step if histogramming
*/
void doublebeam::readFiles(const std::string& filename1,
			   const std::string& filename2, bool histogram,
			   double histogramlogstep) {

  int status;
  fitsfile *fptr;
  char card[FLEN_CARD];
  int nkeys;

  status = 0;
  fptr = NULL;

  fits_open_file(&fptr, filename1.c_str(), READONLY, &status);
  if (fptr == NULL) {
    std::stringstream errstr("");
    errstr << "Error opening input file: " << filename1;
    fits_close_file(fptr,&status);
    throw affineExcept("doublebeam","readFiles",
		       errstr.str(),1);
  }
  if (status) {
    std::stringstream errstr("");
    errstr << "Error opening input file: " << filename1;
    fits_report_error(stderr,status);
    fits_close_file(fptr,&status);
    throw affineExcept("doublebeam","readFiles",
		       errstr.str(),2);
  }

  fits_get_hdrspace(fptr,&nkeys,NULL,&status);
  if (status) {
    std::stringstream errstr("");
    errstr << "Error getting header space for: " << filename1;
    fits_report_error(stderr,status);
    fits_close_file(fptr,&status);
    throw affineExcept("doublebeam","readFiles",
		       errstr.str(),4);
  }

  //Get the pixel scale
  double pixsize1;
  fits_read_key(fptr, TDOUBLE, const_cast<char*>("PIXSIZE"),
		&pixsize1,card,&status);
  if (status) {
    //Try pixscale
    status=0;
    fits_read_key(fptr, TDOUBLE, const_cast<char*>("PIXSCALE"),
		  &pixsize1,card,&status);
    if (status) {
      //Ok, give up
      std::stringstream errstr("");
      errstr << "Unable to find pixel scale information in "
	     << filename1;
      fits_report_error(stderr,status);
      fits_close_file(fptr,&status);
      throw affineExcept("doublebeam","readFiles",
			 errstr.str(),8);
    }
  }

  //Read the actual data
  int naxis;
  long naxes[2];
  fits_get_img_dim(fptr, &naxis, &status);
  fits_get_img_size(fptr, 2, naxes, &status);
  if (status || naxis != 2) {
    std::stringstream errstr("");
    errstr << "ERROR: input DOUBLEBEAM is not 2D";
    fits_close_file(fptr,&status);
    throw affineExcept("doublebeam","readFiles",
		       errstr.str(),8);
  }

  unsigned int n = naxes[0] * naxes[1];
  double* rpixarr1 = new double[n];
  long fpixel[2];
  fpixel[0] = 1;
  for (unsigned int i = 1; i <= static_cast<unsigned int>(naxes[1]); ++i) {
    fpixel[1] = i;
    if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], 0, 
		      rpixarr1 + naxes[0]*(i-1), 0, &status)) 
      break;   /* jump out of loop on error */
  }

  //Clean up
  fits_close_file(fptr,&status);
  if (status) {
    fits_report_error(stderr,status);
    delete[] rpixarr1;
    throw affineExcept("doublebeam","readFiles",
		       "Error closing file",16);
  } 

  //Read the second, just like the first
  status = 0;
  fits_open_file(&fptr, filename2.c_str(), READONLY, &status);
  if (fptr == NULL) {
    fits_close_file(fptr,&status);
    delete[] rpixarr1;
    std::stringstream errstr("");
    errstr << "Error opening input file: " << filename2;
    throw affineExcept("doublebeam","readFiles",
		       errstr.str(),32);
  }
  if (status) {
    fits_report_error(stderr,status);
    fits_close_file(fptr,&status);
    delete[] rpixarr1;
    std::stringstream errstr("");
    errstr << "Error opening input file: " << filename2;
    throw affineExcept("doublebeam","readFiles",
		       errstr.str(),64);
  }
  fits_get_hdrspace(fptr,&nkeys,NULL,&status);
  if (status) {
    fits_report_error(stderr,status);
    fits_close_file(fptr,&status);
    delete[] rpixarr1;
    std::stringstream errstr("");
    errstr << "Error getting header space for: " << filename2;
    throw affineExcept("doublebeam","readFiles",
		       errstr.str(),128);
  }
  double pixsize2;
  fits_read_key(fptr, TDOUBLE, const_cast<char*>("PIXSIZE"),
		&pixsize2,card,&status);
  if (status) {
    status=0;
    fits_read_key(fptr, TDOUBLE, const_cast<char*>("PIXSCALE"),
		  &pixsize2,card,&status);
    if (status) {
      fits_report_error(stderr,status);
      fits_close_file(fptr,&status);
      delete[] rpixarr1;
      std::stringstream errstr("");
      errstr << "Unable to find pixel scale information in "
	     << filename2;
      throw affineExcept("doublebeam","readFiles",
			 errstr.str(),256);
    }
  }
  if (fabs((pixsize1 - pixsize2) / pixsize1 ) > 0.01) {
    delete[] rpixarr1;
    fits_close_file(fptr,&status);
    std::stringstream errstr("");
    errstr << "Pixel scale in " << filename2 << " doesn't match that of "
	   << filename1;
    throw affineExcept("doublebeam","readFiles",
		       errstr.str(),512);
  }
  int naxis2;
  long naxes2[2];
  fits_get_img_dim(fptr, &naxis2, &status);
  fits_get_img_size(fptr, 2, naxes2, &status);
  if (status || naxis2 != 2) {
    fits_close_file(fptr,&status);
    delete[] rpixarr1;
    std::stringstream errstr("");
    errstr << "ERROR: input DOUBLEBEAM is not 2D in " << filename2;
    throw affineExcept("doublebeam","readFiles",
		       errstr.str(),1024);
  }
  if (naxes2[0] != naxes[0]) {

    fits_close_file(fptr,&status);
    delete[] rpixarr1;
    std::stringstream errstr("");
    errstr << "Dimension 1 extent for " << filename2 << " doesn't match "
	      << "that of " << filename1;
    throw affineExcept("doublebeam","readFiles",
		       errstr.str(),2048);
  }
  if (naxes2[1] != naxes[1]) {
    fits_close_file(fptr,&status);
    delete[] rpixarr1;
    std::stringstream errstr("");
    errstr << "Dimension 2 extent for " << filename2 << " doesn't match "
	      << "that of " << filename1;
    throw affineExcept("doublebeam","readFiles",
		       errstr.str(),2048);
  }
  double* rpixarr2 = new double[n];
  for (unsigned int i = 1; i <= static_cast<unsigned int>(naxes[1]); ++i) {
    fpixel[1] = i;
    if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], 0, 
		      rpixarr2 + naxes[0]*(i-1), 0, &status)) 
      break;   /* jump out of loop on error */
  }
  fits_close_file(fptr,&status);
  if (status) {
    fits_report_error(stderr,status);
    delete[] rpixarr1; delete[] rpixarr2;
    std::stringstream errstr("");
    errstr << "Error closing DOUBLEBEAM file " << filename2;
    throw affineExcept("doublebeam","readFiles",
		       errstr.str(),4096);
  } 

  setBeams(n, rpixarr1, rpixarr2, pixsize1, histogram, histogramlogstep);

  delete[] rpixarr1;
  delete[] rpixarr2;
}

/*!
  \returns Effective area of beam in square deg, band 1
*/
double doublebeam::getEffectiveArea1() const {
  bool hasval = hassign[0];
  for (unsigned int i = 1; i < 4; ++i) hasval |= hassign[i];
  if (!hasval) return 0.0;
  double convfac = pixsize / 3600.0;
  return totsm1 * convfac * convfac;
}

/*!
  \returns Effective area of beam in square deg, band 2
*/
double doublebeam::getEffectiveArea2() const {
  bool hasval = hassign[0];
  for (unsigned int i = 1; i < 4; ++i) hasval |= hassign[i];
  if (!hasval) return 0.0;
  double convfac = pixsize / 3600.0;
  return totsm2 * convfac * convfac;
}

/*!
  \returns Geometric mean of effective area of beams
*/
double doublebeam::getEffectiveAreaPixGeoMean() const {
 bool hasval = hassign[0];
  for (unsigned int i = 1; i < 4; ++i) hasval |= hassign[i];
  if (!hasval) return 0.0;
  return sqrt(totsm1*totsm2);
}

/*!
  \returns Minimum area of any beam sign component
*/
double doublebeam::getMinAreaPix() const {
 bool hasval = hassign[0];
  for (unsigned int i = 1; i < 4; ++i) hasval |= hassign[i];
  if (!hasval) return 0.0;
  return totsm1 < totsm2 ? totsm1 : totsm2;
}

/*!
  \param[in] idx Sign component (pp, pn, np, nn)
  \returns Effective area of beam, band 1
*/
double doublebeam::getEffectiveAreaSign1(unsigned int idx) const {
  if (!hassign[idx]) return 0.0;
  double convfac = pixsize / 3600.0;
  double area = tot1[idx] * convfac * convfac;
  return area;
}

/*!
  \param[in] idx Sign component (pp, pn, np, nn)
  \returns Effective area of beam, band 2
*/
double doublebeam::getEffectiveAreaSign2(unsigned int idx) const {
  if (!hassign[idx]) return 0.0;
  double convfac = pixsize / 3600.0;
  double area = tot2[idx] * convfac * convfac;
  return area;
}

/*!
  \param[in] idx Sign component (pp, pn, np, nn)
  \returns Effective area of squared beam, band 1
*/
double doublebeam::getEffectiveAreaSqSign1(unsigned int idx) const {
  if (!hassign[idx]) return 0.0;
  double convfac = pixsize / 3600.0;
  double area = totsq1[idx] * convfac * convfac;
  return area;
}

/*!
  \param[in] idx Sign component (pp, pn, np, nn)
  \returns Effective area of squared beam, band 2
*/
double doublebeam::getEffectiveAreaSqSign2(unsigned int idx) const {
  if (!hassign[idx]) return 0.0;
  double convfac = pixsize / 3600.0;
  double area = totsq2[idx] * convfac * convfac;
  return area;
}

/*!
  \param[in] sgn Index into [pp,pn,np,nn] bits of beam 1
  \returns Maximum beam value in band 1
*/
double doublebeam::getMax1(unsigned int sgn) const {
  if (!hassign[sgn])
    return std::numeric_limits<double>::quiet_NaN();
  double* parr = pixarr1[sgn];
  double max, val;
  max = parr[0];
  for (unsigned int i = 0; i < npix[sgn]; ++i) {
    val = parr[i];
    if (val > max) max = val;
  }
  return max;
}

/*!
  \param[in] sgn Index into [pp,pn,np,nn] bits of beam 2
  \returns Maximum beam value in band 2
*/
double doublebeam::getMax2(unsigned int sgn) const {
  if (!hassign[sgn])
    return std::numeric_limits<double>::quiet_NaN();
  double* parr = pixarr2[sgn];
  double max, val;
  max = parr[0];
  for (unsigned int i = 0; i < npix[sgn]; ++i) {
    val = parr[i];
    if (val > max) max = val;
  }
  return max;
}

/*!
  \param[in] bnd Index into [pp,pn,np,nn] bits of beam 1
  \param[out] min Minimum value
  \param[out] max Maximum value
  
  Returns NaNs if the beam doesn't have that piece
*/
void doublebeam::getMinMax1(unsigned int bnd, double& min,
			    double& max) const {
  if (!hassign[bnd]) {
    min = std::numeric_limits<double>::quiet_NaN();
    max = std::numeric_limits<double>::quiet_NaN();
    return;
  }
  
  double* parr = pixarr1[bnd];
  min = max = parr[0];
  double val;
  for (unsigned int i = 1; i < npix[bnd]; ++i) {
    val = parr[i];
    if (val < min) min = val;
    if (val > max) max = val;
  }

}

/*!
  \param[in] bnd Index into [pp,pn,np,nn] bits of beam 2
  \param[out] min Minimum value
  \param[out] max Maximum value
  
  Returns NaN if the beam doesn't have that piece
*/
void doublebeam::getMinMax2(unsigned int bnd, double& min,
			    double& max) const {
  if (!hassign[bnd]) {
    min = std::numeric_limits<double>::quiet_NaN();
    max = std::numeric_limits<double>::quiet_NaN();
    return;
  }
  
  double* parr = pixarr2[bnd];
  min = max = parr[0];
  double val;
  for (unsigned int i = 1; i < npix[bnd]; ++i) {
    val = parr[i];
    if (val < min) min = val;
    if (val > max) max = val;
  }

}

/*!
  \param[in] comm MPI communicator
  \param[in] dest Destination for messages
*/
void doublebeam::sendSelf(MPI_Comm comm, int dest) const {
  MPI_Send(const_cast<double*>(&pixsize), 1, MPI_DOUBLE, dest, 
	   pofd_mcmc::DOUBLEBEAMSENDPIXSIZE, comm);
  MPI_Send(const_cast<unsigned int*>(npix), 4, MPI_UNSIGNED, 
	   dest, pofd_mcmc::DOUBLEBEAMSENDN, comm);
  MPI_Send(const_cast<bool*>(has_weights), 4, MPI::BOOL, dest, 
	   pofd_mcmc::DOUBLEBEAMSENDHASWEIGHTS, comm);
  
  for (unsigned int i = 0; i < 4; ++i)
    if (npix[i] > 0) {
      MPI_Send(pixarr1[i], npix[i], MPI_DOUBLE, dest,
	       pofd_mcmc::DOUBLEBEAMSENDPIXARR1, comm);
      MPI_Send(pixarr2[i], npix[i], MPI_DOUBLE, dest,
	       pofd_mcmc::DOUBLEBEAMSENDPIXARR2, comm);
      MPI_Send(ipixarr1[i], npix[i], MPI_DOUBLE, dest,
	       pofd_mcmc::DOUBLEBEAMSENDIPIXARR1, comm);
      MPI_Send(ipixarr2[i], npix[i], MPI_DOUBLE, dest,
	       pofd_mcmc::DOUBLEBEAMSENDIPIXARR2, comm);
      if (has_weights[i])
	MPI_Send(weights[i], npix[i], MPI_DOUBLE, dest,
		 pofd_mcmc::DOUBLEBEAMSENDWEIGHTS, comm);
    }
  MPI_Send(const_cast<double*>(tot1), 4, MPI_DOUBLE, dest, 
	   pofd_mcmc::DOUBLEBEAMSENDTOT1, comm);
  MPI_Send(const_cast<double*>(tot2), 4, MPI_DOUBLE, dest, 
	   pofd_mcmc::DOUBLEBEAMSENDTOT2, comm);
  MPI_Send(const_cast<double*>(totsq1), 4, MPI_DOUBLE, dest, 
	   pofd_mcmc::DOUBLEBEAMSENDTOTSQ1, comm);
  MPI_Send(const_cast<double*>(totsq2), 4, MPI_DOUBLE, dest, 
	   pofd_mcmc::DOUBLEBEAMSENDTOTSQ2, comm);
  MPI_Send(const_cast<double*>(&totsm1), 1, MPI_DOUBLE, dest, 
	   pofd_mcmc::DOUBLEBEAMSENDTOTSM1, comm);
  MPI_Send(const_cast<double*>(&totsm2), 1, MPI_DOUBLE, dest, 
	   pofd_mcmc::DOUBLEBEAMSENDTOTSM2, comm);
}

/*!
  \param[in] comm MPI communicator
  \param[in] src Source of messages
*/
void doublebeam::recieveCopy(MPI_Comm comm, int src) {
  MPI_Status Info;

  unsigned int new_n[4];
  MPI_Recv(&pixsize, 1, MPI_DOUBLE, src, pofd_mcmc::DOUBLEBEAMSENDPIXSIZE,
	   comm, &Info);

  //Number of elems
  MPI_Recv(new_n, 4, MPI_UNSIGNED, src, pofd_mcmc::DOUBLEBEAMSENDN,
	   comm, &Info);
  
  //Weight elements
  MPI_Recv(&has_weights, 4, MPI::BOOL, src, 
	   pofd_mcmc::DOUBLEBEAMSENDHASWEIGHTS, comm, &Info);

  //Actual beam elems
  for (unsigned int i = 0; i < 4; ++i) {
    if (new_n[i] != npix[i]) {
      if (pixarr1[i] != NULL) delete[] pixarr1[i];
      if (pixarr2[i] != NULL) delete[] pixarr2[i];
      if (ipixarr1[i] != NULL) delete[] ipixarr1[i];
      if (ipixarr2[i] != NULL) delete[] ipixarr2[i];
      if (weights[i] != NULL) delete[] weights[i];
      if (new_n[i] > 0) {
	pixarr1[i] = new double[ new_n[i] ];
	pixarr2[i] = new double[ new_n[i] ];
	ipixarr1[i] = new double[ new_n[i] ];
	ipixarr2[i] = new double[ new_n[i] ];
	if (has_weights[i]) weights[i] = new double[new_n[i]]; else
	  weights[i] = NULL;
      } else {
	pixarr1[i] = pixarr2[i] = NULL;
	ipixarr1[i] = ipixarr2[i] = NULL;
	weights[i] = NULL;
      }
    } else {
      if (has_weights[i]) {
	if (weights[i] == NULL)
	  weights[i] = new double[ new_n[i] ];
      } else if (weights[i] != NULL) {
	delete[] weights[i];
	weights[i] = NULL;
      }
    }
    
    if (new_n[i] > 0) {
      MPI_Recv(pixarr1[i], new_n[i], MPI_DOUBLE, src,
	       pofd_mcmc::DOUBLEBEAMSENDPIXARR1, comm, &Info);
      MPI_Recv(pixarr2[i], new_n[i], MPI_DOUBLE, src,
	       pofd_mcmc::DOUBLEBEAMSENDPIXARR2, comm, &Info);
      MPI_Recv(ipixarr1[i], new_n[i], MPI_DOUBLE, src,
	       pofd_mcmc::DOUBLEBEAMSENDIPIXARR1, comm, &Info);
      MPI_Recv(ipixarr2[i], new_n[i], MPI_DOUBLE, src,
	       pofd_mcmc::DOUBLEBEAMSENDIPIXARR2, comm, &Info);
      if (has_weights[i])
	MPI_Recv(weights[i], new_n[i], MPI_DOUBLE, src,
		 pofd_mcmc::DOUBLEBEAMSENDWEIGHTS, comm, &Info);
    } 

    npix[i] = new_n[i];
    if (npix[i] > 0) hassign[i]=true; else hassign[i]=false;
  }
  
  //Totals
  MPI_Recv(tot1, 4, MPI_DOUBLE, src, pofd_mcmc::DOUBLEBEAMSENDTOT1,
	   comm, &Info);
  MPI_Recv(tot2, 4, MPI_DOUBLE, src, pofd_mcmc::DOUBLEBEAMSENDTOT2,
	   comm, &Info);
  MPI_Recv(totsq1, 4, MPI_DOUBLE, src, pofd_mcmc::DOUBLEBEAMSENDTOTSQ1,
	   comm, &Info);
  MPI_Recv(totsq2, 4, MPI_DOUBLE, src, pofd_mcmc::DOUBLEBEAMSENDTOTSQ2,
	   comm, &Info);
  MPI_Recv(&totsm1, 1, MPI_DOUBLE, src, pofd_mcmc::DOUBLEBEAMSENDTOTSM1,
	   comm, &Info);
  MPI_Recv(&totsm2, 1, MPI_DOUBLE, src, pofd_mcmc::DOUBLEBEAMSENDTOTSM2,
	   comm, &Info);
}

