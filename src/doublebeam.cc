//doublebeam.cc

#include<limits>
#include<iostream>
#include<cmath>
#include<sstream>
#include<algorithm>
#include<cstring>

#include<fitsio.h>

#include "../include/doublebeam.h"
#include "../include/global_settings.h"
#include "../include/affineExcept.h"

const unsigned int doublebeam::histothresh = 15;

// Note we don't set nbins or minval here.  Why not?
// Because these only really make sense when reading or histogramming
doublebeam::doublebeam() {
  double nan = std::numeric_limits<double>::quiet_NaN();
  for (unsigned int i = 0; i < 4; ++i) hassign[i] = false;
  for (unsigned int i = 0; i < 4; ++i) npix[i] = 0;
  for (unsigned int i = 0; i < 4; ++i) pixarr1[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) pixarr2[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) invpixarr1[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) invpixarr2[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) haslogratio[i] = false;
  for (unsigned int i = 0; i < 4; ++i) logratio[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) ishistogrammed[i] = false;
  for (unsigned int i = 0; i < 4; ++i) nhist[i] = 0;
  for (unsigned int i = 0; i < 4; ++i) binweights[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) binvals1[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) binvals2[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) hasbinlogratio[i] = false;
  for (unsigned int i = 0; i < 4; ++i) binlogratio[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) tot1[i] = nan;
  for (unsigned int i = 0; i < 4; ++i) tot2[i] = nan;
  for (unsigned int i = 0; i < 4; ++i) totsq1[i] = nan;
  for (unsigned int i = 0; i < 4; ++i) totsq2[i] = nan;
  for (unsigned int i = 0; i < 4; ++i) minbm1[i] = nan;
  for (unsigned int i = 0; i < 4; ++i) maxbm1[i] = nan;
  for (unsigned int i = 0; i < 4; ++i) minbm2[i] = nan;
  for (unsigned int i = 0; i < 4; ++i) maxbm2[i] = nan;
}

/*!
  \param[in] filename1   File to read from, beam 1
  \param[in] filename2   File to read from, beam 2
  \param[in] histogram Do beam histogramming?
  \param[in] NBINS Number of histogram bins
  \param[in] MINVAL Minimum bean value to consider
*/
doublebeam::doublebeam(const std::string& filename1,
		       const std::string& filename2,
		       bool histogram, unsigned int NBINS,
		       double MINVAL) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  for (unsigned int i = 0; i < 4; ++i) hassign[i] = false;
  for (unsigned int i = 0; i < 4; ++i) npix[i] = 0;
  for (unsigned int i = 0; i < 4; ++i) pixarr1[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) pixarr2[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) invpixarr1[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) invpixarr2[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) haslogratio[i] = false;
  for (unsigned int i = 0; i < 4; ++i) logratio[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) ishistogrammed[i] = false;
  for (unsigned int i = 0; i < 4; ++i) nhist[i] = 0;
  for (unsigned int i = 0; i < 4; ++i) binweights[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) binvals1[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) binvals2[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) hasbinlogratio[i] = false;
  for (unsigned int i = 0; i < 4; ++i) binlogratio[i] = NULL;
  for (unsigned int i = 0; i < 4; ++i) tot1[i] = nan;
  for (unsigned int i = 0; i < 4; ++i) tot2[i] = nan;
  for (unsigned int i = 0; i < 4; ++i) totsq1[i] = nan;
  for (unsigned int i = 0; i < 4; ++i) totsq2[i] = nan;
  for (unsigned int i = 0; i < 4; ++i) minbm1[i] = nan;
  for (unsigned int i = 0; i < 4; ++i) maxbm1[i] = nan;
  for (unsigned int i = 0; i < 4; ++i) minbm2[i] = nan;
  for (unsigned int i = 0; i < 4; ++i) maxbm2[i] = nan;

  readFiles(filename1, filename2, MINVAL);
  if (histogram) makeHistogram(NBINS);
}

void doublebeam::cleanup() {
  double NaN = std::numeric_limits<double>::quiet_NaN();
  for (unsigned int i = 0; i < 4; ++i) {
    hassign[i] = false;
    npix[i] = 0;
    if (pixarr1[i] != NULL) { delete[] pixarr1[i]; pixarr1[i] = NULL; }
    if (pixarr2[i] != NULL) { delete[] pixarr2[i]; pixarr2[i] = NULL; }
    if (invpixarr1[i] != NULL) { delete[] invpixarr1[i]; invpixarr1[i] = NULL; }
    if (invpixarr2[i] != NULL) { delete[] invpixarr2[i]; invpixarr2[i] = NULL; }
    haslogratio[i] = false;
    if (logratio[i] != NULL) { delete[] logratio[i]; logratio[i] = NULL; }
    ishistogrammed[i] = false;
    nhist[i] = 0;
    if (binweights[i] != NULL) { delete[] binweights[i]; binweights[i] = NULL; }
    if (binvals1[i] != NULL) { delete[] binvals1[i]; binvals1[i] = NULL; }
    if (binvals2[i] != NULL) { delete[] binvals2[i]; binvals2[i] = NULL; }
    hasbinlogratio[i] = false;
    if (binlogratio[i] != NULL) { delete[] binlogratio[i]; binlogratio[i] = NULL; }
    tot1[i] = tot2[i] = totsq1[i] = totsq2[i] = NaN;
    minbm1[i] = maxbm1[i] = minbm2[i] = maxbm2[i] = NaN;
  }
}

/*!
  \param[in] filename1   File to read from, beam 1
  \param[in] filename2   File to read from, beam 2
  \param[in] MINVAL      Minimum absolute value.  Specifically, only
          beam elements with fabs(val) > minval in both beams are kept.
*/
void doublebeam::readFiles(const std::string& filename1,
			   const std::string& filename2,
			   double MINVAL) {

  // Two stages:
  //  1) read in the FITS files
  //  2) Process them into the pixel arrays
  // This function does the FITS reading, then uses setBeams to
  // do step 2, since it's conceivable somebody might want to do
  // the second by hand.

  // Read in the FITS files
  int status;
  fitsfile *fptr;
  char card[FLEN_CARD];
  int nkeys;

  status = 0;
  fptr = NULL;

  // Start with beam 1
  fits_open_file(&fptr, filename1.c_str(), READONLY, &status);
  if (status) {
    std::stringstream errstr("");
    errstr << "Error opening input file: " << filename1;
    fits_report_error(stderr,status);
    fits_close_file(fptr,&status);
    throw affineExcept("doublebeam", "readFiles", errstr.str());
  }

  fits_get_hdrspace(fptr, &nkeys, NULL, &status);
  if (status) {
    std::stringstream errstr("");
    errstr << "Error getting header space for: " << filename1;
    fits_report_error(stderr, status);
    fits_close_file(fptr, &status);
    throw affineExcept("doublebeam", "readFiles", errstr.str());
  }

  //Get the pixel scale.  We will check later to make sure the other beam
  // has the same size.
  double pixsize1;
  fits_read_key(fptr, TDOUBLE, const_cast<char*>("PIXSIZE"),
		&pixsize1, card, &status);
  if (status) { // Nope, didn't find it.  Try pixscale instead
    status=0;
    fits_read_key(fptr, TDOUBLE, const_cast<char*>("PIXSCALE"),
		  &pixsize1, card, &status);
    if (status) {
      //Didn't find that either
      std::stringstream errstr("");
      errstr << "Unable to find pixel scale information in "
	     << filename1;
      fits_report_error(stderr, status);
      fits_close_file(fptr, &status);
      throw affineExcept("doublebeam", "readFiles", errstr.str());
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
    fits_close_file(fptr, &status);
    throw affineExcept("doublebeam", "readFiles", errstr.str());
  }

  unsigned int n = naxes[0] * naxes[1];
  double* rpixarr1 = new double[n]; // Temporary reading array
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
    fits_report_error(stderr, status);
    delete[] rpixarr1;
    throw affineExcept("doublebeam", "readFiles", "Error closing file");
  } 

  //Read the second file
  status = 0;
  fits_open_file(&fptr, filename2.c_str(), READONLY, &status);
  if (status) {
    fits_report_error(stderr, status);
    fits_close_file(fptr, &status);
    delete[] rpixarr1;
    std::stringstream errstr("");
    errstr << "Error opening input file: " << filename2;
    throw affineExcept("doublebeam", "readFiles", errstr.str());
  }
  fits_get_hdrspace(fptr, &nkeys, NULL, &status);
  if (status) {
    fits_report_error(stderr, status);
    fits_close_file(fptr, &status);
    delete[] rpixarr1;
    std::stringstream errstr("");
    errstr << "Error getting header space for: " << filename2;
    throw affineExcept("doublebeam", "readFiles", errstr.str());
  }
  double pixsize2; // Make sure it's the same as band 1
  fits_read_key(fptr, TDOUBLE, const_cast<char*>("PIXSIZE"),
		&pixsize2, card, &status);
  if (status) {
    status=0;
    fits_read_key(fptr, TDOUBLE, const_cast<char*>("PIXSCALE"),
		  &pixsize2, card, &status);
    if (status) {
      fits_report_error(stderr, status);
      fits_close_file(fptr, &status);
      delete[] rpixarr1;
      std::stringstream errstr("");
      errstr << "Unable to find pixel scale information in "
	     << filename2;
      throw affineExcept("doublebeam", "readFiles", errstr.str());
    }
  }

  if (fabs((pixsize1 - pixsize2) / pixsize1 ) > 1e-5) {
    delete[] rpixarr1;
    fits_close_file(fptr, &status);
    std::stringstream errstr("");
    errstr << "Pixel scale in " << filename2 << " doesn't match that of "
	   << filename1 << ": " << pixsize1 << " vs " << pixsize2;
    throw affineExcept("doublebeam", "readFiles", errstr.str());
  }

  int naxis2;
  long naxes2[2];
  fits_get_img_dim(fptr, &naxis2, &status);
  fits_get_img_size(fptr, 2, naxes2, &status);
  if (status || naxis2 != 2) {
    fits_close_file(fptr, &status);
    delete[] rpixarr1;
    std::stringstream errstr("");
    errstr << "ERROR: input DOUBLEBEAM is not 2D in " << filename2;
    throw affineExcept("doublebeam", "readFiles", errstr.str());
  }

  if (naxes2[0] != naxes[0]) {
    fits_close_file(fptr,&status);
    delete[] rpixarr1;
    std::stringstream errstr("");
    errstr << "Dimension 1 extent for " << filename2 << " doesn't match "
	      << "that of " << filename1;
    throw affineExcept("doublebeam", "readFiles", errstr.str());
  }
  if (naxes2[1] != naxes[1]) {
    fits_close_file(fptr, &status);
    delete[] rpixarr1;
    std::stringstream errstr("");
    errstr << "Dimension 2 extent for " << filename2 << " doesn't match "
	      << "that of " << filename1;
    throw affineExcept("doublebeam", "readFiles", errstr.str());
  }
  double* rpixarr2 = new double[n]; // Another temp reading array
  for (unsigned int i = 1; i <= static_cast<unsigned int>(naxes[1]); ++i) {
    fpixel[1] = i;
    if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], 0, 
		      rpixarr2 + naxes[0]*(i-1), 0, &status)) 
      break;   /* jump out of loop on error */
  }
  fits_close_file(fptr,&status);
  if (status) {
    fits_report_error(stderr, status);
    delete[] rpixarr1; delete[] rpixarr2;
    std::stringstream errstr("");
    errstr << "Error closing DOUBLEBEAM file " << filename2;
    throw affineExcept("doublebeam", "readFiles", errstr.str());
  } 

  // Ok, finally all read!  So push into the beam arrays
  setBeams(n, rpixarr1, rpixarr2, pixsize1, MINVAL);

  delete[] rpixarr1;
  delete[] rpixarr2;
}

/*!
  \param[in] n Number of elements in beam1, beam2
  \param[in] beam1 Band 1 beam, size n
  \param[in] beam2 Band 2 beam, size n
  \param[in] PIXSIZE pixel size in arcsec.
  \param[in] MINVAL      Minimum absolute value.  Specifically, only
          beam elements with fabs(val) > minval in both beams are kept.
*/
void doublebeam::setBeams(unsigned int n, const double* const beam1, 
			  const double* const beam2, double PIXSIZE,
			  double MINVAL) {
  // This does the sorting into beam components, as well as allocating
  // internal arrays.  So -- first -- clean up any old info
  cleanup();

  if (n == 0)
    throw affineExcept("doublebeam", "setBeams", "Invalid (zero) n");

  if (MINVAL < 0.0)
    throw affineExcept("doublebeam", "setBeams", "Invalid (negative) MINVAL");

  minval = MINVAL;

  // First, figure out which component is going where
  // Use a temporary array of the sign, where 0, 1, 2, 3 indicate
  // the pp, pn, np, nn component, 4 means it will not be kept.
  unsigned int *component;
  component = new unsigned int[n];
  unsigned int sgn;
  double val1, val2, aval1, aval2;
  for (unsigned int i = 0; i < n; ++i) {
    val1 = beam1[i];
    val2 = beam2[i];
    aval1 = fabs(val1);
    aval2 = fabs(val2);
    if ((aval1 > minval) && (aval2 > minval)) { // Note it's > not >= so 0 works
      if (val1 > 0) {
	if (val2 > 0) sgn=0; else sgn=1;
      } else {
	if (val2 > 0) sgn=2; else sgn=3;
      }
      component[i] = sgn;
      npix[sgn] += 1;
    } else {
      // Not keeping
      component[i] = 4;
    }
  }

  unsigned int totpix = npix[0];
  for (unsigned int i = 1; i < 4; ++i) totpix += npix[i];
  if (totpix == 0)
    throw affineExcept("doublebeam", "setBeams", "No acceptable values found");

  // Now set the beam pieces
  unsigned int idx, curr_n;
  double val;
  double *p1, *p2;
  for (unsigned int i = 0; i < 4; ++i) {
    curr_n = npix[i];
    if (curr_n > 0) {
      hassign[i] = true; // All false from cleanup
      pixarr1[i] = new double[npix[i]];
      pixarr2[i] = new double[npix[i]];
      invpixarr1[i] = new double[npix[i]];
      invpixarr2[i] = new double[npix[i]];

      p1 = pixarr1[i];
      p2 = pixarr2[i];

      // Loop and set values for this component
      idx = 0;
      for (unsigned int j = 0; j < n; ++j)
	if (component[j] == i) {
	  aval1 = fabs(beam1[j]);
	  aval2 = fabs(beam2[j]);
	  p1[idx] = aval1;
	  p2[idx] = aval2;
	  ++idx;
	}
      //assert(idx == curr_n);

      // Do totals
      val = p1[0];
      tot1[i] = val;
      totsq1[i] = val * val;
      for (unsigned int j = 1; j < curr_n; ++j) {
	val = p1[j];
	tot1[i] += val;
	totsq1[i] += val * val;
      }
      val = p2[0];
      tot2[i] = val;
      totsq2[i] = val * val;
      for (unsigned int j = 1; j < curr_n; ++j) {
	val = p2[j];
	tot2[i] += val;
	totsq2[i] += val * val;
      }
      
      // Inverses
      for (unsigned int j = 0; j < curr_n; ++j)
	invpixarr1[i][j] = 1.0 / p1[j];
      for (unsigned int j = 0; j < curr_n; ++j)
	invpixarr2[i][j] = 1.0 / p2[j];
    }
  }

  delete[] component;

  // Find the min/max values for each component
  double curr_min, curr_max;
  for (unsigned int i = 0; i < 4; ++i) 
    if (hassign[i]) {
      curr_n = npix[i];
      p1 = pixarr1[i];
      curr_min = curr_max = p1[0];
      for (unsigned int j = 1; j < curr_n; ++j) {
	val = p1[j];
	if (val > curr_max) curr_max = val;
	else if (val < curr_min) curr_min = val;
      }
      minbm1[i] = curr_min;
      maxbm1[i] = curr_max;
      p2 = pixarr2[i];
      curr_min = curr_max = p2[0];
      for (unsigned int j = 1; j < curr_n; ++j) {
	val = p2[j];
	if (val > curr_max) curr_max = val;
	else if (val < curr_min) curr_min = val;
      }
      minbm2[i] = curr_min;
      maxbm2[i] = curr_max;
    } // min/max bm[12] already NaNed by cleanup

  if (PIXSIZE <= 0.0)
    throw affineExcept("doublebeam", "setBeams",
		       "Invalid (non-positive) pixel size");
  pixsize = PIXSIZE;
}

void doublebeam::makeHistogram(unsigned int NBINS) {
  if (!hasData())
    throw affineExcept("doublebeam", "makeHistogram",
		       "Trying to histogram empty beam");
  if (NBINS == 0)
    throw affineExcept("doublebeam", "makeHistogram",
		       "Invalid (zero) NBINS");
  nbins = NBINS;

  // Clean out any old values
  for (unsigned int i = 0; i < 4; ++i) {
    ishistogrammed[i] = false;
    nhist[i] = 0;
    if (binweights[i] != NULL) { delete[] binweights[i]; binweights[i] = NULL; }
    if (binvals1[i] != NULL) { delete[] binvals1[i]; binvals1[i] = NULL; }
    if (binvals2[i] != NULL) { delete[] binvals2[i]; binvals2[i] = NULL; }
  }

  // Working variables
  unsigned int nbins2 = nbins * nbins;
  unsigned int *tmpwt;
  double *tmphist1, *tmphist2;
  tmpwt = new unsigned int[nbins2];
  tmphist1 = new double[nbins2];
  tmphist2 = new double[nbins2];

  // Main loop
  // We bin in log2 space.  And, even though we are binning
  // the inverse beam, we bin on the non-inverse one then invert,
  // which seems to work better in practice (in that a smaller number
  // of bins produces a R calculation that is closer to the unbinned value).
  const double log2outscale = 0.0014419741739063218; // log2(1.001)
  unsigned int idx1, idx2, totidx, curr_n, utmp;
  double dnbins = static_cast<double>(nbins);
  std::pair<double, double> tmp; // For getting min/max
  double *p1, *p2, *p3; // Temporary convenience pointers
  double ihiststep1, ihiststep2; // Inverse histogram bin size (in log2 space)
  double minbinval1, maxbinval1, minbinval2, maxbinval2; // Min/max vals
  double val1, val2, wtnorm; // Working vars
  for (unsigned int sgn = 0; sgn < 4; ++sgn) 
    if (hassign[sgn] && npix[sgn] >= histothresh) {

      // Get min/max to figure out minimum bin value and step
      tmp = getMinMax1(sgn);
      minbinval1 = log2(tmp.first) - log2outscale;
      maxbinval1 = log2(tmp.second) + log2outscale;
      ihiststep1 = dnbins / (maxbinval1 - minbinval1);
      tmp = getMinMax2(sgn);
      minbinval2 = log2(tmp.first) - log2outscale;
      maxbinval2 = log2(tmp.second) + log2outscale;
      ihiststep2 = dnbins / (maxbinval2 - minbinval2);

      // Zero temporaries
      std::memset(tmpwt, 0, nbins2 * sizeof(unsigned int));
      std::memset(tmphist1, 0, nbins2 * sizeof(double));
      std::memset(tmphist2, 0, nbins2 * sizeof(double));

      p1 = pixarr1[sgn];
      p2 = pixarr2[sgn];
      curr_n = npix[sgn];
      
      for (unsigned int i = 0; i < curr_n; ++i) {
	val1 = p1[i]; // Recall -- already fabsed
	val2 = p2[i]; // also fabsed
	idx1 = static_cast<unsigned int>((log2(val1) - minbinval1) * 
					 ihiststep1);
	// assert(idx1 < nbins);
	idx2 = static_cast<unsigned int>((log2(val2) - minbinval2) * 
					 ihiststep2);
	// assert(idx2 < nbins);
	totidx = idx1 * nbins + idx2;
	tmpwt[totidx] += 1;
	tmphist1[totidx] += val1;
	tmphist2[totidx] += val2;
      }

      // Count the number of non-zero bins
      unsigned int n_nonzero;
      n_nonzero = 0;
      for (unsigned int i = 0; i < nbins2; ++i)
	if (tmpwt[i] > 0) ++n_nonzero; 
      // assert(n_nonzero > 0);
      
      // Allocate the final arrays
      binweights[sgn] = new double[n_nonzero];
      binvals1[sgn] = new double[n_nonzero];
      binvals2[sgn] = new double[n_nonzero];
      
      totidx = 0;
      p1 = binweights[sgn];
      p2 = binvals1[sgn];
      p3 = binvals2[sgn];
      for (unsigned int i = 0; i < nbins2; ++i) {
	utmp = tmpwt[i];
	if (utmp > 0) {
	  wtnorm = static_cast<double>(utmp);
	  p1[totidx] = wtnorm;
	  p2[totidx] = wtnorm / tmphist1[i]; // Invert and normalize
	  p3[totidx] = wtnorm / tmphist2[i];
	  ++totidx;
	}
      }
      ishistogrammed[sgn] = true;
      nhist[sgn] = n_nonzero;
    }

  delete[] tmpwt;
  delete[] tmphist1;
  delete[] tmphist2;
}

/*!
  \returns True if the beam has any data loaded, false otherwise
*/
bool doublebeam::hasData() const {
  for (unsigned int i = 0; i < 4; ++i)
    if (hassign[i]) return true;
  return false;
}

/*!
  \returns Effective area of beam in square deg, band 1
*/
double doublebeam::getEffectiveArea1() const {
  if (!hasData()) return std::numeric_limits<double>::quiet_NaN();
  double convfac = pixsize / 3600.0;
  double totarea = 0.0;
  for (unsigned int i = 0; i < 4; ++i)
    if (hassign[i]) totarea += tot1[i];
  return totarea * convfac * convfac;
}

/*!
  \returns Effective area of beam in square deg, band 2
*/
double doublebeam::getEffectiveArea2() const {
  if (!hasData()) return std::numeric_limits<double>::quiet_NaN();
  double convfac = pixsize / 3600.0;
  double totarea = 0.0;
  for (unsigned int i = 0; i < 4; ++i)
    if (hassign[i]) totarea += tot2[i];
  return totarea * convfac * convfac;
}

/*!
  \param[in] idx Sign component (pp, pn, np, nn)
  \returns Effective area of beam, band 1, in square degrees
*/
double doublebeam::getEffectiveAreaSign1(unsigned int idx) const {
  if (idx >= 4) throw affineExcept("doublebeam", "getEffectiveAreaSign1",
				   "Invalid index");
  
  if (!hassign[idx]) return 0.0;
  double convfac = pixsize / 3600.0;
  double area = tot1[idx] * convfac * convfac;
  return area;
}

/*!
  \param[in] idx Sign component (pp, pn, np, nn)
  \returns Effective area of beam, band 2, in square degrees
*/
double doublebeam::getEffectiveAreaSign2(unsigned int idx) const {
  if (idx >= 4) throw affineExcept("doublebeam", "getEffectiveAreaSign2",
				   "Invalid index");
  
  if (!hassign[idx]) return 0.0;
  double convfac = pixsize / 3600.0;
  double area = tot2[idx] * convfac * convfac;
  return area;
}

/*!
  \param[in] idx Sign component (pp, pn, np, nn)
  \returns Effective area of squared beam, band 1, in square degrees
*/
double doublebeam::getEffectiveAreaSqSign1(unsigned int idx) const {
  if (idx >= 4) throw affineExcept("doublebeam", "getEffectiveAreaSqSign1",
				   "Invalid index");
  
  if (!hassign[idx]) return 0.0;
  double convfac = pixsize / 3600.0;
  double area = totsq1[idx] * convfac * convfac;
  return area;
}

/*!
  \param[in] idx Sign component (pp, pn, np, nn)
  \returns Effective area of squared beam, band 2, in square degrees
*/
double doublebeam::getEffectiveAreaSqSign2(unsigned int idx) const {
  if (idx >= 4) throw affineExcept("doublebeam", "getEffectiveAreaSqSign2",
				   "Invalid index");
  
  if (!hassign[idx]) return 0.0;
  double convfac = pixsize / 3600.0;
  double area = totsq2[idx] * convfac * convfac;
  return area;
}

unsigned int doublebeam::getMaxNPix() const {
  unsigned int n = 0;
  for (unsigned int i = 0; i < 4; ++i)
    if (hassign[i] && npix[i] > n) n = npix[i];
  return n;
}

unsigned int doublebeam::getTotalNPix() const {
  unsigned int n = 0;
  for (unsigned int i = 0; i < 4; ++i)
    if (hassign[i]) n += npix[i];
  return n;
}

/*!
  \param[in] idx Sign component (pp, pn, np, nn)
  \returns Minimum and maximum values for band 1 beam in this sign component
*/
dblpair doublebeam::getMinMax1(unsigned int idx) const {
  if (idx >= 4) throw affineExcept("doublebeam", "getMinMax1",
				   "Invalid index");
  if (!hassign[idx]) 
    return std::make_pair(std::numeric_limits<double>::quiet_NaN(),
			  std::numeric_limits<double>::quiet_NaN());
  return std::make_pair(minbm1[idx], maxbm1[idx]);
}

/*!
  \param[in] idx Sign component (pp, pn, np, nn)
  \returns Minimum and maximum values for band 2 beam in this sign component
*/
dblpair doublebeam::getMinMax2(unsigned int idx) const {
  if (idx >= 4) throw affineExcept("doublebeam", "getMinMax2",
				   "Invalid index");
  if (!hassign[idx]) 
    return std::make_pair(std::numeric_limits<double>::quiet_NaN(),
			  std::numeric_limits<double>::quiet_NaN());
  return std::make_pair(minbm2[idx], maxbm2[idx]);
}

/*
!
  \param[in] idx Sign component (pp, pn, np, nn)
  \returns Access to Log beam1 / beam2.

  Computed if not already set, and saved for later
*/
const double* const doublebeam::getLogRatio(unsigned int idx) const {
  if (idx >= 4) throw affineExcept("doublebeam", "getLogRatio",
				   "Invalid index");
  if (!hassign[idx]) throw affineExcept("doublebeam", "getLogRatio",
					"No data in specified index");
  if (!haslogratio[idx]) {
    // Compute
    unsigned int n = npix[idx];
    if (logratio[idx] != NULL) delete[] logratio[idx];
    logratio[idx] = new double[n];
    const double* p1 = pixarr1[idx];  // beam1
    const double* p2 = invpixarr2[idx];  // 1 / beam2
    double* lg = logratio[idx];
    // Note that zero beam elements are not kept
    for (unsigned int i = 0; i < n; ++i)
      lg[i] = log(p1[i] * p2[i]);
    haslogratio[idx] = true;
  }
  return logratio[idx];
}

/*
!
  \param[in] idx Sign component (pp, pn, np, nn)
  \returns Access to Log beam1 / beam2 of histogrammed beam.

  Computed if not already set, and saved for later
*/
const double* const doublebeam::getBinLogRatio(unsigned int idx) const {
  if (idx >= 4) throw affineExcept("doublebeam", "getBinLogRatio",
				   "Invalid index");
  if (!hassign[idx]) throw affineExcept("doublebeam", "getLogRatio",
					"No data in specified index");
  if (!ishistogrammed[idx]) throw affineExcept("doublebeam", "getLogRatio",
					       "No histogram in that index");
  if (!hasbinlogratio[idx]) {
    // Compute
    unsigned int n = nhist[idx];
    if (binlogratio[idx] != NULL) delete[] binlogratio[idx];
    binlogratio[idx] = new double[n];
    double* lg = binlogratio[idx];
    const double* p1 = binvals1[idx];  // 1 / beam1
    const double* p2 = binvals2[idx];  // 1 / beam2
    // Note that zero beam elements are not kept,
    //  and that binvals stores the inverse beam
    for (unsigned int i = 0; i < n; ++i)
      lg[i] = log(p2[i] / p1[i]); // beam 1 / beam 2
    hasbinlogratio[idx] = true;
  }
  return binlogratio[idx];
}

/*!
  \param[in] comm MPI communicator
  \param[in] dest Destination for messages
*/
void doublebeam::sendSelf(MPI_Comm comm, int dest) const {
  MPI_Send(const_cast<double*>(&pixsize), 1, MPI_DOUBLE, dest, 
	   pofd_mcmc::DOUBLEBEAMSENDPIXSIZE, comm);
  MPI_Send(const_cast<double*>(&minval), 1, MPI_DOUBLE, dest,
	   pofd_mcmc::DOUBLEBEAMSENDMINVAL, comm);

  // Send raw beam
  MPI_Send(const_cast<bool*>(hassign), 4, MPI::BOOL, dest,
	   pofd_mcmc::DOUBLEBEAMSENDHASSIGN, comm);
  MPI_Send(const_cast<double*>(minbm1), 4, MPI_DOUBLE, dest,
	   pofd_mcmc::DOUBLEBEAMSENDMIN1, comm);
  MPI_Send(const_cast<double*>(maxbm1), 4, MPI_DOUBLE, dest,
	   pofd_mcmc::DOUBLEBEAMSENDMAX1, comm);
  MPI_Send(const_cast<double*>(minbm2), 4, MPI_DOUBLE, dest,
	   pofd_mcmc::DOUBLEBEAMSENDMIN2, comm);
  MPI_Send(const_cast<double*>(maxbm2), 4, MPI_DOUBLE, dest,
	   pofd_mcmc::DOUBLEBEAMSENDMAX2, comm);
  MPI_Send(const_cast<unsigned int*>(npix), 4, MPI_UNSIGNED, dest,
	   pofd_mcmc::DOUBLEBEAMSENDNPIX, comm);
  MPI_Send(const_cast<bool*>(haslogratio), 4, MPI::BOOL, dest,
	   pofd_mcmc::DOUBLEBEAMSENDHASLOGRATIO, comm);
  for (unsigned int i = 0; i < 4; ++i)
    if (hassign[i]) {
      MPI_Send(pixarr1[i], npix[i], MPI_DOUBLE, dest,
	       pofd_mcmc::DOUBLEBEAMSENDPIXARR1, comm);
      MPI_Send(pixarr2[i], npix[i], MPI_DOUBLE, dest,
	       pofd_mcmc::DOUBLEBEAMSENDPIXARR2, comm);
      MPI_Send(invpixarr1[i], npix[i], MPI_DOUBLE, dest,
	       pofd_mcmc::DOUBLEBEAMSENDIPIXARR1, comm);
      MPI_Send(invpixarr2[i], npix[i], MPI_DOUBLE, dest,
	       pofd_mcmc::DOUBLEBEAMSENDIPIXARR2, comm);
      if (haslogratio[i])
	MPI_Send(logratio[i], npix[i], MPI_DOUBLE, dest,
		 pofd_mcmc::DOUBLEBEAMSENDLOGRATIO, comm);
    } 
  
  // Histogram
  MPI_Send(const_cast<unsigned int*>(&nbins), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::DOUBLEBEAMSENDNBINS, comm);
  MPI_Send(const_cast<bool*>(ishistogrammed), 4, MPI::BOOL, dest,
	   pofd_mcmc::DOUBLEBEAMSENDISHISTOGRAMMED, comm);
  MPI_Send(const_cast<unsigned int*>(nhist), 4, MPI_UNSIGNED, dest,
	   pofd_mcmc::DOUBLEBEAMSENDNHIST, comm);
  MPI_Send(const_cast<bool*>(hasbinlogratio), 4, MPI::BOOL, dest,
	   pofd_mcmc::DOUBLEBEAMSENDHASBINLOGRATIO, comm);
  for (unsigned int i = 0; i < 4; ++i)
    if (ishistogrammed[i] && (nhist[i] > 0)) {
      MPI_Send(binweights[i], nhist[i], MPI_DOUBLE, dest,
	       pofd_mcmc::DOUBLEBEAMSENDBINWEIGHTS, comm);
      MPI_Send(binvals1[i], nhist[i], MPI_DOUBLE, dest,
	       pofd_mcmc::DOUBLEBEAMSENDBINVALS1, comm);
      MPI_Send(binvals2[i], nhist[i], MPI_DOUBLE, dest,
	       pofd_mcmc::DOUBLEBEAMSENDBINVALS2, comm);
      if (hasbinlogratio[i])
	MPI_Send(binlogratio[i], nhist[i], MPI_DOUBLE, dest,
		 pofd_mcmc::DOUBLEBEAMSENDBINLOGRATIO, comm);
    }

  // Descriptive
  MPI_Send(const_cast<double*>(tot1), 4, MPI_DOUBLE, dest, 
	   pofd_mcmc::DOUBLEBEAMSENDTOT1, comm);
  MPI_Send(const_cast<double*>(tot2), 4, MPI_DOUBLE, dest, 
	   pofd_mcmc::DOUBLEBEAMSENDTOT2, comm);
  MPI_Send(const_cast<double*>(totsq1), 4, MPI_DOUBLE, dest, 
	   pofd_mcmc::DOUBLEBEAMSENDTOTSQ1, comm);
  MPI_Send(const_cast<double*>(totsq2), 4, MPI_DOUBLE, dest, 
	   pofd_mcmc::DOUBLEBEAMSENDTOTSQ2, comm);
}

/*!
  \param[in] comm MPI communicator
  \param[in] src Source of messages
*/
void doublebeam::receiveCopy(MPI_Comm comm, int src) {

  cleanup(); // Get rid of all info.

  MPI_Status Info;
  MPI_Recv(&pixsize, 1, MPI_DOUBLE, src, pofd_mcmc::DOUBLEBEAMSENDPIXSIZE,
	   comm, &Info);
  MPI_Recv(&minval, 1, MPI_DOUBLE, src, pofd_mcmc::DOUBLEBEAMSENDMINVAL, 
	   comm, &Info);

  // Raw beam
  MPI_Recv(hassign, 4, MPI::BOOL, src, pofd_mcmc::DOUBLEBEAMSENDHASSIGN, 
	   comm, &Info);
  MPI_Recv(minbm1, 4, MPI_DOUBLE, src, pofd_mcmc::DOUBLEBEAMSENDMIN1, 
	   comm, &Info);
  MPI_Recv(maxbm1, 4, MPI_DOUBLE, src, pofd_mcmc::DOUBLEBEAMSENDMAX1, 
	   comm, &Info);
  MPI_Recv(minbm2, 4, MPI_DOUBLE, src, pofd_mcmc::DOUBLEBEAMSENDMIN2, 
	   comm, &Info);
  MPI_Recv(maxbm2, 4, MPI_DOUBLE, src, pofd_mcmc::DOUBLEBEAMSENDMAX2, 
	   comm, &Info);
  MPI_Recv(npix, 4, MPI_UNSIGNED, src, pofd_mcmc::DOUBLEBEAMSENDNPIX, 
	   comm, &Info);
  MPI_Recv(haslogratio, 4, MPI::BOOL, src, 
	   pofd_mcmc::DOUBLEBEAMSENDHASLOGRATIO, comm, &Info);

  // Raw beam
  for (unsigned int i = 0; i < 4; ++i)
    if (hassign[i]) {
      pixarr1[i] = new double[npix[i]];
      pixarr2[i] = new double[npix[i]];
      invpixarr1[i] = new double[npix[i]];
      invpixarr2[i] = new double[npix[i]];
      if (haslogratio[i])
	logratio[i] = new double[npix[i]];

      MPI_Recv(pixarr1[i], npix[i], MPI_DOUBLE, src,
	       pofd_mcmc::DOUBLEBEAMSENDPIXARR1, comm, &Info);
      MPI_Recv(pixarr2[i], npix[i], MPI_DOUBLE, src,
	       pofd_mcmc::DOUBLEBEAMSENDPIXARR2, comm, &Info);
      MPI_Recv(invpixarr1[i], npix[i], MPI_DOUBLE, src,
	       pofd_mcmc::DOUBLEBEAMSENDIPIXARR1, comm, &Info);
      MPI_Recv(invpixarr2[i], npix[i], MPI_DOUBLE, src,
	       pofd_mcmc::DOUBLEBEAMSENDIPIXARR2, comm, &Info);
      if (haslogratio[i])
	MPI_Recv(logratio[i], npix[i], MPI_DOUBLE, src,
		 pofd_mcmc::DOUBLEBEAMSENDLOGRATIO, comm, &Info);
    }
  
  // Histogrammed beam
  MPI_Recv(&nbins, 1, MPI_UNSIGNED, src,
	   pofd_mcmc::DOUBLEBEAMSENDNBINS, comm, &Info);
  MPI_Recv(ishistogrammed, 4, MPI::BOOL, src,
	   pofd_mcmc::DOUBLEBEAMSENDISHISTOGRAMMED, comm, &Info);
  MPI_Recv(nhist, 4, MPI_UNSIGNED, src,
	   pofd_mcmc::DOUBLEBEAMSENDNHIST, comm, &Info);
  MPI_Recv(hasbinlogratio, 4, MPI::BOOL, src, 
	   pofd_mcmc::DOUBLEBEAMSENDHASBINLOGRATIO, comm, &Info);
  for (unsigned int i = 0; i < 4; ++i)
    if (ishistogrammed[i] && (nhist[i] > 0)) {
      binweights[i] = new double[nhist[i]];
      binvals1[i] = new double[nhist[i]];
      binvals2[i] = new double[nhist[i]];
      if (hasbinlogratio[i])
	binlogratio[i] = new double[nhist[i]];
      MPI_Recv(binweights[i], nhist[i], MPI_DOUBLE, src,
	       pofd_mcmc::DOUBLEBEAMSENDBINWEIGHTS, comm, &Info);
      MPI_Recv(binvals1[i], nhist[i], MPI_DOUBLE, src,
	       pofd_mcmc::DOUBLEBEAMSENDBINVALS1, comm, &Info);
      MPI_Recv(binvals2[i], nhist[i], MPI_DOUBLE, src,
	       pofd_mcmc::DOUBLEBEAMSENDBINVALS2, comm, &Info);
      if (hasbinlogratio[i])
	MPI_Recv(binlogratio[i], nhist[i], MPI_DOUBLE, src,
		 pofd_mcmc::DOUBLEBEAMSENDBINLOGRATIO, comm, &Info);
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
}

/*!
  \param[inout] os Stream to write to
*/
bool doublebeam::writeToStream(std::ostream& os) const {
  const char* cmpnames[4] = {"pp", "pn", "np", "nn"};
  os << "Beam summary: " << std::endl;
  os << " Pixel size: " << pixsize << " Min kept: " << minval << std::endl;
  os << " Components info: " << std::endl;
  for (unsigned int i = 0; i < 4; ++i)
    if (hassign[i]) {
      os << " Component: " << cmpnames[i];
      os << " raw npix: " << npix[i];
      if (ishistogrammed[i]) 
	os << " hist bins: " << nhist[i];
      os << std::endl;
      os << "   Total1: " << tot1[i] << " Min1: " << minbm1[i]
	 << " Max1: " << maxbm1[i] << std::endl;
      os << "   Total2: " << tot2[i] << " Min2: " << minbm2[i]
	 << " Max2: " << maxbm2[i] << std::endl;
    }
  return true;
}

/*!
  \param[inout] os Stream to write to
  \param[in] b Number counts model
*/
std::ostream& operator<<(std::ostream& os, const doublebeam& b) {
  b.writeToStream(os);
  return os;
}
