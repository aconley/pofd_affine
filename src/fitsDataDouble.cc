#include<limits>
#include<cstring>

#include<fitsio.h>
#include<fftw3.h>

#include "../include/fitsDataDouble.h"
#include "../include/affineExcept.h"
#include "../include/global_settings.h"

fitsDataDouble::fitsDataDouble() {
  n = 0;
  data1 = data2 = nullptr;
  is_binned = false;
  nbins1 = nbins2 = 0;
  bincent01 = bindelta1 = bincent02 = bindelta2 = 0.0;
  binval = nullptr;
  file1 = "";
  file2 = "";
  dataext1 = dataext2 = maskext1 = maskext2 = 0;
  has_mask1 = has_mask2 = false;
}

/*!
  \param[in] filename1  Name of data file, band 1
  \param[in] filename2  Name of data file, band 2
  \param[in] ignoremask Don't consider MASK information from files
  \param[in] meansub Do a mean subtraction

  Like the 1D code, this assumes that any mask is stored in unsigned
  integer form -- important because the FITS standard doesn't fully
  support this.  So the code will not work correctly if your mask
  is stored as integers
*/
fitsDataDouble::fitsDataDouble(const std::string& filename1, 
			       const std::string& filename2,
			       bool ignoremask, bool meansub) {
  n = 0;
  data1 = data2 = nullptr;
  is_binned = false;
  nbins1 = nbins2 = 0;
  bincent01 = bindelta1 = bincent02 = bindelta2 = 0.0;
  binval = nullptr;
  readData(filename1, filename2, ignoremask, meansub);
}

fitsDataDouble::~fitsDataDouble() {
  if (data1 != nullptr) fftw_free(data1);
  if (data2 != nullptr) fftw_free(data2);
  if (binval != nullptr) fftw_free(binval);
}


/*
  \param[in] file Name of file to read
  \param[out] ndata Number of data points
  \param[out] data  On output, holds image pixels
  \param[out] dataext Data extension number
  \param[out] mask  On output, holds mask pixels
  \param[out] maskext Mask extension number
  \param[in] ignore_mask Don't consider mask information from files
  \returns True if a mask was read, false if not

  Internal reading function to avoid mass repetition

  It's up to the caller to not have stuff in data or mask before
  calling, or a memory leak will result
*/
bool fitsDataDouble::readFile(const std::string& file, long& ndata,
			      double*& data, unsigned int& dataext, 
			      unsigned int*& mask, unsigned int& maskext,
			      bool ignore_mask) {
  //Lots of fits reading fun
  fitsfile *fptr;
  int status, naxis, hduok, dataok;
  long nmask, naxes[2], fpixel[2], masknaxes[2];
  bool this_has_mask;

  this_has_mask = false;
  hduok = dataok = 0;
  data = nullptr;
  mask = nullptr;
  ndata = 0;

  //Here we go -- start by trying to open file
  status = 0;
  fits_open_file(&fptr, file.c_str(), READONLY, &status);
  if (status) {
    fits_report_error(stderr, status);
    fits_close_file(fptr, &status);
    throw affineExcept("fitsDataDouble", "readFile", "Problem opening file");
  }

  // Next, try to read in the data.  We do this before the mask
  // because the user may have specified an extension, so we
  // should start there.
  if (file.find("[") == std::string::npos) {
    // The caller did -not- specify an extension, so we
    //  have to try to find the right one.
    int hdutype = IMAGE_HDU;
    fits_movnam_hdu(fptr, hdutype, const_cast<char*>("image"),
		    0, &status);
    if (status == BAD_HDU_NUM) {
      status = 0;
      fits_movnam_hdu(fptr, hdutype, const_cast<char*>("IMAGE"),
		      0, &status);
    }
    if (status == BAD_HDU_NUM) {
      status = 0;
      fits_movabs_hdu(fptr, 1, &hdutype, &status);
      if (hdutype != IMAGE_HDU) {
	fits_close_file(fptr, &status);
	throw affineExcept("fitsDataDouble","readFile",
			   "Couldn't find image");
      }
    }
    if (status) {
      fits_report_error(stderr, status);
      fits_close_file(fptr, &status);
      throw affineExcept("fitsDataDouble", "readFile",
			 "Error locating image data");
    }
  }
  int st;
  fits_get_hdu_num(fptr, &st);
  dataext = static_cast<unsigned int>(st - 1);

  //Check checksum if available
  fits_verify_chksum(fptr, &dataok, &hduok, &status);
  if (dataok == -1) {
    fits_close_file(fptr, &status);
    throw affineExcept("fitsDataDouble", "readFile",
		       "Image data fails checksum");
  }
  if (hduok == -1) {
    fits_close_file(fptr, &status);
    throw affineExcept("fitsDataDouble", "readFile",
		       "Image HDU fails checksum");
  }

  fits_get_img_dim(fptr, &naxis, &status);
  if (status) {
    fits_close_file(fptr, &status);
    throw affineExcept("fitsDataDouble", "readFile",
		       "Error getting image dimensions");
  }
  if (naxis != 2) {
    fits_close_file(fptr, &status);
    throw affineExcept("fitsDataDouble", "readFile",
		       "Image is not 2D");
  }


  fits_get_img_size(fptr, 2, naxes, &status);
  if (status) {
    fits_close_file(fptr, &status);
    throw affineExcept("fitsDataDouble", "readFile",
		       "Error getting image dimensions");
  }

  ndata = static_cast<unsigned int>(naxes[0] * naxes[1]);
  if (ndata == 0) {
    fits_close_file(fptr, &status);
    throw affineExcept("fitsDataDouble", "readFile", "Image is of zero extent");
  }

  //Read, after allocating
  data = (double*) fftw_malloc(sizeof(double) * ndata);
  fpixel[0] = 1; fpixel[1] = 1;
  fits_read_pix(fptr, TDOUBLE, fpixel, ndata, 0, data, 0, &status);
  if (status) {
    fits_report_error(stderr, status);
    fits_close_file(fptr, &status);
    fftw_free(data);
    throw affineExcept("fitsDataDouble", "readFile", "Error reading image");
  }

  // Now try to find and read in the mask, if present
  if (!ignore_mask) {
    fits_movnam_hdu(fptr, IMAGE_HDU, const_cast<char*>("mask"),
		    0, &status);
    if (status == 0) this_has_mask = true; else {
      if (status != BAD_HDU_NUM) {
	fftw_free(data);
	fits_report_error(stderr, status);
	fits_close_file(fptr, &status);
	throw affineExcept("fitsDataDouble", "readFile", 
			   "Error searching for mask");
      } else {
	//Try uppercase
	fits_movnam_hdu(fptr, IMAGE_HDU, const_cast<char*>("MASK"),
			0, &status);
	if (status == 0) this_has_mask = true; else
	  if (status != BAD_HDU_NUM) {
	    fftw_free(data);
	    fits_report_error(stderr, status);
	    fits_close_file(fptr, &status);
	    throw affineExcept("fitsDataDouble", "readFile", 
			       "Error searching for mask");
	  }
      }
    }
    fits_get_hdu_num(fptr, &st);
    maskext = static_cast<unsigned int>(st - 1);
  } else maskext = 0;

  //Read in mask if present
  status = 0;
  if (this_has_mask) {
    //Check checksum if available
    fits_verify_chksum(fptr, &dataok, &hduok, &status);
    if (dataok == -1) {
      fits_close_file(fptr, &status);
      fftw_free(data);
      throw affineExcept("fitsDataDouble", "readFile",
			 "Data in mask extension fails checksum test");
    }
    if (hduok == -1) {
      fits_close_file(fptr, &status);
      fftw_free(data);
      throw affineExcept("fitsDataDouble", "readFile",
			 "Mask extension HDU fails checksum test");
    }
    
    fits_get_img_dim(fptr, &naxis, &status);
    if (status) {
      fits_close_file(fptr, &status);
      fftw_free(data);
      throw affineExcept("fitsDataDouble", "readFile",
			 "Error getting mask dimensions");
    }
    if (naxis != 2) {
      fits_close_file(fptr, &status);
      fftw_free(data);
      throw affineExcept("fitsDataDouble", "readFile",
			 "Mask present but not 2D");
    }

    fits_get_img_size(fptr, 2, masknaxes, &status);
    if (status) {
      fits_close_file(fptr, &status);
      fftw_free(data);
      throw affineExcept("fitsDataDouble", "readFile",
			 "Mask failure getting mask extents");
    }
    
    nmask = masknaxes[0] * masknaxes[1];
    if (nmask == 0) {
      fits_close_file(fptr, &status);
      fftw_free(data);
      throw affineExcept("fitsDataDouble", "readFile",
			 "Mask present zero extent");
    } 

    if (naxes[0] != masknaxes[0]) {
      fits_close_file(fptr, &status);
      fftw_free(data);
      throw affineExcept("fitsDataDouble", "readFile",
			 "Mask does not match image extent along first dimension");
    }
    if (masknaxes[1] != naxes[1]) {
      fits_close_file(fptr, &status);
      fftw_free(data);
      throw affineExcept("fitsDataDouble", "readFile",
			 "Mask does not match image extent along second dimension");
    }
    
    //Read the mask, after allocating
    mask = (unsigned int*) fftw_malloc(sizeof(unsigned int) * nmask);
    fpixel[0] = 1; fpixel[1] = 1;
    fits_read_pix(fptr, TUINT, fpixel, nmask, 0, mask, 0, &status);
    if (status) {
      fits_report_error(stderr, status);
      fits_close_file(fptr, &status);
      fftw_free(data);
      fftw_free(mask);
      throw affineExcept("fitsDataDouble", "readFile",
			 "Error reading mask data");
    }
  }

  //Close out
  fits_close_file(fptr, &status);
  if (status) {
    fits_report_error(stderr, status);
    if (this_has_mask) fftw_free(mask);
    fftw_free(data);
    throw affineExcept("fitsDataDouble", "readFile", "Error closing file");
  }

  return this_has_mask;
}



/*!
  \param[in] filename1  Name of data file, band 1
  \param[in] filename2  Name of data file, band 2
  \param[in] ignore_mask Don't consider MASK information from files
  \param[in] meansub Do a mean subtraction

  The mean subtraction is based on unmasked pixels.
  This will blow away any binning.
*/
void fitsDataDouble::readData(const std::string& filename1,
			      const std::string& filename2,
			      bool ignore_mask, bool meansub) {
  //Free up old arrays
  if (data1 != nullptr) { fftw_free(data1); data1 = nullptr; }
  if (data2 != nullptr) { fftw_free(data2); data2 = nullptr; }
  if (binval != nullptr) { fftw_free(binval); binval = nullptr; }
  is_binned = false;
  has_mask1 = has_mask2 = false;
  maskext1 = maskext2 = 0;
  n = 0;

  file1 = filename1;
  file2 = filename2;

  //Read data into temporary arrays
  long ndat1_, ndat2_;
  double *data1_, *data2_;
  unsigned int *mask1_, *mask2_;
  data1_ = data2_ = nullptr;
  mask1_ = mask2_ = nullptr;
  has_mask1 = readFile(file1, ndat1_, data1_, dataext1, mask1_, maskext1, 
		       ignore_mask);
  has_mask2 = readFile(file2, ndat2_, data2_, dataext2, mask2_, maskext2,
		       ignore_mask);

  //Check to make sure they are the same size, etc.
  if (ndat1_ != ndat2_) {
    if (data1_ != nullptr) fftw_free(data1_);
    if (data2_ != nullptr) fftw_free(data2_);
    if (mask1_ != nullptr) fftw_free(mask1_);
    if (mask2_ != nullptr) fftw_free(mask2_);
    throw affineExcept("fitsDataDouble", "readData",
		       "Data from file1 and file2 not the same size");
  }
  
  //If neither has a mask, we can just copy directly
  if (!(has_mask1 || has_mask2)) {
    n = ndat1_;
    data1 = data1_;
    data2 = data2_;
  } else {
    if (has_mask1 && has_mask2) {
      //Both have masks
      //First, count the number of bits to copy
      unsigned int nkeep = 0;
      for (unsigned int i = 0; i < ndat1_; ++i)
	if ((mask1_[i] & mask2_[i]) == 0) ++nkeep;
      if (nkeep == 0) {
	if (data1_ != nullptr) { fftw_free(data1_); data1_ = nullptr; }
	if (data2_ != nullptr) { fftw_free(data2_); data2_ = nullptr; }
	if (mask1_ != nullptr) { fftw_free(mask1_); mask1_ = nullptr; }
	if (mask2_ != nullptr) { fftw_free(mask2_); mask2_ = nullptr; }
	throw affineExcept("fitsDataDouble","readData",
			   "All data masked");
      }
      //Now actually copy
      data1 = (double*) fftw_malloc(sizeof(double) * nkeep);
      data2 = (double*) fftw_malloc(sizeof(double) * nkeep);
      unsigned int cntr = 0;
      for (unsigned int i = 0; i < ndat1_; ++i)
	if ((mask1_[i] & mask2_[i]) == 0) {
	  data1[cntr] = data1_[i];
	  data2[cntr] = data2_[i];
	  ++cntr;
	}
      n = nkeep;
    } else {
      //Only one has a mask
      unsigned int *msk;
      if (has_mask1) msk = mask1_; else msk=mask2_;
      unsigned int nkeep = 0;
      for (unsigned int i = 0; i < ndat1_; ++i)
	if (msk[i] == 0) ++nkeep;
      if (nkeep == 0) {
	if (data1_ != nullptr) { fftw_free(data1_); data1_ = nullptr; }
	if (data2_ != nullptr) { fftw_free(data2_); data2_ = nullptr; }
	if (mask1_ != nullptr) { fftw_free(mask1_); mask1_ = nullptr; }
	if (mask2_ != nullptr) { fftw_free(mask2_); mask2_ = nullptr; }
	throw affineExcept("fitsDataDouble", "readData",
			   "All data masked");
      }
      data1 = (double*) fftw_malloc(sizeof(double) * nkeep);
      data2 = (double*) fftw_malloc(sizeof(double) * nkeep);
      unsigned int cntr = 0;
      for (unsigned int i = 0; i < ndat1_; ++i)
	if ( msk[i] == 0 ) {
	  data1[cntr] = data1_[i];
	  data2[cntr] = data2_[i];
	  ++cntr;
	}
      n = nkeep;
    }
    fftw_free(data1_);
    data1_ = nullptr;
    fftw_free(data2_);
    data2_ = nullptr;
  }
  if (mask1_ != nullptr) { fftw_free(mask1_); mask1_ = nullptr; }
  if (mask2_ != nullptr) { fftw_free(mask2_); mask2_ = nullptr; }
  if (meansub) meanSubtract();

}

/*!
  \returns Minimum values in each band.  NaN if no data read
*/
dblpair fitsDataDouble::getMin() const {
  dblpair retval;
  if (n == 0) {
    retval.first  = std::numeric_limits<double>::quiet_NaN();
    retval.second = std::numeric_limits<double>::quiet_NaN();
    return retval;
  }

  //There is a shortcut if we are binned, since bincent0
  // is the minimum value by construction
  if (is_binned) {
    retval.first = bincent01;
    retval.second = bincent02;
  } else {
    double minval, cval;
    minval = data1[0];
    for (unsigned int i = 1; i < n; ++i) {
      cval = data1[i];
      if (cval < minval) minval=cval;
    }
    retval.first = minval;
    minval = data2[0];
    for (unsigned int i = 1; i < n; ++i) {
      cval = data2[i];
      if (cval < minval) minval=cval;
    }
    retval.second = minval;
  }
  return retval;
}

/*!
  \returns Maximum values in each band.  NaN if no data read
*/
dblpair fitsDataDouble::getMax() const {
  dblpair retval;
  if (n == 0) {
    retval.first  = std::numeric_limits<double>::quiet_NaN();
    retval.second = std::numeric_limits<double>::quiet_NaN();
    return retval;
  }
  //There is a shortcut if we are binned, since 
  // bincent0 + (nbins - 1)*bindelta
  // is the maximum value by construction
  if (is_binned) {
    retval.first  = bincent01 + static_cast<double>(nbins1-1)*bindelta1;
    retval.second = bincent02 + static_cast<double>(nbins2-1)*bindelta2;
  } else {
    double maxval, cval;
    maxval = data1[0];
    for (unsigned int i = 1; i < n; ++i) {
      cval = data1[i];
      if (cval > maxval) maxval=cval;
    }
    retval.first = maxval;
    maxval = data2[0];
    for (unsigned int i = 1; i < n; ++i) {
      cval = data2[i];
      if (cval > maxval) maxval=cval;
    }
    retval.second = maxval;
  }
  return retval;
}

/*!
  \returns A pair of min/max pairs.
*/
std::pair<dblpair, dblpair>
fitsDataDouble::getMinMax() const {
  double min1, max1, min2, max2;
  if (n == 0) {
    min1 = max1 = min2 = max2 = std::numeric_limits<double>::quiet_NaN();
    return std::make_pair(std::make_pair(min1, max1),
			  std::make_pair(min2, max2));
  }
  
  if (is_binned) {
    min1 = bincent01;
    max1 = bincent01 + static_cast<double>(nbins1 - 1) * bindelta1;
    min2 = bincent02;
    max2 = bincent02 + static_cast<double>(nbins2 - 1) * bindelta2;
    return std::make_pair(std::make_pair(min1, max1),
			  std::make_pair(min2, max2));
  }

  double cval;
  min1 = max1 = data1[0];
  for (unsigned int i = 1; i < n; ++i) {
    cval = data1[i];
    if (cval < min1) min1 = cval;
    if (cval > max1) max1 = cval;
  }
  min2 = max2 = data2[0];
  for (unsigned int i = 1; i < n; ++i) {
    cval = data2[i];
    if (cval < min2) min2 = cval;
    if (cval > max2) max2 = cval;
  }
  return std::make_pair(std::make_pair(min1, max1),
			std::make_pair(min2, max2));
}

/*!
  \returns Mean values in each band.  NaN if no data read
*/
dblpair fitsDataDouble::getMean() const {
  dblpair retval;
  if (n == 0) {
    retval.first  = std::numeric_limits<double>::quiet_NaN();
    retval.second = std::numeric_limits<double>::quiet_NaN();
    return retval;
  }
  double mnval;
  mnval = data1[0];
  for (unsigned int i=1; i < n; ++i) mnval += data1[i];
  retval.first = mnval / static_cast<double>(n);
  mnval = data2[0];
  for (unsigned int i=1; i < n; ++i) mnval += data2[i];
  retval.second = mnval / static_cast<double>(n);
  return retval;
}

/*!
  \returns Mean values in each band before subtraction.  NaN if no data read
*/
dblpair fitsDataDouble::meanSubtract() {
  dblpair mnval = getMean();
  if (n == 0) return mnval;

  double cmn = mnval.first;
  for (unsigned int i = 0; i < n; ++i) data1[i] -= cmn;
  if (is_binned) bincent01 -= cmn;

  cmn = mnval.second;
  for (unsigned int i = 0; i < n; ++i) data2[i] -= cmn;
  if (is_binned) bincent02 -= cmn;

  return mnval;
}

/*!						
  \param[in] NBINS1 Number of bins, band 1
  \param[in] NBINS2 Number of bins, band 2
*/
void fitsDataDouble::applyBinning(unsigned int NBINS1, unsigned int NBINS2) {
  if (n == 0) return;
  if (NBINS1 == 0) throw affineExcept("fitsDataDouble", "applyBinning",
				      "Trying to bin with no bins, band 1");
  if (NBINS2 == 0) throw affineExcept("fitsDataDouble", "applyBinning",
				      "Trying to bin with no bins, band 2");
  if (is_binned && (NBINS1 == nbins1) && (NBINS2 == nbins2)) return;

  if ((NBINS1 != nbins1) || (NBINS2 != nbins2)) {
    if (binval != nullptr) fftw_free(binval);
    binval =
      (unsigned int*) fftw_malloc(sizeof(unsigned int) *NBINS1 * NBINS2);
    nbins1 = NBINS1;
    nbins2 = NBINS2;
  }
  std::memset(binval, 0, nbins1 * nbins2 * sizeof(unsigned int));

  //First, we need the minimum and maximum
  std::pair<dblpair, dblpair> minmax = getMinMax();
  if (std::isnan(minmax.first.first) || std::isinf(minmax.first.first) ||
      std::isnan(minmax.first.second) || std::isinf(minmax.first.second) ||
      std::isnan(minmax.second.first) || std::isinf(minmax.second.first) ||
      std::isnan(minmax.second.second) || std::isinf(minmax.second.second))
    throw affineExcept("fitsDataDouble", "applyBinning",
		       "Non-valid min or max");

  //We want to put the max and min in the center of the top and
  // bottom bin (not at the edges -- a somewhat arbitrary choice)
  bincent01 = minmax.first.first;
  if (nbins1 == 1)
    bindelta1 = 2 * (minmax.first.second - minmax.first.first);
  else
    bindelta1 = (minmax.first.second - minmax.first.first) /
      static_cast<double>(nbins1 - 1);
  bincent02 = minmax.second.first;
  if (nbins2 == 1)
    bindelta2 = 2 * (minmax.second.second - minmax.second.first);
  else
    bindelta2 = (minmax.second.second - minmax.second.first) /
      static_cast<double>(nbins2-1);

  //And... bin
  double ibindelta1 = 1.0 / bindelta1;
  double ibindelta2 = 1.0 / bindelta2;
  unsigned int idx1, idx2;
  for (unsigned int i = 0; i < n; ++i) {
    idx1 = static_cast<unsigned int>((data1[i]-bincent01)*ibindelta1 + 0.5);
    idx2 = static_cast<unsigned int>((data2[i]-bincent02)*ibindelta2 + 0.5);
    binval[idx1 * nbins2 + idx2] += 1;
  }
  is_binned = true;
}

void fitsDataDouble::removeBinning() {
  if (!is_binned) return;
  if (binval != nullptr) { fftw_free(binval); binval = nullptr; }
  nbins1 = nbins2 = 0;
  is_binned = false;
}

/*!
  \returns Number of bins in each band.  Does not imply that binning
  has been applied.
*/
std::pair<unsigned int, unsigned int> fitsDataDouble::getNBins() const {
  return std::make_pair(nbins1, nbins2);
}

/*!
  \returns Center of bottom bin in each band.  Does not imply that binning
  has been applied.
*/
dblpair fitsDataDouble::getBinCent0() const {
  return std::make_pair(bincent01, bincent02);
}

/*!
  \returns Size of bins in each band.  Does not imply that binning
  has been applied.
*/
dblpair fitsDataDouble::getBinDelta() const {
  return std::make_pair(bindelta1, bindelta2);
}

/*!
  \param[in] comm MPI communicator
  \param[in] dest Destination of messages
*/
void fitsDataDouble::sendSelf(MPI_Comm comm, int dest) const {
  MPI_Send(const_cast<unsigned int*>(&n), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::FDDSENDN, comm);
  if (n > 0) {
    MPI_Send(data1, n, MPI_DOUBLE, dest, pofd_mcmc::FDDSENDDATA1, comm);
    MPI_Send(data2, n, MPI_DOUBLE, dest, pofd_mcmc::FDDSENDDATA2, comm);
    MPI_Send(const_cast<bool*>(&is_binned), 1, MPI::BOOL, dest,
	     pofd_mcmc::FDDSENDISBINNED, comm);
    if (is_binned) {
      MPI_Send(const_cast<unsigned int*>(&nbins1), 1, MPI_UNSIGNED, dest,
	       pofd_mcmc::FDDSENDNBINS1, comm);
      MPI_Send(const_cast<unsigned int*>(&nbins2), 1, MPI_UNSIGNED, dest,
	       pofd_mcmc::FDDSENDNBINS2, comm);
      if ((nbins1 > 0) && (nbins2 > 0)) {
	MPI_Send(binval, nbins1*nbins2, MPI_UNSIGNED, dest,
		 pofd_mcmc::FDDSENDBINVAL, comm);
	MPI_Send(const_cast<double*>(&bincent01), 1, MPI_DOUBLE, dest,
		 pofd_mcmc::FDDSENDBINCENT01, comm);
	MPI_Send(const_cast<double*>(&bindelta1), 1, MPI_DOUBLE, dest,
		 pofd_mcmc::FDDSENDBINDELTA1, comm);
	MPI_Send(const_cast<double*>(&bincent02), 1, MPI_DOUBLE, dest,
		 pofd_mcmc::FDDSENDBINCENT02, comm);
	MPI_Send(const_cast<double*>(&bindelta2), 1, MPI_DOUBLE, dest,
		 pofd_mcmc::FDDSENDBINDELTA2, comm);
      }
    }
  }
}

/*!
  \param[in] comm MPI communicator
  \param[in] src Source of messages
*/
void fitsDataDouble::receiveCopy(MPI_Comm comm, int src) {
  unsigned int newn;
  MPI_Status Info;

  MPI_Recv(&newn, 1, MPI_UNSIGNED, src, pofd_mcmc::FDDSENDN, comm, &Info);

  if (newn != n) {
    //Have to change sizes
    if (data1 != nullptr) fftw_free(data1);
    if (data2 != nullptr) fftw_free(data2);
    if (newn > 0) {
      data1 = (double*) fftw_malloc(sizeof(double) * newn);
      data2 = (double*) fftw_malloc(sizeof(double) * newn);
    } else data1 = data2 = nullptr;
    n = newn;
  }
  if (n > 0) {
    //Data to receive
    MPI_Recv(data1, n, MPI_DOUBLE, src, pofd_mcmc::FDDSENDDATA1, comm, &Info);
    MPI_Recv(data2, n, MPI_DOUBLE, src, pofd_mcmc::FDDSENDDATA2, comm, &Info);

    bool recbin;
    MPI_Recv(&recbin, 1, MPI::BOOL, src, pofd_mcmc::FDDSENDISBINNED, 
	     comm, &Info);
    if (recbin) {
      unsigned int newnbins1, newnbins2;
      MPI_Recv(&newnbins1, 1, MPI_UNSIGNED, src, pofd_mcmc::FDDSENDNBINS1, 
	       comm, &Info);
      MPI_Recv(&newnbins2, 1, MPI_UNSIGNED, src, pofd_mcmc::FDDSENDNBINS2, 
	       comm, &Info);
      //Realloc if needed
      if ((newnbins1 != nbins1) || (newnbins2 != nbins2)) {
	if (binval != nullptr) fftw_free(binval);
	if ( (newnbins1 > 0) && (newnbins2 > 0) ) {
	  binval = (unsigned int*) 
	    fftw_malloc( sizeof(unsigned int)*newnbins1*newnbins2 );
	  nbins1 = newnbins1;
	  nbins2 = newnbins2;
	} else binval = nullptr;
      }
      if ((nbins1 > 0) && (nbins2 > 0)) {
	MPI_Recv(binval, nbins1*nbins2, MPI_UNSIGNED, src,
		 pofd_mcmc::FDDSENDBINVAL, comm, &Info);
	MPI_Recv(&bincent01, 1, MPI_DOUBLE, src, pofd_mcmc::FDDSENDBINCENT01,
		 comm, &Info);
	MPI_Recv(&bindelta1, 1, MPI_DOUBLE, src, pofd_mcmc::FDDSENDBINDELTA1,
		 comm, &Info);
	MPI_Recv(&bincent02, 1, MPI_DOUBLE, src, pofd_mcmc::FDDSENDBINCENT02,
		 comm, &Info);
	MPI_Recv(&bindelta2, 1, MPI_DOUBLE, src, pofd_mcmc::FDDSENDBINDELTA2,
		 comm, &Info);
	is_binned = true;
      } else is_binned = false;
    }
  } else {
    if (binval != nullptr) { fftw_free(binval); binval = nullptr; }
    is_binned = false;
    nbins1 = nbins2 = 0;
  }
}
