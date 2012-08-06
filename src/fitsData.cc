#include<limits>

#include <fitsio.h>
#include <fftw3.h>

#include <fitsData.h>
#include <affineExcept.h>
#include <global_settings.h>

fitsData::fitsData() {
  n = 0;
  data=NULL;
  is_binned = false;
  nbins = 0;
  bincent0 = bindelta = 0.0;
  binval = NULL;
}

/*!
  \param[in] file  Name of data file
  \param[in] ignoremask Don't consider MASK information from files
  \param[in] meansub Do a mean subtraction
 */
fitsData::fitsData(const std::string& file, bool ignoremask, bool meansub) {
  n = 0;
  data=NULL;
  is_binned = false;
  nbins = 0;
  bincent0 = bindelta = 0.0;
  binval = NULL;
  readData(file,ignoremask,meansub);
}

fitsData::~fitsData() {
  if (data != NULL) fftw_free(data);
  if (binval != NULL) fftw_free(binval);
}

/*!
  \param[in] file  Name of data file
  \param[in] ignore_mask Don't consider MASK information from files
  \param[in] meansub Do a mean subtraction

  The mean subtraction is based on unmasked pixels.
  This will blow away any binning.
 */
void fitsData::readData(const std::string& file,
			bool ignore_mask, bool meansub) {

  //Free up old arrays
  if (data != NULL) fftw_free(data);
  n = 0;
  if (binval != NULL) fftw_free(binval);
  nbins = 0;
  is_binned = false;

  //Read data into temporary arrays so we can apply mask
  bool has_mask;
  unsigned int ndat_;
  double *data_;
  int *mask;
  ndat_ = 0;
  data_ = NULL;
  mask  = NULL;

   //Lots of fits reading fun
  fitsfile *fptr;
  int status, naxis, hduok, dataok;
  long nmask, naxes[2], fpixel[2], masknaxes[2];
  
  has_mask = false;
  hduok = dataok = 0;

  //Here we go -- start by trying to open file
  status = 0;
  fits_open_file(&fptr, file.c_str(), READONLY, &status);
  if (status) {
    fits_report_error(stderr,status);
    fits_close_file(fptr,&status);
    throw affineExcept("fitsData","readData","Problem opening file",1);
  }

  //Try to find a mask extension
  nmask = 0;
  if ( ! ignore_mask ) {
    fits_movnam_hdu(fptr, IMAGE_HDU, const_cast<char*>("mask"),
		    0, &status);
    if (status == 0) has_mask = true; else {
      if (status != BAD_HDU_NUM) {
	fits_report_error(stderr, status);
	fits_close_file(fptr, &status);
	throw affineExcept("fitsData","readData","Error looking for mask",2);
      } else {
	//Try uppercase
	fits_movnam_hdu(fptr, IMAGE_HDU, const_cast<char*>("MASK"),
			0, &status);
	if (status == 0) has_mask = true; else
	  if (status != BAD_HDU_NUM) {
	    fits_report_error(stderr, status);
	    fits_close_file(fptr, &status);
	    throw affineExcept("fitsData","readData",
			       "Error looking for mask",2);
	  }
      }
    }
  }

  //Read in mask if present
  status = 0;
  if (has_mask) {
    //Check checksum if available
    fits_verify_chksum(fptr, &dataok, &hduok, &status );
    if (dataok == -1) {
      fits_close_file(fptr,&status);
      throw affineExcept("fitsData","readData",
			 "Data in mask extension fails checksum test",4);
    }
    if (hduok == -1) {
      fits_close_file(fptr,&status);
      throw affineExcept("fitsData","readData",
			 "Mask extension HDU fails checksum test",8);
    }
    
    fits_get_img_dim(fptr, &naxis, &status);
    if (status) {
      fits_close_file(fptr,&status);
      throw affineExcept("fitsData","readData",
			 "Error getting mask dimensions",16);
    }
    if (naxis != 2) {
      fits_close_file(fptr,&status);
      throw affineExcept("fitsData","readData",
			 "Mask present but not 2D",32);
    }

    fits_get_img_size(fptr, 2, masknaxes, &status);
    if (status) {
      fits_close_file(fptr,&status);
      throw affineExcept("fitsData","readData",
			 "Mask failure getting mask extents",64);
    }
    
    nmask = masknaxes[0] * masknaxes[1];
    if (nmask == 0) {
      fits_close_file(fptr,&status);
      throw affineExcept("fitsData","readData",
			 "Mask present zero extent",128);
    } 

    //The actual read.  Man fitsio takes a lot of code...
    mask = (int*) fftw_malloc( sizeof(int)*nmask );
    fpixel[0] = 1; fpixel[1] = 1;
    fits_read_pix(fptr, TINT, fpixel, nmask, 0, mask, 0, &status);
    if (status) {
      fits_report_error(stderr,status);
      fits_close_file(fptr,&status);
      fftw_free(mask);
      throw affineExcept("fitsData","readData",
			 "Error reading mask data",256);
    }
  }

  //Now, read the image
  //We look for image or IMAGE.  If we can't find it, just try to
  // go to the primary and hope
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
      if (has_mask) fftw_free(mask);
      throw affineExcept("fitsData","readData",
			 "Couldn't find image",512);
    }
  }
  if (status) {
    fits_report_error(stderr, status);
    fits_close_file(fptr, &status);
    if (has_mask) fftw_free(mask);
    throw affineExcept("fitsData","readData",
		       "Error locating image data",1024);
  }

  //Check checksum if available
  fits_verify_chksum(fptr, &dataok, &hduok, &status );
  if (dataok == -1) {
    fits_close_file(fptr,&status);
    if (has_mask) fftw_free(mask);
    throw affineExcept("fitsData","readData",
		       "Image data fails checksum",2048);
  }
  if (hduok == -1) {
    fits_close_file(fptr,&status);
    if (has_mask) fftw_free(mask);
    throw affineExcept("fitsData","readData",
		       "Image HDU fails checksum",512);
  }

  fits_get_img_dim(fptr, &naxis, &status);
  if (status) {
    fits_close_file(fptr,&status);
    if (has_mask) fftw_free(mask);
    throw affineExcept("fitsData","readData",
		       "Error getting image dimensions",1024);
  }
  if (naxis != 2) {
    fits_close_file(fptr,&status);
    if (has_mask) fftw_free(mask);
    throw affineExcept("fitsData","readData",
		     "Image is not 2D",2048);
  }

  fits_get_img_size(fptr, 2, naxes, &status);
  if (status) {
    fits_close_file(fptr,&status);
    if (has_mask) fftw_free(mask);
    throw affineExcept("fitsData","readData",
		       "Error getting image dimensions",4096);
  }
  if (has_mask) {
    if (naxes[0] != masknaxes[0] ) {
      fits_close_file(fptr,&status);
      fftw_free(mask);
      throw affineExcept("fitsData","readData",
			 "Image does not match mask extent along first dimension",
			 8192);
    }
    if (naxes[1] != masknaxes[1] ) {
      fits_close_file(fptr,&status);
      fftw_free(mask);
      throw affineExcept("fitsData","readData",
			 "Image does not match mask extent along second dimension",
			 16384);
    }
  }

  ndat_ = static_cast<unsigned int>(naxes[0] * naxes[1]);
  if (ndat_ == 0) {
    fits_close_file(fptr,&status);
    if (has_mask) fftw_free(mask);
    throw affineExcept("fitsData","readData",
		       "Image is of zero extent",
		       32768);
  }

  //Read
  data_ = (double*) fftw_malloc(sizeof(double)*ndat_);
  fpixel[0] = 1; fpixel[1] = 1;
  fits_read_pix(fptr, TDOUBLE, fpixel, ndat_, 0, data_, 0, &status);
  if (status) {
    fits_report_error(stderr,status);
    fits_close_file(fptr,&status);
    if (has_mask) fftw_free(mask);
    fftw_free(data_);
    throw affineExcept("fitsData","readData","Error reading image",
		       65536);
  }

  //Close out
  fits_close_file(fptr,&status);
  if (status) {
    fits_report_error(stderr,status);
    if (has_mask) fftw_free(mask);
    fftw_free(data_);
    throw affineExcept("fitsData","readData","Error closing file",
		       131072);
  }

  //If there is no mask, easy cakes
  if ( ! has_mask ) {
    data = data_;
    n = ndat_;
  } else {
    //First, count the number of bits to copy
    unsigned int nkeep = 0;
    for (unsigned int i = 0; i < ndat_; ++i)
      if ( mask[i] == 0 ) ++nkeep;
    if (nkeep == 0) {
      if (data_ != NULL) fftw_free(data_);      
      if (mask != NULL) fftw_free(mask);
      throw affineExcept("fitsData","readData",
			 "All data masked",262144);
    }

    //Now actually copy
    data = (double*) fftw_malloc(sizeof(double)*nkeep);
    unsigned int cntr = 0;
    for (unsigned int i = 0; i < ndat_; ++i)
      if ( mask[i] == 0 ) data[++cntr] = data_[i];
    
    fftw_free(data_);
    n = nkeep;
    fftw_free(mask);
  }

  if (meansub) meanSubtract();
}

double fitsData::getMin() const {
  if (n == 0) return std::numeric_limits<double>::quiet_NaN();
  //There is a shortcut if we are binned, since bincent0
  // is the minimum value by construction
  if (is_binned) return bincent0;
  double minval, cval;
  minval = data[0];
  for (unsigned int i = 1; i < n; ++i) {
    cval = data[i];
    if (cval < minval) minval=cval;
  }
  return minval;
}

double fitsData::getMax() const {
  if (n == 0) return std::numeric_limits<double>::quiet_NaN();
  //There is a shortcut if we are binned, since 
  // bincent0 + (nbins - 1)*bindelta
  // is the maximum value by construction
  if (is_binned) return bincent0 + static_cast<double>(nbins-1)*bindelta;
  double maxval, cval;
  maxval = data[0];
  for (unsigned int i = 1; i < n; ++i) {
    cval = data[i];
    if (cval > maxval) maxval=cval;
  }
  return maxval;
}

/*!
  \param[out] min Minimum flux
  \param[out] max Maximum flux
 */
void fitsData::getMinMax(double& min, double& max) const {
  if (n == 0) {
    min = max = std::numeric_limits<double>::quiet_NaN();
    return;
  }
  
  if (is_binned) {
    min = bincent0;
    max = bincent0 + static_cast<double>(nbins-1)*bindelta;
    return;
  }

  double cval;
  min = max = data[0];
  for (unsigned int i = 1; i < n; ++i) {
    cval = data[i];
    if (cval < min) min=cval;
    if (cval > max) max=cval;
  }
}

double fitsData::getMean() const {
  if (n == 0) return std::numeric_limits<double>::quiet_NaN();
  double mnval;
  mnval = data[0];
  for (unsigned int i=1; i < n; ++i) mnval += data[i];
  mnval /= static_cast<double>(n);
  return mnval;
}

double fitsData::meanSubtract() {
  double mnval = getMean();
  for (unsigned int i = 0; i < n; ++i) data[i] -= mnval;
  
  //Have to adjust binning if we are binned
  if (is_binned) bincent0 -= mnval;

  return mnval;
}

void fitsData::applyBinning(unsigned int NBINS) {
  if (n == 0) return;
  if (NBINS == 0) throw affineExcept("fitsData","applyBinning",
				     "Trying to bin with no bins",1);
  if (is_binned && NBINS==nbins) return;
  if (NBINS != nbins) {
    if (binval != NULL) fftw_free(binval);
    binval = (unsigned int*) fftw_malloc( sizeof(unsigned int)*NBINS );
    nbins = NBINS;
  }
  for (unsigned int i = 0; i < nbins; ++i) binval[i] = 0;

  //First, we need the minimum and maximum
  double minval, maxval;
  getMinMax(minval,maxval);

  //We want to put the max and min in the center of the top and
  // bottom bin (not at the edges -- a somewhat arbitrary choice)
   bincent0 = minval;
   if (nbins == 1)
     bindelta = 2*(maxval-minval);
   else
     bindelta = (maxval-minval)/static_cast<double>(nbins-1);

   //And... bin
   double ibindelta = 1.0/bindelta;
   unsigned int idx;
   for (unsigned int i = 0; i < n; ++i) {
     idx = static_cast<unsigned int>( (data[i]-bincent0)*ibindelta + 0.5 );
     binval[idx] += 1;
   }
   is_binned = true;
}

void fitsData::removeBinning() {
  if (!is_binned) return;
  if (binval != NULL) fftw_free(binval);
  nbins = 0;
  is_binned = false;
}

void fitsData::SendSelf(MPI::Comm& comm, int dest) const {

  comm.Send(&n,1,MPI::UNSIGNED,dest,pofd_mcmc::FDSENDN);
  if (n > 0) 
    comm.Send(data,n,MPI::DOUBLE,dest,pofd_mcmc::FDSENDDATA);
  comm.Send(&is_binned,1,MPI::BOOL,dest,pofd_mcmc::FDSENDISBINNED);
  if (is_binned) {
    comm.Send(&nbins,1,MPI::UNSIGNED,dest,pofd_mcmc::FDSENDNBINS);
    if (nbins > 0) {
      comm.Send(binval,nbins,MPI::UNSIGNED,dest,pofd_mcmc::FDSENDBINVAL);
      comm.Send(&bincent0,1,MPI::DOUBLE,dest,pofd_mcmc::FDSENDBINCENT0);
      comm.Send(&bindelta,1,MPI::DOUBLE,dest,pofd_mcmc::FDSENDBINDELTA);
    }
  }
}

void fitsData::RecieveCopy(MPI::Comm& comm, int src) {
  unsigned int newn;
  comm.Recv(&newn,1,MPI::UNSIGNED,src,pofd_mcmc::FDSENDN);

  if (newn != n) {
    if (data != NULL) fftw_free(data);
    if (newn > 0)
      data = (double*) fftw_malloc( sizeof(double)*newn );
    else data = NULL;
    n = newn;
  }
  if (n > 0) 
    comm.Recv(data,n,MPI::DOUBLE,src,pofd_mcmc::FDSENDDATA);

  bool recbin;
  comm.Recv(&recbin,1,MPI::BOOL,src,pofd_mcmc::FDSENDISBINNED);
  if (recbin) {
    unsigned int newnbins;
    comm.Recv(&newnbins,1,MPI::UNSIGNED,src,pofd_mcmc::FDSENDNBINS);
    if (newnbins != nbins) {
      if (binval != NULL) fftw_free(binval);
      if (newnbins > 0) {
	binval = (unsigned int*) fftw_malloc( sizeof(unsigned int)*newnbins );
	nbins = newnbins;
	comm.Recv(binval,nbins,MPI::UNSIGNED,src,pofd_mcmc::FDSENDBINVAL);
	comm.Recv(&bincent0,1,MPI::DOUBLE,src,pofd_mcmc::FDSENDBINCENT0);
	comm.Recv(&bindelta,1,MPI::DOUBLE,src,pofd_mcmc::FDSENDBINDELTA);
	is_binned = true;
      } else is_binned = false;
    }
  } else if (is_binned) {
    if (binval != NULL) fftw_free(binval);
    binval = NULL;
    is_binned = false;
    nbins = 0;
  }

}
