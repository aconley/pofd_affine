#include<limits>
#include<cmath>
#include<map>
#include<sstream>

#include "../include/calcLikeDouble.h"
#include "../include/affineExcept.h"

const double calcLikeDoubleSingle::bad_like = 1e25;
const double calcLikeDoubleSingle::flux_safety = 1.1;
const double NaN = std::numeric_limits<double>::quiet_NaN();

/*!
  \param[in] NEDGE Number of edge values if doing edge integration
*/
calcLikeDoubleSingle::calcLikeDoubleSingle(unsigned int NEDGE) :
  data_read(false), ndatasets(0), data(NULL), minDataFlux1(NaN),
  maxDataFlux1(NaN), minDataFlux2(NaN), maxDataFlux2(NaN),
  pdfac(NEDGE), minRFlux1(NaN), maxRFlux1(NaN), minRFlux2(NaN),
  maxRFlux2(NaN), sigma_base1(NULL), sigma_base2(NULL), maxsigma_base1(NaN),
  maxsigma_base2(NaN), exp_conf1(0.0), exp_conf2(0.0), 
  like_norm(NULL), like_offset(NULL), has_beam(false), verbose(false) {}

calcLikeDoubleSingle::~calcLikeDoubleSingle() {
  if (data != NULL) delete[] data;
  if (sigma_base1 != NULL) delete[] sigma_base1;
  if (sigma_base2 != NULL) delete[] sigma_base2;
  if (like_norm   != NULL) delete[] like_norm;
  if (like_offset != NULL) delete[] like_offset;
}

void calcLikeDoubleSingle::free() {
  if (data != NULL) { delete[] data; data = NULL; }
  ndatasets = 0;
  data_read = false;
  minDataFlux1 = maxDataFlux1 = minDataFlux2 = maxDataFlux2 = NaN;
  minRFlux1 = maxRFlux1 = minRFlux2 = maxRFlux2 = NaN;

  pd.strict_resize(0, 0);
  pdfac.free();
  
  if (sigma_base1 != NULL) { delete[] sigma_base1; sigma_base1 = NULL; }
  maxsigma_base1 = NaN;
  if (sigma_base2 != NULL) { delete[] sigma_base2; sigma_base2 = NULL; }
  maxsigma_base2 = NaN;
  exp_conf1 = exp_conf2 = 0.0;

  if (like_norm != NULL) { delete[] like_norm; like_norm = NULL; }
  if (like_offset != NULL) { delete[] like_offset; like_offset = NULL; }

  has_beam = false;
  bm.free();
}

/*!
  \param[in] n New number of datasets
*/
void calcLikeDoubleSingle::resize(unsigned int n) {
  if (ndatasets == n) return;  //Don't have to do anything

  if (data != NULL)        delete[] data;
  if (sigma_base1 != NULL) delete[] sigma_base1;
  if (sigma_base2 != NULL) delete[] sigma_base2;
  if (like_offset != NULL) delete[] like_offset;
  if (like_norm   != NULL) delete[] like_norm;

  if (n > 0) {
    data        = new fitsDataDouble[n];
    sigma_base1 = new double[n];
    sigma_base2 = new double[n];
    like_offset = new double[n];
    like_norm   = new double[n];
    for (unsigned int i = 0; i < n; ++i)
      sigma_base1[i] = std::numeric_limits<double>::quiet_NaN();
    for (unsigned int i = 0; i < n; ++i)
      sigma_base2[i] = std::numeric_limits<double>::quiet_NaN();
    for (unsigned int i = 0; i < n; ++i)
      like_offset[i] = std::numeric_limits<double>::quiet_NaN();
    for (unsigned int i = 0; i < n; ++i)
      like_norm[i]   = 1.0;
  } else {
    data        = NULL;
    sigma_base1 = NULL;
    sigma_base2 = NULL;
    like_offset = NULL;
    like_norm   = NULL;
  }
  
  data_read = false;
  minDataFlux1 = maxDataFlux1 = minDataFlux2 = maxDataFlux2 = NaN;
  minRFlux1 = maxRFlux1 = minRFlux2 = maxRFlux2 = NaN;
  maxsigma_base1 = maxsigma_base2 = NaN;
  ndatasets = n;
}

/*!
  \param[in] datafile1 File to read data from, band 1
  \param[in] datafile2 File to read data from, band 2
  \param[in] IGNOREMASK Ignore mask info in files
  \param[in] MEANSUB Subtract mean from data
  \param[in] BINDATA Bin the data
  \param[in] NBINS Number of bins to use

  Special case for only a single data file
*/
void calcLikeDoubleSingle::readDataFromFile(const std::string& datafile1, 
					    const std::string& datafile2, 
					    bool IGNOREMASK, bool MEANSUB,
					    bool BINDATA, unsigned int NBINS) {
  resize(1);

  data[0].readData(datafile1, datafile2, IGNOREMASK, MEANSUB);
  unsigned int nd = data[0].getN();
  if (nd == 0)
    throw affineExcept("calcLikeDoubleSingle", "readDataFromFile",
		       "No unmasked pixels in data");
  if (BINDATA) data[0].applyBinning(NBINS, NBINS);
  
  std::pair<dblpair, dblpair> pr = data[0].getMinMax();
  if (std::isnan(pr.first.first) || std::isinf(pr.first.first) || 
      std::isnan(pr.first.second) || std::isinf(pr.first.second)) {
    std::stringstream errstr;
    errstr << "Problem with data range in file: " << datafile1;
    throw affineExcept("calcLikeDoubleSingle", "readDataFromFile", 
		       errstr.str());
  }
  if (std::isnan(pr.second.first) || std::isinf(pr.second.first) || 
      std::isnan(pr.second.second) || std::isinf(pr.second.second)) {
    std::stringstream errstr;
    errstr << "Problem with data range in file: " << datafile2;
    throw affineExcept("calcLikeDoubleSingle", "readDataFromFile", 
		       errstr.str());
  }

  minDataFlux1 = pr.first.first;
  maxDataFlux1 = pr.first.second;
  minDataFlux1 = pr.second.first;
  maxDataFlux1 = pr.second.second;

  // Can't be computed yet
  like_offset[0] = 0.0;

  data_read = true;
  minRFlux1 = maxRFlux1 = minRFlux2 = maxRFlux2 = NaN; // No longer valid
}

/*!
  \param[in] datafiles1 Files to read data from, band 1
  \param[in] datafiles2 Files to read data from, band 2
  \param[in] IGNOREMASK Ignore mask info in files
  \param[in] MEANSUB Subtract mean from data
  \param[in] BINDATA Bin the data
  \param[in] NBINS Number of bins to use
*/
void calcLikeDoubleSingle::
readDataFromFiles(const std::vector<std::string>& datafiles1, 
		  const std::vector<std::string>& datafiles2, 
		  bool IGNOREMASK, bool MEANSUB, bool BINDATA,
		  unsigned int NBINS) {
  unsigned int n = datafiles1.size();
  if (n == 0)
    throw affineExcept("calcLikeDoubleSingle", "readDataFromFiles",
		       "No data sets");
  if (datafiles2.size() != n)
    throw affineExcept("calcLikeDoubleSingle", "readDataFromFiles",
		       "Different number of datafiles in each band");
  resize(n);

  for (unsigned int i = 0; i < n; ++i) {
    data[i].readData(datafiles1[i], datafiles2[i], IGNOREMASK, MEANSUB);
    unsigned int nd = data[i].getN();
    if (nd == 0) {
      std::stringstream errstr("");
      errstr << "No unmasked pixels in data files: "
	     << datafiles1[i] << " " << datafiles2[i];
      throw affineExcept("calcLikeDoubleSingle", "readDataFromFiles",
			 errstr.str());
    }
    if (BINDATA) data[i].applyBinning(NBINS, NBINS);

    // Can't be computed yet
    like_offset[i] = 0.0;
  }

  // Determine maximum flux for all this data
  // Start with first file
  std::pair<dblpair, dblpair> pr = data[0].getMinMax();
  if (std::isnan(pr.first.first) || std::isinf(pr.first.first) || 
      std::isnan(pr.first.second) || std::isinf(pr.first.second)) {
    std::stringstream errstr;
    errstr << "Problem with data range in file: " << datafiles1[0];
    throw affineExcept("calcLikeDoubleSingle", "readDataFromFile", 
		       errstr.str());
  }
  if (std::isnan(pr.second.first) || std::isinf(pr.second.first) || 
      std::isnan(pr.second.second) || std::isinf(pr.second.second)) {
    std::stringstream errstr;
    errstr << "Problem with data range in file: " << datafiles2[0];
    throw affineExcept("calcLikeDoubleSingle", "readDataFromFile", 
		       errstr.str());
  }
  minDataFlux1 = pr.first.first;
  maxDataFlux1 = pr.first.second;
  minDataFlux1 = pr.second.first;
  maxDataFlux1 = pr.second.second;
  // Now compare to all the others.
  for (unsigned int i = 1; i < n; ++i) {
    pr = data[i].getMinMax();
    if (std::isnan(pr.first.first) || std::isinf(pr.first.first) || 
	std::isnan(pr.first.second) || std::isinf(pr.first.second)) {
      std::stringstream errstr;
      errstr << "Problem with data range in file: " << datafiles1[i];
      throw affineExcept("calcLikeDoubleSingle", "readDataFromFile", 
			 errstr.str());
    }
    if (std::isnan(pr.second.first) || std::isinf(pr.second.first) || 
	std::isnan(pr.second.second) || std::isinf(pr.second.second)) {
      std::stringstream errstr;
      errstr << "Problem with data range in file: " << datafiles2[i];
      throw affineExcept("calcLikeDoubleSingle", "readDataFromFile", 
			 errstr.str());
    }
    if (pr.first.first < minDataFlux1) minDataFlux1 = pr.first.first;
    if (pr.first.second > maxDataFlux1) maxDataFlux1 = pr.first.second;
    if (pr.second.first < minDataFlux2) minDataFlux2 = pr.second.first;
    if (pr.second.second > maxDataFlux2) maxDataFlux2 = pr.second.second;
  }

  data_read = true;
  minRFlux1 = maxRFlux1 = minRFlux2 = maxRFlux2 = NaN; // No longer valid
}

/*!
  \param[in] fl1 FITS filename to read band 1 beam from
  \param[in] fl2 FITS filename to read band 2 beam from
  \param[in] MINVAL Minimum value; specificall, only values where
     both beams satisfy fabs(val) > minval are kept.
  \param[in] histogram Apply beam histogramming
  \param[in] NBINS Number of beam histogram bins to use
*/
void calcLikeDoubleSingle::readBeam(const std::string& fl1, 
				    const std::string& fl2,
				    double MINVAL,
				    bool histogram, 
				    unsigned int NBINS) {
  bm.readFiles(fl1, fl2, MINVAL);
  if (histogram) bm.makeHistogram(NBINS);
  has_beam = true;
}

/*!
  \param[in] nbins Number of bins
*/
void calcLikeDoubleSingle::applyBinning(unsigned int nbins) {
  if ((!data_read) || ndatasets == 0) return;
  //Does nothing if the data is already binned at the same size
  for (unsigned int i = 0; i < ndatasets; ++i)
    data[i].applyBinning(nbins, nbins);
}

void calcLikeDoubleSingle::removeBinning() {
  if ((!data_read) || ndatasets == 0) return;
  for (unsigned int i = 0; i < ndatasets; ++i)
    data[i].removeBinning();
}

/*!
  \param[in] lnorm Vector of likelihood normalization.  Must have
    ndatasets number of elements.
*/
void calcLikeDoubleSingle::setLikeNorm(const std::vector<double>& lnorm) {
  unsigned int n = lnorm.size();
  if (n != ndatasets)
    throw affineExcept("calcLikeDoubleSingle", "setLikeNorm",
		       "like_norm vector not same size as number of data sets");
  if (n == 0) return;
  for (unsigned int i = 0; i < ndatasets; ++i)
    like_norm[i] = lnorm[i];
}

/*!
  \param[in] n Length of lnorm.  Must be equal to the number of datasets.
  \param[in] lnorm Array of likelihood normalization. 
*/
void calcLikeDoubleSingle::setLikeNorm(unsigned int n, 
				       const double* const lnorm) {
  if (n != ndatasets)
    throw affineExcept("calcLikeDoubleSingle", "setLikeNorm",
		       "like_norm array not same size as number of data sets");
  if (n == 0) return;
  for (unsigned int i = 0; i < ndatasets; ++i)
    like_norm[i] = lnorm[i];
}

/*!
  \param[in] s Vector of instrument noise base values, band 1.  Must have
    ndatasets number of elements.
*/
void calcLikeDoubleSingle::setSigmaBase1(const std::vector<double>& s) {
  unsigned int n = s.size();
  if (n != ndatasets)
    throw affineExcept("calcLikeDoubleSingle", "setSigmaBase1",
		       "sigma vectors not same size as number of data sets");
  if (n == 0) return;
  for (unsigned int i = 0; i < ndatasets; ++i)
    sigma_base1[i] = s[i];
  maxsigma_base1 = sigma_base1[0];
  for (unsigned int i = 1; i < ndatasets; ++i)
    if (sigma_base1[i] > maxsigma_base1) maxsigma_base1 = sigma_base1[i];

  if (data_read) {
    double var1, var2, sigprod;
    for (unsigned int i = 0; i < ndatasets; ++i) {
      var1 = sigma_base1[i] * sigma_base1[i] + exp_conf1 * exp_conf1;
      var2 = sigma_base2[i] * sigma_base2[i] + exp_conf2 * exp_conf2;
      sigprod = sqrt(var1 * var2);
      if (sigprod > 0) {
	double ndd = static_cast<double>(data[i].getN());
	like_offset[i] = -ndd * (log(mcmc_affine::two_pi * sigprod) + 1.0);
      }
    }
  }
}

/*!
  \param[in] s Vector of instrument noise base values, band 2.  Must have
    ndatasets number of elements.
*/
void calcLikeDoubleSingle::setSigmaBase2(const std::vector<double>& s) {
  unsigned int n = s.size();
  if (n != ndatasets)
    throw affineExcept("calcLikeDoubleSingle", "setSigmaBase2",
		       "sigma vectors not same size as number of data sets");
  if (n == 0) return;
  for (unsigned int i = 0; i < ndatasets; ++i)
    sigma_base2[i] = s[i];
  maxsigma_base2 = sigma_base2[0];
  for (unsigned int i = 1; i < ndatasets; ++i)
    if (sigma_base2[i] > maxsigma_base2) maxsigma_base2 = sigma_base2[i];

  if (data_read) {
    double var1, var2, sigprod;
    for (unsigned int i = 0; i < ndatasets; ++i) {
      var1 = sigma_base1[i] * sigma_base1[i] + exp_conf1 * exp_conf1;
      var2 = sigma_base2[i] * sigma_base2[i] + exp_conf2 * exp_conf2;
      sigprod = sqrt(var1 * var2);
      if (sigprod > 0) {
	double ndd = static_cast<double>(data[i].getN());
	like_offset[i] = -ndd * (log(mcmc_affine::two_pi * sigprod) + 1.0);
      }
    }
  }
}

/*!
  \param[in] n Length of s.  Must be same as the number of datasets
  \param[in] s Array of instrument noise base values, band 1
*/
void calcLikeDoubleSingle::setSigmaBase1(unsigned int n,const double* const s) {
  if (n != ndatasets)
    throw affineExcept("calcLikeDoubleSingle", "setSigmaBase1",
		       "sigma arrays not same size as number of data sets");
  if (n == 0) return;
  for (unsigned int i = 0; i < ndatasets; ++i)
    sigma_base1[i] = s[i];
  for (unsigned int i = 1; i < ndatasets; ++i)
    if (sigma_base1[i] > maxsigma_base1) maxsigma_base1 = sigma_base1[i];

  if (data_read) {
    double var1, var2, sigprod;
    for (unsigned int i = 0; i < ndatasets; ++i) {
      var1 = sigma_base1[i] * sigma_base1[i] + exp_conf1 * exp_conf1;
      var2 = sigma_base2[i] * sigma_base2[i] + exp_conf2 * exp_conf2;
      sigprod = sqrt(var1 * var2);
      if (sigprod > 0) {
	double ndd = static_cast<double>(data[i].getN());
	like_offset[i] = -ndd * (log(mcmc_affine::two_pi * sigprod) + 1.0);
      }
    }
  }
}

/*!
  \param[in] n Length of s.  Must be same as the number of datasets
  \param[in] s Array of instrument noise base values, band 2
*/
void calcLikeDoubleSingle::setSigmaBase2(unsigned int n,const double* const s) {
  if (n != ndatasets)
    throw affineExcept("calcLikeDoubleSingle", "setSigmaBase2",
		       "sigma arrays not same size as number of data sets");
  if (n == 0) return;
  for (unsigned int i = 0; i < ndatasets; ++i)
    sigma_base2[i] = s[i];
  for (unsigned int i = 2; i < ndatasets; ++i)
    if (sigma_base2[i] > maxsigma_base2) maxsigma_base2 = sigma_base2[i];

  if (data_read) {
    double var1, var2, sigprod;
    for (unsigned int i = 0; i < ndatasets; ++i) {
      var1 = sigma_base1[i] * sigma_base1[i] + exp_conf1 * exp_conf1;
      var2 = sigma_base2[i] * sigma_base2[i] + exp_conf2 * exp_conf2;
      sigprod = sqrt(var1 * var2);
      if (sigprod > 0) {
	double ndd = static_cast<double>(data[i].getN());
	like_offset[i] = -ndd * (log(mcmc_affine::two_pi * sigprod) + 1.0);
      }
    }
  }
}

/*!
  \param[in] model Model to set min/max R range from 

  You must also have loaded in data, set maxsigma_base, read the beam, etc.
  This should not be called during a fit, only before the fit, because
  moving the R ranges around causes numerical jitter in the likelihoods.
*/
void calcLikeDoubleSingle::setRRange(const numberCountsDouble& model) { 
  if (std::isnan(maxsigma_base1) || std::isnan(maxsigma_base2))
    throw affineExcept("calcLikeDoubleSingle", "setRRange", 
		       "Sigma base1 not yet set");
  if (!has_beam)
    throw affineExcept("calcLikeDoubleSingle", "setRRange", "Beam not loaded");
  if (!data_read)
    throw affineExcept("calcLikeDoubleSingle", "setRRange", "Data not loaded");

  // First get the raw range where R is non-zero
  // model validity, beam validity checked in getRRange
  std::pair<dblpair, dblpair> rawrange = model.getRRange(bm);

  // Now try to estimate sigma.  Ideally the user has helped us out here
  // by setting exp_conf to something non-zero
  double sigma1 = sqrt(exp_conf1 * exp_conf1 + maxsigma_base1 * maxsigma_base1);
  double sigma2 = sqrt(exp_conf2 * exp_conf2 + maxsigma_base2 * maxsigma_base2);
  
  // Add some padding to the top of the r range
  minRFlux1 = rawrange.first.first;
  maxRFlux1 = rawrange.first.second + pofd_mcmc::n_zero_pad * sigma1;
  minRFlux2 = rawrange.second.first;
  maxRFlux2 = rawrange.second.second + pofd_mcmc::n_zero_pad * sigma2;
  
  // Make sure this actually covers the data
  if (minDataFlux1 < minRFlux1 && (bm.hasSign(2) || bm.hasSign(3))) 
    minRFlux1 = minDataFlux1;
  if (maxDataFlux1 > maxRFlux1 && (bm.hasSign(0) || bm.hasSign(1)))
    maxRFlux1 = maxDataFlux1;
  if (minDataFlux2 < minRFlux2 && (bm.hasSign(1) || bm.hasSign(3))) 
    minRFlux2 = minDataFlux2;
  if (maxDataFlux2 > maxRFlux2 && (bm.hasSign(0) || bm.hasSign(2)))
    maxRFlux2 = maxDataFlux2;
}


/*
  \param[in] model   Model parameters must already be set
  \param[out] pars_invalid Set to true if parameters are determined to 
                            be invalid
  \param[in] sigmult1 Sigma multiplier in band 1
  \param[in] sigmult2 Sigma multiplier in band 2
  \param[in] fftsize Size of FFT to use
  \param[in] setedge Fill in edges of R with integral average

  \returns Log likelihood of data relative to model

  Any auxilliary parameters (positions of knots, etc.) must already
  be set in model.  The last model parameter is the sigma multiplier,
  the previous ones are the values of the number counts at the knots
*/
double 
calcLikeDoubleSingle::getLogLike(const numberCountsDouble& model,
				 bool& pars_invalid,
				 double sigmult1, double sigmult2,
				 unsigned int fftsize, bool setedge) const {

  if (!data_read)
    throw affineExcept("calcLikeDoubleSingle", "getNegLogLike",
		       "Data not read");
  if (!has_beam)
    throw affineExcept("calcLikeDoubleSingle", "getNegLogLike",
		       "Beam not read");

  //Maximum noise value with multiplier 
  double max_sigma1 = maxsigma_base1 * sigmult1;
  if (max_sigma1 < 0.0) return calcLikeDoubleSingle::bad_like;
  double max_sigma2 = maxsigma_base2 * sigmult2;
  if (max_sigma2 < 0.0) return calcLikeDoubleSingle::bad_like;

  // Initialize P(D)
  // Use the pre-set min/max R ranges
  if (std::isnan(minRFlux1) || std::isnan(maxRFlux1) ||
      std::isnan(minRFlux2) || std::isnan(maxRFlux2))
    throw affineExcept("calcLikeDoubleSingle", "getNegLogLike",
		       "R range not set");
  pars_invalid = ! pdfac.initPD(fftsize, minRFlux1, maxRFlux1, 
				minRFlux2, maxRFlux2, model, 
				bm, setedge);

  if (pars_invalid) return 0.0; // No point in continuing

  //Compute likelihood of each bit of data
  double curr_LogLike, LogLike;
  LogLike = 0.0;
  for (unsigned int i = 0; i < ndatasets; ++i) {
    // Get PD for this particuar set of sigmas
    pdfac.getPD(sigmult1 * sigma_base1[i], sigmult2 * sigma_base2[i],
		pd, true);

    // Get log like
    curr_LogLike = pd.getLogLike(data[i]);

    // Apply beam norm and zero offset factor
    LogLike += (curr_LogLike - like_offset[i]) / like_norm[i];
  }
  
  return LogLike;
}

/*!
  \param[inout] os Stream to write most recent P(D) to
*/
void calcLikeDoubleSingle::writePDToStream(std::ostream& os) const {
  PDDouble cpy(pd);
  cpy.deLog();
  os << cpy;
}

/*!
  \param[in] comm MPI communicator
  \param[in] dest Destination of messages
*/
void calcLikeDoubleSingle::sendSelf(MPI_Comm comm, int dest) const {
  //Data
  MPI_Send(const_cast<unsigned int*>(&ndatasets), 1, 
	   MPI_UNSIGNED, dest, pofd_mcmc::CLDSENDNDATA, comm);
  
  MPI_Send(const_cast<bool*>(&data_read), 1, MPI::BOOL, dest,
	   pofd_mcmc::CLDSENDDATAREAD, comm);
  if (data_read) {
    for (unsigned int i = 0; i < ndatasets; ++i)
      data[i].sendSelf(comm, dest);
    MPI_Send(sigma_base1, ndatasets, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDSIGMABASE1, comm);
    MPI_Send(sigma_base2, ndatasets, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDSIGMABASE2, comm);
    MPI_Send(const_cast<double*>(&maxsigma_base1), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDMAXSIGMABASE1, comm);
    MPI_Send(const_cast<double*>(&maxsigma_base2), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDMAXSIGMABASE2, comm);
    MPI_Send(const_cast<double*>(&exp_conf1), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDEXPCONF1, comm);
    MPI_Send(const_cast<double*>(&exp_conf2), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDEXPCONF2, comm);
    MPI_Send(like_offset, ndatasets, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDLIKEOFFSET, comm);
    MPI_Send(like_norm, ndatasets, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDLIKENORM, comm);
    MPI_Send(const_cast<double*>(&minDataFlux1), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDMAXDATAFLUX1, comm);
    MPI_Send(const_cast<double*>(&minDataFlux1), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDMAXDATAFLUX1, comm);
    MPI_Send(const_cast<double*>(&minDataFlux2), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDMINDATAFLUX2, comm);
    MPI_Send(const_cast<double*>(&maxDataFlux2), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDMAXDATAFLUX2, comm);
    MPI_Send(const_cast<double*>(&minRFlux1), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDMAXRFLUX1, comm);
    MPI_Send(const_cast<double*>(&minRFlux1), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDMAXRFLUX1, comm);
    MPI_Send(const_cast<double*>(&minRFlux2), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDMINRFLUX2, comm);
    MPI_Send(const_cast<double*>(&maxRFlux2), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDMAXRFLUX2, comm);
  }

  //Beam
  MPI_Send(const_cast<bool*>(&has_beam), 1, MPI::BOOL, dest, 
	   pofd_mcmc::CLDSENDHASBEAM, comm);
  if (has_beam) bm.sendSelf(comm, dest);

  //PDFactory
  pdfac.sendSelf(comm, dest);

  //Note we don't send verbose
}

/*!
  \param[in] comm MPI communicator
  \param[in] src Source of messages
*/
void calcLikeDoubleSingle::recieveCopy(MPI_Comm comm, int src) {
  MPI_Status Info;
  //Data
  unsigned int newn;
  MPI_Recv(&newn, 1, MPI_UNSIGNED, src, pofd_mcmc::CLDSENDNDATA,
	   comm, &Info);
  resize(newn);

  bool newread = false;
  MPI_Recv(&newread, 1, MPI::BOOL, src, pofd_mcmc::CLDSENDDATAREAD,
	   comm, &Info);
  if (newread) {
    for (unsigned int i = 0; i < newn; ++i)
      data[i].recieveCopy(comm, src);
    MPI_Recv(sigma_base1, newn, MPI_DOUBLE, src, 
	     pofd_mcmc::CLDSENDSIGMABASE1, comm, &Info);
    MPI_Recv(sigma_base2, newn, MPI_DOUBLE, src, pofd_mcmc::CLDSENDSIGMABASE2,
	     comm, &Info);
    MPI_Recv(&maxsigma_base1, 1, MPI_DOUBLE, src, 
	     pofd_mcmc::CLDSENDMAXSIGMABASE1, comm, &Info);
    MPI_Recv(&maxsigma_base2, 1, MPI_DOUBLE, src,
	     pofd_mcmc::CLDSENDMAXSIGMABASE2, comm, &Info);
    MPI_Recv(&exp_conf1, 1, MPI_DOUBLE, src, pofd_mcmc::CLDSENDEXPCONF1,
	     comm, &Info);
    MPI_Recv(&exp_conf2, 1, MPI_DOUBLE, src, pofd_mcmc::CLDSENDEXPCONF2,
	     comm, &Info);
    MPI_Recv(like_offset, newn, MPI_DOUBLE, src,
	     pofd_mcmc::CLDSENDLIKEOFFSET, comm, &Info);
    MPI_Recv(like_norm, newn, MPI_DOUBLE, src, pofd_mcmc::CLDSENDLIKENORM,
	     comm, &Info);
    MPI_Recv(&minDataFlux1, 1, MPI_DOUBLE, src, pofd_mcmc::CLDSENDMINDATAFLUX1,
	     comm, &Info);
    MPI_Recv(&maxDataFlux1, 1, MPI_DOUBLE, src, pofd_mcmc::CLDSENDMAXDATAFLUX1,
	     comm, &Info);
    MPI_Recv(&minDataFlux2, 1, MPI_DOUBLE, src, pofd_mcmc::CLDSENDMINDATAFLUX2,
	     comm, &Info);
    MPI_Recv(&maxDataFlux2, 1, MPI_DOUBLE, src, pofd_mcmc::CLDSENDMAXDATAFLUX2,
	     comm, &Info);
    MPI_Recv(&minRFlux1, 1, MPI_DOUBLE, src, pofd_mcmc::CLDSENDMINRFLUX1,
	     comm, &Info);
    MPI_Recv(&maxRFlux1, 1, MPI_DOUBLE, src, pofd_mcmc::CLDSENDMAXRFLUX1,
	     comm, &Info);
    MPI_Recv(&minRFlux2, 1, MPI_DOUBLE, src, pofd_mcmc::CLDSENDMINRFLUX2,
	     comm, &Info);
    MPI_Recv(&maxRFlux2, 1, MPI_DOUBLE, src, pofd_mcmc::CLDSENDMAXRFLUX2,
	     comm, &Info);
    data_read = true;
  } else {
    minDataFlux1 = maxDataFlux1 = minDataFlux2 = maxDataFlux2 = NaN;
    minRFlux1 = maxRFlux1 = minRFlux2 = maxRFlux2 = NaN;
    maxsigma_base1 = maxsigma_base2 = NaN;
  }

  //Beam
  bool hsbm;
  MPI_Recv(&hsbm, 1, MPI::BOOL, src, pofd_mcmc::CLDSENDHASBEAM, comm, &Info);
  if (hsbm) {
    bm.recieveCopy(comm, src);
    has_beam = true;
  } else has_beam = false;
  
  //PDFactory
  pdfac.recieveCopy(comm, src);

}


/////////////////////////////////////////////////////////////////

/*!
  \param[in] FFTSIZE Size of FFT transformation (along each edge)
  \param[in] NEDGE Number of elements for edge integration (if used)
  \param[in] EDGEINTEG Use integration to set edge values
  \param[in] BINNED Apply binning to data sets
  \param[in] NBINS Number of bins if binning data set
*/
calcLikeDouble::calcLikeDouble(unsigned int FFTSIZE, unsigned int NEDGE, 
			       bool EDGEINTEG, bool BINNED, 
			       unsigned int NBINS) :
  fftsize(FFTSIZE), nedge(NEDGE), edgeInteg(EDGEINTEG),
  bin_data(BINNED), nbins(NBINS), 
  has_cfirb_prior1(false), cfirb_prior_mean1(0.0), cfirb_prior_sigma1(0.0), 
  has_cfirb_prior2(false), cfirb_prior_mean2(0.0), cfirb_prior_sigma2(0.0), 
  has_sigma_prior1(false), sigma_prior_width1(0.0), has_sigma_prior2(false), 
  sigma_prior_width2(0.0), verbose(false) {

  nbeamsets = 0;
  beamsets  = NULL;
}

calcLikeDouble::~calcLikeDouble() {
  if (beamsets != NULL) delete[] beamsets;
}

/*!
  This doesn't actually deallocate beamsets, just frees the internal
  data storage
*/
void calcLikeDouble::freeData() {
  for (unsigned int i = 0; i < nbeamsets; ++i)
    beamsets[i].free();
}


/*!
  \param[in] filename Name of FFTW wisdom file

  Don't do this before reading in the data files, or it will be overwritten
*/
void calcLikeDouble::addWisdom(const std::string& filename) {
  for (unsigned int i = 0; i < nbeamsets; ++i)
    beamsets[i].addWisdom(filename);
}

/*!
  \param[in] datafiles1 Files to read data from, band 1
  \param[in] datafiles2 Files to read data from, band 2
  \param[in] beamfiles1 Files to read beams from, band 1
  \param[in] beamfiles2 Files to read beams from, band 2
  \param[in] sigmas1 Instrumental sigma, band 1
  \param[in] sigmas2 Instrumental sigma, band 2
  \param[in] like_norms Vector of likelihood normalizations
  \param[in] IGNOREMASK Ignore any mask in data files
  \param[in] MEANSUB Mean subtract the data
  \param[in] MINBEAMVAL Minimum beam value used
  \param[in] HISTOGRAMBEAMS Histogram the beams
  \param[in] NBEAMHIST Number of beam histogram bins to use
  \param[in] EXPCONF1 Expected confusion noise value, band 1
  \param[in] EXPCONF2 Expected confusion noise value, band 2

  Read in a set of data, grouping the inputs by beam and storing
  all of the relevant instrument noise and likelihood normalization values
*/
void calcLikeDouble::
readDataFromFiles(const std::vector<std::string>& datafiles1, 
		  const std::vector<std::string>& datafiles2, 
		  const std::vector<std::string>& beamfiles1,
		  const std::vector<std::string>& beamfiles2,
		  const std::vector<double>& sigmas1,
		  const std::vector<double>& sigmas2,
		  const std::vector<double>& like_norms,
		  bool IGNOREMASK, bool MEANSUB, double MINBEAMVAL, 
		  bool HISTOGRAM, unsigned int NBEAMHIST, 
		  double EXPCONF1, double EXPCONF2) {

  //Make sure they are all the same length
  unsigned int ndat = datafiles1.size();
  if (ndat != datafiles2.size())
    throw affineExcept("calcLikeDouble", "readDataFromFiles",
		       "Datafiles1 and datafile2 not same length");
  if (ndat == 0)
    throw affineExcept("calcLikeDouble", "readDataFromFiles",
		       "No datafiles");
  if (ndat != beamfiles1.size())
    throw affineExcept("calcLikeDouble", "readDataFromFiles",
		       "datafiles and beamfiles1 not same length");
  if (ndat != beamfiles2.size())
    throw affineExcept("calcLikeDouble", "readDataFromFiles",
		       "datafiles and beamfiles2 not same length");
  if (ndat != sigmas1.size())
    throw affineExcept("calcLikeDouble", "readDataFromFiles",
		       "datafiles and sigmas1 not same length");
  if (ndat != sigmas2.size())
    throw affineExcept("calcLikeDouble", "readDataFromFiles",
		       "datafiles and sigmas2 not same length");
  if (ndat != like_norms.size())
    throw affineExcept("calcLikeDouble", "readDataFromFiles",
		       "datafiles and like_norm not same length");
  
  //Use a map.  Could also use a multi-map, but meh
  //Key will be the combination of both beam file names
  std::map< std::string, doublebeam_group > grpmap;
  std::map< std::string, doublebeam_group >::iterator grpmap_it;
  std::string key;
  for (unsigned int i = 0; i < ndat; ++i) {
    key = beamfiles1[i]+":"+beamfiles2[i];
    //See if map already has key
    grpmap_it = grpmap.find(key);
    if (grpmap_it == grpmap.end()) {
      //Previously unknown key
      doublebeam_group newgrp;
      newgrp.n = 1;
      newgrp.datafiles1.push_back(datafiles1[i]);
      newgrp.datafiles2.push_back(datafiles2[i]);
      newgrp.beamfile1 = beamfiles1[i];
      newgrp.beamfile2 = beamfiles2[i];
      newgrp.sigmas1.push_back(sigmas1[i]);
      newgrp.sigmas2.push_back(sigmas2[i]);
      newgrp.like_norms.push_back(like_norms[i]);
      grpmap[key] = newgrp;
    } else {
      //Previously known key
      grpmap_it->second.n += 1;
      grpmap_it->second.datafiles1.push_back(datafiles1[i]);
      grpmap_it->second.datafiles2.push_back(datafiles2[i]);
      grpmap_it->second.sigmas1.push_back(sigmas1[i]);
      grpmap_it->second.sigmas2.push_back(sigmas2[i]);
      grpmap_it->second.like_norms.push_back(like_norms[i]);
    }
  }

  unsigned int newnbeamsets = grpmap.size();
  if (newnbeamsets != nbeamsets) {
    if (beamsets != NULL) delete[] beamsets;
    if (newnbeamsets > 0) beamsets = new calcLikeDoubleSingle[newnbeamsets];
    else beamsets = NULL;
    nbeamsets = newnbeamsets;
  }

  if (nbeamsets > 0) {
    grpmap_it = grpmap.begin();
    for (unsigned int i=0; grpmap_it != grpmap.end(); ++grpmap_it, ++i) {
      beamsets[i].setNEdge(nedge);
      beamsets[i].readDataFromFiles(grpmap_it->second.datafiles1,
				    grpmap_it->second.datafiles2,
				    IGNOREMASK, MEANSUB, bin_data, nbins);
      beamsets[i].readBeam(grpmap_it->second.beamfile1, 
			   grpmap_it->second.beamfile2, 
			   MINBEAMVAL, HISTOGRAM, NBEAMHIST);
      beamsets[i].setExpConf1(EXPCONF1);
      beamsets[i].setExpConf2(EXPCONF2);
      beamsets[i].setSigmaBase1(grpmap_it->second.sigmas1);
      beamsets[i].setSigmaBase2(grpmap_it->second.sigmas2);
      beamsets[i].setLikeNorm(grpmap_it->second.like_norms);
    }
  }
}

/*! 
  \param[in] nedg New edge integration size
 */
void calcLikeDouble::setNEdge(unsigned int nedg) {
  if (nedg == nedge) return;
  if (beamsets != NULL)
    for (unsigned int i = 0; i < nbeamsets; ++i)
      beamsets[i].setNEdge(nedg);
  nedge = nedg;
}

void calcLikeDouble::setBinData() {
  if (bin_data) return;

  //Easy if no data is read.
  if (nbeamsets == 0) {
    bin_data = true;
    return;
  }

  //Now we have to actually bin
  for (unsigned int i = 0; i < nbeamsets; ++i)
    beamsets[i].applyBinning(nbins);
  bin_data = true;

}

void calcLikeDouble::unSetBinData() {
  if (!bin_data) return;

  if (nbeamsets == 0) {
    bin_data = false;
    return;
  }

  for (unsigned int i = 0; i < nbeamsets; ++i)
    beamsets[i].removeBinning();
  bin_data = false;
}

/*!
  \param[in] nbns New number of bins for data binning
*/
void calcLikeDouble::setNBins(unsigned int nbns) {
  if (nbns == nbins) return;
  
  if (nbeamsets == 0) {
    nbins = nbns;
    return;
  }

  if (bin_data) //rebin
    for (unsigned int i = 0; i < nbeamsets; ++i)
      beamsets[i].applyBinning(nbns);

  nbins = nbns;
}

/*!
  \param[in] p Model parameters to use
*/
void calcLikeDouble::setRRanges(const paramSet& p) {
  model.setParams(p);
  for (unsigned int i = 0; i < nbeamsets; ++i)
    beamsets[i].setRRange(model);
}

/*!
  \param[in] mn Mean value of CFIRB prior, band 1
  \param[in] sg Sigma of CFIRB prior, band 1
  
  The prior is assumed Gaussian
*/
void calcLikeDouble::setCFIRBPrior1(double mn, double sg) {
  has_cfirb_prior1   = true;
  cfirb_prior_mean1  = mn;
  cfirb_prior_sigma1 = sg;
}

/*!
  \param[in] mn Mean value of CFIRB prior, band 2
  \param[in] sg Sigma of CFIRB prior, band 2
  
  The prior is assumed Gaussian
*/
void calcLikeDouble::setCFIRBPrior2(double mn, double sg) {
  has_cfirb_prior2   = true;
  cfirb_prior_mean2  = mn;
  cfirb_prior_sigma2 = sg;
}

/*!
  \param[in] ifile Object holding information from model initializaiton file
*/
void calcLikeDouble::setPositions(const initFileDoubleLogNormal& ifile) {
  ifile.getModelPositions(model);
}

/*!
  \param[in] p Parameters to evaluate likelihood for
  \param[out] pars_invalid True if parameters were not valid, False otherwise
  \returns Log likelihood of data relative to model, including priors
*/
double calcLikeDouble::getLogLike(const paramSet& p, bool& pars_invalid) const {
  const double half_log_2pi = 0.918938517570495605469;

  if (nbeamsets == 0) return std::numeric_limits<double>::quiet_NaN();
  unsigned int npar = p.getNParams();
  unsigned int nmodelpars = model.getNParams();
  if (npar < (nmodelpars+2)) //+2 for the sigma multiplier in each band
    throw affineExcept("calcLikeDouble", "getLogLike",
		       "Not enough elements in params");

  //Set the model
  model.setParams(p);

  double LogLike = 0.0;

  //Do the datasets likelihood
  bool pinvalid;
  double sigmult1 = p[nmodelpars]; //Sigma multiplier is after knot values
  double sigmult2 = p[nmodelpars+1];
  pars_invalid = false;
  for (unsigned int i = 0; i < nbeamsets; ++i) {
    LogLike += beamsets[i].getLogLike(model, pinvalid, sigmult1, sigmult2, 
				      fftsize, edgeInteg);
    pars_invalid &= pinvalid;
  }

  if (pars_invalid) return LogLike; //!< Not much point in doing the priors...

  //Add on cfirb prior and sigma prior if needed
  //Only do this once for all data sets so as not to multi-count the prior
  if (has_sigma_prior1) {
    //Assume the mean is always at 1 -- otherwise, the
    // user would have specified different noise level
    double val = (sigmult1-1.0) / sigma_prior_width1;
    LogLike -= half_log_2pi + log(sigma_prior_width1) + 0.5*val*val;
  }
  if (has_sigma_prior2) {
    double val = (sigmult2-1.0) / sigma_prior_width2;
    LogLike -= half_log_2pi + log(sigma_prior_width2) + 0.5*val*val;
  }

  if (has_cfirb_prior1) {
    double s_per_area = model.getFluxPerArea(0);
    double val = (cfirb_prior_mean1 - s_per_area) / cfirb_prior_sigma1;
    LogLike -=  half_log_2pi + log(cfirb_prior_sigma1) + 0.5*val*val;
  }
  if (has_cfirb_prior2) {
    double s_per_area = model.getFluxPerArea(1);
    double val = (cfirb_prior_mean2 - s_per_area) / cfirb_prior_sigma2;
    LogLike -=  half_log_2pi + log(cfirb_prior_sigma2) + 0.5*val*val;
  }

  return LogLike;
}

/*!						
  \param[in] objid HDF5 handle to write information to.  Must already be
    open
*/
void calcLikeDouble::writeToHDF5Handle(hid_t objid) const {
  // Writes some meta information as a sub-group
  hsize_t adims;
  hid_t mems_id, att_id;

  if (H5Iget_ref(objid) < 0)
    throw affineExcept("calcLikeDouble", "writeToHDF5",
		       "Input handle is not valid");

  // FFTSIZE
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate2(objid, "fftsize", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &fftsize);
  H5Aclose(att_id);

  // NEDGE
  att_id = H5Acreate2(objid, "nedge", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &nedge);
  H5Aclose(att_id);

  // Edge integrate
  hbool_t bl;
  att_id = H5Acreate2(objid, "edge_integrate", H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  bl = static_cast<hbool_t>(edgeInteg);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, &bl);
  H5Aclose(att_id);

  // NBEAMSETS
  att_id = H5Acreate2(objid, "nbeamsets", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &nbeamsets);
  H5Aclose(att_id);

  // DATA BINNING
  att_id = H5Acreate2(objid, "bin_data", H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  bl = static_cast<hbool_t>(bin_data);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, &bl);
  H5Aclose(att_id);
  if (bin_data) {
    att_id = H5Acreate2(objid, "ndatabins", H5T_NATIVE_UINT,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_UINT, &nbins);
    H5Aclose(att_id);
  }

  // CFIRB PRIOR1
  att_id = H5Acreate2(objid, "has_cfirb_prior1", H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  bl = static_cast<hbool_t>(has_cfirb_prior1);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, &bl);
  H5Aclose(att_id);
  if (has_cfirb_prior1) {
    att_id = H5Acreate2(objid, "cfirb_prior_mean1", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_DOUBLE, &cfirb_prior_mean1);
    H5Aclose(att_id);
    att_id = H5Acreate2(objid, "cfirb_prior_sigma1", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_DOUBLE, &cfirb_prior_sigma1);
    H5Aclose(att_id);
  }

  // CFIRB PRIOR2
  att_id = H5Acreate2(objid, "has_cfirb_prior2", H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  bl = static_cast<hbool_t>(has_cfirb_prior2);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, &bl);
  H5Aclose(att_id);
  if (has_cfirb_prior2) {
    att_id = H5Acreate2(objid, "cfirb_prior_mean2", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_DOUBLE, &cfirb_prior_mean1);
    H5Aclose(att_id);
    att_id = H5Acreate2(objid, "cfirb_prior_sigma2", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_DOUBLE, &cfirb_prior_sigma2);
    H5Aclose(att_id);
  }

  // SIGMA PRIOR1
  att_id = H5Acreate2(objid, "has_sigma_prior1", H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  bl = static_cast<hbool_t>(has_sigma_prior1);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, &bl);
  H5Aclose(att_id);
  if (has_sigma_prior1) {
    att_id = H5Acreate2(objid, "sigma_prior_width1", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_DOUBLE, &sigma_prior_width1);
    H5Aclose(att_id);
  }

  // SIGMA PRIOR2
  att_id = H5Acreate2(objid, "has_sigma_prior2", H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  bl = static_cast<hbool_t>(has_sigma_prior2);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, &bl);
  H5Aclose(att_id);
  if (has_sigma_prior1) {
    att_id = H5Acreate2(objid, "sigma_prior_width2", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_DOUBLE, &sigma_prior_width2);
    H5Aclose(att_id);
  }

  H5Sclose(mems_id);

  // Model info
  model.writeToHDF5Handle(objid);
}

/*!
  \param[in] comm Communicator
  \param[in] dest Destination of messages
*/
void calcLikeDouble::sendSelf(MPI_Comm comm, int dest) const { 
  //Transform
  MPI_Send(const_cast<unsigned int*>(&fftsize), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::CLDSENDFFTSIZE, comm);
  MPI_Send(const_cast<bool*>(&edgeInteg), 1, MPI::BOOL, dest,
	   pofd_mcmc::CLDSENDEDGEINTEG, comm);
  MPI_Send(const_cast<unsigned int*>(&nedge), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::CLDSENDNEDGE, comm);

  //Data
  MPI_Send(const_cast<unsigned int*>(&nbeamsets), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::CLDSENDNBEAM, comm);
  if (nbeamsets > 0) 
    for (unsigned int i = 0; i < nbeamsets; ++i) {
      MPI_Send(&i, 1, MPI_UNSIGNED, dest, pofd_mcmc::CLDSENDSETNUM, comm);
      beamsets[i].sendSelf(comm, dest);
    }
  
  MPI_Send(const_cast<bool*>(&bin_data), 1, MPI::BOOL, dest,
	   pofd_mcmc::CLDSENDBINDATA, comm);
  MPI_Send(const_cast<unsigned int*>(&nbins), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::CLDSENDNBINS, comm);

  //Model
  model.sendSelf(comm, dest);

  //CFIRB prior
  MPI_Send(const_cast<bool*>(&has_cfirb_prior1), 1, MPI::BOOL, dest,
	   pofd_mcmc::CLDSENDHASCFIRBPRIOR1, comm);
  if (has_cfirb_prior1) {
    MPI_Send(const_cast<double*>(&cfirb_prior_mean1), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDCFIRBPRIORMEAN1, comm);
    MPI_Send(const_cast<double*>(&cfirb_prior_sigma1), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDCFRIBPRIORSIGMA1, comm);
  }
  MPI_Send(const_cast<bool*>(&has_cfirb_prior2), 1, MPI::BOOL, dest,
	   pofd_mcmc::CLDSENDHASCFIRBPRIOR2, comm);
  if (has_cfirb_prior2) {
    MPI_Send(const_cast<double*>(&cfirb_prior_mean2), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDCFIRBPRIORMEAN2, comm);
    MPI_Send(const_cast<double*>(&cfirb_prior_sigma2), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDCFRIBPRIORSIGMA2, comm);
  }

  //Sigma prior
  MPI_Send(const_cast<bool*>(&has_sigma_prior1), 1, MPI::BOOL, dest,
	   pofd_mcmc::CLDSENDHASSIGMAPRIOR1, comm);
  if (has_sigma_prior1)
    MPI_Send(const_cast<double*>(&sigma_prior_width1), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDSIGMAPRIORWIDTH1, comm);
  MPI_Send(const_cast<bool*>(&has_sigma_prior2), 1, MPI::BOOL, dest,
	   pofd_mcmc::CLDSENDHASSIGMAPRIOR2, comm);
  if (has_sigma_prior2)
    MPI_Send(const_cast<double*>(&sigma_prior_width2), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLDSENDSIGMAPRIORWIDTH2, comm);

}

/*!
  \param[in] comm Communicator
  \param[in] src Source of messages
*/
void calcLikeDouble::recieveCopy(MPI_Comm comm, int src) {
  MPI_Status Info;

  //Transform
  MPI_Recv(&fftsize, 1, MPI_UNSIGNED, src, pofd_mcmc::CLDSENDFFTSIZE,
	   comm, &Info);
  MPI_Recv(&edgeInteg, 1, MPI::BOOL, src, pofd_mcmc::CLDSENDEDGEINTEG,
	   comm, &Info);
  MPI_Recv(&nedge, 1, MPI_UNSIGNED, src, pofd_mcmc::CLDSENDNEDGE,
	   comm, &Info);

  //Data
  unsigned int newnbeamsets;
  MPI_Recv(&newnbeamsets, 1, MPI_UNSIGNED, src, pofd_mcmc::CLDSENDNBEAM,
	   comm, &Info);
  if (newnbeamsets != nbeamsets) {
    if (beamsets != NULL) delete[] beamsets;
    if (newnbeamsets > 0) beamsets = new calcLikeDoubleSingle[newnbeamsets];
    else beamsets = NULL;
    nbeamsets = newnbeamsets;
  }
  unsigned int idx;
  for (unsigned int i = 0; i < nbeamsets; ++i) {
    MPI_Recv(&idx, 1, MPI_UNSIGNED, src, pofd_mcmc::CLDSENDSETNUM,
	     comm, &Info);
    beamsets[idx].recieveCopy(comm, src);
  }

  MPI_Recv(&bin_data, 1, MPI::BOOL, src, pofd_mcmc::CLDSENDBINDATA, 
	   comm, &Info);
  MPI_Recv(&nbins, 1, MPI_UNSIGNED, src, pofd_mcmc::CLDSENDNBINS,
	   comm, &Info);

  //Model
  model.recieveCopy(comm, src);

  //CFIRB prior
  MPI_Recv(&has_cfirb_prior1, 1, MPI::BOOL, src, 
	   pofd_mcmc::CLDSENDHASCFIRBPRIOR1, comm, &Info);
  if (has_cfirb_prior1) {
    MPI_Recv(&cfirb_prior_mean1, 1, MPI_DOUBLE, src,
	     pofd_mcmc::CLDSENDCFIRBPRIORMEAN1, comm, &Info);
    MPI_Recv(&cfirb_prior_sigma1, 1, MPI_DOUBLE, src,
	     pofd_mcmc::CLDSENDCFRIBPRIORSIGMA1, comm, &Info);
  }
  MPI_Recv(&has_cfirb_prior2, 1, MPI::BOOL, src, 
	   pofd_mcmc::CLDSENDHASCFIRBPRIOR2, comm, &Info);
  if (has_cfirb_prior2) {
    MPI_Recv(&cfirb_prior_mean2, 1, MPI_DOUBLE, src,
	     pofd_mcmc::CLDSENDCFIRBPRIORMEAN2, comm, &Info);
    MPI_Recv(&cfirb_prior_sigma2, 1, MPI_DOUBLE, src,
	     pofd_mcmc::CLDSENDCFRIBPRIORSIGMA2, comm, &Info);
  }

  //Sigma prior
  MPI_Recv(&has_sigma_prior1, 1, MPI::BOOL, src,
	   pofd_mcmc::CLDSENDHASSIGMAPRIOR1, comm, &Info);
  if (has_sigma_prior1) 
    MPI_Recv(&sigma_prior_width1, 1, MPI_DOUBLE, src,
	     pofd_mcmc::CLDSENDSIGMAPRIORWIDTH1, comm, &Info);
  MPI_Recv(&has_sigma_prior2, 1, MPI::BOOL, src,
	   pofd_mcmc::CLDSENDHASSIGMAPRIOR2, comm, &Info);
  if (has_sigma_prior2) 
    MPI_Recv(&sigma_prior_width2, 1, MPI_DOUBLE, src,
	     pofd_mcmc::CLDSENDSIGMAPRIORWIDTH2, comm, &Info);

}
