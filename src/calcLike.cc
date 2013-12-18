#include<limits>
#include<cmath>
#include<map>
#include<sstream>

#include "../include/calcLike.h"
#include "../include/affineExcept.h"

const double calcLikeSingle::bad_like = 1e25;
const double calcLikeSingle::flux_safety = 1.1;

/*!
  \param[in] NINTERP Number of interpolation points in R evaluation
*/
calcLikeSingle::calcLikeSingle(unsigned int NINTERP) :
  data_read(false), ndatasets(0), data(NULL), maxflux(0.0), 
  pdfac(NINTERP), sigma_base(NULL), maxsigma_base(0.0),
  exp_conf(0.0), like_norm(NULL), like_offset(NULL), has_beam(false), 
  verbose(false) {}

calcLikeSingle::~calcLikeSingle() {
  if (data != NULL) delete[] data;
  if (sigma_base != NULL)  delete[] sigma_base;
  if (like_norm   != NULL) delete[] like_norm;
  if (like_offset != NULL) delete[] like_offset;
}

void calcLikeSingle::free() {
  if (data != NULL) { delete[] data; data = NULL; }
  ndatasets = 0;
  data_read = false;
  maxflux = std::numeric_limits<double>::quiet_NaN();

  pd.strict_resize(0);
  pdfac.free();
  
  if (sigma_base != NULL) { delete[] sigma_base; sigma_base = NULL; }
  maxsigma_base = std::numeric_limits<double>::quiet_NaN();
  exp_conf = 0.0;

  if (like_norm != NULL) { delete[] like_norm; like_norm = NULL; }
  if (like_offset != NULL) { delete[] like_offset; like_offset = NULL; }

  has_beam = false;
  bm.free();
}

/*
  \param[in] n Number of datasets
*/
void calcLikeSingle::resize(unsigned int n) {
  if (ndatasets == n) return;  //Don't have to do anything

  if (data != NULL)        delete[] data;
  if (sigma_base != NULL)  delete[] sigma_base;
  if (like_offset != NULL) delete[] like_offset;
  if (like_norm   != NULL) delete[] like_norm;

  if (n > 0) {
    data        = new fitsData[n];
    sigma_base  = new double[n];
    like_offset = new double[n];
    like_norm   = new double[n];
    for (unsigned int i = 0; i < n; ++i)
      sigma_base[i] = std::numeric_limits<double>::quiet_NaN();
    for (unsigned int i = 0; i < n; ++i)
      like_offset[i] = std::numeric_limits<double>::quiet_NaN();
    for (unsigned int i = 0; i < n; ++i)
      like_norm[i] = 1.0;
  } else {
    data        = NULL;
    sigma_base  = NULL;
    like_offset = NULL;
    like_norm   = NULL;
  }
  
  data_read = false;
  maxflux   = std::numeric_limits<double>::quiet_NaN();
  maxsigma_base = 0.0;
  ndatasets = n;
}

/*!
  \param[in] datafile File to read data from
  \param[in] IGNOREMASK Ignore mask info in files
  \param[in] MEANSUB Subtract mean from data
  \param[in] BINDATA Bin the data
  \param[in] NBINS Number of data bins to use

  Special case for only a single data file
*/
void calcLikeSingle::readDataFromFile(const std::string& datafile, 
				      bool IGNOREMASK, bool MEANSUB,
				      bool BINDATA, unsigned int NBINS) {
  resize(1);

  data[0].readData(datafile, IGNOREMASK, MEANSUB);
  unsigned int nd = data[0].getN();
  if (nd == 0)
    throw affineExcept("calcLikeSingle", "readDataFromFile",
		       "No unmasked pixels in data", 1);
  if (BINDATA) data[0].applyBinning(NBINS);
  
  maxflux = data[0].getMax();
  if (std::isnan(maxflux) || std::isinf(maxflux))
    throw affineExcept("calcLikeSingle", "readDataFromFile",
		       "Problem with maxflux", 2);

  // We can't set the likelihood offset unless sigma base is set.
  like_offset[0] = 0.0;
  data_read = true;
}

/*!
  \param[in] datafiles Files to read data from
  \param[in] IGNOREMASK Ignore mask info in files
  \param[in] MEANSUB Subtract mean from data
  \param[in] BINDATA Bin the data
  \param[in] NBINS Number of bins to use
*/
void calcLikeSingle::
readDataFromFiles(const std::vector<std::string>& datafiles, 
		  bool IGNOREMASK, bool MEANSUB, bool BINDATA,
		  unsigned int NBINS) {
  unsigned int n = datafiles.size();
  if (n == 0)
    throw affineExcept("calcLikeSingle", "readDataFromFile",
		       "No data sets", 1);
  resize(n);

  for (unsigned int i = 0; i < n; ++i) {
    data[i].readData(datafiles[i], IGNOREMASK, MEANSUB);
    unsigned int nd = data[i].getN();
    if (nd == 0) {
      std::stringstream errstr("");
      errstr << "No unmasked pixels in data file: "
	     << datafiles[i];
      throw affineExcept("calcLikeSingle", "readDataFromFile",
			 errstr.str(), 3);
    }
    if (BINDATA) data[i].applyBinning(NBINS);

    // We can't set like_offsets until the sigma base values are set
    like_offset[i] = 0.0;
  }

  //Determine maximum flux (with safety factor) for all this data
  double cmaxflux;
  maxflux = data[0].getMax();
  if (std::isnan(maxflux) || std::isinf(maxflux))
    throw affineExcept("calcLikeSingle", "readDataFromFile",
		       "Problem with maxflux", 4);
  for (unsigned int i = 1; i < n; ++i) {
    cmaxflux = data[i].getMax();
    if (std::isnan(cmaxflux) || std::isinf(cmaxflux))
      throw affineExcept("calcLikeSingle", "readDataFromFile",
			 "Problem with maxflux", 5);
    if (cmaxflux > maxflux) maxflux=cmaxflux;
  }

  data_read = true;
}

/*!
  \param[in] fl Name of FITS file to read beam from
  \param[in] MINVAL Minimum bin value used
  \param[in] histogram Apply beam histogramming
  \param[in] NBINS Number of beam histogram bins to use
*/
void calcLikeSingle::readBeam(const std::string& fl, double MINVAL,
			      bool histogram, unsigned int NBINS) {
  bm.readFile(fl, MINVAL);
  if (histogram) bm.makeHistogram(NBINS);
  has_beam = true;
}

/*!
  \param[in] nbins Number of bins
*/
void calcLikeSingle::applyBinning(unsigned int nbins) {
  if ((!data_read) || ndatasets == 0) return;
  //Does nothing if the data is already binned at the same size
  for (unsigned int i = 0; i < ndatasets; ++i)
    data[i].applyBinning(nbins);
}

void calcLikeSingle::removeBinning() {
  if ((!data_read) || ndatasets == 0) return;
  for (unsigned int i = 0; i < ndatasets; ++i)
    data[i].removeBinning();
}

/*!
  \param[in] lnorm Vector of likelihood normalization values, of length
    the number of datasets
*/
void calcLikeSingle::setLikeNorm(const std::vector<double>& lnorm) {
  unsigned int n = lnorm.size();
  if (n != ndatasets)
    throw affineExcept("calcLikeSingle", "setLikeNorm",
		       "like_norm vector not same size as number of data sets",
		       1);
  if (n == 0) return;
  for (unsigned int i = 0; i < ndatasets; ++i)
    like_norm[i] = lnorm[i];
}

/*!
  \param[in] n Length of lnorm -- must be same as number of data sets
  \param[in] lnorm Array of likelihood normalization values
*/
void calcLikeSingle::setLikeNorm(unsigned int n, 
				 const double* const lnorm) {
  if (n != ndatasets)
    throw affineExcept("calcLikeSingle", "setLikeNorm",
		       "like_norm array not same size as number of data sets",
		       1);
  if (n == 0) return;
  for (unsigned int i = 0; i < ndatasets; ++i)
    like_norm[i] = lnorm[i];
}

/*!
  \param[in] s Vector of instrument sigma base values, of length
    the number of datasets
*/
void calcLikeSingle::setSigmaBase(const std::vector<double>& s) {
  unsigned int n = s.size();
  if (n != ndatasets)
    throw affineExcept("calcLikeSingle", "setSigmaBase",
		       "sigma vectors not same size as number of data sets", 1);
  if (n == 0) return;
  for (unsigned int i = 0; i < ndatasets; ++i)
    sigma_base[i] = s[i];
  maxsigma_base = sigma_base[0];
  for (unsigned int i = 1; i < ndatasets; ++i)
    if (sigma_base[i] > maxsigma_base) maxsigma_base = sigma_base[i];

  // Try to set a reasonable likelihood expectation value
  // In Glenn et al. we used log(N!), but this doesn't seem to work well
  // so instead use a Gaussian expectation for E( \sum_pix ln P )
  if (data_read) {
    double var;
    for (unsigned int i = 0; i < ndatasets; ++i) {
      var = sigma_base[i] * sigma_base[i] + exp_conf * exp_conf;
      if (var > 0) {
	double ndd = static_cast<double>(data[i].getN());
	like_offset[i] = -0.5 * ndd * (log(mcmc_affine::two_pi * var) + 1.0);
      }
    }
  }
}

/*!
  \param[in] n Number of elements in s, must be same as number of datasets
  \param[in] s Array of instrument sigma base values
*/
void calcLikeSingle::setSigmaBase(unsigned int n, const double* const s) {
  if (n != ndatasets)
    throw affineExcept("calcLikeSingle", "setSigmaBase",
		       "sigma arrays not same size as number of data sets", 1);
  if (n == 0) return;
  for (unsigned int i = 0; i < ndatasets; ++i)
    sigma_base[i] = s[i];
  for (unsigned int i = 1; i < ndatasets; ++i)
    if (sigma_base[i] > maxsigma_base) maxsigma_base = sigma_base[i];

  // Try to set a reasonable likelihood expectation value
  // In Glenn et al. we used log(N!), but this doesn't seem to work well
  // so instead use a Gaussian expectation for E( \sum_pix ln P )
  if (data_read) {
    double var;
    for (unsigned int i = 0; i < ndatasets; ++i) {
      var = sigma_base[i] * sigma_base[i] + exp_conf * exp_conf;
      if (var > 0) {
	double ndd = static_cast<double>(data[i].getN());
	like_offset[i] = -0.5 * ndd * (log(mcmc_affine::two_pi * var) + 1.0);
      }
    }
  }
}

/*
  \param[in] model   Model parameters must already be set
  \param[out] pars_invalid Set to true if parameters are determined to 
                            be invalid
  \param[in] sigmul  Sigma multiplier
  \param[in] fftsize Size of FFT to use
  
  \returns The Log Likelihood of the model relative to the data.

  This is the guts of the operation, computing the P(D) and
  finding the log Likelihood

  Any auxilliary parameters (positions of knots, etc.) must already
  be set in model.  The last model parameter is the sigma multiplier,
  the previous ones are the values of the number counts at the knots
*/
double 
calcLikeSingle::getLogLike(const numberCounts& model, bool& pars_invalid,
			   double sigmult, unsigned int fftsize, 
			   bool edgefix) const {

  if (!data_read)
    throw affineExcept("calcLikeSingle", "getNegLogLike",
		       "Data not read", 1);
  if (!has_beam)
    throw affineExcept("calcLikeSingle", "getNegLogLike",
		       "Beam not read", 2);

  //Maximum noise value with multiplier 
  double max_sigma = maxsigma_base * sigmult;
  if (max_sigma < 0.0) return calcLikeSingle::bad_like;

  // Initialize P(D)
  // We have to decide what maxflux to ask for.  This will
  // be the larger of the data maximum flux or the
  // highest knot plus some padding.  This assumes that
  // model.getMaxFlux doesn't change
  double modelmax = model.getMaxFlux() + 2 * maxsigma_base;
  double maxRflux = maxflux > modelmax ? maxflux : modelmax;
  maxRflux *= calcLikeSingle::flux_safety;
  pars_invalid = ! pdfac.initPD(fftsize, max_sigma, maxRflux, model, bm);

  if (pars_invalid) return 0.0; // No point in continuing

  //Compute likelihood of each bit of data
  double curr_LogLike, LogLike;
  LogLike = 0.0;
  for (unsigned int i = 0; i < ndatasets; ++i) {
    // Get PD for this particuar set of sigmas
    pdfac.getPD(sigmult * sigma_base[i], pd, true, edgefix);

    // Get log like
    curr_LogLike = pd.getLogLike(data[i]);

    // Apply beam and zero offset factor
    LogLike += (curr_LogLike - like_offset[i]) / like_norm[i];
  }
  
  return LogLike;
}

/*!
  \param[inout] os Stream to write current P(D) to
*/
void calcLikeSingle::writePDToStream(std::ostream& os) const {
  PD cpy(pd);
  cpy.deLog();
  os << cpy;
}

/*!
  \param[in] comm Communicator
  \param[in] dest Destination of messages
*/
void calcLikeSingle::sendSelf(MPI_Comm comm, int dest) const {
  //Data
  MPI_Send(const_cast<unsigned int*>(&ndatasets), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::CLSENDNDATA, comm);
  
  MPI_Send(const_cast<bool*>(&data_read), 1, MPI::BOOL, dest,
	   pofd_mcmc::CLSENDDATAREAD, comm);
  if (data_read) {
    for (unsigned int i = 0; i < ndatasets; ++i)
      data[i].sendSelf(comm, dest);
    MPI_Send(sigma_base, ndatasets, MPI_DOUBLE, dest,
	     pofd_mcmc::CLSENDSIGMABASE, comm);
    MPI_Send(const_cast<double*>(&maxsigma_base), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLSENDMAXSIGMABASE, comm);
    MPI_Send(const_cast<double*>(&exp_conf), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLSENDEXPCONF, comm);
    MPI_Send(like_offset, ndatasets, MPI_DOUBLE, dest,
	     pofd_mcmc::CLSENDLIKEOFFSET, comm);
    MPI_Send(like_norm, ndatasets, MPI_DOUBLE, dest,
	     pofd_mcmc::CLSENDLIKENORM, comm);
    MPI_Send(const_cast<double*>(&maxflux), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLSENDMAXFLUX, comm);
  }

  //Beam
  MPI_Send(const_cast<bool*>(&has_beam), 1, MPI::BOOL, dest,
	   pofd_mcmc::CLSENDHASBEAM, comm);
  if (has_beam) bm.sendSelf(comm, dest);

  //PDFactory
  pdfac.sendSelf(comm, dest);

  //Note we don't send verbose
}

/*!
  \param[in] comm Communicator
  \param[in] src Source of messages
*/
void calcLikeSingle::recieveCopy(MPI_Comm comm, int src) {
  MPI_Status Info;

  //Data
  unsigned int newn;
  MPI_Recv(&newn, 1, MPI_UNSIGNED, src, pofd_mcmc::CLSENDNDATA, comm, &Info);
  resize(newn);

  bool newread = false;
  MPI_Recv(&newread, 1, MPI::BOOL, src, pofd_mcmc::CLSENDDATAREAD,
	   comm, &Info);
  if (newread) {
    for (unsigned int i = 0; i < newn; ++i)
      data[i].recieveCopy(comm,src);
    MPI_Recv(sigma_base, newn, MPI_DOUBLE, src, pofd_mcmc::CLSENDSIGMABASE,
	     comm, &Info);
    MPI_Recv(&maxsigma_base, 1, MPI_DOUBLE, src, pofd_mcmc::CLSENDMAXSIGMABASE,
	     comm, &Info);
    MPI_Recv(&exp_conf, 1, MPI_DOUBLE, src, pofd_mcmc::CLSENDEXPCONF,
	     comm, &Info);
    MPI_Recv(like_offset, newn, MPI_DOUBLE, src,
	     pofd_mcmc::CLSENDLIKEOFFSET, comm, &Info);
    MPI_Recv(like_norm, newn, MPI_DOUBLE, src, pofd_mcmc::CLSENDLIKENORM,
	     comm, &Info);
    MPI_Recv(&maxflux, 1, MPI_DOUBLE, src, pofd_mcmc::CLSENDMAXFLUX,
	     comm, &Info);
    data_read = true;
  } 

  //Beam
  bool hsbm;
  MPI_Recv(&hsbm, 1, MPI::BOOL, src, pofd_mcmc::CLSENDHASBEAM,
	   comm, &Info);
  if (hsbm) {
    bm.recieveCopy(comm, src);
    has_beam = true;
  } else has_beam = false;
  
  //PDFactory
  pdfac.recieveCopy(comm, src);

}

/////////////////////////////////////////////////////////////////

/*!
  \param[in] FFTSIZE Number of elements in FFT
  \param[in] NINTERP Number of interpolation elements in R calculation
  \param[in] EDGEFIX Apply edge fix
  \param[in] BINNED Bin data
  \param[in] NBINS Number of data bins, if binning
*/
calcLike::calcLike(unsigned int FFTSIZE, unsigned int NINTERP, 
		   bool EDGEFIX, bool BINNED, unsigned int NBINS):
  fftsize(FFTSIZE), ninterp(NINTERP), edgeFix(EDGEFIX), bin_data(BINNED),
  nbins(NBINS), has_cfirb_prior(false), cfirb_prior_mean(0.0),
  cfirb_prior_sigma(0.0), has_sigma_prior(false), sigma_prior_width(0.0),
  verbose(false) {

  nbeamsets = 0;
  beamsets  = NULL;
}

calcLike::~calcLike() {
  if (beamsets != NULL) delete[] beamsets;
}

/*!
  This doesn't actually deallocate beamsets, just frees the internal
  data storage
*/
void calcLike::freeData() {
  for (unsigned int i = 0; i < nbeamsets; ++i)
    beamsets[i].free();
}

/*!
  \param[in] filename Name of FFTW wisdom file

  Don't do this before reading in the data files, or it will be overwritten
*/
void calcLike::addWisdom(const std::string& filename) {
  for (unsigned int i = 0; i < nbeamsets; ++i)
    beamsets[i].addWisdom(filename);
}

/*!
  \param[in] datafiles Vector of data file names
  \param[in] beamfiles Vector of beam file names
  \param[in] sigmas Vector of instrumental sigmas
  \param[in] like_norms Vector of likelihood normalizations
  \param[in] IGNOREMASK Ignore any mask in data files
  \param[in] MEANSUB Mean subtract the data
  \param[in] MINBEAMVAL Minimum beam value used
  \param[in] HISTOGRAMBEAMS Histogram the beams
  \param[in] NBEAMHIST Number of beam histogram bins to use
  \param[in] EXPCONF Expected confusion noise value

  Read in a set of data, grouping the inputs by beam and storing
  all of the relevant instrument noise and likelihood normalization values
*/
void calcLike::readDataFromFiles(const std::vector<std::string>& datafiles, 
				 const std::vector<std::string>& beamfiles,
				 const std::vector<double>& sigmas,
				 const std::vector<double>& like_norms,
				 bool IGNOREMASK, bool MEANSUB,
				 double MINBEAMVAL, bool HISTOGRAM, 
				 unsigned int NBEAMHIST, double EXPCONF) {

  //Make sure they are all the same length
  unsigned int ndat = datafiles.size();
  if (ndat == 0)
    throw affineExcept("calcLike", "readDataFromFiles",
		       "No datafiles", 1);
  if (ndat != beamfiles.size())
    throw affineExcept("calcLike", "readDataFromFiles",
		       "datafiles and beamfiles not same length", 2);
  if (ndat != sigmas.size())
    throw affineExcept("calcLike", "readDataFromFiles",
		       "datafiles and sigma not same length", 3);
  if (ndat != like_norms.size())
    throw affineExcept("calcLike", "readDataFromFiles",
		       "datafiles and like_norm not same length", 4);
  
  //Use a map.  Could also use a multi-map, but meh
  std::map< std::string, beam_group > grpmap;
  std::map< std::string, beam_group >::iterator grpmap_it;
  std::string key;
  for (unsigned int i = 0; i < ndat; ++i) {
    key = beamfiles[i];
    //See if map already has key
    grpmap_it = grpmap.find(key);
    if (grpmap_it == grpmap.end()) {
      //Previously unknown key
      beam_group newgrp;
      newgrp.n = 1;
      newgrp.datafiles.push_back(datafiles[i]);
      newgrp.beamfile = key;
      newgrp.sigmas.push_back(sigmas[i]);
      newgrp.like_norms.push_back(like_norms[i]);
      grpmap[key] = newgrp;
    } else {
      //Previously known key
      grpmap_it->second.n += 1;
      grpmap_it->second.datafiles.push_back(datafiles[i]);
      grpmap_it->second.sigmas.push_back(sigmas[i]);
      grpmap_it->second.like_norms.push_back(like_norms[i]);
    }
  }

  unsigned int newnbeamsets = grpmap.size();
  if (newnbeamsets != nbeamsets) {
    if (beamsets != NULL) delete[] beamsets;
    if (newnbeamsets > 0) beamsets = new calcLikeSingle[newnbeamsets];
    else beamsets = NULL;
    nbeamsets = newnbeamsets;
  }

  if (nbeamsets > 0) {
    grpmap_it = grpmap.begin();
    for (unsigned int i=0; grpmap_it != grpmap.end(); ++grpmap_it, ++i) {
      beamsets[i].setNInterp(ninterp);
      beamsets[i].readDataFromFiles(grpmap_it->second.datafiles,
				    IGNOREMASK, MEANSUB, bin_data, nbins);
      beamsets[i].readBeam(grpmap_it->second.beamfile, MINBEAMVAL,
			   HISTOGRAM, NBEAMHIST);
      beamsets[i].setExpConf(EXPCONF);
      beamsets[i].setSigmaBase(grpmap_it->second.sigmas);
      beamsets[i].setLikeNorm(grpmap_it->second.like_norms);
    }
  }
}

/*! 
  \param[in] nint New interpolation length
*/
void calcLike::setNInterp(unsigned int nint) {
  if (nint == ninterp) return;
  if (beamsets != NULL)
    for (unsigned int i = 0; i < nbeamsets; ++i)
      beamsets[i].setNInterp(nint);
  ninterp = nint;
}

void calcLike::setBinData() {
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

void calcLike::unSetBinData() {
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
  \param[in] nbns Number of bins
*/
void calcLike::setNBins(unsigned int nbns) {
  if (nbns == nbins) return;
  
  if (nbeamsets == 0) {
    nbins = nbns;
    return;
  }

  if (bin_data)
    for (unsigned int i = 0; i < nbeamsets; ++i)
      beamsets[i].applyBinning(nbns);

  nbins = nbns;
}


/*!
  \param[in] mn Mean value of CFIRB prior
  \param[in] sg Sigma of CFIRB prior
  
  The prior is assumed Gaussian
*/
void calcLike::setCFIRBPrior(double mn, double sg) {
  has_cfirb_prior = true;
  cfirb_prior_mean = mn;
  cfirb_prior_sigma = sg;
}

/*!
  \param[in] p Parameters to evaluate
  \param[out] pars_invalid True if parameters were not valid,
   otherwise False
*/
double calcLike::getLogLike(const paramSet& p, bool& pars_invalid) const {
  const double half_log_2pi = 0.918938517570495605469;

  if (nbeamsets == 0) return std::numeric_limits<double>::quiet_NaN();
  unsigned int npar = p.getNParams();
  if (npar < 2) throw affineExcept("calcLike", "getLogLike",
				   "Not enough elements in params", 1);
  unsigned int nknots = model.getNKnots();
  if (nknots == 0)
    throw affineExcept("calcLike", "getLogLike",
		       "Model knot positions not loaded", 2);
  if (nknots > (npar-1))
    throw affineExcept("calcLike", "getLogLike",
		       "Not enough elements in paramSet", 3);

  //Set the model
  model.setParams(p);

  double LogLike = 0.0;

  //Do the datasets likelihood
  pars_invalid = false;
  double sigmult = p[nknots]; //Sigma multiplier is after knot values
  bool pinvalid;
  for (unsigned int i = 0; i < nbeamsets; ++i) {
    LogLike += beamsets[i].getLogLike(model, pinvalid, sigmult, 
				      fftsize, edgeFix);
    pars_invalid &= pinvalid;
  }

  if (pars_invalid) return LogLike; //!< Not much point in doing the priors...

  //Add on cfirb prior and sigma prior if needed
  //Only do this once for all data sets so as not to multi-count the prior
  if (has_sigma_prior) {
    //Assume the mean is always at 1 -- otherwise, the
    // user would have specified different noise level
    double val = (sigmult-1.0) / sigma_prior_width;
    LogLike -= half_log_2pi + log(sigma_prior_width) + 
      0.5*val*val;
  }

  if (has_cfirb_prior) {
    double s_per_area = model.getFluxPerArea();
    double val = (cfirb_prior_mean - s_per_area) / cfirb_prior_sigma;
    LogLike -=  half_log_2pi + log(cfirb_prior_sigma) + 0.5*val*val;
  }
  return LogLike;
}

/*!						
  \param[in] objid HDF5 handle to write information to.  Must already be
    open
*/
void calcLike::writeToHDF5Handle(hid_t objid) const {
  // Writes some meta information
  hsize_t adims;
  hid_t mems_id, att_id;

  if (H5Iget_ref(objid) < 0)
    throw affineExcept("calcLike", "writeToHDF5",
		       "Input handle is not valid", 1);

  // FFTSIZE
  adims = 1;
  mems_id = H5Screate_simple(1, &adims, NULL);
  att_id = H5Acreate2(objid, "fftsize", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &fftsize);
  H5Aclose(att_id);

  // NINTERP
  att_id = H5Acreate2(objid, "ninterp", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &ninterp);
  H5Aclose(att_id);

  // NBEAMSETS
  att_id = H5Acreate2(objid, "nbeamsets", H5T_NATIVE_UINT,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(att_id, H5T_NATIVE_UINT, &nbeamsets);
  H5Aclose(att_id);

  // DATA BINNING
  hbool_t bl;
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

  // CFIRB PRIOR
  att_id = H5Acreate2(objid, "has_cfirb_prior", H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  bl = static_cast<hbool_t>(has_cfirb_prior);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, &bl);
  H5Aclose(att_id);
  if (has_cfirb_prior) {
    att_id = H5Acreate2(objid, "cfirb_prior_mean", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_DOUBLE, &cfirb_prior_mean);
    H5Aclose(att_id);
    att_id = H5Acreate2(objid, "cfirb_prior_sigma", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_DOUBLE, &cfirb_prior_sigma);
    H5Aclose(att_id);
  }

  // SIGMA PRIOR
  att_id = H5Acreate2(objid, "has_sigma_prior", H5T_NATIVE_HBOOL,
		      mems_id, H5P_DEFAULT, H5P_DEFAULT);
  bl = static_cast<hbool_t>(has_sigma_prior);
  H5Awrite(att_id, H5T_NATIVE_HBOOL, &bl);
  H5Aclose(att_id);
  if (has_sigma_prior) {
    att_id = H5Acreate2(objid, "sigma_prior_width", H5T_NATIVE_DOUBLE,
			mems_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(att_id, H5T_NATIVE_DOUBLE, &sigma_prior_width);
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
void calcLike::sendSelf(MPI_Comm comm, int dest) const { 
  //Transform
  MPI_Send(const_cast<unsigned int*>(&fftsize), 1, MPI_UNSIGNED,
	   dest, pofd_mcmc::CLSENDFFTSIZE, comm);
  MPI_Send(const_cast<unsigned int*>(&ninterp), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::CLSENDNINTERP, comm);
  MPI_Send(const_cast<bool*>(&edgeFix), 1, MPI::BOOL, dest,
	   pofd_mcmc::CLSENDEDGEFIX, comm);

  //Data
  MPI_Send(const_cast<unsigned int*>(&nbeamsets), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::CLSENDNBEAM, comm);
  if (nbeamsets > 0) 
    for (unsigned int i = 0; i < nbeamsets; ++i) {
      MPI_Send(const_cast<unsigned int*>(&i), 1, MPI_UNSIGNED, dest,
	       pofd_mcmc::CLSENDSETNUM, comm);
      beamsets[i].sendSelf(comm,dest);
    }

  MPI_Send(const_cast<bool*>(&bin_data), 1, MPI::BOOL, dest,
	   pofd_mcmc::CLSENDBINDATA, comm);
  MPI_Send(const_cast<unsigned int*>(&nbins), 1, MPI_UNSIGNED, dest,
	   pofd_mcmc::CLSENDNBINS, comm);

  //Model
  model.sendSelf(comm,dest);

  //CFIRB prior
  MPI_Send(const_cast<bool*>(&has_cfirb_prior), 1, MPI::BOOL, dest,
	   pofd_mcmc::CLSENDHASCFIRBPRIOR, comm);
  if (has_cfirb_prior) {
    MPI_Send(const_cast<double*>(&cfirb_prior_mean), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLSENDCFIRBPRIORMEAN, comm);
    MPI_Send(const_cast<double*>(&cfirb_prior_sigma), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLSENDCFIRBPRIORSIGMA, comm);
  }

  //Sigma prior
  MPI_Send(const_cast<bool*>(&has_sigma_prior), 1, MPI::BOOL, dest,
	   pofd_mcmc::CLSENDHASSIGMAPRIOR, comm);
  if (has_sigma_prior)
    MPI_Send(const_cast<double*>(&sigma_prior_width), 1, MPI_DOUBLE, dest,
	     pofd_mcmc::CLSENDSIGMAPRIORWIDTH, comm);

}

/*!
  \param[in] comm Communicator
  \param[in] src Source of messages
*/
void calcLike::recieveCopy(MPI_Comm comm, int src) {
  MPI_Status Info;

  //Transform
  MPI_Recv(&fftsize, 1, MPI_UNSIGNED, src, pofd_mcmc::CLSENDFFTSIZE,
	   comm, &Info);
  MPI_Recv(&ninterp, 1, MPI_UNSIGNED, src, pofd_mcmc::CLSENDNINTERP,
	   comm, &Info);
  MPI_Recv(&edgeFix, 1, MPI::BOOL, src, pofd_mcmc::CLSENDEDGEFIX,
	   comm, &Info);

  //Data
  unsigned int newnbeamsets;
  MPI_Recv(&newnbeamsets, 1, MPI_UNSIGNED, src, pofd_mcmc::CLSENDNBEAM,
	   comm, &Info);
  if (newnbeamsets != nbeamsets) {
    if (beamsets != NULL) delete[] beamsets;
    if (newnbeamsets > 0) beamsets = new calcLikeSingle[newnbeamsets];
    else beamsets = NULL;
    nbeamsets = newnbeamsets;
  }
  unsigned int idx;
  for (unsigned int i = 0; i < nbeamsets; ++i) {
    MPI_Recv(&idx, 1, MPI_UNSIGNED, src, pofd_mcmc::CLSENDSETNUM,
	     comm, &Info);
    beamsets[idx].recieveCopy(comm, src);
  }

  MPI_Recv(&bin_data, 1, MPI::BOOL, src, pofd_mcmc::CLSENDBINDATA,
	   comm, &Info);
  MPI_Recv(&nbins, 1, MPI_UNSIGNED, src, pofd_mcmc::CLSENDNBINS,
	   comm, &Info);

  //Model
  model.recieveCopy(comm, src);

  //CFIRB prior
  MPI_Recv(&has_cfirb_prior, 1, MPI::BOOL, src, pofd_mcmc::CLSENDHASCFIRBPRIOR,
	   comm, &Info);
  if (has_cfirb_prior) {
    MPI_Recv(&cfirb_prior_mean, 1, MPI_DOUBLE, src,
	     pofd_mcmc::CLSENDCFIRBPRIORMEAN, comm, &Info);
    MPI_Recv(&cfirb_prior_sigma, 1, MPI_DOUBLE, src,
	     pofd_mcmc::CLSENDCFIRBPRIORSIGMA, comm, &Info);
  }
  //Sigma prior
  MPI_Recv(&has_sigma_prior, 1, MPI::BOOL, src, 
	   pofd_mcmc::CLSENDHASSIGMAPRIOR, comm, &Info);
  if (has_sigma_prior) 
    MPI_Recv(&sigma_prior_width, 1, MPI_DOUBLE, src,
	     pofd_mcmc::CLSENDSIGMAPRIORWIDTH, comm, &Info);
}

