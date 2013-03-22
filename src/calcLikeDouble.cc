#include<limits>
#include<cmath>
#include<map>
#include<sstream>

#include "../include/calcLikeDouble.h"
#include "../include/affineExcept.h"

const double calcLikeDoubleSingle::bad_like = 1e25;
const double calcLikeDoubleSingle::flux_safety = 1.2;

calcLikeDoubleSingle::calcLikeDoubleSingle(unsigned int NEDGE) :
  data_read(false), ndatasets(0), data(NULL), maxflux1(0.0), maxflux2(0.0),
  pdfac(NEDGE), sigma_base1(NULL), sigma_base2(NULL), maxsigma_base1(0.0),
  maxsigma_base2(0.0), like_norm(NULL), like_offset(NULL), 
  has_beam(false), verbose(false) {}

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
  maxflux1 = std::numeric_limits<double>::quiet_NaN();
  maxflux2 = std::numeric_limits<double>::quiet_NaN();

  pd.strict_resize(0, 0);
  pdfac.free();
  
  if (sigma_base1 != NULL) { delete[] sigma_base1; sigma_base1 = NULL; }
  maxsigma_base1 = std::numeric_limits<double>::quiet_NaN();
  if (sigma_base2 != NULL) { delete[] sigma_base2; sigma_base2 = NULL; }
  maxsigma_base2 = std::numeric_limits<double>::quiet_NaN();

  if (like_norm != NULL) { delete[] like_norm; like_norm = NULL; }
  if (like_offset != NULL) { delete[] like_offset; like_offset = NULL; }

  has_beam = false;
  bm.free();
}


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
  maxflux1  = std::numeric_limits<double>::quiet_NaN();
  maxflux2  = std::numeric_limits<double>::quiet_NaN();
  maxsigma_base1 = 0.0;
  maxsigma_base2 = 0.0;
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

  data[0].readData(datafile1,datafile2,IGNOREMASK,MEANSUB);
  unsigned int n = data[0].getN();
  if (n == 0)
    throw affineExcept("calcLikeDoubleSingle","readDataFromFile",
		       "No unmasked pixels in data",1);
  if (BINDATA) data[0].applyBinning(NBINS,NBINS);
  
  std::pair<double,double> pr = data[0].getMax();
  maxflux1 = pr.first; maxflux2 = pr.second;
  if (std::isnan(maxflux1) || std::isinf(maxflux1))
    throw affineExcept("calcLikeDoubleSingle","readDataFromFile",
		       "Problem with maxflux, band 1",2);
  if (std::isnan(maxflux2) || std::isinf(maxflux2))
    throw affineExcept("calcLikeDoubleSingle","readDataFromFile",
		       "Problem with maxflux, band 2",3);
  maxflux1 *= calcLikeDoubleSingle::flux_safety;
  maxflux2 *= calcLikeDoubleSingle::flux_safety;

  //log(N!)
  like_offset[0] = n*log(n) - n + 0.5*log(2.0*M_PI*n) + 
    1.0/(12.0*n) - 1.0/(360.0*n*n*n) + 1/(1260.0*n*n*n*n*n) -
    1.0/(1680.0*pow(n,7));

  data_read = true;
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
    throw affineExcept("calcLikeDoubleSingle","readDataFromFiles",
		       "No data sets",1);
  if (datafiles2.size() != n)
    throw affineExcept("calcLikeDoubleSingle","readDataFromFiles",
		       "Different number of datafiles in each band",2);
  resize(n);

  for (unsigned int i = 0; i < n; ++i) {
    data[i].readData(datafiles1[i],datafiles2[i],IGNOREMASK,MEANSUB);
    unsigned int nd = data[i].getN();
    if (nd == 0) {
      std::stringstream errstr("");
      errstr << "No unmasked pixels in data files: "
	     << datafiles1[i] << " " << datafiles2[i];
      throw affineExcept("calcLikeDoubleSingle","readDataFromFiles",
			 errstr.str(),3);
    }
    if (BINDATA) data[i].applyBinning(NBINS,NBINS);

    //log(N!), a bit beyond the stirling approximation
    like_offset[i] = nd * log(nd) - nd + 
      0.5 * log(2.0 * M_PI * nd) + 
      1.0 / (12.0*nd) - 1.0 / (360.0*nd*nd*nd) + 
      1.0 / (1260.0*nd*nd*nd*nd*nd) -
      1.0 / (1680.0*pow(nd,7));
  }

  //Determine maximum flux (with safety factor) for all this data
  std::pair<double,double> pr;
  pr = data[0].getMax();
  maxflux1 = pr.first; maxflux2 = pr.second;
  if (std::isnan(maxflux1) || std::isinf(maxflux1))
    throw affineExcept("calcLikeDoubleSingle","readDataFromFiles",
		       "Problem with maxflux, band 1",4);
  if (std::isnan(maxflux2) || std::isinf(maxflux2))
    throw affineExcept("calcLikeDoubleSingle","readDataFromFiles",
		       "Problem with maxflux, band 2",5);
  for (unsigned int i = 1; i < n; ++i) {
    pr = data[i].getMax();
    if (std::isnan(pr.first) || std::isinf(pr.first))
      throw affineExcept("calcLikeDoubleSingle","readDataFromFiles",
			 "Problem with maxflux, band 1",6);
    if (std::isnan(pr.second) || std::isinf(pr.second))
      throw affineExcept("calcLikeDoubleSingle","readDataFromFiles",
			 "Problem with maxflux, band 2",7);
    if (pr.first > maxflux1) maxflux1 = pr.first;
    if (pr.second > maxflux2) maxflux2 = pr.second;
  }

  maxflux1 *= calcLikeDoubleSingle::flux_safety;
  maxflux2 *= calcLikeDoubleSingle::flux_safety;

  data_read = true;
}

void calcLikeDoubleSingle::readBeam(const std::string& fl1, 
				    const std::string& fl2,
				    bool histogram, 
				    double histogramlogstep) {
  bm.readFiles(fl1,fl2,histogram,histogramlogstep);
  has_beam = true;
}

/*!
  \param[in] nbins Number of bins
 */
void calcLikeDoubleSingle::applyBinning(unsigned int nbins) {
  if ((!data_read) || ndatasets == 0) return;
  //Does nothing if the data is already binned at the same size
  for (unsigned int i = 0; i < ndatasets; ++i)
    data[i].applyBinning(nbins,nbins);
}

void calcLikeDoubleSingle::removeBinning() {
  if ((!data_read) || ndatasets == 0) return;
  for (unsigned int i = 0; i < ndatasets; ++i)
    data[i].removeBinning();
}


void calcLikeDoubleSingle::setLikeNorm(const std::vector<double>& lnorm) {
  unsigned int n = lnorm.size();
  if (n != ndatasets)
    throw affineExcept("calcLikeDoubleSingle","setLikeNorm",
		       "like_norm vector not same size as number of data sets",
		       1);
  if (n == 0) return;
  for (unsigned int i = 0; i < ndatasets; ++i)
    like_norm[i] = lnorm[i];
}

void calcLikeDoubleSingle::setLikeNorm(unsigned int n, 
			   const double* const lnorm) {
  if (n != ndatasets)
    throw affineExcept("calcLikeDoubleSingle","setLikeNorm",
		       "like_norm array not same size as number of data sets",
		       1);
  if (n == 0) return;
  for (unsigned int i = 0; i < ndatasets; ++i)
    like_norm[i] = lnorm[i];
}

void calcLikeDoubleSingle::setSigmaBase1(const std::vector<double>& s) {
  unsigned int n = s.size();
  if (n != ndatasets)
    throw affineExcept("calcLikeDoubleSingle","setSigmaBase1",
		       "sigma vectors not same size as number of data sets",1);
  if (n == 0) return;
  for (unsigned int i = 0; i < ndatasets; ++i)
    sigma_base1[i] = s[i];
  maxsigma_base1 = sigma_base1[0];
  for (unsigned int i = 1; i < ndatasets; ++i)
    if (sigma_base1[i] > maxsigma_base1) maxsigma_base1 = sigma_base1[i];
}

void calcLikeDoubleSingle::setSigmaBase2(const std::vector<double>& s) {
  unsigned int n = s.size();
  if (n != ndatasets)
    throw affineExcept("calcLikeDoubleSingle","setSigmaBase2",
		       "sigma vectors not same size as number of data sets",1);
  if (n == 0) return;
  for (unsigned int i = 0; i < ndatasets; ++i)
    sigma_base2[i] = s[i];
  maxsigma_base2 = sigma_base2[0];
  for (unsigned int i = 1; i < ndatasets; ++i)
    if (sigma_base2[i] > maxsigma_base2) maxsigma_base2 = sigma_base2[i];
}


void calcLikeDoubleSingle::setSigmaBase1(unsigned int n,const double* const s) {
  if (n != ndatasets)
    throw affineExcept("calcLikeDoubleSingle","setSigmaBase1",
		       "sigma arrays not same size as number of data sets",1);
  if (n == 0) return;
  for (unsigned int i = 0; i < ndatasets; ++i)
    sigma_base1[i] = s[i];
  for (unsigned int i = 1; i < ndatasets; ++i)
    if (sigma_base1[i] > maxsigma_base1) maxsigma_base1 = sigma_base1[i];
}



void calcLikeDoubleSingle::setSigmaBase2(unsigned int n,const double* const s) {
  if (n != ndatasets)
    throw affineExcept("calcLikeDoubleSingle","setSigmaBase2",
		       "sigma arrays not same size as number of data sets",1);
  if (n == 0) return;
  for (unsigned int i = 0; i < ndatasets; ++i)
    sigma_base2[i] = s[i];
  for (unsigned int i = 2; i < ndatasets; ++i)
    if (sigma_base2[i] > maxsigma_base2) maxsigma_base2 = sigma_base2[i];
}

/*
  \param[in] model   Model parameters must already be set
  \param[in] sigmult1 Sigma multiplier in band 1
  \param[in] sigmult2 Sigma multiplier in band 2
  \param[in] fftsize Size of FFT to use
  \param[in] edgefix Apply edge fix to P(D) for wrapping
  \param[in] setedge Fill in edges of R with integral average

  This is the guts of the operation, computing the P(D) and
  finding the -log Likelihood

  Any auxilliary parameters (positions of knots, etc.) must already
  be set in model.  The last model parameter is the sigma multiplier,
  the previous ones are the values of the number counts at the knots
*/
double 
calcLikeDoubleSingle::getLogLike(const numberCountsDouble& model, 
				 double sigmult1, double sigmult2,
				 unsigned int fftsize, bool edgefix, 
				 bool setedge) const {

  if (!data_read)
    throw affineExcept("calcLikeDoubleSingle","getNegLogLike",
		       "Data not read",1);
  if (!has_beam)
    throw affineExcept("calcLikeDoubleSingle","getNegLogLike",
		       "Beam not read",2);

  //Maximum noise value with multiplier 
  double max_sigma1 = maxsigma_base1 * sigmult1;
  if (max_sigma1 < 0.0) return calcLikeDoubleSingle::bad_like;
  double max_sigma2 = maxsigma_base2 * sigmult2;
  if (max_sigma2 < 0.0) return calcLikeDoubleSingle::bad_like;

  //Initialize P(D)
  pdfac.initPD(fftsize,max_sigma1,max_sigma2, maxflux1, maxflux2, 
	       model, bm, setedge);

  //Compute likelihood of each bit of data
  double curr_LogLike, lnorm, LogLike, bmarea;
  bmarea = bm.getEffectiveArea1() * bm.getEffectiveArea2();
  LogLike = 0.0;
  for (unsigned int i = 0; i < ndatasets; ++i) {
    //Get PD for this particuar set of sigmas
    pdfac.getPD(sigmult1*sigma_base1[i], sigmult2*sigma_base2[i],
		 pd, true, edgefix);

    //Get log like
    curr_LogLike = pd.getLogLike(data[i]) - like_offset[i];

    //Beam norm
    lnorm = like_norm[i];
    if (lnorm > 0)
      curr_LogLike /= lnorm * bmarea;
    else if (lnorm < 0)
      curr_LogLike /= fabs(lnorm);

    LogLike += curr_LogLike;
  }
  
  return LogLike;
}

void calcLikeDoubleSingle::writePDToStream(std::ostream& os) const {
  PDDouble cpy(pd);
  cpy.deLog();
  os << cpy;
}

void calcLikeDoubleSingle::sendSelf(MPI::Comm& comm, int dest) const {
  //Data
  comm.Send(&ndatasets,1,MPI::UNSIGNED,dest,pofd_mcmc::CLDSENDNDATA);
  
  comm.Send(&data_read,1,MPI::BOOL,dest,pofd_mcmc::CLDSENDDATAREAD);
  if (data_read) {
    for (unsigned int i = 0; i < ndatasets; ++i)
      data[i].sendSelf(comm,dest);
    comm.Send(sigma_base1,ndatasets,MPI::DOUBLE,dest,
	      pofd_mcmc::CLDSENDSIGMABASE1);
    comm.Send(sigma_base2,ndatasets,MPI::DOUBLE,dest,
	      pofd_mcmc::CLDSENDSIGMABASE2);
    comm.Send(&maxsigma_base1,1,MPI::DOUBLE,dest,
	      pofd_mcmc::CLDSENDMAXSIGMABASE1);
    comm.Send(&maxsigma_base2,1,MPI::DOUBLE,dest,
	      pofd_mcmc::CLDSENDMAXSIGMABASE2);
    comm.Send(like_offset,ndatasets,MPI::DOUBLE,dest,
	      pofd_mcmc::CLDSENDLIKEOFFSET);
    comm.Send(like_norm,ndatasets,MPI::DOUBLE,dest,
	      pofd_mcmc::CLDSENDLIKENORM);
    comm.Send(&maxflux1,1,MPI::DOUBLE,dest,pofd_mcmc::CLDSENDMAXFLUX1);
    comm.Send(&maxflux2,1,MPI::DOUBLE,dest,pofd_mcmc::CLDSENDMAXFLUX2);
  }

  //Beam
  comm.Send(&has_beam,1,MPI::BOOL,dest,pofd_mcmc::CLDSENDHASBEAM);
  if (has_beam) bm.sendSelf(comm,dest);

  //PDFactory
  pdfac.sendSelf(comm,dest);

  //Note we don't send verbose
}


void calcLikeDoubleSingle::recieveCopy(MPI::Comm& comm, int src) {
  //Data
  unsigned int newn;
  comm.Recv(&newn,1,MPI::UNSIGNED,src,pofd_mcmc::CLDSENDNDATA);
  resize(newn);

  bool newread = false;
  comm.Recv(&newread,1,MPI::BOOL,src,pofd_mcmc::CLDSENDDATAREAD);
  if (newread) {
    for (unsigned int i = 0; i < newn; ++i)
      data[i].recieveCopy(comm,src);
    comm.Recv(sigma_base1,newn,MPI::DOUBLE,src,pofd_mcmc::CLDSENDSIGMABASE1);
    comm.Recv(sigma_base2,newn,MPI::DOUBLE,src,pofd_mcmc::CLDSENDSIGMABASE2);
    comm.Recv(&maxsigma_base1,1,MPI::DOUBLE,src,
	      pofd_mcmc::CLDSENDMAXSIGMABASE1);
    comm.Recv(&maxsigma_base2,1,MPI::DOUBLE,src,
	      pofd_mcmc::CLDSENDMAXSIGMABASE2);
    comm.Recv(like_offset,newn,MPI::DOUBLE,src,
	      pofd_mcmc::CLDSENDLIKEOFFSET);
    comm.Recv(like_norm,newn,MPI::DOUBLE,src,pofd_mcmc::CLDSENDLIKENORM);
    comm.Recv(&maxflux1,1,MPI::DOUBLE,src,pofd_mcmc::CLDSENDMAXFLUX1);
    comm.Recv(&maxflux2,1,MPI::DOUBLE,src,pofd_mcmc::CLDSENDMAXFLUX2);
    data_read = true;
  } 

  //Beam
  bool hsbm;
  comm.Recv(&hsbm,1,MPI::BOOL,src,pofd_mcmc::CLDSENDHASBEAM);
  if (hsbm) {
    bm.recieveCopy(comm,src);
    has_beam = true;
  } else has_beam = false;
  
  //PDFactory
  pdfac.recieveCopy(comm,src);

}


/////////////////////////////////////////////////////////////////

calcLikeDouble::calcLikeDouble(unsigned int FFTSIZE, unsigned int NEDGE, 
			       bool EDGEFIX, bool EDGEINTEG,
			       bool BINNED, unsigned int NBINS) :
  fftsize(FFTSIZE), nedge(NEDGE), edgeInteg(EDGEINTEG),
  edgeFix(EDGEFIX), bin_data(BINNED), nbins(NBINS), 
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
  Don't do this before reading in the data files, or it will be overwritten
 */
void calcLikeDouble::addWisdom(const std::string& filename) {
  for (unsigned int i = 0; i < nbeamsets; ++i)
    beamsets[i].addWisdom(filename);
}

/*!
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
		  bool IGNOREMASK, bool MEANSUB,
		  bool HISTOGRAM, double HISTOGRAMLOGSTEP) {

  //Make sure they are all the same length
  unsigned int ndat = datafiles1.size();
  if (ndat != datafiles2.size())
    throw affineExcept("calcLikeDouble","readDataFromFiles",
		       "Datafiles1 and datafile2 not same length",1);
  if (ndat == 0)
    throw affineExcept("calcLikeDouble","readDataFromFiles",
		       "No datafiles",2);
  if (ndat != beamfiles1.size())
    throw affineExcept("calcLikeDouble","readDataFromFiles",
		       "datafiles and beamfiles1 not same length",3);
  if (ndat != beamfiles2.size())
    throw affineExcept("calcLikeDouble","readDataFromFiles",
		       "datafiles and beamfiles2 not same length",4);
  if (ndat != sigmas1.size())
    throw affineExcept("calcLikeDouble","readDataFromFiles",
		       "datafiles and sigmas1 not same length",5);
  if (ndat != sigmas2.size())
    throw affineExcept("calcLikeDouble","readDataFromFiles",
		       "datafiles and sigmas2 not same length",6);
  if (ndat != like_norms.size())
    throw affineExcept("calcLikeDouble","readDataFromFiles",
		       "datafiles and like_norm not same length",7);
  
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
  grpmap_it = grpmap.begin();
  for (unsigned int i=0; grpmap_it != grpmap.end(); ++grpmap_it, ++i) {
    beamsets[i].setNEdge(nedge);
    beamsets[i].readDataFromFiles(grpmap_it->second.datafiles1,
				  grpmap_it->second.datafiles2,
				  IGNOREMASK, MEANSUB, bin_data, nbins);
    beamsets[i].readBeam(grpmap_it->second.beamfile1, 
			 grpmap_it->second.beamfile2, 
			 HISTOGRAM, HISTOGRAMLOGSTEP);
    beamsets[i].setSigmaBase1(grpmap_it->second.sigmas1);
    beamsets[i].setSigmaBase2(grpmap_it->second.sigmas2);
    beamsets[i].setLikeNorm(grpmap_it->second.like_norms);
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

double calcLikeDouble::getLogLike(const paramSet& p) const {
  const double half_log_2pi = 0.918938517570495605469;

  if (nbeamsets == 0) return std::numeric_limits<double>::quiet_NaN();
  unsigned int npar = p.getNParams();
  unsigned int nmodelpars = model.getNParams();
  if (npar < (nmodelpars+2)) //+2 for the sigma multiplier in each band
    throw affineExcept("calcLikeDouble","getLogLike",
		       "Not enough elements in params",1);

  //Set the model
  model.setParams(p);

  double LogLike = 0.0;

  //Do the datasets likelihood
  double sigmult1 = p[nmodelpars]; //Sigma multiplier is after knot values
  double sigmult2 = p[nmodelpars+1]; 
  for (unsigned int i = 0; i < nbeamsets; ++i)
    LogLike += beamsets[i].getLogLike(model,sigmult1,sigmult2,fftsize,
				      edgeFix,edgeInteg);

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
    double save_per_area = model.getMeanFluxPerArea(0);
    double val = (cfirb_prior_mean1 - save_per_area) / cfirb_prior_sigma1;
    LogLike -=  half_log_2pi + log(cfirb_prior_sigma1) + 0.5*val*val;
  }
  if (has_cfirb_prior2) {
    double save_per_area = model.getMeanFluxPerArea(1);
    double val = (cfirb_prior_mean2 - save_per_area) / cfirb_prior_sigma2;
    LogLike -=  half_log_2pi + log(cfirb_prior_sigma2) + 0.5*val*val;
  }

  return LogLike;
}

void calcLikeDouble::sendSelf(MPI::Comm& comm, int dest) const { 
  //Transform
  comm.Send(&fftsize,1,MPI::UNSIGNED,dest,pofd_mcmc::CLDSENDFFTSIZE);
  comm.Send(&edgeFix,1,MPI::BOOL,dest,pofd_mcmc::CLDSENDEDGEFIX);
  comm.Send(&edgeInteg,1,MPI::BOOL,dest,pofd_mcmc::CLDSENDEDGEINTEG);
  comm.Send(&nedge,1,MPI::UNSIGNED,dest,pofd_mcmc::CLDSENDNEDGE);

  //Data
  comm.Send(&nbeamsets,1,MPI::UNSIGNED,dest,pofd_mcmc::CLDSENDNBEAM);
  if (nbeamsets > 0) 
    for (unsigned int i = 0; i < nbeamsets; ++i) {
      comm.Send(&i,1,MPI::UNSIGNED,dest,pofd_mcmc::CLDSENDSETNUM);
      beamsets[i].sendSelf(comm,dest);
    }

  comm.Send(&bin_data,1,MPI::BOOL,dest,pofd_mcmc::CLDSENDBINDATA);
  comm.Send(&nbins,1,MPI::UNSIGNED,dest,pofd_mcmc::CLDSENDNBINS);

  //Model
  model.sendSelf(comm,dest);

  //CFIRB prior
  comm.Send(&has_cfirb_prior1,1,MPI::BOOL,dest,
	    pofd_mcmc::CLDSENDHASCFIRBPRIOR1);
  if (has_cfirb_prior1) {
    comm.Send(&cfirb_prior_mean1,1,MPI::DOUBLE,dest,
	      pofd_mcmc::CLDSENDCFIRBPRIORMEAN1);
    comm.Send(&cfirb_prior_sigma1,1,MPI::DOUBLE,dest,
	      pofd_mcmc::CLDSENDCFRIBPRIORSIGMA1);
  }
  comm.Send(&has_cfirb_prior2,1,MPI::BOOL,dest,
	    pofd_mcmc::CLDSENDHASCFIRBPRIOR2);
  if (has_cfirb_prior2) {
    comm.Send(&cfirb_prior_mean2,1,MPI::DOUBLE,dest,
	      pofd_mcmc::CLDSENDCFIRBPRIORMEAN2);
    comm.Send(&cfirb_prior_sigma2,1,MPI::DOUBLE,dest,
	      pofd_mcmc::CLDSENDCFRIBPRIORSIGMA2);
  }

  //Sigma prior
  comm.Send(&has_sigma_prior1,1,MPI::BOOL,dest,
	    pofd_mcmc::CLDSENDHASSIGMAPRIOR1);
  if (has_sigma_prior1)
    comm.Send(&sigma_prior_width1,1,MPI::DOUBLE,dest,
	      pofd_mcmc::CLDSENDSIGMAPRIORWIDTH1);
  comm.Send(&has_sigma_prior2,1,MPI::BOOL,dest,
	    pofd_mcmc::CLDSENDHASSIGMAPRIOR2);
  if (has_sigma_prior2)
    comm.Send(&sigma_prior_width2,1,MPI::DOUBLE,dest,
	      pofd_mcmc::CLDSENDSIGMAPRIORWIDTH2);

}

void calcLikeDouble::recieveCopy(MPI::Comm& comm, int src) {
  //Transform
  comm.Recv(&fftsize,1,MPI::UNSIGNED,src,pofd_mcmc::CLDSENDFFTSIZE);
  comm.Recv(&edgeFix,1,MPI::BOOL,src,pofd_mcmc::CLDSENDEDGEFIX);
  comm.Recv(&edgeInteg,1,MPI::BOOL,src,pofd_mcmc::CLDSENDEDGEINTEG);
  comm.Recv(&nedge,1,MPI::UNSIGNED,src,pofd_mcmc::CLDSENDNEDGE);

  //Data
  unsigned int newnbeamsets;
  comm.Recv(&newnbeamsets,1,MPI::UNSIGNED,src,pofd_mcmc::CLDSENDNBEAM);
  if (newnbeamsets != nbeamsets) {
    if (beamsets != NULL) delete[] beamsets;
    if (newnbeamsets > 0) beamsets = new calcLikeDoubleSingle[newnbeamsets];
    else beamsets = NULL;
    nbeamsets = newnbeamsets;
  }
  unsigned int idx;
  for (unsigned int i = 0; i < nbeamsets; ++i) {
    comm.Recv(&idx,1,MPI::UNSIGNED,src,pofd_mcmc::CLDSENDSETNUM);
    beamsets[idx].recieveCopy(comm,src);
  }

  comm.Recv(&bin_data,1,MPI::BOOL,src,pofd_mcmc::CLDSENDBINDATA);
  comm.Recv(&nbins,1,MPI::UNSIGNED,src,pofd_mcmc::CLDSENDNBINS);

  //Model
  model.recieveCopy(comm,src);

  //CFIRB prior
  comm.Recv(&has_cfirb_prior1,1,MPI::BOOL,src,pofd_mcmc::CLDSENDHASCFIRBPRIOR1);
  if (has_cfirb_prior1) {
    comm.Recv(&cfirb_prior_mean1,1,MPI::DOUBLE,src,
	      pofd_mcmc::CLDSENDCFIRBPRIORMEAN1);
    comm.Recv(&cfirb_prior_sigma1,1,MPI::DOUBLE,src,
	      pofd_mcmc::CLDSENDCFRIBPRIORSIGMA1);
  }
  comm.Recv(&has_cfirb_prior2,1,MPI::BOOL,src,pofd_mcmc::CLDSENDHASCFIRBPRIOR2);
  if (has_cfirb_prior2) {
    comm.Recv(&cfirb_prior_mean2,1,MPI::DOUBLE,src,
	      pofd_mcmc::CLDSENDCFIRBPRIORMEAN2);
    comm.Recv(&cfirb_prior_sigma2,1,MPI::DOUBLE,src,
	      pofd_mcmc::CLDSENDCFRIBPRIORSIGMA2);
  }

  //Sigma prior
  comm.Recv(&has_sigma_prior1,1,MPI::BOOL,src,pofd_mcmc::CLDSENDHASSIGMAPRIOR1);
  if (has_sigma_prior1) 
    comm.Recv(&sigma_prior_width1,1,MPI::DOUBLE,src,
	      pofd_mcmc::CLDSENDSIGMAPRIORWIDTH1);
  comm.Recv(&has_sigma_prior2,1,MPI::BOOL,src,pofd_mcmc::CLDSENDHASSIGMAPRIOR2);
  if (has_sigma_prior2) 
    comm.Recv(&sigma_prior_width2,1,MPI::DOUBLE,src,
	      pofd_mcmc::CLDSENDSIGMAPRIORWIDTH2);

}
