//calcLikeDouble

#ifndef __calcLikeDouble__
#define __calcLikeDouble__

#include<string>
#include<vector>

#include "../include/PDFactoryDouble.h"
#include "../include/PDDouble.h"
#include "../include/fitsDataDouble.h"
#include "../include/numberCountsDoubleLogNormal.h"
#include "../include/doublebeam.h"
#include "../include/paramSet.h"
#include "../include/utility.h"
#include "../include/global_settings.h"

/*!
  \brief Class for computing likelihood of given model for multiple
  datasets that share the same beam, 2D case.  

  Doesn't include priors.

  Datasets are grouped by their beam because some of the compuation
  can be shared that way.

  \ingroup Likelihoods
*/

class calcLikeDoubleSingle { //Odd name...
 private:
  //Data set.
  static const double flux_safety; //!< maxflux multiplier to ensure margin
  bool data_read; //!< Have we read data?
  unsigned int ndatasets; //!< Number of data sets
  std::string* filenames1; //!< Names of files, band 1
  std::string* filenames2; //!< Names of files, band 2
  fitsDataDouble* data; //!< Actual data sets -- len ndatasets
  double minDataFlux1; //!< Minimum flux density in actual data, band 1
  double maxDataFlux1; //!< Maximum flux density in actual data, band 1
  double minDataFlux2; //!< Minimum flux density in actual data, band 2
  double maxDataFlux2; //!< Maximum flux density in actual data, band 2
  unsigned int* dataext1; //!< Data extension, band 1
  unsigned int* dataext2; //!< Data extension, band 2
  bool *hasmask1; //!< Had mask, band 1
  bool *hasmask2; //!< Had mask, band 2
  unsigned int* maskext1; //!< Mask extension, band 1
  unsigned int* maskext2; //!< Mask extension, band 2
  
  mutable PDDouble pd; //!< Holds PD; convenience variable
  mutable PDFactoryDouble pdfac; //!< Computes PD
  double minRFlux1; //!< Minimum flux for initPD call, band 1
  double maxRFlux1; //!< Maximum flux for initPD call, band 1
  double minRFlux2; //!< Minimum flux for initPD call, band 2
  double maxRFlux2; //!< Maximum flux for initPD call, band 2

  //Noise information
  double* sigma_base1; //!< Base value of sigma, len ndatasets, band 1
  double* sigma_base2; //!< Base value of sigma, len ndatasets, band 2
  double maxsigma_base1; //!< Maximum value of sigma_base, band 1
  double maxsigma_base2; //!< Maximum value of sigma_base, band 2
  double exp_conf1; //!< Expected confusion noise, band 1
  double exp_conf2; //!< Expected confusion noise, band 2

  //Likelihood computation information
  double* like_norm; //!< Likelihood normalization factor for beam area
  double* like_offset; //!< Normalization constant for number of data points

  //Beam
  bool has_beam; //!< Is beam loaded?
  std::string beamfile1; //!< Filename for band 1 beam
  std::string beamfile2; //!< Filename for band 2 beam
  doublebeam bm; //!< Beam for this data

  bool verbose; //!< Output informational messages while running

  void resize(unsigned int); //!< Change number of data sets

 public :
  const static double bad_like; //!< Bad log likelihood value

  /*!\brief Default constructor */
  calcLikeDoubleSingle(unsigned int NEDGE=256); 
  ~calcLikeDoubleSingle(); //!< Destructor

  void free(); //!< Frees memory for large structures, keeps some summary stats

  /*!\brief Add wisdom information */
  void addWisdom(const std::string& filename) { pdfac.addWisdom(filename); }

  void setNEdge(unsigned int n) { pdfac.setNEdge(n); } //!< Set interpolation length
  unsigned int getNEdge() const { return pdfac.getNEdge(); } //!< Get interpolation length

  void setRRange(const numberCountsDouble&); //!< Set range of R evaluations

  void applyBinning(unsigned int); //!< Bin data
  void removeBinning(); //!< Remove binning

  /*! \brief Reads data to compute likelihood for a single file*/
  void readDataFromFile(const std::string&, const std::string&,
			bool IGNOREMASK=false,
			bool MEANSUB=false,bool BINDATA=false,
			unsigned int NBINS=100);

  /*! \brief Reads data to compute likelihood for multiple data files*/
  void readDataFromFiles(const std::vector<std::string>&, 
			 const std::vector<std::string>&, 
			 bool IGNOREMASK=false, bool MEANSUB=false,
			 bool BINDATA=false, unsigned int NBINS=100);

  /*! \brief Read in beam files */
  void readBeam(const std::string&, const std::string&, 
		double MINVAL=1e-6, bool histogram=false, 
		unsigned int NBINS=150);

  /*! \brief Access to beam info */
  const doublebeam& getBeam() const { return bm; }

  void setVerbose() { verbose=true; } //!< Activate verbose mode
  void unSetVerbose() { verbose = false; } //!< Turn off verbose mode

  void setLikeNorm(unsigned int i, double val) { like_norm[i]=val;} //!< Set likelihood normalization factor relative to beam area; set to zero to nor normalize
  void setLikeNorm(const std::vector<double>&); //!< Set likelihood normalization factor
  void setLikeNorm(unsigned int n, const double* const); //!< Set likelihood normalization factor
  double getLikeNorm(unsigned int i) const { return like_norm[i]; } //!< Return likelihood normalization factor

  void setSigmaBase1(const std::vector<double>&); //!< Set base sigma values, band 1
  void setSigmaBase1(unsigned int, const double* const); //!< Set base sigma values, band 1
  /*! \brief Gets sigma base value for band 1 */
  double getSigmaBase1(unsigned int i) const { return sigma_base1[i]; }
  void setSigmaBase2(const std::vector<double>&); //!< Set base sigma values, band 2
  void setSigmaBase2(unsigned int, const double* const); //!< Set base sigma values, band 2
  /*! \brief Gets sigma base value for band 2 */
  double getSigmaBase2(unsigned int i) const { return sigma_base2[i]; }

  void setExpConf1(double v) { exp_conf1 = v;} //!< Set expected confusion noise, band 1
  double getExpConf1() const { return exp_conf1;} //!< Get expected confusion noise, band 1
  void setExpConf2(double v) { exp_conf2 = v;} //!< Set expected confusion noise, band 1
  double getExpConf2() const { return exp_conf2;} //!< Get expected confusion noise, band 1

  unsigned int getNDataSets() const { return ndatasets; } //!< Number of data sets
  unsigned int getNData(unsigned int i) const { return data[i].getN(); } //!< Number of data points in given data set

  /*! \brief Returns \f$\log L\f$ over all data sets.  Model must be set */
  double getLogLike(const numberCountsDouble&, bool&, double sigmult1=1.0, 
		    double sigmult2=1.0, unsigned int fftsize=4096, 
		    bool edgeinteg=true) const;

  void writeToHDF5Handle(hid_t) const; //!< Write dataset info to HDF5 handle
  /*! \brief Make a new group and write to it */
  void writeToNewHDF5Group(hid_t, const std::string& groupname) const;

  /*! \brief MPI copy send operation */
  void sendSelf(MPI_Comm, int dest) const;
  /*! \brief MPI copy receive operation */
  void receiveCopy(MPI_Comm, int dest);
};

////////////////////////////////////////////////////
/*! \brief Utility structure for grouping things by beam */
struct doublebeam_group {
  unsigned int n; //!< Number in group
  std::vector<std::string> datafiles1; //!< Band 1 datafiles (len n)
  std::vector<std::string> datafiles2; //!< Band 2 datafiles (len n)
  std::string beamfile1; //!< Beamfile, band 1, for all elements in group
  std::string beamfile2; //!< Beamfile, band 2, for all elements in group
  std::vector<double> sigmas1; //!< Instrument noise, band 1 (len n)
  std::vector<double> sigmas2; //!< Instrument noise, band 2 (len n)
  std::vector<double> like_norms; //!< Likelihood normalization (len n)
 
  doublebeam_group() { n = 0; } //!< Constructor
};


////////////////////////////////////////////////////

/*!
  \brief Class for computing likelihood of given model, 2D case. 

  Allows for multiple beams and includes priors.

  \ingroup Likelihoods
*/
class calcLikeDouble {
 private :
  // Transform stuff
  unsigned int fftsize; //!< Size of FFT (on each dim)
  unsigned int nedge; //!< R edge size
  bool edgeInteg; //!< Do edge integration

  // Data
  unsigned int nbeamsets; //!< Number of beam sets
  calcLikeDoubleSingle* beamsets; //!< Sets of data grouped by beam
  bool bin_data; //!< Bin the data
  unsigned int nbins; //!< Number of bins (if bin_data)

  // Mean flux per area from model and flux^2; used in cfirb 
  //  and poisson priors, but also as bonus parameters
  mutable paramSet meanParams; //!< params used to compute mean flux values
  mutable double mean_flux_per_area1; //!< Mean flux per area, band 1
  mutable double mean_fluxsq_per_area1; //!< Mean flux^2 per area, band 1
  mutable double mean_flux_per_area2; //!< Mean flux per area, band 2
  mutable double mean_fluxsq_per_area2; //!< Mean flux^2 per area, band 2

  // Priors
  bool has_cfirb_prior1; //!< Are we using CFIRB prior, band 1
  double cfirb_prior_mean1; //!< Value of CFIRB prior mean, band 1
  double cfirb_prior_sigma1; //!< Value of CFIRB prior error, band 1
  bool has_cfirb_prior2; //!< Are we using CFIRB prior, band 2
  double cfirb_prior_mean2; //!< Value of CFIRB prior mean, band 2
  double cfirb_prior_sigma2; //!< Value of CFIRB prior error, band 2
  bool has_sigma_prior1; //!< Sigma multiplier value, band 1
  double sigma_prior_width1; //!< Sigma multiplier width, band 1
  bool has_sigma_prior2; //!< Sigma multiplier value, band 2
  double sigma_prior_width2; //!< Sigma multiplier width, band 2
  bool has_poisson_prior1; //!< Do we have prior on Poisson noise, band 1?
  double poisson_prior_mean1; //!< Value of Poisson prior mean, band 1
  double poisson_prior_sigma1; //!< Value of Poisson prior sigma, band 1
  bool has_poisson_prior2; //!< Do we have prior on Poisson noise, band 2?
  double poisson_prior_mean2; //!< Value of Poisson prior mean, band 2
  double poisson_prior_sigma2; //!< Value of Poisson prior sigma, band 2

  // Regularization
  double regularization_alpha; //!< Regularization multiplier

  // Model
  mutable numberCountsDoubleLogNormal model; //!< Number counts model

  bool verbose; //!< Output informational messages while running

 public:
  /*! \brief Constructor */
  calcLikeDouble(unsigned int FFTSIZE=4096, unsigned int NEDGE=256, 
		 bool EDGEINTEG=true, bool BINNED=false, 
		 unsigned int NBINS=1000);
  ~calcLikeDouble(); //!< Destructor 

  void freeData(); //!< Remove internal data

  void addWisdom(const std::string& filename); //!< Add FFTW wisdom information

  /*! \brief Reads data to compute likelihood for multiple data and beam files*/
  void readDataFromFiles(const std::vector<std::string>&, 
			 const std::vector<std::string>&,
			 const std::vector<std::string>&, 
			 const std::vector<std::string>&,
			 const std::vector<double>&,
			 const std::vector<double>&,
			 const std::vector<double>&,
			 bool IGNOREMASK=false, bool MEANSUB=false,
			 double MINBEAMVAL=1e-6, bool HISTOGRAM=false,
			 unsigned int NBEAMHIST=160,
			 double EXPCONF1=0.0, double EXPCONF2=0.0);
  
  void setVerbose() { verbose=true; } //!< Turn on verbose mode
  void unSetVerbose() { verbose = false; } //!< Turn off verbose mode

  void setFFTSize(unsigned int val) { fftsize = val; } //!< Set FFT size
  unsigned int getFFTSize() const { return fftsize; } //!< Return current FFT size
  void setEdgeInteg() { edgeInteg=true; } //!< Turn on edge integration
  void unSetEdgeInteg() { edgeInteg=false; } //!< Turn off edge integration
  bool getEdgeInteg() const { return edgeInteg; } //!< Are we doing edge integration?
  void setNEdge(unsigned int); //!< Set edge integration size
  unsigned int getNEdge() const { return nedge; } //!< Get edge integration size

  void setRRanges(const paramSet& p); //!< Set R ranges for all datasets

  void setBinData(); //!< Turn on binning of data
  void unSetBinData(); //!< Turn off binning of data
  void setNBins(unsigned int); //!< Change number of bins

  /*! \brief Return number of knots for 1st band spline*/
  unsigned int getNKnots() const { return model.getNKnots(); }
  /*! \brief Return number of knots in color sigma model */
  unsigned int getNSigmas() const { return model.getNSigmas(); }
  /*! \brief Return number of knots in color offset model */
  unsigned int getNOffsets() const { return model.getNOffsets(); }
  /*! \brief Set knot positoins for 1st band spline */
  void setKnotPositions(std::vector<double>& knts) 
  { model.setKnotPositions(knts); }
  /*! \brief Set knot positions for color sigma model */
  void setSigmaPositions(std::vector<double>& knts) 
  { model.setSigmaPositions(knts); }
  /*! \brief Set knot positions for color offset model */
  void setOffsetPositions(std::vector<double>& knts) 
  { model.setOffsetPositions(knts); }
  /*! \brief Set positions of all knots */
  void setPositions(std::vector<double>& K, std::vector<double>& S,
		    std::vector<double>& O) {
    model.setPositions(K, S, O);
  }
  /*! \brief Set positions of all knots */
  void setPositions(const initFileDoubleLogNormal&);

  // Sigma priors
  /*! \brief Activates the sigma prior with width set to value, band 1*/
  void setSigmaPrior1(double val) { 
    has_sigma_prior1 = true; sigma_prior_width1 = val;}
  /*! \brief De-activate the sigma prior, band 1*/
  void unsetSigmaPrior1() { has_sigma_prior1 = false; }
  /*! \brief Activates the sigma prior with width set to value, band 2*/
  void setSigmaPrior2(double val) { 
    has_sigma_prior2 = true; sigma_prior_width2 = val;}
  /*! \brief De-activate the sigma prior, band 2*/
  void unsetSigmaPrior2() { has_sigma_prior2 = false; }

  // CFIRB priors
  /*! \brief Activates the CFIRB prior with the specified values, band 1 */
  void setCFIRBPrior1(double, double);
  /*! \brief De-activate CFIRB prior, band 1 */
  void unsetCFIRBPrior1() { has_cfirb_prior1 = false; }
  /*! \brief Activates the CFIRB prior with the specified values, band 2 */
  void setCFIRBPrior2(double, double);
  /*! \brief De-activate CFIRB prior, band 2 */
  void unsetCFIRBPrior2() { has_cfirb_prior2 = false; }

  // Poisson priors
  /*! \brief Activates the Poisson prior with the specified values, band 1 */
  void setPoissonPrior1( double, double);
  /*! \brief De-activates Poisson prior, band 1 */
  void unsetPoissonPrior1() { has_poisson_prior1 = false; }
  /*! \brief Activates the Poisson prior with the specified values, band 2 */
  void setPoissonPrior2(double, double);
  /*! \brief De-activates Poisson prior, band 2 */
  void unsetPoissonPrior2() { has_poisson_prior2 = false; }

  /*! \brief Set Regularization Alpha*/
  void setRegularizationAlpha(double);

  /*! \brief Get number of beam sets */
  unsigned int getNBeamSets() const { return nbeamsets; }

  /*! \brief Get Log-Likelihood of data for a set of parameters */
  double getLogLike(const paramSet&, bool&) const;

  /*! \brief Fill bonus params (mean fluxes) */
  void fillBonusParams(paramSet& par) const;

  /*! \brief Write to HDF5 handle */
  void writeToHDF5Handle(hid_t) const;

  /*! \brief MPI copy send operation */
  void sendSelf(MPI_Comm, int dest) const;
  /*! \brief MPI copy receive operation */
  void receiveCopy(MPI_Comm, int dest);
};

#endif
