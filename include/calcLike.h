//calcLike

#ifndef __calcLike__
#define __calcLike__

#include<string>
#include<vector>

#include "../include/PDFactory.h"
#include "../include/PD.h"
#include "../include/fitsData.h"
#include "../include/numberCountsKnotsSpline.h"
#include "../include/beam.h"
#include "../include/paramSet.h"
#include "../include/utility.h"
#include "../include/global_settings.h"

/*!
  \brief Class for computing likelihood of given model for multiple
  datasets that share the same beam, 1D case.  

  Doesn't include priors.

  Datasets are grouped by their beam because some of the computation
  can be shared that way.

  \ingroup Likelihoods  
*/

class calcLikeSingle {
 private:
  //Data set
  static const double flux_safety; //!< maxflux multiplier to ensure margin
  bool data_read; //!< Have we read data?
  unsigned int ndatasets; //!< Number of data sets
  fitsData* data; //!< Actual data sets -- len ndatasets
  double minDataFlux; //!< Minimum flux density in actual data 
  double maxDataFlux; //!< Maximum flux density in actual data 

  mutable PD pd; //!< Holds PD; convenience variable
  mutable PDFactory pdfac; //!< Computes PD
  double minRFlux; //!< Minimum flux for initPD call
  double maxRFlux; //!< Maximum flux for initPD call

  //Noise information
  double* sigma_base; //!< Base value of sigma, len ndatasets
  double maxsigma_base; //!< Maximum value of sigma_base
  double exp_conf; //!< Expected confusion noise, in Jy

  //Likelihood computation information
  double* like_norm; //!< Likelihood normalization factor for beam area
  double* like_offset; //!< Normalization constant for number of data points

  //Beam
  bool has_beam; //!< Is beam loaded?
  beam bm; //!< Beam for this data

  bool verbose; //!< Output informational messages while running

  void resize(unsigned int); //!< Change number of data sets

 public :
  const static double bad_like; //!< Bad log likelihood value

  /*!\brief Default constructor */
  calcLikeSingle(unsigned int NINTERP=1024); 
  ~calcLikeSingle(); //!< Destructor

  void free(); //!< Frees all memory

  /*!\brief Add wisdom information */
  void addWisdom(const std::string& filename) { pdfac.addWisdom(filename); }

  void setNInterp(unsigned int n) { pdfac.setNInterp(n); } //!< Set interpolation length
  unsigned int getNInterp() const { return pdfac.getNInterp(); } //!< Get interpolation length

  /*! \brief Reads data to compute likelihood for a single file*/
  void readDataFromFile(const std::string&, bool IGNOREMASK=false,
			bool MEANSUB=false,bool BINDATA=false,
			unsigned int NBINS=1000);

  /*! \brief Reads data to compute likelihood for multiple data files*/
  void readDataFromFiles(const std::vector<std::string>&, 
			 bool IGNOREMASK=false, bool MEANSUB=false,
			 bool BINDATA=false, unsigned int NBINS=1000);

  /*! \brief Read in beam file */
  void readBeam(const std::string&, double MINVAL=1e-5, bool histogram=false,
		unsigned int NBINS=120);

  /*! \brief Access to beam info */
  const beam& getBeam() const { return bm; }

  void setVerbose() { verbose=true; } //!< Activate verbose mode
  void unSetVerbose() { verbose = false; } //!< Turn off verbose mode

  void applyBinning(unsigned int); //!< Bins data
  void removeBinning(); //!< Remove binning of data

  void setRRange(const numberCounts&); //!< Set range of R evaluations

  void setLikeNorm(unsigned int i, double val) { like_norm[i]=val;} //!< Set likelihood normalization factor relative to beam area; set to zero to nor normalize
  void setLikeNorm(const std::vector<double>&); //!< Set likelihood normalization factor
  void setLikeNorm(unsigned int n, const double* const); //!< Set likelihood normalization factor
  double getLikeNorm(unsigned int i) const { return like_norm[i]; } //!< Return likelihood normalization factor

  void setSigmaBase(const std::vector<double>&); //!< Set base sigma values
  void setSigmaBase(unsigned int, const double* const); //!< Set base sigma values
  /*! \brief Get sigma base value */
  double getSigmaBase(unsigned int i) const { return sigma_base[i]; }

  void setExpConf(double v) { exp_conf = v;} //!< Set expected confusion noise
  double getExpConf() const { return exp_conf;} //!< Get expected confusion noise

  unsigned int getNDataSets() const { return ndatasets; } //!< Number of data sets
  unsigned int getNData(unsigned int i) const { return data[i].getN(); } //!< Number of data points in given data set

  /*! \brief Returns \f$\log L\f$ over all data sets.  Model must be set */
  double getLogLike(const numberCounts&, bool& pars_invalid, 
		    double sigmult=1.0, unsigned int fftsize=131072) const;

  void writePDToStream(std::ostream& os) const; //!< Write out computed P(D)
  
  /*! \brief MPI copy send operation */
  void sendSelf(MPI_Comm, int dest) const;
  /*! \brief MPI copy recieve operation */
  void recieveCopy(MPI_Comm, int dest);
};

////////////////////////////////////////////////////
/*!
  \brief Utility structure for grouping things by beam
*/
struct beam_group {
  unsigned int n; //!< Number in group
  std::vector<std::string> datafiles; //!< List of datafiles (len n)
  std::string beamfile; //!< Beam file, shared between all members
  std::vector<double> sigmas; //!< List of instrumental sigmas (len n)
  std::vector<double> like_norms; //!< List of likelihood normalizations (len n)

  beam_group() { n = 0; } //!< Constructor
};


////////////////////////////////////////////////////

/*!
  \brief Class for computing likelihood of given model, 1D case. 

  Allows for multiple beams and includes priors.

  \ingroup Likelihoods
*/
class calcLike {
 private :
  //Transform stuff
  unsigned int fftsize; //!< Size of FFT (on each dim)
  unsigned int ninterp; //!< R interpolation size

  // Data
  unsigned int nbeamsets; //!< Number of beam sets 
  calcLikeSingle* beamsets; //!< Sets of data grouped by beam
  bool bin_data; //!< Bin the data
  unsigned int nbins; //!< Number of bins (if bin_data)

  // Mean flux per area from model and flux^2; used in cfirb 
  //  and poisson priors, but also as bonus parameters
  mutable paramSet meanParams; //!< Mean fluxes computed using this
  mutable double mean_flux_per_area; //!< Mean flux per area
  mutable double mean_fluxsq_per_area; //!< Mean flux^2 per area

  //Priors
  bool has_cfirb_prior; //!< Are we using CFIRB prior
  double cfirb_prior_mean; //!< Value of CFIRB prior mean
  double cfirb_prior_sigma; //!< Value of CFIRB prior error
  bool has_sigma_prior; //!< Sigma multiplier value
  double sigma_prior_width; //!< Sigma multiplier width
  bool has_poisson_prior; //!< Do we have prior on Poisson noise?
  double poisson_prior_mean; //!< Value of Poisson prior mean
  double poisson_prior_sigma; //!< Value of Poisson prior sigma

  // Regularization
  double regularization_alpha; //!< Regularization multiplier

  // Model
  mutable numberCountsKnotsSpline model; //!< Holds number counts model

  bool verbose; //!< Output informational messages while running

 public:
  /*! \brief Constructor */
  calcLike(unsigned int FFTSIZE=262144, unsigned int NINTERP=1024, 
	   bool BINNED=false, unsigned int NBINS=1000);
  ~calcLike(); //!< Destuctor

  void freeData(); //!< Remove internal data

  void addWisdom(const std::string& filename); //!< Add FFTW wisdom

  /*! \brief Reads data to compute likelihood for multiple data and beam files*/
  void readDataFromFiles(const std::vector<std::string>&, 
			 const std::vector<std::string>&,
			 const std::vector<double>&,
			 const std::vector<double>&,
			 bool IGNOREMASK=false, bool MEANSUB=false,
			 double MINBEAMVAL=1e-5, bool HISTOGRAMBEAMS=false, 
			 unsigned int NBEAMHIST=120, double EXPCONF=0.0);
  
  void setVerbose() { verbose=true; } //!< Turn on verbose mode
  void unSetVerbose() { verbose = false; } //!< Turn off verbose mode

  void setFFTSize(unsigned int val) { fftsize = val; } //!< Set FFT size
  unsigned int getFFTSize() const { return fftsize; } //!< Get FFT size

  void setNInterp(unsigned int); //!< Set interpolation size
  unsigned int getNInterp() const { return ninterp; } //!< Get interpolation size
  void setBinData(); //!< Turn on binning of data
  void unSetBinData(); //!< Turn off binning of data
  void setNBins(unsigned int); //!< Change number of bins

  void setRRanges(const paramSet& p); //!< Set R ranges for all datasets

  /*! \brief Get number of knots in model */
  unsigned int getNKnots() const { return model.getNKnots(); } 
  /*! \brief Set positions of knots in model */
  void setKnotPositions(std::vector<double>& knts) 
  { model.setKnotPositions(knts); }
  /*! \brief Set positions of knots in model */
  void setKnotPositions(const initFileKnots& ifile)
  { ifile.getKnotPos(model); }

  // Sigma prior
  /*! \brief Activates the sigma prior with width set to value */
  void setSigmaPrior(double val) { 
    has_sigma_prior=true; sigma_prior_width = val;}
  /*! \brief De-activate the sigma prior */
  void unsetSigmaPrior() { has_sigma_prior = false; }

  // CFIRB prior
  /*! \brief Activates the CFIRB prior with the specified values */
  void setCFIRBPrior(double, double);
  /*! \brief De-activate CFIRB prior */
  void unsetCFIRBPrior() { has_cfirb_prior = false; }
  
  // Poisson prior
  /*! \brief Activates the Poisson prior with the specified values */
  void setPoissonPrior(double, double);
  /*! \brief De-activate Poisson prior */
  void unsetPoissonPrior() { has_poisson_prior = false; }

  /*! \brief Set Regularization Alpha*/
  void setRegularizationAlpha(double);

  /*! \brief Get number of beam sets */
  unsigned int getNBeamSets() const { return nbeamsets; }

  /*! \brief Get Log Likelihood for set of parameters */
  double getLogLike(const paramSet&, bool&) const;

  /*! \brief Fills bonus params (mean flux per area values) */
  void fillBonusParams(paramSet& p) const;

  /*! \brief Write to HDF5 handle */
  void writeToHDF5Handle(hid_t) const;

  /*! \brief MPI copy send operation */
  void sendSelf(MPI_Comm, int dest) const;
  /*! \brief MPI copy send operation */
  void recieveCopy(MPI_Comm, int dest);
};

#endif
