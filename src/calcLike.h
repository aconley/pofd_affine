//calcLike

#ifndef __calcLike__
#define __calcLike__

#include<string>
#include<vector>

#include<PDFactory.h>
#include<PD.h>
#include<fitsData.h>
#include<numberCountsKnotsSpline.h>
#include<beam.h>
#include<paramSet.h>
#include<utility.h>
#include<global_settings.h>

/*!
  \brief Class for computing likelihood of given model for multiple
  datasets that share the same beam.  Doesn't include priors.

  \ingroup Likelihoods

  Datasets are grouped by their beam because some of the compuation
  can be shared that way
*/

class calcLikeSingle {
 private:
  //Data set
  static const double flux_safety; //!< maxflux multiplier to ensure margin
  bool data_read; //!< Have we read data?
  unsigned int ndatasets; //!< Number of data sets
  fitsData* data; //!< Actual data sets -- len ndatasets
  double maxflux; //!< Maximum flux across all data sets for initPD call

  mutable PD pd; //!< Holds PD; convenience variable
  mutable PDFactory pdfac; //!< Computes PD

  //Noise information
  double* sigma_base; //!< Base value of sigma, len ndatasets
  double maxsigma_base; //!< Maximum value of sigma_base

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
  void readBeam(const std::string&, bool histogram=false, 
		double histogramlogstep=0.2);

  /*! \brief Access to beam info */
  const beam& getBeam() const { return bm; }

  void setVerbose() { verbose=true; } //!< Activate verbose mode
  void unSetVerbose() { verbose = false; } //!< Turn off verbose mode

  void applyBinning(unsigned int); //!< Bins data
  void removeBinning(); //!< Remove binning of data

  void setLikeNorm(unsigned int i, double val) { like_norm[i]=val;} //!< Set likelihood normalization factor relative to beam area; set to zero to nor normalize
  void setLikeNorm(const std::vector<double>&); //!< Set likelihood normalization factor
  void setLikeNorm(unsigned int n, const double* const); //!< Set likelihood normalization factor
  double getLikeNorm(unsigned int i) const { return like_norm[i]; } //!< Return likelihood normalization factor

  void setSigmaBase(const std::vector<double>&); //!< Set base sigma values
  void setSigmaBase(unsigned int, const double* const); //!< Set base sigma values
  /*! \brief Get sigma base value */
  double getSigmaBase(unsigned int i) const { return sigma_base[i]; }


  unsigned int getNDataSets() const { return ndatasets; } //!< Number of data sets
  unsigned int getNData(unsigned int i) const { return data[i].getN(); } //!< Number of data points in given data set

  /*! \brief Returns \f$\log L\f$ over all data sets.  Model must be set */
  double getLogLike(const numberCounts&, double sigmult=1.0, 
		    unsigned int fftsize=131072, bool edgefix=true) const;

  void writePDToStream( std::ostream& os ) const; //!< Write out computed P(D)
  
  /*! \brief MPI copy send operation */
  void SendSelf(MPI::Comm&, int dest) const;
  /*! \brief MPI copy recieve operation */
  void RecieveCopy(MPI::Comm&, int dest);
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
  \brief Class for computing likelihood of given model. Allows
  for multiple beams and includes priors.

  \ingroup Likelihoods
*/
class calcLike {
 private :
  //Transform stuff
  unsigned int fftsize; //!< Size of FFT (on each dim)
  unsigned int ninterp; //!< R interpolation size
  bool edgeFix; //!< Apply edge fix

  //Data
  unsigned int nbeamsets; //!< Number of beam sets 
  calcLikeSingle* beamsets; //!< Sets of data grouped by beam
  bool bin_data; //!< Bin the data
  unsigned int nbins; //!< Number of bins (if bin_data)

  //Priors
  bool has_cfirb_prior; //!< Are we using CFIRB prior
  double cfirb_prior_mean; //!< Value of CFIRB prior mean
  double cfirb_prior_sigma; //!< Value of CFIRB prior error
  bool has_sigma_prior; //!< Sigma multiplier value
  double sigma_prior_width; //!< Sigma multiplier width

  //Model
  mutable numberCountsKnotsSpline model; //!< Holds number counts model

  bool verbose; //!< Output informational messages while running

 public:
  /*! \brief Constructor */
  calcLike(unsigned int FFTSIZE=262144, unsigned int NINTERP=1024, 
	   bool EDGEFIX=true, bool BINNED=false, unsigned int NBINS=1000 );
  ~calcLike(); //!< Destuctor

  void addWisdom(const std::string& filename); //!< Add FFTW wisdom

  /*! \brief Reads data to compute likelihood for multiple data and beam files*/
  void readDataFromFiles(const std::vector<std::string>&, 
			 const std::vector<std::string>&,
			 const std::vector<double>&,
			 const std::vector<double>&,
			 bool IGNOREMASK=false, bool MEANSUB=false,
			 bool HISTOGRAM=false, double HISTOGRAMLOGSTEP=0.2);
  
  void setVerbose() { verbose=true; } //!< Turn on verbose mode
  void unSetVerbose() { verbose = false; } //!< Turn off verbose mode

  void setFFTSize(unsigned int val) { fftsize = val; } //!< Set FFT size
  unsigned int getFFTSize() const { return fftsize; } //!< Get FFT size

  void setNInterp(unsigned int); //!< Set interpolation size
  unsigned int getNInterp() const { return ninterp; } //!< Get interpolation size

  void setBinData(); //!< Turn on binning of data
  void unSetBinData(); //!< Turn off binning of data
  void setNBins(unsigned int); //!< Change number of bins

  void setEdgeFix() { edgeFix=true; } //!< Turn on edge fixing
  void unSetEdgeFix() { edgeFix=false; } //!< Turn off edge fixing
  bool getEdgeFix() const { return edgeFix; } //!< Are we edge fixing?

  /*! \brief Get number of knots in model */
  unsigned int getNKnots() const { return model.getNKnots(); } 
  /*! \brief Set positions of knots in model */
  void setKnotPositions(std::vector<double>& knts) 
  { model.setKnotPositions(knts); }
  /*! \brief Set positions of knots in model */
  void setKnotPositions(const initFileKnots& ifile)
  { ifile.getKnotPos(model); }

  //Sigma prior
  /*! \brief Activates the sigma prior with width set to value */
  void setSigmaPrior( double val ) { 
    has_sigma_prior=true; sigma_prior_width = val;}
  /*! \brief De-activate the sigma prior */
  void unsetSigmaPrior() { has_sigma_prior = false; }

  //CFIRB prior
  /*! \brief Activates the CFIRB prior with the specified values */
  void setCFIRBPrior( double, double );
  /*! \brief De-activated CFIRB prior */
  void unsetCFIRBPrior() { has_cfirb_prior = false; }
  
  /*! \brief Get number of beam sets */
  unsigned int getNBeamSets() const { return nbeamsets; }

  /*! \brief Get Log Likelihood for set of parameters */
  double getLogLike(const paramSet&) const;

  /*! \brief MPI copy send operation */
  void SendSelf(MPI::Comm&, int dest) const;
  /*! \brief MPI copy send operation */
  void RecieveCopy(MPI::Comm&, int dest);
};

#endif
