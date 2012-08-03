//calcLikeDouble

#ifndef __calcLikeDouble__
#define __calcLikeDouble__

#include<string>
#include<vector>

#include<PDFactoryDouble.h>
#include<PDDouble.h>
#include<fitsDataDouble.h>
#include<numberCountsDoubleLogNormal.h>
#include<doublebeam.h>
#include<paramSet.h>
#include<utility.h>
#include<global_settings.h>

/*!
  \brief Class for computing likelihood of given 2D model for multiple
  datasets that share the same beam.  Doesn't include priors.

  \ingroup Likelihoods

  Datasets are grouped by their beam because some of the compuation
  can be shared that way
*/

class calcLikeDoubleSingle { //Odd name...
 private:
  //Data set
  static const double flux_safety; //!< maxflux multiplier to ensure margin
  bool data_read; //!< Have we read data?
  unsigned int ndatasets; //!< Number of data sets
  fitsDataDouble* data; //!< Actual data sets -- len ndatasets
  double maxflux1; //!< Maximum flux across all data sets for initPD call, band1
  double maxflux2; //!< Maximum flux across all data sets for initPD call, band2

  mutable PDDouble pd; //!< Holds PD; convenience variable
  mutable PDFactoryDouble pdfac; //!< Computes PD

  //Noise information
  double* sigma_base1; //!< Base value of sigma, len ndatasets, band 1
  double* sigma_base2; //!< Base value of sigma, len ndatasets, band 2
  double maxsigma_base1; //!< Maximum value of sigma_base, band 1
  double maxsigma_base2; //!< Maximum value of sigma_base, band 2

  //Likelihood computation information
  double* like_norm; //!< Likelihood normalization factor for beam area
  double* like_offset; //!< Normalization constant for number of data points

  //Beam
  bool has_beam; //!< Is beam loaded?
  doublebeam bm; //!< Beam for this data

  bool verbose; //!< Output informational messages while running

  void resize(unsigned int); //!< Change number of data sets

 public :
  const static double bad_like; //!< Bad log likelihood value

  /*!\brief Default constructor */
  calcLikeDoubleSingle(unsigned int NEDGE=256); 
  ~calcLikeDoubleSingle(); //!< Destructor

  /*!\brief Add wisdom information */
  void addWisdom(const std::string& filename) { pdfac.addWisdom(filename); }

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
  void readBeam(const std::string&, const std::string&, bool histogram=false, 
		double histogramlogstep=0.2);

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


  unsigned int getNDataSets() const { return ndatasets; } //!< Number of data sets
  unsigned int getNData(unsigned int i) const { return data[i].getN(); } //!< Number of data points in given data set

  /*! \brief Returns \f$\log L\f$ over all data sets.  Model must be set */
  double getLogLike(const numberCountsDouble&, double sigmult1=1.0, 
		    double sigmult2=1.0, unsigned int fftsize=4096, 
		    bool edgefix=true, bool edgeinteg=true) const;

  void writePDToStream( std::ostream& os ) const; //!< Write out computed P(D)

  /*! \brief MPI copy send operation */
  void SendSelf(MPI::Comm&, int dest) const;
  /*! \brief MPI copy recieve operation */
  void RecieveCopy(MPI::Comm&, int dest);
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
  \brief Class for computing likelihood of given model. Allows
  for multiple beams and includes priors.

  \ingroup Likelihoods
*/
class calcLikeDouble {
 private :
  //Transform stuff
  unsigned int fftsize; //!< Size of FFT (on each dim)
  unsigned int nedge; //!< R edge size
  bool edgeInteg; //!< Do edge integration
  bool edgeFix; //!< Apply edge fix

  //Data
  unsigned int nbeamsets; //!< Number of beam sets
  calcLikeDoubleSingle* beamsets; //!< Sets of data grouped by beam
  bool bin_data; //!< Bin the data
  unsigned int nbins; //!< Number of bins (if bin_data)

  //Priors
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

  //Model
  mutable numberCountsDoubleLogNormal model; //!< Number counts model

  bool verbose; //!< Output informational messages while running

 public:
  /*! \brief Constructor */
  calcLikeDouble(unsigned int FFTSIZE=4096, unsigned int NEDGE=256, 
		 bool EDGEFIX=true, bool EDGEINTEG=true,
		 bool BINNED=false, unsigned int NBINS=1000);
  ~calcLikeDouble(); //!< Destructor 

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
			 bool HISTOGRAM=false, double HISTOGRAMLOGSTEP=0.2);
  
  void setVerbose() { verbose=true; } //!< Turn on verbose mode
  void unSetVerbose() { verbose = false; } //!< Turn off verbose mode

  void setFFTSize(unsigned int val) { fftsize = val; } //!< Set FFT size
  unsigned int getFFTSize() const { return fftsize; } //!< Return current FFT size

  void setEdgeFix() { edgeFix=true; } //!< Turn on edge fixing
  void unSetEdgeFix() { edgeFix=false; } //!< Turn off edge fixing
  bool getEdgeFix() const { return edgeFix; } //!< Are we using edge fixing?

  //Note -- you can't set the edge integration size except at construction
  void setEdgeInteg() { edgeInteg=true; } //!< Turn on edge integration
  void unSetEdgeInteg() { edgeInteg=false; } //!< Turn off edge integration
  bool getEdgeInteg() const { return edgeInteg; } //!< Are we doing edge integration?

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
    model.setPositions(K,S,O);
  }
  /*! \brief Set positions of all knots */
  void setPositions(const initFileDoubleLogNormal&);

  //Sigma prior
  /*! \brief Activates the sigma prior with width set to value, band 1*/
  void setSigmaPrior1( double val ) { 
    has_sigma_prior1 = true; sigma_prior_width1 = val;}
  /*! \brief De-activate the sigma prior, band 1*/
  void unsetSigmaPrior1() { has_sigma_prior1 = false; }
  /*! \brief Activates the sigma prior with width set to value, band 2*/
  void setSigmaPrior2( double val ) { 
    has_sigma_prior2 = true; sigma_prior_width2 = val;}
  /*! \brief De-activate the sigma prior, band 2*/
  void unsetSigmaPrior2() { has_sigma_prior2 = false; }

  //CFIRB prior
  /*! \brief Activates the CFIRB prior with the specified values, band 1 */
  void setCFIRBPrior1( double, double );
  /*! \brief De-activated CFIRB prior, band 1 */
  void unsetCFIRBPrior() { has_cfirb_prior1 = false; }
  /*! \brief Activates the CFIRB prior with the specified values, band 2 */
  void setCFIRBPrior2( double, double );
  /*! \brief De-activated CFIRB prior, band 2 */
  void unsetCFIRBPrior2() { has_cfirb_prior2 = false; }

  /*! \brief Get number of beam sets */
  unsigned int getNBeamSets() const { return nbeamsets; }

  /*! \brief Get Log-Likelihood of data for a set of parameters */
  double getLogLike(const paramSet&) const;

  /*! \brief MPI copy send operation */
  void SendSelf(MPI::Comm&, int dest) const;
  /*! \brief MPI copy recieve operation */
  void RecieveCopy(MPI::Comm&, int dest);
};

//////////////////////////////////////////////////
/*!
  \brief Structure to read in spec files for fits

  Lines in the spec file control the likelihood calculations and fits.
  The format is: identifier= values
  The appropriate entries for values depends upon identifier=

  The valid values for identifier= and their associated values are:
   dataset= datafile1 datafile2 inst_sigma1 inst_sigma2 psffile1 \
             psffile2 [like_norm]
               datafile1 and datafile2 are the names of fits data files,
	        in band 1 and band 2, respectively
	       inst_sigma1, inst_sigma2 are the instrumental noise (Gaussian) 
	        in the same units as the data in datafile (e.g., Jy)
		for datafile1 and datafile2
	       psffile1, psffile2 are FITS files giving the beam in each band
	       like_norm is an optional input controlling how the
	        likelihood is normalized for the beam area (def: 1)
   sigmaprior1=  stdev
                  This activates the Gaussian prior on the instrumental noise
		  multiplier in band 1.  It has mean one and standard 
		  deviation stdev
   sigmaprior2=  stdev
                  Same as sigmaprior1, but in band 2
   cfirbprior1=  mean stdev
                  This activates a Gaussian prior on the integrated CFIRB 
		  luminoisty in band 1. It has a mean of mean and a standard 
                  deviation stdev expressed in the same flux and area units as 
		  the model.  So if your model dN/dS is in Jy^-1 deg^-2,
		  then your CFIRB prior should be in Jy deg^-2
   cfirbprior2= mean stdev
                  Same as cfirbprior1 but in band 2
  
  Having multiple lines with dataset= is additive.
  Having multiple lines with sigmaprior= or cfirbprior= results in only
  the values from the final such line being used.

  Including a space after= is important.
*/
struct specFileDouble {
  //dataset
  std::vector<std::string> datafiles1; //!< Data files, band 1
  std::vector<std::string> datafiles2; //!< Data files, band 2
  std::vector<double> sigmas1; //!< Instrumental noise values, band 1
  std::vector<double> sigmas2; //!< Instrumental noise values, band 2
  std::vector<std::string> psffiles1; //!< Beam file, band 1
  std::vector<std::string> psffiles2; //!< Beam file, band 2
  std::vector<double> like_norm; //!< Likelihood normalization

  
  //sigma prior
  bool has_sigprior1; //!< Is the sigma prior on, band 1
  double sigprior_stdev1; //!< Stdev of sigma prior, band 1
  bool has_sigprior2; //!< Is the sigma prior on, band 2
  double sigprior_stdev2; //!< Stdev of sigma prior, band 2

  //cfirb_prior
  bool has_cfirbprior1; //!< Is cfirb prior on, band 1
  double cfirbprior_mean1; //!< Mean value of prior, band 1
  double cfirbprior_stdev1; //!< Stdev of prior, band 1
  bool has_cfirbprior2; //!< Is cfirb prior on, band 2
  double cfirbprior_mean2; //!< Mean value of prior, band 2
  double cfirbprior_stdev2; //!< Stdev of prior, band 2

  specFileDouble(); //!< Default constructor
  specFileDouble(const std::string&); //!< Constructor with file read

  void readFile(const std::string&); //!< Read in file
};



#endif
