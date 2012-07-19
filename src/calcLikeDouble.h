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

  void setVerbose() { verbose=true; }
  void unSetVerbose() { verbose = false; }

  void setLikeNorm(unsigned int i, double val) { like_norm[i]=val;} //!< Set likelihood normalization factor relative to beam area; set to zero to nor normalize
  void setLikeNorm(const std::vector<double>&); //!< Set likelihood normalization factor
  void setLikeNorm(unsigned int n, const double* const); //!< Set likelihood normalization factor
  double getLikeNorm(unsigned int i) const { return like_norm[i]; } //!< Return likelihood normalization factor

  void setSigmaBase1(const std::vector<double>&); //!< Set base sigma values, band 1
  void setSigmaBase1(unsigned int, const double* const); //!< Set base sigma values, band 1
  double getSigmaBase1(unsigned int i) const { return sigma_base1[i]; }
  void setSigmaBase2(const std::vector<double>&); //!< Set base sigma values, band 2
  void setSigmaBase2(unsigned int, const double* const); //!< Set base sigma values, band 2
  double getSigmaBase2(unsigned int i) const { return sigma_base2[i]; }


  unsigned int getNDataSets() const { return ndatasets; } //!< Number of data sets
  unsigned int getNData(unsigned int i) const { return data[i].getN(); } //!< Number of data points in given data set

  /*! \brief Returns \f$\log L\f$ over all data sets.  Model must be set */
  double getLogLike(const numberCountsDouble&, double sigmult=1.0, 
		    unsigned int fftsize=4096, bool edgefix=true,
		    bool setedge=true) const;

  void writePDToStream( std::ostream& os ) const; //!< Write out computed P(D)

  void SendSelf(MPI::Comm&, int dest) const;
  void RecieveCopy(MPI::Comm&, int dest);
};

////////////////////////////////////////////////////
// Utility structure for grouping things by beam
struct doublebeam_group {
  unsigned int n;
  std::vector<std::string> datafiles1;
  std::vector<std::string> datafiles2;
  std::string beamfile1;
  std::string beamfile2;
  std::vector<double> sigmas1;
  std::vector<double> sigmas2;
  std::vector<double> like_norms;

  doublebeam_group() { n = 0; }
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
  bool edgeInterp; //!< Do edge interpolation
  bool edgeFix; //!< Apply edge fix

  //Data
  unsigned int nbeamsets;
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
  bool has_sigma_prior; //!< Sigma multiplier value, both bands
  double sigma_prior_width; //!< Sigma multiplier width, both bands

  //Model
  mutable numberCountsDoubleLogNormal model;

  bool verbose; //!< Output informational messages while running

 public:
  calcLikeDouble(unsigned int FFTSIZE=4096, unsigned int NEDGE=256, 
		 bool EDGEFIX=true, bool EDGEINTERP=true,
		 bool BINNED=false, unsigned int NBINS=1000);
  ~calcLikeDouble();

  void addWisdom(const std::string& filename);

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
  
  void setVerbose() { verbose=true; }
  void unSetVerbose() { verbose = false; }

  void setFFTSize(unsigned int val) { fftsize = val; }
  unsigned int getFFTSize() const { return fftsize; }

  void setEdgeFix() { edgeFix=true; }
  void unSetEdgeFix() { edgeFix=false; }
  bool getEdgeFix() const { return edgeFix; }

  //Note -- you can't set the edge interp size except at construction
  void setEdgeInterp() { edgeInterp=true; }
  void unSetEdgeInterp() { edgeInterp=false; }
  bool getEdgeInterp() const { return edgeInterp; }

  unsigned int getNKnots() const { return model.getNKnots(); }
  unsigned int getNSigmas() const { return model.getNSigmas(); }
  unsigned int getNOffsets() const { return model.getNOffsets(); }
  void setKnotPositions(std::vector<double>& knts) 
  { model.setKnotPositions(knts); }
  void setSigmaPositions(std::vector<double>& knts) 
  { model.setSigmaPositions(knts); }
  void setOffsetPositions(std::vector<double>& knts) 
  { model.setOffsetPositions(knts); }
  void setPositions(std::vector<double>& K, std::vector<double>& S,
		    std::vector<double>& O) {
    model.setPositions(K,S,O);
  }

  //Sigma prior
  /*! \brief Activates the sigma prior with width set to value*/
  void setSigmaPrior( double val ) { 
    has_sigma_prior=true; sigma_prior_width = val;}
  /*! \brief De-activate the sigma prior*/
  void unsetSigmaPrior() { has_sigma_prior = false; }

  //CFIRB prior
  /*! \brief Activates the CFIRB prior with the specified values, band 1 */
  void setCFIRBPrior1( double, double );
  /*! \brief De-activated CFIRB prior, band 1 */
  void unsetCFIRBPrior() { has_cfirb_prior1 = false; }
  /*! \brief Activates the CFIRB prior with the specified values, band 2 */
  void setCFIRBPrior2( double, double );
  /*! \brief De-activated CFIRB prior, band 2 */
  void unsetCFIRBPrior2() { has_cfirb_prior2 = false; }

  unsigned int getNBeamSets() const { return nbeamsets; }

  double getLogLike(const paramSet&) const;

  void SendSelf(MPI::Comm&, int dest) const;
  void RecieveCopy(MPI::Comm&, int dest);
};

#endif
