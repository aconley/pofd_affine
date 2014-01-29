//numberCountsDoubleLogNormal.h

#ifndef __numberCountsDoubleLogNormal__
#define __numberCountsDoubleLogNormal__

#include<gsl/gsl_errno.h>
#include<gsl/gsl_spline.h>
#include<gsl/gsl_interp.h>
#include<gsl/gsl_integration.h>

#include "../include/global_settings.h"
#include "../include/numberCountsDouble.h"
#include "../include/ran.h"

/*!
  \brief Spline number counts model for 2D case.

  The counts are modeled as a spline at one frequency times a Log-Normal 
  colour (flux ratio) model.  

  The full expression is
  \f[
    \frac{dN}{dS_1\, dS_2} = \frac{n_1\left(S_1 \right)}{S_1} 
    \mathrm{L} \left( \frac{S_2}{S_1}; \mu\left(S_1\right),
    \sigma\left(S_1\right) \right)
  \f] 
  where \f$n_1\f$ is just the spline as in the one-dimensional case
  as implemented in numberCountsKnotsSpline.
  The division by \f$S_1\f$ may seem odd, but assures that the
  distribution marginalized over \f$S_2\f$ is just \f$ n_1 \f$.
  \f$L\f$ is just the Log-Normal distribution:
  \f[
    \mathrm{L}\left( x ; \mu, \sigma \right) =
    \frac{1}{\sqrt{2 \pi} \sigma} \frac{1}{x}
    \exp\left[ \frac{ - \left(\log x - \mu\right)^2 }{ \sigma^2 } \right] .
  \f]
  \f$\sigma\f$ and \f$\mu\f$ are stored as splines as functions of
  \f$S_1\f$.  Note that \f$\sigma\f$ and \f$\mu\f$ (the actual parameters 
  of the model)
  are not the mean and square root of the variance of \f$S_2 / S_1\f$,
  but rather the mean and square root of the variance of the 
  (natural) log quantities.  These are related by
  \f[
    \left< S_2/S_1 \right> = \exp \left[ \mu\left(S_1\right) +
       \frac{1}{2}\sigma^2\left(S_1\right) \right]
  \f]
  and
  \f[
    Var[ S_2/S_1 ] = \exp \left( 2 \mu \left( S_1 \right) +
       \sigma^2\left(S_1\right) \right) 
       \left[ \exp \left( \sigma^2\left(S_1\right) \right) - 1 \right] .
  \f]
  Note that the number counts require that \f$S_2/S_1 > 0\f$.
  We also explicitly require that \f$S_1 > 0\f$.  Both are justified
  on physical grounds.
    
  \ingroup Models
*/
class numberCountsDoubleLogNormal : public numberCountsDouble {
 protected :
  
  //Number counts in band one
  unsigned int nknots; //!< Number of knot positions for counts in band one
  double *knots; //!< Locations of knots
  double *logknots; //!< Log of knot positions (stored as log 2)
  double *logknotvals; //!< Log 2 values of differential number counts at knots, band 1
  gsl_interp_accel *acc; //!< Spline lookup accelerator
  gsl_spline *splinelog; //!< Spline in log2/log2 space of values

  //LogNormal for band 2
  unsigned int nsigmaknots; //!< Number of knot positions in sigma
  double *sigmaknots; //!< Positions of sigma knots
  double *sigmavals; //!< Sigma values at positions of knots
  gsl_interp_accel *accsigma; //!< Sigma accelerator
  gsl_interp *sigmainterp; //!< Sigma interpolator
  unsigned int noffsetknots; //!< Number of knot positions in offset
  double *offsetknots; //!< Positions of offset knots
  double *offsetvals; //!< Offset values at positions of knots
  gsl_interp_accel *accoffset; //!< Offset accelerator
  gsl_interp *offsetinterp; //!< Offset interpolator

  //Model validity stuff
  mutable bool knots_valid; //!< Knots ready
  mutable bool sigmas_valid; //!< Sigmas ready
  mutable bool offsets_valid; //!< Offsets ready
  bool knotpos_loaded; //!< Knot positions loaded (not validity checked)
  bool sigmapos_loaded; //!< Sigma positions loaded (not validity checked)
  bool offsetpos_loaded; //!< Offset positions loaded (not validity checked)
  bool knotvals_loaded; //!< Knot values loaded (not validity checked)
  bool sigmavals_loaded; //!< Sigma values loaded (not validity checked)
  bool offsetvals_loaded; //!< Offset values loaded (not validity checked)
  void checkKnotsValid() const; //!< Set knots_valid by checking
  void checkSigmasValid() const; //!< Set sigmas_valid by checking
  void checkOffsetsValid() const; //!< Set offsets_valid by checking

  //Workspace
  gsl_integration_workspace *gsl_work; //!< Integration workspace for QAG
  mutable unsigned int nRWork; //!< Number of elements in R working arrays
  mutable bool* RWorkValid; //!< Has R work been set for a given index
  mutable double* RWork; //!< R working array; for pre-computing R bits
  void setRWorkSize(unsigned int) const; //!< Controls R working arrays

  //Convenience functions for spline computations without input checking
  double getSigmaInner(double) const; //!< Inner sigma computation
  double getOffsetInner(double) const; //!< Inner offset computation
  double getNumberCountsInner(double,double) const; //!< Inner number counts computation

  /*! \brief Get range over which R is expected to be nonzero, unchecked version */
  std::pair<dblpair, dblpair> getRRangeInternal(const doublebeam&) const 
    throw(affineExcept);

  /*! \brief Integrate powers of fluxes over number counts */
  double splineInt(double alpha, double beta) const;
  static const unsigned int nvarr; //!< Number of elements in varr
  void **varr; //!< Internal evil casting array for splineInt

  bool isValidLoaded() const; //!< Are loaded values valid?

 public :
  numberCountsDoubleLogNormal(); //!< Default
  explicit numberCountsDoubleLogNormal(unsigned int, unsigned int, 
				       unsigned int); //!< Constructor

  /*! \brief Constructor */
  numberCountsDoubleLogNormal(const std::vector<float>&, 
			      const std::vector<float>&,
			      const std::vector<float>&); //!< Constructor
  /*! \brief Constructor */
  numberCountsDoubleLogNormal(const std::vector<double>&, 
			      const std::vector<double>&,
			      const std::vector<double>&); //!< Constructor
  /*! \brief Constructor */
  numberCountsDoubleLogNormal(unsigned int, const float* const,
			      unsigned int, const float* const,
			      unsigned int, const float* const);
  /*! \brief Constructor */
  numberCountsDoubleLogNormal(unsigned int, const double* const,
			      unsigned int, const double* const,
			      unsigned int, const double* const);
  /*! \brief Copy constructor */
  numberCountsDoubleLogNormal(const numberCountsDoubleLogNormal& other);
  ~numberCountsDoubleLogNormal(); //!< Destructor

  /*! \brief Copy operator */
  numberCountsDoubleLogNormal& operator=(const numberCountsDoubleLogNormal&);

  bool isValid() const; //!< Are model settings loaded and valid?

  void setNKnots(unsigned int n); //!< Sets number of knots in band 1
  unsigned int getNKnots() const { return nknots; } //!< Number of knots in band 1 spline
  /*! \brief Sets knot positions, band 1, vector version */
  void setKnotPositions(const std::vector<float>&);
  /*! \brief Sets knot positions, band 1, vector version */
  void setKnotPositions(const std::vector<double>&);
  /*! \brief Sets knot positions, band 1, c array version */
  void setKnotPositions(unsigned int, const float* const);
  /*! \brief Sets knot positions, band 1, c array version */
  void setKnotPositions(unsigned int, const double* const);
  /*! \brief Sets number of sigma positions for color model */
  void setNSigmas(unsigned int);
  /*! \brief Gets number of sigma knots in band 2 color model */
  unsigned int getNSigmas() const { return nsigmaknots; }
  /*! \brief Sets knot positions for band 2 sigma model, vector version */
  void setSigmaPositions(const std::vector<float>&);
  /*! \brief Sets knot positions for band 2 sigma model, vector version */
  void setSigmaPositions(const std::vector<double>&);
  /*! \brief Sets knot positions for band 2 sigma model, c array version */
  void setSigmaPositions(unsigned int, const float* const);
  /*! \brief Sets knot positions for band 2 sigma model, c array version */
  void setSigmaPositions(unsigned int, const double* const);
  /*! \brief Sets number of offset positions for color model */
  void setNOffsets(unsigned int);
  /*! \brief Gets number of offset knots in band 2 color model */
  unsigned int getNOffsets() const { return noffsetknots; }
  /*! \brief Sets knot positions for band 2 offset model, vector version */
  void setOffsetPositions(const std::vector<float>&);
  /*! \brief Sets knot positions for band 2 offset model, vector version */
  void setOffsetPositions(const std::vector<double>&);
  /*! \brief Sets knot positions for band 2 offset model, c array version */
  void setOffsetPositions(unsigned int, const float* const);
  /*! \brief Sets knot positions for band 2 offset model, c array version */
  void setOffsetPositions(unsigned int, const double* const);

  //Get individual positions
  double getKnotPosition(unsigned int i) const; //!< Get band 1 model knot position (unchecked)
  double getSigmaPosition(unsigned int i) const; //!< Get sigma knot position (unchecked)
  double getOffsetPosition(unsigned int i) const; //!< Get offset knot position (unchecked)

  /*! \brief Gets total number of parameters */
  unsigned int getNParams() const { return nknots+nsigmaknots+noffsetknots; }
  /*! \brief Gets positions of all knots */
  void getPositions(std::vector<double>&, std::vector<double>&,
		    std::vector<double>&) const;
  /*! \brief Sets positions of all knots */
  void setPositions(const std::vector<float>&,const std::vector<float>&,
		    const std::vector<float>&); 
  /*! \brief Sets positions of all knots */
  void setPositions(const std::vector<double>&,const std::vector<double>&,
		    const std::vector<double>&);
  
  void getParams(paramSet&) const; //!< Get parameters
  void setParams(const paramSet&); //!< Set parameters
 
  /*! \brief Evaluates Sigma */
  double getSigma(double) const;
  /*! \brief Evaluates Offset */
  double getOffset(double) const;
  
  /*! \brief Evaluates number counts model */
  double getNumberCounts(double, double) const;
  
  /*! \brief Minimum flux model is defined for */
  dblpair getMinFlux() const;
  
  /*! \brief Maxium flux model is defined for */
  dblpair getMaxFlux() const;
  
  /*! \brief Get range over which R is expected to be nonzero */
  std::pair<dblpair, dblpair> getRRange(const doublebeam&) const 
    throw(affineExcept);

  double getNS() const; //!< Total number of sources per area
  double getFluxPerArea(unsigned int) const; //!< Flux per unit area (sq deg)
  double getFluxSqPerArea(unsigned int) const; //!< Flux^2 per unit area (sq deg)
  double getFluxPowerPerArea(double p1, double p2) const; //!< Integral of s1 to some power and s2 to some power over number counts
  
  //Routines for getting R
  /*! \brief Get R value at single flux1,flux2 value */
  double getR(double, double, const doublebeam&) const;
  /*! \brief Gets R values, array version */
  void getR(unsigned int n1, const double* const, unsigned int n2,
	    const double* const, const doublebeam&, double*) const;
  
  void writeToHDF5Handle(hid_t) const; //!< Write to HDF5 handle
  bool writeToStream(std::ostream& os) const; //<! Output to stream
  
  void sendSelf(MPI_Comm, int dest) const; //!< Send self
  void recieveCopy(MPI_Comm, int src); //!< Recieve
};

//////////////////////////////

/*!
  \brief A class to read in model specifications from init files, 2D case

  \ingroup Models
 */
class initFileDoubleLogNormal {
 private:
  //These are all c arrays rather than vectors for ease of compatability
  // with MPI.  We also store everything as one vector

  //The first two (knotpos, knotval) are required
  unsigned int nknots; //!< Number of knots in band 1 number counts
  unsigned int nsigmas; //!< Number of knots in color sigma
  unsigned int noffsets; //!< Number of knots in color offset
  unsigned int sigmaidx; //!< Index of first sigma value
  unsigned int offsetidx; //!< Index of first offset value
  double* knotpos; //!< Positions of knots (all, including sigma and offsets)
  double* knotval; //!< Initial value center for knot value (all, inc. sigmas and offsets)

  //These are optional
  bool has_range; //!< Has initial value range
  double* range; //!< Range for initial positions
  bool has_lower_limits; //!< Has some lower limit information
  bool* has_lowlim; //!< Knots have lower limit
  double* lowlim; //!< Value of lower limit
  bool has_upper_limits; //!< Has some upper limit information
  bool* has_uplim; //!< Knots have upper limit
  double* uplim; //!< Value of upper limit

  void checkLimitsDontCross() const; //!< Make sure limits make sense, throw if not
  void checkRange() const; //!< Make sure param ranges make sense, throw if not

 public:
  initFileDoubleLogNormal(); //!< Basic constructor
  /*! \brief Constructor with file read */
  initFileDoubleLogNormal(const std::string&, bool=false, bool=false);
  ~initFileDoubleLogNormal(); //!< Destructor

  unsigned int getNKnots() const { return nknots; } //!< Get number of knots in band 1
  unsigned int getNSigmas() const { return nsigmas; } //!< Get number of knots in color model sigma
  unsigned int getNOffsets() const { return noffsets; } //!< Get number of knots in color model offset
  unsigned int getNTot() const { return nknots+nsigmas+noffsets; } //!< Get total number of model knots

  std::pair<double,double> getKnot(unsigned int idx) const; //!< Get band 1 knot pos and value
  std::pair<double,double> getSigma(unsigned int idx) const; //!< Get color sigma pos and value
  std::pair<double,double> getOffset(unsigned int idx) const; //!< Get color offset pos and value

  void readFile(const std::string&, bool=false, bool=false); //!< Read file

  void getKnotPos(std::vector<double>&) const; //!< Gets the knot positions for band 1
  void getKnotVals(std::vector<double>&) const; //!< Gets the knot values for band 1
  void getSigmaPos(std::vector<double>&) const; //!< Gets the knot positions for color model sigma
  void getSigmaVals(std::vector<double>&) const; //!< Gets the knot values for color model sigma
  void getOffsetPos(std::vector<double>&) const; //!< Gets the knot positions for color model offset
  void getOffsetVals(std::vector<double>&) const; //!< Gets the knot values for color model offset

  void getModelPositions(numberCountsDoubleLogNormal&) const; //!< Sets knot locations in model for all model components
  void getParams(paramSet& p) const; //!< Sets param values to central values
  void generateRandomKnotValues(ran&, paramSet& pnew) const; //!< Seed knot values
  void generateRandomKnotValues(ran&, paramSet& pnew, const paramSet& pcen) const; //!< Seed knot values

  double getKnotRange(unsigned int) const; //!< Get knot range
  bool isKnotFixed(unsigned int) const; //!< Is a knot fixed?

  bool knotHasLowerLimit(unsigned int) const; //!< Does knot have a lower limit
  double getLowerLimit(unsigned int) const; //!< Get knot lower limit

  bool knotHasUpperLimit(unsigned int) const; //!< Does knot have a lower limit
  double getUpperLimit(unsigned int) const; //!< Get knot lower limit

  bool isValid(const paramSet&) const; //!< Checks if parameters are within allowed ranges

  void writeToHDF5Handle(hid_t) const; //!< Write to HDF5 handle

  void sendSelf(MPI_Comm, int dest) const; //!< Send self
  void recieveCopy(MPI_Comm, int src); //!< Recieve
};

#endif
