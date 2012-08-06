//numberCountsDoubleLogNormal.h

#ifndef __numberCountsDoubleLogNormal__
#define __numberCountsDoubleLogNormal__

#include<gsl/gsl_errno.h>
#include<gsl/gsl_spline.h>
#include<gsl/gsl_interp.h>
#include<gsl/gsl_integration.h>

#include<numberCountsDouble.h>
#include<ran.h>

/*!
  \brief Spline number counts model for 2D, where counts
    are modeled as a spline at one frequency plus a Log-Normal colour
    model.  

    The full expression is
    \f[
       \frac{dN}{dS_1\, dS_2} = \frac{n_1\left(S_1 \right)}{S_1} 
            \mathrm{L} \left( \frac{S_2}{S_1}; \mu\left(S_1\right),
            \sigma\left(S_1\right) \right)
    \f] 
    where \f$n_1\f$ is just the spline as in the one-dimensional case.
    The division by \f$S_1\f$ may seem odd, but assures that the
    distribution marginalized over \f$S_2\f$ is just \f$ n_1 \f$.
    \f$\sigma\f$ and \f$\mu\f$ are also strored as splines. \f$L\f$
    is just the Log-Normal form:
    \f[
      \mathrm{L}\left( x ; \mu, \sigma \right) =
       \frac{1}{\sqrt{2 \pi} \sigma} \frac{1}{x}
       \exp\left[ \frac{ - \left(\log x - \mu\right)^2 }{ \sigma^2 } \right] .
    \f]
    Note that \f$\sigma\f$ and \f$\mu\f$ are not the mean and square root
    of the variance of the actual distribution, but rather the mean
    and square root of the variance of the log quantities.  These
    are related by
    \f[
       \left< S_2/S_1 \right> = \exp \left[ \mu\left(S_1\right) +
             \frac{1}{2}\sigma^2\left(S_1\right) \right]
    \f]
    and
    \f[
       Var[ S_2/S_1 ] = \exp \left( 2 \mu \left( S_1 \right) +
           \sigma^2\left(S_1\right) \right) 
          \left[ \exp \left( \sigma^2\left(S_1\right) \right) - 1 \right]
    \f].
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

  bool knotvals_loaded; //!< Knot information loaded

  //Workspace
  gsl_integration_workspace *gsl_work; //!< Integration workspace for QAG
  mutable unsigned int nRWork; //!< Number of elements in R working arrays
  mutable bool* RWorkValid; //!< Valid entries in R working array
  mutable double* RWork1; //!< R working array 1
  mutable double* RWork2; //!< R working array 2
  mutable double* RWork3; //!< R working array 3
  void setRWorkSize( unsigned int ) const; //!< Controls R working arrays

  //Convenience functions for spline computations without input checking
  double getSigmaInner(double) const; //!< Inner sigma computation
  double getOffsetInner(double) const; //!< Inner offset computation
  double getNumberCountsInner(double,double) const; //!< Inner number counts computation

  /*! \brief Integrate powers of fluxes over number counts */
  double splineInt(double alpha, double beta) const;
  static const unsigned int nvarr; //!< Number of elements in varr
  void **varr; //!< Internal evil casting array for splineInt

  //Internal R computation
  double getRInternal(double,double,const doublebeam&,unsigned int) const;
  void getRInternal(unsigned int, const double* const, 
		    unsigned int, const double* const,
		    const doublebeam&,unsigned int, double*) const;

 public :
  numberCountsDoubleLogNormal(); //!< Default
  explicit numberCountsDoubleLogNormal(unsigned int, unsigned int, 
				       unsigned int); //!< Constructor
  numberCountsDoubleLogNormal( const std::vector<double>&, 
			       const std::vector<double>&,
			       const std::vector<double>& ); //!< Constructor
  /*! \brief Constructor */
  numberCountsDoubleLogNormal( unsigned int, const double* const,
			       unsigned int, const double* const,
			       unsigned int, const double* const);
  /*! \brief Copy constructor */
  numberCountsDoubleLogNormal( const numberCountsDoubleLogNormal& other );
  ~numberCountsDoubleLogNormal(); //!< Destructor

  /*! \brief Copy operator */
  numberCountsDoubleLogNormal& operator=(const numberCountsDoubleLogNormal&);

  virtual bool isValid() const; //!< Are model settings loaded and valid?

  void setNKnots(unsigned int n); //!< Sets number of knots in band 1
  unsigned int getNKnots() const { return nknots; } //!< Number of knots in band 1 spline
  /*! \brief Sets knot positions, band 1, vector version */
  void setKnotPositions( const std::vector<double>& );
  /*! \brief Sets knot positions, band 1, c array version */
  void setKnotPositions( unsigned int, const double* const );
  /*! \brief Sets number of sigma positions for color model */
  void setNSigmas(unsigned int);
  /*! \brief Gets number of sigma knots in band 2 color model */
  unsigned int getNSigmas() const { return nsigmaknots; }
  /*! \brief Sets knot positions for band 2 sigma model, vector version */
  void setSigmaPositions( const std::vector<double>& );
  /*! \brief Sets knot positions for band 2 sigma model, c array version */
  void setSigmaPositions( unsigned int, const double* const );
  /*! \brief Sets number of offset positions for color model */
  void setNOffsets(unsigned int);
  /*! \brief Gets number of offset knots in band 2 color model */
  unsigned int getNOffsets() const { return noffsetknots; }
  /*! \brief Sets knot positions for band 2 offset model, vector version */
  void setOffsetPositions( const std::vector<double>& );
  /*! \brief Sets knot positions for band 2 offset model, c array version */
  void setOffsetPositions( unsigned int, const double* const );
  /*! \brief Gets total number of parameters */
  unsigned int getNParams() const { return nknots+nsigmaknots+noffsetknots; }
  /*! \brief Sets positions of all knots */
  void setPositions(const std::vector<double>&,const std::vector<double>&,
		    const std::vector<double>&); //!< Set all positions
  

  void setParams(const paramSet&); //!< Set parameters
 
  /*! \brief Evaluates Sigma */
   double getSigma(double) const;
  /*! \brief Evaluates Offset */
   double getOffset(double) const;

  /*! \brief Evaluates number counts model */
   double getNumberCounts( double, double ) const;

  /*! \brief Minimum flux model is defined for */
   double getMinFlux(unsigned int) const;

  /*! \brief Maxium flux model is defined for */
   double getMaxFlux(unsigned int) const;

   double getNS() const; //!< Total number of sources per area
   double getMeanFluxPerArea(unsigned int) const; //!< Mean flux per unit area (sq deg)
   double getMeanFluxPowerPerArea(double p1,double p2) const; //!< Integral of s1 to some power and s2 to some power over number counts

   //Routines for getting R
   /*! \brief Get R value at single flux1,flux2 value */
   double getR(double, double, const doublebeam&, rtype=BEAMALL) const;
   /*! \brief Gets R values, array version */
   void getR(unsigned int n1,const double* const, unsigned int n2,
	     const double* const, const doublebeam&, double*, 
	     rtype=BEAMALL) const;

   bool writeToStream(std::ostream& os) const; //<! Output

   void SendSelf(MPI::Comm&, int dest) const; //!< Send self
   void RecieveCopy(MPI::Comm&, int src); //!< Recieve
};

//Function to pass to GSL integrator
/*! \brief Evaluates flux1^power1 * exp(const1*mu + const2*sigma^2) dN/dS1 */
double evalPowfNDoubleLogNormal(double,void*); 

//////////////////////////////

/*!
  \brief A class to read in model specifications from init files
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
  double* knotpos; //!< Positions of knots (all)
  double* knotval; //!< Initial value center for knot value (all)

  //These are optional
  bool has_sigma; //!< Has initial value sigma
  double* sigma; //!< Sigma value for initial positions
  bool has_lower_limits; //!< Has some lower limit information
  bool* has_lowlim; //!< Knots have lower limit
  double* lowlim; //!< Value of lower limit
  bool has_upper_limits; //!< Has some upper limit information
  bool* has_uplim; //!< Knots have upper limit
  double* uplim; //!< Value of upper limit

  mutable ran rangen; //!< Random number generator

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

  /*! \brief Set seed of random number generator */
  void setSeed( unsigned long long int seed ) const { rangen.setSeed(seed); }

  void getKnotPos(std::vector<double>&) const; //!< Gets the knot positions for band 1
  void getKnotVals(std::vector<double>&) const; //!< Gets the knot values for band 1
  void getSigmaPos(std::vector<double>&) const; //!< Gets the knot positions for color model sigma
  void getOffsetPos(std::vector<double>&) const; //!< Gets the knot positions for color model offset
  void getModelPositions(numberCountsDoubleLogNormal&) const; //!< Sets knot locations in model for all model components
  void getParams(paramSet& p) const; //!< Sets knot values to central values
  void generateRandomKnotValues(paramSet& p) const; //!< Seed knot values
  
  double getKnotSigma(unsigned int) const; //!< Get knot sigma
  bool isKnotFixed(unsigned int) const; //!< Is a knot fixed?

  bool knotHasLowerLimit(unsigned int) const; //!< Does knot have a lower limit
  double getLowerLimit(unsigned int) const; //!< Get knot lower limit

  bool knotHasUpperLimit(unsigned int) const; //!< Does knot have a lower limit
  double getUpperLimit(unsigned int) const; //!< Get knot lower limit

  bool isValid(const paramSet&) const; //!< Checks if parameters are within allowed ranges

  void SendSelf(MPI::Comm&, int dest) const; //!< Send self
  void RecieveCopy(MPI::Comm&, int src); //!< Recieve
};

#endif
