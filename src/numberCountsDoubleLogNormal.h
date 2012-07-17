//numberCountsDoubleLogNormal.h

#ifndef __numberCountsDoubleLogNormal__
#define __numberCountsDoubleLogNormal__

#include<gsl/gsl_errno.h>
#include<gsl/gsl_spline.h>
#include<gsl/gsl_interp.h>
#include<gsl/gsl_integration.h>

#include<numberCountsDouble.h>


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
    \f$\sigma\f$ and \f$\mu\f$ are also strored as splines.

  \ingroup Models
*/
class numberCountsDoubleLogNormal : public numberCountsDouble {
 protected :
  
  //Number counts in band one
  unsigned int nknots; //!< Number of knot positions for counts in band one
  double *knots; //!< Locations of knots
  double *logknots; //!< Log of knot positions
  double *logknotvals; //!< Log values of differential number counts at knots, band 1
  gsl_interp_accel *acc; //!< Spline lookup accelerator
  gsl_spline *splinelog; //!< Spline in log/log space of logknots.

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
  double getSigmaInner(double) const;
  double getOffsetInner(double) const;
  double getNumberCountsInner(double,double) const;

  /*! \brief Integrate spline * flux1^power1 * EXP( const1*mu + 
    const2*sigma^2 ) */
  double splineInt(double power1, double const1, 
		   double const2) const;
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
				       unsigned int);
  numberCountsDoubleLogNormal( const std::vector<double>&, 
			       const std::vector<double>&,
			       const std::vector<double>& );
  numberCountsDoubleLogNormal( unsigned int, const double* const,
			       unsigned int, const double* const,
			       unsigned int, const double* const);
  numberCountsDoubleLogNormal( const numberCountsDoubleLogNormal& other );
  ~numberCountsDoubleLogNormal(); //!< Destructor

  numberCountsDoubleLogNormal& operator=(const numberCountsDoubleLogNormal&);

  virtual bool isValid() const; //!< Are model settings loaded and valid?

  void setNKnots(unsigned int n); //!< Sets number of knots in band 1
  unsigned int getNKnots() const { return nknots; }
  void setKnotPositions( const std::vector<double>& );
  void setKnotPositions( unsigned int, const double* const );
  void setNSigmas(unsigned int);
  unsigned int getNSigmas() const { return nsigmaknots; }
  void setSigmaPositions( const std::vector<double>& );
  void setSigmaPositions( unsigned int, const double* const );
  void setNOffsets(unsigned int);
  unsigned int getNOffsets() const { return noffsetknots; }
  void setOffsetPositions( const std::vector<double>& );
  void setOffsetPositions( unsigned int, const double* const );

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
   double getR(double, double, const doublebeam&, rtype=BEAMALL) const;
   void getR(unsigned int n1,const double* const, unsigned int n2,
	     const double* const, const doublebeam&, double*, 
	     rtype=BEAMALL) const;

   bool writeToStream(std::ostream& os) const; //<! Output

   void SendSelf(MPI::Comm&, int dest) const; //!< Send self
   void RecieveCopy(MPI::Comm&, int src); //!< Recieve
};

//Function to pass to GSL integrator
/*! \brief Evaluates flux1^power1 * exp(const1*mu + const2*sigma^2) dN/dS1 */
double evalPowfNLogNormal(double,void*); 

#endif
