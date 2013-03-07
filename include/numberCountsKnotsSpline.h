//numberCountsKnotsSpline.h

#ifndef __numberCountsKnotsSpline__
#define __numberCountsKnotsSpline__

#include<gsl/gsl_errno.h>
#include<gsl/gsl_spline.h>
#include<gsl/gsl_integration.h>

#include<numberCountsKnots.h>

/*!
  \brief Spline number counts model
  \ingroup Models
 
  The user inputs knot values in \f$\log_{10}\f$, but they are
  stored internally as base 2 logarithms.  The parameters are
  the knots
 */
class numberCountsKnotsSpline : public numberCountsKnots {
 private :
  double *logknots; //!< Log2 of knot positions

  void setNKnots(unsigned int n); //!< Sets number of knots
  gsl_interp_accel *acc; //!< Spline lookup accelerator
  gsl_spline *splinelog; //!< Spline in log/log space

  gsl_integration_workspace *gsl_work; //!< Integration workspace for QAG

  /*! \brief Integrates flux^power over number counts */
  double splineInt(double power) const; 
  void **varr; //!< Internal evil casting array for integration

  //Internal pos/neg parts
  /*! \brief Get R for positive beam, single flux version */
  double getRPos(double,const beam&) const;
  /*! \brief Get R for positive beam, array of fluxes version */
  void getRPos(unsigned int n,const double* const,
	       const beam&,double*) const;
  /*! \brief Get R for negative beam, single flux version */
  double getRNeg(double,const beam&) const;
  /*! \brief Get R for negative beam, array of fluxes version */
  void getRNeg(unsigned int n,const double* const,
	       const beam&,double*) const;


 public :
  numberCountsKnotsSpline(); //!< Default
  explicit numberCountsKnotsSpline(unsigned int); //!< Constructor with number of knots
  numberCountsKnotsSpline( const std::vector<double>& ); //!< Vector constructor
  numberCountsKnotsSpline( unsigned int, const double* const); //!< C array constructor
  numberCountsKnotsSpline( const numberCountsKnotsSpline& ); //!< Copy constructor
  ~numberCountsKnotsSpline(); //!< Destructor

  /*! \brief copy operator */
  numberCountsKnotsSpline& operator=(const numberCountsKnotsSpline&);
 
  void setKnotPositions(const std::vector<double>&); //!< Set knot positions
  void setKnotPositions(unsigned int, const double* const); //!< Set knot positions

  void setParams(const paramSet&); //!< Set parameters

  double getNumberCounts(double) const; //!< Evaluates number counts model
  double getNS() const; //!< Return the number of sources with \f$S > S_{min}\f$ per square degree

  double getMeanFluxPerArea() const; //!< Mean flux per unit area (sq deg)
  double getMeanFluxSqPerArea() const; //!< Mean flux squared per unit area (sq deg)
  
  /*! \brief Get R at single flux value */
  double getR(double,const beam&, rtype=BEAMPOS) const;
  /*! \brief Get R, array version */
  void getR(unsigned int n,const double* const,
	    const beam&,double*, rtype=BEAMPOS) const;

  virtual void sendSelf(MPI::Comm&, int dest) const; //!< Send self
  virtual void recieveCopy(MPI::Comm&, int src); //!< Recieve
};

#endif
