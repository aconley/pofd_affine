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
  stored internally as natural logarithms.  The parameters are
  the knots
 */
class numberCountsKnotsSpline : public numberCountsKnots {
 private :
  double *logknots; //!< Log of knot positions

  void setNKnots(unsigned int n); //!< Sets number of knots
  gsl_interp_accel *acc; //!< Spline lookup accelerator
  gsl_spline *splinelog; //!< Spline in log/log space

  gsl_integration_workspace *gsl_work; //!< Integration workspace for QAG

  double splineInt(double power) const; //!< Integrate spline * flux^power
  void **varr; //!< Internal evil casting array for integration

  //Internal pos/neg parts
  double getRPos(double,const beam&);
  void getRPos(unsigned int n,const double* const,
	       const beam&,double*);
  double getRNeg(double,const beam&);
  void getRNeg(unsigned int n,const double* const,
	       const beam&,double*);


 public :
  numberCountsKnotsSpline(); //!< Default
  explicit numberCountsKnotsSpline(unsigned int);
  numberCountsKnotsSpline( const std::vector<double>& );
  numberCountsKnotsSpline( unsigned int, const double* const);
  numberCountsKnotsSpline( const numberCountsKnotsSpline& );
  ~numberCountsKnotsSpline(); //!< Destructor

  numberCountsKnotsSpline& operator=(const numberCountsKnotsSpline&);
 
  void setKnotPositions(const std::vector<double>&); //!< Set knot positions
  void setKnotPositions(unsigned int, const double* const); //!< Set knot positions

  void setParams(const paramSet&); //!< Set parameters

  double getNumberCounts(double) const; //!< Evaluates number counts model
  double getNS() const; //!< Return the number of sources with \f$S > S_{min}\f$ per square degree

  double getMeanFluxPerArea() const; //!< Mean flux per unit area (sq deg)
  double getMeanFluxSqPerArea() const; //!< Mean flux squared per unit area (sq deg)
  
  double getR(double,const beam&, rtype=BEAMPOS);
  void getR(unsigned int n,const double* const,
	    const beam&,double*, rtype=BEAMPOS);

  virtual void SendSelf(MPI::Comm&, int dest) const; //!< Send self
  virtual void RecieveCopy(MPI::Comm&, int src); //!< Recieve
};

//Function to pass to GSL integrator
double evalPowfN(double,void*); //!< Evaluates f^pow dN/dS

#endif
