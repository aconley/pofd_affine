//numberCountsKnotsSpline.h

#ifndef __numberCountsKnotsSpline__
#define __numberCountsKnotsSpline__

#include<gsl/gsl_errno.h>
#include<gsl/gsl_spline.h>
#include<gsl/gsl_integration.h>

#include "../include/numberCountsKnots.h"

/*!
  \brief Spline number counts model, 1D case
 
  The user inputs knot values in \f$\log_{10}\f$, but they are
  stored internally as base 2 logarithms.  The parameters are
  the knots

  \ingroup Models
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

 public :
  numberCountsKnotsSpline(); //!< Default
  explicit numberCountsKnotsSpline(unsigned int); //!< Constructor with number of knots
  numberCountsKnotsSpline(const std::vector<float>&); //!< Vector constructor
  numberCountsKnotsSpline(const std::vector<double>&); //!< Vector constructor
  numberCountsKnotsSpline(unsigned int, const float* const); //!< C array constructor
  numberCountsKnotsSpline(unsigned int, const double* const); //!< C array constructor
  numberCountsKnotsSpline(const numberCountsKnotsSpline&); //!< Copy constructor
  virtual ~numberCountsKnotsSpline(); //!< Destructor

  /*! \brief copy operator */
  numberCountsKnotsSpline& operator=(const numberCountsKnotsSpline&);
 
  void setKnotPositions(const std::vector<float>&); //!< Set knot positions
  void setKnotPositions(unsigned int, const float* const); //!< Set knot positions
  void setKnotPositions(const std::vector<double>&); //!< Set knot positions
  void setKnotPositions(unsigned int, const double* const); //!< Set knot positions

  void setParams(const paramSet&); //!< Set parameters

  double getNumberCounts(double) const; //!< Evaluates number counts model
  double getNS() const; //!< Return the number of sources with \f$S > S_{min}\f$ per square degree

  double getFluxPerArea() const; //!< Flux per unit area (sq deg)
  double getFluxSqPerArea() const; //!< Flux squared per unit area (sq deg)
  
  /*! \brief Get R at single flux value */
  double getR(double, const beam&) const;
  /*! \brief Get R, array version */
  void getR(unsigned int n, const double* const,
	    const beam&, double*) const;

  double differenceRegularize(double) const; //!< Tikhonov regularization log Likelihood penalty

  void writeToHDF5Handle(hid_t, bool=false) const; //!< Write to HDF5 handle

  virtual void sendSelf(MPI_Comm, int dest) const; //!< Send self
  virtual void receiveCopy(MPI_Comm, int src); //!< Receive
};

#endif
