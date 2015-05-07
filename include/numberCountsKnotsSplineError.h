#ifndef __numberCountsKnotsSplineError__
#define __numberCountsKnotsSplineError__

#include<string>
#include<iostream>

#include "../include/numberCountsKnotsSpline.h"

/*!
  \brief Class for handling uncertainty information for fits using
  the numberCountsKnotsSpline model.

  Builds statistics about a fit from a fit output file, including
    - Mean fit
    - Best fit
    - Uncertainties at knot values
    - Error snakes

  \ingroup models
*/
class numberCountsKnotsSplineError final {
private:
  static const unsigned int nprob = 3;
  const float plevels[nprob] = {0.683, 0.954, 0.997};

  // Stuff having to do with values and uncertainties at the knots
  unsigned int nparams; //!< Number of model parameters
  std::vector<std::string> paramNames; //!< Names of parameters
  unsigned int nknots; //!< Number of knots
  float *knotPositions; //!< Knot positions, len nknots
  double bestLike; //!< Likelihood of best fit
  float *bestParams; //!< Best fit values of params
  float *meanParams; //!< Mean parameter values
  float *medParams; //!< Median parameter values
  float *uncertaintyParams; //!< Param uncert, nprob x nparams x 2
  
  // Error snake stuff
  unsigned int npoints;  //!< Number of points snake is evaulated at
  float *s; //!< Flux densities of error snake, len npoints
  // These are expressed in log10 -- so the mean is the mean
  //   log10 counts, etc.
  float *mean; //!< Mean model value, len npoints (mean at each point)
  float *median; //!< Median model value, len npoints
  float *snake; //!< Actual error snake, nprob by npoints by 2, flattened

  void setNParams(unsigned int); //!< Set number of params

  mutable numberCountsKnotsSpline model; //!< For doing interpolation
public:
  numberCountsKnotsSplineError(unsigned int);
  ~numberCountsKnotsSplineError();

  // Main processing function
  void build(const std::string&, unsigned int=5, bool=true);

  // Write results
  void writeAsHDF5(const std::string&) const;
};

#endif
