#ifndef __numberCountsDoubleLogNormalError__
#define __numberCountsDoubleLogNormalError__

#include "../include/statsAccumulator.h"
#include "../include/numberCountsDoubleLogNormal.h"

/*!
  \brief Class for handling uncertainty information for fits using
  the numberCountsDoubleLogNormal model.

  Builds statistics about a fit from a fit output file, including
    - Mean fit
    - Best fit
    - Uncertainties at knot values
    - Error snakes

  \ingroup models
*/
class numberCountsDoubleLogNormalStats final : public statsAccumulator {
private:

  // Stuff having to do with values and uncertainties at the knots
  unsigned int nparams; //!< Number of model parameters
  std::vector<std::string> paramNames; //!< Names of parameters
  unsigned int nknots; //!< Number of knots
  unsigned int nsigmas; //!< Number of sigmas
  unsigned int noffsets; //!< Number of offsets
  float *knotPositions; //!< 1D model Knot positions, len nknots
  float *sigmaKnotPositions; //!< Sigma knot positions
  float *offsetKnotPositions; //!< Offset knot positions
  double bestLike; //!< Likelihood of best fit
  float *bestParams; //!< Best fit values of params
  float *meanParams; //!< Mean parameter values
  float *medParams; //!< Median parameter values
  float *uncertaintyParams; //!< Param uncert, nprob x nparams x 2
  
  // Error snake stuff
  unsigned int npoints;  //!< Number of points snake is evaulated at
  // 1D model.  These are expressed in log10 -- so the mean is the mean
  //   log10 counts, etc.
  float *s1D; //!< Flux densities of error snake for 1D model, len npoints
  float *mean1D; //!< Mean 1D model value, len npoints (mean at each point)
  float *median1D; //!< Median 1D model value, len npoints
  float *snake1D; //!< Actual 1D error snake, nprob by npoints by 2, flattened

  // Color model.  There is an option to compute parameters in terms
  //   of f2 / f1 instead of the nominal model parameters.
  //  If so, then sigma is sigma(f2/f1), and offset is <f2/f1>
  bool areRatios; //!< True if f2 / f1 ratios are tabulated
  float *sSigma; //!< Flux densities of sigma error snake, len npoints
  float *meanSigma; //!< Mean color sigma
  float *medianSigma; //!< Median color sigma
  float *snakeSigma; //!< Sigma snake
  float *sOffset; //!< Flux densities of offset error snake, len npoints
  float *meanOffset; //!< Mean color offset
  float *medianOffset; //!< Median color offset
  float *snakeOffset; //!< Sigma snake
 
  void setNParams(unsigned int); //!< Set number of params

  mutable numberCountsDoubleLogNormal model; //!< For doing interpolation
public:
  numberCountsDoubleLogNormalStats(unsigned int, bool=true);
  ~numberCountsDoubleLogNormalStats();

  // Main processing function
  void build(const std::string&, unsigned int=5, 
             bool=true) override;

  // Write results
  void writeAsHDF5(const std::string&) const override;
};

#endif
