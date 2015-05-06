#ifndef __errorSnakeSingle__
#define __errorSnakeSingle__

#include<string>
#include<iostream>

#include "../include/numberCountsKnotsSpline.h"

class errorSnakeSingle final {
private:
  static const unsigned int nprob = 3;
  const float plevels[nprob] = {0.683, 0.954, 0.997};

  unsigned int nknots; //!< Number of knots in model
  float *knots; //!< Location of knots, len nknots
  float *uncertainty; //!< Uncertainty at knot points, nprob by nknots by 2
  
  unsigned int npoints;
  float *s; //!< Flux density values, len npoints
  // These are expressed in log10 point -- so the mean is the mean
  //   log10 counts.
  float *mean; //!< Mean model value, len npoints (mean at each point)
  float *median; //!< Median model value, len npoints
  float *snake; //!< Actual error snake, nprob by npoints by 2, flattened

  mutable numberCountsKnotsSpline model; //!< For doing interpolation
public:
  errorSnakeSingle(unsigned int);
  ~errorSnakeSingle();

  // Main processing function
  void build(const std::string&, unsigned int, bool=true);

  // Write results
  void writeAsHDF5(const std::string&) const;
};

#endif
