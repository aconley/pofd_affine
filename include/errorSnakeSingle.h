#ifndef __errorSnakeSingle__
#define __errorSnakeSingle__

#include<string>
#include<iostream>

#include "../include/numberCountsKnotsSpline.h"

class errorSnakeSingle final {
private:
  static const unsigned int nprob = 3;
  const float plevels[nprob] = {0.683, 0.954, 0.997};
  
  unsigned int npoints;
  float *s; //!< Flux density values, len npoints
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
