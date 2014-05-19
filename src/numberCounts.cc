#include "../include/numberCounts.h"

/*!
  \param[in] alpha Regularization multiplier
  \returns log Likelhood penalty

  The default is to return zero, but subclasses can implement
  this to add a 'smoothing' penalty to the likelihood.  This
  is intended to be used as Tikhonov regularization with a 
  difference-like operator
 */
double numberCounts::differenceRegularize(double alpha) const {
  return 0.0;
}

/*!
  \param[inout] os Stream to write to
  \param[in] b Number counts model
*/
std::ostream& operator<<(std::ostream& os, const numberCounts& b) {
  b.writeToStream(os);
  return os;
}
