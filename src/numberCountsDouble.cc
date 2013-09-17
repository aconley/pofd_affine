#include "../include/numberCountsDouble.h"

/*!
  \param[inout] os Stream to write to
  \param[in] b Number counts model to write
*/
std::ostream& operator<<(std::ostream& os, const numberCountsDouble& b) {
  b.writeToStream(os);
  return os;
}
