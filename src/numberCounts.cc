#include "../include/numberCounts.h"

/*!
  \param[inout] os Stream to write to
  \param[in] b Number counts model
*/
std::ostream& operator<<(std::ostream& os, const numberCounts& b) {
  b.writeToStream(os);
  return os;
}
