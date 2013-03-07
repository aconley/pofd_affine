#include "../include/numberCountsDouble.h"

std::ostream& operator<<(std::ostream& os, const numberCountsDouble& b) {
  b.writeToStream(os);
  return os;
}
