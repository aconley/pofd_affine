#include<numberCounts.h>

std::ostream& operator<<(std::ostream& os, const numberCounts& b) {
  b.writeToStream(os);
  return os;
}
