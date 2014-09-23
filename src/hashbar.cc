#include <iomanip>

#include "../include/hashbar.h"

hashbar::hashbar(unsigned int MAXLEN, unsigned int NSTEPS) :
  maxlen(MAXLEN), currlen(0), nsteps(NSTEPS) {}

void hashbar::update(unsigned int step, std::ostream& os) {
  unsigned int nhash;
  if (step >= nsteps) nhash = maxlen;
  else nhash = step * maxlen / nsteps;
  if (nhash == currlen) return;
  if (nhash == 0) {
    // Empty bar
    os << "  0.0%[";
    for (unsigned int i = 0; i < maxlen; ++i)
      os << " ";
    os << "]\r" << std::flush;
  } else if (nhash < maxlen) {
    // Partially full bar
    float perc = 100.0 * static_cast<float>(step) / nsteps;
    perc = perc > 100.0 ? 100.0 : perc;
    os << std::setw(5) << std::setprecision(1) << std::fixed << perc << "%[";
    for (unsigned int h = 0; h < nhash-1; ++h) os << "=";
    os << ">";
    for (unsigned int h = nhash; h < maxlen; ++h)
      os << " ";
    os << "]\r" << std::flush;
  } else {
    // Full bar, but no newline
    os << "100.0%[";
    for (unsigned int h = 0; h < maxlen; ++h) os << "=";
    os << "]\r" << std::flush;
  }
  currlen = nhash;
}

void hashbar::fill(std::ostream& os) {
  os << "100.0%[";
  for (unsigned int h = 0; h < maxlen; ++h) os << "=";
  os << "]" << std::endl;
}
