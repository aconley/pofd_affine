#include <iomanip>

#include "../include/hashbar.h"

/*!
  \param[in] MAXLEN Maximum number of marks
  \param[in] NSTEPS Number of total steps expected
  \param[in] MARK Mark to draw
 */
hashbar::hashbar(unsigned int MAXLEN, unsigned int NSTEPS,
		 char MARK) :
  maxlen(MAXLEN), currlen(0), nsteps(NSTEPS), mark(MARK) {}

/*!
  \param[in] os Stream to write to

  Writes an empty progress bar, but without a newline.
*/
void hashbar::initialize(std::ostream& os) {
  os << "   0.0%[";
  for (unsigned int i = 0; i < maxlen; ++i)
    os << " ";
  os << "]\r" << std::flush;
}

/*!
  \param[in] step Step we are on (related to nsteps)
  \param[in] os Stream to write to

  Writes the progress bar without a newline.
*/
void hashbar::update(unsigned int step, std::ostream& os) {
  unsigned int nhash;
  if (step >= nsteps) nhash = maxlen;
  else nhash = step * maxlen / nsteps;
  if (nhash == currlen) return;
  if (nhash == 0) {
    // Empty bar
    initialize(os);
  } else if (nhash < maxlen) {
    // Partially full bar
    float perc = 100.0 * static_cast<float>(step) / nsteps;
    perc = perc > 100.0 ? 100.0 : perc;
    std::streamsize sp = os.precision();
    os << " " << std::setw(5) << std::setprecision(1) 
       << std::fixed << perc << "%[";
    os << std::setprecision(sp);
    for (unsigned int h = 0; h < nhash-1; ++h) os << mark;
    os << ">";
    for (unsigned int h = nhash; h < maxlen; ++h)
      os << " ";
    os << "]\r" << std::flush;
  } else {
    // Full bar, but no newline
    os << " 100.0%[";
    for (unsigned int h = 0; h < maxlen; ++h) os << mark;
    os << "]\r" << std::flush;
  }
  currlen = nhash;
}

/*!
  \param[in] os Stream to write to

  Writes a full progress bar with a newline.  This should
  be the last thing you call.
*/
void hashbar::fill(std::ostream& os) {
  os << " 100.0%[";
  for (unsigned int h = 0; h < maxlen; ++h) os << mark;
  os << "]" << std::endl;
}
