// hashbar.h

#ifndef __hashbar__
#define __hashbar__

#include<ostream>

/*!
  \brief Class for drawing hash progress bars
*/
class hashbar {
 private:
  unsigned int maxlen; //!< Maximum length
  unsigned int currlen; //!< Current number of hashes
  unsigned int nsteps; //!< Number of progress steps.
  char mark; //!< Type of mark to make
 public:
  hashbar(unsigned int, unsigned int, char='#'); //!< Constructor
  
  void initialize(std::ostream&); //!< Make empty progress bar, no newline
  void update(unsigned int, std::ostream&); //<! Update progress bar
  void fill(std::ostream&); //!< Fill progress bar and print newline
};

#endif
