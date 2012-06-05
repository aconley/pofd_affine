//ran.h

#ifndef __ran__
#define __ran__

/*!
  \brief Random number generator

  From Numerical Recipes 3
*/

class ran {
 private:
  unsigned long long int u,v,w;
 public:
  ran(unsigned long long int seed=10214L) { setSeed(seed); }
  void setSeed(unsigned long long int);
  unsigned long long int int64();
  double doub() { return 5.42101086242752217E-20 * int64(); }
  unsigned int int32() { return (unsigned int)int64(); }
  unsigned int selectFromRange(unsigned int, unsigned int);
  double gauss(); //!< Generate a Gaussian random number
};

#endif
