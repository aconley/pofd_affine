//ran.h

#ifndef __ran__
#define __ran__

/*!
  \brief Random number generator

  From Numerical Recipes 3
*/

class ran {
 private:
  unsigned long long int u; //!< Internal state variable
  unsigned long long int v; //!< Internal state variable
  unsigned long long int w; //!< Internal state variable
 public:
  /*! \brief Constructor */
  explicit ran(unsigned long long int seed=10214L) { setSeed(seed); }
  void setSeed(unsigned long long int); //!< Set the seed
  unsigned long long int int64(); //!< Get uniform 64 bit integer
  /*! \brief Get uniform float */
  float flt() { return 5.42101086242752217E-20 * int64(); }
  /*! \brief Get uniform double */
  double doub() { return 5.42101086242752217E-20 * int64(); }
  /*! \brief Get uniform 32 bit integer */
  unsigned int int32() { return (unsigned int)int64(); }
  /*! \brief Select integer from range */
  unsigned int selectFromRange(unsigned int, unsigned int);
  /*! \brief Get unit variance, zero mean Gaussian deviate */
  double gauss(); //!< Generate a Gaussian random number
};

#endif
