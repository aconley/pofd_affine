//utility.h

#ifndef __utility__
#define __utility__

#include<vector>
#include<string>

/*!
  \brief Utility functions
*/
namespace utility {

  /*! \brief Output file types */
  enum outfiletype { UNKNOWN=0, TXT=1, FITS=2, HDF5=3 };

  /*! \brief Determine file type from extension */
  outfiletype getOutputFileType(const std::string& str);

  /*! \brief Break input string up into words breaking on spaces */
  void stringwords(const std::string &ins,
                   std::vector<std::string> &words);

  /*! \brief Break input string up into words, also breaking on equals sign */
  void stringwords_eq(const std::string &ins,
                      std::vector<std::string> &words);

  /*! \brief Checks to see if a string should be interpreted as true */
  bool string_true(const std::string &ins);

  /*! \brief Binary search for last value less than specified amount */
  unsigned int binary_search_lt(double value, double* data, unsigned int ndata);

  /*! \brief Binary search for last value less than specified amount */
  unsigned int binary_search_lt(double value, const std::vector<double>&);

  /*! \brief Binary search for first value greater than specified amount */
  unsigned int binary_search_gt(double value, double* data, unsigned int ndata);

  /*! \brief Binary search for first value greater than specified amount */
  unsigned int binary_search_gt(double value,std::vector<double> data);

  /*! \brief Binary search for last value less than or equal to 
     specified amount */
  unsigned int binary_search_lte(double value, double* data, 
                                 unsigned int ndata);

  /*! \brief Binary search for last value greater than specified amount in reverse array*/
  int binary_search_rev(double value, double* data, 
                        unsigned int ndata);

  /*! \brief Binary logarithm (\f$log_{2}\f$) of 32 bit integer */
  unsigned int log2(unsigned int);

  /*! \brief Log of factorial of argument */
  double logfactorial(double);
}

#endif
